//! Packed sparse LU built from the mutable elimination kernel.
//!
//! The factorization records `P * B * Q = L * U`. L is stored by pivot
//! column without its unit diagonal; U stores its diagonal separately and
//! off-diagonal entries by pivot row. Original row/column IDs are retained in
//! the packed streams and mapped to pivot positions by compact lookup arrays.

const std = @import("std");
const sparse_basis = @import("sparse_basis.zig");
const sparse_kernel = @import("sparse_kernel.zig");
const sparse_symbolic = @import("sparse_symbolic.zig");
const sparse_ft = @import("sparse_ft.zig");

pub const SparseLuError = sparse_kernel.KernelError || sparse_symbolic.SymbolicError || sparse_ft.FtError || error{ DimensionMismatch, NumericalFailure };

pub const PivotTrace = struct { rows: []const u32, columns: []const u32 };
pub const OrderingStrategy = enum { automatic, dod_markowitz, highs_kernel };

pub const SparseLU = struct {
    allocator: std.mem.Allocator,
    kernel: sparse_kernel.MutableSparseKernel,
    symbolic: sparse_symbolic.SymbolicWorkspace,
    ft: sparse_ft.SparseForrestTomlin,
    ft_ready: bool = false,
    dimension_capacity: usize = 0,
    factor_capacity: usize = 0,
    dimension: usize = 0,
    l_nonzeros: usize = 0,
    u_nonzeros: usize = 0,
    inserted_fill: usize = 0,
    active_count: usize = 0,
    hyper_views_ready: bool = false,
    ordering_strategy: OrderingStrategy = .automatic,
    selected_ordering: OrderingStrategy = .dod_markowitz,
    peeled_pivots: usize = 0,
    kernel_dimension: usize = 0,
    kernel_nonzeros: usize = 0,
    kernel_maximum_row_count: u32 = 0,
    kernel_maximum_column_count: u32 = 0,
    /// Pivots accepted from the beginning of a validated previous trace.
    trace_replayed_pivots: usize = 0,
    /// Pivots reordered after the first trace validation failure.
    trace_repaired_pivots: usize = 0,

    pivot_rows: []u32 = &.{},
    pivot_columns: []u32 = &.{},
    row_position: []u32 = &.{},
    column_position: []u32 = &.{},
    pivot_values: []f64 = &.{},
    l_starts: []usize = &.{},
    u_starts: []usize = &.{},
    l_rows: []u32 = &.{},
    l_values: []f64 = &.{},
    u_columns: []u32 = &.{},
    u_values: []f64 = &.{},
    work: []f64 = &.{},
    active: []u32 = &.{},
    hyper_output: []u32 = &.{},
    marked: []bool = &.{},
    u_column_starts: []usize = &.{},
    u_column_rows: []u32 = &.{},
    u_column_values: []f64 = &.{},
    l_row_starts: []usize = &.{},
    l_row_columns: []u32 = &.{},
    l_row_values: []f64 = &.{},

    pivot_threshold: f64 = 0.1,
    zero_tolerance: f64 = 1e-14,

    pub fn init(allocator: std.mem.Allocator) SparseLU {
        return .{
            .allocator = allocator,
            .kernel = sparse_kernel.MutableSparseKernel.init(allocator),
            .symbolic = sparse_symbolic.SymbolicWorkspace.init(allocator),
            .ft = sparse_ft.SparseForrestTomlin.init(allocator),
        };
    }

    pub fn deinit(self: *SparseLU) void {
        self.kernel.deinit();
        self.symbolic.deinit();
        self.ft.deinit();
        self.allocator.free(self.pivot_rows);
        self.allocator.free(self.pivot_columns);
        self.allocator.free(self.row_position);
        self.allocator.free(self.column_position);
        self.allocator.free(self.pivot_values);
        self.allocator.free(self.l_starts);
        self.allocator.free(self.u_starts);
        self.allocator.free(self.l_rows);
        self.allocator.free(self.l_values);
        self.allocator.free(self.u_columns);
        self.allocator.free(self.u_values);
        self.allocator.free(self.work);
        self.allocator.free(self.active);
        self.allocator.free(self.hyper_output);
        self.allocator.free(self.marked);
        self.allocator.free(self.u_column_starts);
        self.allocator.free(self.u_column_rows);
        self.allocator.free(self.u_column_values);
        self.allocator.free(self.l_row_starts);
        self.allocator.free(self.l_row_columns);
        self.allocator.free(self.l_row_values);
        self.* = .{
            .allocator = self.allocator,
            .kernel = sparse_kernel.MutableSparseKernel.init(self.allocator),
            .symbolic = sparse_symbolic.SymbolicWorkspace.init(self.allocator),
            .ft = sparse_ft.SparseForrestTomlin.init(self.allocator),
        };
    }

    pub fn factorize(self: *SparseLU, basis: sparse_basis.SparseBasisView) SparseLuError!void {
        return self.factorizeImpl(basis, true, null, false);
    }

    /// Zero-copy trusted reinversion entry for canonical engine-owned basis
    /// CSC. Workspace and factor capacities are retained across calls.
    pub fn factorizeAssumeValid(self: *SparseLU, basis: sparse_basis.SparseBasisView) SparseLuError!void {
        return self.factorizeImpl(basis, false, null, false);
    }

    /// Replay a previously recorded row/column pivot sequence. This is mainly
    /// a diagnostic control for comparing numerical kernels under identical
    /// ordering; normal reinversion should use threshold Markowitz selection.
    pub fn factorizeWithTraceAssumeValid(self: *SparseLU, basis: sparse_basis.SparseBasisView, trace: PivotTrace) SparseLuError!void {
        if (trace.rows.len != basis.dimension or trace.columns.len != basis.dimension) return error.DimensionMismatch;
        return self.factorizeImpl(basis, false, trace, false);
    }

    /// Reuse the valid prefix of a pivot trace from a neighbouring basis.
    /// Once a recorded pivot is missing or fails the numerical threshold, the
    /// remaining suffix is repaired with the configured ordering backend.
    /// This path retains all workspaces and performs no warm allocations.
    pub fn factorizeWithTraceRepairAssumeValid(self: *SparseLU, basis: sparse_basis.SparseBasisView, trace: PivotTrace) SparseLuError!void {
        if (trace.rows.len != basis.dimension or trace.columns.len != basis.dimension) return error.DimensionMismatch;
        return self.factorizeImpl(basis, false, trace, true);
    }

    fn factorizeImpl(self: *SparseLU, basis: sparse_basis.SparseBasisView, comptime validate: bool, trace: ?PivotTrace, repair_trace: bool) SparseLuError!void {
        const n = basis.dimension;
        try self.ensureDimension(n);
        try self.ensureFactorCapacity(@max(basis.nnz(), n));
        self.dimension = n;
        self.l_nonzeros = 0;
        self.u_nonzeros = 0;
        self.inserted_fill = 0;
        self.peeled_pivots = 0;
        self.kernel_dimension = n;
        self.kernel_nonzeros = basis.nnz();
        self.kernel_maximum_row_count = 0;
        self.kernel_maximum_column_count = 0;
        self.trace_replayed_pivots = 0;
        self.trace_repaired_pivots = 0;
        self.ft_ready = false;
        self.selected_ordering = switch (self.ordering_strategy) {
            .automatic => .dod_markowitz,
            else => |forced| forced,
        };
        self.l_starts[0] = 0;
        self.u_starts[0] = 0;

        var first_kernel_pivot: usize = 0;
        const use_pre_kernel_peeling = trace == null and n >= 192 and basis.nnz() / @max(n, 1) > 4;
        if (use_pre_kernel_peeling) {
            const plan = try self.symbolic.planSingletonPrefix(basis);
            first_kernel_pivot = plan.singleton_pivots;
            self.peeled_pivots = first_kernel_pivot;
            @memset(self.row_position[0..n], @as(u32, @intCast(n)));
            @memset(self.column_position[0..n], @as(u32, @intCast(n)));
            for (plan.pivot_rows[0..first_kernel_pivot], plan.pivot_columns[0..first_kernel_pivot], 0..) |row, column, position| {
                self.row_position[row] = @intCast(position);
                self.column_position[column] = @intCast(position);
            }
            for (0..first_kernel_pivot) |position|
                try self.packSymbolicSingleton(basis, plan, position);
            try self.kernel.loadReducedSymbolicAssumeValid(
                basis,
                plan.active_rows,
                plan.active_columns,
                plan.row_counts,
                plan.column_counts,
            );
            const shape = self.kernel.shape();
            self.kernel_dimension = shape.dimension;
            self.kernel_nonzeros = shape.nonzeros;
            self.kernel_maximum_row_count = shape.maximum_row_count;
            self.kernel_maximum_column_count = shape.maximum_column_count;
            self.selected_ordering = self.selectOrdering(shape);
            self.kernel.setMarkowitzSearchBudget(self.markowitzSearchBudget(shape));
        } else {
            if (validate) try self.kernel.load(basis) else try self.kernel.loadAssumeValid(basis);
            if (trace == null or repair_trace) {
                const shape = self.kernel.shape();
                self.kernel_dimension = shape.dimension;
                self.kernel_nonzeros = shape.nonzeros;
                self.kernel_maximum_row_count = shape.maximum_row_count;
                self.kernel_maximum_column_count = shape.maximum_column_count;
                self.selected_ordering = self.selectOrdering(shape);
                self.kernel.setMarkowitzSearchBudget(self.markowitzSearchBudget(shape));
            }
        }

        var repairing = false;
        for (first_kernel_pivot..n) |pivot_index| {
            var choice: sparse_kernel.PivotChoice = undefined;
            if (trace) |recorded| {
                if (!repairing) {
                    const recorded_choice = if (repair_trace)
                        self.kernel.chooseRecordedPivotThreshold(recorded.rows[pivot_index], recorded.columns[pivot_index], self.pivot_threshold)
                    else
                        self.kernel.chooseRecordedPivot(recorded.rows[pivot_index], recorded.columns[pivot_index]);
                    if (recorded_choice) |accepted| {
                        choice = accepted;
                        self.trace_replayed_pivots += 1;
                    } else if (repair_trace) {
                        repairing = true;
                    } else return error.Singular;
                }
            }
            if (trace == null or repairing) {
                choice = switch (self.selected_ordering) {
                    // The backend split is deliberately established before
                    // the HiGHS-style search lands, so dispatch/API changes
                    // can be verified independently from ordering changes.
                    .dod_markowitz => self.kernel.choosePivot(self.pivot_threshold) orelse return error.Singular,
                    .highs_kernel => self.kernel.choosePivotHighs(self.pivot_threshold) orelse return error.Singular,
                    .automatic => unreachable,
                };
                if (repairing) self.trace_repaired_pivots += 1;
            }
            const pivot = try self.kernel.applyPivot(choice, self.zero_tolerance);
            if (trace == null and self.selected_ordering == .highs_kernel and
                self.kernel_dimension >= self.kernel.adaptive_markowitz_minimum_dimension)
                self.kernel.observeMarkowitzPivot();
            try self.ensureFactorCapacity(@max(self.l_nonzeros + pivot.l_rows.len, self.u_nonzeros + pivot.u_columns.len));
            self.pivot_rows[pivot_index] = pivot.pivot_row;
            self.pivot_columns[pivot_index] = pivot.pivot_column;
            self.row_position[pivot.pivot_row] = @intCast(pivot_index);
            self.column_position[pivot.pivot_column] = @intCast(pivot_index);
            self.pivot_values[pivot_index] = pivot.pivot_value;
            @memcpy(self.l_rows[self.l_nonzeros..][0..pivot.l_rows.len], pivot.l_rows);
            @memcpy(self.l_values[self.l_nonzeros..][0..pivot.l_values.len], pivot.l_values);
            self.l_nonzeros += pivot.l_rows.len;
            self.l_starts[pivot_index + 1] = self.l_nonzeros;
            @memcpy(self.u_columns[self.u_nonzeros..][0..pivot.u_columns.len], pivot.u_columns);
            @memcpy(self.u_values[self.u_nonzeros..][0..pivot.u_values.len], pivot.u_values);
            self.u_nonzeros += pivot.u_columns.len;
            self.u_starts[pivot_index + 1] = self.u_nonzeros;
            self.inserted_fill += pivot.inserted_fill;
        }
        self.hyper_views_ready = false;
    }

    fn packSymbolicSingleton(self: *SparseLU, basis: sparse_basis.SparseBasisView, plan: sparse_symbolic.SymbolicPlanView, position: usize) SparseLuError!void {
        const row = plan.pivot_rows[position];
        const column = plan.pivot_columns[position];
        var pivot_value: f64 = 0.0;
        const column_begin: usize = @intCast(basis.starts[column]);
        const column_end: usize = @intCast(basis.starts[column + 1]);
        for (column_begin..column_end) |entry| if (basis.rows[entry].toUsize() == row) {
            pivot_value = basis.values[entry];
            break;
        };
        if (pivot_value == 0.0 or !std.math.isFinite(pivot_value)) return error.Singular;
        self.pivot_rows[position] = row;
        self.pivot_columns[position] = column;
        self.pivot_values[position] = pivot_value;
        for (column_begin..column_end) |entry| {
            const candidate_row: u32 = @intCast(basis.rows[entry].toUsize());
            if (self.row_position[candidate_row] <= position) continue;
            try self.ensureFactorCapacity(self.l_nonzeros + 1);
            self.l_rows[self.l_nonzeros] = candidate_row;
            self.l_values[self.l_nonzeros] = basis.values[entry] / pivot_value;
            self.l_nonzeros += 1;
        }
        self.l_starts[position + 1] = self.l_nonzeros;
        const row_begin: usize = @intCast(plan.row_starts[row]);
        const row_end: usize = @intCast(plan.row_starts[row + 1]);
        for (plan.row_entries[row_begin..row_end]) |entry| {
            const candidate_column = plan.entry_columns[entry];
            if (self.column_position[candidate_column] <= position) continue;
            try self.ensureFactorCapacity(self.u_nonzeros + 1);
            self.u_columns[self.u_nonzeros] = candidate_column;
            self.u_values[self.u_nonzeros] = basis.values[entry];
            self.u_nonzeros += 1;
        }
        self.u_starts[position + 1] = self.u_nonzeros;
    }

    /// Solve `B x = rhs` in place using the published row/column permutations.
    pub fn solve(self: *SparseLU, rhs: []f64) SparseLuError!void {
        if (self.ft_ready) return self.solveFt(rhs, false);
        const n = self.dimension;
        if (rhs.len != n) return error.DimensionMismatch;
        for (0..n) |k| self.work[k] = rhs[self.pivot_rows[k]];

        // Forward solve L y = P b. L columns scatter one completed component
        // into later pivot rows; no sparse gather or temporary index list.
        for (0..n) |k| {
            const solved = self.work[k];
            for (self.l_rows[self.l_starts[k]..self.l_starts[k + 1]], self.l_values[self.l_starts[k]..self.l_starts[k + 1]]) |row, value| {
                const position = self.row_position[row];
                self.work[position] -= value * solved;
            }
        }
        var k = n;
        while (k > 0) {
            k -= 1;
            var value = self.work[k];
            for (self.u_columns[self.u_starts[k]..self.u_starts[k + 1]], self.u_values[self.u_starts[k]..self.u_starts[k + 1]]) |column, coefficient|
                value -= coefficient * self.work[self.column_position[column]];
            value /= self.pivot_values[k];
            if (!std.math.isFinite(value)) return error.NumericalFailure;
            self.work[k] = value;
        }
        for (0..n) |position| rhs[self.pivot_columns[position]] = self.work[position];
    }

    /// FTRAN that retains the partial vector after L and existing FT row
    /// corrections. The next update consumes this captured `aq` spike.
    pub fn solveForUpdate(self: *SparseLU, rhs: []f64) SparseLuError!void {
        try self.ensureFt();
        return self.solveFt(rhs, true);
    }

    fn solveFt(self: *SparseLU, rhs: []f64, comptime capture: bool) SparseLuError!void {
        const n = self.dimension;
        if (rhs.len != n) return error.DimensionMismatch;
        for (0..n) |k| self.work[k] = rhs[self.pivot_rows[k]];
        for (0..n) |k| {
            const solved = self.work[k];
            for (self.l_rows[self.l_starts[k]..self.l_starts[k + 1]], self.l_values[self.l_starts[k]..self.l_starts[k + 1]]) |row, value|
                self.work[self.row_position[row]] -= value * solved;
        }
        for (0..n) |position| rhs[self.pivot_columns[position]] = self.work[position];
        try self.ft.prepareFtran(rhs, capture);
        try self.ft.solveUpper(rhs);
    }

    /// Solve `B^T x = rhs` in place. U^T uses forward scatter and L^T uses a
    /// reverse gather over the same packed factor streams as FTRAN.
    pub fn solveTranspose(self: *SparseLU, rhs: []f64) SparseLuError!void {
        if (self.ft_ready) return self.solveTransposeFt(rhs);
        const n = self.dimension;
        if (rhs.len != n) return error.DimensionMismatch;
        for (0..n) |k| self.work[k] = rhs[self.pivot_columns[k]];
        for (0..n) |k| {
            const solved = self.work[k] / self.pivot_values[k];
            if (!std.math.isFinite(solved)) return error.NumericalFailure;
            self.work[k] = solved;
            for (self.u_columns[self.u_starts[k]..self.u_starts[k + 1]], self.u_values[self.u_starts[k]..self.u_starts[k + 1]]) |column, coefficient|
                self.work[self.column_position[column]] -= coefficient * solved;
        }
        var k = n;
        while (k > 0) {
            k -= 1;
            var solved = self.work[k];
            for (self.l_rows[self.l_starts[k]..self.l_starts[k + 1]], self.l_values[self.l_starts[k]..self.l_starts[k + 1]]) |row, coefficient|
                solved -= coefficient * self.work[self.row_position[row]];
            self.work[k] = solved;
        }
        for (0..n) |position| rhs[self.pivot_rows[position]] = self.work[position];
    }

    fn solveTransposeFt(self: *SparseLU, rhs: []f64) SparseLuError!void {
        const n = self.dimension;
        if (rhs.len != n) return error.DimensionMismatch;
        try self.ft.solveUpperTranspose(rhs);
        for (0..n) |position| self.work[position] = rhs[self.pivot_columns[position]];
        var k = n;
        while (k > 0) {
            k -= 1;
            var solved = self.work[k];
            for (self.l_rows[self.l_starts[k]..self.l_starts[k + 1]], self.l_values[self.l_starts[k]..self.l_starts[k + 1]]) |row, coefficient|
                solved -= coefficient * self.work[self.row_position[row]];
            self.work[k] = solved;
        }
        for (0..n) |position| rhs[self.pivot_rows[position]] = self.work[position];
    }

    /// Form partial BTRAN `ep`, then mutate U/UR with the captured FTRAN
    /// spike and append its FT row correction.
    pub fn applyForrestTomlinUpdate(self: *SparseLU, leaving_column: u32, direction: []const f64, column_scale: f64) SparseLuError!void {
        if (direction.len != self.dimension or leaving_column >= self.dimension or !std.math.isFinite(column_scale)) return error.DimensionMismatch;
        try self.ensureFt();
        try self.ft.captureEp(leaving_column);
        const alpha = direction[leaving_column] * column_scale;
        try self.ft.update(leaving_column, alpha, self.zero_tolerance);
    }

    pub fn mutableUpperView(self: *SparseLU) SparseLuError!sparse_ft.MutableUpperView {
        try self.ensureFt();
        return self.ft.mutableUpperView();
    }

    /// Solve a hyper-sparse FTRAN from the explicitly nonzero entries of
    /// `rhs`. The factor graph is traversed by scatter reachability; unrelated
    /// L/U rows are never visited. Only returned positions in `rhs` are valid.
    pub fn solveHyperSparse(self: *SparseLU, rhs: []f64, input_indices: []const u32, output_indices: []u32) SparseLuError!usize {
        const n = self.dimension;
        if (rhs.len != n or output_indices.len < n) return error.DimensionMismatch;
        if (self.ft_ready and self.ft.hasUpdates()) {
            // FT-aware hyper-sparse FTRAN: L is unaffected by FT updates, so
            // hyper-sparse L avoids visiting unrelated rows. FT corrections
            // and the modified U solve operate on the (denser) result in
            // original column space.
            self.ensureHyperViews();
            self.clearActive();
            for (input_indices) |row| {
                if (row >= n) return error.DimensionMismatch;
                const position = self.row_position[row];
                self.activate(position, rhs[row]);
            }
            var k: usize = 0;
            while (k < n) : (k += 1) if (self.marked[k]) {
                const solved = self.work[k];
                for (self.l_rows[self.l_starts[k]..self.l_starts[k + 1]], self.l_values[self.l_starts[k]..self.l_starts[k + 1]]) |row, value|
                    self.accumulate(self.row_position[row], -value * solved);
            };
            // Scatter hyper-sparse L result to original column space.
            @memset(rhs[0..n], 0.0);
            for (self.active[0..self.active_count]) |position| {
                const column = self.pivot_columns[position];
                rhs[column] = self.work[position];
            }
            // FT corrections and modified U solve on original-space rhs.
            try self.ft.prepareFtran(rhs, false);
            try self.ft.solveUpper(rhs);
            var output_count: usize = 0;
            for (rhs, 0..) |value, index| if (value != 0.0) {
                output_indices[output_count] = @intCast(index);
                output_count += 1;
            };
            return output_count;
        }
        self.ensureHyperViews();
        self.clearActive();
        for (input_indices) |row| {
            if (row >= n) return error.DimensionMismatch;
            const position = self.row_position[row];
            self.activate(position, rhs[row]);
        }
        // L is stored by column and therefore already is its reachability graph.
        var k: usize = 0;
        while (k < n) : (k += 1) if (self.marked[k]) {
            const solved = self.work[k];
            for (self.l_rows[self.l_starts[k]..self.l_starts[k + 1]], self.l_values[self.l_starts[k]..self.l_starts[k + 1]]) |row, value|
                self.accumulate(self.row_position[row], -value * solved);
        };
        // The companion U column view turns reverse substitution into scatter.
        k = n;
        while (k > 0) {
            k -= 1;
            if (!self.marked[k]) continue;
            const solved = self.work[k] / self.pivot_values[k];
            if (!std.math.isFinite(solved)) return error.NumericalFailure;
            self.work[k] = solved;
            for (self.u_column_rows[self.u_column_starts[k]..self.u_column_starts[k + 1]], self.u_column_values[self.u_column_starts[k]..self.u_column_starts[k + 1]]) |row, value|
                self.accumulate(row, -value * solved);
        }
        var count: usize = 0;
        for (self.active[0..self.active_count]) |position| {
            const column = self.pivot_columns[position];
            rhs[column] = self.work[position];
            output_indices[count] = column;
            count += 1;
        }
        return count;
    }

    /// Hyper-sparse BTRAN using U rows followed by the companion L row graph.
    pub fn solveTransposeHyperSparse(self: *SparseLU, rhs: []f64, input_indices: []const u32, output_indices: []u32) SparseLuError!usize {
        const n = self.dimension;
        if (rhs.len != n or output_indices.len < n) return error.DimensionMismatch;
        if (self.ft_ready and self.ft.hasUpdates()) {
            // FT-aware hyper-sparse BTRAN: U^T and corrections in original
            // space (already sparse via FT linked lists), then hyper-sparse
            // L^T in pivot order to avoid visiting unrelated rows.
            try self.ft.solveUpperTranspose(rhs);
            self.ensureHyperViews();
            self.clearActive();
            for (0..n) |position| {
                const column = self.pivot_columns[position];
                if (rhs[column] != 0.0)
                    self.activate(@intCast(position), rhs[column]);
            }
            var k = n;
            while (k > 0) {
                k -= 1;
                if (!self.marked[k]) continue;
                const solved = self.work[k];
                for (self.l_row_columns[self.l_row_starts[k]..self.l_row_starts[k + 1]], self.l_row_values[self.l_row_starts[k]..self.l_row_starts[k + 1]]) |column, value|
                    self.accumulate(column, -value * solved);
            }
            @memset(rhs[0..n], 0.0);
            for (self.active[0..self.active_count]) |position| {
                const row = self.pivot_rows[position];
                rhs[row] = self.work[position];
            }
            var output_count: usize = 0;
            for (rhs, 0..) |value, index| if (value != 0.0) {
                output_indices[output_count] = @intCast(index);
                output_count += 1;
            };
            return output_count;
        }
        self.ensureHyperViews();
        self.clearActive();
        for (input_indices) |column| {
            if (column >= n) return error.DimensionMismatch;
            self.activate(self.column_position[column], rhs[column]);
        }
        for (0..n) |k| if (self.marked[k]) {
            const solved = self.work[k] / self.pivot_values[k];
            if (!std.math.isFinite(solved)) return error.NumericalFailure;
            self.work[k] = solved;
            for (self.u_columns[self.u_starts[k]..self.u_starts[k + 1]], self.u_values[self.u_starts[k]..self.u_starts[k + 1]]) |column, value|
                self.accumulate(self.column_position[column], -value * solved);
        };
        var k = n;
        while (k > 0) {
            k -= 1;
            if (!self.marked[k]) continue;
            const solved = self.work[k];
            for (self.l_row_columns[self.l_row_starts[k]..self.l_row_starts[k + 1]], self.l_row_values[self.l_row_starts[k]..self.l_row_starts[k + 1]]) |column, value|
                self.accumulate(column, -value * solved);
        }
        var count: usize = 0;
        for (self.active[0..self.active_count]) |position| {
            const row = self.pivot_rows[position];
            rhs[row] = self.work[position];
            output_indices[count] = row;
            count += 1;
        }
        return count;
    }

    /// Select hyper-sparse traversal below 12.5% RHS density, otherwise use
    /// the sequential packed kernel. Hyper-sparse now handles FT updates via
    /// an FT-aware L-only traversal, so the FT-update guard is removed: the
    /// L factor is unaffected by Forrest-Tomlin modifications and visiting
    /// only the active rows still saves work on sparse RHS. The returned RHS
    /// is always dense.
    pub fn solveAdaptive(self: *SparseLU, rhs: []f64, input_indices: []const u32, transpose: bool) SparseLuError!void {
        if (input_indices.len * 8 >= self.dimension) return if (transpose) self.solveTranspose(rhs) else self.solve(rhs);
        const count = if (transpose)
            try self.solveTransposeHyperSparse(rhs, input_indices, self.hyper_output)
        else
            try self.solveHyperSparse(rhs, input_indices, self.hyper_output);
        // Hyper methods leave only their output positions defined. Preserve the
        // sparse result while materializing the conventional dense API.
        @memset(self.work[0..self.dimension], 0.0);
        for (self.hyper_output[0..count]) |index| self.work[index] = rhs[index];
        @memcpy(rhs, self.work[0..self.dimension]);
    }

    pub fn factorNonzeros(self: *const SparseLU) usize {
        return self.dimension + self.l_nonzeros + self.u_nonzeros;
    }

    /// Cheap stability signal matching DenseLU's policy boundary. This is a
    /// pivot spread indicator, not a formal condition-number estimate.
    pub fn pivotConditionEstimate(self: *const SparseLU) f64 {
        if (self.dimension == 0) return 1.0;
        var minimum = std.math.inf(f64);
        var maximum: f64 = 0.0;
        for (self.pivot_values[0..self.dimension]) |pivot| {
            const magnitude = @abs(pivot);
            minimum = @min(minimum, magnitude);
            maximum = @max(maximum, magnitude);
        }
        if (minimum <= 0.0 or !std.math.isFinite(minimum) or !std.math.isFinite(maximum)) return std.math.inf(f64);
        return maximum / minimum;
    }

    /// Retained requested bytes for the complete numerical kernel, factors,
    /// permutations and solve workspace. Allocator bookkeeping is excluded.
    pub fn requestedBytes(self: *const SparseLU) usize {
        var total = self.kernel.requestedBytes() + self.symbolic.retainedBytes() + self.ft.retainedBytes();
        total += std.mem.sliceAsBytes(self.pivot_rows).len;
        total += std.mem.sliceAsBytes(self.pivot_columns).len;
        total += std.mem.sliceAsBytes(self.row_position).len;
        total += std.mem.sliceAsBytes(self.column_position).len;
        total += std.mem.sliceAsBytes(self.pivot_values).len;
        total += std.mem.sliceAsBytes(self.l_starts).len;
        total += std.mem.sliceAsBytes(self.u_starts).len;
        total += std.mem.sliceAsBytes(self.l_rows).len;
        total += std.mem.sliceAsBytes(self.l_values).len;
        total += std.mem.sliceAsBytes(self.u_columns).len;
        total += std.mem.sliceAsBytes(self.u_values).len;
        total += std.mem.sliceAsBytes(self.work).len;
        total += std.mem.sliceAsBytes(self.active).len;
        total += std.mem.sliceAsBytes(self.hyper_output).len;
        total += std.mem.sliceAsBytes(self.marked).len;
        total += std.mem.sliceAsBytes(self.u_column_starts).len;
        total += std.mem.sliceAsBytes(self.u_column_rows).len;
        total += std.mem.sliceAsBytes(self.u_column_values).len;
        total += std.mem.sliceAsBytes(self.l_row_starts).len;
        total += std.mem.sliceAsBytes(self.l_row_columns).len;
        total += std.mem.sliceAsBytes(self.l_row_values).len;
        return total;
    }

    fn selectOrdering(self: *const SparseLU, shape: sparse_kernel.KernelShape) OrderingStrategy {
        return switch (self.ordering_strategy) {
            .automatic => if (self.peeled_pivots >= self.dimension - self.peeled_pivots and
                shape.nonzeros / @max(shape.dimension, 1) > 4 and
                shape.maximum_row_count <= 64 and shape.maximum_column_count <= 64)
                .highs_kernel
            else
                .dod_markowitz,
            else => |forced| forced,
        };
    }

    /// Compact high-fill kernels left after a long singleton prefix benefit
    /// from taking the first threshold-safe low-degree frontier candidate.
    /// This DOD cost gate avoids paying for a wider HiGHS-style search whose
    /// additional candidates empirically increase both latency and fill.
    fn markowitzSearchBudget(self: *const SparseLU, shape: sparse_kernel.KernelShape) usize {
        if (self.selected_ordering == .highs_kernel and
            self.peeled_pivots >= self.dimension - self.peeled_pivots and
            shape.dimension <= 128 and shape.nonzeros / @max(shape.dimension, 1) > 4 and
            shape.maximum_row_count <= 64 and shape.maximum_column_count <= 64)
            return 1;
        return 8;
    }

    fn ensureDimension(self: *SparseLU, required: usize) SparseLuError!void {
        if (required <= self.dimension_capacity) return;
        const capacity = grow(self.dimension_capacity, required) catch return error.CapacityOverflow;
        try self.resizeRetained(&self.pivot_rows, capacity);
        try self.resizeRetained(&self.pivot_columns, capacity);
        try self.resizeRetained(&self.row_position, capacity);
        try self.resizeRetained(&self.column_position, capacity);
        try self.resizeRetained(&self.pivot_values, capacity);
        try self.resizeRetained(&self.work, capacity);
        try self.resizeRetained(&self.active, capacity);
        try self.resizeRetained(&self.hyper_output, capacity);
        try self.resizeRetained(&self.marked, capacity);
        self.l_starts = self.allocator.realloc(self.l_starts, capacity + 1) catch return error.OutOfMemory;
        self.u_starts = self.allocator.realloc(self.u_starts, capacity + 1) catch return error.OutOfMemory;
        self.u_column_starts = self.allocator.realloc(self.u_column_starts, capacity + 1) catch return error.OutOfMemory;
        self.l_row_starts = self.allocator.realloc(self.l_row_starts, capacity + 1) catch return error.OutOfMemory;
        self.dimension_capacity = capacity;
    }

    fn ensureFactorCapacity(self: *SparseLU, required: usize) SparseLuError!void {
        if (required <= self.factor_capacity) return;
        const capacity = grow(self.factor_capacity, required) catch return error.CapacityOverflow;
        try self.resizeRetained(&self.l_rows, capacity);
        try self.resizeRetained(&self.l_values, capacity);
        try self.resizeRetained(&self.u_columns, capacity);
        try self.resizeRetained(&self.u_values, capacity);
        try self.resizeRetained(&self.u_column_rows, capacity);
        try self.resizeRetained(&self.u_column_values, capacity);
        try self.resizeRetained(&self.l_row_columns, capacity);
        try self.resizeRetained(&self.l_row_values, capacity);
        self.factor_capacity = capacity;
    }

    fn resizeRetained(self: *SparseLU, slice: anytype, capacity: usize) SparseLuError!void {
        slice.* = self.allocator.realloc(slice.*, capacity) catch return error.OutOfMemory;
    }

    fn buildHyperViews(self: *SparseLU) void {
        const n = self.dimension;
        @memset(self.u_column_starts[0 .. n + 1], 0);
        @memset(self.l_row_starts[0 .. n + 1], 0);
        for (self.u_columns[0..self.u_nonzeros]) |column| self.u_column_starts[self.column_position[column] + 1] += 1;
        for (self.l_rows[0..self.l_nonzeros]) |row| self.l_row_starts[self.row_position[row] + 1] += 1;
        for (0..n) |i| {
            self.u_column_starts[i + 1] += self.u_column_starts[i];
            self.l_row_starts[i + 1] += self.l_row_starts[i];
        }
        for (0..n) |i| self.active[i] = @intCast(self.u_column_starts[i]);
        for (0..n) |row| for (self.u_starts[row]..self.u_starts[row + 1]) |entry| {
            const column = self.column_position[self.u_columns[entry]];
            const output = self.active[column];
            self.u_column_rows[output] = @intCast(row);
            self.u_column_values[output] = self.u_values[entry];
            self.active[column] += 1;
        };
        for (0..n) |i| self.active[i] = @intCast(self.l_row_starts[i]);
        for (0..n) |column| for (self.l_starts[column]..self.l_starts[column + 1]) |entry| {
            const row = self.row_position[self.l_rows[entry]];
            const output = self.active[row];
            self.l_row_columns[output] = @intCast(column);
            self.l_row_values[output] = self.l_values[entry];
            self.active[row] += 1;
        };
        @memset(self.marked[0..n], false);
        self.active_count = 0;
        self.hyper_views_ready = true;
    }

    inline fn ensureHyperViews(self: *SparseLU) void {
        if (!self.hyper_views_ready) self.buildHyperViews();
    }

    fn ensureFt(self: *SparseLU) SparseLuError!void {
        if (self.ft_ready) return;
        self.ensureHyperViews();
        try self.ft.reset(
            self.pivot_columns[0..self.dimension],
            self.pivot_values[0..self.dimension],
            self.u_column_starts[0 .. self.dimension + 1],
            self.u_column_rows[0..self.u_nonzeros],
            self.u_column_values[0..self.u_nonzeros],
        );
        self.ft_ready = true;
    }

    fn clearActive(self: *SparseLU) void {
        for (self.active[0..self.active_count]) |position| self.marked[position] = false;
        self.active_count = 0;
    }
    fn activate(self: *SparseLU, position: u32, value: f64) void {
        if (!self.marked[position]) {
            self.marked[position] = true;
            self.active[self.active_count] = position;
            self.active_count += 1;
            self.work[position] = value;
        } else self.work[position] += value;
    }
    fn accumulate(self: *SparseLU, position: u32, value: f64) void {
        self.activate(position, value);
    }
};

fn grow(current: usize, required: usize) error{Overflow}!usize {
    var capacity = @max(current, 8);
    while (capacity < required) capacity = std.math.add(usize, capacity, capacity / 2 + 8) catch return error.Overflow;
    return capacity;
}

test "packed sparse LU solves FTRAN and BTRAN with fill" {
    const foundation = @import("foundation");
    const RowId = foundation.RowId;
    const Offset = foundation.HUInt;
    const starts = [_]Offset{ 0, 2, 4, 6 };
    const rows = [_]RowId{
        RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(1),
        RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(2),
        RowId.fromUsizeAssumeValid(1), RowId.fromUsizeAssumeValid(2),
    };
    const values = [_]f64{ 2, 1, 1, 2, 3, 4 };
    const basis = sparse_basis.SparseBasisView{ .dimension = 3, .starts = &starts, .rows = &rows, .values = &values };
    var lu = SparseLU.init(std.testing.allocator);
    defer lu.deinit();
    try lu.factorize(basis);
    try std.testing.expect(lu.inserted_fill >= 1);
    var trace_rows: [3]u32 = undefined;
    var trace_columns: [3]u32 = undefined;
    @memcpy(&trace_rows, lu.pivot_rows[0..3]);
    @memcpy(&trace_columns, lu.pivot_columns[0..3]);
    const factor_nonzeros = lu.factorNonzeros();
    try lu.factorizeWithTraceAssumeValid(basis, .{ .rows = &trace_rows, .columns = &trace_columns });
    try std.testing.expectEqual(factor_nonzeros, lu.factorNonzeros());

    var rhs = [_]f64{ 4, 7, 10 };
    const original_rhs = rhs;
    try lu.solve(&rhs);
    for (0..3) |row| {
        var product: f64 = 0.0;
        for (0..3) |column| {
            for (@as(usize, starts[column])..@as(usize, starts[column + 1])) |entry| {
                if (rows[entry].toUsize() == row) product += values[entry] * rhs[column];
            }
        }
        try std.testing.expectApproxEqAbs(original_rhs[row], product, 1e-11);
    }

    var transpose_rhs = [_]f64{ 3, -2, 5 };
    const original_transpose_rhs = transpose_rhs;
    try lu.solveTranspose(&transpose_rhs);
    for (0..3) |column| {
        var product: f64 = 0.0;
        for (@as(usize, starts[column])..@as(usize, starts[column + 1])) |entry|
            product += values[entry] * transpose_rhs[rows[entry].toUsize()];
        try std.testing.expectApproxEqAbs(original_transpose_rhs[column], product, 1e-11);
    }
}

test "validated trace replays prefix and repairs invalid suffix" {
    const foundation = @import("foundation");
    const basis = sparse_basis.SparseBasisView{
        .dimension = 4,
        .starts = &[_]foundation.HUInt{ 0, 2, 4, 6, 8 },
        .rows = &[_]foundation.RowId{
            foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1),
            foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(2),
            foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(3),
            foundation.RowId.fromUsizeAssumeValid(2), foundation.RowId.fromUsizeAssumeValid(3),
        },
        .values = &[_]f64{ 4, 1, 1, 3, 2, 5, 1, 4 },
    };
    var lu = SparseLU.init(std.testing.allocator);
    defer lu.deinit();
    try lu.factorize(basis);
    var trace_rows: [4]u32 = undefined;
    var trace_columns: [4]u32 = undefined;
    @memcpy(&trace_rows, lu.pivot_rows[0..4]);
    @memcpy(&trace_columns, lu.pivot_columns[0..4]);

    // Reusing the first pivot's now-retired row makes the second recorded
    // pivot invalid. The valid prefix remains useful and only the suffix is
    // reordered.
    trace_rows[1] = trace_rows[0];
    trace_columns[1] = trace_columns[0];
    try lu.factorizeWithTraceRepairAssumeValid(basis, .{ .rows = &trace_rows, .columns = &trace_columns });
    try std.testing.expectEqual(@as(usize, 1), lu.trace_replayed_pivots);
    try std.testing.expectEqual(@as(usize, 3), lu.trace_repaired_pivots);

    var rhs = [_]f64{ 3, -2, 5, 7 };
    const original = rhs;
    try lu.solve(&rhs);
    for (0..4) |row| {
        var product: f64 = 0.0;
        for (0..4) |column| for (@as(usize, basis.starts[column])..@as(usize, basis.starts[column + 1])) |entry|
            if (basis.rows[entry].toUsize() == row) {
                product += basis.values[entry] * rhs[column];
            };
        try std.testing.expectApproxEqAbs(original[row], product, 1e-11);
    }
}

test "validated trace rejects a pivot below the current column threshold" {
    const foundation = @import("foundation");
    const basis = sparse_basis.SparseBasisView{
        .dimension = 2,
        .starts = &[_]foundation.HUInt{ 0, 2, 4 },
        .rows = &[_]foundation.RowId{
            foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1),
            foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1),
        },
        .values = &[_]f64{ 1e-3, 1, 2, 3 },
    };
    var kernel = sparse_kernel.MutableSparseKernel.init(std.testing.allocator);
    defer kernel.deinit();
    try kernel.load(basis);
    try std.testing.expect(kernel.chooseRecordedPivot(0, 0) != null);
    try std.testing.expect(kernel.chooseRecordedPivotThreshold(0, 0, 0.1) == null);
    try std.testing.expect(kernel.chooseRecordedPivotThreshold(1, 0, 0.1) != null);
}

test "sparse LU requested bytes account for retained buffers" {
    const foundation = @import("foundation");
    const basis = sparse_basis.SparseBasisView{
        .dimension = 2,
        .starts = &[_]foundation.HUInt{ 0, 1, 2 },
        .rows = &[_]foundation.RowId{
            foundation.RowId.fromUsizeAssumeValid(0),
            foundation.RowId.fromUsizeAssumeValid(1),
        },
        .values = &[_]f64{ 1, 1 },
    };
    var lu = SparseLU.init(std.testing.allocator);
    defer lu.deinit();
    try lu.factorize(basis);
    try std.testing.expect(lu.requestedBytes() >= lu.factorNonzeros() * @sizeOf(f64));
}

test "factor graph hyper sparse FTRAN BTRAN and adaptive dispatch match dense" {
    const foundation = @import("foundation");
    const basis = sparse_basis.SparseBasisView{
        .dimension = 4,
        .starts = &[_]foundation.HUInt{ 0, 2, 4, 6, 8 },
        .rows = &[_]foundation.RowId{
            foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1),
            foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(2),
            foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(3),
            foundation.RowId.fromUsizeAssumeValid(2), foundation.RowId.fromUsizeAssumeValid(3),
        },
        .values = &[_]f64{ 4, 1, 1, 3, 2, 5, 1, 4 },
    };
    var lu = SparseLU.init(std.testing.allocator);
    defer lu.deinit();
    try lu.factorize(basis);
    var dense = [_]f64{ 0, 2, 0, 0 };
    var hyper = dense;
    var output: [4]u32 = undefined;
    try lu.solve(&dense);
    const count = try lu.solveHyperSparse(&hyper, &[_]u32{1}, &output);
    var rebuilt = [_]f64{0} ** 4;
    for (output[0..count]) |index| rebuilt[index] = hyper[index];
    try std.testing.expectEqualSlices(f64, &dense, &rebuilt);

    dense = .{ 0, 0, -3, 0 };
    hyper = dense;
    try lu.solveTranspose(&dense);
    const transpose_count = try lu.solveTransposeHyperSparse(&hyper, &[_]u32{2}, &output);
    rebuilt = .{0} ** 4;
    for (output[0..transpose_count]) |index| rebuilt[index] = hyper[index];
    try std.testing.expectEqualSlices(f64, &dense, &rebuilt);

    hyper = .{ 7, 0, 0, 0 };
    dense = hyper;
    try lu.solve(&dense);
    try lu.solveAdaptive(&hyper, &[_]u32{0}, false);
    try std.testing.expectEqualSlices(f64, &dense, &hyper);
}

test "sparse ordering backends can be forced without changing the public factors" {
    const foundation = @import("foundation");
    const basis = sparse_basis.SparseBasisView{
        .dimension = 3,
        .starts = &[_]foundation.HUInt{ 0, 2, 4, 6 },
        .rows = &[_]foundation.RowId{
            foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1),
            foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(2),
            foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(2),
        },
        .values = &[_]f64{ 2, 1, 1, 2, 3, 4 },
    };
    var lu = SparseLU.init(std.testing.allocator);
    defer lu.deinit();
    inline for (.{ OrderingStrategy.dod_markowitz, OrderingStrategy.highs_kernel }) |strategy| {
        lu.ordering_strategy = strategy;
        try lu.factorize(basis);
        try std.testing.expectEqual(strategy, lu.selected_ordering);
        var rhs = [_]f64{ 4, 7, 10 };
        try lu.solve(&rhs);
        try std.testing.expectApproxEqAbs(@as(f64, 1.375), rhs[0], 1e-12);
    }
}

test "compact peeled high fill kernel gates Markowitz search effort" {
    var lu = SparseLU.init(std.testing.allocator);
    defer lu.deinit();
    lu.dimension = 220;
    lu.peeled_pivots = 130;
    lu.selected_ordering = .highs_kernel;

    const compact = sparse_kernel.KernelShape{
        .dimension = 90,
        .nonzeros = 636,
        .maximum_row_count = 26,
        .maximum_column_count = 23,
    };
    try std.testing.expectEqual(@as(usize, 1), lu.markowitzSearchBudget(compact));

    lu.selected_ordering = .dod_markowitz;
    try std.testing.expectEqual(@as(usize, 8), lu.markowitzSearchBudget(compact));
    lu.selected_ordering = .highs_kernel;
    lu.peeled_pivots = 100;
    try std.testing.expectEqual(@as(usize, 8), lu.markowitzSearchBudget(compact));
}

test "sparse LU matches original products across deterministic sparse bases" {
    const foundation = @import("foundation");
    const RowId = foundation.RowId;
    const Offset = foundation.HUInt;
    var starts: [9]Offset = undefined;
    var rows: [64]RowId = undefined;
    var values: [64]f64 = undefined;

    for (2..9) |n| {
        var nnz: usize = 0;
        for (0..n) |column| {
            starts[column] = @intCast(nnz);
            for (0..n) |row| {
                if (row != column and (row * 7 + column * 11) % 4 != 0) continue;
                rows[nnz] = RowId.fromUsizeAssumeValid(row);
                values[nnz] = if (row == column)
                    8.0 + @as(f64, @floatFromInt(row))
                else
                    @as(f64, @floatFromInt(@as(i32, @intCast((row * 5 + column * 3) % 7)) - 3)) * 0.125 + 0.25;
                nnz += 1;
            }
        }
        starts[n] = @intCast(nnz);
        const basis = sparse_basis.SparseBasisView{ .dimension = n, .starts = starts[0 .. n + 1], .rows = rows[0..nnz], .values = values[0..nnz] };
        var lu = SparseLU.init(std.testing.allocator);
        defer lu.deinit();
        try lu.factorize(basis);

        var rhs: [8]f64 = undefined;
        var original: [8]f64 = undefined;
        for (0..n) |index| rhs[index] = @as(f64, @floatFromInt(index + 1)) * 0.75 - 1.0;
        @memcpy(original[0..n], rhs[0..n]);
        try lu.solve(rhs[0..n]);
        for (0..n) |row| {
            var product: f64 = 0.0;
            for (0..n) |column| {
                for (@as(usize, starts[column])..@as(usize, starts[column + 1])) |entry| {
                    if (rows[entry].toUsize() == row) product += values[entry] * rhs[column];
                }
            }
            try std.testing.expectApproxEqAbs(original[row], product, 1e-10);
        }

        for (0..n) |index| rhs[index] = @as(f64, @floatFromInt(index + 2)) * -0.5 + 0.25;
        @memcpy(original[0..n], rhs[0..n]);
        try lu.solveTranspose(rhs[0..n]);
        for (0..n) |column| {
            var product: f64 = 0.0;
            for (@as(usize, starts[column])..@as(usize, starts[column + 1])) |entry|
                product += values[entry] * rhs[rows[entry].toUsize()];
            try std.testing.expectApproxEqAbs(original[column], product, 1e-10);
        }
    }
}

test "warm reinversion and solves perform no allocator calls" {
    const foundation = @import("foundation");
    const RowId = foundation.RowId;
    const Offset = foundation.HUInt;
    const basis = sparse_basis.SparseBasisView{
        .dimension = 3,
        .starts = &[_]Offset{ 0, 2, 4, 6 },
        .rows = &[_]RowId{
            RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(1),
            RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(2),
            RowId.fromUsizeAssumeValid(1), RowId.fromUsizeAssumeValid(2),
        },
        .values = &[_]f64{ 2, 1, 1, 2, 3, 4 },
    };
    var counted = std.testing.FailingAllocator.init(std.testing.allocator, .{});
    var lu = SparseLU.init(counted.allocator());
    defer lu.deinit();
    try lu.factorize(basis);
    const allocations = counted.allocations;
    const resizes = counted.resize_index;
    try lu.factorizeAssumeValid(basis);
    var rhs = [_]f64{ 1, 2, 3 };
    try lu.solve(&rhs);
    try lu.solveTranspose(&rhs);
    try std.testing.expectEqual(allocations, counted.allocations);
    try std.testing.expectEqual(resizes, counted.resize_index);
}

test "Forrest Tomlin repeated basis replacements match dense reinversion" {
    const foundation = @import("foundation");
    const n = 4;
    const basis = sparse_basis.SparseBasisView{
        .dimension = n,
        .starts = &[_]foundation.HUInt{ 0, 3, 6, 9, 10 },
        .rows = &[_]foundation.RowId{
            foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(3),
            foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(2),
            foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(2), foundation.RowId.fromUsizeAssumeValid(3),
            foundation.RowId.fromUsizeAssumeValid(3),
        },
        .values = &[_]f64{ 6, 1, 1, 1, 5, 1, 1, 4, 1, 5 },
    };
    var sparse = SparseLU.init(std.testing.allocator);
    defer sparse.deinit();
    try sparse.factorize(basis);

    var dense_matrix = [_]f64{
        6, 1, 0, 0,
        1, 5, 1, 0,
        0, 1, 4, 0,
        1, 0, 1, 5,
    };
    const replacements = [_]struct { column: usize, values: [n]f64 }{
        .{ .column = 1, .values = .{ 2, 3, 1, 1 } },
        .{ .column = 3, .values = .{ 1, 2, 1, 4 } },
        .{ .column = 0, .values = .{ 5, 1, 2, 1 } },
        .{ .column = 2, .values = .{ 1, 2, 5, 2 } },
    };
    var oracle = @import("dense_lu.zig").DenseLU.init(std.testing.allocator);
    defer oracle.deinit();

    for (replacements, 0..) |replacement, replacement_index| {
        var aq = replacement.values;
        try sparse.solveForUpdate(&aq);
        const column_scale: f64 = if (replacement_index == 2) -1.0 else 1.0;
        if (column_scale < 0.0) {
            for (&aq) |*value| value.* = -value.*;
        }
        try sparse.applyForrestTomlinUpdate(@intCast(replacement.column), &aq, column_scale);
        for (replacement.values, 0..) |value, row| dense_matrix[row * n + replacement.column] = value;
        try oracle.factorize(n, &dense_matrix);

        var sparse_rhs: [n]f64 = undefined;
        for (&sparse_rhs, 0..) |*value, index| value.* = @as(f64, @floatFromInt(2 + replacement_index * 3 + index));
        var dense_rhs = sparse_rhs;
        try sparse.solve(&sparse_rhs);
        try oracle.solve(&dense_rhs);
        for (sparse_rhs, dense_rhs) |actual, expected| try std.testing.expectApproxEqAbs(expected, actual, 2e-11);

        for (&sparse_rhs, 0..) |*value, index| value.* = @as(f64, @floatFromInt(1 + replacement_index + index * 2));
        dense_rhs = sparse_rhs;
        try sparse.solveTranspose(&sparse_rhs);
        try oracle.solveTranspose(&dense_rhs);
        for (sparse_rhs, dense_rhs) |actual, expected| try std.testing.expectApproxEqAbs(expected, actual, 2e-11);
    }
    const upper = try sparse.mutableUpperView();
    try std.testing.expectEqual(n + replacements.len, upper.pivot_ids.len);
    try std.testing.expectEqual(replacements.len, sparse.ft.update_count);
}

test "Forrest Tomlin update ep excludes historical row corrections" {
    const foundation = @import("foundation");
    const n = 3;
    const basis = sparse_basis.SparseBasisView{
        .dimension = n,
        .starts = &[_]foundation.HUInt{ 0, 1, 3, 5 },
        .rows = &[_]foundation.RowId{
            foundation.RowId.fromUsizeAssumeValid(0),
            foundation.RowId.fromUsizeAssumeValid(0),
            foundation.RowId.fromUsizeAssumeValid(1),
            foundation.RowId.fromUsizeAssumeValid(1),
            foundation.RowId.fromUsizeAssumeValid(2),
        },
        .values = &[_]f64{ 2, 1, 3, 1, 4 },
    };
    var sparse = SparseLU.init(std.testing.allocator);
    defer sparse.deinit();
    try sparse.factorize(basis);
    var dense_matrix = [_]f64{ 2, 1, 0, 0, 3, 1, 0, 0, 4 };
    const replacements = [_]struct { column: usize, values: [n]f64 }{
        .{ .column = 0, .values = .{ 1, 2, 1 } },
        .{ .column = 1, .values = .{ 2, 1, 3 } },
    };
    for (replacements) |replacement| {
        var direction = replacement.values;
        try sparse.solveForUpdate(&direction);
        try sparse.applyForrestTomlinUpdate(@intCast(replacement.column), &direction, 1.0);
        for (replacement.values, 0..) |value, row| dense_matrix[row * n + replacement.column] = value;
    }
    var oracle = @import("dense_lu.zig").DenseLU.init(std.testing.allocator);
    defer oracle.deinit();
    try oracle.factorize(n, &dense_matrix);
    var actual = [_]f64{ 3, -2, 5 };
    var expected = actual;
    try sparse.solve(&actual);
    try oracle.solve(&expected);
    for (actual, expected) |value, reference| try std.testing.expectApproxEqAbs(reference, value, 1e-12);
    actual = .{ -1, 4, 2 };
    expected = actual;
    try sparse.solveTranspose(&actual);
    try oracle.solveTranspose(&expected);
    for (actual, expected) |value, reference| try std.testing.expectApproxEqAbs(reference, value, 1e-12);
}

test "retained Forrest Tomlin workspace performs warm updates without allocation" {
    const foundation = @import("foundation");
    const basis = sparse_basis.SparseBasisView{
        .dimension = 4,
        .starts = &[_]foundation.HUInt{ 0, 1, 2, 3, 4 },
        .rows = &[_]foundation.RowId{
            foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1),
            foundation.RowId.fromUsizeAssumeValid(2), foundation.RowId.fromUsizeAssumeValid(3),
        },
        .values = &[_]f64{ 1, 1, 1, 1 },
    };
    var counted = std.testing.FailingAllocator.init(std.testing.allocator, .{});
    var sparse = SparseLU.init(counted.allocator());
    defer sparse.deinit();
    try sparse.factorize(basis);
    var aq = [_]f64{ 2, 1, 1, 1 };
    try sparse.solveForUpdate(&aq);
    const allocations = counted.allocations;
    const resizes = counted.resize_index;
    try sparse.applyForrestTomlinUpdate(0, &aq, 1);
    aq = .{ 1, 3, 1, 1 };
    try sparse.solveForUpdate(&aq);
    try sparse.applyForrestTomlinUpdate(1, &aq, 1);
    aq = .{ 1, 1, 4, 1 };
    try sparse.solveForUpdate(&aq);
    try sparse.applyForrestTomlinUpdate(2, &aq, 1);
    var rhs = [_]f64{ 1, 2, 3, 4 };
    try sparse.solve(&rhs);
    try sparse.solveTranspose(&rhs);
    try std.testing.expectEqual(allocations, counted.allocations);
    try std.testing.expectEqual(resizes, counted.resize_index);
}

test "Forrest Tomlin deterministic replacement corpus matches dense oracle" {
    const foundation = @import("foundation");
    var starts: [9]foundation.HUInt = undefined;
    var rows: [64]foundation.RowId = undefined;
    var values: [64]f64 = undefined;
    var dense_matrix: [64]f64 = undefined;

    for (2..9) |n| {
        var nnz: usize = 0;
        @memset(dense_matrix[0 .. n * n], 0.0);
        for (0..n) |column| {
            starts[column] = @intCast(nnz);
            for (0..n) |row| {
                const value = if (row == column)
                    9.0 + @as(f64, @floatFromInt(column))
                else
                    @as(f64, @floatFromInt((row * 7 + column * 5) % 5 + 1)) * 0.04;
                rows[nnz] = foundation.RowId.fromUsizeAssumeValid(row);
                values[nnz] = value;
                dense_matrix[row * n + column] = value;
                nnz += 1;
            }
        }
        starts[n] = @intCast(nnz);
        const basis = sparse_basis.SparseBasisView{ .dimension = n, .starts = starts[0 .. n + 1], .rows = rows[0..nnz], .values = values[0..nnz] };
        var sparse = SparseLU.init(std.testing.allocator);
        defer sparse.deinit();
        try sparse.factorize(basis);
        var oracle = @import("dense_lu.zig").DenseLU.init(std.testing.allocator);
        defer oracle.deinit();

        for (0..2 * n) |replacement_index| {
            const column = (replacement_index * 3 + 1) % n;
            var entering: [8]f64 = undefined;
            for (0..n) |row| {
                entering[row] = if (row == column)
                    8.0 + @as(f64, @floatFromInt(replacement_index)) * 0.125
                else
                    @as(f64, @floatFromInt((row * 11 + replacement_index * 3) % 7 + 1)) * 0.03;
            }
            var aq: [8]f64 = undefined;
            @memcpy(aq[0..n], entering[0..n]);
            try sparse.solveForUpdate(aq[0..n]);
            try sparse.applyForrestTomlinUpdate(@intCast(column), aq[0..n], 1);
            for (entering[0..n], 0..) |value, row| dense_matrix[row * n + column] = value;
            try oracle.factorize(n, dense_matrix[0 .. n * n]);

            var actual: [8]f64 = undefined;
            for (actual[0..n], 0..) |*value, index| value.* = @as(f64, @floatFromInt(index + replacement_index + 1)) * 0.375;
            var expected = actual;
            try sparse.solve(actual[0..n]);
            try oracle.solve(expected[0..n]);
            for (actual[0..n], expected[0..n]) |a, e| try std.testing.expectApproxEqAbs(e, a, 1e-10);

            for (actual[0..n], 0..) |*value, index| value.* = @as(f64, @floatFromInt(index * 2 + replacement_index + 2)) * -0.25;
            expected = actual;
            try sparse.solveTranspose(actual[0..n]);
            try oracle.solveTranspose(expected[0..n]);
            for (actual[0..n], expected[0..n]) |a, e| try std.testing.expectApproxEqAbs(e, a, 1e-10);
        }
    }
}
