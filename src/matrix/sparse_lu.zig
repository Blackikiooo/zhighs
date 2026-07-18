//! Packed sparse LU built from the mutable elimination kernel.
//!
//! The factorization records `P * B * Q = L * U`. L is stored by pivot
//! column without its unit diagonal; U stores its diagonal separately and
//! off-diagonal entries by pivot row. Original row/column IDs are retained in
//! the packed streams and mapped to pivot positions by compact lookup arrays.

const std = @import("std");
const sparse_basis = @import("sparse_basis.zig");
const sparse_kernel = @import("sparse_kernel.zig");

pub const SparseLuError = sparse_kernel.KernelError || error{ DimensionMismatch, NumericalFailure };

pub const PivotTrace = struct { rows: []const u32, columns: []const u32 };
pub const OrderingStrategy = enum { automatic, dod_markowitz, highs_kernel };

pub const SparseLU = struct {
    allocator: std.mem.Allocator,
    kernel: sparse_kernel.MutableSparseKernel,
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
        return .{ .allocator = allocator, .kernel = sparse_kernel.MutableSparseKernel.init(allocator) };
    }

    pub fn deinit(self: *SparseLU) void {
        self.kernel.deinit();
        inline for (.{
            "pivot_rows", "pivot_columns", "row_position", "column_position", "pivot_values",
            "l_starts",   "u_starts",      "l_rows",       "l_values",        "u_columns",
            "u_values",   "work",           "active",       "hyper_output", "marked",
            "u_column_starts", "u_column_rows", "u_column_values",
            "l_row_starts", "l_row_columns", "l_row_values",
        }) |field_name| self.allocator.free(@field(self, field_name));
        self.* = .{ .allocator = self.allocator, .kernel = sparse_kernel.MutableSparseKernel.init(self.allocator) };
    }

    pub fn factorize(self: *SparseLU, basis: sparse_basis.SparseBasisView) SparseLuError!void {
        return self.factorizeImpl(basis, true, null);
    }

    /// Zero-copy trusted reinversion entry for canonical engine-owned basis
    /// CSC. Workspace and factor capacities are retained across calls.
    pub fn factorizeAssumeValid(self: *SparseLU, basis: sparse_basis.SparseBasisView) SparseLuError!void {
        return self.factorizeImpl(basis, false, null);
    }

    /// Replay a previously recorded row/column pivot sequence. This is mainly
    /// a diagnostic control for comparing numerical kernels under identical
    /// ordering; normal reinversion should use threshold Markowitz selection.
    pub fn factorizeWithTraceAssumeValid(self: *SparseLU, basis: sparse_basis.SparseBasisView, trace: PivotTrace) SparseLuError!void {
        if (trace.rows.len != basis.dimension or trace.columns.len != basis.dimension) return error.DimensionMismatch;
        return self.factorizeImpl(basis, false, trace);
    }

    fn factorizeImpl(self: *SparseLU, basis: sparse_basis.SparseBasisView, comptime validate: bool, trace: ?PivotTrace) SparseLuError!void {
        const n = basis.dimension;
        try self.ensureDimension(n);
        try self.ensureFactorCapacity(@max(basis.nnz(), n));
        if (validate) try self.kernel.load(basis) else try self.kernel.loadAssumeValid(basis);
        self.dimension = n;
        self.l_nonzeros = 0;
        self.u_nonzeros = 0;
        self.inserted_fill = 0;
        self.peeled_pivots = 0;
        self.kernel_dimension = n;
        self.kernel_nonzeros = basis.nnz();
        self.selected_ordering = self.selectOrdering(basis);
        self.l_starts[0] = 0;
        self.u_starts[0] = 0;

        var peeling = trace == null;
        for (0..n) |pivot_index| {
            const choice = if (trace) |recorded|
                self.kernel.chooseRecordedPivot(recorded.rows[pivot_index], recorded.columns[pivot_index]) orelse return error.Singular
            else if (peeling) peel: {
                if (self.kernel.chooseSingleton()) |singleton| {
                    self.peeled_pivots += 1;
                    break :peel singleton;
                }
                peeling = false;
                self.kernel_dimension = n - pivot_index;
                self.kernel_nonzeros = self.kernel.activeEntries();
                break :peel switch (self.selected_ordering) {
                    .dod_markowitz => self.kernel.choosePivot(self.pivot_threshold) orelse return error.Singular,
                    .highs_kernel => self.kernel.choosePivotHighs(self.pivot_threshold) orelse return error.Singular,
                    .automatic => unreachable,
                };
            } else
                switch (self.selected_ordering) {
                    // The backend split is deliberately established before
                    // the HiGHS-style search lands, so dispatch/API changes
                    // can be verified independently from ordering changes.
                    .dod_markowitz => self.kernel.choosePivot(self.pivot_threshold) orelse return error.Singular,
                    .highs_kernel => self.kernel.choosePivotHighs(self.pivot_threshold) orelse return error.Singular,
                    .automatic => unreachable,
                };
            const pivot = try self.kernel.applyPivot(choice, self.zero_tolerance);
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
        if (peeling) {
            self.kernel_dimension = 0;
            self.kernel_nonzeros = 0;
        }
        self.hyper_views_ready = false;
    }

    /// Solve `B x = rhs` in place using the published row/column permutations.
    pub fn solve(self: *SparseLU, rhs: []f64) SparseLuError!void {
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

    /// Solve `B^T x = rhs` in place. U^T uses forward scatter and L^T uses a
    /// reverse gather over the same packed factor streams as FTRAN.
    pub fn solveTranspose(self: *SparseLU, rhs: []f64) SparseLuError!void {
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

    /// Solve a hyper-sparse FTRAN from the explicitly nonzero entries of
    /// `rhs`. The factor graph is traversed by scatter reachability; unrelated
    /// L/U rows are never visited. Only returned positions in `rhs` are valid.
    pub fn solveHyperSparse(self: *SparseLU, rhs: []f64, input_indices: []const u32, output_indices: []u32) SparseLuError!usize {
        const n = self.dimension;
        if (rhs.len != n or output_indices.len < n) return error.DimensionMismatch;
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
    /// the sequential packed kernel. The returned RHS is always dense.
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

    /// Retained requested bytes for the complete numerical kernel, factors,
    /// permutations and solve workspace. Allocator bookkeeping is excluded.
    pub fn requestedBytes(self: *const SparseLU) usize {
        var total = self.kernel.requestedBytes();
        inline for (.{
            "pivot_rows", "pivot_columns", "row_position", "column_position", "pivot_values",
            "l_starts", "u_starts", "l_rows", "l_values", "u_columns", "u_values", "work", "active", "hyper_output", "marked",
            "u_column_starts", "u_column_rows", "u_column_values", "l_row_starts", "l_row_columns", "l_row_values",
        }) |name| total += @sizeOf(std.meta.Elem(@TypeOf(@field(self, name)))) * @field(self, name).len;
        return total;
    }

    fn selectOrdering(self: *const SparseLU, _: sparse_basis.SparseBasisView) OrderingStrategy {
        return switch (self.ordering_strategy) {
            .automatic => .dod_markowitz,
            else => |forced| forced,
        };
    }

    fn ensureDimension(self: *SparseLU, required: usize) SparseLuError!void {
        if (required <= self.dimension_capacity) return;
        const capacity = grow(self.dimension_capacity, required) catch return error.CapacityOverflow;
        inline for (.{ "pivot_rows", "pivot_columns", "row_position", "column_position", "pivot_values", "work", "active", "hyper_output", "marked" }) |field_name|
            @field(self, field_name) = self.allocator.realloc(@field(self, field_name), capacity) catch return error.OutOfMemory;
        self.l_starts = self.allocator.realloc(self.l_starts, capacity + 1) catch return error.OutOfMemory;
        self.u_starts = self.allocator.realloc(self.u_starts, capacity + 1) catch return error.OutOfMemory;
        self.u_column_starts = self.allocator.realloc(self.u_column_starts, capacity + 1) catch return error.OutOfMemory;
        self.l_row_starts = self.allocator.realloc(self.l_row_starts, capacity + 1) catch return error.OutOfMemory;
        self.dimension_capacity = capacity;
    }

    fn ensureFactorCapacity(self: *SparseLU, required: usize) SparseLuError!void {
        if (required <= self.factor_capacity) return;
        const capacity = grow(self.factor_capacity, required) catch return error.CapacityOverflow;
        inline for (.{ "l_rows", "l_values", "u_columns", "u_values", "u_column_rows", "u_column_values", "l_row_columns", "l_row_values" }) |field_name|
            @field(self, field_name) = self.allocator.realloc(@field(self, field_name), capacity) catch return error.OutOfMemory;
        self.factor_capacity = capacity;
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
    fn accumulate(self: *SparseLU, position: u32, value: f64) void { self.activate(position, value); }
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
