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

pub const SparseLU = struct {
    allocator: std.mem.Allocator,
    kernel: sparse_kernel.MutableSparseKernel,
    dimension_capacity: usize = 0,
    factor_capacity: usize = 0,
    dimension: usize = 0,
    l_nonzeros: usize = 0,
    u_nonzeros: usize = 0,
    inserted_fill: usize = 0,

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
            "u_values",   "work",
        }) |field_name| self.allocator.free(@field(self, field_name));
        self.* = .{ .allocator = self.allocator, .kernel = sparse_kernel.MutableSparseKernel.init(self.allocator) };
    }

    pub fn factorize(self: *SparseLU, basis: sparse_basis.SparseBasisView) SparseLuError!void {
        return self.factorizeImpl(basis, true);
    }

    /// Zero-copy trusted reinversion entry for canonical engine-owned basis
    /// CSC. Workspace and factor capacities are retained across calls.
    pub fn factorizeAssumeValid(self: *SparseLU, basis: sparse_basis.SparseBasisView) SparseLuError!void {
        return self.factorizeImpl(basis, false);
    }

    fn factorizeImpl(self: *SparseLU, basis: sparse_basis.SparseBasisView, comptime validate: bool) SparseLuError!void {
        const n = basis.dimension;
        try self.ensureDimension(n);
        try self.ensureFactorCapacity(@max(basis.nnz(), n));
        if (validate) try self.kernel.load(basis) else try self.kernel.loadAssumeValid(basis);
        self.dimension = n;
        self.l_nonzeros = 0;
        self.u_nonzeros = 0;
        self.inserted_fill = 0;
        self.l_starts[0] = 0;
        self.u_starts[0] = 0;

        for (0..n) |pivot_index| {
            const choice = self.kernel.choosePivot(self.pivot_threshold) orelse return error.Singular;
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

    pub fn factorNonzeros(self: *const SparseLU) usize {
        return self.dimension + self.l_nonzeros + self.u_nonzeros;
    }

    /// Retained requested bytes for the complete numerical kernel, factors,
    /// permutations and solve workspace. Allocator bookkeeping is excluded.
    pub fn requestedBytes(self: *const SparseLU) usize {
        var total = self.kernel.requestedBytes();
        inline for (.{
            "pivot_rows", "pivot_columns", "row_position", "column_position", "pivot_values",
            "l_starts", "u_starts", "l_rows", "l_values", "u_columns", "u_values", "work",
        }) |name| total += @sizeOf(std.meta.Elem(@TypeOf(@field(self, name)))) * @field(self, name).len;
        return total;
    }

    fn ensureDimension(self: *SparseLU, required: usize) SparseLuError!void {
        if (required <= self.dimension_capacity) return;
        const capacity = grow(self.dimension_capacity, required) catch return error.CapacityOverflow;
        inline for (.{ "pivot_rows", "pivot_columns", "row_position", "column_position", "pivot_values", "work" }) |field_name|
            @field(self, field_name) = self.allocator.realloc(@field(self, field_name), capacity) catch return error.OutOfMemory;
        self.l_starts = self.allocator.realloc(self.l_starts, capacity + 1) catch return error.OutOfMemory;
        self.u_starts = self.allocator.realloc(self.u_starts, capacity + 1) catch return error.OutOfMemory;
        self.dimension_capacity = capacity;
    }

    fn ensureFactorCapacity(self: *SparseLU, required: usize) SparseLuError!void {
        if (required <= self.factor_capacity) return;
        const capacity = grow(self.factor_capacity, required) catch return error.CapacityOverflow;
        inline for (.{ "l_rows", "l_values", "u_columns", "u_values" }) |field_name|
            @field(self, field_name) = self.allocator.realloc(@field(self, field_name), capacity) catch return error.OutOfMemory;
        self.factor_capacity = capacity;
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
