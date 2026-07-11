//! Canonical compressed sparse column (CSC) matrix storage.
//!
//! This is the read/compute representation used after construction. Columns
//! are contiguous; row indices in every column are strictly increasing; and
//! duplicate coordinates and explicit zeros are forbidden.

const std = @import("std");
const foundation = @import("foundation");
const sparse_vector = @import("sparse_vector.zig");

const RowId = foundation.RowId;
const ColId = foundation.ColId;

pub const MatrixError = sparse_vector.SparseVectorError || error{
    /// The matrix has a nonzero value at a row index outside its number of rows.
    DimensionTooLarge,
    /// The matrix has a nonzero value at a column index outside its number of columns.
    InvalidColumnStarts,
    /// The matrix has a nonzero value at a row index outside its number of rows.
    InconsistentStorage,
    /// A cached matrix view was built for an older model revision.
    StaleView,
    /// CSR row offsets do not describe the entry arrays correctly.
    InvalidRowStarts,
    /// A row or column scaling factor is non-finite or not strictly positive.
    InvalidScaling,
    /// A numerical transformation produced zero or a non-finite matrix entry.
    NumericalOverflow,
    /// A permutation is not a bijection over the requested dimension.
    InvalidPermutation,
    /// A dynamic-row rollback checkpoint does not belong to current storage.
    InvalidCheckpoint,
};

/// Owning CSC matrix. Column j occupies the half-open range described by
/// col_starts[j] and col_starts[j + 1].
pub const CscMatrix = struct {
    num_rows: usize,
    num_cols: usize,
    col_starts: []usize,
    row_indices: []RowId,
    values: []f64,

    const Self = @This();

    /// Creates a structurally valid all-zero matrix.
    pub fn initZero(allocator: std.mem.Allocator, num_rows: usize, num_cols: usize) (std.mem.Allocator.Error || MatrixError)!Self {
        try validateDimensions(num_rows, num_cols);
        const col_starts = try allocator.alloc(usize, num_cols + 1);
        errdefer allocator.free(col_starts);
        @memset(col_starts, 0);
        const row_indices = try allocator.alloc(RowId, 0);
        errdefer allocator.free(row_indices);
        const values = try allocator.alloc(f64, 0);
        return .{
            .num_rows = num_rows,
            .num_cols = num_cols,
            .col_starts = col_starts,
            .row_indices = row_indices,
            .values = values,
        };
    }

    pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
        allocator.free(self.col_starts);
        allocator.free(self.row_indices);
        allocator.free(self.values);
        self.* = undefined;
    }

    /// Returns the number of nonzero values in the matrix.
    pub inline fn nnz(self: Self) usize {
        return self.values.len;
    }

    /// Performs complete validation at API, parser, or test boundaries.
    /// Matrices emitted by MatrixBuilder.freeze need not be rechecked in every
    /// numerical kernel.
    pub fn validate(self: Self) MatrixError!void {
        try validateDimensions(self.num_rows, self.num_cols);
        if (self.col_starts.len != self.num_cols + 1) return error.InvalidColumnStarts;
        if (self.row_indices.len != self.values.len) return error.InconsistentStorage;
        if (self.col_starts[0] != 0) return error.InvalidColumnStarts;
        if (self.col_starts[self.num_cols] != self.nnz()) return error.InvalidColumnStarts;

        for (0..self.num_cols) |col| {
            const begin = self.col_starts[col];
            const end = self.col_starts[col + 1];
            if (begin > end or end > self.nnz()) return error.InvalidColumnStarts;

            var previous_row: ?usize = null;
            for (begin..end) |position| {
                const row = self.row_indices[position].toUsize();
                const value = self.values[position];
                if (row >= self.num_rows) return error.IndexOutOfBounds;
                if (previous_row) |previous| {
                    if (row <= previous) return error.IndicesNotStrictlyIncreasing;
                }
                if (!std.math.isFinite(value)) return error.NonFiniteValue;
                if (value == 0.0) return error.ExplicitZero;
                previous_row = row;
            }
        }
    }

    /// Returns a checked, non-owning view of one column.
    pub fn column(self: Self, col: ColId) MatrixError!sparse_vector.SparseVectorView(RowId) {
        const col_index = col.toUsize();
        if (col_index >= self.num_cols) return error.IndexOutOfBounds;
        return self.columnAssumeValid(col_index);
    }

    /// Caller guarantees a valid CSC matrix and col_index < num_cols.
    pub inline fn columnAssumeValid(self: Self, col_index: usize) sparse_vector.SparseVectorView(RowId) {
        const begin = self.col_starts[col_index];
        const end = self.col_starts[col_index + 1];
        return .{
            .dimension = self.num_rows,
            .indices = self.row_indices[begin..end],
            .values = self.values[begin..end],
        };
    }

    /// Computes y = A * x, checking only vector dimensions.
    pub fn multiply(self: Self, x: []const f64, y: []f64) MatrixError!void {
        if (x.len != self.num_cols or y.len != self.num_rows) return error.DimensionMismatch;
        self.multiplyAssumeValid(x, y);
    }

    /// Hot-path version: dimensions and CSC invariants must already be valid.
    pub fn multiplyAssumeValid(self: Self, x: []const f64, y: []f64) void {
        @memset(y, 0.0);
        for (0..self.num_cols) |col| {
            const x_value = x[col];
            // Exact zeros are frequent in dense simplex work vectors.
            if (x_value == 0.0) continue;
            for (self.col_starts[col]..self.col_starts[col + 1]) |position| {
                const row = self.row_indices[position].toUsize();
                y[row] += self.values[position] * x_value;
            }
        }
    }

    /// Computes y = transpose(A) * x, checking only vector dimensions.
    pub fn transposeMultiply(self: Self, x: []const f64, y: []f64) MatrixError!void {
        if (x.len != self.num_rows or y.len != self.num_cols) return error.DimensionMismatch;
        self.transposeMultiplyAssumeValid(x, y);
    }

    /// Hot-path transpose multiply with prevalidated dimensions and structure.
    pub fn transposeMultiplyAssumeValid(self: Self, x: []const f64, y: []f64) void {
        for (0..self.num_cols) |col| {
            var sum: f64 = 0.0;
            for (self.col_starts[col]..self.col_starts[col + 1]) |position| {
                sum += self.values[position] * x[self.row_indices[position].toUsize()];
            }
            y[col] = sum;
        }
    }
};

/// Checks whether dimensions can be represented by the strong ID types.
pub fn validateDimensions(num_rows: usize, num_cols: usize) MatrixError!void {
    // CSC needs num_cols + 1 offsets; reject the sole overflowing usize value.
    if (num_cols == std.math.maxInt(usize)) return error.DimensionTooLarge;
    if (num_rows != 0) _ = RowId.fromUsize(num_rows - 1) catch return error.DimensionTooLarge;
    if (num_cols != 0) _ = ColId.fromUsize(num_cols - 1) catch return error.DimensionTooLarge;
}

test "CSC zero matrix is valid and multiplies to zero" {
    var matrix = try CscMatrix.initZero(std.testing.allocator, 3, 2);
    defer matrix.deinit(std.testing.allocator);
    try matrix.validate();
    try std.testing.expectEqual(@as(usize, 0), matrix.nnz());
    var y = [_]f64{ 9.0, 9.0, 9.0 };
    try matrix.multiply(&.{ 4.0, -2.0 }, &y);
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 0.0, 0.0 }, &y);
}

test "CSC column views and both multiplication directions" {
    // Matrix rows: [2,0,-1], [0,4,0], [3,0,5].
    var starts = [_]usize{ 0, 2, 3, 5 };
    var rows = [_]RowId{ try RowId.init(0), try RowId.init(2), try RowId.init(1), try RowId.init(0), try RowId.init(2) };
    var values = [_]f64{ 2.0, 3.0, 4.0, -1.0, 5.0 };
    const matrix: CscMatrix = .{ .num_rows = 3, .num_cols = 3, .col_starts = &starts, .row_indices = &rows, .values = &values };
    try matrix.validate();
    const second = try matrix.column(try ColId.init(1));
    try std.testing.expectEqual(@as(usize, 1), second.nnz());
    try std.testing.expectEqual(@as(usize, 1), second.indices[0].toUsize());
    var ax: [3]f64 = undefined;
    try matrix.multiply(&.{ 2.0, 3.0, -1.0 }, &ax);
    try std.testing.expectEqualSlices(f64, &.{ 5.0, 12.0, 1.0 }, &ax);
    var atx: [3]f64 = undefined;
    try matrix.transposeMultiply(&.{ 1.0, 2.0, 3.0 }, &atx);
    try std.testing.expectEqualSlices(f64, &.{ 11.0, 8.0, 14.0 }, &atx);
}

test "CSC rejects malformed structure and values" {
    var bad_starts = [_]usize{ 1, 1 };
    var no_rows = [_]RowId{};
    var no_values = [_]f64{};
    const bad_start: CscMatrix = .{ .num_rows = 1, .num_cols = 1, .col_starts = &bad_starts, .row_indices = &no_rows, .values = &no_values };
    try std.testing.expectError(error.InvalidColumnStarts, bad_start.validate());

    var starts = [_]usize{ 0, 2 };
    var duplicate_rows = [_]RowId{ try RowId.init(0), try RowId.init(0) };
    var values = [_]f64{ 1.0, 2.0 };
    const duplicate: CscMatrix = .{ .num_rows = 2, .num_cols = 1, .col_starts = &starts, .row_indices = &duplicate_rows, .values = &values };
    try std.testing.expectError(error.IndicesNotStrictlyIncreasing, duplicate.validate());

    var one_start = [_]usize{ 0, 1 };
    var one_row = [_]RowId{try RowId.init(0)};
    var zero = [_]f64{0.0};
    const explicit_zero: CscMatrix = .{ .num_rows = 1, .num_cols = 1, .col_starts = &one_start, .row_indices = &one_row, .values = &zero };
    try std.testing.expectError(error.ExplicitZero, explicit_zero.validate());
}

test "CSC checked access rejects incompatible dimensions" {
    var matrix = try CscMatrix.initZero(std.testing.allocator, 2, 3);
    defer matrix.deinit(std.testing.allocator);
    var y2: [2]f64 = undefined;
    var y3: [3]f64 = undefined;
    try std.testing.expectError(error.DimensionMismatch, matrix.multiply(&.{1.0}, &y2));
    try std.testing.expectError(error.DimensionMismatch, matrix.transposeMultiply(&.{1.0}, &y3));
    try std.testing.expectError(error.IndexOutOfBounds, matrix.column(try ColId.init(3)));
}
