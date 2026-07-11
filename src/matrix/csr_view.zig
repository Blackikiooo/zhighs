//! Rebuildable compressed sparse row (CSR) cache and borrowed views.
//!
//! CSC remains the only authoritative matrix representation. CSR exists solely
//! to accelerate row-oriented algorithms and may be discarded whenever the
//! model revision changes. Keeping this ownership direction explicit prevents
//! two independently mutable copies of the same mathematical matrix.

const std = @import("std");
const foundation = @import("foundation");
const sparse_vector = @import("sparse_vector.zig");
const csc = @import("csc.zig");

const RowId = foundation.RowId;
const ColId = foundation.ColId;

/// Non-owning CSR slices. A view is cheap to copy and must not outlive its cache.
pub const CsrView = struct {
    num_rows: usize,
    num_cols: usize,
    row_starts: []const usize,
    col_indices: []const ColId,
    values: []const f64,

    const Self = @This();

    pub inline fn nnz(self: Self) usize {
        return self.values.len;
    }

    /// Full structural validation for parser, conversion, and test boundaries.
    pub fn validate(self: Self) csc.MatrixError!void {
        try csc.validateDimensions(self.num_rows, self.num_cols);
        if (self.num_rows == std.math.maxInt(usize)) return error.DimensionTooLarge;
        if (self.row_starts.len != self.num_rows + 1) return error.InvalidRowStarts;
        if (self.col_indices.len != self.values.len) return error.InconsistentStorage;
        if (self.row_starts[0] != 0) return error.InvalidRowStarts;
        if (self.row_starts[self.num_rows] != self.nnz()) return error.InvalidRowStarts;

        for (0..self.num_rows) |row_index| {
            const begin = self.row_starts[row_index];
            const end = self.row_starts[row_index + 1];
            if (begin > end or end > self.nnz()) return error.InvalidRowStarts;

            var previous_col: ?usize = null;
            for (begin..end) |position| {
                const col = self.col_indices[position].toUsize();
                const value = self.values[position];
                if (col >= self.num_cols) return error.IndexOutOfBounds;
                if (previous_col) |previous| {
                    if (col <= previous) return error.IndicesNotStrictlyIncreasing;
                }
                if (!std.math.isFinite(value)) return error.NonFiniteValue;
                if (value == 0.0) return error.ExplicitZero;
                previous_col = col;
            }
        }
    }

    /// Returns one checked row as the common canonical sparse-vector view.
    pub fn row(self: Self, row_id: RowId) csc.MatrixError!sparse_vector.SparseVectorView(ColId) {
        const row_index = row_id.toUsize();
        if (row_index >= self.num_rows) return error.IndexOutOfBounds;
        return self.rowAssumeValid(row_index);
    }

    /// Hot-path row access after the cache and row index have been validated.
    pub inline fn rowAssumeValid(self: Self, row_index: usize) sparse_vector.SparseVectorView(ColId) {
        const begin = self.row_starts[row_index];
        const end = self.row_starts[row_index + 1];
        return .{
            .dimension = self.num_cols,
            .indices = self.col_indices[begin..end],
            .values = self.values[begin..end],
        };
    }
};

/// Owning, rebuildable CSR cache.
///
/// The revision is supplied by the owner of the authoritative CSC matrix. The
/// cache cannot infer mutations from raw slices, so revision ownership belongs
/// at Model level where structural changes are already coordinated.
pub const CsrCache = struct {
    source_revision: u64,
    num_rows: usize,
    num_cols: usize,
    row_starts: []usize,
    col_indices: []ColId,
    values: []f64,

    const Self = @This();

    /// Builds CSR in O(rows + nnz), without coordinate sorting.
    ///
    /// CSC columns are visited in ascending order. Appending each encountered
    /// column into its destination row therefore produces strictly increasing
    /// column IDs inside every CSR row automatically.
    pub fn build(allocator: std.mem.Allocator, matrix: csc.CscMatrix, source_revision: u64) (std.mem.Allocator.Error || csc.MatrixError)!Self {
        try matrix.validate();
        return buildAssumeValid(allocator, matrix, source_revision);
    }

    /// Conversion path for canonical CSC matrices from MatrixBuilder.freeze.
    /// It avoids a redundant O(nnz) structural scan while retaining allocation
    /// and dimension-overflow errors. The caller owns the validity proof.
    pub fn buildAssumeValid(allocator: std.mem.Allocator, matrix: csc.CscMatrix, source_revision: u64) (std.mem.Allocator.Error || csc.MatrixError)!Self {
        if (matrix.num_rows == std.math.maxInt(usize)) return error.DimensionTooLarge;

        const row_starts = try allocator.alloc(usize, matrix.num_rows + 1);
        errdefer allocator.free(row_starts);
        @memset(row_starts, 0);

        // Count row entries. Index + 1 leaves row_starts[0] as the zero prefix.
        for (matrix.row_indices) |row_id| row_starts[row_id.toUsize() + 1] += 1;
        for (0..matrix.num_rows) |row_index| row_starts[row_index + 1] += row_starts[row_index];

        const col_indices = try allocator.alloc(ColId, matrix.nnz());
        errdefer allocator.free(col_indices);
        const values = try allocator.alloc(f64, matrix.nnz());
        errdefer allocator.free(values);

        // One temporary cursor per row is the only conversion workspace. It is
        // released immediately; the persistent cache contains exactly CSR data.
        const next = try allocator.dupe(usize, row_starts[0..matrix.num_rows]);
        defer allocator.free(next);

        for (0..matrix.num_cols) |col_index| {
            // Matrix dimensions were validated above, so every column offset is
            // representable and cannot equal the reserved invalid sentinel.
            const col_id = ColId.fromUsize(col_index) catch unreachable;
            for (matrix.col_starts[col_index]..matrix.col_starts[col_index + 1]) |csc_position| {
                const row_index = matrix.row_indices[csc_position].toUsize();
                const csr_position = next[row_index];
                col_indices[csr_position] = col_id;
                values[csr_position] = matrix.values[csc_position];
                next[row_index] += 1;
            }
        }

        return .{
            .source_revision = source_revision,
            .num_rows = matrix.num_rows,
            .num_cols = matrix.num_cols,
            .row_starts = row_starts,
            .col_indices = col_indices,
            .values = values,
        };
    }

    pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
        allocator.free(self.row_starts);
        allocator.free(self.col_indices);
        allocator.free(self.values);
        self.* = undefined;
    }

    pub inline fn isCurrent(self: Self, current_revision: u64) bool {
        return self.source_revision == current_revision;
    }

    /// Checked access for component and model boundaries.
    pub fn view(self: Self, current_revision: u64) csc.MatrixError!CsrView {
        if (!self.isCurrent(current_revision)) return error.StaleView;
        return self.viewAssumeCurrent();
    }

    /// Zero-overhead access for loops whose caller already checked the revision.
    pub inline fn viewAssumeCurrent(self: Self) CsrView {
        return .{
            .num_rows = self.num_rows,
            .num_cols = self.num_cols,
            .row_starts = self.row_starts,
            .col_indices = self.col_indices,
            .values = self.values,
        };
    }

    /// Combines revision and row-bound checks for occasional row access.
    pub fn row(self: Self, current_revision: u64, row_id: RowId) csc.MatrixError!sparse_vector.SparseVectorView(ColId) {
        const current_view = try self.view(current_revision);
        return current_view.row(row_id);
    }
};

test "CSC to CSR conversion preserves rows including empty rows" {
    // Matrix rows: [2,0,-1], [0,0,0], [3,4,5].
    var starts = [_]usize{ 0, 2, 3, 5 };
    var rows = [_]RowId{ try RowId.init(0), try RowId.init(2), try RowId.init(2), try RowId.init(0), try RowId.init(2) };
    var values = [_]f64{ 2.0, 3.0, 4.0, -1.0, 5.0 };
    const matrix: csc.CscMatrix = .{
        .num_rows = 3,
        .num_cols = 3,
        .col_starts = &starts,
        .row_indices = &rows,
        .values = &values,
    };

    var cache = try CsrCache.build(std.testing.allocator, matrix, 17);
    defer cache.deinit(std.testing.allocator);
    const csr = try cache.view(17);
    try csr.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 2, 2, 5 }, csr.row_starts);
    try std.testing.expectEqualSlices(f64, &.{ 2.0, -1.0, 3.0, 4.0, 5.0 }, csr.values);

    const empty_row = try csr.row(try RowId.init(1));
    try std.testing.expectEqual(@as(usize, 0), empty_row.nnz());
    const last_row = try cache.row(17, try RowId.init(2));
    try std.testing.expectEqualSlices(f64, &.{ 3.0, 4.0, 5.0 }, last_row.values);
}

test "CSR cache rejects stale revision and invalid row" {
    var matrix = try csc.CscMatrix.initZero(std.testing.allocator, 2, 3);
    defer matrix.deinit(std.testing.allocator);
    var cache = try CsrCache.build(std.testing.allocator, matrix, 4);
    defer cache.deinit(std.testing.allocator);

    try std.testing.expect(cache.isCurrent(4));
    try std.testing.expect(!cache.isCurrent(5));
    try std.testing.expectError(error.StaleView, cache.view(5));
    try std.testing.expectError(error.StaleView, cache.row(5, try RowId.init(0)));
    try std.testing.expectError(error.IndexOutOfBounds, cache.row(4, try RowId.init(2)));
}

test "CSR validation rejects malformed offsets and unordered columns" {
    var starts = [_]usize{ 0, 2 };
    var duplicate_cols = [_]ColId{ try ColId.init(0), try ColId.init(0) };
    var values = [_]f64{ 1.0, 2.0 };
    const duplicate: CsrView = .{
        .num_rows = 1,
        .num_cols = 2,
        .row_starts = &starts,
        .col_indices = &duplicate_cols,
        .values = &values,
    };
    try std.testing.expectError(error.IndicesNotStrictlyIncreasing, duplicate.validate());

    var bad_starts = [_]usize{ 1, 1 };
    var no_cols = [_]ColId{};
    var no_values = [_]f64{};
    const bad: CsrView = .{
        .num_rows = 1,
        .num_cols = 1,
        .row_starts = &bad_starts,
        .col_indices = &no_cols,
        .values = &no_values,
    };
    try std.testing.expectError(error.InvalidRowStarts, bad.validate());
}

test "random CSC and CSR represent identical coordinates" {
    const builder_module = @import("builder.zig");
    var builder = try builder_module.MatrixBuilder.init(9, 7);
    defer builder.deinit(std.testing.allocator);
    var prng = std.Random.DefaultPrng.init(0xc5c5_5a17);
    const random = prng.random();
    for (0..250) |_| {
        const row = random.intRangeLessThan(usize, 0, 9);
        const col = random.intRangeLessThan(usize, 0, 7);
        const value: f64 = @floatFromInt(random.intRangeAtMost(i8, -4, 4));
        try builder.append(std.testing.allocator, try RowId.fromUsize(row), try ColId.fromUsize(col), value);
    }
    var matrix = try builder.freeze(std.testing.allocator, 0.0);
    defer matrix.deinit(std.testing.allocator);
    var cache = try CsrCache.build(std.testing.allocator, matrix, 0);
    defer cache.deinit(std.testing.allocator);
    const csr = cache.viewAssumeCurrent();
    try csr.validate();
    try std.testing.expectEqual(matrix.nnz(), csr.nnz());

    for (0..matrix.num_rows) |row_index| {
        const csr_row = csr.rowAssumeValid(row_index);
        for (csr_row.indices, csr_row.values) |col_id, value| {
            const column = matrix.columnAssumeValid(col_id.toUsize());
            var found = false;
            for (column.indices, column.values) |row_id, csc_value| {
                if (row_id.toUsize() == row_index) {
                    try std.testing.expectEqual(value, csc_value);
                    found = true;
                    break;
                }
            }
            try std.testing.expect(found);
        }
    }
}
