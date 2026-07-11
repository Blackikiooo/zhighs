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
const memory = @import("memory.zig");

const RowId = foundation.RowId;
const ColId = foundation.ColId;

/// Non-owning CSR slices. A view is cheap to copy and must not outlive its cache.
pub const CsrView = struct {
    num_rows: usize,
    num_cols: usize,
    row_starts: []const foundation.HUInt,
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
            const begin: usize = @intCast(self.row_starts[row_index]);
            const end: usize = @intCast(self.row_starts[row_index + 1]);
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
        const begin: usize = @intCast(self.row_starts[row_index]);
        const end: usize = @intCast(self.row_starts[row_index + 1]);
        return .{
            .dimension = self.num_cols,
            .indices = self.col_indices[begin..end],
            .values = self.values[begin..end],
        };
    }

    /// CSR-native y = A*x. Each output row is an independent contiguous dot
    /// product, avoiding the random output scatter required by CSC.
    pub fn multiply(self: Self, x: []const f64, y: []f64) csc.MatrixError!void {
        if (x.len != self.num_cols or y.len != self.num_rows) return error.DimensionMismatch;
        self.multiplyAssumeValid(x, y);
    }

    pub fn multiplyAssumeValid(self: Self, x: []const f64, y: []f64) void {
        const nrow = self.num_rows;
        const rs = self.row_starts;
        const ci = self.col_indices;
        const vs = self.values;
        var row_idx: usize = 0;
        while (row_idx < nrow) : (row_idx += 1) {
            var sum: f64 = 0.0;
            var pos: usize = @intCast(rs[row_idx]);
            const end: usize = @intCast(rs[row_idx + 1]);
            while (pos < end) : (pos += 1)
                sum = @mulAdd(f64, vs[pos], x[ci[pos].toUsize()], sum);
            y[row_idx] = sum;
        }
    }

    /// CSR-native y = transpose(A)*x. CSC is normally faster for this direction;
    /// this API avoids a format conversion when only CSR is available.
    pub fn transposeMultiply(self: Self, x: []const f64, y: []f64) csc.MatrixError!void {
        if (x.len != self.num_rows or y.len != self.num_cols) return error.DimensionMismatch;
        self.transposeMultiplyAssumeValid(x, y);
    }

    pub fn transposeMultiplyAssumeValid(self: Self, x: []const f64, y: []f64) void {
        const nrow = self.num_rows;
        const rs = self.row_starts;
        const ci = self.col_indices;
        const vs = self.values;
        memory.clearF64(y);
        var i: usize = 0;
        while (i < nrow) : (i += 1) {
            const multiplier = x[i];
            var pos: usize = @intCast(rs[i]);
            const end: usize = @intCast(rs[i + 1]);
            while (pos < end) : (pos += 1) {
                const col = ci[pos].toUsize();
                y[col] = @mulAdd(f64, vs[pos], multiplier, y[col]);
            }
        }
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
    row_starts: []foundation.HUInt,
    col_indices: []ColId,
    values: []f64,
    storage: []align(64) u8,

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
        const next = try allocator.alloc(foundation.HUInt, matrix.num_rows);
        defer allocator.free(next);
        return buildWithScratchAssumeValid(allocator, matrix, source_revision, next);
    }

    pub fn buildWithScratch(allocator: std.mem.Allocator, matrix: csc.CscMatrix, source_revision: u64, row_cursor_scratch: []foundation.HUInt) (std.mem.Allocator.Error || csc.MatrixError)!Self {
        try matrix.validate();
        return buildWithScratchAssumeValid(allocator, matrix, source_revision, row_cursor_scratch);
    }

    /// Reuses caller-owned row cursor storage and removes one temporary
    /// allocation from repeated CSC-to-CSR conversions.
    pub fn buildWithScratchAssumeValid(allocator: std.mem.Allocator, matrix: csc.CscMatrix, source_revision: u64, row_cursor_scratch: []foundation.HUInt) (std.mem.Allocator.Error || csc.MatrixError)!Self {
        if (matrix.num_rows == std.math.maxInt(usize)) return error.DimensionTooLarge;
        if (row_cursor_scratch.len < matrix.num_rows) return error.DimensionMismatch;
        if (matrix.nnz() > std.math.maxInt(foundation.HUInt)) return error.DimensionTooLarge;

        // One allocation both removes allocator traffic and gives deterministic
        // page colors. Each concurrently accessed stream gets a different
        // 64-byte offset within its 4 KiB page, avoiding address-alias stalls.
        const starts_bytes = std.math.mul(usize, matrix.num_rows + 1, @sizeOf(foundation.HUInt)) catch return error.DimensionTooLarge;
        const indices_bytes = std.math.mul(usize, matrix.nnz(), @sizeOf(ColId)) catch return error.DimensionTooLarge;
        const values_bytes = std.math.mul(usize, matrix.nnz(), @sizeOf(f64)) catch return error.DimensionTooLarge;
        const starts_offset: usize = 0;
        const indices_offset = std.mem.alignForward(usize, starts_bytes, 4096) + 64;
        const values_offset = std.mem.alignForward(usize, indices_offset + indices_bytes, 4096) + 128;
        const storage_len = std.math.add(usize, values_offset, values_bytes) catch return error.DimensionTooLarge;
        const storage = try allocator.alignedAlloc(u8, .@"64", storage_len);
        errdefer allocator.free(storage);
        const row_starts_ptr: [*]foundation.HUInt = @ptrCast(@alignCast(storage.ptr + starts_offset));
        const col_indices_ptr: [*]ColId = @ptrCast(@alignCast(storage.ptr + indices_offset));
        const values_ptr: [*]f64 = @ptrCast(@alignCast(storage.ptr + values_offset));
        const row_starts = row_starts_ptr[0 .. matrix.num_rows + 1];
        const col_indices = col_indices_ptr[0..matrix.nnz()];
        const values = values_ptr[0..matrix.nnz()];

        try fillFromCscAssumeValid(matrix, row_starts, col_indices, values, row_cursor_scratch);

        return .{
            .source_revision = source_revision,
            .num_rows = matrix.num_rows,
            .num_cols = matrix.num_cols,
            .row_starts = row_starts,
            .col_indices = col_indices,
            .values = values,
            .storage = storage,
        };
    }

    pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
        allocator.free(self.storage);
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

/// Allocation-free CSC-to-CSR conversion into exact-size caller buffers.
pub fn fillFromCsc(matrix: csc.CscMatrix, row_starts: []foundation.HUInt, col_indices: []ColId, values: []f64, row_cursor_scratch: []foundation.HUInt) csc.MatrixError!void {
    try matrix.validate();
    if (matrix.nnz() > std.math.maxInt(foundation.HUInt)) return error.DimensionTooLarge;
    return fillFromCscAssumeValid(matrix, row_starts, col_indices, values, row_cursor_scratch);
}

pub fn fillFromCscAssumeValid(matrix: csc.CscMatrix, row_starts: []foundation.HUInt, col_indices: []ColId, values: []f64, row_cursor_scratch: []foundation.HUInt) csc.MatrixError!void {
    if (row_starts.len != matrix.num_rows + 1 or col_indices.len != matrix.nnz() or
        values.len != matrix.nnz() or row_cursor_scratch.len < matrix.num_rows)
        return error.DimensionMismatch;
    memory.clearTyped(foundation.HUInt, row_starts);
    for (matrix.row_indices) |row_id| row_starts[row_id.toUsize() + 1] += 1;
    for (0..matrix.num_rows) |row_index| row_starts[row_index + 1] += row_starts[row_index];
    const next = row_cursor_scratch[0..matrix.num_rows];
    @memcpy(next, row_starts[0..matrix.num_rows]);
    for (0..matrix.num_cols) |col_index| {
        const col_id = ColId.fromUsize(col_index) catch unreachable;
        for (matrix.col_starts[col_index]..matrix.col_starts[col_index + 1]) |position| {
            // Load the sequential source value before following the dependent
            // row -> cursor -> destination chain, giving the CPU useful work
            // while the cursor load is outstanding.
            const source_value = matrix.values[position];
            const row = matrix.row_indices[position].toUsize();
            const destination: usize = @intCast(next[row]);
            next[row] += 1;
            col_indices[destination] = col_id;
            values[destination] = source_value;
        }
    }
}

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
    try std.testing.expectEqualSlices(foundation.HUInt, &.{ 0, 2, 2, 5 }, csr.row_starts);
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
    var starts = [_]foundation.HUInt{ 0, 2 };
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

    var bad_starts = [_]foundation.HUInt{ 1, 1 };
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

test "CSR native products match CSC products and scratch build" {
    var starts = [_]usize{ 0, 2, 3 };
    var rows = [_]RowId{ try RowId.init(0), try RowId.init(1), try RowId.init(1) };
    var values = [_]f64{ 2.0, 3.0, 4.0 };
    const matrix: csc.CscMatrix = .{ .num_rows = 2, .num_cols = 2, .col_starts = &starts, .row_indices = &rows, .values = &values };
    var scratch: [2]foundation.HUInt = undefined;
    var cache = try CsrCache.buildWithScratch(std.testing.allocator, matrix, 0, &scratch);
    defer cache.deinit(std.testing.allocator);
    const csr = cache.viewAssumeCurrent();
    var csc_y: [2]f64 = undefined;
    var csr_y: [2]f64 = undefined;
    try matrix.multiply(&.{ 5.0, 6.0 }, &csc_y);
    try csr.multiply(&.{ 5.0, 6.0 }, &csr_y);
    try std.testing.expectEqualSlices(f64, &csc_y, &csr_y);
    try matrix.transposeMultiply(&.{ 5.0, 6.0 }, &csc_y);
    try csr.transposeMultiply(&.{ 5.0, 6.0 }, &csr_y);
    try std.testing.expectEqualSlices(f64, &csc_y, &csr_y);
}
