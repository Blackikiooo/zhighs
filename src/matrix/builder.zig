//! Mutable sparse matrix construction in Structure of Arrays (SoA) form.
//!
//! Construction and computation have different optimal layouts. This builder
//! accepts coordinates in arbitrary order, then freeze sorts, merges duplicates,
//! removes numerical zeros, and emits canonical CSC storage. MultiArrayList
//! owns one allocation while exposing a separate contiguous slice per field,
//! preventing the parallel-array lengths from ever becoming inconsistent.

const std = @import("std");
const foundation = @import("foundation");
const memory = @import("memory.zig");
const sparse_vector = @import("sparse_vector.zig");
const csc = @import("csc.zig");

const RowId = foundation.RowId;
const ColId = foundation.ColId;

const Triplet = struct {
    row: RowId,
    col: ColId,
    value: f64,
    // Makes (col,row,sequence) a total order so unstable PDQ sorting still
    // preserves deterministic duplicate summation order.
    sequence: usize,
};

const TripletList = std.MultiArrayList(Triplet);

pub const MatrixBuilder = struct {
    num_rows: usize,
    num_cols: usize,

    // Official Zig SoA container: row, column and value fields are contiguous,
    // but capacity and logical length are managed atomically as one collection.
    entries: TripletList = .empty,

    const Self = @This();

    pub fn init(num_rows: usize, num_cols: usize) csc.MatrixError!Self {
        try csc.validateDimensions(num_rows, num_cols);
        return .{ .num_rows = num_rows, .num_cols = num_cols };
    }

    pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
        self.entries.deinit(allocator);
        self.* = undefined;
    }

    pub inline fn len(self: Self) usize {
        return self.entries.len;
    }

    pub fn clearRetainingCapacity(self: *Self) void {
        self.entries.clearRetainingCapacity();
    }

    /// Reserves the single MultiArrayList allocation before mutation.
    pub fn reserve(self: *Self, allocator: std.mem.Allocator, additional: usize) std.mem.Allocator.Error!void {
        try self.entries.ensureUnusedCapacity(allocator, additional);
    }

    /// Appends one coordinate. Explicit zeros are accepted here because later
    /// duplicate entries may cancel; freeze removes the final numerical zero.
    pub fn append(self: *Self, allocator: std.mem.Allocator, row: RowId, col: ColId, value: f64) (std.mem.Allocator.Error || csc.MatrixError)!void {
        if (row.toUsize() >= self.num_rows or col.toUsize() >= self.num_cols)
            return error.IndexOutOfBounds;
        if (!std.math.isFinite(value)) return error.NonFiniteValue;

        try self.reserve(allocator, 1);
        self.entries.appendAssumeCapacity(.{ .row = row, .col = col, .value = value, .sequence = self.len() });
    }

    /// Appends a canonical sparse column in one reserved batch.
    pub fn appendColumn(self: *Self, allocator: std.mem.Allocator, col: ColId, vector: sparse_vector.SparseVectorView(RowId)) (std.mem.Allocator.Error || csc.MatrixError)!void {
        if (col.toUsize() >= self.num_cols) return error.IndexOutOfBounds;
        if (vector.dimension != self.num_rows) return error.DimensionMismatch;
        try vector.validate();

        try self.reserve(allocator, vector.nnz());
        for (vector.indices, vector.values) |row, value| {
            self.entries.appendAssumeCapacity(.{ .row = row, .col = col, .value = value, .sequence = self.len() });
        }
    }

    /// Appends a canonical sparse row in one reserved batch.
    pub fn appendRow(self: *Self, allocator: std.mem.Allocator, row: RowId, vector: sparse_vector.SparseVectorView(ColId)) (std.mem.Allocator.Error || csc.MatrixError)!void {
        if (row.toUsize() >= self.num_rows) return error.IndexOutOfBounds;
        if (vector.dimension != self.num_cols) return error.DimensionMismatch;
        try vector.validate();

        try self.reserve(allocator, vector.nnz());
        for (vector.indices, vector.values) |col, value| {
            self.entries.appendAssumeCapacity(.{ .row = row, .col = col, .value = value, .sequence = self.len() });
        }
    }

    /// Canonicalizes the builder and returns an independent owning CSC matrix.
    /// Values with abs(value) <= zero_tolerance are omitted after duplicates are
    /// summed. The builder remains reusable and holds the compacted triplets.
    pub fn freeze(self: *Self, allocator: std.mem.Allocator, zero_tolerance: f64) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
        return self.freezeInternal(allocator, zero_tolerance, true);
    }

    /// Linear-time freeze for callers that already emit nondecreasing
    /// (column,row) coordinates. A linear order check is retained.
    pub fn freezeSorted(self: *Self, allocator: std.mem.Allocator, zero_tolerance: f64) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
        if (!self.isSorted()) return error.TripletsNotSorted;
        return self.freezeInternal(allocator, zero_tolerance, false);
    }

    /// Fastest construction path: caller guarantees nondecreasing coordinates.
    /// Duplicate coordinates may remain adjacent and are merged stably.
    pub fn freezeSortedAssumeValid(self: *Self, allocator: std.mem.Allocator, zero_tolerance: f64) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
        return self.freezeInternal(allocator, zero_tolerance, false);
    }

    fn freezeInternal(self: *Self, allocator: std.mem.Allocator, zero_tolerance: f64, sort_entries: bool) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
        if (!std.math.isFinite(zero_tolerance) or zero_tolerance < 0.0)
            return error.InvalidTolerance;
        const length = self.len();
        // Zig 0.16 MultiArrayList stable context sort is insertion sort. Use
        // PDQ plus sequence as a total tie-break key for O(n log n) behavior
        // without sacrificing deterministic duplicate summation.
        if (sort_entries) self.entries.sortUnstable(SortContext{ .fields = self.entries.slice() });
        var fields = self.entries.slice();
        const rows = fields.items(.row);
        const cols = fields.items(.col);
        const values = fields.items(.value);

        // Merge in place. Read positions never trail the write position, so no
        // temporary triplet allocation is needed after sorting.
        @setFloatMode(.optimized);
        var read: usize = 0;
        var write: usize = 0;
        while (read < length) {
            const row = rows[read];
            const col = cols[read];
            var sum = values[read];
            read += 1;
            while (read < length and sameCoordinate(rows[read], cols[read], row, col)) : (read += 1) {
                sum += values[read];
            }
            if (!std.math.isFinite(sum)) return error.NonFiniteValue;
            if (@abs(sum) <= zero_tolerance) continue;

            fields.set(write, .{ .row = row, .col = col, .value = sum, .sequence = write });
            write += 1;
        }
        self.entries.shrinkRetainingCapacity(write);
        fields = self.entries.slice();
        const compact_rows = fields.items(.row);
        const compact_cols = fields.items(.col);
        const compact_values = fields.items(.value);

        // Allocate each output stream transactionally. Separate errdefer guards
        // make every allocation-failure point leak- and double-free-safe.
        const col_starts = try allocator.alloc(usize, self.num_cols + 1);
        errdefer allocator.free(col_starts);
        @memset(col_starts, 0);
        const row_indices = try allocator.dupe(RowId, compact_rows);
        errdefer allocator.free(row_indices);
        const output_values = try allocator.dupe(f64, compact_values);
        errdefer allocator.free(output_values);

        // Count entries per column, then prefix-sum counts into CSC offsets.
        const ncol = self.num_cols;
        for (compact_cols) |col| col_starts[col.toUsize() + 1] += 1;
        var c: usize = 0;
        while (c < ncol) : (c += 1) col_starts[c + 1] += col_starts[c];
        return .{
            .num_rows = self.num_rows,
            .num_cols = self.num_cols,
            .col_starts = col_starts,
            .row_indices = row_indices,
            .values = output_values,
        };
    }

    fn isSorted(self: Self) bool {
        const fields = self.entries.slice();
        const rows = fields.items(.row);
        const cols = fields.items(.col);
        for (1..self.entries.len) |index| {
            const previous_col = cols[index - 1].toUsize();
            const current_col = cols[index].toUsize();
            if (current_col < previous_col) return false;
            if (current_col == previous_col and rows[index].toUsize() < rows[index - 1].toUsize()) return false;
        }
        return true;
    }
};

inline fn sameCoordinate(a_row: RowId, a_col: ColId, b_row: RowId, b_col: ColId) bool {
    return a_row.toUsize() == b_row.toUsize() and a_col.toUsize() == b_col.toUsize();
}

/// Comparator reads only the coordinate streams. MultiArrayList moves all
/// fields together internally during its stable sort.
const SortContext = struct {
    fields: TripletList.Slice,

    pub fn lessThan(self: @This(), a: usize, b: usize) bool {
        const cols = self.fields.items(.col);
        const a_col = cols[a].toUsize();
        const b_col = cols[b].toUsize();
        if (a_col != b_col) return a_col < b_col;
        const rows = self.fields.items(.row);
        const a_row = rows[a].toUsize();
        const b_row = rows[b].toUsize();
        if (a_row != b_row) return a_row < b_row;
        const sequences = self.fields.items(.sequence);
        return sequences[a] < sequences[b];
    }
};

test "SoA builder sorts, merges duplicates, and preserves empty columns" {
    var builder = try MatrixBuilder.init(3, 4);
    defer builder.deinit(std.testing.allocator);
    try builder.append(std.testing.allocator, try RowId.init(2), try ColId.init(3), 5.0);
    try builder.append(std.testing.allocator, try RowId.init(0), try ColId.init(0), 2.0);
    try builder.append(std.testing.allocator, try RowId.init(1), try ColId.init(3), 4.0);
    try builder.append(std.testing.allocator, try RowId.init(2), try ColId.init(3), -1.0);
    try builder.append(std.testing.allocator, try RowId.init(1), try ColId.init(1), 0.0);

    var matrix = try builder.freeze(std.testing.allocator, 0.0);
    defer matrix.deinit(std.testing.allocator);
    try matrix.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 1, 1, 1, 3 }, matrix.col_starts);
    try std.testing.expectEqualSlices(f64, &.{ 2.0, 4.0, 4.0 }, matrix.values);
    try std.testing.expectEqual(@as(usize, 1), matrix.row_indices[1].toUsize());
    try std.testing.expectEqual(@as(usize, 2), matrix.row_indices[2].toUsize());
}

test "builder tolerance removes cancellation and small residuals" {
    var builder = try MatrixBuilder.init(2, 2);
    defer builder.deinit(std.testing.allocator);
    try builder.append(std.testing.allocator, try RowId.init(0), try ColId.init(0), 3.0);
    try builder.append(std.testing.allocator, try RowId.init(0), try ColId.init(0), -3.0);
    try builder.append(std.testing.allocator, try RowId.init(1), try ColId.init(1), 1e-9);
    var matrix = try builder.freeze(std.testing.allocator, 1e-8);
    defer matrix.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 0), matrix.nnz());
    try matrix.validate();
}

test "builder batch-appends canonical rows and columns" {
    var builder = try MatrixBuilder.init(3, 3);
    defer builder.deinit(std.testing.allocator);
    var column_rows = [_]RowId{ try RowId.init(0), try RowId.init(2) };
    var column_values = [_]f64{ 2.0, 6.0 };
    try builder.appendColumn(std.testing.allocator, try ColId.init(1), .{ .dimension = 3, .indices = &column_rows, .values = &column_values });
    var row_cols = [_]ColId{ try ColId.init(0), try ColId.init(2) };
    var row_values = [_]f64{ 3.0, 7.0 };
    try builder.appendRow(std.testing.allocator, try RowId.init(1), .{ .dimension = 3, .indices = &row_cols, .values = &row_values });
    var matrix = try builder.freeze(std.testing.allocator, 0.0);
    defer matrix.deinit(std.testing.allocator);
    try matrix.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 1, 3, 4 }, matrix.col_starts);
    try std.testing.expectEqualSlices(f64, &.{ 3.0, 2.0, 6.0, 7.0 }, matrix.values);
}

test "builder validates coordinates values tolerance and batch dimensions" {
    var builder = try MatrixBuilder.init(1, 1);
    defer builder.deinit(std.testing.allocator);
    try std.testing.expectError(error.IndexOutOfBounds, builder.append(std.testing.allocator, try RowId.init(1), try ColId.init(0), 1.0));
    try std.testing.expectError(error.NonFiniteValue, builder.append(std.testing.allocator, try RowId.init(0), try ColId.init(0), std.math.nan(f64)));
    try std.testing.expectError(error.InvalidTolerance, builder.freeze(std.testing.allocator, -1.0));
    var no_rows = [_]RowId{};
    var no_values = [_]f64{};
    try std.testing.expectError(error.DimensionMismatch, builder.appendColumn(std.testing.allocator, try ColId.init(0), .{ .dimension = 2, .indices = &no_rows, .values = &no_values }));
}

test "duplicate summation follows insertion order" {
    var builder = try MatrixBuilder.init(1, 1);
    defer builder.deinit(std.testing.allocator);
    const row = try RowId.init(0);
    const col = try ColId.init(0);
    try builder.append(std.testing.allocator, row, col, 1e16);
    try builder.append(std.testing.allocator, row, col, 1.0);
    try builder.append(std.testing.allocator, row, col, -1e16);
    var matrix = try builder.freeze(std.testing.allocator, 0.0);
    defer matrix.deinit(std.testing.allocator);
    // Sequential IEEE-754 evaluation gives zero and canonicalization drops it.
    try std.testing.expectEqual(@as(usize, 0), matrix.nnz());
}

test "randomized CSC multiply agrees with dense reference" {
    const row_count = 8;
    const col_count = 7;
    var prng = std.Random.DefaultPrng.init(0x5a17_c5c5);
    const random = prng.random();

    // Integer-valued samples keep this structural property test exact: any
    // mismatch points to indexing/canonicalization, not tolerance selection.
    var dense = [_]f64{0.0} ** (row_count * col_count);
    var builder = try MatrixBuilder.init(row_count, col_count);
    defer builder.deinit(std.testing.allocator);
    try builder.reserve(std.testing.allocator, 300);

    for (0..300) |_| {
        const row = random.intRangeLessThan(usize, 0, row_count);
        const col = random.intRangeLessThan(usize, 0, col_count);
        const integer = random.intRangeAtMost(i8, -5, 5);
        const value: f64 = @floatFromInt(integer);
        try builder.append(std.testing.allocator, try RowId.fromUsize(row), try ColId.fromUsize(col), value);
        dense[row * col_count + col] += value;
    }

    var matrix = try builder.freeze(std.testing.allocator, 0.0);
    defer matrix.deinit(std.testing.allocator);
    try matrix.validate();

    const x = [_]f64{ 1.0, -2.0, 3.0, 0.0, 4.0, -1.0, 2.0 };
    var expected = [_]f64{0.0} ** row_count;
    for (0..row_count) |row| {
        for (0..col_count) |col| expected[row] += dense[row * col_count + col] * x[col];
    }
    var actual: [row_count]f64 = undefined;
    try matrix.multiply(&x, &actual);
    try std.testing.expectEqualSlices(f64, &expected, &actual);
}

test "sorted freeze avoids sorting and rejects unordered coordinates" {
    var builder = try MatrixBuilder.init(2, 2);
    defer builder.deinit(std.testing.allocator);
    try builder.append(std.testing.allocator, try RowId.init(0), try ColId.init(0), 1.0);
    try builder.append(std.testing.allocator, try RowId.init(1), try ColId.init(0), 2.0);
    try builder.append(std.testing.allocator, try RowId.init(0), try ColId.init(1), 3.0);
    var matrix = try builder.freezeSorted(std.testing.allocator, 0.0);
    defer matrix.deinit(std.testing.allocator);
    try matrix.validate();

    builder.clearRetainingCapacity();
    try builder.append(std.testing.allocator, try RowId.init(0), try ColId.init(1), 1.0);
    try builder.append(std.testing.allocator, try RowId.init(0), try ColId.init(0), 2.0);
    try std.testing.expectError(error.TripletsNotSorted, builder.freezeSorted(std.testing.allocator, 0.0));
}
