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
    /// Callers should `reserve` capacity before appending in a hot loop.
    pub fn append(self: *Self, allocator: std.mem.Allocator, row: RowId, col: ColId, value: f64) (std.mem.Allocator.Error || csc.MatrixError)!void {
        if (row.toUsize() >= self.num_rows or col.toUsize() >= self.num_cols)
            return error.IndexOutOfBounds;
        if (!std.math.isFinite(value)) return error.NonFiniteValue;

        // Inline capacity check: the separate reserve call is measurable overhead
        // in tight loops with pre-reserved capacity (builder_freeze_sorted hot path).
        if (self.entries.len >= self.entries.capacity)
            try self.entries.ensureUnusedCapacity(allocator, 1);
        self.entries.appendAssumeCapacity(.{ .row = row, .col = col, .value = value, .sequence = self.len() });
    }

    /// Unchecked append for benchmark and hot-path callers that pre-validate.
    ///
    /// Safety requirements (caller MUST guarantee all of these):
    /// 1. `reserve` or `ensureUnusedCapacity` has been called with sufficient space.
    /// 2. `row` and `col` are in-bounds for the builder dimensions.
    /// 3. `value` is finite (not NaN, not Inf).
    ///
    /// Violating any precondition is undefined behaviour. Prefer `append` for
    /// general-purpose code; use this only when the capacity check and validation
    /// have been hoisted out of a tight loop and verified by a test or assertion.
    pub fn appendPreReserved(self: *Self, row: RowId, col: ColId, value: f64) void {
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
        return self.freezeInternal(allocator, zero_tolerance, true, true);
    }

    /// Linear-time freeze for callers that already emit nondecreasing
    /// (column,row) coordinates. A linear order check is retained.
    pub fn freezeSorted(self: *Self, allocator: std.mem.Allocator, zero_tolerance: f64) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
        if (!self.isSorted()) return error.TripletsNotSorted;
        return self.freezeInternal(allocator, zero_tolerance, false, true);
    }

    /// Fastest construction path: caller guarantees nondecreasing coordinates.
    /// Duplicate coordinates may remain adjacent and are merged stably.
    pub fn freezeSortedAssumeValid(self: *Self, allocator: std.mem.Allocator, zero_tolerance: f64) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
        return self.freezeInternal(allocator, zero_tolerance, false, true);
    }

    /// Construction-latency variant for callers that do not need compact
    /// offsets immediately. It emits the same canonical CSC values and indices
    /// but skips the secondary HUInt offset array. Later kernels transparently
    /// use `col_starts`; choose the default API when the matrix will be used in
    /// many hot-path products.
    pub fn freezeSortedLeanAssumeValid(self: *Self, allocator: std.mem.Allocator, zero_tolerance: f64) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
        return self.freezeInternal(allocator, zero_tolerance, false, false);
    }

    fn freezeInternal(self: *Self, allocator: std.mem.Allocator, zero_tolerance: f64, sort_entries: bool, comptime include_compact_starts: bool) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
        if (!std.math.isFinite(zero_tolerance) or zero_tolerance < 0.0)
            return error.InvalidTolerance;
        const length = self.len();
        if (sort_entries) self.entries.sortUnstable(SortContext{ .fields = self.entries.slice() });
        var fields = self.entries.slice();
        const rows = fields.items(.row);
        const cols = fields.items(.col);
        const values = fields.items(.value);

        // Temporary column-start counts (freed after prefix-sum and output copy).
        const col_starts = try allocator.alloc(usize, self.num_cols + 1);
        defer allocator.free(col_starts);
        @memset(col_starts, 0);

        // Merge duplicates in-place while counting column occupancy.
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

            col_starts[col.toUsize() + 1] += 1;
            fields.set(write, .{ .row = row, .col = col, .value = sum, .sequence = write });
            write += 1;
        }
        self.entries.shrinkRetainingCapacity(write);
        fields = self.entries.slice();
        const compact_rows = fields.items(.row);
        const compact_values = fields.items(.value);

        // Prefix-sum raw counts into CSC offsets.
        const ncol = self.num_cols;
        var c: usize = 0;
        while (c < ncol) : (c += 1) col_starts[c + 1] += col_starts[c];

        if (comptime include_compact_starts) {
            return self.buildCompact(allocator, write, col_starts, rows, compact_rows, compact_values);
        } else {
            return self.buildLean(allocator, write, col_starts, compact_rows, compact_values);
        }
    }

    /// Build a canonical CSC matrix from merged triplets, allocating all output
    /// arrays in a single page-colored buffer with compact HUInt offsets.
    fn buildCompact(self: *Self, allocator: std.mem.Allocator, write: usize, col_starts: []usize, fields_rows: []RowId, compact_rows: []RowId, compact_values: []f64) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
        _ = fields_rows;
        const starts_bytes = std.math.mul(usize, self.num_cols + 1, @sizeOf(usize)) catch return error.DimensionTooLarge;
        const rows_bytes = std.math.mul(usize, write, @sizeOf(RowId)) catch return error.DimensionTooLarge;
        const values_bytes = std.math.mul(usize, write, @sizeOf(f64)) catch return error.DimensionTooLarge;
        const compact_starts_bytes = std.math.mul(usize, self.num_cols + 1, @sizeOf(foundation.HUInt)) catch return error.DimensionTooLarge;
        const layout = memory.computePageColoredLayout(4, .{ starts_bytes, compact_starts_bytes, rows_bytes, values_bytes }, .{ 0, 64, 128, 192 }) catch return error.DimensionTooLarge;
        const compact_starts_offset = layout.offsets[1];
        const rows_offset = layout.offsets[2];
        const values_offset = layout.offsets[3];
        const storage_len = layout.total;
        const storage = try allocator.alignedAlloc(u8, .@"64", storage_len);
        errdefer allocator.free(storage);
        const output_starts: [*]usize = @ptrCast(@alignCast(storage.ptr));
        const output_rows: [*]RowId = @ptrCast(@alignCast(storage.ptr + rows_offset));
        const output_values: [*]f64 = @ptrCast(@alignCast(storage.ptr + values_offset));
        @memcpy(output_starts[0 .. self.num_cols + 1], col_starts);
        @memcpy(output_rows[0..write], compact_rows);
        @memcpy(output_values[0..write], compact_values);
        const compact_starts: []foundation.HUInt = block: {
            const ptr: [*]foundation.HUInt = @ptrCast(@alignCast(storage.ptr + compact_starts_offset));
            const result = ptr[0 .. self.num_cols + 1];
            for (result, col_starts) |*dest, s| dest.* = @intCast(s);
            break :block result;
        };
        return .{
            .num_rows = self.num_rows,
            .num_cols = self.num_cols,
            .col_starts = output_starts[0 .. self.num_cols + 1],
            .row_indices = output_rows[0..write],
            .values = output_values[0..write],
            .storage = storage,
            .compact_col_starts = compact_starts,
        };
    }

    /// Lightweight output builder for callers that do not need compact offsets.
    /// No page coloring — the three output arrays are packed consecutively with
    /// natural alignment padding so each array starts on its preferred boundary.
    fn buildLean(self: *Self, allocator: std.mem.Allocator, write: usize, col_starts: []usize, compact_rows: []RowId, compact_values: []f64) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
        const starts_bytes = std.math.mul(usize, self.num_cols + 1, @sizeOf(usize)) catch return error.DimensionTooLarge;
        const rows_bytes = std.math.mul(usize, write, @sizeOf(RowId)) catch return error.DimensionTooLarge;
        const values_bytes = std.math.mul(usize, write, @sizeOf(f64)) catch return error.DimensionTooLarge;
        const layout = memory.computeLayout(3, .{ starts_bytes, rows_bytes, values_bytes }, .{ @alignOf(usize), @alignOf(RowId), @alignOf(f64) }) catch return error.DimensionTooLarge;
        const rows_offset = layout.offsets[1];
        const values_offset = layout.offsets[2];
        const storage_len = layout.total;
        const storage = try allocator.alignedAlloc(u8, .@"64", storage_len);
        errdefer allocator.free(storage);
        const output_starts = @as([*]usize, @ptrCast(@alignCast(storage.ptr)))[0 .. self.num_cols + 1];
        const row_indices = @as([*]RowId, @ptrCast(@alignCast(storage.ptr + rows_offset)))[0..write];
        const output_values = @as([*]f64, @ptrCast(@alignCast(storage.ptr + values_offset)))[0..write];
        @memcpy(output_starts, col_starts);
        @memcpy(row_indices, compact_rows);
        @memcpy(output_values, compact_values);
        return .{
            .num_rows = self.num_rows,
            .num_cols = self.num_cols,
            .col_starts = output_starts,
            .row_indices = row_indices,
            .values = output_values,
            .storage = storage,
            .compact_col_starts = null,
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

test "lean sorted freeze preserves canonical data without compact offsets" {
    var builder = try MatrixBuilder.init(3, 2);
    defer builder.deinit(std.testing.allocator);
    try builder.append(std.testing.allocator, try RowId.init(0), try ColId.init(0), 2.0);
    try builder.append(std.testing.allocator, try RowId.init(2), try ColId.init(1), 5.0);

    var matrix = try builder.freezeSortedLeanAssumeValid(std.testing.allocator, 0.0);
    defer matrix.deinit(std.testing.allocator);
    try matrix.validate();
    try std.testing.expect(matrix.compact_col_starts == null);
    try std.testing.expectEqualSlices(usize, &.{ 0, 1, 2 }, matrix.col_starts);
    try std.testing.expectEqualSlices(f64, &.{ 2.0, 5.0 }, matrix.values);
}

/// Experimental API: two-pass freeze from sorted arrays with duplicate merge.
///
/// This API may be corrected, renamed, or replaced after integration with the
/// LP/presolve layers. Do not treat its signature or allocation layout as stable.
///
/// Freeze directly from pre-sorted arrays, bypassing the builder's append
/// mechanism.  Useful for callers that already have sorted (col, row) data
/// and want to skip the MultiArrayList indirection.
///
/// The arrays must be sorted by (col, row). Adjacent duplicates are merged
/// stably. The output is a canonical CSC matrix.
/// Prefer `freezeCanonicalIntoAssumeValid` + `CscBuildBuffers` for hot paths
/// where the caller can own the buffer lifecycle.
pub fn freezeFromSortedArraysAssumeValid(allocator: std.mem.Allocator, num_rows: usize, num_cols: usize, sorted_rows: []const RowId, sorted_cols: []const ColId, sorted_values: []const f64, zero_tolerance: f64, comptime include_compact_starts: bool) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    if (!std.math.isFinite(zero_tolerance) or zero_tolerance < 0.0)
        return error.InvalidTolerance;
    try csc.validateDimensions(num_rows, num_cols);
    const length = sorted_rows.len;
    if (sorted_cols.len != length or sorted_values.len != length)
        return error.InconsistentStorage;

    // Allocate temporary column-start counts.
    const col_starts = try allocator.alloc(usize, num_cols + 1);
    defer allocator.free(col_starts);
    @memset(col_starts, 0);

    // Pass 1: merge + count columns (no data storage — only col_starts)
    var read1: usize = 0;
    while (read1 < length) {
        const row = sorted_rows[read1];
        const col = sorted_cols[read1];
        var sum = sorted_values[read1];
        read1 += 1;
        while (read1 < length and sameCoordinate(sorted_rows[read1], sorted_cols[read1], row, col)) : (read1 += 1) {
            sum += sorted_values[read1];
        }
        if (!std.math.isFinite(sum)) return error.NonFiniteValue;
        if (@abs(sum) <= zero_tolerance) continue;
        col_starts[col.toUsize() + 1] += 1;
    }

    // Prefix-sum raw counts into CSC offsets.
    var c: usize = 0;
    while (c < num_cols) : (c += 1) col_starts[c + 1] += col_starts[c];
    const final_nnz = col_starts[num_cols];

    // Allocate output storage now that final_nnz is known.
    // Use a cursor array for pass 2 — borrow the first num_cols elements
    // of col_starts as cursor (col_starts is prefix-sum, cursor is copy).
    // After allocation, col_starts is memcpy'd to output and then reused as cursor.
    const cursor = try allocator.alloc(foundation.HUInt, num_cols);
    defer allocator.free(cursor);

    if (comptime include_compact_starts) {
        const starts_bytes = std.math.mul(usize, num_cols + 1, @sizeOf(usize)) catch return error.DimensionTooLarge;
        const rows_bytes = std.math.mul(usize, final_nnz, @sizeOf(RowId)) catch return error.DimensionTooLarge;
        const values_bytes = std.math.mul(usize, final_nnz, @sizeOf(f64)) catch return error.DimensionTooLarge;
        const compact_starts_bytes = std.math.mul(usize, num_cols + 1, @sizeOf(foundation.HUInt)) catch return error.DimensionTooLarge;
        const layout = memory.computePageColoredLayout(4, .{ starts_bytes, compact_starts_bytes, rows_bytes, values_bytes }, .{ 0, 64, 128, 192 }) catch return error.DimensionTooLarge;
        const storage_len = layout.total;
        const storage = try allocator.alignedAlloc(u8, .@"64", storage_len);
        errdefer allocator.free(storage);
        const output_starts: [*]usize = @ptrCast(@alignCast(storage.ptr));
        const output_rows: [*]RowId = @ptrCast(@alignCast(storage.ptr + layout.offsets[2]));
        const output_values: [*]f64 = @ptrCast(@alignCast(storage.ptr + layout.offsets[3]));
        @memcpy(output_starts[0 .. num_cols + 1], col_starts);

        // Pass 2: merge again, write directly to output.
        for (cursor[0..num_cols], col_starts[0..num_cols]) |*dest, s| dest.* = @intCast(s);
        var read2: usize = 0;
        while (read2 < length) {
            const row = sorted_rows[read2];
            const col = sorted_cols[read2];
            var sum = sorted_values[read2];
            read2 += 1;
            while (read2 < length and sameCoordinate(sorted_rows[read2], sorted_cols[read2], row, col)) : (read2 += 1) {
                sum += sorted_values[read2];
            }
            if (@abs(sum) <= zero_tolerance) continue;
            const dest: usize = @intCast(cursor[col.toUsize()]);
            cursor[col.toUsize()] += 1;
            output_rows[dest] = row;
            output_values[dest] = sum;
        }

        const compact_starts: []foundation.HUInt = block: {
            const ptr: [*]foundation.HUInt = @ptrCast(@alignCast(storage.ptr + layout.offsets[1]));
            const result = ptr[0 .. num_cols + 1];
            for (result, col_starts) |*dest, s| dest.* = @intCast(s);
            break :block result;
        };
        return .{ .num_rows = num_rows, .num_cols = num_cols, .col_starts = output_starts[0 .. num_cols + 1], .row_indices = output_rows[0..final_nnz], .values = output_values[0..final_nnz], .storage = storage, .compact_col_starts = compact_starts };
    } else {
        const starts_bytes = std.math.mul(usize, num_cols + 1, @sizeOf(usize)) catch return error.DimensionTooLarge;
        const rows_bytes = std.math.mul(usize, final_nnz, @sizeOf(RowId)) catch return error.DimensionTooLarge;
        const values_bytes = std.math.mul(usize, final_nnz, @sizeOf(f64)) catch return error.DimensionTooLarge;
        const layout = memory.computeLayout(3, .{ starts_bytes, rows_bytes, values_bytes }, .{ @alignOf(usize), @alignOf(RowId), @alignOf(f64) }) catch return error.DimensionTooLarge;
        const storage = try allocator.alignedAlloc(u8, .@"64", layout.total);
        errdefer allocator.free(storage);
        const output_starts = @as([*]usize, @ptrCast(@alignCast(storage.ptr)))[0 .. num_cols + 1];
        const output_rows = @as([*]RowId, @ptrCast(@alignCast(storage.ptr + layout.offsets[1])))[0..final_nnz];
        const output_values = @as([*]f64, @ptrCast(@alignCast(storage.ptr + layout.offsets[2])))[0..final_nnz];
        @memcpy(output_starts, col_starts);

        // Pass 2: merge again, write directly to output.
        for (cursor[0..num_cols], col_starts[0..num_cols]) |*dest, s| dest.* = @intCast(s);
        var read2: usize = 0;
        while (read2 < length) {
            const row = sorted_rows[read2];
            const col = sorted_cols[read2];
            var sum = sorted_values[read2];
            read2 += 1;
            while (read2 < length and sameCoordinate(sorted_rows[read2], sorted_cols[read2], row, col)) : (read2 += 1) {
                sum += sorted_values[read2];
            }
            if (@abs(sum) <= zero_tolerance) continue;
            const dest: usize = @intCast(cursor[col.toUsize()]);
            cursor[col.toUsize()] += 1;
            output_rows[dest] = row;
            output_values[dest] = sum;
        }

        return .{ .num_rows = num_rows, .num_cols = num_cols, .col_starts = output_starts, .row_indices = output_rows, .values = output_values, .storage = storage, .compact_col_starts = null };
    }
}

/// Experimental API: owning freeze for canonical (no-duplicate) arrays.
///
/// This API may be corrected, renamed, or replaced by the reusable-buffer path.
/// Do not treat its signature or allocation layout as stable.
///
/// Fast-path freeze for canonical arrays.  Caller guarantees:
/// - Sorted by (col, row)
/// - No duplicate coordinates
/// - No explicit zeros (abs(value) > zero_tolerance)
/// - All values are finite
///
/// Because no merging is needed, this is a single-pass count + prefix +
/// memcpy with zero branch-miss pressure from merge logic.
/// Prefer `freezeCanonicalIntoAssumeValid` + `CscBuildBuffers` for hot paths.
pub fn freezeFromCanonicalArraysAssumeValid(allocator: std.mem.Allocator, num_rows: usize, num_cols: usize, sorted_rows: []const RowId, sorted_cols: []const ColId, sorted_values: []const f64, zero_tolerance: f64, comptime include_compact_starts: bool) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    _ = zero_tolerance; // caller guarantees canonical — no zeros to filter
    try csc.validateDimensions(num_rows, num_cols);
    const nnz = sorted_rows.len;
    if (sorted_cols.len != nnz or sorted_values.len != nnz)
        return error.InconsistentStorage;
    if (nnz == 0) {
        const starts = try allocator.alloc(usize, num_cols + 1);
        @memset(starts, 0);
        return .{ .num_rows = num_rows, .num_cols = num_cols, .col_starts = starts, .row_indices = &.{}, .values = &.{}, .storage = null, .compact_col_starts = null };
    }

    // Pass 1: count entries per column (no merge needed — input is canonical)
    const col_starts = try allocator.alloc(usize, num_cols + 1);
    defer allocator.free(col_starts);
    @memset(col_starts, 0);
    for (sorted_cols) |col| col_starts[col.toUsize() + 1] += 1;

    // Prefix-sum counts into CSC offsets
    var c: usize = 0;
    while (c < num_cols) : (c += 1) col_starts[c + 1] += col_starts[c];

    if (comptime include_compact_starts) {
        const starts_bytes = std.math.mul(usize, num_cols + 1, @sizeOf(usize)) catch return error.DimensionTooLarge;
        const rows_bytes = std.math.mul(usize, nnz, @sizeOf(RowId)) catch return error.DimensionTooLarge;
        const values_bytes = std.math.mul(usize, nnz, @sizeOf(f64)) catch return error.DimensionTooLarge;
        const compact_starts_bytes = std.math.mul(usize, num_cols + 1, @sizeOf(foundation.HUInt)) catch return error.DimensionTooLarge;
        const layout = memory.computePageColoredLayout(4, .{ starts_bytes, compact_starts_bytes, rows_bytes, values_bytes }, .{ 0, 64, 128, 192 }) catch return error.DimensionTooLarge;
        const storage_len = layout.total;
        const storage = try allocator.alignedAlloc(u8, .@"64", storage_len);
        errdefer allocator.free(storage);
        const output_starts: [*]usize = @ptrCast(@alignCast(storage.ptr));
        const output_rows: [*]RowId = @ptrCast(@alignCast(storage.ptr + layout.offsets[2]));
        const output_values: [*]f64 = @ptrCast(@alignCast(storage.ptr + layout.offsets[3]));
        @memcpy(output_starts[0 .. num_cols + 1], col_starts);
        @memcpy(output_rows[0..nnz], sorted_rows);
        @memcpy(output_values[0..nnz], sorted_values);
        const compact_starts: []foundation.HUInt = block: {
            const ptr: [*]foundation.HUInt = @ptrCast(@alignCast(storage.ptr + layout.offsets[1]));
            const result = ptr[0 .. num_cols + 1];
            for (result, col_starts) |*dest, s| dest.* = @intCast(s);
            break :block result;
        };
        return .{ .num_rows = num_rows, .num_cols = num_cols, .col_starts = output_starts[0 .. num_cols + 1], .row_indices = output_rows[0..nnz], .values = output_values[0..nnz], .storage = storage, .compact_col_starts = compact_starts };
    } else {
        const starts_bytes = std.math.mul(usize, num_cols + 1, @sizeOf(usize)) catch return error.DimensionTooLarge;
        const rows_bytes = std.math.mul(usize, nnz, @sizeOf(RowId)) catch return error.DimensionTooLarge;
        const values_bytes = std.math.mul(usize, nnz, @sizeOf(f64)) catch return error.DimensionTooLarge;
        const layout = memory.computeLayout(3, .{ starts_bytes, rows_bytes, values_bytes }, .{ @alignOf(usize), @alignOf(RowId), @alignOf(f64) }) catch return error.DimensionTooLarge;
        const storage = try allocator.alignedAlloc(u8, .@"64", layout.total);
        errdefer allocator.free(storage);
        const output_starts = @as([*]usize, @ptrCast(@alignCast(storage.ptr)))[0 .. num_cols + 1];
        const output_rows = @as([*]RowId, @ptrCast(@alignCast(storage.ptr + layout.offsets[1])))[0..nnz];
        const output_values = @as([*]f64, @ptrCast(@alignCast(storage.ptr + layout.offsets[2])))[0..nnz];
        @memcpy(output_starts, col_starts);
        @memcpy(output_rows, sorted_rows);
        @memcpy(output_values, sorted_values);
        return .{ .num_rows = num_rows, .num_cols = num_cols, .col_starts = output_starts, .row_indices = output_rows, .values = output_values, .storage = storage, .compact_col_starts = null };
    }
}

test "freezeFromCanonical handles tridiagonal with no duplicates" {
    const rows = [_]RowId{ try RowId.init(0), try RowId.init(1), try RowId.init(2) };
    const cols = [_]ColId{ try ColId.init(0), try ColId.init(1), try ColId.init(2) };
    const values = [_]f64{ 1.0, 2.0, 3.0 };
    var matrix = try freezeFromCanonicalArraysAssumeValid(std.testing.allocator, 3, 3, &rows, &cols, &values, 0.0, false);
    defer matrix.deinit(std.testing.allocator);
    try matrix.validate();
    try std.testing.expectEqual(@as(usize, 3), matrix.nnz());
    try std.testing.expectEqualSlices(usize, &.{ 0, 1, 2, 3 }, matrix.col_starts);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 2.0, 3.0 }, matrix.values);
}

/// Stable API: caller-owned reusable buffers for CSC construction.
///
/// The caller allocates and manages the lifetime of all three slices.
/// Write functions accept a `*CscBuildBuffers` to update `col_starts` in
/// place and return a `CscView` whose `row_indices` and `values` lengths
/// reflect the actual nnz written.
pub const CscBuildBuffers = struct {
    col_starts: []usize,
    row_indices: []RowId,
    values: []f64,

    pub fn initCapacity(allocator: std.mem.Allocator, num_cols: usize, nnz_capacity: usize) (std.mem.Allocator.Error || csc.MatrixError)!CscBuildBuffers {
        if (num_cols == std.math.maxInt(usize)) return error.DimensionTooLarge;
        const starts = try allocator.alloc(usize, num_cols + 1);
        errdefer allocator.free(starts);
        const indices = try allocator.alloc(RowId, nnz_capacity);
        errdefer allocator.free(indices);
        const vals = try allocator.alloc(f64, nnz_capacity);
        return .{ .col_starts = starts, .row_indices = indices, .values = vals };
    }

    pub fn deinit(self: *CscBuildBuffers, allocator: std.mem.Allocator) void {
        allocator.free(self.col_starts);
        allocator.free(self.row_indices);
        allocator.free(self.values);
        self.* = undefined;
    }
};

/// Stable API: caller-owned reusable CSC freeze.
///
/// Write canonical (sorted, unique, nonzero, finite) coordinate triples into
/// caller-owned `CscBuildBuffers`.  Returns a borrowed `CscView` pointing into
/// the buffers; the caller retains ownership and must keep the buffers alive
/// while the view is in use.
///
/// Capacity: `buffers.row_indices.len` and `buffers.values.len` must be >= the
/// final nnz, and `buffers.col_starts.len` must be >= num_cols + 1.
/// Returns `error.BufferTooSmall` if any buffer is undersized.
///
/// Only `buffers.col_starts[0..num_cols+1]` is written; the remainder of the
/// buffer (if any) is untouched.  Callers may pre-touch pages by zero-filling
/// the full capacity before the first call, but the function itself does not
/// force any pre-touch policy.
///
/// Contract: the input arrays (`rows`, `cols`, `vals`) must NOT overlap with
/// `buffers.row_indices` or `buffers.values`.  Overlap with `buffers.col_starts`
/// is permitted only if the caller guarantees that `cols` is disjoint from the
/// written region of `col_starts` — otherwise the count phase may be corrupted
/// by the subsequent prefix-sum pass.
pub fn freezeCanonicalIntoAssumeValid(buffers: *CscBuildBuffers, num_rows: usize, num_cols: usize, rows: []const RowId, cols: []const ColId, vals: []const f64) csc.MatrixError!csc.CscView {
    const nnz = rows.len;
    if (cols.len != nnz or vals.len != nnz) return error.InconsistentStorage;
    if (num_cols == std.math.maxInt(usize)) return error.DimensionTooLarge;
    if (buffers.col_starts.len < num_cols + 1) return error.BufferTooSmall;
    if (buffers.row_indices.len < nnz) return error.BufferTooSmall;
    if (buffers.values.len < nnz) return error.BufferTooSmall;

    // Count + prefix-sum into caller's col_starts — only the first num_cols+1 slots
    @memset(buffers.col_starts[0 .. num_cols + 1], 0);
    for (cols) |col| buffers.col_starts[col.toUsize() + 1] += 1;
    var c: usize = 0;
    while (c < num_cols) : (c += 1) buffers.col_starts[c + 1] += buffers.col_starts[c];

    @memcpy(buffers.row_indices[0..nnz], rows);
    @memcpy(buffers.values[0..nnz], vals);

    return .{ .num_rows = num_rows, .num_cols = num_cols, .col_starts = buffers.col_starts[0 .. num_cols + 1], .row_indices = buffers.row_indices[0..nnz], .values = buffers.values[0..nnz] };
}

test "freezeCanonicalIntoAssumeValid writes correct view" {
    var bufs = try CscBuildBuffers.initCapacity(std.testing.allocator, 3, 10);
    defer bufs.deinit(std.testing.allocator);
    const rows = [_]RowId{ try RowId.init(0), try RowId.init(1), try RowId.init(2) };
    const cols = [_]ColId{ try ColId.init(0), try ColId.init(1), try ColId.init(2) };
    const vals = [_]f64{ 1.0, 2.0, 3.0 };
    const view = try freezeCanonicalIntoAssumeValid(&bufs, 3, 3, &rows, &cols, &vals);
    try view.validate();
    try std.testing.expectEqual(@as(usize, 3), view.nnz());
    try std.testing.expectEqualSlices(usize, &.{ 0, 1, 2, 3 }, view.col_starts);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 2.0, 3.0 }, view.values);
}

test "freezeCanonicalIntoAssumeValid rejects undersized buffer" {
    var bufs = try CscBuildBuffers.initCapacity(std.testing.allocator, 3, 1);
    defer bufs.deinit(std.testing.allocator);
    const rows = [_]RowId{ try RowId.init(0), try RowId.init(1) };
    const cols = [_]ColId{ try ColId.init(0), try ColId.init(1) };
    const vals = [_]f64{ 1.0, 2.0 };
    try std.testing.expectError(error.BufferTooSmall, freezeCanonicalIntoAssumeValid(&bufs, 3, 2, &rows, &cols, &vals));
}

test "freezeCanonicalIntoAssumeValid large buffer small matrix" {
    // Buffer allocated for 8 columns and 100 nnz, but only 3 columns + 2 entries used.
    var bufs = try CscBuildBuffers.initCapacity(std.testing.allocator, 8, 100);
    defer bufs.deinit(std.testing.allocator);
    // Pre-fill with sentinel to detect overwrites beyond valid region
    @memset(bufs.col_starts, 0xAA);
    @memset(std.mem.sliceAsBytes(bufs.row_indices), 0xBB);
    @memset(bufs.values, -1.0);
    const rows = [_]RowId{ try RowId.init(0), try RowId.init(2) };
    const cols = [_]ColId{ try ColId.init(0), try ColId.init(0) };
    const vals = [_]f64{ 1.0, 2.0 };
    const view = try freezeCanonicalIntoAssumeValid(&bufs, 3, 3, &rows, &cols, &vals);
    try view.validate();
    try std.testing.expectEqual(@as(usize, 2), view.nnz());
    try std.testing.expectEqualSlices(usize, &.{ 0, 2, 2, 2 }, view.col_starts);
    // Verify untouched starts region beyond num_cols+1 is unchanged
    for (4..bufs.col_starts.len) |k| try std.testing.expectEqual(@as(u8, 0xAA), @as(u8, @truncate(bufs.col_starts[k])));
    // Verify untouched row_indices/values beyond nnz are unchanged
    for (2..bufs.row_indices.len) |k| try std.testing.expect(bufs.row_indices[k].toUsize() != 0);
    for (2..bufs.values.len) |k| try std.testing.expectEqual(@as(f64, -1.0), bufs.values[k]);
}

test "freezeCanonicalIntoAssumeValid matches owning freeze" {
    var prng = std.Random.DefaultPrng.init(0xC5C5_2026);
    const random = prng.random();
    const num_rows: usize = 15;
    const num_cols: usize = 12;
    const max_nnz: usize = 200;

    // Use fixed arrays (worst-case size)
    var rows_buf: [200]RowId = undefined;
    var cols_buf: [200]ColId = undefined;
    var vals_buf: [200]f64 = undefined;
    var n: usize = 0;

    // Track for duplicate removal via simple linear scan
    for (0..max_nnz) |_| {
        const row = random.intRangeLessThan(usize, 0, num_rows);
        const col = random.intRangeLessThan(usize, 0, num_cols);
        const value: f64 = @floatFromInt(random.intRangeAtMost(i16, -10, 10));
        if (value == 0.0) continue;
        // Check for duplicates in already-added entries
        var dup = false;
        for (0..n) |k| {
            if (rows_buf[k].toUsize() == row and cols_buf[k].toUsize() == col) {
                dup = true;
                break;
            }
        }
        if (dup) continue;
        if (n >= 200) break;
        rows_buf[n] = try RowId.fromUsize(row);
        cols_buf[n] = try ColId.fromUsize(col);
        vals_buf[n] = value;
        n += 1;
    }

    // Sort by (col, row) — bubble sort for simplicity
    for (0..n) |i| {
        for (i + 1..n) |j| {
            if (cols_buf[i].toUsize() > cols_buf[j].toUsize() or
                (cols_buf[i].toUsize() == cols_buf[j].toUsize() and
                    rows_buf[i].toUsize() > rows_buf[j].toUsize()))
            {
                const tmp_r = rows_buf[i];
                rows_buf[i] = rows_buf[j];
                rows_buf[j] = tmp_r;
                const tmp_c = cols_buf[i];
                cols_buf[i] = cols_buf[j];
                cols_buf[j] = tmp_c;
                const tmp_v = vals_buf[i];
                vals_buf[i] = vals_buf[j];
                vals_buf[j] = tmp_v;
            }
        }
    }

    const rows = rows_buf[0..n];
    const cols = cols_buf[0..n];
    const vals = vals_buf[0..n];

    // Build via owning API
    var owning = try freezeFromCanonicalArraysAssumeValid(std.testing.allocator, num_rows, num_cols, rows, cols, vals, 0.0, false);
    defer owning.deinit(std.testing.allocator);

    // Build via reusable API
    var bufs = try CscBuildBuffers.initCapacity(std.testing.allocator, num_cols, n);
    defer bufs.deinit(std.testing.allocator);
    const view = try freezeCanonicalIntoAssumeValid(&bufs, num_rows, num_cols, rows, cols, vals);

    // Verify identical
    try std.testing.expectEqual(owning.nnz(), view.nnz());
    try std.testing.expectEqualSlices(usize, owning.col_starts, view.col_starts);
    try std.testing.expectEqualSlices(RowId, owning.row_indices, view.row_indices);
    try std.testing.expectEqualSlices(f64, owning.values, view.values);
}

test "freezeFromSortedArrays merges duplicates" {
    const rows = [_]RowId{ try RowId.init(0), try RowId.init(0), try RowId.init(1) };
    const cols = [_]ColId{ try ColId.init(0), try ColId.init(0), try ColId.init(0) };
    const values = [_]f64{ 2.0, -2.0, 3.0 };
    var matrix = try freezeFromSortedArraysAssumeValid(std.testing.allocator, 3, 2, &rows, &cols, &values, 0.0, false);
    defer matrix.deinit(std.testing.allocator);
    try matrix.validate();
    try std.testing.expectEqual(@as(usize, 1), matrix.nnz());
    try std.testing.expectEqual(@as(f64, 3.0), matrix.values[0]);
}

test "freezeFromSortedArrays handles odd nnz" {
    const rows = [_]RowId{ try RowId.init(0), try RowId.init(1), try RowId.init(2) };
    const cols = [_]ColId{ try ColId.init(0), try ColId.init(1), try ColId.init(1) };
    const values = [_]f64{ 1.0, 2.0, 3.0 };
    var matrix = try freezeFromSortedArraysAssumeValid(std.testing.allocator, 3, 2, &rows, &cols, &values, 0.0, false);
    defer matrix.deinit(std.testing.allocator);
    try matrix.validate();
    try std.testing.expectEqual(@as(usize, 3), matrix.nnz());
    try std.testing.expectEqualSlices(usize, &.{ 0, 1, 3 }, matrix.col_starts);
}

test "freezeFromSortedArrays handles empty matrix (zero nnz)" {
    const rows = [_]RowId{};
    const cols = [_]ColId{};
    const values = [_]f64{};
    var matrix = try freezeFromSortedArraysAssumeValid(std.testing.allocator, 3, 2, &rows, &cols, &values, 0.0, false);
    defer matrix.deinit(std.testing.allocator);
    try matrix.validate();
    try std.testing.expectEqual(@as(usize, 0), matrix.nnz());
    try std.testing.expectEqualSlices(usize, &.{ 0, 0, 0 }, matrix.col_starts);
}

test "freezeFromSortedArrays handles empty rows and columns" {
    const rows = [_]RowId{ try RowId.init(0), try RowId.init(2) };
    const cols = [_]ColId{ try ColId.init(0), try ColId.init(0) };
    const values = [_]f64{ 1.0, 2.0 };
    var matrix = try freezeFromSortedArraysAssumeValid(std.testing.allocator, 4, 3, &rows, &cols, &values, 0.0, false);
    defer matrix.deinit(std.testing.allocator);
    try matrix.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 2, 2, 2 }, matrix.col_starts);
}
