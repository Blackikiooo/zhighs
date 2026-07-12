//! Explicit CSC transpose construction.
//!
//! The conversion uses counting plus prefix sums. It is O(rows + nnz), and the
//! output is canonical without sorting because source columns are visited in
//! increasing order.

const std = @import("std");
const foundation = @import("foundation");
const csc = @import("csc.zig");
const memory = @import("memory.zig");

pub const TransposeBuffers = struct {
    storage: []align(64) u8,
    starts: []usize,
    compact_starts: []foundation.HUInt,
    rows: []foundation.RowId,
    values: []f64,
    cursor: []foundation.HUInt,

    pub fn init(allocator: std.mem.Allocator, num_cols: usize, nnz: usize) (std.mem.Allocator.Error || csc.MatrixError)!TransposeBuffers {
        if (num_cols == std.math.maxInt(usize)) return error.DimensionTooLarge;
        const starts_bytes = std.math.mul(usize, num_cols + 1, @sizeOf(usize)) catch return error.DimensionTooLarge;
        const compact_bytes = std.math.mul(usize, num_cols + 1, @sizeOf(foundation.HUInt)) catch return error.DimensionTooLarge;
        const rows_bytes = std.math.mul(usize, nnz, @sizeOf(foundation.RowId)) catch return error.DimensionTooLarge;
        const values_bytes = std.math.mul(usize, nnz, @sizeOf(f64)) catch return error.DimensionTooLarge;
        const cursor_bytes = std.math.mul(usize, num_cols, @sizeOf(foundation.HUInt)) catch return error.DimensionTooLarge;
        const layout = memory.computePageColoredLayout(5, .{ starts_bytes, compact_bytes, rows_bytes, values_bytes, cursor_bytes }, .{ 0, 64, 128, 192, 256 }) catch return error.DimensionTooLarge;
        const storage = try allocator.alignedAlloc(u8, .@"64", layout.total);
        const starts_ptr: [*]usize = @ptrCast(@alignCast(storage.ptr + layout.offsets[0]));
        const compact_ptr: [*]foundation.HUInt = @ptrCast(@alignCast(storage.ptr + layout.offsets[1]));
        const rows_ptr: [*]foundation.RowId = @ptrCast(@alignCast(storage.ptr + layout.offsets[2]));
        const values_ptr: [*]f64 = @ptrCast(@alignCast(storage.ptr + layout.offsets[3]));
        const cursor_ptr: [*]foundation.HUInt = @ptrCast(@alignCast(storage.ptr + layout.offsets[4]));
        return .{
            .storage = storage,
            .starts = starts_ptr[0 .. num_cols + 1],
            .compact_starts = compact_ptr[0 .. num_cols + 1],
            .rows = rows_ptr[0..nnz],
            .values = values_ptr[0..nnz],
            .cursor = cursor_ptr[0..num_cols],
        };
    }

    pub fn deinit(self: *TransposeBuffers, allocator: std.mem.Allocator) void {
        allocator.free(self.storage);
        self.* = undefined;
    }
};

/// Checked transpose for matrices entering from an untrusted boundary.
pub fn transpose(allocator: std.mem.Allocator, matrix: csc.CscMatrix) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    try matrix.validate();
    return transposeAssumeValid(allocator, matrix);
}

/// Transposes a canonical CSC matrix without repeating structural validation.
/// The result carries compact HUInt offsets for downstream CSR/transpose speed.
pub fn transposeAssumeValid(allocator: std.mem.Allocator, matrix: csc.CscMatrix) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    if (matrix.num_rows == std.math.maxInt(usize)) return error.DimensionTooLarge;
    var buffers = try TransposeBuffers.init(allocator, matrix.num_rows, matrix.nnz());
    errdefer buffers.deinit(allocator);
    try transposeIntoAssumeValid(matrix, buffers.starts, buffers.rows, buffers.values, buffers.cursor);
    for (buffers.compact_starts, buffers.starts) |*destination, start| destination.* = @intCast(start);

    return .{
        .num_rows = matrix.num_cols,
        .num_cols = matrix.num_rows,
        .col_starts = buffers.starts,
        .row_indices = buffers.rows,
        .values = buffers.values,
        .storage = buffers.storage,
        .compact_col_starts = buffers.compact_starts,
    };
}

/// Fast transpose without compact HUInt offsets.
/// Uses byte-clear instead of volatile SIMD to keep starts zeroing fast.
/// Cursor scratch space is allocated and freed within the call so the returned
/// matrix does not retain dead scratch bytes.
pub fn transposeLeanAssumeValid(allocator: std.mem.Allocator, matrix: csc.CscMatrix) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    if (matrix.num_rows == std.math.maxInt(usize)) return error.DimensionTooLarge;
    const nnz = matrix.nnz();
    const new_num_cols = matrix.num_rows;
    const starts_bytes = std.math.mul(usize, new_num_cols + 1, @sizeOf(usize)) catch return error.DimensionTooLarge;
    const rows_bytes = std.math.mul(usize, nnz, @sizeOf(foundation.RowId)) catch return error.DimensionTooLarge;
    const values_bytes = std.math.mul(usize, nnz, @sizeOf(f64)) catch return error.DimensionTooLarge;
    const layout = memory.computeLayout(3, .{ starts_bytes, rows_bytes, values_bytes }, .{ @alignOf(usize), @alignOf(foundation.RowId), @alignOf(f64) }) catch return error.DimensionTooLarge;
    const rows_offset = layout.offsets[1];
    const values_offset = layout.offsets[2];
    const storage_len = layout.total;
    const storage = try allocator.alignedAlloc(u8, .@"64", storage_len);
    errdefer allocator.free(storage);
    const starts = @as([*]usize, @ptrCast(@alignCast(storage.ptr)))[0 .. new_num_cols + 1];
    const rows = @as([*]foundation.RowId, @ptrCast(@alignCast(storage.ptr + rows_offset)))[0..nnz];
    const out_values = @as([*]f64, @ptrCast(@alignCast(storage.ptr + values_offset)))[0..nnz];

    // Cursor is scratch — allocate separately so the returned matrix doesn't
    // permanently retain dead bytes.
    const cursor = try allocator.alloc(foundation.HUInt, new_num_cols);
    defer allocator.free(cursor);

    @memset(std.mem.sliceAsBytes(starts), 0);
    for (matrix.row_indices) |row_id| starts[row_id.toUsize() + 1] += 1;
    for (0..new_num_cols) |row| starts[row + 1] += starts[row];
    for (cursor, starts[0..new_num_cols]) |*dest, s| dest.* = @intCast(s);
    if (matrix.compact_col_starts) |source_starts| {
        fillTransposeEntries(foundation.HUInt, matrix.num_cols, source_starts, matrix.row_indices, matrix.values, cursor, rows, out_values);
    } else {
        fillTransposeEntries(usize, matrix.num_cols, matrix.col_starts, matrix.row_indices, matrix.values, cursor, rows, out_values);
    }
    return .{
        .num_rows = matrix.num_cols,
        .num_cols = matrix.num_rows,
        .col_starts = starts,
        .row_indices = rows,
        .values = out_values,
        .storage = storage,
        .compact_col_starts = null,
    };
}

/// Owning transpose using HUInt (4-byte) internal starts for reduced memory
/// traffic.  Otherwise identical to `transposeLeanAssumeValid`.
pub fn transposeLeanAssumeValidCompact(allocator: std.mem.Allocator, matrix: csc.CscMatrix) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    if (matrix.num_rows == std.math.maxInt(usize)) return error.DimensionTooLarge;
    const nnz = matrix.nnz();
    const new_num_cols = matrix.num_rows;
    const starts_compact = try allocator.alloc(foundation.HUInt, new_num_cols + 1);
    defer allocator.free(starts_compact);

    const starts_bytes = std.math.mul(usize, new_num_cols + 1, @sizeOf(usize)) catch return error.DimensionTooLarge;
    const rows_bytes = std.math.mul(usize, nnz, @sizeOf(foundation.RowId)) catch return error.DimensionTooLarge;
    const values_bytes = std.math.mul(usize, nnz, @sizeOf(f64)) catch return error.DimensionTooLarge;
    const layout = memory.computeLayout(3, .{ starts_bytes, rows_bytes, values_bytes }, .{ @alignOf(usize), @alignOf(foundation.RowId), @alignOf(f64) }) catch return error.DimensionTooLarge;
    const storage = try allocator.alignedAlloc(u8, .@"64", layout.total);
    errdefer allocator.free(storage);
    const starts = @as([*]usize, @ptrCast(@alignCast(storage.ptr)))[0 .. new_num_cols + 1];
    const rows = @as([*]foundation.RowId, @ptrCast(@alignCast(storage.ptr + layout.offsets[1])))[0..nnz];
    const out_values = @as([*]f64, @ptrCast(@alignCast(storage.ptr + layout.offsets[2])))[0..nnz];

    const cursor = try allocator.alloc(foundation.HUInt, new_num_cols);
    defer allocator.free(cursor);

    // Histogram + prefix on HUInt (4-byte, matching C++)
    memory.clearTyped(foundation.HUInt, starts_compact);
    for (matrix.row_indices) |row_id| starts_compact[row_id.toUsize() + 1] += 1;
    for (0..new_num_cols) |row| starts_compact[row + 1] += starts_compact[row];

    // Convert HUInt → usize for output, copy cursor
    for (starts[0..new_num_cols], starts_compact[0..new_num_cols]) |*dest, s| dest.* = @intCast(s);
    starts[new_num_cols] = @intCast(starts_compact[new_num_cols]);
    for (cursor, starts_compact[0..new_num_cols]) |*dest, s| dest.* = s;

    if (matrix.compact_col_starts) |source_starts| {
        fillTransposeEntries(foundation.HUInt, matrix.num_cols, source_starts, matrix.row_indices, matrix.values, cursor, rows, out_values);
    } else {
        fillTransposeEntries(usize, matrix.num_cols, matrix.col_starts, matrix.row_indices, matrix.values, cursor, rows, out_values);
    }
    return .{ .num_rows = matrix.num_cols, .num_cols = matrix.num_rows, .col_starts = starts, .row_indices = rows, .values = out_values, .storage = storage, .compact_col_starts = null };
}

pub fn transposeInto(matrix: csc.CscMatrix, starts: []usize, rows: []foundation.RowId, values: []f64, cursor_scratch: []foundation.HUInt) csc.MatrixError!void {
    try matrix.validate();
    if (matrix.nnz() > std.math.maxInt(foundation.HUInt)) return error.DimensionTooLarge;
    return transposeIntoAssumeValid(matrix, starts, rows, values, cursor_scratch);
}

/// Allocation-free explicit transpose into exact-size caller buffers.
pub fn transposeIntoAssumeValid(matrix: csc.CscMatrix, starts: []usize, rows: []foundation.RowId, values: []f64, cursor_scratch: []foundation.HUInt) csc.MatrixError!void {
    if (starts.len != matrix.num_rows + 1 or rows.len != matrix.nnz() or
        values.len != matrix.nnz() or cursor_scratch.len < matrix.num_rows)
        return error.DimensionMismatch;
    memory.clearUsize(starts);
    for (matrix.row_indices) |row_id| starts[row_id.toUsize() + 1] += 1;
    for (0..matrix.num_rows) |row| starts[row + 1] += starts[row];
    const next = cursor_scratch[0..matrix.num_rows];
    for (next, starts[0..matrix.num_rows]) |*destination, start| destination.* = @intCast(start);
    if (matrix.compact_col_starts) |source_starts|
        return fillTransposeEntries(foundation.HUInt, matrix.num_cols, source_starts, matrix.row_indices, matrix.values, next, rows, values);
    return fillTransposeEntries(usize, matrix.num_cols, matrix.col_starts, matrix.row_indices, matrix.values, next, rows, values);
}

/// Like `transposeIntoAssumeValid` but uses `[]HUInt` (4-byte) for the internal
/// histogram and prefix-sum phases, matching C++ `HighsInt` data width.
/// `starts_compact` is the work array; `starts_out` receives the final `usize`
/// offsets.  Callers that already own a `[]HUInt` scratch region (e.g. the
/// `compact_starts` field of `TransposeBuffers`) can reuse it here.
pub fn transposeIntoAssumeValidCompact(matrix: csc.CscMatrix, starts_compact: []foundation.HUInt, starts_out: []usize, rows: []foundation.RowId, values: []f64, cursor_scratch: []foundation.HUInt) csc.MatrixError!void {
    if (starts_compact.len != matrix.num_rows + 1 or starts_out.len != matrix.num_rows + 1 or
        rows.len != matrix.nnz() or values.len != matrix.nnz() or
        cursor_scratch.len < matrix.num_rows)
        return error.DimensionMismatch;

    // Histogram on HUInt (4 bytes per write, matches C++)
    memory.clearTyped(foundation.HUInt, starts_compact);
    for (matrix.row_indices) |row_id| starts_compact[row_id.toUsize() + 1] += 1;

    // Prefix sum on HUInt
    for (0..matrix.num_rows) |row| starts_compact[row + 1] += starts_compact[row];

    // Copy cursor (same type, no cast)
    const next = cursor_scratch[0..matrix.num_rows];
    @memcpy(next, starts_compact[0..matrix.num_rows]);

    // Convert HUInt → usize for output
    for (starts_out[0 .. matrix.num_rows + 1], starts_compact[0 .. matrix.num_rows + 1]) |*dest, src| dest.* = @intCast(src);

    if (matrix.compact_col_starts) |source_starts|
        return fillTransposeEntries(foundation.HUInt, matrix.num_cols, source_starts, matrix.row_indices, matrix.values, next, rows, values);
    return fillTransposeEntries(usize, matrix.num_cols, matrix.col_starts, matrix.row_indices, matrix.values, next, rows, values);
}

noinline fn fillTransposeEntries(
    comptime Offset: type,
    num_cols: usize,
    source_starts: []const Offset,
    source_rows: []const foundation.RowId,
    source_values: []const f64,
    next: []foundation.HUInt,
    rows: []foundation.RowId,
    values: []f64,
) void {
    for (0..num_cols) |source_col| {
        const target_row = foundation.RowId.fromUsize(source_col) catch unreachable;
        const begin: usize = @intCast(source_starts[source_col]);
        const end: usize = @intCast(source_starts[source_col + 1]);
        for (begin..end) |position| {
            const source_value = source_values[position];
            const target_col = source_rows[position].toUsize();
            const destination: usize = @intCast(next[target_col]);
            next[target_col] += 1;
            rows[destination] = target_row;
            values[destination] = source_value;
        }
    }
}

test "explicit transpose is canonical and swaps dimensions" {
    var starts = [_]usize{ 0, 2, 3, 5 };
    var rows = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(2), try foundation.RowId.init(1), try foundation.RowId.init(0), try foundation.RowId.init(2) };
    var values = [_]f64{ 2.0, 3.0, 4.0, -1.0, 5.0 };
    const matrix: csc.CscMatrix = .{ .num_rows = 3, .num_cols = 3, .col_starts = &starts, .row_indices = &rows, .values = &values };

    var result = try transpose(std.testing.allocator, matrix);
    defer result.deinit(std.testing.allocator);
    try result.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 2, 3, 5 }, result.col_starts);
    try std.testing.expectEqualSlices(f64, &.{ 2.0, -1.0, 4.0, 3.0, 5.0 }, result.values);
}

test "transposing twice reproduces canonical CSC exactly" {
    const builder_module = @import("builder.zig");
    var builder = try builder_module.MatrixBuilder.init(5, 7);
    defer builder.deinit(std.testing.allocator);
    var prng = std.Random.DefaultPrng.init(0x7a11_2026);
    const random = prng.random();
    for (0..120) |_| {
        const row = random.intRangeLessThan(usize, 0, 5);
        const col = random.intRangeLessThan(usize, 0, 7);
        const value: f64 = @floatFromInt(random.intRangeAtMost(i8, -3, 3));
        try builder.append(std.testing.allocator, try foundation.RowId.fromUsize(row), try foundation.ColId.fromUsize(col), value);
    }
    var original = try builder.freeze(std.testing.allocator, 0.0);
    defer original.deinit(std.testing.allocator);
    var once = try transposeAssumeValid(std.testing.allocator, original);
    defer once.deinit(std.testing.allocator);
    var twice = try transposeAssumeValid(std.testing.allocator, once);
    defer twice.deinit(std.testing.allocator);

    try std.testing.expectEqualSlices(usize, original.col_starts, twice.col_starts);
    try std.testing.expectEqualSlices(foundation.RowId, original.row_indices, twice.row_indices);
    try std.testing.expectEqualSlices(f64, original.values, twice.values);
}

test "transposeLean handles odd nnz" {
    var starts = [_]usize{ 0, 1, 3, 4 };
    var rows = [_]foundation.RowId{ try foundation.RowId.init(1), try foundation.RowId.init(0), try foundation.RowId.init(2), try foundation.RowId.init(1) };
    var values = [_]f64{ 3.0, 4.0, 5.0, 6.0 };
    const matrix: csc.CscMatrix = .{ .num_rows = 3, .num_cols = 3, .col_starts = &starts, .row_indices = &rows, .values = &values };

    var result = try transposeLeanAssumeValid(std.testing.allocator, matrix);
    defer result.deinit(std.testing.allocator);
    try result.validate();
    try std.testing.expectEqual(@as(usize, 4), result.nnz());
    try std.testing.expect(result.compact_col_starts == null);
    // Verify values survived: (1,0)=3→(0,1)=3, (0,1)=4→(1,0)=4, (2,1)=5→(1,2)=5, (1,2)=6→(2,1)=6
    var expected_starts = [_]usize{ 0, 1, 3, 4 };
    try std.testing.expectEqualSlices(usize, &expected_starts, result.col_starts);
}

test "transposeLean handles empty matrix (zero nnz)" {
    var starts = [_]usize{0} ** 4; // [0,0,0,0]
    var rows = [_]foundation.RowId{};
    var values = [_]f64{};
    const matrix: csc.CscMatrix = .{ .num_rows = 2, .num_cols = 3, .col_starts = &starts, .row_indices = &rows, .values = &values };

    var result = try transposeLeanAssumeValid(std.testing.allocator, matrix);
    defer result.deinit(std.testing.allocator);
    try result.validate();
    try std.testing.expectEqual(@as(usize, 0), result.nnz());
    try std.testing.expect(result.compact_col_starts == null);
    // Transposed empty 2×3 → 3×2, all col_starts are zero
    try std.testing.expectEqualSlices(usize, &.{ 0, 0, 0 }, result.col_starts);
}

test "transposeLean handles single-column odd-nnz matrix" {
    var starts = [_]usize{ 0, 3 };
    var rows = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(2), try foundation.RowId.init(4) };
    var values = [_]f64{ 1.0, 2.0, 3.0 };
    const matrix: csc.CscMatrix = .{ .num_rows = 5, .num_cols = 1, .col_starts = &starts, .row_indices = &rows, .values = &values };

    var result = try transposeLeanAssumeValid(std.testing.allocator, matrix);
    defer result.deinit(std.testing.allocator);
    try result.validate();
    try std.testing.expectEqual(@as(usize, 3), result.nnz());
    try std.testing.expect(result.compact_col_starts == null);
}
