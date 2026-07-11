//! Explicit CSC transpose construction.
//!
//! The conversion uses counting plus prefix sums. It is O(rows + nnz), and the
//! output is canonical without sorting because source columns are visited in
//! increasing order.

const std = @import("std");
const foundation = @import("foundation");
const csc = @import("csc.zig");

/// Checked transpose for matrices entering from an untrusted boundary.
pub fn transpose(allocator: std.mem.Allocator, matrix: csc.CscMatrix) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    try matrix.validate();
    return transposeAssumeValid(allocator, matrix);
}

/// Transposes a canonical CSC matrix without repeating structural validation.
pub fn transposeAssumeValid(allocator: std.mem.Allocator, matrix: csc.CscMatrix) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    if (matrix.num_rows == std.math.maxInt(usize)) return error.DimensionTooLarge;

    const transposed_starts = try allocator.alloc(usize, matrix.num_rows + 1);
    errdefer allocator.free(transposed_starts);
    const transposed_rows = try allocator.alloc(foundation.RowId, matrix.nnz());
    errdefer allocator.free(transposed_rows);
    const transposed_values = try allocator.alloc(f64, matrix.nnz());
    errdefer allocator.free(transposed_values);
    const next = try allocator.alloc(usize, matrix.num_rows);
    defer allocator.free(next);
    try transposeIntoAssumeValid(&matrix, transposed_starts, transposed_rows, transposed_values, next);

    return .{
        .num_rows = matrix.num_cols,
        .num_cols = matrix.num_rows,
        .col_starts = transposed_starts,
        .row_indices = transposed_rows,
        .values = transposed_values,
    };
}

pub fn transposeInto(matrix: csc.CscMatrix, starts: []usize, rows: []foundation.RowId, values: []f64, cursor_scratch: []usize) csc.MatrixError!void {
    try matrix.validate();
    return transposeIntoAssumeValid(&matrix, starts, rows, values, cursor_scratch);
}

/// Allocation-free explicit transpose into exact-size caller buffers.
pub fn transposeIntoAssumeValid(matrix: *const csc.CscMatrix, starts: []usize, rows: []foundation.RowId, values: []f64, cursor_scratch: []usize) csc.MatrixError!void {
    if (starts.len != matrix.num_rows + 1 or rows.len != matrix.nnz() or
        values.len != matrix.nnz() or cursor_scratch.len < matrix.num_rows)
        return error.DimensionMismatch;
    @memset(starts, 0);
    for (matrix.row_indices) |row_id| starts[row_id.toUsize() + 1] += 1;
    for (0..matrix.num_rows) |row| starts[row + 1] += starts[row];
    const next = cursor_scratch[0..matrix.num_rows];
    @memcpy(next, starts[0..matrix.num_rows]);
    for (0..matrix.num_cols) |source_col| {
        const target_row = foundation.RowId.fromUsize(source_col) catch unreachable;
        for (matrix.col_starts[source_col]..matrix.col_starts[source_col + 1]) |position| {
            const target_col = matrix.row_indices[position].toUsize();
            const destination = next[target_col];
            rows[destination] = target_row;
            values[destination] = matrix.values[position];
            next[target_col] += 1;
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
