//! Structural sparse matrix edits that return a new canonical CSC matrix.
//!
//! Functions do not mutate the source. This gives callers strong failure
//! safety: an allocation or validation error leaves the authoritative matrix
//! untouched, and Model can swap revisions only after successful construction.

const std = @import("std");
const foundation = @import("foundation");
const sparse_vector = @import("sparse_vector.zig");
const slice = @import("slice.zig");
const csc = @import("csc.zig");

pub fn appendColumns(
    allocator: std.mem.Allocator,
    matrix: csc.CscMatrix,
    columns: []const sparse_vector.SparseVectorView(foundation.RowId),
) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    try matrix.validate();
    const new_num_cols = std.math.add(usize, matrix.num_cols, columns.len) catch return error.DimensionTooLarge;
    try csc.validateDimensions(matrix.num_rows, new_num_cols);
    var added_nnz: usize = 0;
    for (columns) |column| {
        if (column.dimension != matrix.num_rows) return error.DimensionMismatch;
        try column.validate();
        added_nnz = std.math.add(usize, added_nnz, column.nnz()) catch return error.DimensionTooLarge;
    }
    const total_nnz = std.math.add(usize, matrix.nnz(), added_nnz) catch return error.DimensionTooLarge;

    var result = try csc.CscMatrix.initPackedUninitialized(allocator, matrix.num_rows, new_num_cols, total_nnz);
    errdefer result.deinit(allocator);
    const starts = result.col_starts;
    const rows = result.row_indices;
    const values = result.values;
    @memcpy(starts[0 .. matrix.num_cols + 1], matrix.col_starts);
    var running = matrix.nnz();
    for (columns, 0..) |column, offset| {
        running += column.nnz();
        starts[matrix.num_cols + offset + 1] = running;
    }
    @memcpy(rows[0..matrix.nnz()], matrix.row_indices);
    @memcpy(values[0..matrix.nnz()], matrix.values);
    var destination = matrix.nnz();
    for (columns) |column| {
        @memcpy(rows[destination..][0..column.nnz()], column.indices);
        @memcpy(values[destination..][0..column.nnz()], column.values);
        destination += column.nnz();
    }
    return result;
}

pub fn appendRows(
    allocator: std.mem.Allocator,
    matrix: csc.CscMatrix,
    rows_to_add: []const sparse_vector.SparseVectorView(foundation.ColId),
) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    try matrix.validate();
    const new_num_rows = std.math.add(usize, matrix.num_rows, rows_to_add.len) catch return error.DimensionTooLarge;
    try csc.validateDimensions(new_num_rows, matrix.num_cols);
    var added_nnz: usize = 0;
    for (rows_to_add) |row| {
        if (row.dimension != matrix.num_cols) return error.DimensionMismatch;
        try row.validate();
        added_nnz = std.math.add(usize, added_nnz, row.nnz()) catch return error.DimensionTooLarge;
    }
    const total_nnz = std.math.add(usize, matrix.nnz(), added_nnz) catch return error.DimensionTooLarge;

    var result = try csc.CscMatrix.initPackedUninitialized(allocator, new_num_rows, matrix.num_cols, total_nnz);
    errdefer result.deinit(allocator);
    const starts = result.col_starts;
    const rows = result.row_indices;
    const values = result.values;
    @memset(starts, 0);
    for (0..matrix.num_cols) |col|
        starts[col + 1] = matrix.col_starts[col + 1] - matrix.col_starts[col];
    for (rows_to_add) |row| {
        for (row.indices) |col_id| starts[col_id.toUsize() + 1] += 1;
    }
    for (0..matrix.num_cols) |col| starts[col + 1] += starts[col];

    // Existing row IDs are smaller than every appended row ID. Copying the old
    // columns first and then visiting new rows in order preserves row sorting.
    for (0..matrix.num_cols) |col| {
        const begin = matrix.col_starts[col];
        const end = matrix.col_starts[col + 1];
        const count = end - begin;
        const destination = starts[col];
        @memcpy(rows[destination..][0..count], matrix.row_indices[begin..end]);
        @memcpy(values[destination..][0..count], matrix.values[begin..end]);
        starts[col] += count;
    }
    for (rows_to_add, 0..) |row, added_row| {
        const row_id = foundation.RowId.fromUsize(matrix.num_rows + added_row) catch unreachable;
        for (row.indices, row.values) |col_id, value| {
            const col = col_id.toUsize();
            const destination = starts[col];
            rows[destination] = row_id;
            values[destination] = value;
            starts[col] += 1;
        }
    }

    var col = matrix.num_cols;
    while (col > 0) {
        starts[col] = starts[col - 1];
        col -= 1;
    }
    starts[0] = 0;
    return result;
}

/// Appends canonical CSR-like row streams beneath a CSC matrix without first
/// materializing an array of per-row SparseVectorView descriptors.
pub fn appendRowsFromCsr(
    allocator: std.mem.Allocator,
    matrix: csc.CscMatrix,
    row_starts: []const usize,
    col_indices: []const foundation.ColId,
    stream_values: []const f64,
) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    try matrix.validate();
    try validateCsrRowStreams(matrix.num_cols, row_starts, col_indices, stream_values);
    return appendRowsFromCsrAssumeValid(allocator, matrix, row_starts, col_indices, stream_values);
}

/// Trusted allocation path for validated canonical CSR-like row streams.
/// `starts[0..num_cols]` doubles as the scatter cursor and is restored after
/// filling, avoiding a separate O(num_cols) cursor allocation.
pub fn appendRowsFromCsrAssumeValid(
    allocator: std.mem.Allocator,
    matrix: csc.CscMatrix,
    row_starts: []const usize,
    col_indices: []const foundation.ColId,
    stream_values: []const f64,
) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    std.debug.assert(row_starts.len > 0);
    std.debug.assert(row_starts[0] == 0);
    std.debug.assert(row_starts[row_starts.len - 1] == col_indices.len);
    std.debug.assert(col_indices.len == stream_values.len);

    const added_rows = row_starts.len - 1;
    const new_num_rows = std.math.add(usize, matrix.num_rows, added_rows) catch return error.DimensionTooLarge;
    try csc.validateDimensions(new_num_rows, matrix.num_cols);
    const total_nnz = std.math.add(usize, matrix.nnz(), stream_values.len) catch return error.DimensionTooLarge;

    var result = try csc.CscMatrix.initPackedUninitialized(allocator, new_num_rows, matrix.num_cols, total_nnz);
    errdefer result.deinit(allocator);
    const starts = result.col_starts;
    const rows = result.row_indices;
    const values = result.values;
    @memset(starts, 0);
    for (0..matrix.num_cols) |col|
        starts[col + 1] = matrix.col_starts[col + 1] - matrix.col_starts[col];
    for (col_indices) |col_id| starts[col_id.toUsize() + 1] += 1;
    for (0..matrix.num_cols) |col| starts[col + 1] += starts[col];

    // Existing row IDs precede every appended row ID. Fill them first, then
    // visit appended rows in increasing order to preserve canonical sorting.
    for (0..matrix.num_cols) |col| {
        const begin = matrix.col_starts[col];
        const end = matrix.col_starts[col + 1];
        const count = end - begin;
        const destination = starts[col];
        @memcpy(rows[destination..][0..count], matrix.row_indices[begin..end]);
        @memcpy(values[destination..][0..count], matrix.values[begin..end]);
        starts[col] += count;
    }
    for (0..added_rows) |added_row| {
        const row_id = foundation.RowId.fromUsize(matrix.num_rows + added_row) catch unreachable;
        for (row_starts[added_row]..row_starts[added_row + 1]) |position| {
            const col = col_indices[position].toUsize();
            const destination = starts[col];
            rows[destination] = row_id;
            values[destination] = stream_values[position];
            starts[col] += 1;
        }
    }

    // Cursor col now holds the original final offset of column col+1. Shift
    // those ends right by one slot to recover canonical column starts.
    var col = matrix.num_cols;
    while (col > 0) {
        starts[col] = starts[col - 1];
        col -= 1;
    }
    starts[0] = 0;

    return result;
}

fn validateCsrRowStreams(
    num_cols: usize,
    row_starts: []const usize,
    col_indices: []const foundation.ColId,
    values: []const f64,
) csc.MatrixError!void {
    if (row_starts.len == 0 or row_starts[0] != 0) return error.InvalidRowStarts;
    if (col_indices.len != values.len) return error.InconsistentStorage;
    if (row_starts[row_starts.len - 1] != values.len) return error.InvalidRowStarts;
    for (0..row_starts.len - 1) |row| {
        const begin = row_starts[row];
        const end = row_starts[row + 1];
        if (begin > end or end > values.len) return error.InvalidRowStarts;
        var previous_col: ?usize = null;
        for (begin..end) |position| {
            const col = col_indices[position].toUsize();
            if (col >= num_cols) return error.IndexOutOfBounds;
            if (previous_col) |previous| if (col <= previous) return error.IndicesNotStrictlyIncreasing;
            if (!std.math.isFinite(values[position])) return error.NonFiniteValue;
            if (values[position] == 0.0) return error.ExplicitZero;
            previous_col = col;
        }
    }
}

pub fn deleteColumns(allocator: std.mem.Allocator, matrix: csc.CscMatrix, deleted: []const foundation.ColId) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    try matrix.validate();
    try validateDeleted(foundation.ColId, deleted, matrix.num_cols);
    const kept = try allocator.alloc(foundation.ColId, matrix.num_cols - deleted.len);
    defer allocator.free(kept);
    var deleted_position: usize = 0;
    var kept_position: usize = 0;
    for (0..matrix.num_cols) |col| {
        if (deleted_position < deleted.len and deleted[deleted_position].toUsize() == col) {
            deleted_position += 1;
        } else {
            kept[kept_position] = foundation.ColId.fromUsize(col) catch unreachable;
            kept_position += 1;
        }
    }
    return slice.extractColumnsAssumeValid(allocator, matrix, kept);
}

pub fn deleteRows(allocator: std.mem.Allocator, matrix: csc.CscMatrix, deleted: []const foundation.RowId) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    try matrix.validate();
    try validateDeleted(foundation.RowId, deleted, matrix.num_rows);
    const kept = try allocator.alloc(foundation.RowId, matrix.num_rows - deleted.len);
    defer allocator.free(kept);
    var deleted_position: usize = 0;
    var kept_position: usize = 0;
    for (0..matrix.num_rows) |row| {
        if (deleted_position < deleted.len and deleted[deleted_position].toUsize() == row) {
            deleted_position += 1;
        } else {
            kept[kept_position] = foundation.RowId.fromUsize(row) catch unreachable;
            kept_position += 1;
        }
    }
    return slice.extractRowsAssumeValid(allocator, matrix, kept);
}

fn validateDeleted(comptime Id: type, deleted: []const Id, dimension: usize) csc.MatrixError!void {
    var previous: ?usize = null;
    for (deleted) |id| {
        const index = id.toUsize();
        if (index >= dimension) return error.IndexOutOfBounds;
        if (previous) |old| if (index <= old) return error.IndicesNotStrictlyIncreasing;
        previous = index;
    }
}

test "append canonical columns and rows" {
    var matrix = try csc.CscMatrix.initZero(std.testing.allocator, 2, 2);
    defer matrix.deinit(std.testing.allocator);
    var column_rows = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(1) };
    var column_values = [_]f64{ 2.0, 3.0 };
    const columns = [_]sparse_vector.SparseVectorView(foundation.RowId){.{ .dimension = 2, .indices = &column_rows, .values = &column_values }};
    var with_column = try appendColumns(std.testing.allocator, matrix, &columns);
    defer with_column.deinit(std.testing.allocator);
    try with_column.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 0, 0, 2 }, with_column.col_starts);

    var row_cols = [_]foundation.ColId{ try foundation.ColId.init(0), try foundation.ColId.init(2) };
    var row_values = [_]f64{ 4.0, 5.0 };
    const new_rows = [_]sparse_vector.SparseVectorView(foundation.ColId){.{ .dimension = 3, .indices = &row_cols, .values = &row_values }};
    var completed = try appendRows(std.testing.allocator, with_column, &new_rows);
    defer completed.deinit(std.testing.allocator);
    try completed.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 1, 1, 4 }, completed.col_starts);
    try std.testing.expectEqualSlices(f64, &.{ 4.0, 2.0, 3.0, 5.0 }, completed.values);
}

test "append CSR row streams directly including an empty row" {
    const builder_module = @import("builder.zig");
    var builder = try builder_module.MatrixBuilder.init(2, 3);
    defer builder.deinit(std.testing.allocator);
    try builder.append(std.testing.allocator, try foundation.RowId.init(0), try foundation.ColId.init(0), 1.0);
    try builder.append(std.testing.allocator, try foundation.RowId.init(1), try foundation.ColId.init(2), 2.0);
    var base = try builder.freezeSortedLeanAssumeValid(std.testing.allocator, 0.0);
    defer base.deinit(std.testing.allocator);

    const row_starts = [_]usize{ 0, 2, 2, 3 };
    const cols = [_]foundation.ColId{
        try foundation.ColId.init(0),
        try foundation.ColId.init(2),
        try foundation.ColId.init(1),
    };
    const stream_values = [_]f64{ 3.0, 4.0, 5.0 };
    var result = try appendRowsFromCsr(std.testing.allocator, base, &row_starts, &cols, &stream_values);
    defer result.deinit(std.testing.allocator);
    try result.validate();
    try std.testing.expectEqual(@as(usize, 5), result.num_rows);
    try std.testing.expectEqualSlices(usize, &.{ 0, 2, 3, 5 }, result.col_starts);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 3.0, 5.0, 2.0, 4.0 }, result.values);
    try std.testing.expectEqual(@as(usize, 2), result.row_indices[1].toUsize());
    try std.testing.expectEqual(@as(usize, 4), result.row_indices[2].toUsize());
}

test "append CSR row streams rejects malformed canonical data" {
    var base = try csc.CscMatrix.initZero(std.testing.allocator, 1, 2);
    defer base.deinit(std.testing.allocator);
    const duplicate_cols = [_]foundation.ColId{ try foundation.ColId.init(0), try foundation.ColId.init(0) };
    try std.testing.expectError(
        error.IndicesNotStrictlyIncreasing,
        appendRowsFromCsr(std.testing.allocator, base, &.{ 0, 2 }, &duplicate_cols, &.{ 1.0, 2.0 }),
    );
    try std.testing.expectError(
        error.InvalidRowStarts,
        appendRowsFromCsr(std.testing.allocator, base, &.{ 1, 1 }, &.{}, &.{}),
    );
}

test "delete rows and columns remaps dimensions canonically" {
    const builder_module = @import("builder.zig");
    var builder = try builder_module.MatrixBuilder.init(3, 3);
    defer builder.deinit(std.testing.allocator);
    try builder.append(std.testing.allocator, try foundation.RowId.init(0), try foundation.ColId.init(0), 1.0);
    try builder.append(std.testing.allocator, try foundation.RowId.init(1), try foundation.ColId.init(1), 2.0);
    try builder.append(std.testing.allocator, try foundation.RowId.init(2), try foundation.ColId.init(2), 3.0);
    var matrix = try builder.freeze(std.testing.allocator, 0.0);
    defer matrix.deinit(std.testing.allocator);
    var without_row = try deleteRows(std.testing.allocator, matrix, &.{try foundation.RowId.init(1)});
    defer without_row.deinit(std.testing.allocator);
    var result = try deleteColumns(std.testing.allocator, without_row, &.{try foundation.ColId.init(0)});
    defer result.deinit(std.testing.allocator);
    try result.validate();
    try std.testing.expectEqual(@as(usize, 2), result.num_rows);
    try std.testing.expectEqual(@as(usize, 2), result.num_cols);
    try std.testing.expectEqualSlices(f64, &.{3.0}, result.values);
    try std.testing.expectEqual(@as(usize, 1), result.row_indices[0].toUsize());
}

test "structural edits validate vectors and deletion sets" {
    var matrix = try csc.CscMatrix.initZero(std.testing.allocator, 2, 2);
    defer matrix.deinit(std.testing.allocator);
    var no_ids = [_]foundation.RowId{};
    var no_values = [_]f64{};
    const wrong_columns = [_]sparse_vector.SparseVectorView(foundation.RowId){.{ .dimension = 3, .indices = &no_ids, .values = &no_values }};
    try std.testing.expectError(error.DimensionMismatch, appendColumns(std.testing.allocator, matrix, &wrong_columns));
    const duplicate = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(0) };
    try std.testing.expectError(error.IndicesNotStrictlyIncreasing, deleteRows(std.testing.allocator, matrix, &duplicate));
}
