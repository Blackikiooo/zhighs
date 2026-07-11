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

    const starts = try allocator.alloc(usize, new_num_cols + 1);
    errdefer allocator.free(starts);
    @memcpy(starts[0 .. matrix.num_cols + 1], matrix.col_starts);
    var running = matrix.nnz();
    for (columns, 0..) |column, offset| {
        running += column.nnz();
        starts[matrix.num_cols + offset + 1] = running;
    }
    const rows = try allocator.alloc(foundation.RowId, total_nnz);
    errdefer allocator.free(rows);
    const values = try allocator.alloc(f64, total_nnz);
    errdefer allocator.free(values);
    @memcpy(rows[0..matrix.nnz()], matrix.row_indices);
    @memcpy(values[0..matrix.nnz()], matrix.values);
    var destination = matrix.nnz();
    for (columns) |column| {
        @memcpy(rows[destination..][0..column.nnz()], column.indices);
        @memcpy(values[destination..][0..column.nnz()], column.values);
        destination += column.nnz();
    }
    return .{ .num_rows = matrix.num_rows, .num_cols = new_num_cols, .col_starts = starts, .row_indices = rows, .values = values };
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

    const starts = try allocator.alloc(usize, matrix.num_cols + 1);
    errdefer allocator.free(starts);
    @memset(starts, 0);
    for (0..matrix.num_cols) |col|
        starts[col + 1] = matrix.col_starts[col + 1] - matrix.col_starts[col];
    for (rows_to_add) |row| {
        for (row.indices) |col_id| starts[col_id.toUsize() + 1] += 1;
    }
    for (0..matrix.num_cols) |col| starts[col + 1] += starts[col];

    const rows = try allocator.alloc(foundation.RowId, total_nnz);
    errdefer allocator.free(rows);
    const values = try allocator.alloc(f64, total_nnz);
    errdefer allocator.free(values);
    const next = try allocator.dupe(usize, starts[0..matrix.num_cols]);
    defer allocator.free(next);

    // Existing row IDs are smaller than every appended row ID. Copying the old
    // columns first and then visiting new rows in order preserves row sorting.
    for (0..matrix.num_cols) |col| {
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position| {
            const destination = next[col];
            rows[destination] = matrix.row_indices[position];
            values[destination] = matrix.values[position];
            next[col] += 1;
        }
    }
    for (rows_to_add, 0..) |row, added_row| {
        const row_id = foundation.RowId.fromUsize(matrix.num_rows + added_row) catch unreachable;
        for (row.indices, row.values) |col_id, value| {
            const col = col_id.toUsize();
            const destination = next[col];
            rows[destination] = row_id;
            values[destination] = value;
            next[col] += 1;
        }
    }
    return .{ .num_rows = new_num_rows, .num_cols = matrix.num_cols, .col_starts = starts, .row_indices = rows, .values = values };
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
