//! Canonical row and column extraction from CSC matrices.
//!
//! Selected IDs must be strictly increasing. This contract makes the output
//! ordering deterministic and lets extraction preserve canonical order without
//! sorting or hash tables.

const std = @import("std");
const foundation = @import("foundation");
const csc = @import("csc.zig");
const transform_buffers = @import("transform_buffers.zig");

const CscTransformBuffers = transform_buffers.CscTransformBuffers;

/// Extracts the half-open contiguous column interval [from_col, to_col).
pub fn extractColumnRange(allocator: std.mem.Allocator, matrix: csc.CscMatrix, from_col: usize, to_col: usize) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    try matrix.validate();
    if (from_col > to_col or to_col > matrix.num_cols) return error.IndexOutOfBounds;
    const output_cols = to_col - from_col;
    const source_begin = matrix.col_starts[from_col];
    const source_end = matrix.col_starts[to_col];
    var result = try csc.CscMatrix.initPackedUninitialized(allocator, matrix.num_rows, output_cols, source_end - source_begin);
    errdefer result.deinit(allocator);
    const starts = result.col_starts;
    for (0..output_cols + 1) |index| starts[index] = matrix.col_starts[from_col + index] - source_begin;
    @memcpy(result.row_indices, matrix.row_indices[source_begin..source_end]);
    @memcpy(result.values, matrix.values[source_begin..source_end]);
    return result;
}

/// Writes a contiguous column range into reusable caller-owned storage.
/// The returned view is invalidated by the next write to `buffers` or deinit.
pub fn extractColumnRangeInto(buffers: *CscTransformBuffers, matrix: csc.CscMatrix, from_col: usize, to_col: usize) csc.MatrixError!csc.CscView {
    try matrix.validate();
    if (from_col > to_col or to_col > matrix.num_cols) return error.IndexOutOfBounds;
    const output_cols = to_col - from_col;
    const source_begin = matrix.col_starts[from_col];
    const source_end = matrix.col_starts[to_col];
    const output_nnz = source_end - source_begin;
    try buffers.requireCapacity(output_cols, output_nnz, 0);
    for (0..output_cols + 1) |index|
        buffers.col_starts[index] = matrix.col_starts[from_col + index] - source_begin;
    @memcpy(buffers.row_indices[0..output_nnz], matrix.row_indices[source_begin..source_end]);
    @memcpy(buffers.values[0..output_nnz], matrix.values[source_begin..source_end]);
    return buffers.viewAssumeValid(matrix.num_rows, output_cols, output_nnz);
}

/// Validate a selection and return an owning CSC matrix of those columns.
pub fn extractColumns(allocator: std.mem.Allocator, matrix: csc.CscMatrix, selected: []const foundation.ColId) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    try matrix.validate();
    try validateSelection(foundation.ColId, selected, matrix.num_cols);
    return extractColumnsAssumeValid(allocator, matrix, selected);
}

/// Trusted owning column extraction; selected IDs must be valid and ordered.
pub fn extractColumnsAssumeValid(allocator: std.mem.Allocator, matrix: csc.CscMatrix, selected: []const foundation.ColId) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    var output_nnz: usize = 0;
    for (selected) |source_col| {
        const col = source_col.toUsize();
        output_nnz += matrix.col_starts[col + 1] - matrix.col_starts[col];
    }

    var result = try csc.CscMatrix.initPackedUninitialized(allocator, matrix.num_rows, selected.len, output_nnz);
    errdefer result.deinit(allocator);
    const starts = result.col_starts;
    const rows = result.row_indices;
    const values = result.values;
    starts[0] = 0;
    var destination: usize = 0;
    for (selected, 0..) |source_col, target_col| {
        const col = source_col.toUsize();
        const begin = matrix.col_starts[col];
        const end = matrix.col_starts[col + 1];
        const count = end - begin;
        @memcpy(rows[destination..][0..count], matrix.row_indices[begin..end]);
        @memcpy(values[destination..][0..count], matrix.values[begin..end]);
        destination += count;
        starts[target_col + 1] = destination;
    }
    return result;
}

/// Validate and extract selected columns into retained caller buffers.
pub fn extractColumnsInto(buffers: *CscTransformBuffers, matrix: csc.CscMatrix, selected: []const foundation.ColId) csc.MatrixError!csc.CscView {
    try matrix.validate();
    try validateSelection(foundation.ColId, selected, matrix.num_cols);
    return extractColumnsIntoAssumeValid(buffers, matrix, selected);
}

/// Trusted allocation-free selected-column extraction.
pub fn extractColumnsIntoAssumeValid(buffers: *CscTransformBuffers, matrix: csc.CscMatrix, selected: []const foundation.ColId) csc.MatrixError!csc.CscView {
    var output_nnz: usize = 0;
    for (selected) |source_col| {
        const col = source_col.toUsize();
        output_nnz = std.math.add(usize, output_nnz, matrix.col_starts[col + 1] - matrix.col_starts[col]) catch return error.DimensionTooLarge;
    }
    try buffers.requireCapacity(selected.len, output_nnz, 0);

    buffers.col_starts[0] = 0;
    var destination: usize = 0;
    for (selected, 0..) |source_col, target_col| {
        const col = source_col.toUsize();
        const begin = matrix.col_starts[col];
        const end = matrix.col_starts[col + 1];
        const count = end - begin;
        @memcpy(buffers.row_indices[destination..][0..count], matrix.row_indices[begin..end]);
        @memcpy(buffers.values[destination..][0..count], matrix.values[begin..end]);
        destination += count;
        buffers.col_starts[target_col + 1] = destination;
    }
    return buffers.viewAssumeValid(matrix.num_rows, selected.len, output_nnz);
}

/// Validate a row selection and return an owning remapped CSC submatrix.
pub fn extractRows(allocator: std.mem.Allocator, matrix: csc.CscMatrix, selected: []const foundation.RowId) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    try matrix.validate();
    try validateSelection(foundation.RowId, selected, matrix.num_rows);
    return extractRowsAssumeValid(allocator, matrix, selected);
}

/// Trusted owning row extraction; selected IDs must be valid and ordered.
pub fn extractRowsAssumeValid(allocator: std.mem.Allocator, matrix: csc.CscMatrix, selected: []const foundation.RowId) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    const missing = std.math.maxInt(usize);
    const row_map = try allocator.alloc(usize, matrix.num_rows);
    defer allocator.free(row_map);
    @memset(row_map, missing);
    for (selected, 0..) |source_row, target_row| row_map[source_row.toUsize()] = target_row;

    var output_nnz: usize = 0;
    for (0..matrix.num_cols) |col| {
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position| {
            if (row_map[matrix.row_indices[position].toUsize()] != missing) output_nnz += 1;
        }
    }

    var result = try csc.CscMatrix.initPackedUninitialized(allocator, selected.len, matrix.num_cols, output_nnz);
    errdefer result.deinit(allocator);
    var destination: usize = 0;
    for (0..matrix.num_cols) |col| {
        result.col_starts[col] = destination;
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position| {
            const target_row = row_map[matrix.row_indices[position].toUsize()];
            if (target_row == missing) continue;
            result.row_indices[destination] = foundation.RowId.fromUsize(target_row) catch unreachable;
            result.values[destination] = matrix.values[position];
            destination += 1;
        }
    }
    result.col_starts[matrix.num_cols] = destination;
    return result;
}

/// Validate and extract selected rows into retained buffers.
pub fn extractRowsInto(buffers: *CscTransformBuffers, matrix: csc.CscMatrix, selected: []const foundation.RowId) csc.MatrixError!csc.CscView {
    try matrix.validate();
    try validateSelection(foundation.RowId, selected, matrix.num_rows);
    return extractRowsIntoAssumeValid(buffers, matrix, selected);
}

/// Trusted allocation-free selected-row extraction and row-ID remapping.
pub fn extractRowsIntoAssumeValid(buffers: *CscTransformBuffers, matrix: csc.CscMatrix, selected: []const foundation.RowId) csc.MatrixError!csc.CscView {
    try buffers.requireCapacity(matrix.num_cols, 0, matrix.num_rows);
    const missing = std.math.maxInt(usize);
    const row_map = buffers.index_scratch[0..matrix.num_rows];
    @memset(row_map, missing);
    for (selected, 0..) |source_row, target_row| row_map[source_row.toUsize()] = target_row;

    var output_nnz: usize = 0;
    for (matrix.row_indices) |row| {
        if (row_map[row.toUsize()] != missing)
            output_nnz = std.math.add(usize, output_nnz, 1) catch return error.DimensionTooLarge;
    }
    try buffers.requireCapacity(matrix.num_cols, output_nnz, matrix.num_rows);

    var destination: usize = 0;
    for (0..matrix.num_cols) |col| {
        buffers.col_starts[col] = destination;
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position| {
            const target_row = row_map[matrix.row_indices[position].toUsize()];
            if (target_row == missing) continue;
            buffers.row_indices[destination] = foundation.RowId.fromUsizeAssumeValid(target_row);
            buffers.values[destination] = matrix.values[position];
            destination += 1;
        }
    }
    buffers.col_starts[matrix.num_cols] = destination;
    return buffers.viewAssumeValid(selected.len, matrix.num_cols, output_nnz);
}

/// Require a strictly increasing in-range row or column selection.
fn validateSelection(comptime Id: type, selected: []const Id, dimension: usize) csc.MatrixError!void {
    var previous: ?usize = null;
    for (selected) |id| {
        const index = id.toUsize();
        if (index >= dimension) return error.IndexOutOfBounds;
        if (previous) |old| {
            if (index <= old) return error.IndicesNotStrictlyIncreasing;
        }
        previous = index;
    }
}

test "column extraction preserves requested canonical columns" {
    var starts = [_]usize{ 0, 2, 3, 5 };
    var rows = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(2), try foundation.RowId.init(1), try foundation.RowId.init(0), try foundation.RowId.init(2) };
    var values = [_]f64{ 2.0, 3.0, 4.0, -1.0, 5.0 };
    const matrix = csc.CscMatrix.initBorrowedAssumeValid(3, 3, &starts, &rows, &values);
    const selected = [_]foundation.ColId{ try foundation.ColId.init(0), try foundation.ColId.init(2) };
    var result = try extractColumns(std.testing.allocator, matrix, &selected);
    defer result.deinit(std.testing.allocator);
    try result.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 2, 4 }, result.col_starts);
    try std.testing.expectEqualSlices(f64, &.{ 2.0, 3.0, -1.0, 5.0 }, result.values);
}

test "row extraction remaps row IDs and preserves columns" {
    var starts = [_]usize{ 0, 2, 3, 5 };
    var rows = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(2), try foundation.RowId.init(1), try foundation.RowId.init(0), try foundation.RowId.init(2) };
    var values = [_]f64{ 2.0, 3.0, 4.0, -1.0, 5.0 };
    const matrix = csc.CscMatrix.initBorrowedAssumeValid(3, 3, &starts, &rows, &values);
    const selected = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(2) };
    var result = try extractRows(std.testing.allocator, matrix, &selected);
    defer result.deinit(std.testing.allocator);
    try result.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 2, 2, 4 }, result.col_starts);
    try std.testing.expectEqualSlices(f64, &.{ 2.0, 3.0, -1.0, 5.0 }, result.values);
    try std.testing.expectEqual(@as(usize, 0), result.row_indices[2].toUsize());
    try std.testing.expectEqual(@as(usize, 1), result.row_indices[3].toUsize());
}

test "empty selections and invalid selections are handled" {
    var matrix = try csc.CscMatrix.initZero(std.testing.allocator, 2, 2);
    defer matrix.deinit(std.testing.allocator);
    var no_columns = try extractColumns(std.testing.allocator, matrix, &.{});
    defer no_columns.deinit(std.testing.allocator);
    try no_columns.validate();
    var no_rows = try extractRows(std.testing.allocator, matrix, &.{});
    defer no_rows.deinit(std.testing.allocator);
    try no_rows.validate();

    const duplicate = [_]foundation.ColId{ try foundation.ColId.init(0), try foundation.ColId.init(0) };
    try std.testing.expectError(error.IndicesNotStrictlyIncreasing, extractColumns(std.testing.allocator, matrix, &duplicate));
    const out_of_bounds = [_]foundation.RowId{try foundation.RowId.init(2)};
    try std.testing.expectError(error.IndexOutOfBounds, extractRows(std.testing.allocator, matrix, &out_of_bounds));
}

test "contiguous column range copies one CSC span" {
    var starts = [_]usize{ 0, 1, 1, 3 };
    var rows = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(0), try foundation.RowId.init(1) };
    var values = [_]f64{ 1.0, 2.0, 3.0 };
    const matrix = csc.CscMatrix.initBorrowedAssumeValid(2, 3, &starts, &rows, &values);
    var result = try extractColumnRange(std.testing.allocator, matrix, 1, 3);
    defer result.deinit(std.testing.allocator);
    try result.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 0, 2 }, result.col_starts);
    try std.testing.expectEqualSlices(f64, &.{ 2.0, 3.0 }, result.values);
    try std.testing.expectError(error.IndexOutOfBounds, extractColumnRange(std.testing.allocator, matrix, 2, 1));
}

test "slice Into APIs reuse one capacity and return borrowed views" {
    var starts = [_]usize{ 0, 2, 3, 5 };
    var rows = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(2), try foundation.RowId.init(1), try foundation.RowId.init(0), try foundation.RowId.init(2) };
    var values = [_]f64{ 2.0, 3.0, 4.0, -1.0, 5.0 };
    const matrix = csc.CscMatrix.initBorrowedAssumeValid(3, 3, &starts, &rows, &values);
    var buffers = try CscTransformBuffers.initCapacity(std.testing.allocator, 3, 5, 3);
    defer buffers.deinit(std.testing.allocator);
    const starts_ptr = buffers.col_starts.ptr;
    const rows_ptr = buffers.row_indices.ptr;
    const values_ptr = buffers.values.ptr;

    const columns = [_]foundation.ColId{ try foundation.ColId.init(0), try foundation.ColId.init(2) };
    const column_view = try extractColumnsInto(&buffers, matrix, &columns);
    try column_view.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 2, 4 }, column_view.col_starts);
    try std.testing.expectEqualSlices(f64, &.{ 2.0, 3.0, -1.0, 5.0 }, column_view.values);

    const selected_rows = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(2) };
    const row_view = try extractRowsInto(&buffers, matrix, &selected_rows);
    try row_view.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 2, 2, 4 }, row_view.col_starts);
    try std.testing.expectEqual(starts_ptr, buffers.col_starts.ptr);
    try std.testing.expectEqual(rows_ptr, buffers.row_indices.ptr);
    try std.testing.expectEqual(values_ptr, buffers.values.ptr);

    const range_view = try extractColumnRangeInto(&buffers, matrix, 1, 3);
    try range_view.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 1, 3 }, range_view.col_starts);
    var too_small = try CscTransformBuffers.initCapacity(std.testing.allocator, 1, 1, 1);
    defer too_small.deinit(std.testing.allocator);
    try std.testing.expectError(error.BufferTooSmall, extractColumnsInto(&too_small, matrix, &columns));
}
