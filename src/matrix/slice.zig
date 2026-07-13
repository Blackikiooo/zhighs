//! Canonical row and column extraction from CSC matrices.
//!
//! Selected IDs must be strictly increasing. This contract makes the output
//! ordering deterministic and lets extraction preserve canonical order without
//! sorting or hash tables.

const std = @import("std");
const foundation = @import("foundation");
const csc = @import("csc.zig");

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

pub fn extractColumns(allocator: std.mem.Allocator, matrix: csc.CscMatrix, selected: []const foundation.ColId) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    try matrix.validate();
    try validateSelection(foundation.ColId, selected, matrix.num_cols);
    return extractColumnsAssumeValid(allocator, matrix, selected);
}

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

pub fn extractRows(allocator: std.mem.Allocator, matrix: csc.CscMatrix, selected: []const foundation.RowId) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    try matrix.validate();
    try validateSelection(foundation.RowId, selected, matrix.num_rows);
    return extractRowsAssumeValid(allocator, matrix, selected);
}

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
    const matrix: csc.CscMatrix = .{ .num_rows = 3, .num_cols = 3, .col_starts = &starts, .row_indices = &rows, .values = &values };
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
    const matrix: csc.CscMatrix = .{ .num_rows = 3, .num_cols = 3, .col_starts = &starts, .row_indices = &rows, .values = &values };
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
    const matrix: csc.CscMatrix = .{ .num_rows = 2, .num_cols = 3, .col_starts = &starts, .row_indices = &rows, .values = &values };
    var result = try extractColumnRange(std.testing.allocator, matrix, 1, 3);
    defer result.deinit(std.testing.allocator);
    try result.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 0, 2 }, result.col_starts);
    try std.testing.expectEqualSlices(f64, &.{ 2.0, 3.0 }, result.values);
    try std.testing.expectError(error.IndexOutOfBounds, extractColumnRange(std.testing.allocator, matrix, 2, 1));
}
