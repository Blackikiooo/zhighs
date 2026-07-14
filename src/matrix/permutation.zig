//! Row and column permutation of canonical CSC matrices.
//!
//! A permutation slice maps old index to new strongly typed ID. Checked entry
//! points verify a full bijection. Column movement is direct; row permutation
//! requires sorting entries inside each destination column to restore canonical
//! increasing row order.

const std = @import("std");
const foundation = @import("foundation");
const csc = @import("csc.zig");
const transform_buffers = @import("transform_buffers.zig");

const CscTransformBuffers = transform_buffers.CscTransformBuffers;

pub fn permute(
    allocator: std.mem.Allocator,
    matrix: csc.CscMatrix,
    row_old_to_new: []const foundation.RowId,
    col_old_to_new: []const foundation.ColId,
) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    try matrix.validate();
    try validatePermutation(foundation.RowId, allocator, row_old_to_new, matrix.num_rows);
    try validatePermutation(foundation.ColId, allocator, col_old_to_new, matrix.num_cols);
    return permuteAssumeValid(allocator, matrix, row_old_to_new, col_old_to_new);
}

pub fn permuteAssumeValid(
    allocator: std.mem.Allocator,
    matrix: csc.CscMatrix,
    row_old_to_new: []const foundation.RowId,
    col_old_to_new: []const foundation.ColId,
) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
    var result = try csc.CscMatrix.initPackedUninitialized(allocator, matrix.num_rows, matrix.num_cols, matrix.nnz());
    errdefer result.deinit(allocator);
    const starts = result.col_starts;
    const rows = result.row_indices;
    const values = result.values;
    @memset(starts, 0);
    for (0..matrix.num_cols) |old_col| {
        const new_col = col_old_to_new[old_col].toUsize();
        starts[new_col + 1] = matrix.col_starts[old_col + 1] - matrix.col_starts[old_col];
    }
    for (0..matrix.num_cols) |col| starts[col + 1] += starts[col];

    for (0..matrix.num_cols) |old_col| {
        const new_col = col_old_to_new[old_col].toUsize();
        var destination = starts[new_col];
        for (matrix.col_starts[old_col]..matrix.col_starts[old_col + 1]) |position| {
            rows[destination] = row_old_to_new[matrix.row_indices[position].toUsize()];
            values[destination] = matrix.values[position];
            destination += 1;
        }
    }

    const context: EntrySortContext = .{ .rows = rows, .values = values };
    for (0..matrix.num_cols) |col|
        std.sort.pdqContext(starts[col], starts[col + 1], context);

    return result;
}

pub fn permuteInto(
    buffers: *CscTransformBuffers,
    matrix: csc.CscMatrix,
    row_old_to_new: []const foundation.RowId,
    col_old_to_new: []const foundation.ColId,
) csc.MatrixError!csc.CscView {
    try matrix.validate();
    try buffers.requireCapacity(matrix.num_cols, matrix.nnz(), @max(matrix.num_rows, matrix.num_cols));
    try validatePermutationWithScratch(foundation.RowId, row_old_to_new, matrix.num_rows, buffers.index_scratch);
    try validatePermutationWithScratch(foundation.ColId, col_old_to_new, matrix.num_cols, buffers.index_scratch);
    return permuteIntoAssumeValid(buffers, matrix, row_old_to_new, col_old_to_new);
}

pub fn permuteIntoAssumeValid(
    buffers: *CscTransformBuffers,
    matrix: csc.CscMatrix,
    row_old_to_new: []const foundation.RowId,
    col_old_to_new: []const foundation.ColId,
) csc.MatrixError!csc.CscView {
    try buffers.requireCapacity(matrix.num_cols, matrix.nnz(), 0);
    const starts = buffers.col_starts[0 .. matrix.num_cols + 1];
    const rows = buffers.row_indices[0..matrix.nnz()];
    const values = buffers.values[0..matrix.nnz()];
    @memset(starts, 0);
    for (0..matrix.num_cols) |old_col| {
        const new_col = col_old_to_new[old_col].toUsize();
        starts[new_col + 1] = matrix.col_starts[old_col + 1] - matrix.col_starts[old_col];
    }
    for (0..matrix.num_cols) |col| starts[col + 1] += starts[col];

    for (0..matrix.num_cols) |old_col| {
        const new_col = col_old_to_new[old_col].toUsize();
        var destination = starts[new_col];
        for (matrix.col_starts[old_col]..matrix.col_starts[old_col + 1]) |position| {
            rows[destination] = row_old_to_new[matrix.row_indices[position].toUsize()];
            values[destination] = matrix.values[position];
            destination += 1;
        }
    }

    const context: EntrySortContext = .{ .rows = rows, .values = values };
    for (0..matrix.num_cols) |col|
        std.sort.pdqContext(starts[col], starts[col + 1], context);
    return buffers.viewAssumeValid(matrix.num_rows, matrix.num_cols, matrix.nnz());
}

fn validatePermutation(comptime Id: type, allocator: std.mem.Allocator, old_to_new: []const Id, dimension: usize) (std.mem.Allocator.Error || csc.MatrixError)!void {
    if (old_to_new.len != dimension) return error.InvalidPermutation;
    const seen = try allocator.alloc(bool, dimension);
    defer allocator.free(seen);
    @memset(seen, false);
    for (old_to_new) |new_id| {
        const new_index = new_id.toUsize();
        if (new_index >= dimension or seen[new_index]) return error.InvalidPermutation;
        seen[new_index] = true;
    }
}

fn validatePermutationWithScratch(comptime Id: type, old_to_new: []const Id, dimension: usize, scratch: []usize) csc.MatrixError!void {
    if (old_to_new.len != dimension) return error.InvalidPermutation;
    if (scratch.len < dimension) return error.BufferTooSmall;
    const seen = scratch[0..dimension];
    @memset(seen, 0);
    for (old_to_new) |new_id| {
        const new_index = new_id.toUsize();
        if (new_index >= dimension or seen[new_index] != 0) return error.InvalidPermutation;
        seen[new_index] = 1;
    }
}

const EntrySortContext = struct {
    rows: []foundation.RowId,
    values: []f64,

    pub fn lessThan(self: @This(), a: usize, b: usize) bool {
        return self.rows[a].toUsize() < self.rows[b].toUsize();
    }

    pub fn swap(self: @This(), a: usize, b: usize) void {
        std.mem.swap(foundation.RowId, &self.rows[a], &self.rows[b]);
        std.mem.swap(f64, &self.values[a], &self.values[b]);
    }
};

test "row and column permutation preserves values at mapped coordinates" {
    // Original rows: [2,0,-1], [0,4,0], [3,0,5].
    var starts = [_]usize{ 0, 2, 3, 5 };
    var rows = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(2), try foundation.RowId.init(1), try foundation.RowId.init(0), try foundation.RowId.init(2) };
    var values = [_]f64{ 2.0, 3.0, 4.0, -1.0, 5.0 };
    const matrix = csc.CscMatrix.initBorrowedAssumeValid(3, 3, &starts, &rows, &values);
    const row_map = [_]foundation.RowId{ try foundation.RowId.init(2), try foundation.RowId.init(0), try foundation.RowId.init(1) };
    const col_map = [_]foundation.ColId{ try foundation.ColId.init(1), try foundation.ColId.init(2), try foundation.ColId.init(0) };
    var result = try permute(std.testing.allocator, matrix, &row_map, &col_map);
    defer result.deinit(std.testing.allocator);
    try result.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 2, 4, 5 }, result.col_starts);
    // New column 0 is old column 2: rows old 0->2 and old 2->1, sorted 1,2.
    try std.testing.expectEqualSlices(f64, &.{ 5.0, -1.0, 3.0, 2.0, 4.0 }, result.values);
}

test "permutation must be a dimension-sized bijection" {
    var matrix = try csc.CscMatrix.initZero(std.testing.allocator, 2, 2);
    defer matrix.deinit(std.testing.allocator);
    const valid_rows = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(1) };
    const duplicate_cols = [_]foundation.ColId{ try foundation.ColId.init(0), try foundation.ColId.init(0) };
    try std.testing.expectError(error.InvalidPermutation, permute(std.testing.allocator, matrix, &valid_rows, &duplicate_cols));
    const short_rows = [_]foundation.RowId{try foundation.RowId.init(0)};
    const valid_cols = [_]foundation.ColId{ try foundation.ColId.init(0), try foundation.ColId.init(1) };
    try std.testing.expectError(error.InvalidPermutation, permute(std.testing.allocator, matrix, &short_rows, &valid_cols));
}

test "permutation followed by inverse reproduces matrix" {
    const builder_module = @import("builder.zig");
    var builder = try builder_module.MatrixBuilder.init(3, 4);
    defer builder.deinit(std.testing.allocator);
    try builder.append(std.testing.allocator, try foundation.RowId.init(0), try foundation.ColId.init(1), 2.0);
    try builder.append(std.testing.allocator, try foundation.RowId.init(2), try foundation.ColId.init(3), -7.0);
    try builder.append(std.testing.allocator, try foundation.RowId.init(1), try foundation.ColId.init(0), 5.0);
    var original = try builder.freeze(std.testing.allocator, 0.0);
    defer original.deinit(std.testing.allocator);
    const rows_forward = [_]foundation.RowId{ try foundation.RowId.init(2), try foundation.RowId.init(0), try foundation.RowId.init(1) };
    const rows_inverse = [_]foundation.RowId{ try foundation.RowId.init(1), try foundation.RowId.init(2), try foundation.RowId.init(0) };
    const cols_forward = [_]foundation.ColId{ try foundation.ColId.init(2), try foundation.ColId.init(0), try foundation.ColId.init(3), try foundation.ColId.init(1) };
    const cols_inverse = [_]foundation.ColId{ try foundation.ColId.init(1), try foundation.ColId.init(3), try foundation.ColId.init(0), try foundation.ColId.init(2) };
    var changed = try permute(std.testing.allocator, original, &rows_forward, &cols_forward);
    defer changed.deinit(std.testing.allocator);
    var restored = try permute(std.testing.allocator, changed, &rows_inverse, &cols_inverse);
    defer restored.deinit(std.testing.allocator);
    try std.testing.expectEqualSlices(usize, original.col_starts, restored.col_starts);
    try std.testing.expectEqualSlices(foundation.RowId, original.row_indices, restored.row_indices);
    try std.testing.expectEqualSlices(f64, original.values, restored.values);
}

test "permutation Into reuses aligned output and validation scratch" {
    var starts = [_]usize{ 0, 2, 3, 5 };
    var rows = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(2), try foundation.RowId.init(1), try foundation.RowId.init(0), try foundation.RowId.init(2) };
    var values = [_]f64{ 2.0, 3.0, 4.0, -1.0, 5.0 };
    const matrix = csc.CscMatrix.initBorrowedAssumeValid(3, 3, &starts, &rows, &values);
    const row_map = [_]foundation.RowId{ try foundation.RowId.init(2), try foundation.RowId.init(0), try foundation.RowId.init(1) };
    const col_map = [_]foundation.ColId{ try foundation.ColId.init(1), try foundation.ColId.init(2), try foundation.ColId.init(0) };
    var buffers = try CscTransformBuffers.initCapacity(std.testing.allocator, 3, 5, 3);
    defer buffers.deinit(std.testing.allocator);
    const view = try permuteInto(&buffers, matrix, &row_map, &col_map);
    try view.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 2, 4, 5 }, view.col_starts);
    try std.testing.expectEqualSlices(f64, &.{ 5.0, -1.0, 3.0, 2.0, 4.0 }, view.values);

    const duplicate_cols = [_]foundation.ColId{ try foundation.ColId.init(0), try foundation.ColId.init(0), try foundation.ColId.init(2) };
    try std.testing.expectError(error.InvalidPermutation, permuteInto(&buffers, matrix, &row_map, &duplicate_cols));
}
