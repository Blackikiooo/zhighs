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
const transform_buffers = @import("transform_buffers.zig");

const CscTransformBuffers = transform_buffers.CscTransformBuffers;

/// Validate and append a canonical CSC column block, returning a new owner.
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

/// Validate both blocks and append columns into reusable caller buffers.
pub fn appendColumnsInto(
    buffers: *CscTransformBuffers,
    matrix: csc.CscMatrix,
    columns: []const sparse_vector.SparseVectorView(foundation.RowId),
) csc.MatrixError!csc.CscView {
    try matrix.validate();
    for (columns) |column| {
        if (column.dimension != matrix.num_rows) return error.DimensionMismatch;
        try column.validate();
    }
    return appendColumnsIntoAssumeValid(buffers, matrix, columns);
}

/// Trusted allocation-free column append into sufficiently large buffers.
pub fn appendColumnsIntoAssumeValid(
    buffers: *CscTransformBuffers,
    matrix: csc.CscMatrix,
    columns: []const sparse_vector.SparseVectorView(foundation.RowId),
) csc.MatrixError!csc.CscView {
    const new_num_cols = std.math.add(usize, matrix.num_cols, columns.len) catch return error.DimensionTooLarge;
    try csc.validateDimensions(matrix.num_rows, new_num_cols);
    var total_nnz = matrix.nnz();
    for (columns) |column|
        total_nnz = std.math.add(usize, total_nnz, column.nnz()) catch return error.DimensionTooLarge;
    try buffers.requireCapacity(new_num_cols, total_nnz, 0);

    @memcpy(buffers.col_starts[0 .. matrix.num_cols + 1], matrix.col_starts);
    @memcpy(buffers.row_indices[0..matrix.nnz()], matrix.row_indices);
    @memcpy(buffers.values[0..matrix.nnz()], matrix.values);
    var destination = matrix.nnz();
    for (columns, 0..) |column, offset| {
        @memcpy(buffers.row_indices[destination..][0..column.nnz()], column.indices);
        @memcpy(buffers.values[destination..][0..column.nnz()], column.values);
        destination += column.nnz();
        buffers.col_starts[matrix.num_cols + offset + 1] = destination;
    }
    return buffers.viewAssumeValid(matrix.num_rows, new_num_cols, total_nnz);
}

/// Validate and append a canonical CSC row block, returning a new owner.
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

/// Validate both blocks and append rows into reusable caller buffers.
pub fn appendRowsInto(
    buffers: *CscTransformBuffers,
    matrix: csc.CscMatrix,
    rows_to_add: []const sparse_vector.SparseVectorView(foundation.ColId),
) csc.MatrixError!csc.CscView {
    try matrix.validate();
    for (rows_to_add) |row| {
        if (row.dimension != matrix.num_cols) return error.DimensionMismatch;
        try row.validate();
    }
    return appendRowsIntoAssumeValid(buffers, matrix, rows_to_add);
}

/// Trusted allocation-free row append into sufficiently large buffers.
pub fn appendRowsIntoAssumeValid(
    buffers: *CscTransformBuffers,
    matrix: csc.CscMatrix,
    rows_to_add: []const sparse_vector.SparseVectorView(foundation.ColId),
) csc.MatrixError!csc.CscView {
    const new_num_rows = std.math.add(usize, matrix.num_rows, rows_to_add.len) catch return error.DimensionTooLarge;
    try csc.validateDimensions(new_num_rows, matrix.num_cols);
    var total_nnz = matrix.nnz();
    for (rows_to_add) |row|
        total_nnz = std.math.add(usize, total_nnz, row.nnz()) catch return error.DimensionTooLarge;
    try buffers.requireCapacity(matrix.num_cols, total_nnz, 0);

    const starts = buffers.col_starts[0 .. matrix.num_cols + 1];
    @memset(starts, 0);
    for (0..matrix.num_cols) |col|
        starts[col + 1] = matrix.col_starts[col + 1] - matrix.col_starts[col];
    for (rows_to_add) |row| for (row.indices) |col_id| {
        starts[col_id.toUsize() + 1] += 1;
    };
    for (0..matrix.num_cols) |col| starts[col + 1] += starts[col];

    for (0..matrix.num_cols) |col| {
        const begin = matrix.col_starts[col];
        const end = matrix.col_starts[col + 1];
        const count = end - begin;
        const destination = starts[col];
        @memcpy(buffers.row_indices[destination..][0..count], matrix.row_indices[begin..end]);
        @memcpy(buffers.values[destination..][0..count], matrix.values[begin..end]);
        starts[col] += count;
    }
    for (rows_to_add, 0..) |row, added_row| {
        const row_id = foundation.RowId.fromUsizeAssumeValid(matrix.num_rows + added_row);
        for (row.indices, row.values) |col_id, value| {
            const col = col_id.toUsize();
            const destination = starts[col];
            buffers.row_indices[destination] = row_id;
            buffers.values[destination] = value;
            starts[col] += 1;
        }
    }
    restoreStartsAfterCursorUse(starts);
    return buffers.viewAssumeValid(new_num_rows, matrix.num_cols, total_nnz);
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

/// Validate a CSR row block and append it below CSC into caller buffers.
pub fn appendRowsFromCsrInto(
    buffers: *CscTransformBuffers,
    matrix: csc.CscMatrix,
    row_starts: []const usize,
    col_indices: []const foundation.ColId,
    stream_values: []const f64,
) csc.MatrixError!csc.CscView {
    try matrix.validate();
    try validateCsrRowStreams(matrix.num_cols, row_starts, col_indices, stream_values);
    return appendRowsFromCsrIntoAssumeValid(buffers, matrix, row_starts, col_indices, stream_values);
}

/// Trusted CSR-to-CSC row append using prevalidated exact-size buffers.
pub fn appendRowsFromCsrIntoAssumeValid(
    buffers: *CscTransformBuffers,
    matrix: csc.CscMatrix,
    row_starts: []const usize,
    col_indices: []const foundation.ColId,
    stream_values: []const f64,
) csc.MatrixError!csc.CscView {
    const added_rows = row_starts.len - 1;
    const new_num_rows = std.math.add(usize, matrix.num_rows, added_rows) catch return error.DimensionTooLarge;
    try csc.validateDimensions(new_num_rows, matrix.num_cols);
    const total_nnz = std.math.add(usize, matrix.nnz(), stream_values.len) catch return error.DimensionTooLarge;
    try buffers.requireCapacity(matrix.num_cols, total_nnz, 0);

    const starts = buffers.col_starts[0 .. matrix.num_cols + 1];
    @memset(starts, 0);
    for (0..matrix.num_cols) |col|
        starts[col + 1] = matrix.col_starts[col + 1] - matrix.col_starts[col];
    for (col_indices) |col_id| starts[col_id.toUsize() + 1] += 1;
    for (0..matrix.num_cols) |col| starts[col + 1] += starts[col];

    for (0..matrix.num_cols) |col| {
        const begin = matrix.col_starts[col];
        const end = matrix.col_starts[col + 1];
        const count = end - begin;
        const destination = starts[col];
        @memcpy(buffers.row_indices[destination..][0..count], matrix.row_indices[begin..end]);
        @memcpy(buffers.values[destination..][0..count], matrix.values[begin..end]);
        starts[col] += count;
    }
    for (0..added_rows) |added_row| {
        const row_id = foundation.RowId.fromUsizeAssumeValid(matrix.num_rows + added_row);
        for (row_starts[added_row]..row_starts[added_row + 1]) |position| {
            const col = col_indices[position].toUsize();
            const destination = starts[col];
            buffers.row_indices[destination] = row_id;
            buffers.values[destination] = stream_values[position];
            starts[col] += 1;
        }
    }
    restoreStartsAfterCursorUse(starts);
    return buffers.viewAssumeValid(new_num_rows, matrix.num_cols, total_nnz);
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

/// Validate canonical CSR offsets, IDs, values, and requested row dimension.
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

/// Return an owning matrix with the checked, strictly increasing columns removed.
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

/// Delete checked columns into reusable output and index-scratch buffers.
pub fn deleteColumnsInto(buffers: *CscTransformBuffers, matrix: csc.CscMatrix, deleted: []const foundation.ColId) csc.MatrixError!csc.CscView {
    try matrix.validate();
    try validateDeleted(foundation.ColId, deleted, matrix.num_cols);
    const output_cols = matrix.num_cols - deleted.len;
    var output_nnz = matrix.nnz();
    for (deleted) |col_id| {
        const col = col_id.toUsize();
        output_nnz -= matrix.col_starts[col + 1] - matrix.col_starts[col];
    }
    try buffers.requireCapacity(output_cols, output_nnz, 0);

    var deleted_pos: usize = 0;
    var output_col: usize = 0;
    var destination: usize = 0;
    buffers.col_starts[0] = 0;
    for (0..matrix.num_cols) |col| {
        if (deleted_pos < deleted.len and deleted[deleted_pos].toUsize() == col) {
            deleted_pos += 1;
            continue;
        }
        const begin = matrix.col_starts[col];
        const end = matrix.col_starts[col + 1];
        const count = end - begin;
        @memcpy(buffers.row_indices[destination..][0..count], matrix.row_indices[begin..end]);
        @memcpy(buffers.values[destination..][0..count], matrix.values[begin..end]);
        destination += count;
        output_col += 1;
        buffers.col_starts[output_col] = destination;
    }
    return buffers.viewAssumeValid(matrix.num_rows, output_cols, output_nnz);
}

/// Return an owning matrix with checked rows removed and remaining IDs remapped.
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

/// Delete checked rows into reusable buffers, using scratch for the ID remap.
pub fn deleteRowsInto(buffers: *CscTransformBuffers, matrix: csc.CscMatrix, deleted: []const foundation.RowId) csc.MatrixError!csc.CscView {
    try matrix.validate();
    try validateDeleted(foundation.RowId, deleted, matrix.num_rows);
    try buffers.requireCapacity(matrix.num_cols, 0, matrix.num_rows);
    const missing = std.math.maxInt(usize);
    const row_map = buffers.index_scratch[0..matrix.num_rows];
    var deleted_pos: usize = 0;
    var next_row: usize = 0;
    for (0..matrix.num_rows) |row| {
        if (deleted_pos < deleted.len and deleted[deleted_pos].toUsize() == row) {
            row_map[row] = missing;
            deleted_pos += 1;
        } else {
            row_map[row] = next_row;
            next_row += 1;
        }
    }
    var output_nnz: usize = 0;
    for (matrix.row_indices) |row| if (row_map[row.toUsize()] != missing) {
        output_nnz += 1;
    };
    try buffers.requireCapacity(matrix.num_cols, output_nnz, matrix.num_rows);

    var destination: usize = 0;
    for (0..matrix.num_cols) |col| {
        buffers.col_starts[col] = destination;
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position| {
            const new_row = row_map[matrix.row_indices[position].toUsize()];
            if (new_row == missing) continue;
            buffers.row_indices[destination] = foundation.RowId.fromUsizeAssumeValid(new_row);
            buffers.values[destination] = matrix.values[position];
            destination += 1;
        }
    }
    buffers.col_starts[matrix.num_cols] = destination;
    return buffers.viewAssumeValid(next_row, matrix.num_cols, output_nnz);
}

/// Shift cursor-mutated CSC end offsets back into canonical start offsets.
fn restoreStartsAfterCursorUse(starts: []usize) void {
    var col = starts.len - 1;
    while (col > 0) {
        starts[col] = starts[col - 1];
        col -= 1;
    }
    starts[0] = 0;
}

/// Require a strictly increasing in-range deletion list.
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
    const columns = [_]sparse_vector.SparseVectorView(foundation.RowId){sparse_vector.SparseVectorView(foundation.RowId).initAssumeValid(2, &column_rows, &column_values)};
    var with_column = try appendColumns(std.testing.allocator, matrix, &columns);
    defer with_column.deinit(std.testing.allocator);
    try with_column.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 0, 0, 2 }, with_column.col_starts);

    var row_cols = [_]foundation.ColId{ try foundation.ColId.init(0), try foundation.ColId.init(2) };
    var row_values = [_]f64{ 4.0, 5.0 };
    const new_rows = [_]sparse_vector.SparseVectorView(foundation.ColId){sparse_vector.SparseVectorView(foundation.ColId).initAssumeValid(3, &row_cols, &row_values)};
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

test "edit Into APIs append and delete with reusable capacity" {
    const builder_module = @import("builder.zig");
    var builder = try builder_module.MatrixBuilder.init(2, 2);
    defer builder.deinit(std.testing.allocator);
    try builder.append(std.testing.allocator, try foundation.RowId.init(0), try foundation.ColId.init(0), 1.0);
    try builder.append(std.testing.allocator, try foundation.RowId.init(1), try foundation.ColId.init(1), 2.0);
    var base = try builder.freezeSortedLeanAssumeValid(std.testing.allocator, 0.0);
    defer base.deinit(std.testing.allocator);
    var buffers = try CscTransformBuffers.initCapacity(std.testing.allocator, 4, 8, 4);
    defer buffers.deinit(std.testing.allocator);

    const column_rows = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(1) };
    const column_values = [_]f64{ 3.0, 4.0 };
    const columns = [_]sparse_vector.SparseVectorView(foundation.RowId){
        sparse_vector.SparseVectorView(foundation.RowId).initAssumeValid(2, &column_rows, &column_values),
    };
    const with_column = try appendColumnsInto(&buffers, base, &columns);
    try with_column.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 1, 2, 4 }, with_column.col_starts);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 2.0, 3.0, 4.0 }, with_column.values);

    // The previous view is invalid after this call; use the owning base again.
    const row_starts = [_]usize{ 0, 2 };
    const row_cols = [_]foundation.ColId{ try foundation.ColId.init(0), try foundation.ColId.init(1) };
    const row_values = [_]f64{ 5.0, 6.0 };
    const with_row = try appendRowsFromCsrInto(&buffers, base, &row_starts, &row_cols, &row_values);
    try with_row.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 2, 4 }, with_row.col_starts);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 5.0, 2.0, 6.0 }, with_row.values);

    const deleted_cols = [_]foundation.ColId{try foundation.ColId.init(0)};
    const without_column = try deleteColumnsInto(&buffers, base, &deleted_cols);
    try without_column.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 1 }, without_column.col_starts);
    try std.testing.expectEqualSlices(f64, &.{2.0}, without_column.values);

    const deleted_rows = [_]foundation.RowId{try foundation.RowId.init(0)};
    const without_row = try deleteRowsInto(&buffers, base, &deleted_rows);
    try without_row.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 0, 1 }, without_row.col_starts);
    try std.testing.expectEqual(@as(usize, 0), without_row.row_indices[0].toUsize());

    var too_small = try CscTransformBuffers.initCapacity(std.testing.allocator, 1, 1, 1);
    defer too_small.deinit(std.testing.allocator);
    try std.testing.expectError(error.BufferTooSmall, appendColumnsInto(&too_small, base, &columns));
    try std.testing.expectError(error.BufferTooSmall, deleteRowsInto(&too_small, base, &deleted_rows));
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
    const wrong_columns = [_]sparse_vector.SparseVectorView(foundation.RowId){sparse_vector.SparseVectorView(foundation.RowId).initAssumeValid(3, &no_ids, &no_values)};
    try std.testing.expectError(error.DimensionMismatch, appendColumns(std.testing.allocator, matrix, &wrong_columns));
    const duplicate = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(0) };
    try std.testing.expectError(error.IndicesNotStrictlyIncreasing, deleteRows(std.testing.allocator, matrix, &duplicate));
}
