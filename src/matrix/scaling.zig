//! Row and column scaling for canonical sparse matrices.
//!
//! Positive finite factors preserve the sparsity pattern. Checked mutation is
//! transactional with respect to numerical errors: it preflights every result
//! before changing the matrix, so failure never leaves partially scaled data.

const std = @import("std");
const foundation = @import("foundation");
const csc = @import("csc.zig");
const memory = @import("memory.zig");

pub const ScalingView = struct {
    row: []const f64,
    col: []const f64,

    pub fn validate(self: @This(), num_rows: usize, num_cols: usize) csc.MatrixError!void {
        if (self.row.len != num_rows or self.col.len != num_cols) return error.DimensionMismatch;
        for (self.row) |factor| if (!validFactor(factor)) return error.InvalidScaling;
        for (self.col) |factor| if (!validFactor(factor)) return error.InvalidScaling;
    }
};

/// Applies A[i,j] := row[i] * A[i,j] * col[j].
pub fn apply(matrix: *csc.CscMatrix, factors: ScalingView) csc.MatrixError!void {
    try factors.validate(matrix.num_rows, matrix.num_cols);
    try preflight(matrix.*, factors, false);
    applyAssumeValid(matrix, factors);
}

/// Hot path for validated factors whose products are known representable.
pub fn applyAssumeValid(matrix: *csc.CscMatrix, factors: ScalingView) void {
    const ncol = matrix.num_cols;
    const starts = matrix.col_starts;
    const ri = matrix.row_indices;
    const vs = matrix.values;
    var col: usize = 0;
    while (col < ncol) : (col += 1) {
        const col_factor = factors.col[col];
        var position = starts[col];
        const end = starts[col + 1];
        while (position < end) : (position += 1) {
            const row: usize = @intFromEnum(ri[position]);
            vs[position] *= factors.row[row] * col_factor;
        }
    }
}

/// Reverses apply using division by the same positive factors.
pub fn remove(matrix: *csc.CscMatrix, factors: ScalingView) csc.MatrixError!void {
    try factors.validate(matrix.num_rows, matrix.num_cols);
    try preflight(matrix.*, factors, true);
    removeAssumeValid(matrix, factors);
}

pub fn removeAssumeValid(matrix: *csc.CscMatrix, factors: ScalingView) void {
    for (0..matrix.num_cols) |col| {
        const col_factor = factors.col[col];
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position| {
            const row = matrix.row_indices[position].toUsize();
            matrix.values[position] /= factors.row[row] * col_factor;
        }
    }
}

pub fn scaleColumn(matrix: *csc.CscMatrix, col: foundation.ColId, factor: f64) csc.MatrixError!void {
    const col_index = col.toUsize();
    if (col_index >= matrix.num_cols) return error.IndexOutOfBounds;
    if (!validFactor(factor)) return error.InvalidScaling;
    const begin = matrix.col_starts[col_index];
    const end = matrix.col_starts[col_index + 1];
    for (matrix.values[begin..end]) |value| {
        const result = value * factor;
        if (!std.math.isFinite(result) or result == 0.0) return error.NumericalOverflow;
    }
    for (matrix.values[begin..end]) |*value| value.* *= factor;
}

pub fn scaleRow(matrix: *csc.CscMatrix, row: foundation.RowId, factor: f64) csc.MatrixError!void {
    const row_index = row.toUsize();
    if (row_index >= matrix.num_rows) return error.IndexOutOfBounds;
    if (!validFactor(factor)) return error.InvalidScaling;
    for (matrix.row_indices, matrix.values) |row_id, value| {
        if (row_id.toUsize() != row_index) continue;
        const result = value * factor;
        if (!std.math.isFinite(result) or result == 0.0) return error.NumericalOverflow;
    }
    for (matrix.row_indices, matrix.values) |row_id, *value| {
        if (row_id.toUsize() == row_index) value.* *= factor;
    }
}

pub fn applyColumnFactors(matrix: *csc.CscMatrix, factors: []const f64) csc.MatrixError!void {
    if (factors.len != matrix.num_cols) return error.DimensionMismatch;
    // Validate and preflight directly, avoiding allocation of unit row factors.
    for (factors) |factor| if (!validFactor(factor)) return error.InvalidScaling;
    for (0..matrix.num_cols) |col| {
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position| {
            const result = matrix.values[position] * factors[col];
            if (!std.math.isFinite(result) or result == 0.0) return error.NumericalOverflow;
        }
    }
    for (0..matrix.num_cols) |col| {
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position| {
            matrix.values[position] *= factors[col];
        }
    }
}

pub fn applyRowFactors(matrix: *csc.CscMatrix, factors: []const f64) csc.MatrixError!void {
    if (factors.len != matrix.num_rows) return error.DimensionMismatch;
    for (factors) |factor| if (!validFactor(factor)) return error.InvalidScaling;
    for (matrix.row_indices, matrix.values) |row_id, value| {
        const result = value * factors[row_id.toUsize()];
        if (!std.math.isFinite(result) or result == 0.0) return error.NumericalOverflow;
    }
    for (matrix.row_indices, matrix.values) |row_id, *value| value.* *= factors[row_id.toUsize()];
}

/// Computes simple max-equilibration factors into caller-owned arrays.
///
/// Columns are normalized first. Row maxima are then measured after column
/// scaling and normalized. Empty rows/columns receive factor 1.
pub fn computeMaxEquilibration(matrix: csc.CscMatrix, row_factors: []f64, col_factors: []f64) csc.MatrixError!void {
    if (row_factors.len != matrix.num_rows or col_factors.len != matrix.num_cols)
        return error.DimensionMismatch;

    for (0..matrix.num_cols) |col| {
        var maximum: f64 = 0.0;
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position|
            maximum = @max(maximum, @abs(matrix.values[position]));
        col_factors[col] = if (maximum == 0.0) 1.0 else 1.0 / maximum;
    }

    memory.clearF64Fast(row_factors);
    for (0..matrix.num_cols) |col| {
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position| {
            const row = matrix.row_indices[position].toUsize();
            row_factors[row] = @max(row_factors[row], @abs(matrix.values[position]) * col_factors[col]);
        }
    }
    for (row_factors) |*maximum|
        maximum.* = if (maximum.* == 0.0) 1.0 else 1.0 / maximum.*;
}

/// HiGHS-style column factors: reciprocal maxima rounded to the nearest power
/// of two and clamped to +/- max_exponent.
pub fn computePowerOfTwoColumnFactors(matrix: csc.CscMatrix, max_exponent: u10, output: []f64) csc.MatrixError!void {
    if (output.len != matrix.num_cols) return error.DimensionMismatch;
    for (0..matrix.num_cols) |col| {
        var maximum: f64 = 0.0;
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position|
            maximum = @max(maximum, @abs(matrix.values[position]));
        output[col] = powerOfTwoReciprocal(maximum, max_exponent);
    }
}

/// Row counterpart. Caller may apply column factors first to match HiGHS's
/// sequential col-then-row equilibration policy.
pub fn computePowerOfTwoRowFactors(matrix: csc.CscMatrix, max_exponent: u10, output: []f64) csc.MatrixError!void {
    if (output.len != matrix.num_rows) return error.DimensionMismatch;
    memory.clearF64Fast(output);
    for (matrix.row_indices, matrix.values) |row_id, value| {
        const row = row_id.toUsize();
        output[row] = @max(output[row], @abs(value));
    }
    for (output) |*maximum| maximum.* = powerOfTwoReciprocal(maximum.*, max_exponent);
}

fn validFactor(factor: f64) bool {
    return std.math.isFinite(factor) and factor > 0.0;
}

fn powerOfTwoReciprocal(maximum: f64, max_exponent: u10) f64 {
    if (maximum == 0.0) return 1.0;
    const limit: f64 = @floatFromInt(max_exponent);
    const exponent = std.math.clamp(@round(-@log2(maximum)), -limit, limit);
    return @exp2(exponent);
}

fn preflight(matrix: csc.CscMatrix, factors: ScalingView, inverse: bool) csc.MatrixError!void {
    for (0..matrix.num_cols) |col| {
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position| {
            const row = matrix.row_indices[position].toUsize();
            const combined = factors.row[row] * factors.col[col];
            if (!std.math.isFinite(combined) or combined == 0.0) return error.NumericalOverflow;
            const result = if (inverse) matrix.values[position] / combined else matrix.values[position] * combined;
            if (!std.math.isFinite(result) or result == 0.0) return error.NumericalOverflow;
        }
    }
}

test "scaling apply and remove are reversible" {
    var starts = [_]usize{ 0, 2, 3 };
    var rows = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(1), try foundation.RowId.init(1) };
    var values = [_]f64{ 2.0, -3.0, 4.0 };
    var matrix = csc.CscMatrix.initBorrowedAssumeValid(2, 2, &starts, &rows, &values);
    const original = values;
    const factors: ScalingView = .{ .row = &.{ 2.0, 0.5 }, .col = &.{ 4.0, 0.25 } };
    try apply(&matrix, factors);
    try std.testing.expectEqualSlices(f64, &.{ 16.0, -6.0, 0.5 }, matrix.values);
    try remove(&matrix, factors);
    try std.testing.expectEqualSlices(f64, &original, matrix.values);
    try matrix.validate();
}

test "checked scaling failure does not partially modify values" {
    var starts = [_]usize{ 0, 2 };
    var rows = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(1) };
    var values = [_]f64{ 2.0, 3.0 };
    var matrix = csc.CscMatrix.initBorrowedAssumeValid(2, 1, &starts, &rows, &values);
    try std.testing.expectError(error.NumericalOverflow, apply(&matrix, .{ .row = &.{ 1.0, 1e308 }, .col = &.{2.0} }));
    try std.testing.expectEqualSlices(f64, &.{ 2.0, 3.0 }, matrix.values);
    try std.testing.expectError(error.InvalidScaling, apply(&matrix, .{ .row = &.{ 1.0, 0.0 }, .col = &.{1.0} }));
}

test "max equilibration handles empty rows and columns" {
    var starts = [_]usize{ 0, 2, 2 };
    var rows = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(2) };
    var values = [_]f64{ 2.0, -8.0 };
    const matrix = csc.CscMatrix.initBorrowedAssumeValid(3, 2, &starts, &rows, &values);
    var row_factors: [3]f64 = undefined;
    var col_factors: [2]f64 = undefined;
    try computeMaxEquilibration(matrix, &row_factors, &col_factors);
    try std.testing.expectEqualSlices(f64, &.{ 4.0, 1.0, 1.0 }, &row_factors);
    try std.testing.expectEqualSlices(f64, &.{ 0.125, 1.0 }, &col_factors);
}

test "single row and column scaling are checked and transactional" {
    var starts = [_]usize{ 0, 2, 3 };
    var rows = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(1), try foundation.RowId.init(1) };
    var values = [_]f64{ 2.0, 3.0, 4.0 };
    var matrix = csc.CscMatrix.initBorrowedAssumeValid(2, 2, &starts, &rows, &values);
    try scaleColumn(&matrix, try foundation.ColId.init(0), 2.0);
    try std.testing.expectEqualSlices(f64, &.{ 4.0, 6.0, 4.0 }, matrix.values);
    try scaleRow(&matrix, try foundation.RowId.init(1), 0.5);
    try std.testing.expectEqualSlices(f64, &.{ 4.0, 3.0, 2.0 }, matrix.values);
    try std.testing.expectError(error.NumericalOverflow, scaleColumn(&matrix, try foundation.ColId.init(0), 1e308));
    try std.testing.expectEqualSlices(f64, &.{ 4.0, 3.0, 2.0 }, matrix.values);
}

test "power-of-two scaling matches bounded nearest exponent policy" {
    var starts = [_]usize{ 0, 1, 2, 2 };
    var rows = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(1) };
    var values = [_]f64{ 3.0, 0.01 };
    var matrix = csc.CscMatrix.initBorrowedAssumeValid(3, 3, &starts, &rows, &values);
    var columns: [3]f64 = undefined;
    try computePowerOfTwoColumnFactors(matrix, 4, &columns);
    try std.testing.expectEqualSlices(f64, &.{ 0.25, 16.0, 1.0 }, &columns);
    try applyColumnFactors(&matrix, &columns);
    var row_factors: [3]f64 = undefined;
    try computePowerOfTwoRowFactors(matrix, 4, &row_factors);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 8.0, 1.0 }, &row_factors);
    try applyRowFactors(&matrix, &row_factors);
    try matrix.validate();
}
