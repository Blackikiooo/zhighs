//! Read-only sparse matrix statistics and norms.
//!
//! Kernels accept caller-owned work buffers where necessary. This keeps memory
//! allocation out of iterative solver loops and makes scratch lifetime explicit.

const std = @import("std");
const foundation = @import("foundation");
const csc = @import("csc.zig");

pub const AbsoluteRange = struct { min: f64, max: f64 };

pub const ValueAssessment = struct {
    range: ?AbsoluteRange,
    small_count: usize,
    large_count: usize,
};

/// Exact canonical-storage equality. Canonical CSC has a unique ordering, so
/// storage equality is also mathematical equality under exact f64 semantics.
pub fn eql(lhs: csc.CscMatrix, rhs: csc.CscMatrix) bool {
    return lhs.num_rows == rhs.num_rows and lhs.num_cols == rhs.num_cols and
        std.mem.eql(usize, lhs.col_starts, rhs.col_starts) and
        std.mem.eql(foundation.RowId, lhs.row_indices, rhs.row_indices) and
        std.mem.eql(f64, lhs.values, rhs.values);
}

/// Minimum and maximum absolute nonzero value, or null for an empty matrix.
pub fn absoluteRange(matrix: csc.CscMatrix) ?AbsoluteRange {
    if (matrix.values.len == 0) return null;
    var result: AbsoluteRange = .{ .min = @abs(matrix.values[0]), .max = @abs(matrix.values[0]) };
    for (matrix.values[1..]) |value| {
        const magnitude = @abs(value);
        result.min = @min(result.min, magnitude);
        result.max = @max(result.max, magnitude);
    }
    return result;
}

/// Counts values at or below small_limit and at or above large_limit.
pub fn assessValues(matrix: csc.CscMatrix, small_limit: f64, large_limit: f64) csc.MatrixError!ValueAssessment {
    if (!std.math.isFinite(small_limit) or !std.math.isFinite(large_limit) or
        small_limit < 0.0 or large_limit < small_limit) return error.InvalidTolerance;
    var small_count: usize = 0;
    var large_count: usize = 0;
    for (matrix.values) |value| {
        const magnitude = @abs(value);
        if (magnitude <= small_limit) small_count += 1;
        if (magnitude >= large_limit) large_count += 1;
    }
    return .{ .range = absoluteRange(matrix), .small_count = small_count, .large_count = large_count };
}

pub fn hasLargeValue(matrix: csc.CscMatrix, large_limit: f64) csc.MatrixError!bool {
    if (!std.math.isFinite(large_limit) or large_limit < 0.0) return error.InvalidTolerance;
    for (matrix.values) |value| if (@abs(value) >= large_limit) return true;
    return false;
}

/// Largest absolute matrix entry. Returns zero for an empty matrix.
pub fn maxAbs(matrix: csc.CscMatrix) f64 {
    var result: f64 = 0.0;
    for (matrix.values) |value| result = @max(result, @abs(value));
    return result;
}

/// Writes the absolute sum of every column.
pub fn columnOneNorms(matrix: csc.CscMatrix, output: []f64) csc.MatrixError!void {
    if (output.len != matrix.num_cols) return error.DimensionMismatch;
    columnOneNormsAssumeValid(matrix, output);
}

pub fn columnOneNormsAssumeValid(matrix: csc.CscMatrix, output: []f64) void {
    for (0..matrix.num_cols) |col| {
        var sum: f64 = 0.0;
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position|
            sum += @abs(matrix.values[position]);
        output[col] = sum;
    }
}

/// Matrix one-norm: maximum absolute column sum.
pub fn oneNorm(matrix: csc.CscMatrix) f64 {
    var result: f64 = 0.0;
    for (0..matrix.num_cols) |col| {
        var sum: f64 = 0.0;
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position|
            sum += @abs(matrix.values[position]);
        result = @max(result, sum);
    }
    return result;
}

/// Writes absolute row sums into caller-provided scratch/output storage.
pub fn rowOneNorms(matrix: csc.CscMatrix, output: []f64) csc.MatrixError!void {
    if (output.len != matrix.num_rows) return error.DimensionMismatch;
    rowOneNormsAssumeValid(matrix, output);
}

pub fn rowOneNormsAssumeValid(matrix: csc.CscMatrix, output: []f64) void {
    @memset(output, 0.0);
    for (matrix.row_indices, matrix.values) |row, value|
        output[row.toUsize()] += @abs(value);
}

/// Matrix infinity-norm using caller-owned row scratch.
pub fn infinityNorm(matrix: csc.CscMatrix, row_scratch: []f64) csc.MatrixError!f64 {
    try rowOneNorms(matrix, row_scratch);
    var result: f64 = 0.0;
    for (row_scratch) |sum| result = @max(result, sum);
    return result;
}

/// Numerically stable Frobenius norm.
///
/// The scaled sum-of-squares algorithm avoids overflow from value*value and
/// underflow when entries have very different magnitudes.
pub fn frobeniusNorm(matrix: csc.CscMatrix) f64 {
    var scale: f64 = 0.0;
    var scaled_squares: f64 = 1.0;
    for (matrix.values) |value| {
        const magnitude = @abs(value);
        if (magnitude == 0.0) continue;
        if (scale < magnitude) {
            const ratio = scale / magnitude;
            scaled_squares = 1.0 + scaled_squares * ratio * ratio;
            scale = magnitude;
        } else {
            const ratio = magnitude / scale;
            scaled_squares += ratio * ratio;
        }
    }
    return if (scale == 0.0) 0.0 else scale * @sqrt(scaled_squares);
}

/// Computes y += alpha * A * x without clearing y.
pub fn addProduct(matrix: csc.CscMatrix, alpha: f64, x: []const f64, y: []f64) csc.MatrixError!void {
    if (x.len != matrix.num_cols or y.len != matrix.num_rows) return error.DimensionMismatch;
    if (!std.math.isFinite(alpha)) return error.NonFiniteValue;
    addProductAssumeValid(matrix, alpha, x, y);
}

pub fn addProductAssumeValid(matrix: csc.CscMatrix, alpha: f64, x: []const f64, y: []f64) void {
    if (alpha == 0.0) return;
    for (0..matrix.num_cols) |col| {
        const multiplier = alpha * x[col];
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position|
            y[matrix.row_indices[position].toUsize()] += multiplier * matrix.values[position];
    }
}

pub fn addProductSkippingZeros(matrix: csc.CscMatrix, alpha: f64, x: []const f64, y: []f64) csc.MatrixError!void {
    if (x.len != matrix.num_cols or y.len != matrix.num_rows) return error.DimensionMismatch;
    if (!std.math.isFinite(alpha)) return error.NonFiniteValue;
    addProductSkippingZerosAssumeValid(matrix, alpha, x, y);
}

pub fn addProductSkippingZerosAssumeValid(matrix: csc.CscMatrix, alpha: f64, x: []const f64, y: []f64) void {
    if (alpha == 0.0) return;
    for (0..matrix.num_cols) |col| {
        const multiplier = alpha * x[col];
        if (multiplier == 0.0) continue;
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position|
            y[matrix.row_indices[position].toUsize()] += multiplier * matrix.values[position];
    }
}

/// Computes y += alpha * transpose(A) * x without clearing y.
pub fn addTransposeProduct(matrix: csc.CscMatrix, alpha: f64, x: []const f64, y: []f64) csc.MatrixError!void {
    if (x.len != matrix.num_rows or y.len != matrix.num_cols) return error.DimensionMismatch;
    if (!std.math.isFinite(alpha)) return error.NonFiniteValue;
    addTransposeProductAssumeValid(matrix, alpha, x, y);
}

pub fn addTransposeProductAssumeValid(matrix: csc.CscMatrix, alpha: f64, x: []const f64, y: []f64) void {
    if (alpha == 0.0) return;
    for (0..matrix.num_cols) |col| {
        var sum: f64 = 0.0;
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position|
            sum += matrix.values[position] * x[matrix.row_indices[position].toUsize()];
        y[col] += alpha * sum;
    }
}

/// High-precision y = A*x using caller-owned HCD row scratch.
pub fn multiplyHighPrecision(matrix: csc.CscMatrix, x: []const f64, y: []f64, scratch: []foundation.HCD) csc.MatrixError!void {
    if (x.len != matrix.num_cols or y.len != matrix.num_rows or scratch.len != matrix.num_rows)
        return error.DimensionMismatch;
    multiplyHighPrecisionAssumeValid(matrix, x, y, scratch);
}

pub fn multiplyHighPrecisionAssumeValid(matrix: csc.CscMatrix, x: []const f64, y: []f64, scratch: []foundation.HCD) void {
    for (scratch) |*sum| sum.* = foundation.HCD.initWithHD(0.0);
    for (0..matrix.num_cols) |col| {
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position| {
            const row = matrix.row_indices[position].toUsize();
            scratch[row] = scratch[row].addHCD(foundation.HCD.twoProduct(matrix.values[position], x[col]));
        }
    }
    for (y, scratch) |*result, sum| result.* = sum.toHD();
}

/// HiGHS-compatible fast HCD accumulation: residuals are accumulated without
/// per-entry renormalization. Faster, but the robust variant above is preferred
/// when a row may contain very many terms or extreme dynamic range.
pub fn multiplyHighPrecisionFastAssumeValid(matrix: csc.CscMatrix, x: []const f64, y: []f64, scratch: []foundation.HCD) void {
    for (scratch) |*sum| sum.* = foundation.HCD.initWithHD(0.0);
    for (0..matrix.num_cols) |col| {
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position| {
            const row = matrix.row_indices[position].toUsize();
            scratch[row] = scratch[row].addHCDFast(foundation.HCD.twoProduct(matrix.values[position], x[col]));
        }
    }
    for (y, scratch) |*result, sum| result.* = sum.toHD();
}

/// Matches HiGHS productQuad semantics: multiplication rounds to f64, while
/// accumulation uses an HCD residual. This is faster than exact two-product.
pub fn multiplyCompensatedAssumeValid(matrix: csc.CscMatrix, x: []const f64, y: []f64, scratch: []foundation.HCD) void {
    for (scratch) |*sum| sum.* = foundation.HCD.initWithHD(0.0);
    for (0..matrix.num_cols) |col| {
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position| {
            const row = matrix.row_indices[position].toUsize();
            scratch[row] = scratch[row].addHDFast(matrix.values[position] * x[col]);
        }
    }
    for (y, scratch) |*result, sum| result.* = sum.toHD();
}

/// High-precision y = transpose(A)*x; CSC needs only one HCD accumulator.
pub fn transposeMultiplyHighPrecision(matrix: csc.CscMatrix, x: []const f64, y: []f64) csc.MatrixError!void {
    if (x.len != matrix.num_rows or y.len != matrix.num_cols) return error.DimensionMismatch;
    transposeMultiplyHighPrecisionAssumeValid(matrix, x, y);
}

pub fn transposeMultiplyHighPrecisionAssumeValid(matrix: csc.CscMatrix, x: []const f64, y: []f64) void {
    for (0..matrix.num_cols) |col| {
        var sum = foundation.HCD.initWithHD(0.0);
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position| {
            const row = matrix.row_indices[position].toUsize();
            sum = sum.addHCD(foundation.HCD.twoProduct(matrix.values[position], x[row]));
        }
        y[col] = sum.toHD();
    }
}

pub fn transposeMultiplyHighPrecisionFastAssumeValid(matrix: csc.CscMatrix, x: []const f64, y: []f64) void {
    for (0..matrix.num_cols) |col| {
        var sum = foundation.HCD.initWithHD(0.0);
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position| {
            const row = matrix.row_indices[position].toUsize();
            sum = sum.addHCDFast(foundation.HCD.twoProduct(matrix.values[position], x[row]));
        }
        y[col] = sum.toHD();
    }
}

pub fn transposeMultiplyCompensatedAssumeValid(matrix: csc.CscMatrix, x: []const f64, y: []f64) void {
    for (0..matrix.num_cols) |col| {
        var sum = foundation.HCD.initWithHD(0.0);
        for (matrix.col_starts[col]..matrix.col_starts[col + 1]) |position| {
            const row = matrix.row_indices[position].toUsize();
            sum = sum.addHDFast(matrix.values[position] * x[row]);
        }
        y[col] = sum.toHD();
    }
}

test "matrix norms and absolute statistics" {
    var starts = [_]usize{ 0, 2, 3 };
    var rows = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(1), try foundation.RowId.init(1) };
    var values = [_]f64{ -3.0, 4.0, -12.0 };
    const matrix: csc.CscMatrix = .{ .num_rows = 2, .num_cols = 2, .col_starts = &starts, .row_indices = &rows, .values = &values };

    try std.testing.expectEqual(@as(f64, 12.0), maxAbs(matrix));
    try std.testing.expectEqual(@as(f64, 12.0), oneNorm(matrix));
    var columns: [2]f64 = undefined;
    try columnOneNorms(matrix, &columns);
    try std.testing.expectEqualSlices(f64, &.{ 7.0, 12.0 }, &columns);
    var rows_scratch: [2]f64 = undefined;
    try std.testing.expectEqual(@as(f64, 16.0), try infinityNorm(matrix, &rows_scratch));
    try std.testing.expectEqualSlices(f64, &.{ 3.0, 16.0 }, &rows_scratch);
    try std.testing.expectEqual(@as(f64, 13.0), frobeniusNorm(matrix));
}

test "stable Frobenius norm avoids intermediate overflow" {
    var starts = [_]usize{ 0, 2 };
    var rows = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(1) };
    var values = [_]f64{ 1e200, 1e200 };
    const matrix: csc.CscMatrix = .{ .num_rows = 2, .num_cols = 1, .col_starts = &starts, .row_indices = &rows, .values = &values };
    const norm = frobeniusNorm(matrix);
    try std.testing.expect(std.math.isFinite(norm));
    try std.testing.expectApproxEqRel(@sqrt(2.0) * 1e200, norm, 1e-15);
}

test "norm output dimensions are checked" {
    var matrix = try csc.CscMatrix.initZero(std.testing.allocator, 2, 3);
    defer matrix.deinit(std.testing.allocator);
    var wrong: [1]f64 = undefined;
    try std.testing.expectError(error.DimensionMismatch, columnOneNorms(matrix, &wrong));
    try std.testing.expectError(error.DimensionMismatch, rowOneNorms(matrix, &wrong));
    try std.testing.expectError(error.DimensionMismatch, infinityNorm(matrix, &wrong));
}

test "range assessment equality and alpha products" {
    const RowId = foundation.RowId;
    var starts = [_]usize{ 0, 2, 3 };
    var rows = [_]RowId{ try RowId.init(0), try RowId.init(1), try RowId.init(1) };
    var values = [_]f64{ -0.25, 4.0, 20.0 };
    const matrix: csc.CscMatrix = .{ .num_rows = 2, .num_cols = 2, .col_starts = &starts, .row_indices = &rows, .values = &values };
    try std.testing.expect(eql(matrix, matrix));
    const range = absoluteRange(matrix).?;
    try std.testing.expectEqual(@as(f64, 0.25), range.min);
    try std.testing.expectEqual(@as(f64, 20.0), range.max);
    const assessment = try assessValues(matrix, 0.5, 10.0);
    try std.testing.expectEqual(@as(usize, 1), assessment.small_count);
    try std.testing.expectEqual(@as(usize, 1), assessment.large_count);
    try std.testing.expect(try hasLargeValue(matrix, 20.0));

    var y = [_]f64{ 1.0, 2.0 };
    try addProduct(matrix, 2.0, &.{ 2.0, 1.0 }, &y);
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 58.0 }, &y);
    var yt = [_]f64{ 1.0, 2.0 };
    try addTransposeProduct(matrix, 0.5, &.{ 2.0, 3.0 }, &yt);
    try std.testing.expectEqualSlices(f64, &.{ 6.75, 32.0 }, &yt);
}

test "high precision products retain small residual" {
    const RowId = foundation.RowId;
    var starts = [_]usize{ 0, 1, 2, 3 };
    var rows = [_]RowId{ try RowId.init(0), try RowId.init(0), try RowId.init(0) };
    var values = [_]f64{ 1e16, 1.0, -1e16 };
    const matrix: csc.CscMatrix = .{ .num_rows = 1, .num_cols = 3, .col_starts = &starts, .row_indices = &rows, .values = &values };
    var y: [1]f64 = undefined;
    var scratch: [1]foundation.HCD = undefined;
    try multiplyHighPrecision(matrix, &.{ 1.0, 1.0, 1.0 }, &y, &scratch);
    try std.testing.expectEqual(@as(f64, 1.0), y[0]);
}
