//! Exact signatures and collision-safe comparison for duplicate LP rows.

const std = @import("std");

fn canonicalBits(value: f64) u64 {
    return if (value == 0.0) 0 else @bitCast(value);
}

fn mix(seed: u64, value: u64) u64 {
    var bits = value;
    return std.hash.Wyhash.hash(seed, std.mem.asBytes(&bits));
}

/// Hash a canonical row. Hash equality is only a candidate filter; callers
/// must use `exactlyEqual` before removing a row.
pub fn hash(
    lower: f64,
    upper: f64,
    columns: []const u32,
    values: []const f64,
) u64 {
    std.debug.assert(columns.len == values.len);
    var result = mix(0, canonicalBits(lower));
    result = mix(result, canonicalBits(upper));
    result = mix(result, columns.len);
    for (columns, values) |column, value| {
        result = mix(result, column);
        result = mix(result, canonicalBits(value));
    }
    return result;
}

/// Confirm exact equality after a signature match. Signed zero bounds are
/// equivalent; matrix storage cannot contain explicit zero coefficients.
pub fn exactlyEqual(
    left_lower: f64,
    left_upper: f64,
    left_columns: []const u32,
    left_values: []const f64,
    right_lower: f64,
    right_upper: f64,
    right_columns: []const u32,
    right_values: []const f64,
) bool {
    if (left_lower != right_lower or left_upper != right_upper or
        !std.mem.eql(u32, left_columns, right_columns) or left_values.len != right_values.len)
        return false;
    for (left_values, right_values) |left, right| {
        if (left != right) return false;
    }
    return true;
}

/// Hash only the column pattern. Rows with equal hashes may still differ in
/// values and bounds; callers must compare the pattern and values explicitly.
pub fn patternHash(columns: []const u32) u64 {
    var result = mix(0, columns.len);
    for (columns) |column| result = mix(result, column);
    return result;
}

/// Return `scale` when `right = scale * left` exactly for every coefficient.
/// The rows must share the same ordered column pattern and be nonempty.
pub fn parallelScale(
    left_columns: []const u32,
    left_values: []const f64,
    right_columns: []const u32,
    right_values: []const f64,
) ?f64 {
    if (left_columns.len == 0 or left_columns.len != right_columns.len or left_values.len != right_values.len)
        return null;
    if (!std.mem.eql(u32, left_columns, right_columns)) return null;
    const scale = right_values[0] / left_values[0];
    if (scale == 0.0 or !std.math.isFinite(scale)) return null;
    for (left_values, right_values) |left, right| {
        if (!sameBits(scale * left, right)) return null;
    }
    return scale;
}

/// Confirm that scaling the left interval by `scale` reproduces the right
/// interval exactly. A negative scale swaps the interval sides.
pub fn boundsMatchScale(
    scale: f64,
    left_lower: f64,
    left_upper: f64,
    right_lower: f64,
    right_upper: f64,
) bool {
    if (!std.math.isFinite(scale) or scale == 0.0) return false;
    if (scale > 0.0)
        return scaledBoundEqual(scale, left_lower, right_lower) and
            scaledBoundEqual(scale, left_upper, right_upper);
    return scaledBoundEqual(scale, left_upper, right_lower) and
        scaledBoundEqual(scale, left_lower, right_upper);
}

fn sameBits(left: f64, right: f64) bool {
    return canonicalBits(left) == canonicalBits(right);
}

fn scaledBoundEqual(scale: f64, bound: f64, target: f64) bool {
    if (std.math.isNan(bound) or std.math.isNan(target)) return false;
    return sameBits(scale * bound, target);
}

test "column patterns and exact scaled rows" {
    const columns = [_]u32{ 0, 3 };
    const left = [_]f64{ 1.0, -2.0 };
    const right = [_]f64{ 2.0, -4.0 };
    try std.testing.expectEqual(patternHash(&columns), patternHash(&columns));
    try std.testing.expectEqual(@as(?f64, 2.0), parallelScale(&columns, &left, &columns, &right));
    try std.testing.expect(boundsMatchScale(2.0, -1.0, 3.0, -2.0, 6.0));
    try std.testing.expect(!boundsMatchScale(2.0, -1.0, 3.0, -2.0, 5.0));
    try std.testing.expect(boundsMatchScale(-2.0, -1.0, 3.0, -6.0, 2.0));
}

test "near-proportional rows are not removed" {
    const columns = [_]u32{ 0, 1 };
    const left = [_]f64{ 1.0, 3.0 };
    const right = [_]f64{ 2.0, 6.000000000000001 };
    try std.testing.expectEqual(@as(?f64, null), parallelScale(&columns, &left, &columns, &right));
}
test "duplicate row signatures preserve exact row semantics" {
    const columns = [_]u32{ 0, 3 };
    const values = [_]f64{ 2.0, -4.0 };
    try std.testing.expectEqual(
        hash(-0.0, 5.0, &columns, &values),
        hash(0.0, 5.0, &columns, &values),
    );
    try std.testing.expect(exactlyEqual(
        -0.0,
        5.0,
        &columns,
        &values,
        0.0,
        5.0,
        &columns,
        &values,
    ));
    try std.testing.expect(!exactlyEqual(
        0.0,
        5.0,
        &columns,
        &values,
        0.0,
        5.0,
        &columns,
        &.{ 2.0, -3.0 },
    ));
}
