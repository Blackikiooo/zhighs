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
