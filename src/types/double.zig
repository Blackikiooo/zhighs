const std = @import("std");

pub const HighsDouble = f64;

pub const HighsCDouble = struct {
    high: HighsDouble,
    low: HighsDouble,

    const Self = @This();
};

fn _twoSum(x: *HighsDouble, y: *HighsDouble, a: HighsDouble, b: HighsDouble) void {
    x.* = a + b;
    const z = x.* - a;
    y.* = (a - (x.* - z)) + (b - z);
}

test "two-sum-test" {
    var a: HighsDouble = 0;
    var b: HighsDouble = 0;

    _twoSum(&a, &b, 1e16, 1.0);
    const expected = a + b;

    try std.testing.expectEqual(a, 1e16);
    try std.testing.expectEqual(b, 1.0);
    try std.testing.expect(1e16 == expected);
}
