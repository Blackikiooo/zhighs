const std = @import("std");

pub const HighsDouble = f64;

pub const HighsCDouble = struct {
    high: HighsDouble,
    low: HighsDouble,

    const Self = @This();

    const default = Self{
        .high = 0.0,
        .low = 0.0,
    };

    pub fn twoSum(a: HighsDouble, b: HighsDouble) Self {
        @setFloatMode(.strict);
        const x = a + b;
        const z = x - a;
        const y = (a - (x - z)) + (b - z);

        return Self{
            .high = x,
            .low = y,
        };
    }

    pub fn split(a: HighsDouble) Self {
        const factor = (1 << 27) + 1; // 2^27 + 1
        const t = factor * a;
        const high = t - (t - a);
        const low = a - high;

        return Self{
            .high = high,
            .low = low,
        };
    }

    pub fn init(high: HighsDouble) Self {
        return Self{
            .high = high,
            .low = 0.0,
        };
    }

    pub fn toHighsDouble(self: Self) HighsDouble {
        return self.high + self.low;
    }

    pub fn addDoubleAssign(self: *Self, o: HighsDouble) void {
        const sum = HighsCDouble.twoSum(self.high, o);
        self.high = sum.high;
        self.low += sum.low;
    }

    pub fn addCDoubleAssign(self: *Self, o: HighsCDouble) void {
        const sum = HighsCDouble.twoSum(self.high, o.high);
        self.high = sum.high;
        self.low += sum.low + o.low;
    }

    pub fn cmp(self: Self, o: Self) std.math.Order {
        const a = self.toHighsDouble(); // mind that this may introduce rounding errors.
        const b = o.toHighsDouble();
        if (a < b) return .lt;
        if (a > b) return .gt;
        return .eq;
    }
};

test "two-sum-test" {
    const res = HighsCDouble.twoSum(1e16, 1.0);
    const a = res.high;
    const b = res.low;
    const expected = a + b;

    try std.testing.expectEqual(a, 1e16);
    try std.testing.expectEqual(b, 1.0);
    try std.testing.expect(1e16 == expected);
}

test "cmp-test" {
    const a = HighsCDouble.init(1.0);
    const b = HighsCDouble.init(2.0);
    try std.testing.expectEqual(a.cmp(b), .lt);
    try std.testing.expectEqual(b.cmp(a), .gt);
    try std.testing.expectEqual(a.cmp(a), .eq);

    var d = HighsCDouble.init(1e16);
    d.addDoubleAssign(1);
    try std.testing.expectEqual(HighsCDouble.init(1e16).cmp(d), .eq);
}
