const std = @import("std");
const builtin = @import("builtin");

/// Double precision.
pub const HD = f64;

const has_fma = switch (builtin.cpu.arch) {
    .x86, .x86_64 => builtin.cpu.has(.x86, .fma),
    else => false,
};

/// Double-double precision.
pub const HCD = struct {
    high: HD,
    low: HD,

    const Self = @This();

    const default = Self{
        .high = 0.0,
        .low = 0.0,
    };

    /// performs an exact transformation such that x + y = a + b
    /// and x = double(a + b). The operation uses 6 flops (addition/subtraction).
    pub fn twoSum(a: HD, b: HD) Self {
        @setFloatMode(.strict);

        const x = a + b;
        const z = x - a;
        const y = (a - (x - z)) + (b - z);

        return Self{
            .high = x,
            .low = y,
        };
    }

    /// splits a 53 bit double precision number into two 26 bit parts
    /// such that x + y = a holds exactly
    pub fn split(a: HD) Self {
        @setFloatMode(.strict);

        const factor: comptime_int = (1 << 27) + 1; // 2^27 + 1
        const t = factor * a;
        const high = t - (t - a);
        const err = a - high;

        return Self{
            .high = high,
            .low = err,
        };
    }

    /// performs an exact transformation such that x + y = a * b
    /// and x = double(a * b). The operation uses 10 flops for
    /// addition/subtraction and 7 flops for multiplication.
    pub fn twoProduct(a: HD, b: HD) Self {
        @setFloatMode(.strict);

        const x = a * b;

        const err = if (comptime has_fma) blk: {
            break :blk @mulAdd(HD, a, b, -x);
        } else blk: {
            const a_split = split(a);
            const b_split = split(b);
            break :blk (a_split.low * b_split.low) - (((x - a_split.high * b_split.high) - a_split.low * b_split.high) - a_split.high * b_split.low);
        };

        return Self{
            .high = x,
            .low = err,
        };
    }

    pub fn initWithHD(value: HD) Self {
        return Self{
            .high = value,
            .low = 0.0,
        };
    }

    /// If you need to convert this to a double, use this function,
    /// but be aware that this may introduce rounding errors again.
    pub inline fn toHD(self: Self) HD {
        return self.high + self.low;
    }

    /// This will make `.low` more and more bigger, and may introduce more and more rounding errors,
    /// so you will take account of the renorm outsiede because invoke `twoSum` is expensive.
    pub fn addHDAssign(self: *Self, o: HD) void {
        const sum = HCD.twoSum(self.high, o);
        self.high = sum.high;
        self.low += sum.low;
    }
    /// The same as `addHDAssign`, but the parameter is a `HCD`.
    pub fn addHCDAssign(self: *Self, o: HCD) void {
        const sum = HCD.twoSum(self.high, o.high);
        self.high = sum.high;
        self.low += sum.low + o.low;
    }

    pub fn minusHDAssign(self: *Self, o: HD) void {
        const sum = HCD.twoSum(self.high, -o);
        self.high = sum.high;
        self.low += sum.low;
    }

    pub fn minusHCDAssign(self: *Self, o: HCD) void {
        const sum = HCD.twoSum(self.high, -o.high);
        self.high = sum.high;
        self.low += sum.low - o.low;
    }

    pub fn cmp(self: Self, o: Self) std.math.Order {
        const a = self.toHD(); // mind that this may introduce rounding errors.
        const b = o.toHD();
        if (a < b) return .lt;
        if (a > b) return .gt;
        return .eq;
    }
};

test "two-sum-test" {
    const res = HCD.twoSum(1e16, 1.0);
    const a = res.high;
    const b = res.low;
    const expected = a + b;

    try std.testing.expectEqual(a, 1e16);
    try std.testing.expectEqual(b, 1.0);
    try std.testing.expect(1e16 == expected);
}

test "cmp-test" {
    const a = HCD.initWithHD(1.0);
    const b = HCD.initWithHD(2.0);
    try std.testing.expectEqual(a.cmp(b), .lt);
    try std.testing.expectEqual(b.cmp(a), .gt);
    try std.testing.expectEqual(a.cmp(a), .eq);

    var d = HCD.initWithHD(1e16);
    d.addHDAssign(1);
    try std.testing.expectEqual(HCD.initWithHD(1e16).cmp(d), .eq);
}

test "two-product" {
    const res1 = HCD.twoProduct(1e16, 1.0);
    try std.testing.expectEqual(res1.toHD(), 1e16);

    const res2 = HCD.twoProduct(1e16, 0.3);
    try std.testing.expectEqual(res2.toHD(), 3e15);
}

test "split" {
    const res = HCD.split(1e16);
    try std.testing.expect(res.toHD() == 1e16);
}
