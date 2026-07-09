const std = @import("std");
const builtin = @import("builtin");

/// Double precision.
pub const HD = f64;

const has_fma = switch (builtin.cpu.arch) {
    // x86 系列
    .x86, .x86_64 => builtin.cpu.has(.x86, .fma),
    // AArch64 / ARM32
    .aarch64 => builtin.cpu.has(.aarch64, .fp_armv8),
    .arm => builtin.cpu.has(.arm, .fp_armv8),
    // RISC-V 带浮点扩展就有FMA
    .riscv64, .riscv32 => builtin.cpu.has(.riscv, .d),
    // PowerPC
    .powerpc, .powerpcle, .powerpc64, .powerpc64le => true,
    // 龙芯
    .loongarch64 => true,
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
    /// and x = double(a * b).
    pub inline fn twoProduct(a: HD, b: HD) Self {
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

    inline fn quickTwoSum(a: HD, b: HD) Self {
        @setFloatMode(.strict);

        const high = a + b;
        const low = b - (high - a);

        return Self{
            .high = high,
            .low = low,
        };
    }

    inline fn quickThreeSum(a: HD, b: HD, c: HD) Self {
        @setFloatMode(.strict);

        const sum = a + b;
        const err = b - (sum - a);
        return quickTwoSum(sum, err + c);
    }

    pub inline fn initWithHD(value: HD) Self {
        return Self{
            .high = value,
            .low = 0.0,
        };
    }

    pub inline fn init(high: HD, low: HD) Self {
        return Self{
            .high = high,
            .low = low,
        };
    }

    /// use const x = self is clone semantics, but use clone will more clear.
    pub inline fn clone(self: Self) Self {
        return Self{
            .high = self.high,
            .low = self.low,
        };
    }

    /// If you need to convert this to a double, use this function,
    /// but be aware that this may introduce rounding errors again.
    pub inline fn toHD(self: Self) HD {
        return self.high + self.low;
    }

    pub fn addHD(self: Self, o: HD) HCD {
        var res: HCD = self.clone();
        const sum = twoSum(res.high, o);
        res.high = sum.high;
        res.low += sum.low;
        return res;
    }

    pub fn addHCD(self: Self, o: HCD) HCD {
        var res: HCD = self.clone();
        const sum = twoSum(res.high, o.high);
        res.high = sum.high;
        res.low += sum.low + o.low;
        return res;
    }

    pub fn minusHD(self: Self, o: HD) HCD {
        var res: HCD = self.clone();
        const sum = twoSum(res.high, -o);
        res.high = sum.high;
        res.low += sum.low;
        return res;
    }

    pub fn minusHCD(self: Self, o: HCD) HCD {
        var res: HCD = self.clone();
        const sum = twoSum(res.high, -o.high);
        res.high = sum.high;
        res.low += sum.low - o.low;
        return res;
    }

    pub inline fn multiplyHD(self: Self, o: HD) HCD {
        const product = twoProduct(self.high, o);
        return quickTwoSum(product.high, product.low + self.low * o);
    }

    pub inline fn multiplyHCD(self: Self, o: HCD) HCD {
        if (o.low == 0.0) {
            return self.multiplyHD(o.high);
        }

        const product = twoProduct(self.high, o.high);
        const cross_products = product.low + self.high * o.low + self.low * o.high;
        return quickTwoSum(product.high, cross_products);
    }

    pub inline fn divideHD(self: Self, o: HD) HCD {
        @setFloatMode(.strict);

        if (o == 0.0) {
            @branchHint(.unlikely);
            return Self.initWithHD(self.toHD() / o);
        }

        const q1 = self.high / o;
        const product = Self.twoProduct(q1, o);
        const q2 = (((self.high - product.high) - product.low) + self.low) / o;

        return quickTwoSum(q1, q2);
    }

    pub inline fn divideHCD(self: Self, o: HCD) HCD {
        @setFloatMode(.strict);

        if (o.high == 0.0 and o.low == 0.0) {
            @branchHint(.unlikely);
            return Self.initWithHD(self.toHD() / o.toHD());
        }
        if (o.low == 0.0) {
            return self.divideHD(o.high);
        }

        const q1 = self.high / o.high;
        const p1 = Self.twoProduct(q1, o.high);
        const p2 = q1 * o.low;
        const r1 = (((self.high - p1.high) - p1.low) + self.low) - p2;

        const q2 = r1 / o.high;
        const p3 = Self.twoProduct(q2, o.high);
        const p4 = q2 * o.low;
        const r2 = ((r1 - p3.high) - p3.low) - p4;

        const q3 = r2 / o.high;
        return quickThreeSum(q1, q2, q3);
    }

    /// The same as '+=' operator c++, but the parameter is a `HD`.
    /// This will make `.low` more and more bigger, and may introduce more and more rounding errors,
    /// so you will take account of the renorm outsiede because invoke `twoSum` is expensive.
    pub fn addHDAssign(self: *Self, o: HD) void {
        const sum = HCD.twoSum(self.high, o);
        self.high = sum.high;
        self.low += sum.low;
    }

    /// The same as '+=' operator in c++, but the parameter is a `HCD`.
    /// The same as `addHDAssign`, but the parameter is a `HCD`.
    pub fn addHCDAssign(self: *Self, o: HCD) void {
        const sum = HCD.twoSum(self.high, o.high);
        self.high = sum.high;
        self.low += sum.low + o.low;
    }

    /// The same as '-=' operator c++, but the parameter is a `HD`.
    pub fn minusHDAssign(self: *Self, o: HD) void {
        const sum = HCD.twoSum(self.high, -o);
        self.high = sum.high;
        self.low += sum.low;
    }

    /// The same as '-=' operator c++, but the parameter is a `HCD`.
    pub fn minusHCDAssign(self: *Self, o: HCD) void {
        const sum = HCD.twoSum(self.high, -o.high);
        self.high = sum.high;
        self.low += sum.low - o.low;
    }

    /// The same as '*=' operator c++, but the parameter is a `HD`.
    pub fn multiplyHDAssign(self: *Self, o: HD) void {
        self.* = self.multiplyHD(o);
    }

    /// The same as '*=' operator c++, but the parameter is a `HCD`.
    pub fn multiplyHCDAssign(self: *Self, o: HCD) void {
        self.* = self.multiplyHCD(o);
    }

    /// The same as '/=' operator c++, but the parameter is a `HD`.
    pub fn divideHDAssign(self: *Self, o: HD) void {
        self.* = self.divideHD(o);
    }

    /// The same as '/=' operator c++, but the parameter is a `HCD`.
    pub fn divideHCDAssign(self: *Self, o: HCD) void {
        self.* = self.divideHCD(o);
    }

    pub fn cmp(self: Self, o: Self) std.math.Order {
        if (self.high < o.high) return .lt;
        if (self.high > o.high) return .gt;
        if (self.low < o.low) return .lt;
        if (self.low > o.low) return .gt;
        return .eq;
    }
};

test "twoSum stores rounded-away addend in low" {
    const res = HCD.twoSum(1e16, 1.0);

    try std.testing.expectEqual(@as(HD, 1e16), res.high);
    try std.testing.expectEqual(@as(HD, 1.0), res.low);
}

test "twoSum cancels equal opposite values" {
    const res = HCD.twoSum(1.0, -1.0);

    try std.testing.expectEqual(@as(HD, 0.0), res.high);
    try std.testing.expectEqual(@as(HD, 0.0), res.low);
}

test "split preserves zero" {
    const res = HCD.split(0.0);

    try std.testing.expectEqual(@as(HD, 0.0), res.high);
    try std.testing.expectEqual(@as(HD, 0.0), res.low);
}

test "split recombines positive value" {
    const value: HD = 1e16;
    const res = HCD.split(value);

    try std.testing.expectEqual(value, res.toHD());
}

test "split recombines negative value" {
    const value: HD = -1e16;
    const res = HCD.split(value);

    try std.testing.expectEqual(value, res.toHD());
}

test "twoProduct stores exact product error in low" {
    const a: HD = 1e16;
    const b: HD = 0.3;
    const high = a * b;
    const res = HCD.twoProduct(a, b);

    try std.testing.expectEqual(high, res.high);
    try std.testing.expectEqual(@mulAdd(HD, a, b, -high), res.low);
}

test "twoProduct returns zero for zero operand" {
    const res = HCD.twoProduct(123.0, 0.0);

    try std.testing.expectEqual(@as(HD, 0.0), res.high);
    try std.testing.expectEqual(@as(HD, 0.0), res.low);
}

test "initWithHD sets low to zero" {
    const res = HCD.initWithHD(2.5);

    try std.testing.expectEqual(@as(HD, 2.5), res.high);
    try std.testing.expectEqual(@as(HD, 0.0), res.low);
}

test "init stores high and low" {
    const res = HCD.init(2.5, 0.25);

    try std.testing.expectEqual(@as(HD, 2.5), res.high);
    try std.testing.expectEqual(@as(HD, 0.25), res.low);
}

test "clone copies high and low" {
    const value = HCD.init(2.5, 0.25);
    const res = value.clone();

    try std.testing.expectEqual(value.high, res.high);
    try std.testing.expectEqual(value.low, res.low);
}

test "toHD returns rounded sum" {
    const value = HCD.init(2.5, 0.25);

    try std.testing.expectEqual(@as(HD, 2.75), value.toHD());
}

test "addHD captures small addend in low" {
    const res = HCD.initWithHD(1e16).addHD(1.0);

    try std.testing.expectEqual(@as(HD, 1e16), res.high);
    try std.testing.expectEqual(@as(HD, 1.0), res.low);
}

test "addHCD adds high and low parts" {
    const res = HCD.init(1e16, 1.0).addHCD(HCD.init(1.0, 0.25));

    try std.testing.expectEqual(@as(HD, 1e16), res.high);
    try std.testing.expectEqual(@as(HD, 2.25), res.low);
}

test "minusHD captures small subtrahend in low" {
    const res = HCD.initWithHD(1e16).minusHD(1.0);

    try std.testing.expectEqual(@as(HD, 1e16), res.high);
    try std.testing.expectEqual(@as(HD, -1.0), res.low);
}

test "minusHCD subtracts high and low parts" {
    const res = HCD.init(1e16, 1.0).minusHCD(HCD.init(1.0, 0.25));

    try std.testing.expectEqual(@as(HD, 1e16), res.high);
    try std.testing.expectEqual(@as(HD, -0.25), res.low);
}

test "multiplyHD includes low part" {
    const res = HCD.init(1e16, 1.0).multiplyHD(2.0);

    try std.testing.expectEqual(@as(HD, 2e16), res.high);
    try std.testing.expectEqual(@as(HD, 2.0), res.low);
}

test "multiplyHD by zero returns zero" {
    const res = HCD.init(4.0, 0.5).multiplyHD(0.0);

    try std.testing.expectEqual(@as(HD, 0.0), res.high);
    try std.testing.expectEqual(@as(HD, 0.0), res.low);
}

test "multiplyHD by one preserves value" {
    const res = HCD.init(4.0, 0.5).multiplyHD(1.0);

    try std.testing.expectEqual(@as(HD, 4.5), res.toHD());
}

test "multiplyHD handles negative multiplier" {
    const res = HCD.init(4.0, 0.5).multiplyHD(-2.0);

    try std.testing.expectEqual(@as(HD, -9.0), res.toHD());
}

test "multiplyHCD by zero returns zero" {
    const res = HCD.init(4.0, 0.5).multiplyHCD(HCD.initWithHD(0.0));

    try std.testing.expectEqual(@as(HD, 0.0), res.high);
    try std.testing.expectEqual(@as(HD, 0.0), res.low);
}

test "multiplyHCD by one preserves value" {
    const res = HCD.init(4.0, 0.5).multiplyHCD(HCD.initWithHD(1.0));

    try std.testing.expectEqual(@as(HD, 4.5), res.toHD());
}

test "multiplyHCD uses multiplier low fast path result" {
    const res = HCD.init(1e16, 1.0).multiplyHCD(HCD.initWithHD(2.0));

    try std.testing.expectEqual(@as(HD, 2e16), res.high);
    try std.testing.expectEqual(@as(HD, 2.0), res.low);
}

test "multiplyHCD includes cross products" {
    const res = HCD.init(1.0, 1e-16).multiplyHCD(HCD.init(2.0, 1e-16));

    try std.testing.expectEqual(@as(HD, 2.0000000000000004), res.toHD());
}

test "multiplyHCD handles negative multiplier" {
    const res = HCD.init(4.0, 0.5).multiplyHCD(HCD.initWithHD(-2.0));

    try std.testing.expectEqual(@as(HD, -9.0), res.toHD());
}

test "divideHD divides by one" {
    const res = HCD.init(4.0, 0.5).divideHD(1.0);

    try std.testing.expectEqual(@as(HD, 4.5), res.toHD());
}

test "divideHD divides zero numerator" {
    const res = HCD.initWithHD(0.0).divideHD(2.0);

    try std.testing.expectEqual(@as(HD, 0.0), res.high);
    try std.testing.expectEqual(@as(HD, 0.0), res.low);
}

test "divideHD keeps low part contribution" {
    const res = HCD.init(1e16, 2.0).divideHD(2.0);

    try std.testing.expectEqual(@as(HD, 5000000000000001.0), res.toHD());
}

test "divideHD handles negative divisor" {
    const res = HCD.init(4.0, 0.5).divideHD(-2.0);

    try std.testing.expectEqual(@as(HD, -2.25), res.toHD());
}

test "divideHD returns infinity for nonzero divided by zero" {
    const res = HCD.initWithHD(1.0).divideHD(0.0);

    try std.testing.expect(std.math.isPositiveInf(res.high));
    try std.testing.expectEqual(@as(HD, 0.0), res.low);
}

test "divideHD returns nan for zero divided by zero" {
    const res = HCD.initWithHD(0.0).divideHD(0.0);

    try std.testing.expect(std.math.isNan(res.high));
    try std.testing.expectEqual(@as(HD, 0.0), res.low);
}

test "divideHCD divides by one" {
    const res = HCD.init(4.0, 0.5).divideHCD(HCD.initWithHD(1.0));

    try std.testing.expectEqual(@as(HD, 4.5), res.toHD());
}

test "divideHCD divides zero numerator" {
    const res = HCD.initWithHD(0.0).divideHCD(HCD.init(2.0, 0.25));

    try std.testing.expectEqual(@as(HD, 0.0), res.high);
    try std.testing.expectEqual(@as(HD, 0.0), res.low);
}

test "divideHCD accounts for denominator low part" {
    const res = HCD.init(1.0, 3e-16).divideHCD(HCD.init(1.0, 1e-16));

    try std.testing.expectEqual(@as(HD, 1.0000000000000002), res.toHD());
}

test "divideHCD handles negative divisor" {
    const res = HCD.init(4.0, 0.5).divideHCD(HCD.initWithHD(-2.0));

    try std.testing.expectEqual(@as(HD, -2.25), res.toHD());
}

test "divideHCD returns infinity for nonzero divided by zero" {
    const res = HCD.initWithHD(1.0).divideHCD(HCD.initWithHD(0.0));

    try std.testing.expect(std.math.isPositiveInf(res.high));
    try std.testing.expectEqual(@as(HD, 0.0), res.low);
}

test "addHDAssign mutates receiver" {
    var value = HCD.initWithHD(1e16);

    value.addHDAssign(1.0);

    try std.testing.expectEqual(@as(HD, 1e16), value.high);
    try std.testing.expectEqual(@as(HD, 1.0), value.low);
}

test "addHCDAssign mutates receiver" {
    var value = HCD.init(1e16, 1.0);

    value.addHCDAssign(HCD.init(1.0, 0.25));

    try std.testing.expectEqual(@as(HD, 1e16), value.high);
    try std.testing.expectEqual(@as(HD, 2.25), value.low);
}

test "minusHDAssign mutates receiver" {
    var value = HCD.initWithHD(1e16);

    value.minusHDAssign(1.0);

    try std.testing.expectEqual(@as(HD, 1e16), value.high);
    try std.testing.expectEqual(@as(HD, -1.0), value.low);
}

test "minusHCDAssign mutates receiver" {
    var value = HCD.init(1e16, 1.0);

    value.minusHCDAssign(HCD.init(1.0, 0.25));

    try std.testing.expectEqual(@as(HD, 1e16), value.high);
    try std.testing.expectEqual(@as(HD, -0.25), value.low);
}

test "multiplyHDAssign includes low part" {
    var value = HCD.init(1e16, 1.0);

    value.multiplyHDAssign(2.0);

    try std.testing.expectEqual(@as(HD, 2e16), value.high);
    try std.testing.expectEqual(@as(HD, 2.0), value.low);
}

test "multiplyHCDAssign includes cross products" {
    var value = HCD.init(1.0, 1e-16);

    value.multiplyHCDAssign(HCD.init(2.0, 1e-16));

    try std.testing.expectEqual(@as(HD, 2.0000000000000004), value.toHD());
}

test "divideHDAssign mutates receiver" {
    var value = HCD.init(4.0, 0.5);

    value.divideHDAssign(2.0);

    try std.testing.expectEqual(@as(HD, 2.25), value.toHD());
}

test "divideHCDAssign mutates receiver" {
    var value = HCD.init(4.0, 0.5);

    value.divideHCDAssign(HCD.initWithHD(2.0));

    try std.testing.expectEqual(@as(HD, 2.25), value.toHD());
}

test "cmp returns lt for smaller value" {
    try std.testing.expectEqual(std.math.Order.lt, HCD.initWithHD(1.0).cmp(HCD.initWithHD(2.0)));
}

test "cmp returns gt for greater value" {
    try std.testing.expectEqual(std.math.Order.gt, HCD.initWithHD(2.0).cmp(HCD.initWithHD(1.0)));
}

test "cmp returns eq for equal value" {
    try std.testing.expectEqual(std.math.Order.eq, HCD.initWithHD(1.0).cmp(HCD.initWithHD(1.0)));
}

test "cmp compares low part when high parts match" {
    try std.testing.expectEqual(std.math.Order.gt, HCD.init(1.0, 1e-16).cmp(HCD.initWithHD(1.0)));
}
