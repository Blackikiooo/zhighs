const std = @import("std");
const builtin = @import("builtin");

/// Double precision.
pub const HD = f64;

const has_fma = switch (builtin.cpu.arch) {
    // x86 系列
    .x86, .x86_64 => builtin.cpu.has(.x86, .fma),
    // AArch64 / ARM32
    .aarch64 => builtin.cpu.has(.aarch64, .fp_armv8),
    .arm => builtin.cpu.hasAny(.arm, &.{ .fp_armv8, .vfp4 }),
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

    /// Error-free transform for `a + b`.
    /// Algorithm: Knuth two-sum, returning `high = round(a + b)` and `low`
    /// such that `high + low` equals the exact real sum when no overflow occurs.
    /// Example: `const r = HCD.twoSum(1e16, 1.0);`
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

    /// Split one f64 into high and low parts.
    /// Algorithm: Dekker split with factor `2^27 + 1`; used by `twoProduct`
    /// when hardware FMA is unavailable. Like most Dekker split routines, this
    /// is intended for normal finite inputs; subnormal inputs may lose the
    /// exact split property.
    /// Example: `const parts = HCD.split(1e16);`
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

    /// Error-free transform for `a * b`.
    /// Algorithm: FMA residual when available, otherwise Dekker product using
    /// `split`; returns `high = round(a * b)` and product error in `low`.
    /// Example: `const p = HCD.twoProduct(1e16, 0.3);`
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

    /// Fast renormalization for a dominant high part plus a small correction.
    /// Algorithm: quick two-sum; valid when `abs(a) >= abs(b)`.
    /// Example: `const r = HCD.quickTwoSum(1e16, 1.0);`
    inline fn quickTwoSum(a: HD, b: HD) Self {
        @setFloatMode(.strict);

        const high = a + b;
        const low = b - (high - a);

        return Self{
            .high = high,
            .low = low,
        };
    }

    /// Fast renormalization of three terms where `a` is the dominant term.
    /// Algorithm: quick two-sum on `a + b`, then folds `c` into the residual.
    /// Example: `const r = HCD.quickThreeSum(q1, q2, q3);`
    inline fn quickThreeSum(a: HD, b: HD, c: HD) Self {
        @setFloatMode(.strict);

        const sum = a + b;
        const err = b - (sum - a);
        return quickTwoSum(sum, err + c);
    }

    /// Full two-term renormalization.
    /// Algorithm: `twoSum(a, b)`; use when `b` may not be tiny relative to `a`.
    /// Example: `const r = HCD.renorm2(value.high, value.low);`
    inline fn renorm2(a: HD, b: HD) Self {
        @setFloatMode(.strict);

        return twoSum(a, b);
    }

    /// Construct an HCD from one f64.
    /// Algorithm: exact embedding with `low = 0`.
    /// Example: `const x = HCD.initWithHD(1.0);`
    pub inline fn initWithHD(value: HD) Self {
        return Self{
            .high = value,
            .low = 0.0,
        };
    }

    /// Construct an HCD from explicit high and low parts.
    /// Algorithm: direct storage; no renormalization is performed.
    /// Example: `const x = HCD.init(1e16, 1.0);`
    pub inline fn init(high: HD, low: HD) Self {
        return Self{
            .high = high,
            .low = low,
        };
    }

    /// Copy this value.
    /// Algorithm: direct field copy; equivalent to value assignment, but clearer.
    /// Example: `const y = x.clone();`
    pub inline fn clone(self: Self) Self {
        return Self{
            .high = self.high,
            .low = self.low,
        };
    }

    /// Convert to one f64.
    /// Algorithm: rounded `high + low`; this intentionally loses extra precision.
    /// Example: `const y: HD = x.toHD();`
    pub inline fn toHD(self: Self) HD {
        return self.high + self.low;
    }

    /// Return a normalized copy of this HCD.
    /// Algorithm: full two-term renormalization of `high` and `low`.
    /// Example: `const y = x.renorm();`
    pub inline fn renorm(self: Self) HCD {
        return renorm2(self.high, self.low);
    }

    /// Normalize this HCD in place.
    /// Algorithm: assign `self.renorm()` back to `self`.
    /// Example: `x.renormAssign();`
    pub inline fn renormAssign(self: *Self) void {
        self.* = self.renorm();
    }

    /// Add one f64 without final renormalization.
    /// Algorithm: `twoSum(high, o)` and accumulate the residual into `low`.
    /// Example: `x.addHDAssignFast(1.0); x.renormAssign();`
    pub inline fn addHDFast(self: Self, o: HD) HCD {
        const sum = twoSum(self.high, o);
        return Self.init(sum.high, self.low + sum.low);
    }

    /// Add one f64 without final renormalization when `self.high` dominates.
    /// Algorithm: quick two-sum on `high + o`, then accumulate the residual.
    /// Requires: `abs(self.high) >= abs(o)`.
    /// Example: `x.addHDOrderedAssignFast(delta); x.renormAssign();`
    pub inline fn addHDOrderedFast(self: Self, o: HD) HCD {
        const sum = quickTwoSum(self.high, o);
        return Self.init(sum.high, self.low + sum.low);
    }

    /// Add one f64 and return a normalized result.
    /// Algorithm: `addHDFast(o).renorm()`; safer, but costs an extra `twoSum`.
    /// Example: `const y = x.addHD(1.0);`
    pub inline fn addHD(self: Self, o: HD) HCD {
        return self.addHDFast(o).renorm();
    }

    /// Add one f64 and return a normalized result when `self.high` dominates.
    /// Algorithm: `addHDOrderedFast(o).renorm()`.
    /// Requires: `abs(self.high) >= abs(o)`.
    /// Example: `const y = x.addHDOrdered(delta);`
    pub inline fn addHDOrdered(self: Self, o: HD) HCD {
        return self.addHDOrderedFast(o).renorm();
    }

    /// Add another HCD without final renormalization.
    /// Algorithm: `twoSum(high, o.high)` and accumulate both low parts.
    /// Example: `x.addHCDAssignFast(y); x.renormAssign();`
    pub inline fn addHCDFast(self: Self, o: HCD) HCD {
        const sum = twoSum(self.high, o.high);
        return Self.init(sum.high, self.low + sum.low + o.low);
    }

    /// Add another HCD without final renormalization when `self.high` dominates.
    /// Algorithm: quick two-sum on the high parts, then accumulate low parts.
    /// Requires: `abs(self.high) >= abs(o.high)`.
    /// Example: `x.addHCDOrderedAssignFast(delta); x.renormAssign();`
    pub inline fn addHCDOrderedFast(self: Self, o: HCD) HCD {
        const sum = quickTwoSum(self.high, o.high);
        return Self.init(sum.high, self.low + sum.low + o.low);
    }

    /// Add another HCD and return a normalized result.
    /// Algorithm: `addHCDFast(o).renorm()`; use when you do not want to manage
    /// delayed renormalization manually.
    /// Example: `const z = x.addHCD(y);`
    pub inline fn addHCD(self: Self, o: HCD) HCD {
        return self.addHCDFast(o).renorm();
    }

    /// Add another HCD and return a normalized result when `self.high` dominates.
    /// Algorithm: `addHCDOrderedFast(o).renorm()`.
    /// Requires: `abs(self.high) >= abs(o.high)`.
    /// Example: `const z = x.addHCDOrdered(delta);`
    pub inline fn addHCDOrdered(self: Self, o: HCD) HCD {
        return self.addHCDOrderedFast(o).renorm();
    }

    /// Subtract one f64 without final renormalization.
    /// Algorithm: `twoSum(high, -o)` and accumulate the residual into `low`.
    /// Example: `x.minusHDAssignFast(1.0); x.renormAssign();`
    pub inline fn minusHDFast(self: Self, o: HD) HCD {
        const sum = twoSum(self.high, -o);
        return Self.init(sum.high, self.low + sum.low);
    }

    /// Subtract one f64 without final renormalization when `self.high` dominates.
    /// Algorithm: quick two-sum on `high - o`, then accumulate the residual.
    /// Requires: `abs(self.high) >= abs(o)`.
    /// Example: `x.minusHDOrderedAssignFast(delta); x.renormAssign();`
    pub inline fn minusHDOrderedFast(self: Self, o: HD) HCD {
        const sum = quickTwoSum(self.high, -o);
        return Self.init(sum.high, self.low + sum.low);
    }

    /// Subtract one f64 and return a normalized result.
    /// Algorithm: `minusHDFast(o).renorm()`.
    /// Example: `const y = x.minusHD(1.0);`
    pub inline fn minusHD(self: Self, o: HD) HCD {
        return self.minusHDFast(o).renorm();
    }

    /// Subtract one f64 and return a normalized result when `self.high` dominates.
    /// Algorithm: `minusHDOrderedFast(o).renorm()`.
    /// Requires: `abs(self.high) >= abs(o)`.
    /// Example: `const y = x.minusHDOrdered(delta);`
    pub inline fn minusHDOrdered(self: Self, o: HD) HCD {
        return self.minusHDOrderedFast(o).renorm();
    }

    /// Subtract another HCD without final renormalization.
    /// Algorithm: `twoSum(high, -o.high)` and subtract `o.low` from the low
    /// accumulator.
    /// Example: `x.minusHCDAssignFast(y); x.renormAssign();`
    pub inline fn minusHCDFast(self: Self, o: HCD) HCD {
        const sum = twoSum(self.high, -o.high);
        return Self.init(sum.high, self.low + sum.low - o.low);
    }

    /// Subtract another HCD without final renormalization when `self.high`
    /// dominates.
    /// Algorithm: quick two-sum on the high parts, then accumulate low parts.
    /// Requires: `abs(self.high) >= abs(o.high)`.
    /// Example: `x.minusHCDOrderedAssignFast(delta); x.renormAssign();`
    pub inline fn minusHCDOrderedFast(self: Self, o: HCD) HCD {
        const sum = quickTwoSum(self.high, -o.high);
        return Self.init(sum.high, self.low + sum.low - o.low);
    }

    /// Subtract another HCD and return a normalized result.
    /// Algorithm: `minusHCDFast(o).renorm()`.
    /// Example: `const z = x.minusHCD(y);`
    pub inline fn minusHCD(self: Self, o: HCD) HCD {
        return self.minusHCDFast(o).renorm();
    }

    /// Subtract another HCD and return a normalized result when `self.high`
    /// dominates.
    /// Algorithm: `minusHCDOrderedFast(o).renorm()`.
    /// Requires: `abs(self.high) >= abs(o.high)`.
    /// Example: `const z = x.minusHCDOrdered(delta);`
    pub inline fn minusHCDOrdered(self: Self, o: HCD) HCD {
        return self.minusHCDOrderedFast(o).renorm();
    }

    /// Multiply by one f64.
    /// Algorithm: double-double by double fast product: `twoProduct(high, o)`
    /// plus the cross term `low * o`, then `quickTwoSum`.
    /// Example: `const y = x.multiplyHD(2.0);`
    pub inline fn multiplyHD(self: Self, o: HD) HCD {
        const product = twoProduct(self.high, o);
        return quickTwoSum(product.high, product.low + self.low * o);
    }

    /// Multiply by another HCD.
    /// Algorithm: fast double-double product using `twoProduct(high, o.high)`
    /// and cross terms `high * o.low + low * o.high`; ignores `low * o.low`.
    /// Example: `const z = x.multiplyHCD(y);`
    pub inline fn multiplyHCD(self: Self, o: HCD) HCD {
        if (o.low == 0.0) {
            return self.multiplyHD(o.high);
        }

        const product = twoProduct(self.high, o.high);
        const cross_products = product.low + self.high * o.low + self.low * o.high;
        return quickTwoSum(product.high, cross_products);
    }

    /// Divide by one f64.
    /// Algorithm: one quotient estimate plus one residual correction:
    /// `q1 = high / o`, correct with `(self - q1 * o) / o`, then quick sum.
    /// Example: `const y = x.divideHD(2.0);`
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

    /// Divide by another HCD with one quotient correction.
    /// Algorithm: fast double-double division: estimate `q1 = high / den.high`,
    /// compute the expanded residual of `self - q1 * den`, then return
    /// `quickTwoSum(q1, residual / den.high)`.
    /// This is faster than `divideHCD`, but skips the second correction term.
    /// Recommendation: use this explicitly in solver hot paths when one
    /// correction is enough for the caller's tolerance; keep `divideHCD` as the
    /// default when stability matters more than throughput.
    /// Example: `const z = x.divideHCDFast(y);`
    pub inline fn divideHCDFast(self: Self, o: HCD) HCD {
        @setFloatMode(.strict);

        var den = o;
        if (den.high == 0.0) {
            @branchHint(.unlikely);
            den = den.renorm();
            if (den.high == 0.0 and den.low == 0.0) {
                return Self.initWithHD(self.toHD() / den.toHD());
            }
        }
        if (den.low == 0.0) {
            return self.divideHD(den.high);
        }

        const q1 = self.high / den.high;
        const p1 = Self.twoProduct(q1, den.high);
        const p2 = q1 * den.low;
        const r1 = (((self.high - p1.high) - p1.low) + self.low) - p2;
        const q2 = r1 / den.high;

        return quickTwoSum(q1, q2);
    }

    /// Divide by another HCD.
    /// Algorithm: two Newton-style quotient correction steps using expanded
    /// double-double residual arithmetic; renormalizes rare `o.high == 0`
    /// inputs and falls back to `divideHD` if `o.low = 0`.
    /// Recommendation: use this as the default safe division path. In measured
    /// hot paths with sufficient error budget, choose `divideHCDFast` manually.
    /// Example: `const z = x.divideHCD(y);`
    pub inline fn divideHCD(self: Self, o: HCD) HCD {
        @setFloatMode(.strict);

        var den = o;
        if (den.high == 0.0) {
            @branchHint(.unlikely);
            den = den.renorm();
            if (den.high == 0.0 and den.low == 0.0) {
                return Self.initWithHD(self.toHD() / den.toHD());
            }
        }
        if (den.low == 0.0) {
            return self.divideHD(den.high);
        }

        const q1 = self.high / den.high;
        const p1 = Self.twoProduct(q1, den.high);
        const p2 = q1 * den.low;
        const r1 = (((self.high - p1.high) - p1.low) + self.low) - p2;

        const q2 = r1 / den.high;
        const p3 = Self.twoProduct(q2, den.high);
        const p4 = q2 * den.low;
        const r2 = ((r1 - p3.high) - p3.low) - p4;

        const q3 = r2 / den.high;
        return quickThreeSum(q1, q2, q3);
    }

    /// Add one f64 in place with automatic renormalization.
    /// Algorithm: assign `addHD(o)` back to `self`.
    /// Example: `x.addHDAssign(1.0);`
    pub inline fn addHDAssign(self: *Self, o: HD) void {
        self.* = self.addHD(o);
    }

    /// Add one f64 in place without final renormalization.
    /// Algorithm: assign `addHDFast(o)` back to `self`.
    /// Example: `x.addHDAssignFast(1.0); x.renormAssign();`
    pub inline fn addHDAssignFast(self: *Self, o: HD) void {
        self.* = self.addHDFast(o);
    }

    /// Add one f64 in place with automatic renormalization when `self.high`
    /// dominates.
    /// Algorithm: assign `addHDOrdered(o)` back to `self`.
    /// Requires: `abs(self.high) >= abs(o)`.
    /// Example: `x.addHDOrderedAssign(delta);`
    pub inline fn addHDOrderedAssign(self: *Self, o: HD) void {
        self.* = self.addHDOrdered(o);
    }

    /// Add one f64 in place without final renormalization when `self.high`
    /// dominates.
    /// Algorithm: assign `addHDOrderedFast(o)` back to `self`.
    /// Requires: `abs(self.high) >= abs(o)`.
    /// Example: `x.addHDOrderedAssignFast(delta); x.renormAssign();`
    pub inline fn addHDOrderedAssignFast(self: *Self, o: HD) void {
        self.* = self.addHDOrderedFast(o);
    }

    /// Add another HCD in place with automatic renormalization.
    /// Algorithm: assign `addHCD(o)` back to `self`.
    /// Example: `x.addHCDAssign(y);`
    pub inline fn addHCDAssign(self: *Self, o: HCD) void {
        self.* = self.addHCD(o);
    }

    /// Add another HCD in place without final renormalization.
    /// Algorithm: assign `addHCDFast(o)` back to `self`.
    /// Example: `x.addHCDAssignFast(y); x.renormAssign();`
    pub inline fn addHCDAssignFast(self: *Self, o: HCD) void {
        self.* = self.addHCDFast(o);
    }

    /// Add another HCD in place with automatic renormalization when `self.high`
    /// dominates.
    /// Algorithm: assign `addHCDOrdered(o)` back to `self`.
    /// Requires: `abs(self.high) >= abs(o.high)`.
    /// Example: `x.addHCDOrderedAssign(delta);`
    pub inline fn addHCDOrderedAssign(self: *Self, o: HCD) void {
        self.* = self.addHCDOrdered(o);
    }

    /// Add another HCD in place without final renormalization when `self.high`
    /// dominates.
    /// Algorithm: assign `addHCDOrderedFast(o)` back to `self`.
    /// Requires: `abs(self.high) >= abs(o.high)`.
    /// Example: `x.addHCDOrderedAssignFast(delta); x.renormAssign();`
    pub inline fn addHCDOrderedAssignFast(self: *Self, o: HCD) void {
        self.* = self.addHCDOrderedFast(o);
    }

    /// Subtract one f64 in place with automatic renormalization.
    /// Algorithm: assign `minusHD(o)` back to `self`.
    /// Example: `x.minusHDAssign(1.0);`
    pub inline fn minusHDAssign(self: *Self, o: HD) void {
        self.* = self.minusHD(o);
    }

    /// Subtract one f64 in place without final renormalization.
    /// Algorithm: assign `minusHDFast(o)` back to `self`.
    /// Example: `x.minusHDAssignFast(1.0); x.renormAssign();`
    pub inline fn minusHDAssignFast(self: *Self, o: HD) void {
        self.* = self.minusHDFast(o);
    }

    /// Subtract one f64 in place with automatic renormalization when `self.high`
    /// dominates.
    /// Algorithm: assign `minusHDOrdered(o)` back to `self`.
    /// Requires: `abs(self.high) >= abs(o)`.
    /// Example: `x.minusHDOrderedAssign(delta);`
    pub inline fn minusHDOrderedAssign(self: *Self, o: HD) void {
        self.* = self.minusHDOrdered(o);
    }

    /// Subtract one f64 in place without final renormalization when `self.high`
    /// dominates.
    /// Algorithm: assign `minusHDOrderedFast(o)` back to `self`.
    /// Requires: `abs(self.high) >= abs(o)`.
    /// Example: `x.minusHDOrderedAssignFast(delta); x.renormAssign();`
    pub inline fn minusHDOrderedAssignFast(self: *Self, o: HD) void {
        self.* = self.minusHDOrderedFast(o);
    }

    /// Subtract another HCD in place with automatic renormalization.
    /// Algorithm: assign `minusHCD(o)` back to `self`.
    /// Example: `x.minusHCDAssign(y);`
    pub inline fn minusHCDAssign(self: *Self, o: HCD) void {
        self.* = self.minusHCD(o);
    }

    /// Subtract another HCD in place without final renormalization.
    /// Algorithm: assign `minusHCDFast(o)` back to `self`.
    /// Example: `x.minusHCDAssignFast(y); x.renormAssign();`
    pub inline fn minusHCDAssignFast(self: *Self, o: HCD) void {
        self.* = self.minusHCDFast(o);
    }

    /// Subtract another HCD in place with automatic renormalization when
    /// `self.high` dominates.
    /// Algorithm: assign `minusHCDOrdered(o)` back to `self`.
    /// Requires: `abs(self.high) >= abs(o.high)`.
    /// Example: `x.minusHCDOrderedAssign(delta);`
    pub inline fn minusHCDOrderedAssign(self: *Self, o: HCD) void {
        self.* = self.minusHCDOrdered(o);
    }

    /// Subtract another HCD in place without final renormalization when
    /// `self.high` dominates.
    /// Algorithm: assign `minusHCDOrderedFast(o)` back to `self`.
    /// Requires: `abs(self.high) >= abs(o.high)`.
    /// Example: `x.minusHCDOrderedAssignFast(delta); x.renormAssign();`
    pub inline fn minusHCDOrderedAssignFast(self: *Self, o: HCD) void {
        self.* = self.minusHCDOrderedFast(o);
    }

    /// Multiply by one f64 in place.
    /// Algorithm: assign `multiplyHD(o)` back to `self`.
    /// Example: `x.multiplyHDAssign(2.0);`
    pub inline fn multiplyHDAssign(self: *Self, o: HD) void {
        self.* = self.multiplyHD(o);
    }

    /// Multiply by another HCD in place.
    /// Algorithm: assign `multiplyHCD(o)` back to `self`.
    /// Example: `x.multiplyHCDAssign(y);`
    pub inline fn multiplyHCDAssign(self: *Self, o: HCD) void {
        self.* = self.multiplyHCD(o);
    }

    /// Divide by one f64 in place.
    /// Algorithm: assign `divideHD(o)` back to `self`.
    /// Example: `x.divideHDAssign(2.0);`
    pub inline fn divideHDAssign(self: *Self, o: HD) void {
        self.* = self.divideHD(o);
    }

    /// Divide by another HCD in place.
    /// Algorithm: assign `divideHCD(o)` back to `self`.
    /// Recommendation: default assign path for stability; use
    /// `divideHCDAssignFast` only when the hot path has enough error budget.
    /// Example: `x.divideHCDAssign(y);`
    pub inline fn divideHCDAssign(self: *Self, o: HCD) void {
        self.* = self.divideHCD(o);
    }

    /// Divide by another HCD in place with one quotient correction.
    /// Algorithm: assign `divideHCDFast(o)` back to `self`.
    /// Recommendation: hot-path assign variant; caller chooses it explicitly
    /// when one correction is enough for the surrounding algorithm.
    /// Example: `x.divideHCDAssignFast(y);`
    pub inline fn divideHCDAssignFast(self: *Self, o: HCD) void {
        self.* = self.divideHCDFast(o);
    }

    /// Compare two already-normalized HCD values.
    /// Algorithm: lexicographic compare on `(high, low)`; faster but assumes
    /// both values have been normalized.
    /// Example: `const order = x.renorm().cmpFast(y.renorm());`
    pub inline fn cmpFast(self: Self, o: Self) std.math.Order {
        if (self.high < o.high) return .lt;
        if (self.high > o.high) return .gt;
        if (self.low < o.low) return .lt;
        if (self.low > o.low) return .gt;
        return .eq;
    }

    /// Compare two HCD values safely.
    /// Algorithm: renormalize both operands, then use lexicographic compare.
    /// Example: `const order = x.cmp(y);`
    pub inline fn cmp(self: Self, o: Self) std.math.Order {
        return self.renorm().cmpFast(o.renorm());
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

test "renorm folds large low part into high" {
    const res = HCD.init(1e16, 2.25).renorm();

    try std.testing.expectEqual(@as(HD, 10000000000000002.0), res.high);
    try std.testing.expectEqual(@as(HD, 0.25), res.low);
}

test "renormAssign folds large low part into high" {
    var value = HCD.init(1e16, 2.25);

    value.renormAssign();

    try std.testing.expectEqual(@as(HD, 10000000000000002.0), value.high);
    try std.testing.expectEqual(@as(HD, 0.25), value.low);
}

test "addHD captures small addend in low" {
    const res = HCD.initWithHD(1e16).addHD(1.0);

    try std.testing.expectEqual(@as(HD, 1e16), res.high);
    try std.testing.expectEqual(@as(HD, 1.0), res.low);
}

test "addHDFast captures small addend in low" {
    const res = HCD.initWithHD(1e16).addHDFast(1.0);

    try std.testing.expectEqual(@as(HD, 1e16), res.high);
    try std.testing.expectEqual(@as(HD, 1.0), res.low);
}

test "addHDOrderedFast captures ordered small addend" {
    const res = HCD.initWithHD(1e16).addHDOrderedFast(1.0);

    try std.testing.expectEqual(@as(HD, 1e16), res.high);
    try std.testing.expectEqual(@as(HD, 1.0), res.low);
}

test "addHDOrdered renormalizes ordered small addend" {
    const res = HCD.init(1e16, 1.0).addHDOrdered(1.25);

    try std.testing.expectEqual(@as(HD, 10000000000000002.0), res.high);
    try std.testing.expectEqual(@as(HD, 0.25), res.low);
}

test "addHCD adds high and low parts" {
    const res = HCD.init(1e16, 1.0).addHCD(HCD.init(1.0, 0.25));

    try std.testing.expectEqual(@as(HD, 10000000000000002.0), res.high);
    try std.testing.expectEqual(@as(HD, 0.25), res.low);
}

test "addHCDFast leaves renormalization to caller" {
    const res = HCD.init(1e16, 1.0).addHCDFast(HCD.init(1.0, 0.25));

    try std.testing.expectEqual(@as(HD, 1e16), res.high);
    try std.testing.expectEqual(@as(HD, 2.25), res.low);
}

test "addHCDFast can be manually renormalized" {
    const res = HCD.init(1e16, 1.0).addHCDFast(HCD.init(1.0, 1.25)).renorm();

    try std.testing.expectEqual(@as(HD, 10000000000000004.0), res.high);
    try std.testing.expectEqual(@as(HD, -0.75), res.low);
}

test "addHCDOrderedFast leaves ordered renormalization to caller" {
    const res = HCD.init(1e16, 1.0).addHCDOrderedFast(HCD.init(1.0, 0.25));

    try std.testing.expectEqual(@as(HD, 1e16), res.high);
    try std.testing.expectEqual(@as(HD, 2.25), res.low);
}

test "addHCDOrdered renormalizes ordered HCD addend" {
    const res = HCD.init(1e16, 1.0).addHCDOrdered(HCD.init(1.0, 1.25));

    try std.testing.expectEqual(@as(HD, 10000000000000004.0), res.high);
    try std.testing.expectEqual(@as(HD, -0.75), res.low);
}

test "addHCD handles high cancellation with larger low part" {
    const res = HCD.init(1e16, 5.0).addHCD(HCD.init(-1e16, 7.0));

    try std.testing.expectEqual(@as(HD, 12.0), res.high);
    try std.testing.expectEqual(@as(HD, 0.0), res.low);
}

test "minusHD captures small subtrahend in low" {
    const res = HCD.initWithHD(1e16).minusHD(1.0);

    try std.testing.expectEqual(@as(HD, 1e16), res.high);
    try std.testing.expectEqual(@as(HD, -1.0), res.low);
}

test "minusHDFast captures small subtrahend in low" {
    const res = HCD.initWithHD(1e16).minusHDFast(1.0);

    try std.testing.expectEqual(@as(HD, 1e16), res.high);
    try std.testing.expectEqual(@as(HD, -1.0), res.low);
}

test "minusHDOrderedFast captures ordered small subtrahend" {
    const res = HCD.initWithHD(1e16).minusHDOrderedFast(1.0);

    try std.testing.expectEqual(@as(HD, 1e16), res.high);
    try std.testing.expectEqual(@as(HD, -1.0), res.low);
}

test "minusHDOrdered renormalizes ordered small subtrahend" {
    const res = HCD.init(1e16, 1.0).minusHDOrdered(1.25);

    try std.testing.expectEqual(@as(HD, 1e16), res.high);
    try std.testing.expectEqual(@as(HD, -0.25), res.low);
}

test "minusHCD subtracts high and low parts" {
    const res = HCD.init(1e16, 1.0).minusHCD(HCD.init(1.0, 0.25));

    try std.testing.expectEqual(@as(HD, 1e16), res.high);
    try std.testing.expectEqual(@as(HD, -0.25), res.low);
}

test "minusHCDFast leaves renormalization to caller" {
    const res = HCD.init(1e16, 1.0).minusHCDFast(HCD.init(1.0, 0.25));

    try std.testing.expectEqual(@as(HD, 1e16), res.high);
    try std.testing.expectEqual(@as(HD, -0.25), res.low);
}

test "minusHCDOrderedFast leaves ordered renormalization to caller" {
    const res = HCD.init(1e16, 1.0).minusHCDOrderedFast(HCD.init(1.0, 0.25));

    try std.testing.expectEqual(@as(HD, 1e16), res.high);
    try std.testing.expectEqual(@as(HD, -0.25), res.low);
}

test "minusHCDOrdered renormalizes ordered HCD subtrahend" {
    const res = HCD.init(1e16, 1.0).minusHCDOrdered(HCD.init(1.0, 1.25));

    try std.testing.expectEqual(@as(HD, 9999999999999998.0), res.high);
    try std.testing.expectEqual(@as(HD, 0.75), res.low);
}

test "minusHCD handles high cancellation with larger low part" {
    const res = HCD.init(1e16, 5.0).minusHCD(HCD.init(1e16, -7.0));

    try std.testing.expectEqual(@as(HD, 12.0), res.high);
    try std.testing.expectEqual(@as(HD, 0.0), res.low);
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

test "divideHCD handles denominator stored in low part" {
    const res = HCD.initWithHD(4.0).divideHCD(HCD.init(0.0, 2.0));

    try std.testing.expectEqual(@as(HD, 2.0), res.toHD());
}

test "divideHCD accounts for denominator low part" {
    const res = HCD.init(1.0, 3e-16).divideHCD(HCD.init(1.0, 1e-16));

    try std.testing.expectEqual(@as(HD, 1.0000000000000002), res.toHD());
}

test "divideHCDFast accounts for denominator low part" {
    const res = HCD.init(1.0, 3e-16).divideHCDFast(HCD.init(1.0, 1e-16));

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

test "addHDAssignFast mutates receiver" {
    var value = HCD.initWithHD(1e16);

    value.addHDAssignFast(1.0);

    try std.testing.expectEqual(@as(HD, 1e16), value.high);
    try std.testing.expectEqual(@as(HD, 1.0), value.low);
}

test "addHDOrderedAssign mutates and renormalizes receiver" {
    var value = HCD.init(1e16, 1.0);

    value.addHDOrderedAssign(1.25);

    try std.testing.expectEqual(@as(HD, 10000000000000002.0), value.high);
    try std.testing.expectEqual(@as(HD, 0.25), value.low);
}

test "addHDOrderedAssignFast leaves receiver unrenormalized" {
    var value = HCD.init(1e16, 1.0);

    value.addHDOrderedAssignFast(1.0);

    try std.testing.expectEqual(@as(HD, 1e16), value.high);
    try std.testing.expectEqual(@as(HD, 2.0), value.low);
}

test "addHCDAssign mutates receiver" {
    var value = HCD.init(1e16, 1.0);

    value.addHCDAssign(HCD.init(1.0, 0.25));

    try std.testing.expectEqual(@as(HD, 10000000000000002.0), value.high);
    try std.testing.expectEqual(@as(HD, 0.25), value.low);
}

test "addHCDAssign renormalizes accumulated low part" {
    var value = HCD.init(1e16, 1.0);

    value.addHCDAssign(HCD.init(1.0, 1.25));

    try std.testing.expectEqual(@as(HD, 10000000000000004.0), value.high);
    try std.testing.expectEqual(@as(HD, -0.75), value.low);
}

test "addHCDAssignFast leaves caller-controlled renormalization" {
    var value = HCD.init(1e16, 1.0);

    value.addHCDAssignFast(HCD.init(1.0, 1.25));

    try std.testing.expectEqual(@as(HD, 1e16), value.high);
    try std.testing.expectEqual(@as(HD, 3.25), value.low);
}

test "addHCDOrderedAssign mutates and renormalizes receiver" {
    var value = HCD.init(1e16, 1.0);

    value.addHCDOrderedAssign(HCD.init(1.0, 1.25));

    try std.testing.expectEqual(@as(HD, 10000000000000004.0), value.high);
    try std.testing.expectEqual(@as(HD, -0.75), value.low);
}

test "addHCDOrderedAssignFast leaves receiver unrenormalized" {
    var value = HCD.init(1e16, 1.0);

    value.addHCDOrderedAssignFast(HCD.init(1.0, 1.25));

    try std.testing.expectEqual(@as(HD, 1e16), value.high);
    try std.testing.expectEqual(@as(HD, 3.25), value.low);
}

test "minusHDAssign mutates receiver" {
    var value = HCD.initWithHD(1e16);

    value.minusHDAssign(1.0);

    try std.testing.expectEqual(@as(HD, 1e16), value.high);
    try std.testing.expectEqual(@as(HD, -1.0), value.low);
}

test "minusHDAssignFast mutates receiver" {
    var value = HCD.initWithHD(1e16);

    value.minusHDAssignFast(1.0);

    try std.testing.expectEqual(@as(HD, 1e16), value.high);
    try std.testing.expectEqual(@as(HD, -1.0), value.low);
}

test "minusHDOrderedAssign mutates and renormalizes receiver" {
    var value = HCD.init(1e16, 1.0);

    value.minusHDOrderedAssign(1.25);

    try std.testing.expectEqual(@as(HD, 1e16), value.high);
    try std.testing.expectEqual(@as(HD, -0.25), value.low);
}

test "minusHDOrderedAssignFast leaves receiver unrenormalized" {
    var value = HCD.init(1e16, 1.0);

    value.minusHDOrderedAssignFast(1.0);

    try std.testing.expectEqual(@as(HD, 1e16), value.high);
    try std.testing.expectEqual(@as(HD, 0.0), value.low);
}

test "minusHCDAssign mutates receiver" {
    var value = HCD.init(1e16, 1.0);

    value.minusHCDAssign(HCD.init(1.0, 0.25));

    try std.testing.expectEqual(@as(HD, 1e16), value.high);
    try std.testing.expectEqual(@as(HD, -0.25), value.low);
}

test "minusHCDAssignFast leaves caller-controlled renormalization" {
    var value = HCD.init(1e16, 1.0);

    value.minusHCDAssignFast(HCD.init(1.0, 0.25));

    try std.testing.expectEqual(@as(HD, 1e16), value.high);
    try std.testing.expectEqual(@as(HD, -0.25), value.low);
}

test "minusHCDOrderedAssign mutates and renormalizes receiver" {
    var value = HCD.init(1e16, 1.0);

    value.minusHCDOrderedAssign(HCD.init(1.0, 1.25));

    try std.testing.expectEqual(@as(HD, 9999999999999998.0), value.high);
    try std.testing.expectEqual(@as(HD, 0.75), value.low);
}

test "minusHCDOrderedAssignFast leaves receiver unrenormalized" {
    var value = HCD.init(1e16, 1.0);

    value.minusHCDOrderedAssignFast(HCD.init(1.0, 1.25));

    try std.testing.expectEqual(@as(HD, 1e16), value.high);
    try std.testing.expectEqual(@as(HD, -1.25), value.low);
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

test "divideHCDAssignFast mutates receiver" {
    var value = HCD.init(1.0, 3e-16);

    value.divideHCDAssignFast(HCD.init(1.0, 1e-16));

    try std.testing.expectEqual(@as(HD, 1.0000000000000002), value.toHD());
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

test "cmpFast compares normalized operands" {
    try std.testing.expectEqual(std.math.Order.gt, HCD.init(2.0, 0.0).cmpFast(HCD.init(1.0, 0.5)));
}

test "cmp compares low part when high parts match" {
    try std.testing.expectEqual(std.math.Order.gt, HCD.init(1.0, 1e-16).cmp(HCD.initWithHD(1.0)));
}

test "cmp renormalizes noncanonical operands" {
    const a = HCD.init(1e16, 3.25);
    const b = HCD.init(10000000000000004.0, -0.75);

    try std.testing.expectEqual(std.math.Order.eq, a.cmp(b));
}
