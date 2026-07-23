//! Exact seed-zero random vectors used by the pinned HiGHS simplex.
//!
//! HiGHS consumes one column shuffle, one all-variable shuffle, then one
//! fraction per all-variable index. Cost perturbation and CHUZC4 tie-breaking
//! must share this sequence or an otherwise source-equivalent ratio test can
//! follow a different pivot path.

const std = @import("std");

const constants = [_]u64{
    0xc8497d2a400d9551, 0x80c8963be3e4c2f3, 0x042d8680e260ae5b, 0x8a183895eeac1536,
    0xa94e9c75f80ad6de, 0x7e92251dec62835e, 0x07294165cb671455, 0x89b0f6212b0a4292,
    0x31900011b96bf554, 0xa44540f8eee2094f, 0xce7ffd372e4c64fc, 0x51c9d471bfe6a10f,
    0x758c2a674483826f, 0xf91a20abe63f8b02, 0xc2a069024a1fcc6f, 0xd5bb18b70c5dbd59,
    0xd510adac6d1ae289, 0x571d069b23050a79, 0x60873b8872933e06, 0x780481cc19670350,
    0x7a48551760216885, 0xb5d68b918231e6ca, 0xa7e5571699aa5274, 0x7b6d309b2cfdcf01,
    0x04e77c3d474daeff, 0x4dbf099fd7247031, 0x5d70dca901130beb, 0x9f8b5f0df4182499,
    0x293a74c9686092da, 0xd09bdab6840f52b3, 0xc05d47f3ab302263, 0x6b79e62b884b65d6,
    0xa581106fc980c34d, 0xf081b7145ea2293e, 0xfb27243dd7c3f5ad, 0x5211bf8860ea667f,
    0x9455e65cb2385e7f, 0x0dfaf6731b449b33, 0x4ec98b3c6f5e68c7, 0x007bfd4a42ae936b,
    0x65c93061f8674518, 0x640816f17127c5d1, 0x6dd4bab17b7c3a74, 0x34d9268c256fa1ba,
    0x0b4d0c6b5b50d7f4, 0x30aa965bc9fadaff, 0xc0ac1d0c2771404d, 0xc5e64509abb76ef2,
    0xd606b11990624a36, 0x0d3f05d242ce2fb7, 0x469a803cb276fe32, 0xa4a44d177a3e23f4,
    0xb9d9a120dcc1ca03, 0x2e15af8165234a2e, 0x10609ba2720573d4, 0xaa4191b60368d1d5,
    0x333dd2300bc57762, 0xdf6ec48f79fb402f, 0x5ed20fcef1b734fa, 0x4c94924ec8be21ee,
    0x5abe6ad9d131e631, 0xbe10136a522e602d, 0x53671115c340e779, 0x9f392fe43e2144da,
};

/// Hash one 64-bit generator state as two 32-bit words.
fn pairHash(comptime index: usize, a: u32, b: u32) u64 {
    return (@as(u64, a) +% constants[2 * index]) *%
        (@as(u64, b) +% constants[2 * index + 1]);
}

/// Deterministic stream reproducing the pinned HiGHS `HighsRandom` sequence.
pub const RandomStream = struct {
    /// Nonzero xorshift state.
    state: u64 = 0,

    /// Seed and mix the stream; seed zero is the simplex compatibility path.
    pub fn init(seed: u32) RandomStream {
        var result = RandomStream{ .state = seed };
        while (true) {
            result.state = pairHash(0, @truncate(result.state), @truncate(result.state >> 32));
            result.state ^= pairHash(1, @truncate(result.state >> 32), seed) >> 32;
            if (result.state != 0) return result;
        }
    }

    /// Advance the underlying xorshift64 state once.
    fn advance(self: *RandomStream) void {
        self.state ^= self.state >> 12;
        self.state ^= self.state << 25;
        self.state ^= self.state >> 27;
    }

    /// Extract one hashed candidate using the selected HiGHS hash lane.
    fn candidate(self: *const RandomStream, comptime index: usize, shift: u6) u32 {
        return @truncate(pairHash(index, @truncate(self.state), @truncate(self.state >> 32)) >> shift);
    }

    /// Draw uniformly from `[0, supremum)` using rejection sampling.
    fn drawUniform(self: *RandomStream, supremum: u32) u32 {
        if (supremum <= 1) return 0;
        const bits: u6 = @intCast(32 - @clz(supremum - 1));
        const shift: u6 = @intCast(64 - @as(u7, bits));
        const indices = [_]u6{ 0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31 };
        while (true) {
            self.advance();
            inline for (indices) |index| {
                const value = self.candidate(index, shift);
                if (value < supremum) return value;
            }
        }
    }

    /// Public machine-sized uniform draw. `supremum <= 1` returns zero.
    pub fn integer(self: *RandomStream, supremum: usize) usize {
        if (supremum <= 1) return 0;
        std.debug.assert(supremum <= std.math.maxInt(u32));
        return self.drawUniform(@intCast(supremum));
    }

    /// Advance exactly as Fisher–Yates would without retaining its output.
    fn consumeShuffle(self: *RandomStream, count: usize) void {
        var index = count;
        while (index > 1) : (index -= 1) _ = self.drawUniform(@intCast(index));
    }

    /// Shuffle `values` in place with the pinned Fisher–Yates sequence.
    fn shuffle(self: *RandomStream, values: []u32) void {
        var index = values.len;
        while (index > 1) : (index -= 1) {
            const position: usize = self.drawUniform(@intCast(index));
            std.mem.swap(u32, &values[position], &values[index - 1]);
        }
    }

    /// Return the next strictly positive binary fraction below one.
    pub fn fraction(self: *RandomStream) f64 {
        self.advance();
        const low: u32 = @truncate(self.state);
        const high: u32 = @truncate(self.state >> 32);
        const output = (pairHash(0, low, high) >> 12) ^ (pairHash(1, low, high) >> 38);
        return @as(f64, @floatFromInt(1 + output)) * 2.2204460492503125e-16;
    }

    /// State immediately after `HEkk::initialiseSimplexLpRandomVectors`.
    /// A cold solve calls that routine once from `initialiseEkk` and again
    /// from `initialiseForSolve`, without reinitialising `random_` between
    /// calls. `HEkkDual::correctDualInfeasibilities` continues after the
    /// second call.
    pub fn afterVectorInitialization(num_cols: usize, num_tot: usize) RandomStream {
        var random = RandomStream.init(0);
        for (0..2) |_| {
            random.consumeShuffle(num_cols);
            random.consumeShuffle(num_tot);
            for (0..num_tot) |_| _ = random.fraction();
        }
        return random;
    }
};

/// Reproduce the second `HEkk::initialiseSimplexLpRandomVectors` call used by
/// a cold solve for random_seed=0. The first call advances the same RNG and
/// its vectors are overwritten before dual simplex starts.
pub fn initializeVectors(random_values: []f64, permutation: []u32, num_cols: usize) void {
    std.debug.assert(random_values.len == permutation.len);
    var random = RandomStream.init(0);
    random.consumeShuffle(num_cols);
    random.consumeShuffle(permutation.len);
    for (permutation) |_| _ = random.fraction();

    random.consumeShuffle(num_cols);
    for (permutation, 0..) |*value, index| value.* = @intCast(index);
    random.shuffle(permutation);
    for (random_values) |*value| value.* = random.fraction();
}

test "HiGHS random vectors are deterministic and bounded" {
    var first_values: [7]f64 = undefined;
    var first_permutation: [7]u32 = undefined;
    var second_values: [7]f64 = undefined;
    var second_permutation: [7]u32 = undefined;
    initializeVectors(&first_values, &first_permutation, 3);
    initializeVectors(&second_values, &second_permutation, 3);
    try std.testing.expectEqualSlices(f64, &first_values, &second_values);
    try std.testing.expectEqualSlices(u32, &first_permutation, &second_permutation);
    for (first_values) |value| try std.testing.expect(value > 0.0 and value < 1.0);
}

test "cold-solve second random vectors match pinned HiGHS" {
    var values: [1088]f64 = undefined;
    var permutation: [1088]u32 = undefined;
    initializeVectors(&values, &permutation, 688);
    try std.testing.expectEqual(@as(u32, 308), permutation[0]);
    try std.testing.expectEqual(@as(u32, 1058), permutation[1]);
    try std.testing.expectEqual(@as(f64, 0.20631502720009948), values[0]);
    try std.testing.expectEqual(@as(f64, 0.7775873865060616), values[1]);
}
