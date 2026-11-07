//! Compile-time target layout policy for matrix hot paths.
//!
//! Role: layout, alignment, vector-width, and unroll constants.
//! Shape: none; this module holds no matrix data.
//!
//! All decisions are derived from Zig's `builtin` target description.  The
//! resulting constants disappear during compilation, so sparse kernels do not
//! pay for architecture checks or runtime dispatch inside their inner loops.
//!
//! Model-level numeric policy does not belong here: the shared LP coefficient
//! floor lives in `foundation.tolerance`.

const std = @import("std");
const builtin = @import("builtin");

/// Cache-line size in bytes for owning matrix workspaces.
///
/// x86, non-Apple AArch64, ARM, RISC-V, s390x and LoongArch cores use 64-byte
/// lines.  Apple Silicon and POWER use 128-byte lines: on those targets two
/// unrelated owning blocks packed into one line can invalidate each other when
/// a solver keeps several matrices live in the same iteration.
pub const cache_line_bytes: comptime_int = switch (builtin.cpu.arch) {
    .aarch64, .aarch64_be => if (isAppleCpu()) 128 else 64,
    .powerpc, .powerpcle, .powerpc64, .powerpc64le => 128,
    else => 64,
};

/// True for any Apple-designed core (M1 and later report `apple_*` models).
fn isAppleCpu() bool {
    return std.mem.startsWith(u8, builtin.cpu.model.name, "apple");
}

/// Number of future sparse columns to prefetch during basis assembly.
/// Sparse columns are pointer-chasing workloads: a modest lead hides the
/// starts/index/value latency without evicting the column currently copied.
pub const sparse_column_prefetch_distance: comptime_int = switch (builtin.cpu.arch) {
    .x86, .x86_64 => 8,
    .aarch64, .aarch64_be => 6,
    else => 4,
};

/// Number of `T` elements in one native vector register for the requested
/// build target (never the host, so cross-compilation stays correct).
///
/// `std.simd.suggestVectorLength` already encodes the preferences that matter
/// for these kernels, so it is used as the single source of truth:
///   * x86 honours `prefer_256_bit` / `prefer_128_bit`.  Client Tiger Lake and
///     Ice Lake cores expose `avx512f` but down-clock on 512-bit traffic, so
///     they get 256-bit vectors instead of the widest encoding;
///   * AArch64 distinguishes NEON (128-bit) from SVE (256-bit default), which
///     a plain `16 / @sizeOf(T)` fallback cannot express;
///   * POWER picks VSX, RISC-V picks its `Zvl*` vector length, s390x and
///     LoongArch pick their own register widths.
///
/// The result is always a power of two so callers may use it as an alignment
/// stride, and it is never zero.
pub fn vectorLanes(comptime T: type) comptime_int {
    comptime {
        const suggested = std.simd.suggestVectorLength(T) orelse fallbackLanes(T);
        return @max(1, std.math.floorPowerOfTwo(usize, suggested));
    }
}

/// Vector register width in bytes for `T`.  Layout code that wants to keep
/// several streams on distinct vector boundaries can compare this against
/// `cache_line_bytes`.
pub fn vectorBytes(comptime T: type) comptime_int {
    return vectorLanes(T) * @sizeOf(T);
}

/// Last-resort 128-bit vector for targets Zig does not model explicitly.
fn fallbackLanes(comptime T: type) comptime_int {
    const element_bits: usize = @max(8, @bitSizeOf(T));
    return @max(1, 128 / element_bits);
}

/// Tuned unroll factor for vector store loops (e.g. memory clear kernels).
///
/// Larger values benefit cores with deep store buffers / write combining;
/// smaller values avoid uop-cache pressure on narrower pipelines.
/// Current policy:
///   x86/x86_64 — AMD Zen family ("znver*") → 4,
///                wide-vector parts without a 256-bit preference → 3,
///                everything else → 2
///   AArch64    — Apple Silicon → 4, server/device cores → 3
///   PowerPC    — 3
///   others     — 2 (conservative default)
///
/// Zen detection must match `znver1` … `znver5`: Zig's model names never
/// contain the substring `zen`, so a plain `"zen"` search silently disables
/// the wide unroll on every AMD part.
pub fn unrollFactor() comptime_int {
    const cpu = builtin.cpu;
    return switch (cpu.arch) {
        .x86, .x86_64 => if (std.mem.indexOf(u8, cpu.model.name, "znver") != null)
            4
        else if (std.Target.x86.featureSetHas(cpu.features, .avx512f) and
            !std.Target.x86.featureSetHas(cpu.features, .prefer_256_bit))
            3
        else
            2,
        .aarch64, .aarch64_be => if (isAppleCpu()) 4 else 3,
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => 3,
        else => 2,
    };
}

test "target policy exposes power-of-two cache and vector widths" {
    try std.testing.expect(std.math.isPowerOfTwo(cache_line_bytes));

    inline for (.{ f64, f32, u8, u16, u32, u64, usize, i32 }) |T| {
        const lanes = vectorLanes(T);
        try std.testing.expect(lanes >= 1);
        try std.testing.expect(std.math.isPowerOfTwo(lanes));
        // Zig never suggests more than one 1024-bit register (Hexagon HVX).
        try std.testing.expect(vectorBytes(T) <= 128);
        // A vector must hold at least one element of the requested type.
        try std.testing.expect(vectorBytes(T) >= @sizeOf(T));
    }
}

test "target policy unroll factor stays in the tuned range" {
    const unroll = unrollFactor();
    try std.testing.expect(unroll >= 1);
    try std.testing.expect(unroll <= 8);
    // Every tuned value must keep the store loop's cursor arithmetic exact.
    try std.testing.expect(std.math.isPowerOfTwo(unroll) or unroll == 3);
}

test "vector lanes agree with the Zig suggestion for f64" {
    const suggested = std.simd.suggestVectorLength(f64) orelse return;
    try std.testing.expectEqual(
        @max(1, std.math.floorPowerOfTwo(usize, suggested)),
        vectorLanes(f64),
    );
}
