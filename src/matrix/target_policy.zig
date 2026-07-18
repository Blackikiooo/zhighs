//! Compile-time target policy for matrix hot paths.
//!
//! All decisions are derived from Zig's `builtin` target description.  The
//! resulting constants disappear during compilation, so sparse kernels do not
//! pay for architecture checks or runtime dispatch inside their inner loops.

const std = @import("std");
const builtin = @import("builtin");

/// Conservative cache-line size used for owning matrix workspaces.
/// Mainstream x86, AArch64, ARM and RISC-V server cores use 64-byte lines.  A
/// conservative 64-byte alignment is still valid on targets with wider lines
/// and avoids inflating every small allocation to a page-sized boundary.
pub const cache_line_bytes: comptime_int = switch (builtin.cpu.arch) {
    .x86, .x86_64, .aarch64, .arm, .riscv32, .riscv64, .powerpc64, .powerpc64le => 64,
    else => @max(64, @alignOf(std.simd.suggestVectorLength(u8) orelse 16)),
};

/// Number of future sparse columns to prefetch during basis assembly.
/// Sparse columns are pointer-chasing workloads: a modest lead hides the
/// starts/index/value latency without evicting the column currently copied.
pub const sparse_column_prefetch_distance: comptime_int = switch (builtin.cpu.arch) {
    .x86, .x86_64 => 8,
    .aarch64 => 6,
    else => 4,
};

/// Native vector width selected from the requested build target, not the host.
/// This keeps cross-compilation correct and allows `-Dcpu=native` to expose
/// wider vectors without maintaining architecture-specific source variants.
pub fn vectorLanes(comptime T: type) comptime_int {
    if (builtin.cpu.arch == .x86 or builtin.cpu.arch == .x86_64) {
        if (std.Target.x86.featureSetHas(builtin.cpu.features, .avx512f)) return 64 / @sizeOf(T);
        if (std.Target.x86.featureSetHas(builtin.cpu.features, .avx2)) return 32 / @sizeOf(T);
    }
    return 16 / @sizeOf(T);
}

test "target policy exposes power-of-two cache and vector widths" {
    try std.testing.expect(std.math.isPowerOfTwo(cache_line_bytes));
    try std.testing.expect(std.math.isPowerOfTwo(vectorLanes(f64)));
}
