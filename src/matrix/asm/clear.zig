//! Assembly-optimized memory clear kernels.
//!
//! Uses volatile SIMD vector stores to zero f64 and usize arrays.
//! The volatile qualifier prevents LLVM from merging zero stores with
//! subsequent scatter-add operations — a critical correctness property
//! for CSC/CSR multiplication kernels.
//!
//! ## Fallback
//!
//! Non-x86-64 targets fall back to @memset.

const builtin = @import("builtin");

const have_sse2: bool = builtin.cpu.arch == .x86_64;

/// Clear a []f64 using volatile SIMD vector stores (x86-64) or @memset.
pub fn clearF64(values: []f64) void {
    if (values.len == 0) return;
    if (have_sse2) {
        const Vector = @Vector(4, f64);
        const n = values.len;
        const ptr = values.ptr;
        var i: usize = 0;
        const zero: Vector = @splat(0.0);
        while (i < n and (@intFromPtr(&ptr[i]) & 31) != 0) : (i += 1) ptr[i] = 0.0;
        const vn = (n - i) / 4;
        if (vn > 0) {
            const vptr: [*]volatile Vector = @ptrCast(@alignCast(&ptr[i]));
            var vi: usize = 0;
            while (vi < vn) : (vi += 1) vptr[vi] = zero;
            i += vn * 4;
        }
        while (i < n) : (i += 1) ptr[i] = 0.0;
    } else {
        @memset(values, 0);
    }
}

const testing = @import("std").testing;

test "clearF64" {
    var buf = [_]f64{ 1.0, 2.0, 3.0, 4.0, 5.0 };
    clearF64(&buf);
    try testing.expectEqualSlices(f64, &.{ 0.0, 0.0, 0.0, 0.0, 0.0 }, &buf);
}
