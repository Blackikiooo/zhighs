//! Small memory kernels shared by sparse matrix hot paths.

const std = @import("std");

/// Generic vectorized clear for any integer or float type.
/// Uses @Vector(4, T) for aligned stores, volatile when specified.
/// This comptime function eliminates 3× copy-paste between f64 and usize variants.
fn clearGeneric(comptime T: type, values: []T, comptime volatile_stores: bool) void {
    if (values.len == 0) return;
    const Vector = @Vector(4, T);
    const align_mask = @alignOf(Vector) - 1;

    // Align: handle leading unaligned elements one at a time
    var i: usize = 0;
    while (i < values.len and (@intFromPtr(&values[i]) & align_mask) != 0) : (i += 1)
        values[i] = 0;

    // Vector-aligned stores (4 elements per iteration)
    const remaining = values[i..];
    const vn = remaining.len / 4;
    const zero: Vector = @splat(0);
    if (vn > 0) {
        if (comptime volatile_stores) {
            const vptr: [*]volatile Vector = @ptrCast(@alignCast(remaining.ptr));
            var vi: usize = 0;
            while (vi + 4 <= vn) : (vi += 4) {
                vptr[vi] = zero;
                vptr[vi + 1] = zero;
                vptr[vi + 2] = zero;
                vptr[vi + 3] = zero;
            }
            while (vi < vn) : (vi += 1) vptr[vi] = zero;
        } else {
            const vptr: [*]Vector = @ptrCast(@alignCast(remaining.ptr));
            var vi: usize = 0;
            while (vi + 4 <= vn) : (vi += 4) {
                vptr[vi] = zero;
                vptr[vi + 1] = zero;
                vptr[vi + 2] = zero;
                vptr[vi + 3] = zero;
            }
            while (vi < vn) : (vi += 1) vptr[vi] = zero;
        }
        i += vn * 4;
    }

    // Tail elements (past last full vector)
    if (comptime volatile_stores) {
        const tptr: [*]volatile T = @ptrCast(values.ptr + i);
        var ti: usize = i;
        while (ti < values.len) : (ti += 1) tptr[ti - i] = 0;
    } else {
        while (i < values.len) : (i += 1) values[i] = 0;
    }
}

/// Typed volatile SIMD clear for compact offset/index work arrays.
pub inline fn clearTyped(comptime T: type, values: []T) void {
    clearGeneric(T, values, true);
}

/// Explicit volatile SIMD vector stores for f64 arrays.
/// Volatile prevents LLVM from merging zero stores with subsequent
/// scatter-add operations — essential for correctness in CSC/CSR kernels.
pub inline fn clearF64(values: []f64) void {
    clearGeneric(f64, values, true);
}

/// Non-volatile clear for hot-path kernels where the compiler barrier
/// is provided externally (via asm "memory" clobber or @memset).
/// Benchmarking shows this causes severe regressions in CSC/CSR scatter
/// kernels because LLVM incorrectly elides zero stores. Reserved for
/// experimental use; most callers should use clearF64 (volatile).
pub inline fn clearF64Fast(values: []f64) void {
    clearGeneric(f64, values, false);
}

/// Byte-level clear via sliceAsBytes. Uses @memset which generates
/// rep stosb or SIMD stores on LLVM backends. Useful when the byte
/// representation is needed rather than typed vector stores.
pub inline fn clearF64Bytes(values: []f64) void {
    @memset(std.mem.sliceAsBytes(values), 0);
}

/// Volatile SIMD vector stores for usize arrays.
pub inline fn clearUsize(values: []usize) void {
    clearGeneric(usize, values, true);
}

/// Non-volatile SIMD vector stores for usize arrays.
/// Use this when the compiler barrier from volatile is not needed,
/// e.g. clearing row_starts before a counting pass.
pub inline fn clearUsizeFast(values: []usize) void {
    clearGeneric(usize, values, false);
}

test "clear kernels handle short unaligned slices" {
    var floats = [_]f64{ 1.0, 2.0, 3.0, 4.0, 5.0 };
    clearF64(floats[1..4]);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 0.0, 0.0, 0.0, 5.0 }, &floats);
    var integers = [_]usize{ 1, 2, 3, 4, 5 };
    clearUsize(integers[1..4]);
    try std.testing.expectEqualSlices(usize, &.{ 1, 0, 0, 0, 5 }, &integers);
}

test "clearGeneric matches clearF64 for f64" {
    var a: [64]f64 = undefined;
    var b: [64]f64 = undefined;
    for (&a, &b, 0..) |*va, *vb, i| {
        va.* = @floatFromInt(i + 1);
        vb.* = @floatFromInt(i + 1);
    }
    clearF64(&a);
    clearGeneric(f64, &b, true);
    try std.testing.expectEqualSlices(f64, &a, &b);
}

test "clearGeneric handles edge cases" {
    // Single element
    var single: [1]f64 = [_]f64{42.0};
    clearGeneric(f64, &single, true);
    try std.testing.expectEqual(@as(f64, 0.0), single[0]);

    // Exactly 4 elements (one vector)
    var four: [4]f64 = [_]f64{ 1, 2, 3, 4 };
    clearGeneric(f64, &four, true);
    try std.testing.expectEqual(@as(f64, 0.0), four[0] + four[1] + four[2] + four[3]);

    // Empty
    var empty: [0]f64 = undefined;
    clearGeneric(f64, &empty, true);
}
