//! Small memory kernels shared by sparse matrix hot paths.

const std = @import("std");

/// Clears f64 storage through its byte representation. Keeping this operation
/// centralized makes target-specific replacement possible if Zig's generated
/// memset regresses. IEEE-754 +0.0 has the all-zero bit pattern.
pub inline fn clearF64Bytes(values: []f64) void {
    @memset(std.mem.sliceAsBytes(values), 0);
}

/// Explicit vector stores avoid a slow scalar memset lowering observed with
/// Zig 0.16 on some x86-64 targets. Unaligned vector pointers are intentional.
pub inline fn clearF64(values: []f64) void {
    const Vector = @Vector(4, f64);
    var scalar_index: usize = 0;
    while (scalar_index < values.len and (@intFromPtr(&values[scalar_index]) & (@alignOf(Vector) - 1)) != 0) : (scalar_index += 1)
        values[scalar_index] = 0.0;

    const remaining = values[scalar_index..];
    const vector_count = remaining.len / 4;
    const zero: Vector = @splat(0.0);
    if (vector_count != 0) {
        const vectors: [*]volatile Vector = @ptrCast(@alignCast(remaining.ptr));
        for (0..vector_count) |index| vectors[index] = zero;
    }
    const tail: [*]volatile f64 = @ptrCast(remaining[vector_count * 4 ..].ptr);
    for (0..remaining.len - vector_count * 4) |index| tail[index] = 0.0;
}

pub inline fn clearUsize(values: []usize) void {
    const Vector = @Vector(4, usize);
    var scalar_index: usize = 0;
    while (scalar_index < values.len and (@intFromPtr(&values[scalar_index]) & (@alignOf(Vector) - 1)) != 0) : (scalar_index += 1)
        values[scalar_index] = 0;
    const remaining = values[scalar_index..];
    const vector_count = remaining.len / 4;
    const zero: Vector = @splat(0);
    if (vector_count != 0) {
        const vectors: [*]volatile Vector = @ptrCast(@alignCast(remaining.ptr));
        for (0..vector_count) |index| vectors[index] = zero;
    }
    const tail: [*]volatile usize = @ptrCast(remaining[vector_count * 4 ..].ptr);
    for (0..remaining.len - vector_count * 4) |index| tail[index] = 0;
}

test "clear kernels handle short unaligned slices" {
    var floats = [_]f64{ 1.0, 2.0, 3.0, 4.0, 5.0 };
    clearF64(floats[1..4]);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 0.0, 0.0, 0.0, 5.0 }, &floats);
    var integers = [_]usize{ 1, 2, 3, 4, 5 };
    clearUsize(integers[1..4]);
    try std.testing.expectEqualSlices(usize, &.{ 1, 0, 0, 0, 5 }, &integers);
}
