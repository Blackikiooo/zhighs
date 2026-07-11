//! Sparse matrix microbenchmark. Run with:
//! zig build bench-matrix -Doptimize=ReleaseFast

const std = @import("std");
const zhighs = @import("zhighs");

const dimension: usize = 50_000;
const multiply_repeats: usize = 200;
const conversion_repeats: usize = 5;

fn nowNs() i128 {
    var ts: std.posix.timespec = undefined;
    switch (std.posix.errno(std.posix.system.clock_gettime(.MONOTONIC, &ts))) {
        .SUCCESS => {},
        else => unreachable,
    }
    return @as(i128, ts.sec) * std.time.ns_per_s + ts.nsec;
}

fn elapsed(start: i128, end: i128) u64 {
    return @intCast(end - start);
}

pub fn main() !void {
    const allocator = std.heap.page_allocator;
    var builder = try zhighs.matrix.MatrixBuilder.init(dimension, dimension);
    defer builder.deinit(allocator);
    try builder.reserve(allocator, dimension * 3 - 2);
    for (0..dimension) |index| {
        const row = try zhighs.RowId.fromUsize(index);
        try builder.append(allocator, row, try zhighs.ColId.fromUsize(index), 4.0);
        if (index != 0) try builder.append(allocator, row, try zhighs.ColId.fromUsize(index - 1), -1.0);
        if (index + 1 < dimension) try builder.append(allocator, row, try zhighs.ColId.fromUsize(index + 1), -1.0);
    }
    var matrix = try builder.freeze(allocator, 0.0);
    defer matrix.deinit(allocator);
    const x = try allocator.alloc(f64, dimension);
    defer allocator.free(x);
    const y = try allocator.alloc(f64, dimension);
    defer allocator.free(y);
    @memset(x, 1.0);

    var checksum: f64 = 0.0;
    var start = nowNs();
    for (0..multiply_repeats) |_| {
        matrix.multiplyAssumeValid(x, y);
        checksum += y[dimension / 2];
    }
    var end = nowNs();
    std.mem.doNotOptimizeAway(checksum);
    std.debug.print("csc_ax,{d},{d},{d:.3}\n", .{ dimension, elapsed(start, end), checksum });

    checksum = 0.0;
    start = nowNs();
    for (0..multiply_repeats) |_| {
        matrix.transposeMultiplyAssumeValid(x, y);
        checksum += y[dimension / 2];
    }
    end = nowNs();
    std.mem.doNotOptimizeAway(checksum);
    std.debug.print("csc_atx,{d},{d},{d:.3}\n", .{ dimension, elapsed(start, end), checksum });

    start = nowNs();
    for (0..conversion_repeats) |revision| {
        var cache = try zhighs.matrix.CsrCache.buildAssumeValid(allocator, matrix, revision);
        std.mem.doNotOptimizeAway(cache.values.ptr);
        cache.deinit(allocator);
    }
    end = nowNs();
    std.debug.print("csc_to_csr,{d},{d},{d}\n", .{ dimension, elapsed(start, end), conversion_repeats });

    start = nowNs();
    for (0..conversion_repeats) |_| {
        var transposed = try zhighs.matrix.transposeAssumeValid(allocator, matrix);
        std.mem.doNotOptimizeAway(transposed.values.ptr);
        transposed.deinit(allocator);
    }
    end = nowNs();
    std.debug.print("transpose,{d},{d},{d}\n", .{ dimension, elapsed(start, end), conversion_repeats });
}
