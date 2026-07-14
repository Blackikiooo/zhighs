//! Single-variable experiment for manual @Vector column scaling.
//! Compares the ordinary scalar Zig loop (which LLVM may auto-vectorize) with
//! an explicit four-lane f64 implementation over identical CSC column spans.

const std = @import("std");

inline fn clobber(pointer: anytype) void {
    asm volatile (""
        :
        : [pointer] "r" (pointer),
        : .{ .memory = true });
}

fn nowNs() i128 {
    var ts: std.posix.timespec = undefined;
    switch (std.posix.errno(std.posix.system.clock_gettime(.MONOTONIC, &ts))) {
        .SUCCESS => {},
        else => unreachable,
    }
    return @as(i128, ts.sec) * std.time.ns_per_s + ts.nsec;
}

noinline fn scaleScalar(starts: []const u32, values: []f64, factors: []const f64) void {
    for (factors, 0..) |factor, col| {
        var position: usize = starts[col];
        const end: usize = starts[col + 1];
        while (position < end) : (position += 1) values[position] *= factor;
    }
}

noinline fn scaleVector(starts: []const u32, values: []f64, factors: []const f64) void {
    const Vec = @Vector(4, f64);
    for (factors, 0..) |factor, col| {
        var position: usize = starts[col];
        const end: usize = starts[col + 1];
        const factor_vector: Vec = @splat(factor);
        while (position + 4 <= end) : (position += 4) {
            const pointer: *align(1) Vec = @ptrCast(values.ptr + position);
            pointer.* *= factor_vector;
        }
        while (position < end) : (position += 1) values[position] *= factor;
    }
}

fn checksum(values: []const f64) f64 {
    var result: f64 = 0.0;
    for (values) |value| result += value;
    return result;
}

const Range = struct { min: f64, max: f64 };

noinline fn rangeScalar(values: []const f64) Range {
    var result: Range = .{ .min = std.math.inf(f64), .max = 0.0 };
    for (values) |value| {
        const magnitude = @abs(value);
        result.min = @min(result.min, magnitude);
        result.max = @max(result.max, magnitude);
    }
    return result;
}

noinline fn rangeVector(values: []const f64) Range {
    const Vec = @Vector(4, f64);
    var minima: Vec = @splat(std.math.inf(f64));
    var maxima: Vec = @splat(0.0);
    var position: usize = 0;
    while (position + 4 <= values.len) : (position += 4) {
        const pointer: *align(1) const Vec = @ptrCast(values.ptr + position);
        const magnitudes = @abs(pointer.*);
        minima = @min(minima, magnitudes);
        maxima = @max(maxima, magnitudes);
    }
    var result: Range = .{ .min = @reduce(.Min, minima), .max = @reduce(.Max, maxima) };
    while (position < values.len) : (position += 1) {
        const magnitude = @abs(values[position]);
        result.min = @min(result.min, magnitude);
        result.max = @max(result.max, magnitude);
    }
    return result;
}

fn benchmarkRange(allocator: std.mem.Allocator, length: usize) !void {
    const values = try allocator.alloc(f64, length);
    defer allocator.free(values);
    for (values, 0..) |*value, index| {
        const magnitude = 0.001 + @as(f64, @floatFromInt(index % 1009)) * 0.125;
        value.* = if (index & 1 == 0) magnitude else -magnitude;
    }
    const repeats = std.math.clamp(200_000_000 / @max(length, 1), 20, 200_000);
    var scalar_result: Range = undefined;
    var started = nowNs();
    for (0..repeats) |_| {
        scalar_result = rangeScalar(values);
        clobber(&scalar_result);
    }
    const scalar_ns: u64 = @intCast(@divTrunc(nowNs() - started, repeats));
    var vector_result: Range = undefined;
    started = nowNs();
    for (0..repeats) |_| {
        vector_result = rangeVector(values);
        clobber(&vector_result);
    }
    const vector_ns: u64 = @intCast(@divTrunc(nowNs() - started, repeats));
    std.debug.print("absolute_range\tlength={d}\tscalar_ns={d}\tvector_ns={d}\tratio={d:.4}\tresults={d:.6},{d:.6}/{d:.6},{d:.6}\n", .{
        length,
        scalar_ns,
        vector_ns,
        @as(f64, @floatFromInt(vector_ns)) / @as(f64, @floatFromInt(scalar_ns)),
        scalar_result.min,
        scalar_result.max,
        vector_result.min,
        vector_result.max,
    });
}

fn benchmark(allocator: std.mem.Allocator, nnz_per_col: usize) !void {
    const columns: usize = 16_384;
    const nnz = columns * nnz_per_col;
    const starts = try allocator.alloc(u32, columns + 1);
    defer allocator.free(starts);
    const scalar_values = try allocator.alloc(f64, nnz);
    defer allocator.free(scalar_values);
    const vector_values = try allocator.alloc(f64, nnz);
    defer allocator.free(vector_values);
    const factors = try allocator.alloc(f64, columns);
    defer allocator.free(factors);

    for (starts, 0..) |*start, index| start.* = @intCast(index * nnz_per_col);
    for (scalar_values, vector_values, 0..) |*scalar, *vector, index| {
        const value = 1.0 + @as(f64, @floatFromInt(index % 31)) * 0.015625;
        scalar.* = value;
        vector.* = value;
    }
    for (factors, 0..) |*factor, index| factor.* = if (index & 1 == 0) 0.999 else 1.001;

    const repeats = std.math.clamp(200_000_000 / @max(nnz, 1), 20, 20_000);
    scaleScalar(starts, scalar_values, factors);
    var started = nowNs();
    for (0..repeats) |_| {
        scaleScalar(starts, scalar_values, factors);
        clobber(scalar_values.ptr);
    }
    const scalar_ns: u64 = @intCast(@divTrunc(nowNs() - started, repeats));

    scaleVector(starts, vector_values, factors);
    started = nowNs();
    for (0..repeats) |_| {
        scaleVector(starts, vector_values, factors);
        clobber(vector_values.ptr);
    }
    const vector_ns: u64 = @intCast(@divTrunc(nowNs() - started, repeats));

    std.debug.print("column_scale\tnnz_per_col={d}\tnnz={d}\tscalar_ns={d}\tvector_ns={d}\tratio={d:.4}\tchecksums={d:.6}/{d:.6}\n", .{
        nnz_per_col,
        nnz,
        scalar_ns,
        vector_ns,
        @as(f64, @floatFromInt(vector_ns)) / @as(f64, @floatFromInt(scalar_ns)),
        checksum(scalar_values),
        checksum(vector_values),
    });
}

pub fn main() !void {
    const allocator = std.heap.smp_allocator;
    for ([_]usize{ 3, 8, 16, 32, 128 }) |nnz_per_col| try benchmark(allocator, nnz_per_col);
    for ([_]usize{ 31, 4096, 131_072, 2_097_152 }) |length| try benchmarkRange(allocator, length);
}
