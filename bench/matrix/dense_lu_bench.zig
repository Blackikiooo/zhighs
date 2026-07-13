//! Independent dense-basis FTRAN/BTRAN benchmark.
//!
//! Default ReleaseFast sweep:
//!   zig build bench-dense-lu -Doptimize=ReleaseFast -Dcpu=native
//!
//! Single-kernel mode for Linux perf:
//!   ZHIGHS_LU_KERNEL=btran ZHIGHS_LU_DIMENSION=512 \
//!     ZHIGHS_LU_REPEATS=1000 zig-out/bin/dense-lu-bench

const std = @import("std");
const zhighs = @import("zhighs");

const Kernel = enum { ftran, btran };
const default_dimensions = [_]usize{ 32, 64, 128, 256, 512 };

fn nowNs() i128 {
    var ts: std.posix.timespec = undefined;
    switch (std.posix.errno(std.posix.system.clock_gettime(.MONOTONIC, &ts))) {
        .SUCCESS => {},
        else => unreachable,
    }
    return @as(i128, ts.sec) * std.time.ns_per_s + ts.nsec;
}

inline fn clobber(pointer: anytype) void {
    asm volatile (""
        :
        : [pointer] "r" (pointer),
        : .{ .memory = true });
}

fn defaultRepeats(n: usize) usize {
    const coefficient_visits: usize = 64 * 1024 * 1024;
    const per_solve = n * n;
    return std.math.clamp(coefficient_visits / @max(per_solve, 1), 8, 20_000);
}

/// Dense, strictly diagonally dominant matrix with a cyclic row permutation.
/// The permutation forces partial pivoting while keeping every benchmark case
/// well-conditioned enough to measure solve traversal rather than failures.
fn fillPermutedBasis(data: []f64, n: usize) void {
    @memset(data, 0.0);
    for (0..n) |physical_row| {
        const source_row = (physical_row + 1) % n;
        var off_diagonal_sum: f64 = 0.0;
        for (0..n) |col| {
            if (col == source_row) continue;
            const code: i32 = @intCast((source_row * 17 + col * 29 + 11) % 23);
            const value = @as(f64, @floatFromInt(code - 11)) * 0.0078125;
            data[physical_row * n + col] = value;
            off_diagonal_sum += @abs(value);
        }
        data[physical_row * n + source_row] = off_diagonal_sum + 1.0;
    }
}

fn runKernel(lu: *zhighs.matrix.DenseLU, allocator: std.mem.Allocator, n: usize, kernel: Kernel, repeats: usize) !void {
    const seed = try allocator.alloc(f64, n);
    defer allocator.free(seed);
    const rhs = try allocator.alloc(f64, n);
    defer allocator.free(rhs);
    for (seed, 0..) |*value, index|
        value.* = @as(f64, @floatFromInt(@as(i64, @intCast(index % 31)) - 15)) * 0.03125 + 1.0;

    @memcpy(rhs, seed);
    switch (kernel) {
        .ftran => try lu.solve(rhs),
        .btran => try lu.solveTranspose(rhs),
    }

    const start = nowNs();
    var checksum: f64 = 0.0;
    for (0..repeats) |iteration| {
        @memcpy(rhs, seed);
        switch (kernel) {
            .ftran => try lu.solve(rhs),
            .btran => try lu.solveTranspose(rhs),
        }
        checksum += rhs[iteration % n];
        clobber(rhs.ptr);
    }
    const total_ns: u64 = @intCast(nowNs() - start);
    const ns_per_solve = @as(f64, @floatFromInt(total_ns)) / @as(f64, @floatFromInt(repeats));
    std.debug.print("zig,{s},{d},{d},{d},{d:.3},{d},{d:.17}\n", .{
        @tagName(kernel),
        n,
        repeats,
        total_ns,
        ns_per_solve,
        n * n * @sizeOf(f64),
        checksum,
    });
}

fn runDimension(allocator: std.mem.Allocator, n: usize, selected_kernel: ?Kernel, repeats_override: ?usize) !void {
    const coefficient_count = std.math.mul(usize, n, n) catch return error.InvalidDimension;
    const basis = try allocator.alloc(f64, coefficient_count);
    fillPermutedBasis(basis, n);
    var lu = zhighs.matrix.DenseLU.init(allocator);
    defer lu.deinit();
    try lu.factorizeOwned(n, basis);

    const repeats = repeats_override orelse defaultRepeats(n);
    if (selected_kernel) |kernel| {
        try runKernel(&lu, allocator, n, kernel, repeats);
    } else {
        try runKernel(&lu, allocator, n, .ftran, repeats);
        try runKernel(&lu, allocator, n, .btran, repeats);
    }
}

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;
    const selected_kernel: ?Kernel = if (init.environ_map.get("ZHIGHS_LU_KERNEL")) |name|
        if (std.mem.eql(u8, name, "ftran")) .ftran else if (std.mem.eql(u8, name, "btran")) .btran else return error.InvalidKernel
    else
        null;
    const selected_dimension: ?usize = if (init.environ_map.get("ZHIGHS_LU_DIMENSION")) |text|
        try std.fmt.parseUnsigned(usize, text, 10)
    else
        null;
    const repeats_override: ?usize = if (init.environ_map.get("ZHIGHS_LU_REPEATS")) |text|
        try std.fmt.parseUnsigned(usize, text, 10)
    else
        null;
    if (repeats_override) |repeats|
        if (repeats == 0) return error.InvalidRepeats;

    std.debug.print("implementation,kernel,dimension,repeats,total_ns,ns_per_solve,lu_bytes,checksum\n", .{});
    if (selected_dimension) |dimension| {
        if (dimension == 0) return error.InvalidDimension;
        try runDimension(allocator, dimension, selected_kernel, repeats_override);
    } else {
        for (default_dimensions) |dimension|
            try runDimension(allocator, dimension, selected_kernel, repeats_override);
    }
}
