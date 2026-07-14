//! Allocator experiments for matrix-specific lifetimes.
//! Measures short-lived HUInt scratch and session-scoped builder/output data.

const std = @import("std");
const zhighs = @import("zhighs");

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

fn scratchHeap(allocator: std.mem.Allocator, count: usize, repeats: usize) !u64 {
    const started = nowNs();
    for (0..repeats) |_| {
        const scratch = try allocator.alloc(zhighs.HUInt, count);
        @memset(scratch, 0);
        clobber(scratch.ptr);
        allocator.free(scratch);
    }
    return @intCast(@divTrunc(nowNs() - started, repeats));
}

fn scratchStackFallback(allocator: std.mem.Allocator, count: usize, repeats: usize) !u64 {
    const started = nowNs();
    for (0..repeats) |_| {
        var state = std.heap.stackFallback(4096, allocator);
        const scratch_allocator = state.get();
        const scratch = try scratch_allocator.alloc(zhighs.HUInt, count);
        @memset(scratch, 0);
        clobber(scratch.ptr);
        scratch_allocator.free(scratch);
    }
    return @intCast(@divTrunc(nowNs() - started, repeats));
}

fn fillBuilder(builder: *zhighs.matrix.MatrixBuilder, allocator: std.mem.Allocator, dimension: usize) !void {
    try builder.reserve(allocator, dimension * 3 - 2);
    for (0..dimension) |col| {
        const col_id = zhighs.ColId.fromUsizeAssumeValid(col);
        if (col != 0) builder.appendPreReserved(zhighs.RowId.fromUsizeAssumeValid(col - 1), col_id, -1.0);
        builder.appendPreReserved(zhighs.RowId.fromUsizeAssumeValid(col), col_id, 4.0);
        if (col + 1 < dimension) builder.appendPreReserved(zhighs.RowId.fromUsizeAssumeValid(col + 1), col_id, -1.0);
    }
}

fn buildSmp(allocator: std.mem.Allocator, dimension: usize, repeats: usize) !u64 {
    const started = nowNs();
    for (0..repeats) |_| {
        var builder = try zhighs.matrix.MatrixBuilder.init(dimension, dimension);
        try fillBuilder(&builder, allocator, dimension);
        var matrix = try builder.freezeSortedAssumeValid(allocator, 0.0);
        clobber(matrix.values.ptr);
        matrix.deinit(allocator);
        builder.deinit(allocator);
    }
    return @intCast(@divTrunc(nowNs() - started, repeats));
}

fn buildArena(child: std.mem.Allocator, dimension: usize, repeats: usize) !struct { ns: u64, capacity: usize } {
    var arena = std.heap.ArenaAllocator.init(child);
    defer arena.deinit();
    const allocator = arena.allocator();
    // Warm once so retain_capacity measures the intended session reuse case.
    {
        var builder = try zhighs.matrix.MatrixBuilder.init(dimension, dimension);
        try fillBuilder(&builder, allocator, dimension);
        const matrix = try builder.freezeSortedAssumeValid(allocator, 0.0);
        clobber(matrix.values.ptr);
        _ = arena.reset(.retain_capacity);
    }
    const started = nowNs();
    for (0..repeats) |_| {
        var builder = try zhighs.matrix.MatrixBuilder.init(dimension, dimension);
        try fillBuilder(&builder, allocator, dimension);
        const matrix = try builder.freezeSortedAssumeValid(allocator, 0.0);
        clobber(matrix.values.ptr);
        _ = arena.reset(.retain_capacity);
    }
    return .{ .ns = @intCast(@divTrunc(nowNs() - started, repeats)), .capacity = arena.queryCapacity() };
}

pub fn main() !void {
    const allocator = std.heap.smp_allocator;
    for ([_]usize{ 64, 256, 1024, 4096 }) |count| {
        const repeats: usize = 100_000;
        const heap_ns = try scratchHeap(allocator, count, repeats);
        const stack_ns = try scratchStackFallback(allocator, count, repeats);
        std.debug.print("scratch\tcount={d}\theap_ns={d}\tstack_fallback_ns={d}\tratio={d:.4}\n", .{
            count,
            heap_ns,
            stack_ns,
            @as(f64, @floatFromInt(stack_ns)) / @as(f64, @floatFromInt(heap_ns)),
        });
    }
    for ([_]usize{ 64, 512, 4096 }) |dimension| {
        const repeats: usize = if (dimension < 4096) 4000 else 500;
        const smp_ns = try buildSmp(allocator, dimension, repeats);
        const arena_result = try buildArena(allocator, dimension, repeats);
        const page_repeats: usize = if (dimension < 4096) 200 else 40;
        const page_ns = try buildSmp(std.heap.page_allocator, dimension, page_repeats);
        std.debug.print("sorted_build\tdimension={d}\tsmp_ns={d}\tpage_ns={d}\tarena_ns={d}\tarena_ratio={d:.4}\tpage_ratio={d:.4}\tarena_capacity={d}\n", .{
            dimension,
            smp_ns,
            page_ns,
            arena_result.ns,
            @as(f64, @floatFromInt(arena_result.ns)) / @as(f64, @floatFromInt(smp_ns)),
            @as(f64, @floatFromInt(page_ns)) / @as(f64, @floatFromInt(smp_ns)),
            arena_result.capacity,
        });
    }
}
