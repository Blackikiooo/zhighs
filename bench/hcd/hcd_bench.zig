const std = @import("std");
const Foundation = @import("foundation");

const HCD = Foundation.HCD;
const HD = Foundation.HD;

const n: usize = 4096;
const repeats: usize = 20_000;

fn makeHD(i: usize) HD {
    const x: u64 = @as(u64, i) *% 6364136223846793005 +% 1442695040888963407;
    const mant = @as(HD, @floatFromInt((x >> 12) & ((@as(u64, 1) << 40) - 1))) / @as(HD, @floatFromInt(@as(u64, 1) << 40));
    const sign: HD = if ((x & 1) == 0) 1.0 else -1.0;
    const scale: HD = @as(HD, @floatFromInt((x >> 52) & 15)) * 0.0625 + 0.5;
    return sign * (mant + 0.125) * scale;
}

fn makeSmallHD(i: usize) HD {
    return makeHD(i) * 1e-12;
}

fn fillHD(values: []HD, comptime small: bool) void {
    for (values, 0..) |*v, i| {
        v.* = if (small) makeSmallHD(i + 1) else makeHD(i + 1);
    }
}

fn fillHCD(values: []HCD, comptime small: bool) void {
    for (values, 0..) |*v, i| {
        const hi = if (small) makeSmallHD(i + 1) else makeHD(i + 1);
        v.* = HCD.init(hi, makeSmallHD(i + 10_001));
    }
}

fn elapsedNs(start: i128, end: i128) u64 {
    return @intCast(end - start);
}

fn nowNs() i128 {
    var ts: std.posix.timespec = undefined;
    switch (std.posix.errno(std.posix.system.clock_gettime(.MONOTONIC, &ts))) {
        .SUCCESS => {},
        else => unreachable,
    }
    return @as(i128, ts.sec) * std.time.ns_per_s + ts.nsec;
}

fn printResult(name: []const u8, ns: u64, checksum: HD) void {
    const ops = @as(HD, @floatFromInt(n * repeats));
    const ns_per_op = @as(HD, @floatFromInt(ns)) / ops;
    std.debug.print("{s},{d},{d:.6},{d:.17}\n", .{ name, ns, ns_per_op, checksum });
}

fn benchAddHDFast(values: []const HD) !void {
    var x = HCD.init(1.0, 1e-16);
    const start = nowNs();
    for (0..repeats) |_| {
        for (values) |v| x = x.addHDFast(v);
    }
    const end = nowNs();
    std.mem.doNotOptimizeAway(x);
    printResult("zig.add_hd_fast", elapsedNs(start, end), x.toHD());
}

fn benchAddHDOrderedFast(values: []const HD) !void {
    var x = HCD.init(1.0e8, 1e-16);
    const start = nowNs();
    for (0..repeats) |_| {
        for (values) |v| x = x.addHDOrderedFast(v);
    }
    const end = nowNs();
    std.mem.doNotOptimizeAway(x);
    printResult("zig.add_hd_ordered_fast", elapsedNs(start, end), x.toHD());
}

fn benchAddHCDFast(values: []const HCD) !void {
    var x = HCD.init(1.0, 1e-16);
    const start = nowNs();
    for (0..repeats) |_| {
        for (values) |v| x = x.addHCDFast(v);
    }
    const end = nowNs();
    std.mem.doNotOptimizeAway(x);
    printResult("zig.add_hcd_fast", elapsedNs(start, end), x.toHD());
}

fn benchMultiplyHD(values: []const HD) !void {
    var checksum: HD = 0.0;
    const start = nowNs();
    for (0..repeats) |r| {
        var x = HCD.init(1.0000001 + @as(HD, @floatFromInt(r)) * 1e-18, 1e-16);
        for (values) |v| x = x.multiplyHD(1.0 + @abs(v) * 1e-12);
        checksum += x.toHD();
    }
    const end = nowNs();
    std.mem.doNotOptimizeAway(checksum);
    printResult("zig.multiply_hd", elapsedNs(start, end), checksum);
}

fn benchMultiplyHCD(values: []const HCD) !void {
    var checksum: HD = 0.0;
    const start = nowNs();
    for (0..repeats) |r| {
        var x = HCD.init(1.0000001 + @as(HD, @floatFromInt(r)) * 1e-18, 1e-16);
        for (values) |v| x = x.multiplyHCD(HCD.init(1.0 + @abs(v.high) * 1e-12, @abs(v.low) * 1e-12));
        checksum += x.toHD();
    }
    const end = nowNs();
    std.mem.doNotOptimizeAway(checksum);
    printResult("zig.multiply_hcd", elapsedNs(start, end), checksum);
}

fn benchDivideHD(values: []const HD) !void {
    var checksum: HD = 0.0;
    const start = nowNs();
    for (0..repeats) |r| {
        var x = HCD.init(1.0000001 + @as(HD, @floatFromInt(r)) * 1e-18, 1e-16);
        for (values) |v| x = x.divideHD(1.0 + @abs(v) * 1e-12);
        checksum += x.toHD();
    }
    const end = nowNs();
    std.mem.doNotOptimizeAway(checksum);
    printResult("zig.divide_hd", elapsedNs(start, end), checksum);
}

fn benchDivideHCD(values: []const HCD) !void {
    var checksum: HD = 0.0;
    const start = nowNs();
    for (0..repeats) |r| {
        var x = HCD.init(1.0000001 + @as(HD, @floatFromInt(r)) * 1e-18, 1e-16);
        for (values) |v| x = x.divideHCD(HCD.init(1.0 + @abs(v.high) * 1e-12, @abs(v.low) * 1e-12));
        checksum += x.toHD();
    }
    const end = nowNs();
    std.mem.doNotOptimizeAway(checksum);
    printResult("zig.divide_hcd", elapsedNs(start, end), checksum);
}

fn benchDivideHCDFast(values: []const HCD) !void {
    var checksum: HD = 0.0;
    const start = nowNs();
    for (0..repeats) |r| {
        var x = HCD.init(1.0000001 + @as(HD, @floatFromInt(r)) * 1e-18, 1e-16);
        for (values) |v| x = x.divideHCDFast(HCD.init(1.0 + @abs(v.high) * 1e-12, @abs(v.low) * 1e-12));
        checksum += x.toHD();
    }
    const end = nowNs();
    std.mem.doNotOptimizeAway(checksum);
    printResult("zig.divide_hcd_fast", elapsedNs(start, end), checksum);
}

pub fn main() !void {
    var hd_values: [n]HD = undefined;
    var small_hd_values: [n]HD = undefined;
    var hcd_values: [n]HCD = undefined;

    fillHD(&hd_values, false);
    fillHD(&small_hd_values, true);
    fillHCD(&hcd_values, false);

    std.debug.print("name,total_ns,ns_per_op,checksum\n", .{});
    try benchAddHDFast(&hd_values);
    try benchAddHDOrderedFast(&small_hd_values);
    try benchAddHCDFast(&hcd_values);
    try benchMultiplyHD(&hd_values);
    try benchMultiplyHCD(&hcd_values);
    try benchDivideHD(&hd_values);
    try benchDivideHCD(&hcd_values);
    try benchDivideHCDFast(&hcd_values);
}
