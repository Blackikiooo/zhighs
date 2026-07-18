//! Retained-workspace Forrest--Tomlin update and updated-solve benchmark.

const std = @import("std");
const zhighs = @import("zhighs");

fn nowNs() i128 {
    var ts: std.posix.timespec = undefined;
    _ = std.posix.system.clock_gettime(.MONOTONIC, &ts);
    return @as(i128, ts.sec) * std.time.ns_per_s + ts.nsec;
}

const ChainTiming = struct { aq_ns: u64 = 0, update_ns: u64 = 0, checksum: f64 = 0 };

fn runChain(lu: *zhighs.matrix.SparseLU, aq: []f64, chain_length: usize, timed: bool) !ChainTiming {
    const n = aq.len;
    var result = ChainTiming{};
    for (0..chain_length) |update_index| {
        const column = (update_index * 17 + 3) % n;
        @memset(aq, 0.0);
        aq[column] = 4.0 + @as(f64, @floatFromInt(column % 7)) * 0.125;
        aq[if (column == 0) n - 1 else column - 1] = -0.5;
        aq[if (column + 1 == n) 0 else column + 1] = -0.5;
        aq[(column + 7) % n] += 0.01;

        const aq_started = nowNs();
        try lu.solveForUpdate(aq);
        const update_started = nowNs();
        try lu.applyForrestTomlinUpdate(@intCast(column), aq, 1.0);
        const finished = nowNs();
        if (timed) {
            result.aq_ns += @intCast(update_started - aq_started);
            result.update_ns += @intCast(finished - update_started);
        }
        result.checksum += aq[column];
    }
    return result;
}

pub fn main(init: std.process.Init) !void {
    var args = std.process.Args.iterate(init.minimal.args);
    _ = args.next();
    const n = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 512;
    const chain_length = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 32;
    const solve_repeats = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 1001;
    if (n < 16 or chain_length == 0 or chain_length > n or solve_repeats == 0 or args.next() != null) return error.InvalidArguments;
    const allocator = std.heap.c_allocator;
    const starts = try allocator.alloc(zhighs.foundation.HUInt, n + 1);
    defer allocator.free(starts);
    const rows = try allocator.alloc(zhighs.foundation.RowId, n * 3);
    defer allocator.free(rows);
    const values = try allocator.alloc(f64, n * 3);
    defer allocator.free(values);
    var output: usize = 0;
    for (0..n) |column| {
        starts[column] = @intCast(output);
        var candidates = [_]usize{ if (column == 0) n - 1 else column - 1, column, if (column + 1 == n) 0 else column + 1 };
        std.mem.sort(usize, &candidates, {}, std.sort.asc(usize));
        for (candidates) |row| {
            rows[output] = zhighs.foundation.RowId.fromUsizeAssumeValid(row);
            values[output] = if (row == column) 4.0 + @as(f64, @floatFromInt(column % 7)) * 0.125 else -0.5;
            output += 1;
        }
    }
    starts[n] = @intCast(output);
    const basis = zhighs.matrix.SparseBasisView{ .dimension = n, .starts = starts, .rows = rows, .values = values };
    var lu = zhighs.matrix.SparseLU.init(allocator);
    defer lu.deinit();
    const aq = try allocator.alloc(f64, n);
    defer allocator.free(aq);
    const rhs = try allocator.alloc(f64, n);
    defer allocator.free(rhs);

    // First chain grows every retained FT stream outside the timed run.
    try lu.factorizeAssumeValid(basis);
    _ = try runChain(&lu, aq, chain_length, false);
    try lu.factorizeAssumeValid(basis);
    const timing = try runChain(&lu, aq, chain_length, true);

    var checksum = timing.checksum;
    const solve_started = nowNs();
    for (0..solve_repeats) |repeat| {
        for (rhs, 0..) |*value, index| value.* = 1.0 + @as(f64, @floatFromInt((index + repeat) % 13)) * 0.0625;
        try lu.solve(rhs);
        checksum += rhs[repeat % n];
        for (rhs, 0..) |*value, index| value.* = 0.5 + @as(f64, @floatFromInt((index * 3 + repeat) % 11)) * 0.125;
        try lu.solveTranspose(rhs);
        checksum += rhs[(repeat * 7) % n];
    }
    const solve_ns: u64 = @intCast(nowNs() - solve_started);
    std.debug.print("dimension,chain,aq_ns_per_update,ft_update_ns,updated_ftran_btran_ns,requested_bytes,checksum\n", .{});
    std.debug.print("{d},{d},{d},{d},{d},{d},{d:.17}\n", .{
        n,
        chain_length,
        timing.aq_ns / chain_length,
        timing.update_ns / chain_length,
        solve_ns / solve_repeats,
        lu.requestedBytes(),
        checksum,
    });
}
