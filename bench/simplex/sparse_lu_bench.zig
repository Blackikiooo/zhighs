//! Reusable sparse INVERT/FTRAN/BTRAN benchmark for perf and HiGHS parity work.

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

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;
    var args = std.process.Args.iterate(init.minimal.args);
    _ = args.next();
    const n = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 512;
    const repeats = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 21;
    if (n < 3 or n > std.math.maxInt(u32) or repeats == 0 or args.next() != null) return error.InvalidArguments;

    const starts = try allocator.alloc(zhighs.foundation.HUInt, n + 1);
    defer allocator.free(starts);
    const rows = try allocator.alloc(zhighs.foundation.RowId, n * 3);
    defer allocator.free(rows);
    const values = try allocator.alloc(f64, n * 3);
    defer allocator.free(values);
    var output: usize = 0;
    for (0..n) |column| {
        starts[column] = @intCast(output);
        const candidates = [_]usize{ if (column == 0) n - 1 else column - 1, column, if (column + 1 == n) 0 else column + 1 };
        var sorted = candidates;
        std.mem.sort(usize, &sorted, {}, std.sort.asc(usize));
        for (sorted) |row| {
            rows[output] = zhighs.foundation.RowId.fromUsizeAssumeValid(row);
            values[output] = if (row == column) 4.0 + @as(f64, @floatFromInt(column % 7)) * 0.125 else -0.5;
            output += 1;
        }
    }
    starts[n] = @intCast(output);
    const basis = zhighs.matrix.SparseBasisView{ .dimension = n, .starts = starts, .rows = rows, .values = values };
    var lu = zhighs.matrix.SparseLU.init(allocator);
    defer lu.deinit();
    try lu.factorize(basis);
    const samples = try allocator.alloc(u64, repeats);
    defer allocator.free(samples);
    for (samples) |*sample| {
        const started = nowNs();
        try lu.factorize(basis);
        sample.* = @intCast(nowNs() - started);
    }
    std.mem.sort(u64, samples, {}, std.sort.asc(u64));
    const factor_ns = samples[samples.len / 2];
    const rhs = try allocator.alloc(f64, n);
    defer allocator.free(rhs);
    const seed = try allocator.alloc(f64, n);
    defer allocator.free(seed);
    for (seed, 0..) |*value, index| value.* = @as(f64, @floatFromInt(index % 17)) * 0.125 + 1.0;
    const started = nowNs();
    var checksum: f64 = 0.0;
    for (0..repeats * 10) |_| {
        @memcpy(rhs, seed);
        try lu.solve(rhs);
        checksum += rhs[0];
        @memcpy(rhs, seed);
        try lu.solveTranspose(rhs);
        checksum += rhs[n - 1];
    }
    const solves_ns: u64 = @intCast(nowNs() - started);
    std.debug.print("implementation,dimension,basis_nnz,factor_nnz,fill,median_invert_ns,ftran_btran_ns,checksum\n", .{});
    std.debug.print("zhighs,{d},{d},{d},{d},{d},{d},{d:.17}\n", .{
        n, basis.nnz(), lu.factorNonzeros(), lu.inserted_fill, factor_ns, solves_ns / (repeats * 10), checksum,
    });
}
