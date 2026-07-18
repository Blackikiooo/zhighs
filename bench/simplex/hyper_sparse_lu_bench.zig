//! Dense versus factor-graph-reachable sparse solve benchmark.

const std = @import("std");
const zhighs = @import("zhighs");

fn nowNs() i128 {
    var ts: std.posix.timespec = undefined;
    _ = std.posix.system.clock_gettime(.MONOTONIC, &ts);
    return @as(i128, ts.sec) * std.time.ns_per_s + ts.nsec;
}

pub fn main(init: std.process.Init) !void {
    var args = std.process.Args.iterate(init.minimal.args);
    _ = args.next();
    const n = if (args.next()) |value| try std.fmt.parseUnsigned(usize, value, 10) else 4096;
    const block = if (args.next()) |value| try std.fmt.parseUnsigned(usize, value, 10) else 32;
    const repeats = if (args.next()) |value| try std.fmt.parseUnsigned(usize, value, 10) else 10001;
    if (n == 0 or block == 0 or n % block != 0 or repeats == 0 or args.next() != null) return error.InvalidArguments;
    const allocator = std.heap.c_allocator;
    const starts = try allocator.alloc(zhighs.foundation.HUInt, n + 1);
    defer allocator.free(starts);
    const rows = try allocator.alloc(zhighs.foundation.RowId, n * 3);
    defer allocator.free(rows);
    const values = try allocator.alloc(f64, n * 3);
    defer allocator.free(values);
    var nnz: usize = 0;
    for (0..n) |column| {
        starts[column] = @intCast(nnz);
        const local = column % block;
        if (local > 0) {
            rows[nnz] = zhighs.foundation.RowId.fromUsizeAssumeValid(column - 1);
            values[nnz] = -0.25;
            nnz += 1;
        }
        rows[nnz] = zhighs.foundation.RowId.fromUsizeAssumeValid(column);
        values[nnz] = 2;
        nnz += 1;
        if (local + 1 < block) {
            rows[nnz] = zhighs.foundation.RowId.fromUsizeAssumeValid(column + 1);
            values[nnz] = -0.25;
            nnz += 1;
        }
    }
    starts[n] = @intCast(nnz);
    const basis = zhighs.matrix.SparseBasisView{ .dimension = n, .starts = starts, .rows = rows[0..nnz], .values = values[0..nnz] };
    var lu = zhighs.matrix.SparseLU.init(allocator);
    defer lu.deinit();
    try lu.factorizeAssumeValid(basis);
    const rhs = try allocator.alloc(f64, n);
    defer allocator.free(rhs);
    const output = try allocator.alloc(u32, n);
    defer allocator.free(output);
    const input = [_]u32{0};
    @memset(rhs, 0);
    rhs[0] = 1;
    _ = try lu.solveHyperSparse(rhs, &input, output); // build companion views outside timing

    const Modes = enum { dense_ftran, hyper_ftran, adaptive_ftran, dense_btran, hyper_btran, adaptive_btran };
    inline for (std.meta.tags(Modes)) |mode| {
        var checksum: f64 = 0;
        const began = nowNs();
        for (0..repeats) |_| {
            if (mode != .hyper_ftran and mode != .hyper_btran) @memset(rhs, 0);
            rhs[0] = 1;
            switch (mode) {
                .dense_ftran => try lu.solve(rhs),
                .hyper_ftran => _ = try lu.solveHyperSparse(rhs, &input, output),
                .adaptive_ftran => try lu.solveAdaptive(rhs, &input, false),
                .dense_btran => try lu.solveTranspose(rhs),
                .hyper_btran => _ = try lu.solveTransposeHyperSparse(rhs, &input, output),
                .adaptive_btran => try lu.solveAdaptive(rhs, &input, true),
            }
            checksum += rhs[0];
        }
        const elapsed: u64 = @intCast(nowNs() - began);
        std.debug.print("{s},{d},{d},{d},{d},{d:.17}\n", .{ @tagName(mode), n, block, repeats, elapsed / repeats, checksum });
    }
}
