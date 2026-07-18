//! Factor and validate a deterministic `ZHIGHS_BASIS_V1` CSC interchange file.

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
    const path = args.next() orelse return error.InvalidArguments;
    const repeats = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 21;
    if (repeats == 0 or args.next() != null) return error.InvalidArguments;
    const allocator = std.heap.c_allocator;
    const file = try std.Io.Dir.cwd().openFile(init.io, path, .{});
    defer file.close(init.io);
    var buffer: [64 * 1024]u8 = undefined;
    var reader = file.reader(init.io, &buffer);
    const content = try reader.interface.allocRemaining(allocator, .limited(2 * 1024 * 1024 * 1024));
    defer allocator.free(content);
    var fields = std.mem.tokenizeAny(u8, content, " \t\r\n");
    if (!std.mem.eql(u8, fields.next() orelse return error.InvalidBasisFile, "ZHIGHS_BASIS_V1")) return error.InvalidBasisFile;
    const n = try std.fmt.parseUnsigned(usize, fields.next() orelse return error.InvalidBasisFile, 10);
    const nnz = try std.fmt.parseUnsigned(usize, fields.next() orelse return error.InvalidBasisFile, 10);
    const starts = try allocator.alloc(zhighs.foundation.HUInt, n + 1);
    defer allocator.free(starts);
    const rows = try allocator.alloc(zhighs.foundation.RowId, nnz);
    defer allocator.free(rows);
    const values = try allocator.alloc(f64, nnz);
    defer allocator.free(values);
    var output: usize = 0;
    for (0..n) |column| {
        starts[column] = @intCast(output);
        const count = try std.fmt.parseUnsigned(usize, fields.next() orelse return error.InvalidBasisFile, 10);
        for (0..count) |_| {
            rows[output] = try zhighs.foundation.RowId.fromUsize(try std.fmt.parseUnsigned(usize, fields.next() orelse return error.InvalidBasisFile, 10));
            values[output] = try std.fmt.parseFloat(f64, fields.next() orelse return error.InvalidBasisFile);
            output += 1;
        }
    }
    if (output != nnz or fields.next() != null) return error.InvalidBasisFile;
    starts[n] = @intCast(nnz);
    const basis = zhighs.matrix.SparseBasisView{ .dimension = n, .starts = starts, .rows = rows, .values = values };
    var lu = zhighs.matrix.SparseLU.init(allocator);
    defer lu.deinit();
    try lu.factorize(basis);
    const trace_rows = try allocator.dupe(u32, lu.pivot_rows[0..n]);
    defer allocator.free(trace_rows);
    const trace_columns = try allocator.dupe(u32, lu.pivot_columns[0..n]);
    defer allocator.free(trace_columns);
    const samples = try allocator.alloc(u64, repeats);
    defer allocator.free(samples);
    for (samples) |*sample| {
        const began = nowNs();
        try lu.factorizeAssumeValid(basis);
        sample.* = @intCast(nowNs() - began);
    }
    std.mem.sort(u64, samples, {}, std.sort.asc(u64));
    const replay_samples = try allocator.alloc(u64, repeats);
    defer allocator.free(replay_samples);
    for (replay_samples) |*sample| {
        const began = nowNs();
        try lu.factorizeWithTraceAssumeValid(basis, .{ .rows = trace_rows, .columns = trace_columns });
        sample.* = @intCast(nowNs() - began);
    }
    std.mem.sort(u64, replay_samples, {}, std.sort.asc(u64));
    const rhs = try allocator.alloc(f64, n);
    defer allocator.free(rhs);
    for (rhs, 0..) |*value, index| value.* = 1.0 + @as(f64, @floatFromInt(index % 13)) * 0.125;
    const original = try allocator.dupe(f64, rhs);
    defer allocator.free(original);
    try lu.solve(rhs);
    var residual: f64 = 0;
    for (0..n) |row| {
        var product: f64 = 0;
        for (0..n) |column| for (@as(usize, starts[column])..@as(usize, starts[column + 1])) |entry|
            if (rows[entry].toUsize() == row) { product += values[entry] * rhs[column]; };
        residual = @max(residual, @abs(product - original[row]));
    }
    std.debug.print("zhighs-basis,{s},{d},{d},{d},{d},{d},{d},{e},{d}\n", .{
        path, n, nnz, lu.factorNonzeros(), lu.inserted_fill, samples[repeats / 2], replay_samples[repeats / 2], residual, lu.requestedBytes(),
    });
}
