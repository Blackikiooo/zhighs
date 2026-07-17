//! End-to-end zhighs model reader benchmark.
//!
//! Measures file open + read + parse + canonical CSC construction. Model
//! destruction and output are outside the timed region. Output is one TSV row
//! compatible with `highs_parser_bench.cpp`.
//!
//! Usage: io-parser-bench MODEL [iterations=7] [warmups=2]
//!                              [input=automatic|buffered|mmap]
//!                              [names=keep|drop]

const std = @import("std");
const model_io = @import("io");

fn nowNs() i128 {
    var timestamp: std.posix.timespec = undefined;
    switch (std.posix.errno(std.posix.system.clock_gettime(.MONOTONIC, &timestamp))) {
        .SUCCESS => {},
        else => unreachable,
    }
    return @as(i128, timestamp.sec) * std.time.ns_per_s + timestamp.nsec;
}

fn peakRssKb() usize {
    const value = std.posix.getrusage(std.posix.rusage.SELF).maxrss;
    return if (value > 0) @intCast(value) else 0;
}

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;
    const io = init.io;
    var args = std.process.Args.iterate(init.minimal.args);
    _ = args.next();
    const path = args.next() orelse return error.MissingModelPath;
    const iterations = if (args.next()) |value| try std.fmt.parseUnsigned(usize, value, 10) else 7;
    const warmups = if (args.next()) |value| try std.fmt.parseUnsigned(usize, value, 10) else 2;
    const input_name = args.next() orelse "automatic";
    const input_mode: model_io.InputMode = if (std.mem.eql(u8, input_name, "automatic"))
        .automatic
    else if (std.mem.eql(u8, input_name, "buffered"))
        .buffered
    else if (std.mem.eql(u8, input_name, "mmap"))
        .memory_map
    else
        return error.InvalidInputMode;
    const names_name = args.next() orelse "keep";
    const keep_names = if (std.mem.eql(u8, names_name, "keep"))
        true
    else if (std.mem.eql(u8, names_name, "drop"))
        false
    else
        return error.InvalidArguments;
    if (iterations == 0 or args.next() != null) return error.InvalidArguments;

    const stat = try std.Io.Dir.cwd().statFile(io, path, .{});
    var samples = try allocator.alloc(u64, iterations);
    defer allocator.free(samples);
    var rows: usize = 0;
    var columns: usize = 0;
    var nonzeros: usize = 0;
    var checksum: f64 = 0.0;
    for (0..warmups + iterations) |run| {
        const started = nowNs();
        var parsed = try model_io.readFile(io, allocator, path, .{ .input_mode = input_mode, .keep_names = keep_names });
        const elapsed: u64 = @intCast(nowNs() - started);
        rows = parsed.row_lower.len;
        columns = parsed.col_cost.len;
        nonzeros = parsed.matrix.nnz();
        checksum = parsed.objective_offset;
        for (parsed.col_cost) |value| checksum += value;
        for (parsed.matrix.values) |value| checksum += value;
        if (run >= warmups) samples[run - warmups] = elapsed;
        parsed.deinit();
    }
    std.mem.sort(u64, samples, {}, std.sort.asc(u64));
    const best = samples[0];
    const median = samples[samples.len / 2];
    const mib = @as(f64, @floatFromInt(stat.size)) / (1024.0 * 1024.0);
    const throughput = mib / (@as(f64, @floatFromInt(median)) / std.time.ns_per_s);
    const stdout = std.Io.File.stdout();
    const implementation = if (!keep_names)
        "zhighs-no-names"
    else switch (input_mode) {
        .automatic => "zhighs",
        .buffered => "zhighs-buffered",
        .memory_map => "zhighs-mmap",
    };
    const line = try std.fmt.allocPrint(allocator, "{s}\t{s}\t{d}\t{d}\t{d}\t{d}\t{d:.3}\t{d:.3}\t{d}\t{d:.17}\n", .{
        implementation, path, rows, columns, nonzeros, best, @as(f64, @floatFromInt(median)) / 1e6, throughput, peakRssKb(), checksum,
    });
    defer allocator.free(line);
    try stdout.writeStreamingAll(io, line);
}
