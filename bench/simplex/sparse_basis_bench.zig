//! Sparse basis-assembly microbenchmark.
//!
//! Build with `zig build build-bench-sparse-basis -Doptimize=ReleaseFast
//! -Dcpu=native`.  The installed binary accepts optional `dimension`,
//! `entries_per_column`, and `repetitions` arguments, which makes the same
//! leaf kernel usable under `perf stat` and disassembly inspection.

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
    const dimension = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 100_000;
    const entries_per_column = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 8;
    const repetitions = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 101;
    if (dimension == 0 or entries_per_column == 0 or entries_per_column > dimension or repetitions == 0 or args.next() != null)
        return error.InvalidArguments;

    const nnz = try std.math.mul(usize, dimension, entries_per_column);
    const starts = try allocator.alloc(usize, dimension + 1);
    defer allocator.free(starts);
    const rows = try allocator.alloc(zhighs.foundation.RowId, nnz);
    defer allocator.free(rows);
    const values = try allocator.alloc(f64, nnz);
    defer allocator.free(values);
    const basic_index = try allocator.alloc(u32, dimension);
    defer allocator.free(basic_index);
    const row_scale = try allocator.alloc(f64, dimension);
    defer allocator.free(row_scale);
    const artificial_sign = try allocator.alloc(f64, dimension);
    defer allocator.free(artificial_sign);

    for (0..dimension) |column| {
        starts[column] = column * entries_per_column;
        basic_index[column] = @intCast(column);
        row_scale[column] = if (column & 1 == 0) 1.0 else -1.0;
        artificial_sign[column] = 0.0;
        // The stride is coprime to common powers of two; sorting the small
        // local set preserves canonical CSC without hiding random-row traffic.
        const column_rows = rows[starts[column]..][0..entries_per_column];
        for (column_rows, 0..) |*row, entry| row.* = zhighs.foundation.RowId.fromUsizeAssumeValid((column * 17 + entry * 7919) % dimension);
        std.mem.sort(zhighs.foundation.RowId, column_rows, {}, struct {
            fn lessThan(_: void, lhs: zhighs.foundation.RowId, rhs: zhighs.foundation.RowId) bool {
                return lhs.toUsize() < rhs.toUsize();
            }
        }.lessThan);
        for (values[starts[column]..][0..entries_per_column], 0..) |*value, entry|
            value.* = 0.25 + @as(f64, @floatFromInt(entry + 1)) * 0.125;
    }
    starts[dimension] = nnz;
    const matrix = zhighs.matrix.CscView.initAssumeValid(dimension, dimension, starts, rows, values);
    var buffers = zhighs.matrix.SparseBasisBuffers.init(allocator);
    defer buffers.deinit();
    _ = try buffers.assemble(matrix, basic_index, row_scale, artificial_sign);

    const samples = try allocator.alloc(u64, repetitions);
    defer allocator.free(samples);
    var checksum: f64 = 0.0;
    for (samples) |*sample| {
        const started = nowNs();
        const basis = try buffers.assemble(matrix, basic_index, row_scale, artificial_sign);
        sample.* = @intCast(nowNs() - started);
        checksum += basis.values[(basis.nnz() / 2) % basis.nnz()];
        std.mem.doNotOptimizeAway(basis.values.ptr);
    }
    std.mem.sort(u64, samples, {}, std.sort.asc(u64));
    const median = samples[samples.len / 2];
    const bytes = nnz * (@sizeOf(zhighs.foundation.RowId) + @sizeOf(f64)) + (dimension + 1) * @sizeOf(zhighs.foundation.HUInt);
    std.debug.print("implementation,dimension,nnz,repetitions,median_ns,GiB_per_s,retained_bytes,checksum\n", .{});
    std.debug.print("zhighs,{d},{d},{d},{d},{d:.3},{d},{d:.17}\n", .{
        dimension,
        nnz,
        repetitions,
        median,
        (@as(f64, @floatFromInt(bytes)) / @as(f64, @floatFromInt(median))) *
            (@as(f64, std.time.ns_per_s) / (1024.0 * 1024.0 * 1024.0)),
        buffers.starts_capacity * @sizeOf(zhighs.foundation.HUInt) + buffers.entry_capacity * (@sizeOf(zhighs.foundation.RowId) + @sizeOf(f64)),
        checksum,
    });
}
