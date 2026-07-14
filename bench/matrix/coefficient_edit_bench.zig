//! Batched model coefficient-edit throughput benchmark.
//! Run with: zig build bench-coefficient-edits -Doptimize=ReleaseFast

const std = @import("std");
const zhighs = @import("zhighs");

const dimension: usize = 4_096;
const value_batch_size: usize = 4_096;
const structural_batch_size: usize = 2_048;
const value_repeats: usize = 100;
const structural_repeats: usize = 40;
const scalar_batch_size: usize = 4_096;
const scalar_repeats: usize = 100;
const small_scalar_repeats: usize = 10_000;

fn nowNs() i128 {
    var ts: std.posix.timespec = undefined;
    switch (std.posix.errno(std.posix.system.clock_gettime(.MONOTONIC, &ts))) {
        .SUCCESS => {},
        else => unreachable,
    }
    return @as(i128, ts.sec) * std.time.ns_per_s + ts.nsec;
}

fn report(name: []const u8, batch_size: usize, repeats: usize, start: i128, checksum: f64) void {
    const total_ns: u64 = @intCast(nowNs() - start);
    const edit_count = batch_size * repeats;
    const edits_per_second = @as(f64, @floatFromInt(edit_count)) *
        @as(f64, @floatFromInt(std.time.ns_per_s)) /
        @as(f64, @floatFromInt(total_ns));
    std.debug.print("{s},{d},{d},{d},{d},{d:.3},{d:.17}\n", .{
        name,
        dimension,
        batch_size,
        repeats,
        total_ns,
        edits_per_second,
        checksum,
    });
}

fn buildTridiagonalModel(allocator: std.mem.Allocator, env: *zhighs.model.Env) !zhighs.model.Model {
    var model = try zhighs.model.Model.init(allocator, env, "coefficient_edit_bench");
    errdefer model.deinit();

    for (0..dimension) |_|
        try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, null);
    for (0..dimension) |col| {
        var indices: [3]usize = undefined;
        var values: [3]f64 = undefined;
        var count: usize = 0;
        if (col > 0) {
            indices[count] = col - 1;
            values[count] = -1.0;
            count += 1;
        }
        indices[count] = col;
        values[count] = 4.0;
        count += 1;
        if (col + 1 < dimension) {
            indices[count] = col + 1;
            values[count] = -1.0;
            count += 1;
        }
        try model.addVar(count, indices[0..count], values[0..count], 0.0, 0.0, 1.0, .continuous, null);
    }
    try model.updateModel();
    return model;
}

pub fn main() !void {
    const allocator = std.heap.smp_allocator;
    var env = try zhighs.model.Env.initSimple(allocator);
    defer env.deinit();
    var model = try buildTridiagonalModel(allocator, &env);
    defer model.deinit();

    std.debug.print("kernel,dimension,batch_size,repeats,total_ns,edits_per_second,checksum\n", .{});

    // Warm the retained pending queue and coefficient scratch allocations.
    for (0..value_batch_size) |index|
        try model.chgCoeff(index, index, 5.0);
    try model.updateModel();

    var start = nowNs();
    for (0..value_repeats) |repeat| {
        const value: f64 = if (repeat & 1 == 0) 6.0 else 5.0;
        for (0..value_batch_size) |index|
            try model.chgCoeff(index, index, value);
        try model.updateModel();
    }
    report("existing_values", value_batch_size, value_repeats, start, model.matrix.csc().values[0]);

    // Each batch changes half absent off-band positions and half existing
    // diagonal positions. Alternating batches restores the previous shape,
    // so every iteration performs one structural canonical merge.
    start = nowNs();
    for (0..structural_repeats) |repeat| {
        const insert = repeat & 1 == 0;
        for (0..structural_batch_size / 2) |index| {
            const col = index * 2;
            const far_row = (col + dimension / 2) % dimension;
            try model.chgCoeff(far_row, col, if (insert) @as(f64, 2.0) else 0.0);
            try model.chgCoeff(col, col, if (insert) @as(f64, 0.0) else 5.0);
        }
        try model.updateModel();
    }
    report("structural_mixed", structural_batch_size, structural_repeats, start, model.matrix.csc().values[0]);

    // Four writes per target exercise scalar-stream last-write-wins. The
    // final value alternates so every batch remains observable.
    start = nowNs();
    for (0..scalar_repeats) |repeat| {
        const final_value: f64 = if (repeat & 1 == 0) 3.0 else 4.0;
        for (0..scalar_batch_size / 4) |index| {
            try model.chgObj(index, 1.0);
            try model.chgObj(index, 2.0);
            try model.chgObj(index, 5.0);
            try model.chgObj(index, final_value);
        }
        try model.updateModel();
    }
    report("scalar_last_write_wins", scalar_batch_size, scalar_repeats, start, model.var_obj[0]);

    const small_batch_sizes = [_]usize{ 1, 4, 8, 16 };
    const small_batch_names = [_][]const u8{ "scalar_batch_1", "scalar_batch_4", "scalar_batch_8", "scalar_batch_16" };
    for (small_batch_sizes, small_batch_names) |batch_size, name| {
        start = nowNs();
        for (0..small_scalar_repeats) |repeat| {
            const value: f64 = if (repeat & 1 == 0) 6.0 else 7.0;
            for (0..batch_size) |index| try model.chgObj(index, value);
            try model.updateModel();
        }
        report(name, batch_size, small_scalar_repeats, start, model.var_obj[0]);
    }

    // Measure only the structural flush: queue construction is outside the
    // timer so this tracks plan flattening plus the single packed CSC rebuild.
    var append_model = try zhighs.model.Model.init(allocator, &env, "structural_append_bench");
    defer append_model.deinit();
    for (0..dimension) |_|
        try append_model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, null);
    for (0..dimension) |col| {
        const row = [_]usize{col};
        const value = [_]f64{4.0};
        try append_model.addVar(1, &row, &value, 0.0, 0.0, 1.0, .continuous, null);
        try append_model.chgObj(col, 2.0);
    }
    start = nowNs();
    try append_model.updateModel();
    report("append_rows_columns_folded_scalars", dimension * 3, 1, start, append_model.var_obj[0]);
}
