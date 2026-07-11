//! Matrix benchmark compatible with the HiGHS C++ comparison harness.
//! Always run with: zig build bench-matrix -Doptimize=ReleaseFast

const std = @import("std");
const zhighs = @import("zhighs");

const dimension: usize = 50_000;
const nnz: usize = dimension * 3 - 2;
const product_repeats: usize = 200;
const quad_repeats: usize = 100;
const transform_repeats: usize = 100;
const accumulator_repeats: usize = 200;

fn nowNs() i128 {
    var ts: std.posix.timespec = undefined;
    switch (std.posix.errno(std.posix.system.clock_gettime(.MONOTONIC, &ts))) {
        .SUCCESS => {},
        else => unreachable,
    }
    return @as(i128, ts.sec) * std.time.ns_per_s + ts.nsec;
}

fn report(name: []const u8, repeats: usize, start: i128, result_checksum: f64) void {
    const total: u64 = @intCast(nowNs() - start);
    const per_repeat = @as(f64, @floatFromInt(total)) / @as(f64, @floatFromInt(repeats));
    std.debug.print("zig,{s},{d},{d},{d},{d},{d:.3},{d:.17}\n", .{ name, dimension, nnz, repeats, total, per_repeat, result_checksum });
}

fn checksum(values: []const f64) f64 {
    var sum: f64 = 0.0;
    for (values) |value| sum += value;
    return sum;
}

inline fn clobber(values: []f64) void {
    clobberPtr(values.ptr);
}

inline fn clobberPtr(pointer: anytype) void {
    asm volatile (""
        :
        : [pointer] "r" (pointer),
        : .{ .memory = true });
}

fn fillSorted(builder: *zhighs.matrix.MatrixBuilder, allocator: std.mem.Allocator) !void {
    for (0..dimension) |col| {
        const col_id = try zhighs.ColId.fromUsize(col);
        if (col != 0) try builder.append(allocator, try zhighs.RowId.fromUsize(col - 1), col_id, -1.0);
        try builder.append(allocator, try zhighs.RowId.fromUsize(col), col_id, 4.0);
        if (col + 1 < dimension) try builder.append(allocator, try zhighs.RowId.fromUsize(col + 1), col_id, -1.0);
    }
}

fn fillUnsorted(builder: *zhighs.matrix.MatrixBuilder, allocator: std.mem.Allocator) !void {
    var col = dimension;
    while (col != 0) {
        col -= 1;
        const col_id = try zhighs.ColId.fromUsize(col);
        if (col + 1 < dimension) try builder.append(allocator, try zhighs.RowId.fromUsize(col + 1), col_id, -1.0);
        try builder.append(allocator, try zhighs.RowId.fromUsize(col), col_id, 4.0);
        if (col != 0) try builder.append(allocator, try zhighs.RowId.fromUsize(col - 1), col_id, -1.0);
    }
}

pub fn main() !void {
    // Fair counterpart to C++ malloc/vector allocation; page_allocator would
    // add mmap/munmap cost to every transformation benchmark.
    const allocator = std.heap.smp_allocator;
    var builder = try zhighs.matrix.MatrixBuilder.init(dimension, dimension);
    defer builder.deinit(allocator);
    try builder.reserve(allocator, nnz);
    try fillSorted(&builder, allocator);
    var matrix = try builder.freezeSortedAssumeValid(allocator, 0.0);
    defer matrix.deinit(allocator);

    const cursor_storage = try allocator.alloc(zhighs.HUInt, dimension + 48);
    defer allocator.free(cursor_storage);
    const cursor = cursor_storage[48..][0..dimension];
    var csr_cache = try zhighs.matrix.CsrCache.buildWithScratchAssumeValid(allocator, matrix, 0, cursor);
    defer csr_cache.deinit(allocator);
    const csr = csr_cache.viewAssumeCurrent();
    const csr_starts_storage = try allocator.alloc(zhighs.HUInt, dimension + 1 + 16);
    defer allocator.free(csr_starts_storage);
    const reusable_csr_starts = csr_starts_storage[16..][0 .. dimension + 1];
    const transpose_starts_storage = try allocator.alloc(usize, dimension + 1 + 8);
    defer allocator.free(transpose_starts_storage);
    const reusable_transpose_starts = transpose_starts_storage[8..][0 .. dimension + 1];
    const cols_storage = try allocator.alloc(zhighs.ColId, nnz + 32);
    defer allocator.free(cols_storage);
    const reusable_cols = cols_storage[32..][0..nnz];
    const rows_storage = try allocator.alloc(zhighs.RowId, nnz + 16);
    defer allocator.free(rows_storage);
    const reusable_rows = rows_storage[16..][0..nnz];
    const values_storage = try allocator.alloc(f64, nnz + 24);
    defer allocator.free(values_storage);
    const reusable_values = values_storage[24..][0..nnz];

    const dense_x_storage = try allocator.alloc(f64, dimension + 8);
    defer allocator.free(dense_x_storage);
    const dense_x = dense_x_storage[8..][0..dimension];
    const sparse_x = try allocator.alloc(f64, dimension);
    defer allocator.free(sparse_x);
    const sparse_ids = try allocator.alloc(zhighs.ColId, dimension / 20);
    defer allocator.free(sparse_ids);
    const sparse_values = try allocator.alloc(f64, dimension / 20);
    defer allocator.free(sparse_values);
    const y_storage = try allocator.alloc(f64, dimension + 40);
    defer allocator.free(y_storage);
    const y = y_storage[40..][0..dimension];
    const hcd_storage = try allocator.alloc(zhighs.HCD, dimension + 12);
    defer allocator.free(hcd_storage);
    const hcd_scratch = hcd_storage[12..][0..dimension];
    const row_scale = try allocator.alloc(f64, dimension);
    defer allocator.free(row_scale);
    const col_scale = try allocator.alloc(f64, dimension);
    defer allocator.free(col_scale);
    @memset(dense_x, 1.0);
    for (sparse_x, 0..) |*value, index| value.* = if (index % 20 == 0) 1.0 else 0.0;
    for (sparse_ids, sparse_values, 0..) |*id, *value, index| {
        id.* = try zhighs.ColId.fromUsize(index * 20);
        value.* = 1.0;
    }
    const sparse_view: zhighs.matrix.SparseVectorView(zhighs.ColId) = .{ .dimension = dimension, .indices = sparse_ids, .values = sparse_values };
    @memset(row_scale, 1.0);
    @memset(col_scale, 1.0);

    std.debug.print("implementation,kernel,dimension,nnz,repeats,total_ns,ns_per_repeat,checksum\n", .{});

    var result_checksum: f64 = 0.0;
    var start = nowNs();
    for (0..product_repeats) |_| clobber(y);
    result_checksum = checksum(y);
    report("barrier_only", product_repeats, start, result_checksum);

    start = nowNs();
    for (0..product_repeats) |_| {
        @memset(std.mem.sliceAsBytes(y), 0);
        clobber(y);
    }
    result_checksum = checksum(y);
    report("clear_output_bytes", product_repeats, start, result_checksum);

    start = nowNs();
    for (0..product_repeats) |_| {
        zhighs.matrix.clearF64(y);
        clobber(y);
    }
    result_checksum = checksum(y);
    report("clear_output", product_repeats, start, result_checksum);

    start = nowNs();
    for (0..product_repeats) |_| {
        matrix.multiplyAssumeValid(dense_x, y);
        clobber(y);
    }
    result_checksum = checksum(y);
    report("csc_ax_dense", product_repeats, start, result_checksum);

    start = nowNs();
    for (0..product_repeats) |_| {
        matrix.multiplySkippingZerosAssumeValid(sparse_x, y);
        clobber(y);
    }
    result_checksum = checksum(y);
    report("csc_ax_sparse_skip", product_repeats, start, result_checksum);

    start = nowNs();
    for (0..product_repeats) |_| {
        matrix.multiplySparseAssumeValid(sparse_view, y);
        clobber(y);
    }
    result_checksum = checksum(y);
    report("csc_ax_sparse_view", product_repeats, start, result_checksum);

    @memset(y, 0.0);
    start = nowNs();
    for (0..product_repeats) |_| {
        matrix.addSparseProductAssumeValid(sparse_view, y);
        clobber(y);
    }
    result_checksum = checksum(y);
    report("csc_sparse_add_no_clear", product_repeats, start, result_checksum);

    start = nowNs();
    for (0..product_repeats) |_| {
        csr.multiplyAssumeValid(dense_x, y);
        clobber(y);
    }
    result_checksum = checksum(y);
    report("csr_ax_dense", product_repeats, start, result_checksum);

    start = nowNs();
    for (0..product_repeats) |_| {
        matrix.transposeMultiplyAssumeValid(dense_x, y);
        clobber(y);
    }
    result_checksum = checksum(y);
    report("csc_atx_dense", product_repeats, start, result_checksum);

    start = nowNs();
    for (0..product_repeats) |_| {
        csr.transposeMultiplyAssumeValid(dense_x, y);
        clobber(y);
    }
    result_checksum = checksum(y);
    report("csr_atx_dense", product_repeats, start, result_checksum);

    @memset(y, 0.0);
    start = nowNs();
    for (0..product_repeats) |_| {
        zhighs.matrix.addProductAssumeValid(matrix, 1.0, dense_x, y);
        clobber(y);
    }
    result_checksum = checksum(y);
    report("alpha_ax_plus_y", product_repeats, start, result_checksum);

    start = nowNs();
    for (0..quad_repeats) |_| {
        zhighs.matrix.multiplyHighPrecisionAssumeValid(matrix, dense_x, y, hcd_scratch);
        clobber(y);
    }
    result_checksum = checksum(y);
    report("product_exact_robust", quad_repeats, start, result_checksum);

    start = nowNs();
    for (0..quad_repeats) |_| {
        zhighs.matrix.multiplyHighPrecisionFastAssumeValid(matrix, dense_x, y, hcd_scratch);
        clobber(y);
    }
    result_checksum = checksum(y);
    report("product_exact_fast", quad_repeats, start, result_checksum);

    start = nowNs();
    for (0..quad_repeats) |_| {
        zhighs.matrix.multiplyCompensatedAssumeValid(matrix, dense_x, y, hcd_scratch);
        clobber(y);
    }
    result_checksum = checksum(y);
    report("product_quad", quad_repeats, start, result_checksum);

    start = nowNs();
    for (0..product_repeats) |_| zhighs.matrix.applyScalingAssumeValid(&matrix, .{ .row = row_scale, .col = col_scale });
    result_checksum = checksum(matrix.values);
    report("apply_scale", product_repeats, start, result_checksum);

    start = nowNs();
    for (0..transform_repeats) |revision| {
        var cache = try zhighs.matrix.CsrCache.buildWithScratchAssumeValid(allocator, matrix, revision, cursor);
        clobberPtr(cache.values.ptr);
        cache.deinit(allocator);
    }
    report("csc_to_csr_scratch", transform_repeats, start, matrix.values[0]);

    start = nowNs();
    for (0..transform_repeats) |_| {
        try zhighs.matrix.fillCsrFromCscAssumeValid(matrix, reusable_csr_starts, reusable_cols, reusable_values, cursor);
        clobberPtr(reusable_values.ptr);
    }
    report("csc_to_csr_into", transform_repeats, start, reusable_values[0]);

    start = nowNs();
    for (0..transform_repeats) |_| {
        var transposed = try zhighs.matrix.transposeAssumeValid(allocator, matrix);
        clobberPtr(transposed.values.ptr);
        transposed.deinit(allocator);
    }
    report("transpose", transform_repeats, start, matrix.values[0]);

    start = nowNs();
    for (0..transform_repeats) |_| {
        try zhighs.matrix.transposeIntoAssumeValid(matrix, reusable_transpose_starts, reusable_rows, reusable_values, cursor);
        clobberPtr(reusable_values.ptr);
    }
    report("transpose_into", transform_repeats, start, reusable_values[0]);

    builder.clearRetainingCapacity();
    try fillSorted(&builder, allocator);
    start = nowNs();
    var sorted_matrix = try builder.freezeSortedAssumeValid(allocator, 0.0);
    const sorted_elapsed_checksum = sorted_matrix.values[sorted_matrix.values.len / 2];
    std.mem.doNotOptimizeAway(sorted_elapsed_checksum);
    report("builder_freeze_sorted", 1, start, sorted_elapsed_checksum);
    sorted_matrix.deinit(allocator);

    builder.clearRetainingCapacity();
    try fillUnsorted(&builder, allocator);
    start = nowNs();
    var general_matrix = try builder.freeze(allocator, 0.0);
    const general_checksum = general_matrix.values[general_matrix.values.len / 2];
    std.mem.doNotOptimizeAway(general_checksum);
    report("builder_freeze_general", 1, start, general_checksum);
    general_matrix.deinit(allocator);

    var accumulator = try zhighs.matrix.SparseAccumulator(zhighs.RowId).initWithCapacity(allocator, dimension, dimension);
    defer accumulator.deinit(allocator);
    start = nowNs();
    for (0..accumulator_repeats) |_| {
        accumulator.clear();
        for (0..dimension) |index| {
            const id = zhighs.RowId.fromUsizeAssumeValid(index);
            accumulator.addAssumeValid(id, 1.0);
            accumulator.addAssumeValid(id, -0.5);
        }
        clobberPtr(&accumulator);
    }
    result_checksum = accumulator.get(zhighs.RowId.fromUsizeAssumeValid(dimension / 2));
    report("sparse_accumulate", accumulator_repeats, start, result_checksum);
}
