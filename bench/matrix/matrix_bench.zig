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
const appended_rows: usize = 1_024;
const appended_nnz: usize = appended_rows * 3;

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
    var store_matrix = try matrix.clone(allocator);
    var matrix_store = zhighs.matrix.MatrixStore.initAssumeValid(store_matrix);
    store_matrix = undefined;
    defer matrix_store.deinit(allocator);

    var csr_buffers = try zhighs.matrix.CsrBuffers.init(allocator, dimension, nnz);
    defer csr_buffers.deinit(allocator);
    var transpose_buffers = try zhighs.matrix.TransposeBuffers.init(allocator, dimension, nnz);
    defer transpose_buffers.deinit(allocator);
    var transform_buffers = try zhighs.matrix.CscTransformBuffers.initCapacity(allocator, dimension, nnz + appended_nnz, dimension);
    defer transform_buffers.deinit(allocator);
    var csr_cache = try zhighs.matrix.CsrCache.buildWithScratchAssumeValid(allocator, matrix, 0, csr_buffers.cursor);
    defer csr_cache.deinit(allocator);
    const csr = csr_cache.viewAssumeCurrent();

    const dense_x_storage = try allocator.alloc(f64, dimension + 8);
    defer allocator.free(dense_x_storage);
    const dense_x = dense_x_storage[8..][0..dimension];
    const sparse_x = try allocator.alloc(f64, dimension);
    defer allocator.free(sparse_x);
    const sparse_ids = try allocator.alloc(zhighs.ColId, dimension / 20);
    defer allocator.free(sparse_ids);
    const sparse_values = try allocator.alloc(f64, dimension / 20);
    defer allocator.free(sparse_values);
    const y = try allocator.alloc(f64, dimension);
    defer allocator.free(y);
    const hcd_storage = try allocator.alloc(zhighs.HCD, dimension + 12);
    defer allocator.free(hcd_storage);
    const hcd_scratch = hcd_storage[12..][0..dimension];
    const row_scale = try allocator.alloc(f64, dimension);
    defer allocator.free(row_scale);
    const col_scale = try allocator.alloc(f64, dimension);
    defer allocator.free(col_scale);
    const row_permutation = try allocator.alloc(zhighs.RowId, dimension);
    defer allocator.free(row_permutation);
    const col_permutation = try allocator.alloc(zhighs.ColId, dimension);
    defer allocator.free(col_permutation);
    const appended_row_starts = try allocator.alloc(usize, appended_rows + 1);
    defer allocator.free(appended_row_starts);
    const appended_col_indices = try allocator.alloc(zhighs.ColId, appended_nnz);
    defer allocator.free(appended_col_indices);
    const appended_values = try allocator.alloc(f64, appended_nnz);
    defer allocator.free(appended_values);
    @memset(dense_x, 1.0);
    for (sparse_x, 0..) |*value, index| value.* = if (index % 20 == 0) 1.0 else 0.0;
    for (sparse_ids, sparse_values, 0..) |*id, *value, index| {
        id.* = try zhighs.ColId.fromUsize(index * 20);
        value.* = 1.0;
    }
    const sparse_view: zhighs.matrix.SparseVectorView(zhighs.ColId) = .{ .dimension = dimension, .indices = sparse_ids, .values = sparse_values };
    @memset(row_scale, 1.0);
    @memset(col_scale, 1.0);
    for (row_permutation, col_permutation, 0..) |*row, *col, index| {
        row.* = try zhighs.RowId.fromUsize(dimension - 1 - index);
        col.* = try zhighs.ColId.fromUsize(dimension - 1 - index);
    }
    for (0..appended_rows) |row| {
        appended_row_starts[row] = row * 3;
        for (0..3) |offset| {
            appended_col_indices[row * 3 + offset] = try zhighs.ColId.fromUsize(row * 3 + offset);
            appended_values[row * 3 + offset] = @floatFromInt(offset + 1);
        }
    }
    appended_row_starts[appended_rows] = appended_nnz;

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
        var cache = try zhighs.matrix.CsrCache.buildWithScratchAssumeValid(allocator, matrix, revision, csr_buffers.cursor);
        clobberPtr(cache.values.ptr);
        cache.deinit(allocator);
    }
    report("csc_to_csr_scratch", transform_repeats, start, matrix.values[0]);

    start = nowNs();
    for (0..transform_repeats) |_| {
        try zhighs.matrix.fillCsrFromCscAssumeValid(matrix, csr_buffers.row_starts, csr_buffers.col_indices, csr_buffers.values, csr_buffers.cursor);
        clobberPtr(csr_buffers.values.ptr);
    }
    report("csc_to_csr_into", transform_repeats, start, csr_buffers.values[0]);

    _ = try matrix_store.csr(allocator);
    start = nowNs();
    for (0..transform_repeats) |repeat| {
        const replacement_value: f64 = if (repeat & 1 == 0) 5.0 else 4.0;
        try matrix_store.updateValuesAtPositionsAssumeValid(allocator, &.{0}, &.{replacement_value});
        const rebuilt = try matrix_store.csr(allocator);
        clobberPtr(rebuilt.values.ptr);
    }
    report("matrix_store_csr_rebuild", transform_repeats, start, matrix_store.csc().values[0]);

    start = nowNs();
    for (0..transform_repeats) |_| {
        var transposed = try zhighs.matrix.transposeAssumeValid(allocator, matrix);
        clobberPtr(transposed.values.ptr);
        transposed.deinit(allocator);
    }
    report("transpose", transform_repeats, start, matrix.values[0]);

    start = nowNs();
    for (0..transform_repeats) |_| {
        try zhighs.matrix.transposeIntoAssumeValid(matrix, transpose_buffers.starts, transpose_buffers.rows, transpose_buffers.values, transpose_buffers.cursor);
        clobberPtr(transpose_buffers.values.ptr);
    }
    report("transpose_into", transform_repeats, start, transpose_buffers.values[0]);

    start = nowNs();
    for (0..transform_repeats) |_| {
        var extracted = try zhighs.matrix.extractColumnRange(allocator, matrix, dimension / 4, 3 * dimension / 4);
        clobberPtr(extracted.values.ptr);
        extracted.deinit(allocator);
    }
    report("extract_column_range", transform_repeats, start, matrix.values[0]);

    start = nowNs();
    for (0..transform_repeats) |_| {
        const view = try zhighs.matrix.extractColumnRangeInto(&transform_buffers, matrix, dimension / 4, 3 * dimension / 4);
        clobberPtr(view.values.ptr);
    }
    report("extract_column_range_into", transform_repeats, start, transform_buffers.values[0]);

    start = nowNs();
    for (0..transform_repeats) |_| {
        var permuted = try zhighs.matrix.permute(allocator, matrix, row_permutation, col_permutation);
        clobberPtr(permuted.values.ptr);
        permuted.deinit(allocator);
    }
    report("permute", transform_repeats, start, matrix.values[0]);

    start = nowNs();
    for (0..transform_repeats) |_| {
        const view = try zhighs.matrix.permuteInto(&transform_buffers, matrix, row_permutation, col_permutation);
        clobberPtr(view.values.ptr);
    }
    report("permute_into", transform_repeats, start, transform_buffers.values[0]);

    start = nowNs();
    for (0..transform_repeats) |_| {
        var appended = try zhighs.matrix.appendRowsFromCsr(allocator, matrix, appended_row_starts, appended_col_indices, appended_values);
        clobberPtr(appended.values.ptr);
        appended.deinit(allocator);
    }
    report("append_rows_csr", transform_repeats, start, matrix.values[0]);

    start = nowNs();
    for (0..transform_repeats) |_| {
        const view = try zhighs.matrix.appendRowsFromCsrInto(&transform_buffers, matrix, appended_row_starts, appended_col_indices, appended_values);
        clobberPtr(view.values.ptr);
    }
    report("append_rows_csr_into", transform_repeats, start, transform_buffers.values[0]);

    builder.clearRetainingCapacity();
    try fillSorted(&builder, allocator);
    start = nowNs();
    // Match the C++ reference's single-offset representation. Applications
    // that prioritize repeated matrix products should use the default compact
    // freeze path and accept its additional construction work.
    var sorted_matrix = try builder.freezeSortedLeanAssumeValid(allocator, 0.0);
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
