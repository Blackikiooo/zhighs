//! Dedicated profiling harness for perf analysis.
//! Build: zig build-exe bench/matrix/perf_profile.zig -I src/ -I lib/ -Doptimize=ReleaseFast -Dcpu=native --name perf-profile
//! Usage: perf record ./perf-profile <kernel-name>

const std = @import("std");
const zhighs = @import("zhighs");

const dimension: usize = 50_000;
const nnz: usize = dimension * 3 - 2;

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

fn benchCscToCsrInto(allocator: std.mem.Allocator, matrix: *const zhighs.matrix.CscMatrix, cursor: []zhighs.HUInt) void {
    const starts_storage = allocator.alloc(zhighs.HUInt, dimension + 1 + 16) catch unreachable;
    defer allocator.free(starts_storage);
    const reusable_starts = starts_storage[16..][0 .. dimension + 1];
    const cols_storage = allocator.alloc(zhighs.ColId, nnz + 32) catch unreachable;
    defer allocator.free(cols_storage);
    const reusable_cols = cols_storage[32..][0..nnz];
    const values_storage = allocator.alloc(f64, nnz + 24) catch unreachable;
    defer allocator.free(values_storage);
    const reusable_values = values_storage[24..][0..nnz];
    for (0..500) |_| {
        zhighs.matrix.fillCsrFromCscAssumeValid(matrix.*, reusable_starts, reusable_cols, reusable_values, cursor) catch unreachable;
        clobberPtr(reusable_values.ptr);
    }
}

fn benchSparseAccumulate(allocator: std.mem.Allocator) void {
    var accumulator = zhighs.matrix.SparseAccumulator(zhighs.RowId).initWithCapacity(allocator, dimension, dimension) catch unreachable;
    defer accumulator.deinit(allocator);
    for (0..500) |_| {
        accumulator.clear();
        for (0..dimension) |index| {
            const id = zhighs.RowId.fromUsizeAssumeValid(index);
            accumulator.addAssumeValid(id, 1.0);
            accumulator.addAssumeValid(id, -0.5);
        }
        clobberPtr(&accumulator);
    }
    _ = accumulator.get(zhighs.RowId.fromUsizeAssumeValid(dimension / 2));
}

fn benchCsrTransposeMultiply(allocator: std.mem.Allocator, csr: *const zhighs.matrix.CsrView, x: []const f64) void {
    const y_storage = allocator.alloc(f64, dimension + 40) catch unreachable;
    defer allocator.free(y_storage);
    const y = y_storage[40..][0..dimension];
    for (0..500) |_| {
        csr.transposeMultiplyAssumeValid(x, y);
        clobberPtr(y.ptr);
    }
}

fn benchCscTransposeMultiply(allocator: std.mem.Allocator, matrix: *const zhighs.matrix.CscMatrix, x: []const f64) void {
    const y = allocator.alloc(f64, dimension) catch unreachable;
    defer allocator.free(y);
    for (0..500) |_| {
        matrix.transposeMultiplyAssumeValid(x, y);
        clobberPtr(y.ptr);
    }
}

fn benchCscAxDense(allocator: std.mem.Allocator, matrix: *const zhighs.matrix.CscMatrix, x: []const f64) void {
    const y = allocator.alloc(f64, dimension) catch unreachable;
    defer allocator.free(y);
    for (0..500) |_| {
        matrix.multiplyAssumeValid(x, y);
        clobberPtr(y.ptr);
    }
}

fn benchCscAxSkippingZeros(allocator: std.mem.Allocator, matrix: *const zhighs.matrix.CscMatrix) void {
    const x = allocator.alloc(f64, dimension) catch unreachable;
    defer allocator.free(x);
    @memset(x, 0.0);
    var index: usize = 0;
    while (index < dimension) : (index += 20) x[index] = 1.0;
    const y = allocator.alloc(f64, dimension) catch unreachable;
    defer allocator.free(y);
    for (0..2000) |_| {
        matrix.multiplySkippingZerosAssumeValid(x, y);
        clobberPtr(y.ptr);
    }
}

fn benchCscSparseInput(allocator: std.mem.Allocator, matrix: *const zhighs.matrix.CscMatrix, add_only: bool) void {
    const count = dimension / 20;
    const indices = allocator.alloc(zhighs.ColId, count) catch unreachable;
    defer allocator.free(indices);
    const values = allocator.alloc(f64, count) catch unreachable;
    defer allocator.free(values);
    for (indices, values, 0..) |*id, *value, i| {
        id.* = zhighs.ColId.fromUsizeAssumeValid(i * 20);
        value.* = 1.0;
    }
    const y_storage = allocator.alloc(f64, dimension + 40) catch unreachable;
    defer allocator.free(y_storage);
    const y = y_storage[40..][0..dimension];
    @memset(y, 0.0);
    const view: zhighs.matrix.SparseVectorView(zhighs.ColId) = .{
        .dimension = dimension,
        .indices = indices,
        .values = values,
    };
    for (0..4000) |_| {
        if (add_only)
            matrix.addSparseProductAssumeValid(view, y)
        else
            matrix.multiplySparseAssumeValid(view, y);
        clobberPtr(y.ptr);
    }
}

fn benchAlphaAxPy(allocator: std.mem.Allocator, matrix: *const zhighs.matrix.CscMatrix, x: []const f64) void {
    const y = allocator.alloc(f64, dimension) catch unreachable;
    defer allocator.free(y);
    @memset(y, 0);
    for (0..500) |_| {
        zhighs.matrix.addProductAssumeValid(matrix.*, 1.0, x, y);
        clobberPtr(y.ptr);
    }
}

fn benchScaling(allocator: std.mem.Allocator, matrix: *zhighs.matrix.CscMatrix) void {
    const row_scale = allocator.alloc(f64, dimension) catch unreachable;
    defer allocator.free(row_scale);
    const col_scale = allocator.alloc(f64, dimension) catch unreachable;
    defer allocator.free(col_scale);
    @memset(row_scale, 1.0);
    @memset(col_scale, 1.0);
    for (0..500) |_| {
        zhighs.matrix.applyScalingAssumeValid(matrix, .{ .row = row_scale, .col = col_scale });
        clobberPtr(matrix.values.ptr);
    }
}

fn benchProductQuad(allocator: std.mem.Allocator, matrix: *const zhighs.matrix.CscMatrix, x: []const f64) void {
    const y_storage = allocator.alloc(f64, dimension + 40) catch unreachable;
    defer allocator.free(y_storage);
    const y = y_storage[40..][0..dimension];
    const scratch_storage = allocator.alloc(zhighs.HCD, dimension + 12) catch unreachable;
    defer allocator.free(scratch_storage);
    const scratch = scratch_storage[12..][0..dimension];
    for (0..500) |_| {
        zhighs.matrix.multiplyCompensatedAssumeValid(matrix.*, x, y, scratch);
        clobberPtr(y.ptr);
    }
}

fn benchTransposeInto(allocator: std.mem.Allocator, matrix: *const zhighs.matrix.CscMatrix, cursor: []zhighs.HUInt) void {
    const starts = allocator.alloc(usize, matrix.num_cols + 1) catch unreachable;
    defer allocator.free(starts);
    const rows = allocator.alloc(zhighs.RowId, nnz) catch unreachable;
    defer allocator.free(rows);
    const values = allocator.alloc(f64, nnz) catch unreachable;
    defer allocator.free(values);
    for (0..500) |_| {
        zhighs.matrix.transposeIntoAssumeValid(matrix.*, starts, rows, values, cursor) catch unreachable;
        clobberPtr(values.ptr);
    }
}

fn benchBuilderFreezeGeneral(allocator: std.mem.Allocator) void {
    var builder = zhighs.matrix.MatrixBuilder.init(dimension, dimension) catch unreachable;
    defer builder.deinit(allocator);
    builder.reserve(allocator, nnz) catch unreachable;
    for (0..100) |_| {
        builder.clearRetainingCapacity();
        var col = dimension;
        while (col != 0) {
            col -= 1;
            const col_id = zhighs.ColId.fromUsize(col) catch unreachable;
            if (col + 1 < dimension) builder.append(allocator, zhighs.RowId.fromUsize(col + 1) catch unreachable, col_id, -1.0) catch unreachable;
            builder.append(allocator, zhighs.RowId.fromUsize(col) catch unreachable, col_id, 4.0) catch unreachable;
            if (col != 0) builder.append(allocator, zhighs.RowId.fromUsize(col - 1) catch unreachable, col_id, -1.0) catch unreachable;
        }
        var matrix = builder.freeze(allocator, 0.0) catch unreachable;
        matrix.deinit(allocator);
    }
}

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;

    // Build matrix once, reuse for all kernels
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

    const dense_x = try allocator.alloc(f64, dimension);
    defer allocator.free(dense_x);
    @memset(dense_x, 1.0);

    const requested = init.environ_map.get("ZHIGHS_PERF_KERNEL") orelse "sparse_accumulate";
    if (std.mem.eql(u8, requested, "csc_to_csr_into")) return benchCscToCsrInto(allocator, &matrix, cursor);
    if (std.mem.eql(u8, requested, "sparse_accumulate")) return benchSparseAccumulate(allocator);
    if (std.mem.eql(u8, requested, "csr_atx")) return benchCsrTransposeMultiply(allocator, &csr, dense_x);
    if (std.mem.eql(u8, requested, "csc_atx")) return benchCscTransposeMultiply(allocator, &matrix, dense_x);
    if (std.mem.eql(u8, requested, "csc_ax")) return benchCscAxDense(allocator, &matrix, dense_x);
    if (std.mem.eql(u8, requested, "csc_ax_skip")) return benchCscAxSkippingZeros(allocator, &matrix);
    if (std.mem.eql(u8, requested, "csc_ax_sparse_view")) return benchCscSparseInput(allocator, &matrix, false);
    if (std.mem.eql(u8, requested, "csc_sparse_add")) return benchCscSparseInput(allocator, &matrix, true);
    if (std.mem.eql(u8, requested, "alpha_ax_py")) return benchAlphaAxPy(allocator, &matrix, dense_x);
    if (std.mem.eql(u8, requested, "product_quad")) return benchProductQuad(allocator, &matrix, dense_x);
    if (std.mem.eql(u8, requested, "scaling")) return benchScaling(allocator, &matrix);
    if (std.mem.eql(u8, requested, "transpose_into")) return benchTransposeInto(allocator, &matrix, cursor);
    if (std.mem.eql(u8, requested, "builder_general")) return benchBuilderFreezeGeneral(allocator);
    std.debug.print("unknown ZHIGHS_PERF_KERNEL={s}\n", .{requested});
    return error.InvalidKernel;
}
