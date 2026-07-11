//! Dedicated profiling harness for perf analysis.
//! Build: zig build-exe bench/matrix/perf_profile.zig -I src/ -I lib/ -Doptimize=ReleaseFast -Dcpu=native --name perf-profile
//! Usage: perf record ./perf-profile <kernel-name>

const std = @import("std");
const zhighs = @import("zhighs");

const dimension: usize = 50_000;
const nnz: usize = dimension * 3 - 2;

fn fillSorted(builder: *zhighs.matrix.MatrixBuilder, allocator: std.mem.Allocator) !void {
    for (0..dimension) |col| {
        const col_id = try zhighs.ColId.fromUsize(col);
        if (col != 0) try builder.append(allocator, try zhighs.RowId.fromUsize(col - 1), col_id, -1.0);
        try builder.append(allocator, try zhighs.RowId.fromUsize(col), col_id, 4.0);
        if (col + 1 < dimension) try builder.append(allocator, try zhighs.RowId.fromUsize(col + 1), col_id, -1.0);
    }
}

fn benchCscToCsrInto(allocator: std.mem.Allocator, matrix: *const zhighs.matrix.CscMatrix, cursor: []usize) void {
    const reusable_starts = allocator.alloc(usize, dimension + 1) catch unreachable;
    defer allocator.free(reusable_starts);
    const reusable_cols = allocator.alloc(zhighs.ColId, nnz) catch unreachable;
    defer allocator.free(reusable_cols);
    const reusable_values = allocator.alloc(f64, nnz) catch unreachable;
    defer allocator.free(reusable_values);
    for (0..500) |_| {
        zhighs.matrix.fillCsrFromCscAssumeValid(matrix, reusable_starts, reusable_cols, reusable_values, cursor) catch unreachable;
        std.mem.doNotOptimizeAway(reusable_values.ptr);
    }
}

fn benchSparseAccumulate(allocator: std.mem.Allocator) void {
    var accumulator = zhighs.matrix.SparseAccumulator(zhighs.RowId).init(allocator, dimension) catch unreachable;
    defer accumulator.deinit(allocator);
    accumulator.reserve(allocator, dimension) catch unreachable;
    for (0..500) |_| {
        accumulator.clear();
        for (0..dimension) |index| {
            const id = zhighs.RowId.fromUsizeAssumeValid(index);
            accumulator.addAssumeValid(id, 1.0);
            accumulator.addAssumeValid(id, -0.5);
        }
        std.mem.doNotOptimizeAway(&accumulator);
    }
    _ = accumulator.get(zhighs.RowId.fromUsizeAssumeValid(dimension / 2));
}

fn benchCsrTransposeMultiply(allocator: std.mem.Allocator, csr: *const zhighs.matrix.CsrView, x: []const f64) void {
    const y = allocator.alloc(f64, dimension) catch unreachable;
    defer allocator.free(y);
    for (0..500) |_| {
        @memset(y, 0);
        csr.transposeMultiplyAssumeValid(x, y);
        std.mem.doNotOptimizeAway(y.ptr);
    }
}

fn benchCscTransposeMultiply(allocator: std.mem.Allocator, matrix: *const zhighs.matrix.CscMatrix, x: []const f64) void {
    const y = allocator.alloc(f64, dimension) catch unreachable;
    defer allocator.free(y);
    for (0..500) |_| {
        @memset(y, 0);
        matrix.transposeMultiplyAssumeValid(x, y);
        std.mem.doNotOptimizeAway(y.ptr);
    }
}

fn benchCscAxDense(allocator: std.mem.Allocator, matrix: *const zhighs.matrix.CscMatrix, x: []const f64) void {
    const y = allocator.alloc(f64, dimension) catch unreachable;
    defer allocator.free(y);
    for (0..500) |_| {
        @memset(y, 0);
        matrix.multiplyAssumeValid(x, y);
        std.mem.doNotOptimizeAway(y.ptr);
    }
}

fn benchAlphaAxPy(allocator: std.mem.Allocator, matrix: *const zhighs.matrix.CscMatrix, x: []const f64) void {
    const y = allocator.alloc(f64, dimension) catch unreachable;
    defer allocator.free(y);
    for (0..500) |_| {
        @memset(y, 0);
        zhighs.matrix.addProductAssumeValid(matrix.*, 1.0, x, y);
        std.mem.doNotOptimizeAway(y.ptr);
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
        std.mem.doNotOptimizeAway(matrix.values.ptr);
    }
}

fn benchTransposeInto(allocator: std.mem.Allocator, matrix: *const zhighs.matrix.CscMatrix, cursor: []usize) void {
    const starts = allocator.alloc(usize, matrix.num_cols + 1) catch unreachable;
    defer allocator.free(starts);
    const rows = allocator.alloc(zhighs.RowId, nnz) catch unreachable;
    defer allocator.free(rows);
    const values = allocator.alloc(f64, nnz) catch unreachable;
    defer allocator.free(values);
    for (0..500) |_| {
        zhighs.matrix.transposeIntoAssumeValid(matrix, starts, rows, values, cursor) catch unreachable;
        std.mem.doNotOptimizeAway(values.ptr);
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

/// Change this to profile a specific kernel:
///   .sparse_accumulate, .csc_to_csr_into, .csr_atx, .csc_atx,
///   .csc_ax, .alpha_ax_py, .scaling, .transpose_into, .builder_general
const KernelToProfile: enum {
    all,
    sparse_accumulate,
    csc_to_csr_into,
    csr_atx,
    csc_atx,
    csc_ax,
    alpha_ax_py,
    scaling,
    transpose_into,
    builder_general,
} = .sparse_accumulate;

pub fn main() !void {
    const allocator = std.heap.smp_allocator;

    // Build matrix once, reuse for all kernels
    var builder = try zhighs.matrix.MatrixBuilder.init(dimension, dimension);
    defer builder.deinit(allocator);
    try builder.reserve(allocator, nnz);
    try fillSorted(&builder, allocator);
    var matrix = try builder.freezeSortedAssumeValid(allocator, 0.0);
    defer matrix.deinit(allocator);

    const cursor = try allocator.alloc(usize, dimension);
    defer allocator.free(cursor);
    var csr_cache = try zhighs.matrix.CsrCache.buildWithScratchAssumeValid(allocator, matrix, 0, cursor);
    defer csr_cache.deinit(allocator);
    const csr = csr_cache.viewAssumeCurrent();

    const dense_x = try allocator.alloc(f64, dimension);
    defer allocator.free(dense_x);
    @memset(dense_x, 1.0);

    const ProfileType = @TypeOf(KernelToProfile);
    inline for (comptime std.meta.tags(ProfileType)) |k| {
        if (k == .all) continue;
        if (KernelToProfile != .all and k != KernelToProfile) continue;

        if (k == .csc_to_csr_into) {
            benchCscToCsrInto(allocator, &matrix, cursor);
        }
        if (k == .sparse_accumulate) {
            benchSparseAccumulate(allocator);
        }
        if (k == .csr_atx) {
            benchCsrTransposeMultiply(allocator, &csr, dense_x);
        }
        if (k == .csc_atx) {
            benchCscTransposeMultiply(allocator, &matrix, dense_x);
        }
        if (k == .csc_ax) {
            benchCscAxDense(allocator, &matrix, dense_x);
        }
        if (k == .alpha_ax_py) {
            benchAlphaAxPy(allocator, &matrix, dense_x);
        }
        if (k == .scaling) {
            benchScaling(allocator, &matrix);
        }
        if (k == .transpose_into) {
            benchTransposeInto(allocator, &matrix, cursor);
        }
        if (k == .builder_general) {
            benchBuilderFreezeGeneral(allocator);
        }
    }
}
