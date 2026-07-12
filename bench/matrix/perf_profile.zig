//! Dedicated profiling harness for perf analysis.
//! Build: zig build build-perf-profile -Doptimize=ReleaseFast -Dcpu=native
//! Usage:
//!   ZHIGHS_PERF_KERNEL=csc_to_csr_into zig-out/bin/perf-profile
//!
//! Each process measures exactly one kernel. Kernel names and repeat counts
//! are shared with highs_perf_profile.cpp for fair comparison.

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

fn nowNs() i128 {
    var ts: std.posix.timespec = undefined;
    switch (std.posix.errno(std.posix.system.clock_gettime(.MONOTONIC, &ts))) {
        .SUCCESS => {},
        else => unreachable,
    }
    return @as(i128, ts.sec) * std.time.ns_per_s + ts.nsec;
}

fn checksum(values: []const f64) f64 {
    var sum: f64 = 0.0;
    for (values) |v| sum += v;
    return sum;
}

/// Structural hash of a sparse matrix: starts + indices + values.
/// FNV-1a 64-bit. All integer fields are normalized to u64 so the hash
/// matches the C++ reference regardless of the native storage width.
fn structuralHash(starts: []const usize, indices: anytype, values: []const f64) u64 {
    var h: u64 = 0xcbf29ce484222325;
    for (starts) |s| {
        const v: u64 = @intCast(s);
        const bytes = std.mem.asBytes(&v);
        for (bytes) |b| { h ^= b; h *%= 0x100000001b3; }
    }
    for (indices) |idx| {
        const v: u64 = @intCast(idx.toUsize());
        const bytes = std.mem.asBytes(&v);
        for (bytes) |b| { h ^= b; h *%= 0x100000001b3; }
    }
    for (values) |v| {
        const bytes = std.mem.asBytes(&v);
        for (bytes) |b| { h ^= b; h *%= 0x100000001b3; }
    }
    return h;
}

/// Structural hash for CSR output (row_starts as HUInt), FNV-1a 64-bit.
/// Normalizes HUInt to u64 for cross-implementation comparison.
fn structuralHashCsr(row_starts: []const zhighs.HUInt, col_indices: anytype, values: []const f64) u64 {
    var h: u64 = 0xcbf29ce484222325;
    for (row_starts) |s| {
        const v: u64 = @intCast(s);
        const bytes = std.mem.asBytes(&v);
        for (bytes) |b| { h ^= b; h *%= 0x100000001b3; }
    }
    for (col_indices) |idx| {
        const v: u64 = @intCast(idx.toUsize());
        const bytes = std.mem.asBytes(&v);
        for (bytes) |b| { h ^= b; h *%= 0x100000001b3; }
    }
    for (values) |v| {
        const bytes = std.mem.asBytes(&v);
        for (bytes) |b| { h ^= b; h *%= 0x100000001b3; }
    }
    return h;
}

/// Build the tridiagonal test matrix and reusable workspace once.
const TestHarness = struct {
    matrix: zhighs.matrix.CscMatrix,
    csr_cache: zhighs.matrix.CsrCache,
    csr_view: zhighs.matrix.CsrView,
    csr_buffers: zhighs.matrix.CsrBuffers,
    transpose_buffers: zhighs.matrix.TransposeBuffers,
    dense_x: []f64,
    sparse_x: []f64,
    sparse_ids: []zhighs.ColId,
    sparse_values: []f64,
    sparse_view: zhighs.matrix.SparseVectorView(zhighs.ColId),
    y: []f64,
    hcd_scratch: []zhighs.HCD,
    row_scale: []f64,
    col_scale: []f64,
    cursor: []zhighs.HUInt,
    allocator: std.mem.Allocator,

    fn init(allocator: std.mem.Allocator) !TestHarness {
        var builder = try zhighs.matrix.MatrixBuilder.init(dimension, dimension);
        defer builder.deinit(allocator);
        try builder.reserve(allocator, nnz);
        for (0..dimension) |col| {
            const col_id = try zhighs.ColId.fromUsize(col);
            if (col != 0) try builder.append(allocator, try zhighs.RowId.fromUsize(col - 1), col_id, -1.0);
            try builder.append(allocator, try zhighs.RowId.fromUsize(col), col_id, 4.0);
            if (col + 1 < dimension) try builder.append(allocator, try zhighs.RowId.fromUsize(col + 1), col_id, -1.0);
        }
        const matrix = try builder.freezeSortedAssumeValid(allocator, 0.0);

        const cursor_storage = try allocator.alloc(zhighs.HUInt, dimension + 48);
        const cursor = cursor_storage[48..][0..dimension];
        const csr_buffers = try zhighs.matrix.CsrBuffers.init(allocator, dimension, nnz);
        const transpose_buffers = try zhighs.matrix.TransposeBuffers.init(allocator, dimension, nnz);
        const csr_cache = try zhighs.matrix.CsrCache.buildWithScratchAssumeValid(allocator, matrix, 0, csr_buffers.cursor);
        const csr_view = csr_cache.viewAssumeCurrent();

        const dense_x = try allocator.alloc(f64, dimension);
        @memset(dense_x, 1.0);
        const sparse_x = try allocator.alloc(f64, dimension);
        for (sparse_x, 0..) |*v, i| v.* = if (i % 20 == 0) 1.0 else 0.0;
        const sparse_count = dimension / 20;
        const sparse_ids = try allocator.alloc(zhighs.ColId, sparse_count);
        const sparse_values_arr = try allocator.alloc(f64, sparse_count);
        for (sparse_ids, sparse_values_arr, 0..) |*id, *v, i| {
            id.* = try zhighs.ColId.fromUsize(i * 20);
            v.* = 1.0;
        }
        const sparse_view: zhighs.matrix.SparseVectorView(zhighs.ColId) = .{
            .dimension = dimension,
            .indices = sparse_ids,
            .values = sparse_values_arr,
        };
        const y = try allocator.alloc(f64, dimension);
        const hcd_storage = try allocator.alloc(zhighs.HCD, dimension + 12);
        const hcd_scratch = hcd_storage[12..][0..dimension];
        const row_scale = try allocator.alloc(f64, dimension);
        const col_scale = try allocator.alloc(f64, dimension);
        @memset(row_scale, 1.0);
        @memset(col_scale, 1.0);

        return TestHarness{
            .matrix = matrix,
            .csr_cache = csr_cache,
            .csr_view = csr_view,
            .csr_buffers = csr_buffers,
            .transpose_buffers = transpose_buffers,
            .dense_x = dense_x,
            .sparse_x = sparse_x,
            .sparse_ids = sparse_ids,
            .sparse_values = sparse_values_arr,
            .sparse_view = sparse_view,
            .y = y,
            .hcd_scratch = hcd_scratch,
            .row_scale = row_scale,
            .col_scale = col_scale,
            .cursor = cursor,
            .allocator = allocator,
        };
    }

    fn deinit(self: *TestHarness) void {
        const a = self.allocator;
        self.matrix.deinit(a);
        self.csr_cache.deinit(a);
        self.csr_buffers.deinit(a);
        self.transpose_buffers.deinit(a);
        a.free(self.dense_x);
        a.free(self.sparse_x);
        a.free(self.sparse_ids);
        a.free(self.sparse_values);
        a.free(self.y);
        a.free(self.hcd_scratch);
        a.free(self.row_scale);
        a.free(self.col_scale);
    }
};

/// Select repeats so each kernel runs roughly 300ms–2s.
fn warmRepeats(kernel: []const u8) usize {
    return if (std.mem.eql(u8, kernel, "csc_to_csr_into") or std.mem.eql(u8, kernel, "csc_to_csr_owning") or std.mem.eql(u8, kernel, "transpose_into") or std.mem.eql(u8, kernel, "transpose_owning") or std.mem.eql(u8, kernel, "builder_freeze_general"))
        100
    else if (std.mem.eql(u8, kernel, "builder_freeze_sorted") or std.mem.eql(u8, kernel, "builder_freeze_prepopulated") or std.mem.eql(u8, kernel, "builder_freeze_canonical") or std.mem.eql(u8, kernel, "builder_freeze_reusable"))
        1000
    else if (std.mem.eql(u8, kernel, "clear_output"))
        100000
    else if (std.mem.eql(u8, kernel, "csc_ax_sparse_view") or std.mem.eql(u8, kernel, "csc_sparse_add_no_clear"))
        4000
    else if (std.mem.eql(u8, kernel, "csc_ax_sparse_skip"))
        2000
    else
        500;
}

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;
    var h = try TestHarness.init(allocator);
    defer h.deinit();

    const requested = init.environ_map.get("ZHIGHS_PERF_KERNEL") orelse {
        std.debug.print("Usage: ZHIGHS_PERF_KERNEL=<name> {s}\n", .{@src().file});
        std.debug.print("Available kernels:\n", .{});
        for (kernels) |k| std.debug.print("  {s}\n", .{k});
        return error.NoKernelSpecified;
    };

    const repeats = warmRepeats(requested);
    var result_checksum: f64 = 0.0;
    var result_struct_hash: u64 = 0;
    const start = nowNs();

    if (std.mem.eql(u8, requested, "clear_output")) {
        for (0..repeats) |_| {
            zhighs.matrix.clearF64(h.y);
            clobberPtr(h.y.ptr);
        }
        result_checksum = checksum(h.y);
    } else if (std.mem.eql(u8, requested, "csc_ax_dense")) {
        for (0..repeats) |_| {
            h.matrix.multiplyAssumeValid(h.dense_x, h.y);
            clobberPtr(h.y.ptr);
        }
        result_checksum = checksum(h.y);
    } else if (std.mem.eql(u8, requested, "csc_ax_sparse_skip")) {
        for (0..repeats) |_| {
            h.matrix.multiplySkippingZerosAssumeValid(h.sparse_x, h.y);
            clobberPtr(h.y.ptr);
        }
        result_checksum = checksum(h.y);
    } else if (std.mem.eql(u8, requested, "csc_ax_sparse_view")) {
        for (0..repeats) |_| {
            h.matrix.multiplySparseAssumeValid(h.sparse_view, h.y);
            clobberPtr(h.y.ptr);
        }
        result_checksum = checksum(h.y);
    } else if (std.mem.eql(u8, requested, "csc_sparse_add_no_clear")) {
        @memset(h.y, 0.0);
        for (0..repeats) |_| {
            h.matrix.addSparseProductAssumeValid(h.sparse_view, h.y);
            clobberPtr(h.y.ptr);
        }
        result_checksum = checksum(h.y);
    } else if (std.mem.eql(u8, requested, "csr_ax_dense")) {
        for (0..repeats) |_| {
            h.csr_view.multiplyAssumeValid(h.dense_x, h.y);
            clobberPtr(h.y.ptr);
        }
        result_checksum = checksum(h.y);
    } else if (std.mem.eql(u8, requested, "csc_atx_dense")) {
        for (0..repeats) |_| {
            h.matrix.transposeMultiplyAssumeValid(h.dense_x, h.y);
            clobberPtr(h.y.ptr);
        }
        result_checksum = checksum(h.y);
    } else if (std.mem.eql(u8, requested, "csr_atx_dense")) {
        for (0..repeats) |_| {
            h.csr_view.transposeMultiplyAssumeValid(h.dense_x, h.y);
            clobberPtr(h.y.ptr);
        }
        result_checksum = checksum(h.y);
    } else if (std.mem.eql(u8, requested, "alpha_ax_plus_y")) {
        @memset(h.y, 0.0);
        for (0..repeats) |_| {
            zhighs.matrix.addProductAssumeValid(h.matrix, 1.0, h.dense_x, h.y);
            clobberPtr(h.y.ptr);
        }
        result_checksum = checksum(h.y);
    } else if (std.mem.eql(u8, requested, "product_quad")) {
        for (0..repeats) |_| {
            zhighs.matrix.multiplyCompensatedAssumeValid(h.matrix, h.dense_x, h.y, h.hcd_scratch);
            clobberPtr(h.y.ptr);
        }
        result_checksum = checksum(h.y);
    } else if (std.mem.eql(u8, requested, "apply_scale")) {
        for (0..repeats) |_| {
            zhighs.matrix.applyScalingAssumeValid(&h.matrix, .{ .row = h.row_scale, .col = h.col_scale });
            clobberPtr(h.matrix.values.ptr);
        }
        result_checksum = checksum(h.matrix.values);
    } else if (std.mem.eql(u8, requested, "csc_to_csr_into")) {
        for (0..repeats) |_| {
            try zhighs.matrix.fillCsrFromCscAssumeValid(h.matrix, h.csr_buffers.row_starts, h.csr_buffers.col_indices, h.csr_buffers.values, h.csr_buffers.cursor);
            clobberPtr(h.csr_buffers.values.ptr);
        }
        result_checksum = checksum(h.csr_buffers.values);
        result_struct_hash = structuralHashCsr(h.csr_buffers.row_starts, h.csr_buffers.col_indices, h.csr_buffers.values);
    } else if (std.mem.eql(u8, requested, "csc_to_csr_owning")) {
        // Match C++ semantics: copy the CSC matrix (like HighsSparseMatrix fresh = csc),
        // then build CSR from the copy (like fresh.ensureRowwise()).
        for (0..repeats) |_| {
            const starts_copy = try h.allocator.dupe(usize, h.matrix.col_starts);
            defer h.allocator.free(starts_copy);
            const indices_copy = try h.allocator.dupe(zhighs.RowId, h.matrix.row_indices);
            defer h.allocator.free(indices_copy);
            const values_copy = try h.allocator.dupe(f64, h.matrix.values);
            defer h.allocator.free(values_copy);
            const matrix_copy = zhighs.matrix.CscMatrix{
                .num_rows = h.matrix.num_rows,
                .num_cols = h.matrix.num_cols,
                .col_starts = starts_copy,
                .row_indices = indices_copy,
                .values = values_copy,
                .compact_col_starts = null,
            };
            var cache = try zhighs.matrix.CsrCache.buildWithScratchAssumeValid(h.allocator, matrix_copy, 0, h.csr_buffers.cursor);
            clobberPtr(cache.values.ptr);
            cache.deinit(h.allocator);
        }
        result_checksum = checksum(h.csr_view.values);
        // Build one extra for struct hash (must be on the owning CSR, not the reusable buffers)
        {
            const starts_copy = try h.allocator.dupe(usize, h.matrix.col_starts);
            defer h.allocator.free(starts_copy);
            const indices_copy = try h.allocator.dupe(zhighs.RowId, h.matrix.row_indices);
            defer h.allocator.free(indices_copy);
            const values_copy = try h.allocator.dupe(f64, h.matrix.values);
            defer h.allocator.free(values_copy);
            const matrix_copy = zhighs.matrix.CscMatrix{
                .num_rows = h.matrix.num_rows,
                .num_cols = h.matrix.num_cols,
                .col_starts = starts_copy,
                .row_indices = indices_copy,
                .values = values_copy,
                .compact_col_starts = null,
            };
            var cache = try zhighs.matrix.CsrCache.buildWithScratchAssumeValid(h.allocator, matrix_copy, 0, h.csr_buffers.cursor);
            defer cache.deinit(h.allocator);
            result_struct_hash = structuralHashCsr(cache.row_starts, cache.col_indices, cache.values);
        }
    } else if (std.mem.eql(u8, requested, "transpose_into")) {
        for (0..repeats) |_| {
            try zhighs.matrix.transposeIntoAssumeValid(h.matrix, h.transpose_buffers.starts, h.transpose_buffers.rows, h.transpose_buffers.values, h.transpose_buffers.cursor);
            clobberPtr(h.transpose_buffers.values.ptr);
        }
        result_checksum = checksum(h.transpose_buffers.values);
        result_struct_hash = structuralHash(h.transpose_buffers.starts, h.transpose_buffers.rows, h.transpose_buffers.values);
    } else if (std.mem.eql(u8, requested, "transpose_owning")) {
        for (0..repeats) |_| {
            var t = try zhighs.matrix.transposeLeanAssumeValid(h.allocator, h.matrix);
            clobberPtr(t.values.ptr);
            t.deinit(h.allocator);
        }
        // Build one extra for checksum + structural hash
        {
            var t = try zhighs.matrix.transposeAssumeValid(h.allocator, h.matrix);
            result_checksum = checksum(t.values);
            result_struct_hash = structuralHash(t.col_starts, t.row_indices, t.values);
            t.deinit(h.allocator);
        }
    } else if (std.mem.eql(u8, requested, "builder_freeze_sorted")) {
        for (0..repeats) |_| {
            var b = try zhighs.matrix.MatrixBuilder.init(dimension, dimension);
            defer b.deinit(h.allocator);
            try b.reserve(h.allocator, nnz);
            for (0..dimension) |col| {
                const col_id = zhighs.ColId.fromUsizeAssumeValid(col);
                if (col != 0) b.appendPreReserved(zhighs.RowId.fromUsizeAssumeValid(col - 1), col_id, -1.0);
                b.appendPreReserved(zhighs.RowId.fromUsizeAssumeValid(col), col_id, 4.0);
                if (col + 1 < dimension) b.appendPreReserved(zhighs.RowId.fromUsizeAssumeValid(col + 1), col_id, -1.0);
            }
            var m = try b.freezeSortedLeanAssumeValid(h.allocator, 0.0);
            clobberPtr(m.values.ptr);
            m.deinit(h.allocator);
        }
        // Build one extra for checksum
        {
            var b = try zhighs.matrix.MatrixBuilder.init(dimension, dimension);
            defer b.deinit(h.allocator);
            try b.reserve(h.allocator, nnz);
            for (0..dimension) |col| {
                const col_id = zhighs.ColId.fromUsizeAssumeValid(col);
                if (col != 0) b.appendPreReserved(zhighs.RowId.fromUsizeAssumeValid(col - 1), col_id, -1.0);
                b.appendPreReserved(zhighs.RowId.fromUsizeAssumeValid(col), col_id, 4.0);
                if (col + 1 < dimension) b.appendPreReserved(zhighs.RowId.fromUsizeAssumeValid(col + 1), col_id, -1.0);
            }
            var m = try b.freezeSortedLeanAssumeValid(h.allocator, 0.0);
            result_checksum = checksum(m.values);
            result_struct_hash = structuralHash(m.col_starts, m.row_indices, m.values);
            m.deinit(h.allocator);
        }
    } else if (std.mem.eql(u8, requested, "builder_freeze_prepopulated")) {
        // Only the freeze step — arrays are pre-built outside the timed loop.
        var prepop_rows = try h.allocator.alloc(zhighs.RowId, nnz);
        defer h.allocator.free(prepop_rows);
        var prepop_cols = try h.allocator.alloc(zhighs.ColId, nnz);
        defer h.allocator.free(prepop_cols);
        var prepop_values = try h.allocator.alloc(f64, nnz);
        defer h.allocator.free(prepop_values);
        {
            var idx: usize = 0;
            for (0..dimension) |col| {
                const col_id = zhighs.ColId.fromUsizeAssumeValid(col);
                if (col != 0) {
                    prepop_rows[idx] = zhighs.RowId.fromUsizeAssumeValid(col - 1);
                    prepop_cols[idx] = col_id;
                    prepop_values[idx] = -1.0;
                    idx += 1;
                }
                prepop_rows[idx] = zhighs.RowId.fromUsizeAssumeValid(col);
                prepop_cols[idx] = col_id;
                prepop_values[idx] = 4.0;
                idx += 1;
                if (col + 1 < dimension) {
                    prepop_rows[idx] = zhighs.RowId.fromUsizeAssumeValid(col + 1);
                    prepop_cols[idx] = col_id;
                    prepop_values[idx] = -1.0;
                    idx += 1;
                }
            }
        }
        for (0..repeats) |_| {
            var m = try zhighs.matrix.freezeFromSortedArraysAssumeValid(h.allocator, dimension, dimension, prepop_rows, prepop_cols, prepop_values, 0.0, false);
            clobberPtr(m.values.ptr);
            m.deinit(h.allocator);
        }
        // Build one extra for checksum + struct hash
        {
            var m = try zhighs.matrix.freezeFromSortedArraysAssumeValid(h.allocator, dimension, dimension, prepop_rows, prepop_cols, prepop_values, 0.0, false);
            result_checksum = checksum(m.values);
            result_struct_hash = structuralHash(m.col_starts, m.row_indices, m.values);
            m.deinit(h.allocator);
        }
    } else if (std.mem.eql(u8, requested, "builder_freeze_canonical")) {
        // Canonical input: no duplicates, no zeros, all finite — single-pass count+memcpy.
        var canon_rows = try h.allocator.alloc(zhighs.RowId, nnz);
        defer h.allocator.free(canon_rows);
        var canon_cols = try h.allocator.alloc(zhighs.ColId, nnz);
        defer h.allocator.free(canon_cols);
        var canon_values = try h.allocator.alloc(f64, nnz);
        defer h.allocator.free(canon_values);
        {
            var idx: usize = 0;
            for (0..dimension) |col| {
                const col_id = zhighs.ColId.fromUsizeAssumeValid(col);
                if (col != 0) {
                    canon_rows[idx] = zhighs.RowId.fromUsizeAssumeValid(col - 1);
                    canon_cols[idx] = col_id;
                    canon_values[idx] = -1.0;
                    idx += 1;
                }
                canon_rows[idx] = zhighs.RowId.fromUsizeAssumeValid(col);
                canon_cols[idx] = col_id;
                canon_values[idx] = 4.0;
                idx += 1;
                if (col + 1 < dimension) {
                    canon_rows[idx] = zhighs.RowId.fromUsizeAssumeValid(col + 1);
                    canon_cols[idx] = col_id;
                    canon_values[idx] = -1.0;
                    idx += 1;
                }
            }
        }
        for (0..repeats) |_| {
            var m = try zhighs.matrix.freezeFromCanonicalArraysAssumeValid(h.allocator, dimension, dimension, canon_rows, canon_cols, canon_values, 0.0, false);
            clobberPtr(m.values.ptr);
            m.deinit(h.allocator);
        }
        {
            var m = try zhighs.matrix.freezeFromCanonicalArraysAssumeValid(h.allocator, dimension, dimension, canon_rows, canon_cols, canon_values, 0.0, false);
            result_checksum = checksum(m.values);
            result_struct_hash = structuralHash(m.col_starts, m.row_indices, m.values);
            m.deinit(h.allocator);
        }
    } else if (std.mem.eql(u8, requested, "builder_freeze_reusable")) {
        // Uses formal CscBuildBuffers API: allocate once, pre-touch, reuse across iterations.
        var bufs = try zhighs.matrix.CscBuildBuffers.initCapacity(h.allocator, dimension, nnz);
        defer bufs.deinit(h.allocator);
        // Pre-touch: write to every page to fault it in (caller policy, not forced by API)
        @memset(bufs.col_starts, 0);
        @memset(std.mem.sliceAsBytes(bufs.row_indices), 0);
        @memset(bufs.values, 0);

        // Pre-build canonical arrays once
        var prep_rows = try h.allocator.alloc(zhighs.RowId, nnz);
        defer h.allocator.free(prep_rows);
        var prep_cols = try h.allocator.alloc(zhighs.ColId, nnz);
        defer h.allocator.free(prep_cols);
        var prep_vals = try h.allocator.alloc(f64, nnz);
        defer h.allocator.free(prep_vals);
        {
            var idx: usize = 0;
            for (0..dimension) |col| {
                const col_id = zhighs.ColId.fromUsizeAssumeValid(col);
                if (col != 0) { prep_rows[idx] = zhighs.RowId.fromUsizeAssumeValid(col - 1); prep_cols[idx] = col_id; prep_vals[idx] = -1.0; idx += 1; }
                prep_rows[idx] = zhighs.RowId.fromUsizeAssumeValid(col); prep_cols[idx] = col_id; prep_vals[idx] = 4.0; idx += 1;
                if (col + 1 < dimension) { prep_rows[idx] = zhighs.RowId.fromUsizeAssumeValid(col + 1); prep_cols[idx] = col_id; prep_vals[idx] = -1.0; idx += 1; }
            }
        }

        for (0..repeats) |_| {
            const view = zhighs.matrix.freezeCanonicalIntoAssumeValid(&bufs, dimension, dimension, prep_rows, prep_cols, prep_vals) catch unreachable;
            clobberPtr(view.values.ptr);
        }
        result_checksum = checksum(bufs.values[0..nnz]);
        result_struct_hash = structuralHash(bufs.col_starts[0 .. dimension + 1], bufs.row_indices[0..nnz], bufs.values[0..nnz]);
    } else if (std.mem.eql(u8, requested, "builder_freeze_general")) {
        for (0..repeats) |_| {
            var b = try zhighs.matrix.MatrixBuilder.init(dimension, dimension);
            defer b.deinit(h.allocator);
            try b.reserve(h.allocator, nnz);
            var col = dimension;
            while (col != 0) {
                col -= 1;
                const col_id = zhighs.ColId.fromUsizeAssumeValid(col);
                if (col + 1 < dimension) b.appendPreReserved(zhighs.RowId.fromUsizeAssumeValid(col + 1), col_id, -1.0);
                b.appendPreReserved(zhighs.RowId.fromUsizeAssumeValid(col), col_id, 4.0);
                if (col != 0) b.appendPreReserved(zhighs.RowId.fromUsizeAssumeValid(col - 1), col_id, -1.0);
            }
            var m = try b.freeze(h.allocator, 0.0);
            clobberPtr(m.values.ptr);
            m.deinit(h.allocator);
        }
        // Build one extra for checksum
        {
            var b = try zhighs.matrix.MatrixBuilder.init(dimension, dimension);
            defer b.deinit(h.allocator);
            try b.reserve(h.allocator, nnz);
            var col = dimension;
            while (col != 0) {
                col -= 1;
                const col_id = zhighs.ColId.fromUsizeAssumeValid(col);
                if (col + 1 < dimension) b.appendPreReserved(zhighs.RowId.fromUsizeAssumeValid(col + 1), col_id, -1.0);
                b.appendPreReserved(zhighs.RowId.fromUsizeAssumeValid(col), col_id, 4.0);
                if (col != 0) b.appendPreReserved(zhighs.RowId.fromUsizeAssumeValid(col - 1), col_id, -1.0);
            }
            var m = try b.freeze(h.allocator, 0.0);
            result_checksum = checksum(m.values);
            result_struct_hash = structuralHash(m.col_starts, m.row_indices, m.values);
            m.deinit(h.allocator);
        }
    } else if (std.mem.eql(u8, requested, "sparse_accumulate")) {
        var accumulator = zhighs.matrix.SparseAccumulator(zhighs.RowId).initWithCapacity(h.allocator, dimension, dimension) catch unreachable;
        defer accumulator.deinit(h.allocator);
        const diag = init.environ_map.get("ZHIGHS_PERF_DIAG") != null;
        if (diag) {
            // Control experiment: output addresses + page offsets for correlation
            const dv_ptr = @intFromPtr(accumulator.dense_values.ptr);
            const ap_ptr = @intFromPtr(accumulator.active_ptr);
            std.debug.print("diag,dense_ptr,{d},{d}\n", .{ dv_ptr, dv_ptr & 4095 });
            std.debug.print("diag,active_ptr,{d},{d}\n", .{ ap_ptr, ap_ptr & 4095 });
            // Three separate batches to test within-process consistency
            for (0..3) |batch| {
                const bstart = nowNs();
                for (0..repeats) |_| {
                    accumulator.clear();
                    for (0..dimension) |index| {
                        const id = zhighs.RowId.fromUsizeAssumeValid(index);
                        accumulator.addAssumeValid(id, 1.0);
                        accumulator.addAssumeValid(id, -0.5);
                    }
                    clobberPtr(&accumulator);
                }
                const btotal: u64 = @intCast(nowNs() - bstart);
                const bper = @as(f64, @floatFromInt(btotal)) / @as(f64, @floatFromInt(repeats));
                std.debug.print("diag,batch,{d},{d:.3}\n", .{ batch, bper });
            }
            // Barrier variant: asm volatile between the two adds
            const bstart2 = nowNs();
            for (0..repeats) |_| {
                accumulator.clear();
                for (0..dimension) |index| {
                    const id = zhighs.RowId.fromUsizeAssumeValid(index);
                    accumulator.addAssumeValid(id, 1.0);
                    asm volatile ("" ::: .{ .memory = true });
                    accumulator.addAssumeValid(id, -0.5);
                }
                clobberPtr(&accumulator);
            }
            const btotal2: u64 = @intCast(nowNs() - bstart2);
            const bper2 = @as(f64, @floatFromInt(btotal2)) / @as(f64, @floatFromInt(repeats));
            std.debug.print("diag,barrier,{d:.3}\n", .{bper2});
        } else {
            for (0..repeats) |_| {
                accumulator.clear();
                for (0..dimension) |index| {
                    const id = zhighs.RowId.fromUsizeAssumeValid(index);
                    accumulator.addAssumeValid(id, 1.0);
                    accumulator.addAssumeValid(id, -0.5);
                }
                clobberPtr(&accumulator);
            }
        }
        result_checksum = accumulator.get(zhighs.RowId.fromUsizeAssumeValid(dimension / 2));
    } else {
        std.debug.print("unknown ZHIGHS_PERF_KERNEL={s}\n", .{requested});
        return error.InvalidKernel;
    }

    const total: u64 = @intCast(nowNs() - start);
    const per_repeat = @as(f64, @floatFromInt(total)) / @as(f64, @floatFromInt(repeats));
    std.debug.print("zig,{s},{d},{d:.3},{d:.17},{d}\n", .{ requested, total, per_repeat, result_checksum, result_struct_hash });
}

const kernels = [_][]const u8{
    "clear_output",
    "csc_ax_dense",
    "csc_ax_sparse_skip",
    "csc_ax_sparse_view",
    "csc_sparse_add_no_clear",
    "csr_ax_dense",
    "csc_atx_dense",
    "csr_atx_dense",
    "alpha_ax_plus_y",
    "product_quad",
    "apply_scale",
    "csc_to_csr_into",
    "csc_to_csr_owning",
    "transpose_into",
    "transpose_owning",
    "builder_freeze_sorted",
    "builder_freeze_prepopulated",
    "builder_freeze_canonical",
    "builder_freeze_reusable",
    "builder_freeze_general",
    "sparse_accumulate",
};
