//! Assembly-optimized CSR multiply kernels.
//!
//! ## CSR A^T x (scatter-add)
//!
//! The inner loop is a scatter-add: y[col] += vs[pos] * mult where col is
//! loaded from col_indices[pos]. LLVM's alias analysis conservatively assumes
//! that a store to y[col] may alias the values, col_indices, or row_starts
//! pointers — forcing a reload of every pointer after each scatter iteration.
//! Hand-tuned assembly keeps all pointers in registers for the full duration.
//!
//! Uses SSE2 instructions (movsd/mulsd/addsd, no VEX prefix) for maximum
//! compatibility across all x86-64 CPUs.
//!
//! ## Fallback
//!
//! Non-x86-64 targets fall back to a plain Zig loop.

const builtin = @import("builtin");
const std = @import("std");

const have_sse2: bool = builtin.cpu.arch == .x86_64;

/// CSR A^T x: y += A^T * x  (scatter-add)
///
/// Parameters are raw addresses via @intFromPtr to ensure direct register input.
pub fn csrTransposeMultiply(
    nrow: usize,
    rs_addr: usize,
    ci_addr: usize,
    vs_addr: usize,
    xp_addr: usize,
    yp_addr: usize,
) void {
    if (nrow == 0) return;

    if (have_sse2) {
        var ri: usize = 0;
        var pp: usize = 0;
        asm volatile (
            \\1:
            \\cmpq %[nrow], %[ri]
            \\jae 3f
            \\movsd (%[xp], %[ri], 8), %%xmm6
            \\movq (%[rs], %[ri], 8), %[pos]
            \\movq 8(%[rs], %[ri], 8), %%rcx
            \\2:
            \\cmpq %%rcx, %[pos]
            \\jae 0f
            \\movl (%[ci], %[pos], 4), %%r8d
            \\movsd (%[vs], %[pos], 8), %%xmm0
            \\mulsd %%xmm6, %%xmm0
            \\addsd (%[yp], %%r8, 8), %%xmm0
            \\movsd %%xmm0, (%[yp], %%r8, 8)
            \\incq %[pos]
            \\jmp 2b
            \\0:
            \\incq %[ri]
            \\jmp 1b
            \\3:
            : [ri] "+r"(ri), [pos] "+r"(pp)
            : [rs] "r"(rs_addr), [ci] "r"(ci_addr), [vs] "r"(vs_addr),
              [xp] "r"(xp_addr), [yp] "r"(yp_addr), [nrow] "r"(nrow)
            : .{ .rcx = true, .r8 = true, .xmm0 = true, .xmm6 = true, .cc = true, .memory = true }
        );
    } else {
        var i: usize = 0;
        while (i < nrow) : (i += 1) {
            const rs_s = @as([*]const usize, @ptrFromInt(rs_addr));
            const ci_s = @as([*]const u32, @ptrFromInt(ci_addr));
            const vs_s = @as([*]const f64, @ptrFromInt(vs_addr));
            const xp_s = @as([*]const f64, @ptrFromInt(xp_addr));
            const yp_s = @as([*]f64, @ptrFromInt(yp_addr));
            const mult = xp_s[i];
            var pos = rs_s[i];
            const end = rs_s[i + 1];
            while (pos < end) : (pos += 1) {
                yp_s[ci_s[pos]] += vs_s[pos] * mult;
            }
        }
    }
}

const testing = std.testing;

test "csrTransposeMultiply matches CSC reference" {
    const foundation = @import("foundation");
    const RowId = foundation.RowId;
    const csc = @import("../csc.zig");
    const csr_cache = @import("../csr_view.zig");

    const ci = [_]u32{ 0, 2, 1, 0, 2 };
    const vs = [_]f64{ 2.0, 3.0, 4.0, -1.0, 5.0 };
    const x = [_]f64{ 1.0, 2.0, 3.0 };
    const matrix: csc.CscMatrix = .{
        .num_rows = 3, .num_cols = 3,
        .col_starts = @constCast(@as(*const [4]usize, &[_]usize{0, 2, 3, 5})),
        .row_indices = @constCast(@as(*const [5]RowId, @ptrCast(&ci))),
        .values = @constCast(&vs),
    };
    var cursor = [_]usize{0} ** 3;
    var cache = try csr_cache.CsrCache.buildWithScratchAssumeValid(testing.allocator, matrix, 0, &cursor);
    defer cache.deinit(testing.allocator);
    const csr = cache.viewAssumeCurrent();

    var ref_y: [3]f64 = undefined;
    @memset(&ref_y, 0);
    csr.transposeMultiplyAssumeValid(&x, &ref_y);

    var asm_y: [3]f64 = undefined;
    @memset(&asm_y, 0);
    csrTransposeMultiply(csr.num_rows,
        @intFromPtr(csr.row_starts.ptr), @intFromPtr(csr.col_indices.ptr),
        @intFromPtr(csr.values.ptr), @intFromPtr(&x), @intFromPtr(&asm_y));
    try testing.expectEqualSlices(f64, &ref_y, &asm_y);
}

test "csrTransposeMultiply with zero rows" {
    var y: [3]f64 = undefined;
    @memset(&y, 0);
    csrTransposeMultiply(0, 0, 0, 0, 0, @intFromPtr(&y));
}

test "csrTransposeMultiply large regular pattern" {
    const dimension: usize = 100;
    const nnz = dimension * 3 - 2;
    const foundation = @import("foundation");
    const RowId = foundation.RowId;
    const ColId = foundation.ColId;
    const builder_module = @import("../builder.zig");
    const csr_cache = @import("../csr_view.zig");

    var builder = try builder_module.MatrixBuilder.init(dimension, dimension);
    defer builder.deinit(testing.allocator);
    try builder.reserve(testing.allocator, nnz);
    for (0..dimension) |col| {
        const col_id = try ColId.fromUsize(col);
        if (col != 0) try builder.append(testing.allocator, try RowId.fromUsize(col - 1), col_id, -1.0);
        try builder.append(testing.allocator, try RowId.fromUsize(col), col_id, 4.0);
        if (col + 1 < dimension) try builder.append(testing.allocator, try RowId.fromUsize(col + 1), col_id, -1.0);
    }
    var matrix = try builder.freezeSortedAssumeValid(testing.allocator, 0.0);
    defer matrix.deinit(testing.allocator);

    const cursor = try testing.allocator.alloc(usize, dimension);
    defer testing.allocator.free(cursor);
    var cache = try csr_cache.CsrCache.buildWithScratchAssumeValid(testing.allocator, matrix, 0, cursor);
    defer cache.deinit(testing.allocator);
    const csr = cache.viewAssumeCurrent();

    const x = try testing.allocator.alloc(f64, dimension);
    defer testing.allocator.free(x);
    for (x, 0..) |*v, i| v.* = @floatFromInt(i % 7);

    const ref_y = try testing.allocator.alloc(f64, dimension);
    defer testing.allocator.free(ref_y);
    @memset(ref_y, 0);
    csr.transposeMultiplyAssumeValid(x, ref_y);

    const asm_y = try testing.allocator.alloc(f64, dimension);
    defer testing.allocator.free(asm_y);
    @memset(asm_y, 0);
    csrTransposeMultiply(csr.num_rows,
        @intFromPtr(csr.row_starts.ptr), @intFromPtr(csr.col_indices.ptr),
        @intFromPtr(csr.values.ptr), @intFromPtr(x.ptr), @intFromPtr(asm_y.ptr));
    try testing.expectEqualSlices(f64, ref_y, asm_y);
}
