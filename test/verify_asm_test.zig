const std = @import("std");
const builtin = @import("builtin");

fn zigTransposeMultiply(nrow: usize, rs: [*]const usize, ci: [*]const u32, vs: [*]const f64, x: [*]const f64, y: [*]f64) void {
    var i: usize = 0;
    while (i < nrow) : (i += 1) {
        const mult = x[i];
        var pos = rs[i];
        const end = rs[i + 1];
        while (pos < end) : (pos += 1)
            y[ci[pos]] += vs[pos] * mult;
    }
}
fn asmTransposeMultiply(nrow: usize, rs: usize, ci: usize, vs: usize, xp: usize, yp: usize) void {
    if (builtin.cpu.arch != .x86_64) return;
    var ri: usize = 0;
    var pp: usize = 0;
    asm volatile (
        \\cmpq %[nrow], %[ri]
        \\jae 3f
        \\movsd (%[xp], %[ri], 8), %%xmm6
        \\movq (%[rs], %[ri], 8), %[pos]
        \\movq 8(%[rs], %[ri], 8), %%rcx
        \\1:
        \\cmpq %%rcx, %[pos]
        \\jae 0f
        \\movl (%[ci], %[pos], 4), %%r8d
        \\movsd (%[vs], %[pos], 8), %%xmm0
        \\mulsd %%xmm6, %%xmm0
        \\addsd (%[yp], %%r8, 8), %%xmm0
        \\movsd %%xmm0, (%[yp], %%r8, 8)
        \\incq %[pos]
        \\jmp 1b
        \\0:
        \\incq %[ri]
        \\jmp 1b
        \\3:
        : [ri] "+r"(ri), [pos] "+r"(pp)
        : [rs] "r"(rs), [ci] "r"(ci), [vs] "r"(vs),
          [xp] "r"(xp), [yp] "r"(yp), [nrow] "r"(nrow)
        : .{ .rcx = true, .r8 = true, .xmm0 = true, .xmm6 = true, .cc = true, .memory = true }
    );
}

pub fn main() !void {
    const dim: usize = 100;
    var rs: [101]usize = undefined;
    var ci: [298]u32 = undefined;
    var vs: [298]f64 = undefined;
    var pos: usize = 0;
    for (0..dim) |row| {
        rs[row] = pos;
        if (row > 0) { ci[pos] = @intCast(row-1); vs[pos] = -1.0; pos += 1; }
        ci[pos] = @intCast(row); vs[pos] = 4.0; pos += 1;
        if (row + 1 < dim) { ci[pos] = @intCast(row+1); vs[pos] = -1.0; pos += 1; }
    }
    rs[dim] = pos;
    var x: [100]f64 = undefined;
    for (&x, 0..) |*v, i| v.* = @floatFromInt((i * 7) % 13);
    var y_zig: [100]f64 = undefined; @memset(&y_zig, 0);
    var y_asm: [100]f64 = undefined; @memset(&y_asm, 0);

    zigTransposeMultiply(dim, &rs, &ci, &vs, &x, &y_zig);
    @memset(&y_asm, 0);
    asmTransposeMultiply(dim, @intFromPtr(&rs), @intFromPtr(&ci),
        @intFromPtr(&vs), @intFromPtr(&x), @intFromPtr(&y_asm));

    var max_diff: f64 = 0;
    inline for (&y_zig, &y_asm) |z, a| { max_diff = @max(max_diff, @abs(z - a)); }
    std.debug.print("Max diff: {d:.6}\n", .{max_diff});
    if (max_diff < 1e-12) { std.debug.print("PASS\n", .{}); } else { std.debug.print("FAIL\n", .{}); }
}
