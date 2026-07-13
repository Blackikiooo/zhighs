//! Cold-vs-warm LP session microbenchmark.
//!
//! Measures first solve, objective reoptimization, RHS reoptimization, and
//! bound reoptimization through the public Model lifecycle. Build/run with
//! `zig build bench-simplex -Doptimize=ReleaseFast`.

const std = @import("std");
const zhighs = @import("zhighs");

pub fn main() !void {
    const allocator = std.heap.page_allocator;
    const io = std.Io.Threaded.global_single_threaded.io();
    var env = try zhighs.model.Env.initSimple(allocator);
    defer env.deinit();
    var model = try zhighs.model.Model.init(allocator, &env, "session_bench");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .greater_equal, 1.0, null);
    try model.addVar(1, &[_]usize{0}, &[_]f64{1.0}, 2.0, 0.0, std.math.inf(f64), .continuous, null);
    try model.addVar(1, &[_]usize{0}, &[_]f64{1.0}, 1.0, 0.0, 1.0, .continuous, null);

    var lap_start = nowNs(io);
    try model.optimize();
    const cold_ns = lapNs(io, &lap_start);

    const repetitions: usize = 10_000;
    for (0..repetitions) |iteration| {
        try model.setDblAttrElement(.rhs, 0, if (iteration & 1 == 0) 0.5 else 1.0);
        try model.optimize();
    }
    const rhs_ns = lapNs(io, &lap_start);
    for (0..repetitions) |iteration| {
        try model.setDblAttrElement(.obj, 0, if (iteration & 1 == 0) 1.5 else 2.0);
        try model.optimize();
    }
    const objective_ns = lapNs(io, &lap_start);
    for (0..repetitions) |iteration| {
        try model.setDblAttrElement(.ub, 1, if (iteration & 1 == 0) 0.75 else 1.0);
        try model.optimize();
    }
    const bounds_ns = lapNs(io, &lap_start);

    std.debug.print(
        "cold={d}ns rhs={d}ns/op objective={d}ns/op bounds={d}ns/op\n",
        .{ cold_ns, rhs_ns / repetitions, objective_ns / repetitions, bounds_ns / repetitions },
    );

    try benchmarkKernels(allocator, io);
}

fn benchmarkKernels(allocator: std.mem.Allocator, io: std.Io) !void {
    const dimension: usize = 64;
    const repetitions: usize = 20_000;
    var factorization = zhighs.lp.simplex.factorization.Factorization.init(allocator);
    defer factorization.deinit();
    try factorization.factorizeIdentity(dimension);
    const rhs = try allocator.alloc(f64, dimension);
    defer allocator.free(rhs);
    @memset(rhs, 1.0);

    var lap_start = nowNs(io);
    for (0..repetitions) |_| {
        try factorization.solve(rhs);
        std.mem.doNotOptimizeAway(rhs.ptr);
    }
    const ftran_ns = lapNs(io, &lap_start);
    for (0..repetitions) |_| {
        try factorization.solveTranspose(rhs);
        std.mem.doNotOptimizeAway(rhs.ptr);
    }
    const btran_ns = lapNs(io, &lap_start);

    const Status = zhighs.lp.simplex.basis.BasisStatus;
    var tableau: [dimension]f64 = undefined;
    var reduced: [dimension]f64 = undefined;
    var status: [dimension]Status = undefined;
    var lower: [dimension]f64 = undefined;
    var upper: [dimension]f64 = undefined;
    var primal: [dimension]f64 = undefined;
    var ratios: [dimension]f64 = undefined;
    var directions: [dimension]f64 = undefined;
    var candidates: [dimension]u32 = undefined;
    for (0..dimension) |i| {
        tableau[i] = -1.0 - @as(f64, @floatFromInt(i % 7));
        reduced[i] = @as(f64, @floatFromInt(i));
        status[i] = .at_lower;
        lower[i] = 0.0;
        upper[i] = if (i % 3 == 0) 1.0 else std.math.inf(f64);
        primal[i] = 0.0;
    }
    const ratio = zhighs.lp.simplex.ratio_test.RatioTest{};
    for (0..repetitions) |_| {
        const choice = ratio.chooseDualEntering(
            &tableau,
            &reduced,
            &status,
            &lower,
            &upper,
            &primal,
            .at_lower,
            16.0,
            &ratios,
            &directions,
            &candidates,
        );
        std.mem.doNotOptimizeAway(choice);
    }
    const bfrt_ns = lapNs(io, &lap_start);
    std.debug.print(
        "ftran64={d}ns/op btran64={d}ns/op bfrt64={d}ns/op\n",
        .{ ftran_ns / repetitions, btran_ns / repetitions, bfrt_ns / repetitions },
    );
}

fn nowNs(io: std.Io) i96 {
    return std.Io.Clock.awake.now(io).nanoseconds;
}

fn lapNs(io: std.Io, start: *i96) u64 {
    const current = nowNs(io);
    const elapsed: u64 = @intCast(@max(current - start.*, 0));
    start.* = current;
    return elapsed;
}
