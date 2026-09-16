//! Reproducible dense/profile comparison on a path-structured LP.
//! max sum(x), x_i+x_(i+1) <= 2, x >= 0. For even n the optimum is n.
const std = @import("std");
const z = @import("zhighs");

pub fn main(init: std.process.Init) !void {
    const a = init.gpa;
    const io = init.io;
    const stdout = std.Io.File.stdout();
    try stdout.writeStreamingAll(io, "n\tbackend\titerations\tobjective\tprimal_violation\tnormal_entries\tmedian_ms\tassembly_ms\tfactor_ms\tdirection_ms\n");
    for ([_]usize{ 128, 512 }) |n| {
        const m = n - 1;
        const starts = try a.alloc(usize, n + 1);
        defer a.free(starts);
        const rows = try a.alloc(z.foundation.RowId, 2 * m);
        defer a.free(rows);
        const values = try a.alloc(f64, 2 * m);
        defer a.free(values);
        @memset(values, 1);
        const cost = try a.alloc(f64, n);
        defer a.free(cost);
        @memset(cost, 1);
        const lower = try a.alloc(f64, n);
        defer a.free(lower);
        @memset(lower, 0);
        const upper = try a.alloc(f64, n);
        defer a.free(upper);
        @memset(upper, std.math.inf(f64));
        const rl = try a.alloc(f64, m);
        defer a.free(rl);
        @memset(rl, -std.math.inf(f64));
        const ru = try a.alloc(f64, m);
        defer a.free(ru);
        @memset(ru, 2);
        var pos: usize = 0;
        for (0..n) |j| {
            starts[j] = pos;
            if (j > 0) {
                rows[pos] = z.foundation.RowId.fromUsizeAssumeValid(j - 1);
                pos += 1;
            }
            if (j < m) {
                rows[pos] = z.foundation.RowId.fromUsizeAssumeValid(j);
                pos += 1;
            }
        }
        starts[n] = pos;
        const problem = z.lp.ProblemView{
            .num_rows = m,
            .num_cols = n,
            .col_cost = cost,
            .col_lower = lower,
            .col_upper = upper,
            .row_lower = rl,
            .row_upper = ru,
            .matrix = z.matrix.CscView.initAssumeValid(m, n, starts, rows, values),
            .objective_sense = .maximize,
            .objective_offset = 0,
        };
        for ([_]bool{ true, false }) |dense| {
            var times: [3]u64 = undefined;
            var representative: z.ipm.engine.Stats = .{};
            var objective: f64 = 0;
            var violation: f64 = 0;
            var iterations: usize = 0;
            // One warm-up, then three complete solves including transformation/allocation.
            for (0..4) |repeat| {
                var stats: z.ipm.engine.Stats = .{};
                const started = std.Io.Clock.awake.now(io).nanoseconds;
                var result = try z.ipm.engine.solveWithStats(a, problem, .{ .dense_normal = dense }, io, &stats);
                defer result.deinit();
                const elapsed: u64 = @intCast(std.Io.Clock.awake.now(io).nanoseconds - started);
                if (result.status != .optimal) return error.NonOptimal;
                objective = result.objective_value;
                if (@abs(objective - @as(f64, @floatFromInt(n))) > 1e-5) return error.WrongObjective;
                violation = 0;
                for (0..m) |i| violation = @max(violation, result.primal[i] + result.primal[i + 1] - 2);
                for (result.primal) |value| {
                    if (!std.math.isFinite(value)) return error.NonFinite;
                    violation = @max(violation, -value);
                }
                if (violation > 1e-7) return error.Infeasible;
                if (repeat > 0) times[repeat - 1] = elapsed;
                representative = stats;
                iterations = result.iterations;
            }
            std.mem.sort(u64, &times, {}, std.sort.asc(u64));
            const s = representative;
            const line = try std.fmt.allocPrint(a, "{d}\t{s}\t{d}\t{d:.8}\t{e}\t{d}\t{d:.3}\t{d:.3}\t{d:.3}\t{d:.3}\n", .{
                n,            if (dense) "dense" else "profile", iterations,      objective,          violation, s.normal_entries,
                ms(times[1]), ms(s.assembly_ns),                 ms(s.factor_ns), ms(s.direction_ns),
            });
            defer a.free(line);
            try stdout.writeStreamingAll(io, line);
        }
    }
}

fn ms(ns: u64) f64 {
    return @as(f64, @floatFromInt(ns)) / std.time.ns_per_ms;
}
