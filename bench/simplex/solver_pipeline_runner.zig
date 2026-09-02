//! MPS acceptance runner for the complete solver -> presolve -> simplex ->
//! postsolve pipeline. Output is one stable TSV row per invocation.

const std = @import("std");
const zhighs = @import("zhighs");

fn nowNs(io_context: std.Io) i96 {
    return std.Io.Clock.awake.now(io_context).nanoseconds;
}

fn violation(value: f64, lower: f64, upper: f64) f64 {
    return @max(@max(lower - value, value - upper), 0.0);
}

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;
    const io_context = init.io;
    var args = try init.minimal.args.iterateAllocator(allocator);
    defer args.deinit();
    _ = args.next();
    const path = args.next() orelse return error.MissingModelPath;
    const max_iterations = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 1_000_000;
    const phase_one_strategy: zhighs.lp.simplex.engine.PhaseOneStrategy = if (args.next()) |text|
        if (std.mem.eql(u8, text, "dual")) .dual else if (std.mem.eql(u8, text, "auto")) .automatic else if (std.mem.eql(u8, text, "primal")) .primal else return error.InvalidArguments
    else
        .primal;
    if (args.next() != null) return error.InvalidArguments;

    var model = try zhighs.io.readFile(io_context, allocator, path, .{ .keep_names = false });
    defer model.deinit();
    const problem = zhighs.solver.LpProblemView{
        .num_rows = model.row_lower.len,
        .num_cols = model.col_cost.len,
        .col_cost = model.col_cost,
        .col_lower = model.col_lower,
        .col_upper = model.col_upper,
        .row_lower = model.row_lower,
        .row_upper = model.row_upper,
        .matrix = model.matrix.view(),
        .objective_sense = if (model.objective_sense == .minimize) .minimize else .maximize,
        .objective_offset = model.objective_offset,
    };

    const started = nowNs(io_context);
    var result = zhighs.solver.solveLp(allocator, problem, .{
        .max_iterations = max_iterations,
        .phase_one_strategy = phase_one_strategy,
    });
    defer result.deinit();
    const solve_ns: u64 = @intCast(@max(nowNs(io_context) - started, 0));
    const view = result.resultView(problem) orelse {
        const line = try std.fmt.allocPrint(allocator, "zhighs-solver\t{s}\t{s}\tunavailable\t{d}\n", .{ path, @tagName(result.status), solve_ns });
        defer allocator.free(line);
        try std.Io.File.stdout().writeStreamingAll(io_context, line);
        return;
    };

    var primal_residual = std.math.inf(f64);
    var dual_residual = std.math.inf(f64);
    var ray_residual = std.math.inf(f64);
    var ray_objective: f64 = 0.0;
    const activity = try allocator.alloc(f64, problem.num_rows);
    defer allocator.free(activity);

    if (view.primal_valid) {
        primal_residual = 0.0;
        for (view.primal, problem.col_lower, problem.col_upper) |value, lower, upper|
            primal_residual = @max(primal_residual, violation(value, lower, upper));
        @memset(activity, 0.0);
        for (0..problem.num_cols) |column| {
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient|
                activity[row.toUsize()] += coefficient * view.primal[column];
        }
        for (activity, problem.row_lower, problem.row_upper) |value, lower, upper|
            primal_residual = @max(primal_residual, violation(value, lower, upper));
    }
    if (view.dual_valid and view.basis_valid) {
        dual_residual = 0.0;
        for (view.reduced_cost, view.structural_status) |cost, status| {
            const current = switch (status) {
                .at_lower => @max(-cost, 0.0),
                .at_upper => @max(cost, 0.0),
                .basic, .free, .superbasic => @abs(cost),
                .fixed => 0.0,
            };
            dual_residual = @max(dual_residual, current);
        }
    }
    if (view.unbounded_ray.len == problem.num_cols) {
        ray_residual = 0.0;
        @memset(activity, 0.0);
        for (0..problem.num_cols) |column| {
            const direction = view.unbounded_ray[column];
            ray_objective += problem.col_cost[column] * direction;
            if (std.math.isFinite(problem.col_upper[column])) ray_residual = @max(ray_residual, direction);
            if (std.math.isFinite(problem.col_lower[column])) ray_residual = @max(ray_residual, -direction);
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient|
                activity[row.toUsize()] += coefficient * direction;
        }
        for (activity, problem.row_lower, problem.row_upper) |direction, lower, upper| {
            if (std.math.isFinite(upper)) ray_residual = @max(ray_residual, direction);
            if (std.math.isFinite(lower)) ray_residual = @max(ray_residual, -direction);
        }
        ray_residual = @max(ray_residual, 0.0);
    }

    const line = try std.fmt.allocPrint(
        allocator,
        "zhighs-solver\t{s}\t{s}\t{d:.17}\t{d}\t{e:.6}\t{e:.6}\t{e:.6}\t{e:.6}\t{}\t{}\t{}\t{}\t{d}\n",
        .{ path, @tagName(view.status), view.objective_value, view.iterations, primal_residual, dual_residual, ray_residual, ray_objective, view.primal_valid, view.dual_valid, view.basis_valid, result.presolved_problem != null, solve_ns },
    );
    defer allocator.free(line);
    try std.Io.File.stdout().writeStreamingAll(io_context, line);
}
