//! End-to-end MPS parser + revised-simplex acceptance runner.
//!
//! Output is a single TSV row so corpus scripts can compare status, objective,
//! iteration count, residuals, factorization lifecycle, and elapsed solve time.

const std = @import("std");
const zhighs = @import("zhighs");

fn nowNs() i128 {
    var timestamp: std.posix.timespec = undefined;
    switch (std.posix.errno(std.posix.system.clock_gettime(.MONOTONIC, &timestamp))) {
        .SUCCESS => {},
        else => unreachable,
    }
    return @as(i128, timestamp.sec) * std.time.ns_per_s + timestamp.nsec;
}

fn violation(value: f64, lower: f64, upper: f64) f64 {
    return @max(@max(lower - value, value - upper), 0.0);
}

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;
    const io_context = init.io;
    var args = std.process.Args.iterate(init.minimal.args);
    _ = args.next();
    const path = args.next() orelse return error.MissingModelPath;
    const max_iterations = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 1_000_000;
    const max_updates = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 100;
    const trace_enabled = if (args.next()) |text| std.mem.eql(u8, text, "trace") else false;
    if (args.next() != null) return error.InvalidArguments;

    const started = nowNs();
    var model = try zhighs.io.readFile(io_context, allocator, path, .{ .keep_names = false });
    defer model.deinit();
    const parsed_ns: u64 = @intCast(nowNs() - started);
    const problem = zhighs.lp.ProblemView{
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
    var engine = zhighs.lp.simplex.engine.SimplexEngine.init(allocator);
    defer engine.deinit();
    engine.numerical.max_update_count = max_updates;
    const trace: []zhighs.lp.simplex.engine.PivotTraceEvent = if (trace_enabled)
        try allocator.alloc(zhighs.lp.simplex.engine.PivotTraceEvent, @min(max_iterations, 100_000))
    else
        &.{};
    defer if (trace_enabled) allocator.free(trace);
    const solve_started = nowNs();
    const status = engine.solveProblem(problem, .{ .max_iterations = max_iterations, .pivot_trace = trace });
    const solve_ns: u64 = @intCast(nowNs() - solve_started);
    if (trace_enabled) for (trace[0..engine.pivot_trace_count]) |event| {
        const trace_line = try std.fmt.allocPrint(allocator, "pivot\t{d}\t{d}\t{d}\t{d}\t{d:.17}\t{d:.17}\t{d}\t{e:.6}\t{e:.6}\n", .{
            event.iteration,          event.entering_column, event.leaving_column, event.leaving_row,
            event.pivot,              event.step,            event.update_count,   event.ftran_relative_residual,
            event.condition_estimate,
        });
        defer allocator.free(trace_line);
        try std.Io.File.stderr().writeStreamingAll(io_context, trace_line);
    };

    var primal_residual: f64 = std.math.inf(f64);
    var dual_residual: f64 = std.math.inf(f64);
    if (engine.solutionView(problem, status)) |solution| {
        primal_residual = 0.0;
        dual_residual = 0.0;
        for (solution.primal, model.col_lower, model.col_upper) |value, lower, upper|
            primal_residual = @max(primal_residual, violation(value, lower, upper));
        const activity = try allocator.alloc(f64, problem.num_rows);
        defer allocator.free(activity);
        @memset(activity, 0.0);
        for (0..problem.num_cols) |column| {
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient|
                activity[row.toUsize()] += coefficient * solution.primal[column];
        }
        for (activity, model.row_lower, model.row_upper) |value, lower, upper|
            primal_residual = @max(primal_residual, violation(value, lower, upper));
        if (engine.basis) |basis| for (solution.reduced_cost, basis.col_status[0..problem.num_cols]) |cost, column_status| {
            const current = switch (column_status) {
                .at_lower => @max(-cost, 0.0),
                .at_upper => @max(cost, 0.0),
                .basic, .free, .superbasic => @abs(cost),
                .fixed => 0.0,
            };
            dual_residual = @max(dual_residual, current);
        };
    }

    const stats = engine.factorization.stats;
    const reinversions = stats.update_limit_reinversions + stats.update_growth_reinversions + stats.solve_residual_reinversions;
    const total_ns = parsed_ns + solve_ns;
    const line = try std.fmt.allocPrint(
        allocator,
        "zhighs\t{s}\t{s}\t{d:.17}\t{d}\t{e:.6}\t{e:.6}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\n",
        .{ path, @tagName(status), engine.objective_value, engine.iterations, primal_residual, dual_residual, stats.factorizations, reinversions, stats.update_limit_reinversions, stats.update_growth_reinversions, stats.ft_updates, engine.factorization.update_count, parsed_ns, solve_ns, total_ns },
    );
    defer allocator.free(line);
    try std.Io.File.stdout().writeStreamingAll(io_context, line);
}
