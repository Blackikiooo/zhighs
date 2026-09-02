//! Cold/warm acceptance runner for the persistent full LP solver session.

const std = @import("std");
const zhighs = @import("zhighs");

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;
    const io_context = init.io;
    var args = try init.minimal.args.iterateAllocator(allocator);
    defer args.deinit();
    _ = args.next();
    const path = args.next() orelse return error.MissingModelPath;
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
    const revision = zhighs.solver.LpRevisionStamp{ .structure = 1, .matrix_values = 1, .bounds = 1, .objective = 1 };
    const control = zhighs.solver.LpSolveControl{ .phase_one_strategy = .dual };
    var session = zhighs.solver.LpSolveSession.init(allocator);
    defer session.deinit();

    const cold_status = session.solve(problem, control, revision);
    const cold_view = session.resultView(problem, cold_status) orelse return error.ResultUnavailable;
    const cold_objective = cold_view.objective_value;
    const cold_iterations = cold_view.iterations;
    const presolved = session.presolved_problem != null;

    const warm_status = session.solve(problem, control, revision);
    const warm_view = session.resultView(problem, warm_status) orelse return error.ResultUnavailable;
    const line = try std.fmt.allocPrint(
        allocator,
        "zhighs-session\t{s}\t{s}\t{s}\t{d:.17}\t{d:.17}\t{d}\t{d}\t{d}\t{d}\t{}\n",
        .{ path, @tagName(cold_status), @tagName(warm_status), cold_objective, warm_view.objective_value, cold_iterations, warm_view.iterations, session.cold_solves, session.reoptimizations, presolved },
    );
    defer allocator.free(line);
    try std.Io.File.stdout().writeStreamingAll(io_context, line);
}
