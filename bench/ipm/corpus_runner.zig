//! One real-model solve per process. JSON preserves full primal/dual vectors
//! for independent validation against a separately parsed reference model.
const std = @import("std");
const z = @import("zhighs");

pub fn main(init: std.process.Init) !void {
    const a = init.gpa;
    var args = try init.minimal.args.iterateAllocator(a);
    defer args.deinit();
    _ = args.next();
    const path = args.next() orelse return error.MissingModel;
    const mode = args.next() orelse "sparse";
    const refinement_steps = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 0;
    const max_iterations = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 200;
    const feasibility_tolerance = if (args.next()) |text| try std.fmt.parseFloat(f64, text) else 1e-7;
    const primal_correction_steps = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 5;
    const max_regularization = if (args.next()) |text| try std.fmt.parseFloat(f64, text) else 1e-4;
    const equilibrate = if (args.next()) |text| std.mem.eql(u8, text, "scaled") else false;
    const ordering_text = args.next() orelse "automatic";
    const row_ordering: z.ipm.options.RowOrdering = if (std.mem.eql(u8, ordering_text, "automatic"))
        .automatic
    else if (std.mem.eql(u8, ordering_text, "natural"))
        .natural
    else if (std.mem.eql(u8, ordering_text, "rcm"))
        .reverse_cuthill_mckee
    else
        .minimum_degree;
    const started = std.Io.Clock.awake.now(init.io).nanoseconds;
    var input = try z.io.readFile(init.io, a, path, .{ .keep_names = false });
    defer input.deinit();
    const parsed = std.Io.Clock.awake.now(init.io).nanoseconds;
    const problem = z.lp.ProblemView{
        .num_rows = input.row_lower.len,
        .num_cols = input.col_cost.len,
        .col_cost = input.col_cost,
        .col_lower = input.col_lower,
        .col_upper = input.col_upper,
        .row_lower = input.row_lower,
        .row_upper = input.row_upper,
        .matrix = input.matrix.view(),
        .objective_sense = if (input.objective_sense == .minimize) .minimize else .maximize,
        .objective_offset = input.objective_offset,
    };
    var stats: z.ipm.Stats = .{};
    var result = try z.ipm.solveWithStats(a, problem, .{
        .dense_normal = std.mem.eql(u8, mode, "dense"),
        .equilibrate = equilibrate,
        .row_ordering = row_ordering,
        .max_iterations = max_iterations,
        .refinement_steps = refinement_steps,
        .primal_correction_steps = primal_correction_steps,
        .tolerances = .{
            .primal_feasibility = feasibility_tolerance,
            .dual_feasibility = feasibility_tolerance,
            .relative_gap = feasibility_tolerance,
            .max_regularization = max_regularization,
        },
        .time_limit_ns = 5 * std.time.ns_per_s,
    }, init.io, &stats);
    defer result.deinit();
    const finished = std.Io.Clock.awake.now(init.io).nanoseconds;
    const output = try std.json.Stringify.valueAlloc(a, .{
        .status = @tagName(result.status),
        .objective = if (result.status == .optimal) @as(?f64, result.objective_value) else null,
        .iterations = result.iterations,
        .rows = problem.num_rows,
        .cols = problem.num_cols,
        .nnz = problem.matrix.values.len,
        .parse_ns = parsed - started,
        .solve_ns = finished - parsed,
        .stats = stats,
        .primal = result.primal,
        .dual = result.dual,
        .reduced_cost = result.reduced_cost,
    }, .{});
    defer a.free(output);
    try std.Io.File.stdout().writeStreamingAll(init.io, output);
}
