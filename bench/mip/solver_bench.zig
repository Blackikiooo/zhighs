//! Repeatable end-to-end MIP benchmark and profiling entry point.
//!
//! Parses once, then reports aggregate branch-and-bound timings and structural
//! counters as TSV. Build with ReleaseFast and run the installed binary under
//! samply, perf, VTune, or another sampling profiler.

const std = @import("std");
const zhighs = @import("zhighs");

fn nowNs(io: std.Io) i128 {
    return std.Io.Clock.now(.awake, io).nanoseconds;
}

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;
    const io_context = init.io;
    var args = try init.minimal.args.iterateAllocator(allocator);
    defer args.deinit();
    _ = args.next();
    const path = args.next() orelse return error.MissingModelPath;
    const repetitions = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 7;
    if (repetitions == 0 or args.next() != null) return error.InvalidArguments;

    var model = try zhighs.io.readFile(io_context, allocator, path, .{ .keep_names = false });
    defer model.deinit();
    const integrality = try allocator.alloc(zhighs.mip.Integrality, model.col_type.len);
    defer allocator.free(integrality);
    for (model.col_type, integrality) |source, *target| target.* = switch (source) {
        .continuous => .continuous,
        .binary => .binary,
        .integer => .integer,
        .semi_continuous => .semi_continuous,
        .semi_integer => .semi_integer,
    };
    const problem = zhighs.mip.ProblemView{
        .relaxation = .{
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
        },
        .integrality = integrality,
    };

    var lp_elapsed_ns: u128 = 0;
    for (0..repetitions) |_| {
        const started = nowNs(io_context);
        var result = zhighs.solver.solveLp(allocator, problem.relaxation, .{});
        lp_elapsed_ns += @intCast(nowNs(io_context) - started);
        if (result.status != .optimal) return error.UnexpectedRelaxationStatus;
        result.deinit();
    }

    var mip_elapsed_ns: u128 = 0;
    var last: zhighs.mip.SolveResult = undefined;
    for (0..repetitions) |run| {
        const started = nowNs(io_context);
        var result = try zhighs.mip.solve(allocator, problem, .{});
        mip_elapsed_ns += @intCast(nowNs(io_context) - started);
        if (run + 1 == repetitions) {
            last = result;
        } else {
            result.deinit();
        }
    }
    defer last.deinit();

    const line = try std.fmt.allocPrint(
        allocator,
        "zhighs-mip\t{s}\tstatus={s}\truns={d}\tlp_mean_ns={d}\tmip_mean_ns={d}\tobjective={d:.17}\tbound={d:.17}\tnodes={d}\tlp_iterations={d}\tlp_work={d}\tcold_solves={d}\treoptimizations={d}\n",
        .{ path, @tagName(last.status), repetitions, lp_elapsed_ns / repetitions, mip_elapsed_ns / repetitions, last.objective_value, last.objective_bound, last.nodes, last.lp_iterations, last.lp_work, last.lp_cold_solves, last.lp_reoptimizations },
    );
    defer allocator.free(line);
    try std.Io.File.stdout().writeStreamingAll(io_context, line);
}
