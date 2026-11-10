//! Presolve ablation benchmark: parse once, then compare cold solve pipelines
//! with no presolve, legacy rules only, each new rule alone, and every rule.
//! Output is TSV so a driver script can aggregate medians and reduction counts.

const std = @import("std");
const zhighs = @import("zhighs");

const Mode = enum {
    off,
    legacy,
    implied,
    parallel,
    full,

    fn options(self: Mode) zhighs.presolve.PresolveOptions {
        return switch (self) {
            .off => .{},
            .legacy => .{ .implied_fixing = false, .parallel_rows = false },
            .implied => .{ .implied_fixing = true, .parallel_rows = false },
            .parallel => .{ .implied_fixing = false, .parallel_rows = true },
            .full => .{ .implied_fixing = true, .parallel_rows = true },
        };
    }

    fn parse(text: []const u8) ?Mode {
        inline for (std.meta.fields(Mode)) |field| {
            if (std.mem.eql(u8, text, field.name)) return @enumFromInt(field.value);
        }
        return null;
    }
};

const Sample = struct {
    status: zhighs.lp.SolveStatus,
    objective: f64,
    iterations: usize,
    rows: usize,
    cols: usize,
    nnz: usize,
    presolve_ns: u128,
    solve_ns: u128,
    postsolve_ns: u128,
    violation: f64,
    column_substitutions: u32,
    free_singletons: u32,
    singleton_rows: u32,
    redundant_rows: u32,
    duplicate_rows: u32,
    implied_fixed: u32,
    parallel_rows: u32,
};

fn nowNs(io: std.Io) i128 {
    return std.Io.Clock.now(.awake, io).nanoseconds;
}

fn maxViolation(allocator: std.mem.Allocator, problem: *const zhighs.lp.ProblemView, primal: []const f64) f64 {
    var worst: f64 = 0.0;
    for (primal, problem.col_lower, problem.col_upper) |value, lower, upper| {
        const scale = @max(1.0, @abs(value));
        if (value < lower) worst = @max(worst, (lower - value) / scale);
        if (value > upper) worst = @max(worst, (value - upper) / scale);
    }
    const activity = allocator.alloc(f64, problem.num_rows) catch return worst;
    defer allocator.free(activity);
    @memset(activity, 0.0);
    for (0..problem.num_cols) |column| {
        const begin = problem.matrix.col_starts[column];
        const end = problem.matrix.col_starts[column + 1];
        for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row_id, coefficient|
            activity[row_id.toUsize()] += coefficient * primal[column];
    }
    for (activity, problem.row_lower, problem.row_upper) |value, lower, upper| {
        const scale = @max(1.0, @abs(value));
        if (value < lower) worst = @max(worst, (lower - value) / scale);
        if (value > upper) worst = @max(worst, (value - upper) / scale);
    }
    return worst;
}


fn runOnce(
    io: std.Io,
    allocator: std.mem.Allocator,
    problem: *const zhighs.lp.ProblemView,
    nature: ?[]const zhighs.presolve.ColumnNature,
    mode: Mode,
) !Sample {
    var engine = zhighs.lp.simplex.engine.SimplexEngine.init(allocator);
    defer engine.deinit();
    var control = zhighs.solver.LpSolveControl{};

    var sample = Sample{
        .status = .numerical_failure,
        .objective = 0.0,
        .iterations = 0,
        .rows = problem.num_rows,
        .cols = problem.num_cols,
        .nnz = problem.matrix.values.len,
        .presolve_ns = 0,
        .solve_ns = 0,
        .postsolve_ns = 0,
        .violation = 0.0,
        .column_substitutions = 0,
        .free_singletons = 0,
        .singleton_rows = 0,
        .redundant_rows = 0,
        .duplicate_rows = 0,
        .implied_fixed = 0,
        .parallel_rows = 0,
    };
    if (mode == .off) {
        const started = nowNs(io);
        sample.status = engine.solveProblem(problem, &control);
        sample.solve_ns = @intCast(nowNs(io) - started);
        const solution = engine.solutionView(problem, sample.status) orelse return error.MissingSolution;
        sample.objective = solution.objective_value;
        sample.iterations = solution.iterations;
        sample.violation = maxViolation(allocator, problem, solution.primal);
        return sample;
    }

    const options = mode.options();
    const presolve_started = nowNs(io);
    var pipeline = if (nature) |kinds|
        try zhighs.presolve.PresolvePipeline.initMipWithOptions(allocator, problem, kinds, options)
    else
        try zhighs.presolve.PresolvePipeline.initWithOptions(allocator, problem, options);
    sample.presolve_ns = @intCast(nowNs(io) - presolve_started);
    defer pipeline.deinit();

    sample.column_substitutions = pipeline.column_substitution_count;
    sample.free_singletons = pipeline.free_singleton_column_count;
    sample.singleton_rows = pipeline.singleton_row_count;
    sample.redundant_rows = pipeline.redundant_row_count;
    sample.duplicate_rows = pipeline.duplicate_row_count;
    sample.implied_fixed = pipeline.implied_fixed_column_count;
    sample.parallel_rows = pipeline.parallel_row_count;

    const reduced = if (pipeline.wasApplied()) pipeline.reducedView() else problem.*;
    sample.rows = reduced.num_rows;
    sample.cols = reduced.num_cols;
    sample.nnz = reduced.matrix.values.len;

    const solve_started = nowNs(io);
    sample.status = engine.solveProblem(&reduced, &control);
    sample.solve_ns = @intCast(nowNs(io) - solve_started);
    const solution = engine.solutionView(&reduced, sample.status) orelse return error.MissingSolution;
    sample.iterations = solution.iterations;

    if (!pipeline.wasApplied()) {
        sample.objective = solution.objective_value;
        sample.violation = maxViolation(allocator, problem, solution.primal);
        return sample;
    }

    const postsolve_started = nowNs(io);
    var recovered = try pipeline.postsolveAll(solution, allocator);
    sample.postsolve_ns = @intCast(nowNs(io) - postsolve_started);
    defer recovered.deinit();
    sample.objective = recovered.objective_value;
    sample.violation = maxViolation(allocator, problem, recovered.primal);
    return sample;
}



pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;
    const io_context = init.io;
    var args = try init.minimal.args.iterateAllocator(allocator);
    defer args.deinit();
    _ = args.next();
    const path = args.next() orelse return error.MissingModelPath;
    const mode_text = args.next() orelse return error.MissingMode;
    const repetitions = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 3;
    if (repetitions == 0 or args.next() != null) return error.InvalidArguments;
    const mode = Mode.parse(mode_text) orelse return error.InvalidMode;

    var model = try zhighs.io.readFile(io_context, allocator, path, .{ .keep_names = false });
    defer model.deinit();
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

    var nature_storage: ?[]zhighs.presolve.ColumnNature = null;
    defer if (nature_storage) |slice| allocator.free(slice);
    const nature: ?[]const zhighs.presolve.ColumnNature = blk: {
        for (model.col_type) |kind| {
            if (kind != .continuous) break;
        } else break :blk null;
        const slice = try allocator.alloc(zhighs.presolve.ColumnNature, model.col_type.len);
        for (model.col_type, slice) |kind, *target| target.* = switch (kind) {
            .continuous, .semi_continuous, .semi_integer => .continuous,
            .binary => .binary,
            .integer => .discrete,
        };
        nature_storage = slice;
        break :blk slice;
    };

    // Cold reference solve without presolve; every optimal mode must reproduce
    // it. The raw simplex may fail on a model that presolve makes trivial, so a
    // non-optimal reference is recorded and the requested mode still runs.
    const reference = try runOnce(io_context, allocator, &problem, nature, .off);
    const reference_optimal = reference.status == .optimal;

    var presolve_total: u128 = 0;
    var solve_total: u128 = 0;
    var postsolve_total: u128 = 0;
    var iterations_total: usize = 0;
    var last: Sample = undefined;
    for (0..repetitions) |_| {
        last = try runOnce(io_context, allocator, &problem, nature, mode);
        presolve_total += last.presolve_ns;
        solve_total += last.solve_ns;
        postsolve_total += last.postsolve_ns;
        iterations_total += last.iterations;
        if (last.status != .optimal) break;
        if (reference_optimal) {
            const scale = @max(1.0, @abs(reference.objective));
            if (@abs(last.objective - reference.objective) > 1e-6 * scale) {
                const line = try std.fmt.allocPrint(allocator, "presolve-bench\t{s}\tmode={s}\tstatus={s}\tobjective_mismatch=1\texpected={d:.17}\tactual={d:.17}\n", .{ path, mode_text, @tagName(last.status), reference.objective, last.objective });
                defer allocator.free(line);
                try std.Io.File.stdout().writeStreamingAll(io_context, line);
                return error.ObjectiveMismatch;
            }
        }
    }

    const reps: u128 = repetitions;
    const total_ns = presolve_total + solve_total + postsolve_total;
    const line = try std.fmt.allocPrint(
        allocator,
        "presolve-bench\t{s}\tmode={s}\tstatus={s}\treps={d}\trows0={d}\tcols0={d}\tnnz0={d}\trows={d}\tcols={d}\tnnz={d}\t" ++
            "subst={d}\tfree_singletons={d}\tsingleton_rows={d}\tredundant_rows={d}\tduplicate_rows={d}\timplied_fixed={d}\tparallel_rows={d}\t" ++
            "presolve_mean_ns={d}\tsolve_mean_ns={d}\tpostsolve_mean_ns={d}\ttotal_mean_ns={d}\titerations={d}\tobjective={d:.17}\tviolation={e:.3}\n",
        .{
            path,
            mode_text,
            @tagName(last.status),
            repetitions,
            problem.num_rows,
            problem.num_cols,
            problem.matrix.values.len,
            last.rows,
            last.cols,
            last.nnz,
            last.column_substitutions,
            last.free_singletons,
            last.singleton_rows,
            last.redundant_rows,
            last.duplicate_rows,
            last.implied_fixed,
            last.parallel_rows,
            presolve_total / reps,
            solve_total / reps,
            postsolve_total / reps,
            total_ns / reps,
            iterations_total / reps,
            last.objective,
            last.violation,
        },
    );
    defer allocator.free(line);
    try std.Io.File.stdout().writeStreamingAll(io_context, line);
}
