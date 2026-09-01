//! Top-level solve orchestration.
//!
//! Owns the incumbent model and coordinates validation, presolve, solver
//! selection, postsolve, callbacks, limits, and final status publication.

const std = @import("std");
const foundation = @import("foundation");
const lp = @import("lp");
const matrix = @import("matrix");
const presolve = @import("presolve");

pub const LpSolveControl = lp.simplex.engine.SolveControl;
pub const LpSolveStatus = lp.simplex.solution.SolveStatus;
pub const LpProblemView = lp.ProblemView;
pub const LpObjectiveSense = lp.ObjectiveSense;
pub const LpBasisStatus = lp.simplex.basis.BasisStatus;
pub const LpBasisView = lp.BasisView;
pub const LpBasisSnapshot = lp.BasisSnapshot;
pub const LpProgressEventView = lp.simplex.engine.ProgressEventView;
pub const LpCallbackAction = lp.simplex.engine.CallbackAction;
pub const LpPricingRule = lp.simplex.pricing.PricingRule;

pub const LpRevisionStamp = struct {
    structure: u64,
    matrix_values: u64,
    bounds: u64,
    objective: u64,
};

/// Borrowed LP result projected from an owning `LpSolveResult`.
///
/// Every slice is invalidated by `LpSolveResult.deinit`.
pub const LpSolveResultView = struct {
    status: LpSolveStatus,
    primal: []const f64,
    dual: []const f64,
    reduced_cost: []const f64,
    unbounded_ray: []const f64,
    infeasibility_ray: []const f64,
    structural_status: []const LpBasisStatus,
    logical_status: []const LpBasisStatus,
    basic_index: []const u32,
    row_scale: []const f64,
    objective_value: f64,
    iterations: usize,
};

/// Owning lifetime boundary around one simplex invocation.
pub const LpSolveResult = struct {
    engine: lp.simplex.engine.SimplexEngine,
    status: LpSolveStatus,
    presolved_problem: ?presolve.PresolvedProblem = null,
    projection: ?PostsolveProjection = null,

    pub fn deinit(self: *LpSolveResult) void {
        if (self.projection) |*value| value.deinit();
        if (self.presolved_problem) |*value| value.deinit();
        self.engine.deinit();
        self.* = undefined;
    }

    pub fn resultView(self: *const LpSolveResult, problem: LpProblemView) ?LpSolveResultView {
        if (self.projection) |*projection| return projection.view(self.engine.iterations);
        if (self.presolved_problem != null) return null;
        const solution = self.engine.solutionView(problem, self.status) orelse return null;
        const basis = if (self.engine.basis) |*value| value else return null;
        return .{
            .status = self.status,
            .primal = solution.primal,
            .dual = solution.dual,
            .reduced_cost = solution.reduced_cost,
            .unbounded_ray = solution.unbounded_ray,
            .infeasibility_ray = solution.infeasibility_ray,
            .structural_status = basis.col_status[0..problem.num_cols],
            .logical_status = basis.col_status[problem.num_cols..][0..problem.num_rows],
            .basic_index = basis.basic_index,
            .row_scale = basis.row_scale,
            .objective_value = solution.objective_value,
            .iterations = solution.iterations,
        };
    }
};

/// Execute the LP backend while borrowing all immutable problem arrays.
pub fn solveLp(allocator: std.mem.Allocator, problem: LpProblemView, control: LpSolveControl) LpSolveResult {
    var engine = lp.simplex.engine.SimplexEngine.init(allocator);
    if (!hasFixedColumns(problem)) {
        const status = engine.solveProblem(problem, control);
        return .{ .engine = engine, .status = status };
    }

    var pp = presolve.presolve(allocator, problem) catch
        return .{ .engine = engine, .status = .numerical_failure };
    var reduced_control = control;
    // An original-coordinate basis cannot be consumed by the reduced model.
    // A mapped warm-start is added at the solver boundary in P6.3.
    reduced_control.initial_basis = null;
    const reduced_problem = pp.reducedView();
    const status = engine.solveProblem(reduced_problem, reduced_control);
    const projection = PostsolveProjection.init(allocator, problem, &pp, &engine, status) catch null;
    return .{
        .engine = engine,
        .status = status,
        .presolved_problem = pp,
        .projection = projection,
    };
}

/// Persistent LP solve lifetime. It owns all mutable numerical state but never
/// retains borrowed problem arrays between calls.
pub const LpSolveSession = struct {
    allocator: std.mem.Allocator,
    engine: lp.simplex.engine.SimplexEngine,
    presolved_problem: ?presolve.PresolvedProblem = null,
    projection: ?PostsolveProjection = null,
    last_revision: ?LpRevisionStamp = null,
    last_rows: usize = 0,
    last_cols: usize = 0,
    cold_solves: usize = 0,
    reoptimizations: usize = 0,

    pub fn init(allocator: std.mem.Allocator) LpSolveSession {
        return .{ .allocator = allocator, .engine = lp.simplex.engine.SimplexEngine.init(allocator) };
    }

    pub fn deinit(self: *LpSolveSession) void {
        self.clearPresolveState();
        self.engine.deinit();
        self.* = undefined;
    }

    pub fn solve(self: *LpSolveSession, problem: LpProblemView, control: LpSolveControl, revision: LpRevisionStamp) LpSolveStatus {
        if (self.projection) |*value| value.deinit();
        self.projection = null;

        if (!hasFixedColumns(problem)) {
            if (self.presolved_problem) |*value| value.deinit();
            self.presolved_problem = null;
            return self.solveReduced(problem, control, revision, false);
        }

        const next_pp = presolve.presolve(self.allocator, problem) catch {
            self.last_revision = revision;
            return .numerical_failure;
        };
        const mapping_changed = if (self.presolved_problem) |*previous|
            !std.mem.eql(u32, previous.reduced_to_original_col, next_pp.reduced_to_original_col)
        else
            true;
        if (self.presolved_problem) |*value| value.deinit();
        self.presolved_problem = next_pp;

        var reduced_control = control;
        reduced_control.initial_basis = null;
        const reduced_problem = self.presolved_problem.?.reducedView();
        const status = self.solveReduced(reduced_problem, reduced_control, revision, mapping_changed);
        self.projection = PostsolveProjection.init(
            self.allocator,
            problem,
            &self.presolved_problem.?,
            &self.engine,
            status,
        ) catch {
            return .numerical_failure;
        };
        return status;
    }

    fn solveReduced(self: *LpSolveSession, problem: LpProblemView, control: LpSolveControl, revision: LpRevisionStamp, force_cold: bool) LpSolveStatus {
        const reusable = !force_cold and if (self.last_revision) |last|
            last.structure == revision.structure and last.matrix_values == revision.matrix_values and
                self.last_rows == problem.num_rows and self.last_cols == problem.num_cols
        else
            false;
        const status = if (reusable) blk: {
            self.reoptimizations += 1;
            break :blk self.engine.reoptimizeProblem(problem, control);
        } else blk: {
            self.cold_solves += 1;
            break :blk self.engine.solveProblem(problem, control);
        };
        self.last_revision = revision;
        self.last_rows = problem.num_rows;
        self.last_cols = problem.num_cols;
        return status;
    }

    pub fn resultView(self: *const LpSolveSession, problem: LpProblemView, status: LpSolveStatus) ?LpSolveResultView {
        if (self.projection) |*projection| return projection.view(self.engine.iterations);
        if (self.presolved_problem != null) return null;
        const solution = self.engine.solutionView(problem, status) orelse return null;
        const basis = if (self.engine.basis) |*value| value else return null;
        return .{
            .status = status,
            .primal = solution.primal,
            .dual = solution.dual,
            .reduced_cost = solution.reduced_cost,
            .unbounded_ray = solution.unbounded_ray,
            .infeasibility_ray = solution.infeasibility_ray,
            .structural_status = basis.col_status[0..problem.num_cols],
            .logical_status = basis.col_status[problem.num_cols..][0..problem.num_rows],
            .basic_index = basis.basic_index,
            .row_scale = basis.row_scale,
            .objective_value = solution.objective_value,
            .iterations = solution.iterations,
        };
    }

    fn clearPresolveState(self: *LpSolveSession) void {
        if (self.projection) |*value| value.deinit();
        if (self.presolved_problem) |*value| value.deinit();
        self.projection = null;
        self.presolved_problem = null;
    }
};

fn hasFixedColumns(problem: LpProblemView) bool {
    for (problem.col_lower, problem.col_upper) |lower, upper| {
        if (lower == upper and std.math.isFinite(lower)) return true;
    }
    return false;
}

/// Owning original-coordinate result produced after solving a reduced LP.
const PostsolveProjection = struct {
    allocator: std.mem.Allocator,
    solution: presolve.PresolveResult,
    structural_status: []LpBasisStatus,
    logical_status: []LpBasisStatus,
    basic_index: []u32,
    row_scale: []f64,

    fn init(
        allocator: std.mem.Allocator,
        original_problem: LpProblemView,
        pp: *const presolve.PresolvedProblem,
        engine: *const lp.simplex.engine.SimplexEngine,
        status: LpSolveStatus,
    ) !PostsolveProjection {
        const reduced_problem = pp.reducedView();
        const reduced_solution = engine.solutionView(reduced_problem, status) orelse return error.SolutionUnavailable;
        var result = PostsolveProjection{
            .allocator = allocator,
            .solution = try presolve.postsolve(pp, reduced_solution, allocator),
            .structural_status = &.{},
            .logical_status = &.{},
            .basic_index = &.{},
            .row_scale = &.{},
        };
        errdefer result.deinit();
        result.structural_status = try allocator.alloc(LpBasisStatus, original_problem.num_cols);
        result.logical_status = try allocator.alloc(LpBasisStatus, original_problem.num_rows);
        result.basic_index = try allocator.alloc(u32, original_problem.num_rows);
        result.row_scale = try allocator.alloc(f64, original_problem.num_rows);

        const basis = if (engine.basis) |*value| value else return error.BasisUnavailable;
        @memset(result.structural_status, .fixed);
        for (pp.reduced_to_original_col, 0..) |original_col, reduced_col| {
            result.structural_status[original_col] = basis.col_status[reduced_col];
        }
        @memcpy(result.logical_status, basis.col_status[pp.num_cols..][0..pp.num_rows]);
        @memcpy(result.row_scale, basis.row_scale[0..pp.num_rows]);
        for (basis.basic_index[0..pp.num_rows], 0..) |reduced_col, row| {
            result.basic_index[row] = if (reduced_col < pp.num_cols)
                pp.reduced_to_original_col[reduced_col]
            else if (reduced_col < pp.num_cols + pp.num_rows)
                original_problem.num_cols + (reduced_col - pp.num_cols)
            else
                std.math.maxInt(u32);
        }
        return result;
    }

    fn deinit(self: *PostsolveProjection) void {
        self.solution.deinit();
        self.allocator.free(self.structural_status);
        self.allocator.free(self.logical_status);
        self.allocator.free(self.basic_index);
        self.allocator.free(self.row_scale);
        self.* = undefined;
    }

    fn view(self: *const PostsolveProjection, iterations: usize) LpSolveResultView {
        return .{
            .status = self.solution.status,
            .primal = self.solution.primal,
            .dual = self.solution.dual,
            .reduced_cost = self.solution.reduced_cost,
            .unbounded_ray = self.solution.unbounded_ray,
            .infeasibility_ray = self.solution.infeasibility_ray,
            .structural_status = self.structural_status,
            .logical_status = self.logical_status,
            .basic_index = self.basic_index,
            .row_scale = self.row_scale,
            .objective_value = self.solution.objective_value,
            .iterations = iterations,
        };
    }
};

fn fixedColumnProblem(fixed_first: bool) LpProblemView {
    return .{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 1.0, 2.0 },
        .col_lower = if (fixed_first) &.{ 2.0, 0.0 } else &.{ 0.0, 3.0 },
        .col_upper = if (fixed_first) &.{ 2.0, std.math.inf(f64) } else &.{ std.math.inf(f64), 3.0 },
        .row_lower = &.{5.0},
        .row_upper = &.{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(
            1,
            2,
            &[_]usize{ 0, 1, 2 },
            &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) },
            &.{ 1.0, 1.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
}

test "solveLp returns fixed-column solutions and basis in original coordinates" {
    const problem = fixedColumnProblem(false);
    var result = solveLp(std.testing.allocator, problem, .{});
    defer result.deinit();
    try std.testing.expectEqual(LpSolveStatus.optimal, result.status);
    const view = result.resultView(problem).?;
    try std.testing.expectEqualSlices(f64, &.{ 2.0, 3.0 }, view.primal);
    try std.testing.expectApproxEqAbs(@as(f64, 8.0), view.objective_value, 1e-9);
    try std.testing.expectEqual(LpBasisStatus.fixed, view.structural_status[1]);
    try std.testing.expectEqual(@as(usize, 2), view.structural_status.len);
    try std.testing.expectEqual(@as(usize, 1), view.logical_status.len);
    try std.testing.expect(view.basic_index[0] < problem.num_cols + problem.num_rows);
}

test "session invalidates reuse when the fixed-column mapping changes" {
    var session = LpSolveSession.init(std.testing.allocator);
    defer session.deinit();
    const revision = LpRevisionStamp{ .structure = 1, .matrix_values = 1, .bounds = 1, .objective = 1 };

    const first = fixedColumnProblem(false);
    try std.testing.expectEqual(LpSolveStatus.optimal, session.solve(first, .{}, revision));
    try std.testing.expectEqual(@as(usize, 1), session.cold_solves);
    try std.testing.expectEqualSlices(f64, &.{ 2.0, 3.0 }, session.resultView(first, .optimal).?.primal);

    const second = fixedColumnProblem(true);
    const changed_revision = LpRevisionStamp{ .structure = 1, .matrix_values = 1, .bounds = 2, .objective = 1 };
    try std.testing.expectEqual(LpSolveStatus.optimal, session.solve(second, .{}, changed_revision));
    try std.testing.expectEqual(@as(usize, 2), session.cold_solves);
    try std.testing.expectEqualSlices(f64, &.{ 2.0, 3.0 }, session.resultView(second, .optimal).?.primal);
}

test {
    std.testing.refAllDecls(@This());
}
