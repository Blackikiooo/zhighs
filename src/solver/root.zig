//! Top-level solve orchestration.
//!
//! Owns the incumbent model and coordinates validation, presolve, solver
//! selection, postsolve, callbacks, limits, and final status publication.

const std = @import("std");
const lp = @import("lp");

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

    pub fn deinit(self: *LpSolveResult) void {
        self.engine.deinit();
        self.* = undefined;
    }

    pub fn resultView(self: *const LpSolveResult, problem: LpProblemView) ?LpSolveResultView {
        const solution = self.engine.solutionView(problem, self.status) orelse return null;
        const basis = if (self.engine.basis) |*value| value else return null;
        return .{
            .status = self.status,
            .primal = solution.primal,
            .dual = solution.dual,
            .reduced_cost = solution.reduced_cost,
            .unbounded_ray = solution.unbounded_ray,
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
    const status = engine.solveProblem(problem, control);
    return .{ .engine = engine, .status = status };
}

/// Persistent LP solve lifetime. It owns all mutable numerical state but never
/// retains borrowed problem arrays between calls.
pub const LpSolveSession = struct {
    engine: lp.simplex.engine.SimplexEngine,
    last_revision: ?LpRevisionStamp = null,
    last_rows: usize = 0,
    last_cols: usize = 0,
    cold_solves: usize = 0,
    reoptimizations: usize = 0,

    pub fn init(allocator: std.mem.Allocator) LpSolveSession {
        return .{ .engine = lp.simplex.engine.SimplexEngine.init(allocator) };
    }

    pub fn deinit(self: *LpSolveSession) void {
        self.engine.deinit();
        self.* = undefined;
    }

    pub fn solve(self: *LpSolveSession, problem: LpProblemView, control: LpSolveControl, revision: LpRevisionStamp) LpSolveStatus {
        const reusable = if (self.last_revision) |last|
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
        const solution = self.engine.solutionView(problem, status) orelse return null;
        const basis = if (self.engine.basis) |*value| value else return null;
        return .{
            .status = status,
            .primal = solution.primal,
            .dual = solution.dual,
            .reduced_cost = solution.reduced_cost,
            .structural_status = basis.col_status[0..problem.num_cols],
            .logical_status = basis.col_status[problem.num_cols..][0..problem.num_rows],
            .basic_index = basis.basic_index,
            .row_scale = basis.row_scale,
            .objective_value = solution.objective_value,
            .iterations = solution.iterations,
        };
    }
};

test {
    std.testing.refAllDecls(@This());
}
