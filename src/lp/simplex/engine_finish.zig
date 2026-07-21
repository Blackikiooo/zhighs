//! Terminal-state handling and solution validation for `SimplexEngine`.
//!
//! ## Responsibility
//!
//! Owns optimal/unbounded finalization, unbounded-ray and optimal-solution
//! validation, and the borrowed `SolutionView` publication.

const std = @import("std");
const basis_module = @import("basis.zig");
const basis_snapshot_module = @import("basis_snapshot.zig");
const factorization_module = @import("factorization.zig");
const pricing_module = @import("pricing.zig");
const ratio_module = @import("ratio_test.zig");
const numerical_module = @import("numerical.zig");
const dual_phase_one_module = @import("dual_phase_one.zig");
const crash_module = @import("crash.zig");
const degeneracy_module = @import("degeneracy.zig");
const pricing_workspace_module = @import("pricing_workspace.zig");
const problem_module = @import("problem.zig");
const solution_module = @import("solution.zig");
const foundation = @import("foundation");
const matrix = @import("matrix");
const SimplexEngine = @import("engine.zig").SimplexEngine;
const SolveStatus = @import("engine.zig").SolveStatus;

pub fn finishOptimal(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
    if (self.recomputeBasicValues(problem) != .optimal) return .numerical_failure;
    if (self.recomputeReducedCosts(problem) != .optimal) return .numerical_failure;
    if (self.validateOptimalSolution(problem) != .optimal) return .numerical_failure;
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    self.objective_value = problem.objective_offset;
    for (problem.col_cost, basis.primal[0..problem.num_cols], basis.column_scale[0..problem.num_cols]) |cost, value, scale|
        self.objective_value += cost * value * scale;
    if (!std.math.isFinite(self.objective_value)) return .numerical_failure;
    return .optimal;
}

/// Construct and validate an original-coordinate primal ray from the
/// signed entering direction and current FTRAN column. No finite iterate
/// is published as evidence: unboundedness is accepted only when variable
/// bounds, row recession directions, and objective improvement all hold.
pub fn finishUnbounded(self: *SimplexEngine, problem: problem_module.ProblemView, entering_col: usize, entering_direction: f64) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    @memset(basis.unbounded_ray, 0.0);
    if (entering_col < problem.num_cols)
        basis.unbounded_ray[entering_col] = entering_direction * basis.column_scale[entering_col];
    for (basis.basic_index, basis.pivot_direction) |column_u32, direction| {
        const column: usize = @intCast(column_u32);
        if (column < problem.num_cols)
            basis.unbounded_ray[column] = -direction * basis.column_scale[column];
    }
    if (!self.validateUnboundedRay(problem, basis.unbounded_ray)) return .numerical_failure;
    self.unbounded_ray_valid = true;
    return .unbounded;
}

pub fn validateUnboundedRay(self: *SimplexEngine, problem: problem_module.ProblemView, ray: []const f64) bool {
    const basis = if (self.basis) |*value| value else return false;
    if (ray.len != problem.num_cols) return false;
    var ray_max: f64 = 0.0;
    for (ray, problem.col_lower, problem.col_upper) |direction, lower, upper| {
        if (!std.math.isFinite(direction)) return false;
        ray_max = @max(ray_max, @abs(direction));
        if (direction > self.numerical.zero_tolerance and std.math.isFinite(upper)) return false;
        if (direction < -self.numerical.zero_tolerance and std.math.isFinite(lower)) return false;
    }
    if (ray_max <= self.numerical.zero_tolerance) return false;

    @memset(basis.rhs_work, 0.0);
    @memset(basis.residual_work, 0.0);
    for (0..problem.num_cols) |column| {
        const direction = ray[column];
        if (direction == 0.0) continue;
        const begin = problem.matrix.col_starts[column];
        const end = problem.matrix.col_starts[column + 1];
        for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
            if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
            const row_index = row.toUsize();
            const product = coefficient * direction;
            basis.rhs_work[row_index] += product;
            basis.residual_work[row_index] += @abs(product);
        }
    }
    for (basis.rhs_work, basis.residual_work, problem.row_lower, problem.row_upper) |direction, magnitude, lower, upper| {
        const tolerance = self.numerical.residual_tolerance * @max(1.0, magnitude);
        if (std.math.isFinite(upper) and direction > tolerance) return false;
        if (std.math.isFinite(lower) and direction < -tolerance) return false;
    }

    var objective_direction: f64 = 0.0;
    var objective_magnitude: f64 = 0.0;
    for (problem.col_cost, ray) |cost, direction| {
        objective_direction += cost * direction;
        objective_magnitude += @abs(cost * direction);
    }
    const objective_tolerance = self.numerical.dual_tolerance * @max(1.0, objective_magnitude);
    return if (problem.objective_sense == .minimize)
        objective_direction < -objective_tolerance
    else
        objective_direction > objective_tolerance;
}

/// Allocation-free KKT feasibility check used before publishing an
/// optimal result. This catches accumulated update drift and incorrect
/// bound/status transitions at the solver boundary.
pub fn validateOptimalSolution(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    const primal_tolerance = self.numerical.primal_tolerance;
    const dual_tolerance = self.numerical.dual_tolerance;

    for (basis.primal[0..problem.num_cols], basis.column_scale[0..problem.num_cols], problem.col_lower, problem.col_upper) |internal, scale, lower, upper| {
        const value = internal * scale;
        if (!std.math.isFinite(value) or value < lower - primal_tolerance or value > upper + primal_tolerance)
            return .numerical_failure;
    }

    @memset(basis.rhs_work, 0.0);
    for (0..problem.num_cols) |column| {
        const value = basis.primal[column] * basis.column_scale[column];
        if (value == 0.0) continue;
        const begin = problem.matrix.col_starts[column];
        const end = problem.matrix.col_starts[column + 1];
        for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
            if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
            basis.rhs_work[row.toUsize()] += coefficient * value;
        }
    }
    for (basis.rhs_work, problem.row_lower, problem.row_upper) |activity, lower, upper| {
        if (!std.math.isFinite(activity) or activity < lower - primal_tolerance or activity > upper + primal_tolerance)
            return .numerical_failure;
    }

    for (basis.reduced_cost[0 .. problem.num_cols + problem.num_rows], basis.col_status[0 .. problem.num_cols + problem.num_rows], 0..) |internal_reduced, status, column| {
        const reduced = if (column < problem.num_cols)
            internal_reduced / (basis.column_scale[column] * self.objective_scale)
        else
            internal_reduced / self.objective_scale;
        if (!std.math.isFinite(reduced)) return .numerical_failure;
        const infeasible = switch (status) {
            .at_lower => reduced < -dual_tolerance,
            .at_upper => reduced > dual_tolerance,
            .free, .superbasic, .basic => @abs(reduced) > dual_tolerance,
            .fixed => false,
        };
        if (infeasible) return .numerical_failure;
    }

    const artificial_begin = problem.num_cols + problem.num_rows;
    for (basis.primal[artificial_begin..]) |value| {
        if (!std.math.isFinite(value) or @abs(value) > primal_tolerance) return .numerical_failure;
    }
    return .optimal;
}

/// Borrow the current engine-owned solution arrays without copying.
pub fn solutionView(self: *const SimplexEngine, problem: problem_module.ProblemView, status: SolveStatus) ?solution_module.SolutionView {
    const basis = if (self.basis) |*value| value else return null;
    for (basis.published_primal, basis.primal[0..problem.num_cols], basis.column_scale[0..problem.num_cols]) |*published, internal, scale|
        published.* = internal * scale;
    for (basis.published_dual, basis.dual) |*published, internal|
        published.* = internal / self.objective_scale;
    for (basis.published_reduced_cost, basis.reduced_cost[0..problem.num_cols], basis.column_scale[0..problem.num_cols]) |*published, internal, scale|
        published.* = internal / (scale * self.objective_scale);
    return .{
        .status = status,
        .primal = basis.published_primal,
        .dual = basis.published_dual,
        .reduced_cost = basis.published_reduced_cost,
        .unbounded_ray = if (self.unbounded_ray_valid) basis.unbounded_ray else &.{},
        .infeasibility_ray = if (self.infeasibility_ray_valid) basis.infeasibility_ray else &.{},
        .objective_value = self.objective_value,
        .iterations = self.iterations,
    };
}

test "engine detects an unbounded improving structural column" {
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{-1},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{-std.math.inf(f64)},
        .row_upper = &[_]f64{1},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &[_]usize{ 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    const status = engine.solveProblem(problem, .{});
    try std.testing.expectEqual(SolveStatus.unbounded, status);
    const solution = engine.solutionView(problem, status).?;
    try std.testing.expectEqual(@as(usize, 1), solution.unbounded_ray.len);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), solution.unbounded_ray[0], 1e-12);
}
