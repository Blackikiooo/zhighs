//! Pivot mechanics, reduced-cost maintenance, and factorization lifecycle
//! for `SimplexEngine`.
//!
//! ## Responsibility
//!
//! Owns direction computation (FTRAN) with residual checks, the pivot state
//! transition, rank-one and exact reduced-cost updates, refactorization
//! decisions, and basic-value recomputation.

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
const RebuildReason = @import("engine.zig").RebuildReason;

/// Compute the revised-simplex entering direction `B^-1 A_j` into the
/// engine-owned SoA workspace without allocating or copying model data.
pub fn computeDirection(self: *SimplexEngine, problem: problem_module.ProblemView, entering_col: usize) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    self.direction_requires_reinversion = false;
    if (self.fillInternalColumn(problem, entering_col, basis.pivot_direction) != .optimal) {
        self.failure_site = .direction_column;
        return .numerical_failure;
    }
    self.factorization.solveForUpdate(basis.pivot_direction) catch {
        self.failure_site = .direction_solve;
        return .numerical_failure;
    };
    if (!self.directionResidualAcceptable(problem, entering_col)) {
        self.factorization.recordReinversion(.solve_residual);
        if (self.refactorizeBasis(problem, .solve_residual) != .optimal) {
            self.failure_site = .direction_refactor;
            return .numerical_failure;
        }
        if (self.fillInternalColumn(problem, entering_col, basis.pivot_direction) != .optimal) return .numerical_failure;
        self.factorization.solveForUpdate(basis.pivot_direction) catch return .numerical_failure;
        var refinement_step: usize = 0;
        while (!self.directionResidualAcceptable(problem, entering_col)) : (refinement_step += 1) {
            if (refinement_step >= self.numerical.max_refinement_steps) {
                self.failure_site = .direction_refinement;
                return .numerical_failure;
            }
            self.factorization.solve(basis.residual_work) catch {
                self.failure_site = .direction_refinement;
                return .numerical_failure;
            };
            for (basis.pivot_direction, basis.residual_work) |*value, correction| value.* += correction;
            self.numerical.refinement_count += 1;
            // The captured partial aq predates refinement. Force a fresh
            // base factorization after this pivot instead of publishing an
            // FT update from inconsistent capture state.
            self.direction_requires_reinversion = true;
        }
    }
    self.observeAqDensity(basis.pivot_direction);
    self.dual_work.changed_row_count = self.factorization.gatherFtranResultIndices(
        basis.pivot_direction,
        self.dual_work.changed_row_index,
    );
    return .optimal;
}

/// Validate an updated FTRAN before its direction reaches ratio testing.
/// A drifting FT chain can otherwise select an invalid leaving row and
/// permanently corrupt the basis before the periodic growth gate fires.
pub fn directionResidualAcceptable(self: *SimplexEngine, problem: problem_module.ProblemView, entering_col: usize) bool {
    const basis = if (self.basis) |*value| value else return false;
    if (self.fillInternalColumn(problem, entering_col, basis.residual_work) != .optimal) return false;
    for (basis.rhs_work, basis.residual_work) |*magnitude, rhs| magnitude.* = @abs(rhs);
    self.subtractBasisProductWithMagnitude(
        problem,
        basis.pivot_direction,
        basis.residual_work,
        basis.rhs_work,
    ) catch return false;
    var residual_max: f64 = 0.0;
    var equation_scale: f64 = 0.0;
    for (basis.residual_work) |value| residual_max = @max(residual_max, @abs(value));
    for (basis.rhs_work) |value| equation_scale = @max(equation_scale, value);
    // Measure normwise backward error. Dividing only by |a_q| rejects a
    // backward-stable solve whenever a poorly scaled basis produces a
    // large x whose B*x terms cancel. |a_q| + |B|*|x| is the natural
    // scale of the equations and is also what iterative refinement can
    // meaningfully improve.
    const relative = residual_max / @max(1.0, equation_scale);
    self.numerical.last_ftran_relative_residual = relative;
    self.numerical.max_ftran_relative_residual = @max(self.numerical.max_ftran_relative_residual, relative);
    return std.math.isFinite(relative) and relative <= self.numerical.residual_tolerance;
}

/// Apply one primal pivot and rebuild the dense basis factorization in
/// existing storage. This is allocation-free; update factorizations can
/// replace the reinversion later without changing the state transition.
pub fn performPivot(self: *SimplexEngine, problem: problem_module.ProblemView, entering_col: usize, entering_direction: f64, leaving_row: usize, leaving_bound: basis_module.BasisStatus, step: f64) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    if (entering_col >= basis.col_status.len or leaving_row >= problem.num_rows) return .numerical_failure;
    const leaving_col = basis.basic_index[leaving_row];
    const pivot = basis.pivot_direction[leaving_row];
    if (self.pricing.rule == .steepest_edge and self.dual_row_index == @as(u32, @intCast(leaving_row))) {
        if (self.updateDualSteepestEdgeWeights(leaving_row, pivot) == .optimal) {
            self.stats.dual_dse_updates += 1;
            self.dual_dse_updates_since_start += 1;
            if (self.active_dual_edge_weight_strategy == .steepest_devex and
                self.active_dual_dse_update_budget != 0 and
                self.dual_dse_updates_since_start >= self.active_dual_dse_update_budget)
                self.switchDualDseToDevex(.budget);
        } else if (self.active_dual_edge_weight_strategy == .steepest_devex) {
            self.switchDualDseToDevex(.invalid);
        } else {
            self.dual_edge_weights_valid = false;
        }
    } else {
        self.dual_edge_weights_valid = false;
        if (self.algorithm == .dual_revised and self.pricing.rule == .devex and
            self.active_dual_edge_weight_strategy == .steepest_devex)
            self.updateDualDevexWeights(leaving_row, pivot)
        else if (self.algorithm == .dual_revised or self.active_devex_strategy == .legacy)
            self.updateLegacyDevexWeights(entering_col, leaving_row, pivot);
    }
    self.dual_row_index = null;
    var forced_rebuild_reason: ?RebuildReason = null;
    const update_succeeded = blk: {
        if (self.direction_requires_reinversion) {
            forced_rebuild_reason = .direction_refinement;
            break :blk false;
        }
        if (self.fresh_factorization_pivots_remaining != 0) {
            forced_rebuild_reason = .fresh_mode;
            break :blk false;
        }
        self.factorization.update(.{
            .leaving_row = @intCast(leaving_row),
            .entering_col = @intCast(entering_col),
            .direction = basis.pivot_direction,
            .column_scale = entering_direction,
        }) catch |err| {
            const failure_kind = self.factorization.recordUpdateFailure(err);
            if (self.stats.first_update_failure_kind == null) {
                self.stats.first_update_failure_kind = failure_kind;
                self.stats.first_update_failure_iteration = self.iterations;
                self.stats.first_update_failure_entering = @intCast(entering_col);
                self.stats.first_update_failure_leaving_row = @intCast(leaving_row);
            }
            forced_rebuild_reason = if (self.cleanup_active) .cleanup else .update_rejected;
            break :blk false;
        };
        break :blk true;
    };
    self.direction_requires_reinversion = false;
    self.numerical.observePivot(pivot);
    for (basis.basic_value, basis.pivot_direction) |*value, direction| value.* -= step * direction;
    const entering_value = basis.primal[entering_col] + entering_direction * step;
    basis.basic_value[leaving_row] = entering_value;
    basis.basic_lower[leaving_row] = basis.col_lower[entering_col];
    basis.basic_upper[leaving_row] = basis.col_upper[entering_col];
    basis.primal[entering_col] = entering_value;
    basis.primal[leaving_col] = if (leaving_bound == .at_upper) basis.col_upper[leaving_col] else basis.col_lower[leaving_col];
    basis.applyPivot(leaving_row, entering_col, leaving_bound) catch return .numerical_failure;
    // A committed basis change invalidates rebuild-derived solution state.
    // If the factor update below triggers reinversion, refactorizeBasis marks
    // it fresh again and the dual controller completes the rebuild pipeline.
    self.dual_control.has_fresh_rebuild = false;
    self.dual_work.notePivot(entering_col, leaving_col, leaving_bound);
    self.dual_work.updatePrimalInfeasibilityList(
        basis,
        self.dual_work.changed_row_index[0..self.dual_work.changed_row_count],
        self.numerical.primal_tolerance,
    );
    if (self.devex_reset_after_pivot) {
        self.initializePrimalDevexFramework();
        self.devex_reset_after_pivot = false;
    }
    for (basis.basic_index, basis.basic_value) |global_col, value| basis.primal[global_col] = value;
    if (self.pivot_trace_count < self.active_pivot_trace.len) {
        self.active_pivot_trace[self.pivot_trace_count] = .{
            .phase = self.current_phase,
            .iteration = self.iterations,
            .entering_column = @intCast(entering_col),
            .leaving_column = leaving_col,
            .leaving_row = @intCast(leaving_row),
            .pivot = pivot,
            .step = step,
            .update_count = self.factorization.update_count,
            .ftran_relative_residual = self.numerical.last_ftran_relative_residual,
            .condition_estimate = self.numerical.pivot_condition_estimate,
            .bound_flip_count = 0,
        };
        self.pivot_trace_count += 1;
    }
    const reinversion_reason = self.factorization.reinversionReason(self.numerical.max_update_count);
    if (reinversion_reason) |reason| self.factorization.recordReinversion(reason);
    if (!update_succeeded or reinversion_reason != null or self.numerical.needsRefactor()) {
        const rebuild_reason = forced_rebuild_reason orelse if (reinversion_reason) |reason| switch (reason) {
            .update_limit => RebuildReason.update_limit,
            .update_growth => RebuildReason.update_growth,
            .solve_residual => RebuildReason.solve_residual,
            .small_pivot => RebuildReason.small_pivot,
        } else RebuildReason.numerical_policy;
        const status = self.refactorizeBasis(problem, rebuild_reason);
        if (status == .optimal and rebuild_reason == .fresh_mode)
            self.fresh_factorization_pivots_remaining -|= 1;
        if (status != .optimal and self.failure_site == .none) self.failure_site = .pivot_update;
        return status;
    }
    return .optimal;
}

/// Recompute reduced costs from persistent dual `workCost + workShift`.
///
/// This is the exact reprice used by dual Phase I/II rebuilds; it overwrites
/// the complete internal reduced-cost vector.
pub fn recomputeReducedCostsFromWork(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    const costs = self.dual_work.work_cost;
    const shifts = self.dual_work.work_shift;
    const original_count = problem.num_cols + problem.num_rows;
    if (costs.len < original_count or shifts.len < original_count) return .numerical_failure;
    for (basis.basic_index, 0..) |column, row| basis.dual[row] = if (column < original_count) costs[column] + shifts[column] else 0.0;
    self.factorization.solveTranspose(basis.dual) catch return .numerical_failure;
    const pricing_started = self.statisticsTimestamp();
    defer self.recordPricingElapsed(pricing_started);
    for (0..problem.num_cols) |column| {
        var reduced = costs[column] + shifts[column];
        const begin = problem.matrix.col_starts[column];
        const end = problem.matrix.col_starts[column + 1];
        for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
            if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
            const row_index = row.toUsize();
            reduced -= basis.row_scale[row_index] * coefficient * basis.column_scale[column] * basis.dual[row_index];
        }
        basis.reduced_cost[column] = reduced;
    }
    for (0..problem.num_rows) |row| basis.reduced_cost[problem.num_cols + row] = costs[problem.num_cols + row] + shifts[problem.num_cols + row] - basis.dual[row];
    @memset(basis.reduced_cost[original_count..], 0.0);
    self.reduced_cost_update_count = 0;
    return .optimal;
}

/// Exact original-objective reprice with drift observation. This is the
/// safety boundary for dual rank-1 updates: exact values replace the
/// incrementally maintained vector before feasibility is revalidated.
pub fn recomputeReducedCostsWithDrift(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    const count = basis.reduced_cost.len;
    if (self.degeneracy.reduced_cost_snapshot.len < count) return .numerical_failure;
    @memcpy(self.degeneracy.reduced_cost_snapshot[0..count], basis.reduced_cost);
    if (self.recomputeReducedCosts(problem) != .optimal) return .numerical_failure;
    var maximum_drift: f64 = 0.0;
    for (basis.reduced_cost, self.degeneracy.reduced_cost_snapshot[0..count]) |exact, updated| {
        maximum_drift = @max(maximum_drift, @abs(exact - updated) / @max(1.0, @abs(exact)));
    }
    self.maximum_reduced_cost_drift = @max(self.maximum_reduced_cost_drift, maximum_drift);
    self.exact_reprices += 1;
    self.stats.dual_exact_reprices += 1;
    if (self.active_adaptive_reprice) {
        if (maximum_drift > self.numerical.dual_tolerance * 10.0)
            self.reduced_cost_refresh_period = @max(self.reduced_cost_refresh_period / 2, 1)
        else if (maximum_drift < self.numerical.dual_tolerance * 0.1)
            self.reduced_cost_refresh_period = @min(self.reduced_cost_refresh_period * 2, 8);
    }
    return .optimal;
}

/// Update all reduced costs after a primal basis replacement using the
/// pivotal row of the old basis: `r' = r - (r_q / alpha_pq) * alpha_p`.
/// One BTRAN and one CSC scan replace rebuilding `c_B`, solving for the
/// full dual vector, and repricing from scratch. A periodic exact refresh
/// bounds accumulated roundoff while keeping the hot path allocation-free.
pub fn updateReducedCostsAfterPrimalPivot(
    self: *SimplexEngine,
    problem: problem_module.ProblemView,
    entering_col: usize,
    entering_direction: f64,
    leaving_row: u32,
) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    const row: usize = @intCast(leaving_row);
    if (row >= problem.num_rows or entering_col >= basis.reduced_cost.len) return .numerical_failure;
    const actual_pivot = basis.pivot_direction[row] * entering_direction;
    if (!std.math.isFinite(actual_pivot) or @abs(actual_pivot) <= self.numerical.pivot_tolerance)
        return .numerical_failure;
    const theta = basis.reduced_cost[entering_col] / actual_pivot;
    if (!std.math.isFinite(theta)) return .numerical_failure;

    @memset(basis.residual_work, 0.0);
    basis.residual_work[row] = 1.0;
    self.factorization.solveTranspose(basis.residual_work) catch return .numerical_failure;
    const row_pricing = self.selectRowPricing(basis.residual_work);
    const pricing_started = self.statisticsTimestamp();
    defer if (row_pricing)
        self.recordRowPricingElapsed(pricing_started)
    else
        self.recordPricingElapsed(pricing_started);
    if (row_pricing) {
        @memset(basis.tableau[0..problem.num_cols], 0.0);
        for (basis.residual_work, 0..) |dual_value, row_index| {
            if (@abs(dual_value) <= self.numerical.zero_tolerance) continue;
            const scaled_dual = dual_value * basis.row_scale[row_index];
            for (self.pricing_row_view.rowColumns(row_index), self.pricing_row_view.rowValues(row_index)) |column_u32, coefficient| {
                const column: usize = @intCast(column_u32);
                basis.tableau[column] += scaled_dual * coefficient * basis.column_scale[column];
            }
        }
        for (basis.reduced_cost[0..problem.num_cols], basis.tableau[0..problem.num_cols]) |*reduced, tableau_entry|
            reduced.* -= theta * tableau_entry;
    } else {
        for (0..problem.num_cols) |column| {
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            var tableau_entry: f64 = 0.0;
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |matrix_row, coefficient| {
                if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                const row_index = matrix_row.toUsize();
                tableau_entry += basis.residual_work[row_index] * basis.row_scale[row_index] * coefficient * basis.column_scale[column];
            }
            basis.tableau[column] = tableau_entry;
            basis.reduced_cost[column] -= theta * tableau_entry;
        }
    }
    const artificial_begin = problem.num_cols + problem.num_rows;
    for (0..problem.num_rows) |logical_row| {
        const multiplier = theta * basis.residual_work[logical_row];
        basis.tableau[problem.num_cols + logical_row] = basis.residual_work[logical_row];
        basis.tableau[artificial_begin + logical_row] = basis.residual_work[logical_row] * basis.artificial_sign[logical_row];
        basis.reduced_cost[problem.num_cols + logical_row] -= multiplier;
        basis.reduced_cost[artificial_begin + logical_row] -= multiplier * basis.artificial_sign[logical_row];
    }
    basis.reduced_cost[entering_col] = 0.0;
    for (basis.reduced_cost) |value| if (!std.math.isFinite(value)) return .numerical_failure;
    if (self.pricing.rule == .devex and self.active_devex_strategy == .framework and
        self.updatePrimalDevexFramework(entering_col, row, actual_pivot) != .optimal)
        return .numerical_failure;
    self.reduced_cost_update_count += 1;
    return .optimal;
}

/// Recompute reduced costs exactly from the original scaled objective.
pub fn recomputeReducedCosts(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    const maximize = problem.objective_sense == .maximize;
    for (basis.basic_index, 0..) |global_col, row| {
        basis.dual[row] = if (global_col < problem.num_cols)
            (if (maximize) -problem.col_cost[global_col] else problem.col_cost[global_col]) * self.objective_scale * basis.column_scale[global_col]
        else
            0.0;
    }
    self.factorization.solveTranspose(basis.dual) catch return .numerical_failure;
    const row_pricing = self.selectRowPricing(basis.dual);
    const pricing_started = self.statisticsTimestamp();
    defer if (row_pricing)
        self.recordRowPricingElapsed(pricing_started)
    else
        self.recordPricingElapsed(pricing_started);
    for (0..problem.num_cols) |col|
        basis.reduced_cost[col] = (if (maximize) -problem.col_cost[col] else problem.col_cost[col]) * self.objective_scale * basis.column_scale[col];
    if (row_pricing) {
        for (basis.dual, 0..) |dual_value, row_index| {
            if (@abs(dual_value) <= self.numerical.zero_tolerance) continue;
            const scaled_dual = dual_value * basis.row_scale[row_index];
            for (self.pricing_row_view.rowColumns(row_index), self.pricing_row_view.rowValues(row_index)) |column_u32, coefficient| {
                const column: usize = @intCast(column_u32);
                basis.reduced_cost[column] -= scaled_dual * coefficient * basis.column_scale[column];
            }
        }
    } else {
        for (0..problem.num_cols) |col| {
            const begin = problem.matrix.col_starts[col];
            const end = problem.matrix.col_starts[col + 1];
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, value| {
                if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(value)) continue;
                const row_index = row.toUsize();
                basis.reduced_cost[col] -= basis.row_scale[row_index] * value * basis.column_scale[col] * basis.dual[row_index];
            }
        }
    }
    for (0..problem.num_rows) |row| basis.reduced_cost[problem.num_cols + row] = -basis.dual[row];
    @memset(basis.reduced_cost[problem.num_cols + problem.num_rows ..], 0.0);
    self.reduced_cost_update_count = 0;
    return .optimal;
}

/// Rebuild the current basis factors and record the classified reason.
///
/// On success this also completes all factor-dependent state refresh performed
/// by `finishRefactorization`.
pub fn refactorizeBasis(self: *SimplexEngine, problem: problem_module.ProblemView, reason: RebuildReason) SolveStatus {
    _ = self.dual_control.beginRebuild();
    const rebuild_started = self.statisticsTimestamp();
    defer {
        self.stats.rebuild_calls += 1;
        self.stats.rebuild_ns = std.math.add(u64, self.stats.rebuild_ns, self.elapsedSince(rebuild_started)) catch std.math.maxInt(u64);
    }
    switch (reason) {
        .phase_one_setup => self.stats.rebuild_phase_one_setup += 1,
        .solve_residual => self.stats.rebuild_solve_residual += 1,
        .small_pivot => self.stats.rebuild_small_pivot += 1,
        .update_limit => self.stats.rebuild_update_limit += 1,
        .update_growth => self.stats.rebuild_update_growth += 1,
        .direction_refinement => self.stats.rebuild_direction_refinement += 1,
        .fresh_mode => self.stats.rebuild_fresh_mode += 1,
        .update_rejected => self.stats.rebuild_update_rejected += 1,
        .cleanup => self.stats.rebuild_cleanup += 1,
        .edge_weight_reset => self.stats.rebuild_edge_weight_reset += 1,
        .numerical_policy => self.stats.rebuild_numerical_policy += 1,
    }
    self.factorizeCurrentBasis(problem) catch {
        self.failure_site = .pivot_factorization;
        return .numerical_failure;
    };
    const status = self.finishRefactorization();
    if (status == .optimal) self.dual_control.finishRebuild();
    return status;
}

/// Assemble and factorize the matrix columns named by the current basis head.
pub fn factorizeCurrentBasis(self: *SimplexEngine, problem: problem_module.ProblemView) factorization_module.FactorizationError!void {
    const basis = if (self.basis) |*value| value else return error.DimensionMismatch;
    return self.factorization.factorizeBasis(
        problem.matrix,
        basis.basic_index,
        basis.row_scale,
        basis.column_scale[0..problem.num_cols],
        basis.artificial_sign,
    );
}

/// Reset update-dependent numerical state and rebuild exact edge weights.
pub fn finishRefactorization(self: *SimplexEngine) SolveStatus {
    self.numerical.markRefactorized();
    self.observeFactorizationStability();
    self.dual_edge_weights_valid = false;
    self.dual_row_index = null;
    if (self.pricing.rule == .steepest_edge and self.ensureExactDualEdgeWeights() != .optimal) {
        self.failure_site = .pivot_edge_weights;
        return .numerical_failure;
    }
    return .optimal;
}

/// Feed the current factor pivot spread into the numerical warning policy.
pub fn observeFactorizationStability(self: *SimplexEngine) void {
    self.numerical.pivot_condition_estimate = self.factorization.pivotConditionEstimate();
    if (!std.math.isFinite(self.numerical.pivot_condition_estimate) or self.numerical.pivot_condition_estimate > 1e12)
        self.numerical.numerical_warning = true;
}

/// Recompute the primal basic solution from the immutable problem and
/// current nonbasic values. Uses engine-owned workspace and performs no
/// allocation; this limits drift from long Eta update chains.
pub fn recomputeBasicValues(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
    return self.recomputeBasicValuesImpl(problem, true);
}

/// Recompute basic values without rejecting primal bound violations.
pub fn recomputeBasicValuesUnchecked(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
    return self.recomputeBasicValuesImpl(problem, false);
}

/// Reconstruct `B^-1 rhs`, optionally rejecting basic values outside bounds.
pub fn recomputeBasicValuesImpl(self: *SimplexEngine, problem: problem_module.ProblemView, enforce_bounds: bool) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    @memcpy(basis.residual_work, basis.row_rhs);
    for (basis.col_status, basis.primal, 0..) |status, value, column| {
        if (status == .basic or value == 0.0) continue;
        if (column < problem.num_cols) {
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
                if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                const row_index = row.toUsize();
                basis.residual_work[row_index] -= basis.row_scale[row_index] * coefficient * basis.column_scale[column] * value;
            }
        } else if (column < problem.num_cols + problem.num_rows) {
            basis.residual_work[column - problem.num_cols] -= value;
        } else {
            const row = column - problem.num_cols - problem.num_rows;
            basis.residual_work[row] -= basis.artificial_sign[row] * value;
        }
    }
    @memcpy(basis.rhs_work, basis.residual_work);
    self.factorization.solve(basis.rhs_work) catch return .numerical_failure;

    var refinement_step: usize = 0;
    while (true) {
        @memcpy(basis.pivot_direction, basis.residual_work);
        self.subtractBasisProduct(problem, basis.rhs_work, basis.pivot_direction) catch return .numerical_failure;
        var residual_max: f64 = 0.0;
        var rhs_max: f64 = 0.0;
        for (basis.pivot_direction, basis.residual_work) |residual, rhs_value| {
            residual_max = @max(residual_max, @abs(residual));
            rhs_max = @max(rhs_max, @abs(rhs_value));
        }
        self.numerical.observeResidual(residual_max, rhs_max);
        if (!std.math.isFinite(residual_max)) return .numerical_failure;
        if (self.numerical.last_relative_residual <= self.numerical.residual_tolerance or
            refinement_step >= self.numerical.max_refinement_steps)
            break;
        self.factorization.solve(basis.pivot_direction) catch return .numerical_failure;
        for (basis.rhs_work, basis.pivot_direction) |*value, correction| value.* += correction;
        refinement_step += 1;
        self.numerical.refinement_count += 1;
    }
    for (basis.basic_index, basis.basic_value, basis.rhs_work, 0..) |column, *value, recomputed, row| {
        value.* = recomputed;
        basis.primal[column] = recomputed;
        if (enforce_bounds and (recomputed < basis.basic_lower[row] - self.numerical.primal_tolerance or
            recomputed > basis.basic_upper[row] + self.numerical.primal_tolerance))
            return .numerical_failure;
    }
    return .optimal;
}

/// Compute `residual -= B * x` directly from the borrowed CSC model and
/// basis membership. No dense basis copy or temporary column is required.
pub fn subtractBasisProduct(self: *SimplexEngine, problem: problem_module.ProblemView, x: []const f64, residual: []f64) !void {
    const basis = if (self.basis) |*value| value else return error.NumericalFailure;
    if (x.len != problem.num_rows or residual.len != problem.num_rows) return error.NumericalFailure;
    for (basis.basic_index, x) |global_col, value| {
        if (value == 0.0) continue;
        const column: usize = @intCast(global_col);
        if (column < problem.num_cols) {
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
                if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                const row_index = row.toUsize();
                if (row_index >= residual.len) return error.NumericalFailure;
                residual[row_index] -= basis.row_scale[row_index] * coefficient * basis.column_scale[column] * value;
            }
        } else if (column < problem.num_cols + problem.num_rows) {
            residual[column - problem.num_cols] -= value;
        } else if (column < problem.num_cols + 2 * problem.num_rows) {
            const row = column - problem.num_cols - problem.num_rows;
            residual[row] -= basis.artificial_sign[row] * value;
        } else return error.NumericalFailure;
    }
}

/// Subtract `B*x` and accumulate `|B|*|x|` in the same sparse traversal.
/// The magnitude output lets FTRAN use a scale-aware backward-error test
/// without constructing a dense basis or making another matrix pass.
pub fn subtractBasisProductWithMagnitude(
    self: *SimplexEngine,
    problem: problem_module.ProblemView,
    x: []const f64,
    residual: []f64,
    magnitude: []f64,
) !void {
    const basis = if (self.basis) |*value| value else return error.NumericalFailure;
    if (x.len != problem.num_rows or residual.len != problem.num_rows or magnitude.len != problem.num_rows)
        return error.NumericalFailure;
    for (basis.basic_index, x) |global_col, value| {
        if (!std.math.isFinite(value)) return error.NumericalFailure;
        const column: usize = @intCast(global_col);
        if (column < problem.num_cols) {
            const start = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            for (problem.matrix.row_indices[start..end], problem.matrix.values[start..end]) |row, coefficient| {
                if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                const row_index = row.toUsize();
                if (row_index >= residual.len) return error.NumericalFailure;
                const product = basis.row_scale[row_index] * coefficient * basis.column_scale[column] * value;
                residual[row_index] -= product;
                magnitude[row_index] += @abs(product);
            }
        } else if (column < problem.num_cols + problem.num_rows) {
            const row = column - problem.num_cols;
            residual[row] -= value;
            magnitude[row] += @abs(value);
        } else {
            const row = column - problem.num_cols - problem.num_rows;
            if (row >= residual.len) return error.NumericalFailure;
            const product = basis.artificial_sign[row] * value;
            residual[row] -= product;
            magnitude[row] += @abs(product);
        }
    }
}
