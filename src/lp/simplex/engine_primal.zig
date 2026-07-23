//! Primal revised-simplex loop, primal Phase I, and primal Devex weights for
//! `SimplexEngine`.
//!
//! ## Responsibility
//!
//! Owns the primal Phase-I/Phase-II iteration loops, leaving-variable
//! selection (including the Bland fallback), primal entering pricing, the
//! primal Devex weight frameworks, and artificial-basis cleanup.

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
const PhaseOneStrategy = @import("engine.zig").PhaseOneStrategy;
const SolveStatus = @import("engine.zig").SolveStatus;
const DegeneracyTraceEvent = @import("engine.zig").DegeneracyTraceEvent;
const PrimalLeavingResult = @import("engine.zig").PrimalLeavingResult;
const SolveControl = @import("engine.zig").SolveControl;

/// Choose the leaving basic row for the currently materialized pivot
/// direction. Returns `unbounded` when no positive pivot coefficient
/// limits the entering variable.
pub fn chooseLeaving(self: *SimplexEngine) PrimalLeavingResult {
    return self.chooseLeavingWithPolicy(true);
}

/// Apply the primal ratio test to the current FTRAN direction.
///
/// For every basic row this converts the signed direction into a nonnegative
/// limiting direction and the distance to the bound that would be hit. The
/// active degeneracy policy then chooses Harris, rank-perturbed, or (when
/// permitted) Bland ordering. `allow_bland` is disabled for the generalized
/// bounded Phase-I problem, where the assumptions behind the standard-form
/// Bland proof do not hold.
pub fn chooseLeavingWithPolicy(self: *SimplexEngine, allow_bland: bool) PrimalLeavingResult {
    const basis = if (self.basis) |*value| value else return .{ .status = .numerical_failure };
    for (basis.basic_margin, basis.ratio_direction, basis.basic_value, basis.basic_lower, basis.basic_upper, basis.pivot_direction) |*margin, *ratio_direction, value, lower, upper, direction| {
        if (direction > self.numerical.zero_tolerance) {
            margin.* = value - lower;
            ratio_direction.* = direction;
        } else if (direction < -self.numerical.zero_tolerance) {
            margin.* = upper - value;
            ratio_direction.* = -direction;
        } else {
            margin.* = 0.0;
            ratio_direction.* = 0.0;
        }
    }
    // Explicit baseline retains the proven Bland fallback (notably
    // feasible starts such as kb2). Automatic mode keeps Harris during its
    // 256-pivot hysteresis window, then uses bounded perturbation; this
    // avoids changing short local ties into weak Bland pivots.
    if (allow_bland and self.numerical.anti_cycling_active and
        self.active_degeneracy_strategy == .baseline)
        return self.chooseLeavingBland();
    const choice = if (self.degeneracy.active)
        self.ratio_test.chooseLeavingPerturbed(
            basis.ratio_direction,
            basis.basic_margin,
            self.degeneracy.row_rank[0..basis.basic_margin.len],
            self.numerical.primal_tolerance,
        )
    else
        self.ratio_test.chooseLeaving(basis.ratio_direction, basis.basic_margin);
    if (choice.row == null) return .{ .status = .unbounded };
    const row = choice.row.?;
    const bound: basis_module.BasisStatus = if (basis.pivot_direction[@intCast(row)] > 0) .at_lower else .at_upper;
    return .{ .status = .optimal, .row = row, .step = choice.step, .bound = bound };
}

/// Choose the limiting basic variable with Bland-compatible tie breaking.
///
/// The minimum nonnegative step wins; ties within the configured perturbation
/// tolerance are resolved by the global basic-column index rather than row
/// position. A scale-aware pivot threshold excludes coefficients too small to
/// support a numerically meaningful basis exchange.
pub fn chooseLeavingBland(self: *SimplexEngine) PrimalLeavingResult {
    const basis = if (self.basis) |*value| value else return .{ .status = .numerical_failure };
    var best_row: ?u32 = null;
    var best_column: u32 = std.math.maxInt(u32);
    var best_step = std.math.inf(f64);
    const tie_tolerance = self.numerical.perturbation;
    var maximum_direction: f64 = 0.0;
    for (basis.ratio_direction) |direction| maximum_direction = @max(maximum_direction, direction);
    const pivot_tolerance = @max(self.numerical.zero_tolerance, self.ratio_test.tolerance * maximum_direction);
    for (basis.ratio_direction, basis.basic_margin, basis.basic_index, 0..) |direction, margin, basic_column, row| {
        if (direction <= pivot_tolerance) continue;
        const step = @max(margin / direction, 0.0);
        if (!std.math.isFinite(step)) continue;
        if (step < best_step - tie_tolerance or
            (@abs(step - best_step) <= tie_tolerance and basic_column < best_column))
        {
            best_step = step;
            best_column = basic_column;
            best_row = @intCast(row);
        }
    }
    const row = best_row orelse return .{ .status = .unbounded };
    const bound: basis_module.BasisStatus = if (basis.pivot_direction[row] > 0.0) .at_lower else .at_upper;
    return .{ .status = .optimal, .row = row, .step = best_step, .bound = bound };
}

/// Forward-error guard for a pivotal entry extracted from an updated
/// FTRAN. `sqrt(epsilon)` is the boundary where cancellation can leave
/// fewer than half the significand bits. The unit floor also catches an
/// entire entering column that is small before column equilibration.
pub fn pivotNeedsFreshFactorization(self: *const SimplexEngine, leaving_row: u32) bool {
    const basis = if (self.basis) |*value| value else return true;
    const row: usize = @intCast(leaving_row);
    if (row >= basis.pivot_direction.len) return true;
    var maximum: f64 = 0.0;
    for (basis.pivot_direction) |value| maximum = @max(maximum, @abs(value));
    return @abs(basis.pivot_direction[row]) <= @sqrt(std.math.floatEps(f64)) * @max(1.0, maximum);
}

/// Lightweight Devex reference update. It deliberately uses the already
/// hot FTRAN direction and does not allocate. Exact steepest-edge weights
/// can replace this policy without changing basis storage or pivot code.
pub fn updateLegacyDevexWeights(self: *SimplexEngine, entering_col: usize, leaving_row: usize, pivot: f64) void {
    const basis = if (self.basis) |*value| value else return;
    if (entering_col >= basis.col_edge_weight.len or leaving_row >= basis.row_edge_weight.len) return;
    var norm_squared: f64 = 1.0;
    for (basis.pivot_direction) |value| norm_squared += value * value;
    if (!std.math.isFinite(norm_squared)) norm_squared = 1.0;
    const entering_weight = @max(1.0, norm_squared);
    basis.col_edge_weight[entering_col] = entering_weight;
    basis.row_edge_weight[leaving_row] = if (@abs(pivot) > self.numerical.pivot_tolerance)
        @max(1.0, entering_weight / (pivot * pivot))
    else
        entering_weight;

    if (self.pricing.devex_reset_period != 0 and self.pricing.iterations % self.pricing.devex_reset_period == 0) {
        @memset(basis.col_edge_weight, 1.0);
        @memset(basis.row_edge_weight, 1.0);
    }
}

/// Freeze the current nonbasic set as a primal Devex reference framework.
/// Subsequent pivots retain these bits until accumulated weight error
/// triggers a deterministic reset after the combinatorial pivot commits.
pub fn initializePrimalDevexFramework(self: *SimplexEngine) void {
    const basis = if (self.basis) |*value| value else return;
    @memset(basis.col_edge_weight, 1.0);
    for (basis.devex_reference, basis.col_status) |*reference, status|
        reference.* = @intFromBool(status != .basic);
    self.devex_framework_iterations = 0;
    self.devex_bad_weight_count = 0;
    self.stats.devex_frameworks += 1;
}

/// Full primal Devex reference recurrence. The pivotal norm is evaluated
/// from the frozen reference set and the hot FTRAN direction, then every
/// nonbasic weight touched by the complete old-basis tableau row is raised
/// to its framework lower bound. No candidate list or allocation is used.
pub fn updatePrimalDevexFramework(
    self: *SimplexEngine,
    entering_col: usize,
    leaving_row: usize,
    pivot: f64,
) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    if (entering_col >= basis.col_edge_weight.len or leaving_row >= basis.basic_index.len or
        @abs(pivot) <= self.numerical.pivot_tolerance)
        return .numerical_failure;
    var pivotal_norm: f64 = @floatFromInt(basis.devex_reference[entering_col]);
    for (basis.basic_index, basis.pivot_direction) |basic_column, direction| {
        if (basis.devex_reference[basic_column] != 0) pivotal_norm += direction * direction;
    }
    if (!std.math.isFinite(pivotal_norm) or pivotal_norm <= 0.0) return .numerical_failure;
    if (basis.col_edge_weight[entering_col] > 3.0 * pivotal_norm) {
        self.devex_bad_weight_count += 1;
        self.stats.devex_bad_weights += 1;
    }
    const new_pivotal_weight = pivotal_norm / (pivot * pivot);
    if (!std.math.isFinite(new_pivotal_weight)) return .numerical_failure;
    for (basis.col_edge_weight, basis.devex_reference, basis.tableau, basis.col_status) |*weight, reference, alpha, status| {
        if (status == .basic) continue;
        const candidate = new_pivotal_weight * alpha * alpha + @as(f64, @floatFromInt(reference));
        if (!std.math.isFinite(candidate)) return .numerical_failure;
        weight.* = @max(weight.*, candidate);
    }
    const leaving_col: usize = @intCast(basis.basic_index[leaving_row]);
    basis.col_edge_weight[leaving_col] = @max(1.0, new_pivotal_weight);
    basis.col_edge_weight[entering_col] = 1.0;
    self.devex_framework_iterations += 1;
    self.stats.devex_framework_updates += 1;
    if (self.devex_bad_weight_count > 3) self.devex_reset_after_pivot = true;
    return .optimal;
}

/// Run primal revised-simplex Phase II on the installed feasible basis.
///
/// The loop performs reduced-cost pricing, FTRAN, the primal ratio test,
/// optional bound exchange, and basis pivots. It also owns Devex framework
/// lifetime, degeneracy/taboo retries, exact-reprice cadence, and the fresh
/// reinversion required before an apparent optimum can be certified.
/// Terminal statuses therefore refer to the original objective and a
/// numerically checked basis, not merely an incrementally updated work state.
pub fn solvePrimal(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
    const phase_started = self.statisticsTimestamp();
    const iteration_started = self.iterations;
    const saved_pricing_rule = self.pricing.rule;
    if (self.active_primal_pricing_strategy == .partial) self.pricing.rule = .partial;
    defer self.pricing.rule = saved_pricing_rule;
    defer self.recordPhaseElapsed(.phase_two, phase_started);
    defer self.recordPhaseIterations(.phase_two, iteration_started);
    if (self.recomputeReducedCosts(problem) != .optimal) {
        self.failure_site = .reduced_cost;
        return .numerical_failure;
    }
    if (self.pricing.rule == .devex and self.active_devex_strategy == .framework)
        self.initializePrimalDevexFramework();
    while (self.iterations < control.max_iterations) : (self.iterations += 1) {
        if (self.beginIteration(problem, control, .phase_two)) |status| return status;
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        const original_cols = problem.num_cols + problem.num_rows;
        const pricing_tolerance = self.scaledDualTolerance(problem);
        if (self.prepareDegeneracyPolicy(basis.basic_value.len, original_cols) != .optimal)
            return .numerical_failure;
        var entering_candidate = if (self.active_degeneracy_strategy == .baseline)
            self.choosePrimalEnteringTimed(original_cols, pricing_tolerance)
        else
            self.choosePrimalEnteringWeightedTimed(original_cols, pricing_tolerance, true);
        if (entering_candidate == null and self.active_degeneracy_strategy == .perturbation_taboo and
            self.degeneracy.active)
        {
            self.degeneracy.invalidateTaboo();
            self.stats.taboo_retries += 1;
            entering_candidate = self.choosePrimalEnteringWeightedTimed(original_cols, pricing_tolerance, false);
        }
        const entering = entering_candidate orelse {
            // Never certify dual feasibility from an updated BTRAN chain.
            // Reinvert the unchanged basis once and price again; a stale
            // transpose solve can otherwise publish a false optimum.
            if (self.factorization.update_count != 0) {
                self.factorization.recordReinversion(.solve_residual);
                if (self.refactorizeBasis(problem, .solve_residual) != .optimal) return .numerical_failure;
                if (self.recomputeReducedCosts(problem) != .optimal) {
                    self.failure_site = .reduced_cost;
                    return .numerical_failure;
                }
                continue;
            }
            const status = self.finishOptimal(problem);
            if (status == .numerical_failure) {
                if (self.active_degeneracy_strategy != .baseline) {
                    self.stats.cold_restart_solves += 1;
                    return self.restartSolveWithoutPerturbation(problem, control);
                }
                // finishOptimal still fails at baseline degeneracy.
                // When the Devex framework changes the pivot path into a
                // basis that fails the residual check, fall back to the
                // well-tested legacy pricing with a clean cold start.
                if (self.active_devex_strategy == .framework) {
                    self.stats.cold_restart_solves += 1;
                    var legacy_control = control;
                    legacy_control.initial_basis = null;
                    legacy_control.degeneracy_strategy = .baseline;
                    legacy_control.devex_strategy = .legacy;
                    return self.solveProblem(problem, legacy_control);
                }
                self.failure_site = .optimality_check;
            }
            return status;
        };
        if (self.computeDirection(problem, entering.column) != .optimal) return .numerical_failure;
        if (entering.direction < 0) {
            for (basis.pivot_direction) |*value| value.* = -value.*;
        }
        var leaving = self.chooseLeaving();
        var small_pivot_retry = false;
        // Backward error bounds the full solve but not a tiny pivotal
        // component's forward error. Recompute suspicious pivots from a
        // fresh basis before mutating basis membership.
        if (leaving.status == .optimal and self.pivotNeedsFreshFactorization(leaving.row.?)) {
            small_pivot_retry = true;
            if (self.factorization.update_count != 0) {
                // Once a basis exposes a component at the forward-accuracy
                // boundary, occasional FT updates can steer later ratio
                // tests onto a different numerical path. Recompute from a
                // fresh basis before deciding whether the pivot is real.
                self.fresh_factorization_pivots_remaining = self.numerical.fresh_factorization_recovery_pivots;
                self.factorization.recordReinversion(.small_pivot);
                if (self.refactorizeBasis(problem, .small_pivot) != .optimal) return .numerical_failure;
                if (self.computeDirection(problem, entering.column) != .optimal) return .numerical_failure;
                if (entering.direction < 0) {
                    for (basis.pivot_direction) |*value| value.* = -value.*;
                }
                leaving = self.chooseLeaving();
            }
            // A warning that survives fresh FTRAN is not update drift.
            // Keep Bland's entering choice, but let Harris replace the
            // numerically weak leaving row before basis membership is
            // mutated. This prevents a known near-rank pivot from being
            // handed to reinversion as an already committed basis.
            if (leaving.status == .optimal and self.pivotNeedsFreshFactorization(leaving.row.?))
                leaving = self.chooseLeavingWithPolicy(false);
        }
        const own_step = if (entering.direction > 0)
            basis.col_upper[entering.column] - basis.primal[entering.column]
        else
            basis.primal[entering.column] - basis.col_lower[entering.column];
        if (std.math.isFinite(own_step) and (leaving.status == .unbounded or own_step < leaving.step)) {
            for (basis.basic_value, basis.pivot_direction) |*value, direction| value.* -= own_step * direction;
            basis.primal[entering.column] = if (entering.direction > 0) basis.col_upper[entering.column] else basis.col_lower[entering.column];
            basis.col_status[entering.column] = if (entering.direction > 0) .at_upper else .at_lower;
            self.observeIterationStep(
                own_step,
                entering.column,
                entering.direction,
                null,
                0,
                null,
                @abs(basis.reduced_cost[entering.column] * own_step),
                small_pivot_retry,
                0,
            );
        } else {
            if (leaving.status == .unbounded) {
                const status = self.finishUnbounded(problem, entering.column, entering.direction);
                if (status == .numerical_failure and self.active_degeneracy_strategy != .baseline) {
                    self.stats.cold_restart_solves += 1;
                    return self.restartSolveWithoutPerturbation(problem, control);
                }
                return status;
            }
            if (leaving.status != .optimal) return leaving.status;
            const leaving_column = basis.basic_index[@intCast(leaving.row.?)];
            const objective_change = @abs(basis.reduced_cost[entering.column] * leaving.step);
            const ratio_tie_count = if ((self.statistics_io != null or self.active_degeneracy_trace.len != 0) and
                leaving.step <= self.numerical.primal_tolerance)
                self.countPrimalRatioTies(leaving.step)
            else
                0;
            if (self.updateReducedCostsAfterPrimalPivot(problem, entering.column, entering.direction, leaving.row.?) != .optimal)
                return .numerical_failure;
            if (self.performPivot(problem, entering.column, entering.direction, leaving.row.?, leaving.bound, leaving.step) != .optimal) {
                if (self.failure_site == .none) self.failure_site = .pivot_update;
                return .numerical_failure;
            }
            self.observeIterationStep(
                leaving.step,
                entering.column,
                entering.direction,
                leaving_column,
                ratio_tie_count,
                own_step,
                objective_change,
                small_pivot_retry,
                0,
            );
        }
        if ((self.fresh_factorization_pivots_remaining != 0 or self.reduced_cost_update_count >= self.reduced_cost_refresh_period) and
            self.recomputeReducedCosts(problem) != .optimal)
        {
            self.failure_site = .reduced_cost;
            return .numerical_failure;
        }
    }
    return .iteration_limit;
}

/// Convert the user-space dual-feasibility tolerance to scaled model space.
///
/// Reduced costs include objective and column scaling, so pricing must use the
/// smallest structural column scale to avoid rejecting a violation that would
/// exceed tolerance after unscaling. Machine epsilon is the absolute floor.
pub fn scaledDualTolerance(self: *const SimplexEngine, problem: problem_module.ProblemView) f64 {
    const basis = if (self.basis) |*value| value else return self.numerical.dual_tolerance;
    var minimum_column_scale: f64 = 1.0;
    for (basis.column_scale[0..problem.num_cols]) |scale| minimum_column_scale = @min(minimum_column_scale, scale);
    return @max(std.math.floatEps(f64), self.numerical.dual_tolerance * self.objective_scale * minimum_column_scale);
}

/// Price primal entering candidates and charge the scan to pricing statistics.
///
/// This baseline dispatcher selects Bland, partial multiple pricing, or the
/// configured weighted full scan. It also records reduced-cost density for
/// later diagnostics. `column_count` restricts pricing to the active prefix,
/// excluding Phase-I artificials during ordinary Phase II.
pub fn choosePrimalEnteringTimed(self: *SimplexEngine, column_count: usize, tolerance: f64) ?pricing_module.EnteringChoice {
    const started = self.statisticsTimestamp();
    defer self.recordPricingElapsed(started);
    const basis = if (self.basis) |*value| value else return null;
    self.observePricingDensity(basis.reduced_cost[0..column_count]);
    self.stats.dense_pricing_dispatches += 1;
    return if (self.numerical.anti_cycling_active)
        self.pricing.choosePrimalEnteringBland(
            basis.reduced_cost[0..column_count],
            basis.col_status[0..column_count],
            tolerance,
        )
    else if (self.pricing.rule == .partial)
        self.pricing.choosePrimalEnteringMultiple(
            basis.reduced_cost[0..column_count],
            basis.col_status[0..column_count],
            basis.col_edge_weight[0..column_count],
            basis.flip_columns[0..column_count],
            tolerance,
        )
    else
        self.pricing.choosePrimalEnteringWeighted(
            basis.reduced_cost[0..column_count],
            basis.col_status[0..column_count],
            basis.col_edge_weight[0..column_count],
            tolerance,
        );
}

/// Price an entering column under the active perturbation/taboo policy.
///
/// When degeneracy handling is active, stable column ranks break score ties
/// and the optional taboo horizon suppresses recently rejected columns.
/// Callers must retry with `respect_taboo == false` before treating an empty
/// result as a feasibility certificate. Outside that mode this has the same
/// weighted/partial dispatch semantics as `choosePrimalEnteringTimed`.
pub fn choosePrimalEnteringWeightedTimed(
    self: *SimplexEngine,
    column_count: usize,
    tolerance: f64,
    respect_taboo: bool,
) ?pricing_module.EnteringChoice {
    const started = self.statisticsTimestamp();
    defer self.recordPricingElapsed(started);
    const basis = if (self.basis) |*value| value else return null;
    self.observePricingDensity(basis.reduced_cost[0..column_count]);
    self.stats.dense_pricing_dispatches += 1;
    return if (self.degeneracy.active)
        self.pricing.choosePrimalEnteringPerturbed(
            basis.reduced_cost[0..column_count],
            basis.col_status[0..column_count],
            basis.col_edge_weight[0..column_count],
            self.degeneracy.column_rank[0..column_count],
            if (respect_taboo and self.active_degeneracy_strategy == .perturbation_taboo)
                self.degeneracy.taboo_until[0..column_count]
            else
                &.{},
            self.iterations,
            tolerance,
        )
    else if (self.pricing.rule == .partial)
        self.pricing.choosePrimalEnteringMultiple(
            basis.reduced_cost[0..column_count],
            basis.col_status[0..column_count],
            basis.col_edge_weight[0..column_count],
            basis.flip_columns[0..column_count],
            tolerance,
        )
    else
        self.pricing.choosePrimalEnteringWeighted(
            basis.reduced_cost[0..column_count],
            basis.col_status[0..column_count],
            basis.col_edge_weight[0..column_count],
            tolerance,
        );
}

/// Deterministic cold-start selector. Automatic mode remains opt-in until
/// corpus timing establishes a benefit over the frozen primal baseline.
pub fn chooseColdPhaseOneStrategy(self: *const SimplexEngine, problem: problem_module.ProblemView) PhaseOneStrategy {
    const basis = if (self.basis) |*value| value else return .primal;
    var primal_count: usize = 0;
    var primal_sum: f64 = 0.0;
    for (basis.basic_value, basis.basic_lower, basis.basic_upper) |value, lower, upper| {
        const amount = @max(@max(lower - value, value - upper), 0.0);
        if (amount > self.numerical.primal_tolerance) {
            primal_count += 1;
            primal_sum += amount;
        }
    }
    var dual_count: usize = 0;
    var free_count: usize = 0;
    var boxed_count: usize = 0;
    for (0..problem.num_cols + problem.num_rows) |column| {
        const lower = basis.col_lower[column];
        const upper = basis.col_upper[column];
        if (!std.math.isFinite(lower) and !std.math.isFinite(upper)) free_count += 1;
        if (std.math.isFinite(lower) and std.math.isFinite(upper) and lower != upper) boxed_count += 1;
        const reduced = basis.reduced_cost[column];
        const infeasible = switch (basis.col_status[column]) {
            .at_lower => reduced < -self.numerical.dual_tolerance,
            .at_upper => reduced > self.numerical.dual_tolerance,
            .free, .superbasic => @abs(reduced) > self.numerical.dual_tolerance,
            else => false,
        };
        if (infeasible) dual_count += 1;
    }
    var coefficient_min = std.math.inf(f64);
    var coefficient_max: f64 = 0.0;
    for (0..problem.num_cols) |column| {
        const begin = problem.matrix.col_starts[column];
        const end = problem.matrix.col_starts[column + 1];
        for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
            if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
            const magnitude = @abs(coefficient * basis.row_scale[row.toUsize()] * basis.column_scale[column]);
            coefficient_min = @min(coefficient_min, magnitude);
            coefficient_max = @max(coefficient_max, magnitude);
        }
    }
    var cost_min = std.math.inf(f64);
    var cost_max: f64 = 0.0;
    for (problem.col_cost) |cost| {
        const magnitude = @abs(cost);
        if (magnitude == 0.0) continue;
        cost_min = @min(cost_min, magnitude);
        cost_max = @max(cost_max, magnitude);
    }
    const coefficient_range = if (coefficient_min < std.math.inf(f64)) coefficient_max / coefficient_min else 1.0;
    const cost_range = if (cost_min < std.math.inf(f64)) cost_max / cost_min else 1.0;
    // Dual Phase I is favoured for many violated logical rows, provided
    // free columns do not dominate its amplified infeasibility objective.
    // Wide sparse models amortize BTRAN row pricing particularly well;
    // highly scaled models retain the more defensive primal path.
    const width_score = if (problem.num_rows == 0) problem.num_cols else problem.num_cols / problem.num_rows;
    const scale_risk = (if (coefficient_range > 1e12) problem.num_rows else 0) +
        (if (cost_range > 1e12) problem.num_rows else 0);
    const risk_score = free_count * 8 + boxed_count + scale_risk;
    const benefit_score = primal_count * 4 + @min(dual_count, problem.num_cols) + width_score;
    if (primal_count > 0 and primal_sum > self.numerical.primal_tolerance and
        problem.num_cols >= problem.num_rows * 2 and benefit_score > risk_score)
        return .dual;
    return .primal;
}

/// Minimize the sum of artificial variables to obtain primal feasibility.
///
/// Phase I uses its own objective, pricing selection, and reduced-cost refresh
/// cadence while sharing the main pivot machinery. An apparent infeasibility
/// or failed artificial cleanup reached through perturbed ordering is never
/// published directly: the logical starting epoch is restored and rerun with
/// the baseline policy first. Success leaves artificials fixed at zero and
/// prepares the surviving original/logical basis for Phase II.
pub fn solvePhaseOne(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
    const phase_started = self.statisticsTimestamp();
    const iteration_started = self.iterations;
    const saved_pricing_rule = self.pricing.rule;
    const phase_two_refresh_period = self.reduced_cost_refresh_period;
    self.pricing.rule = switch (control.phase_one_pricing) {
        .inherit => if (self.active_primal_pricing_strategy == .partial) .partial else saved_pricing_rule,
        .dantzig => .dantzig,
        .devex => .devex,
        .steepest_edge => .steepest_edge,
    };
    defer self.pricing.rule = saved_pricing_rule;
    // Phase-I drift observations do not justify relaxing the original
    // objective's refresh cadence across the phase boundary.
    defer self.reduced_cost_refresh_period = phase_two_refresh_period;
    defer self.recordPhaseElapsed(.phase_one, phase_started);
    defer self.recordPhaseIterations(.phase_one, iteration_started);
    if (self.pricing.rule == .devex and self.active_devex_strategy == .framework)
        self.initializePrimalDevexFramework();
    if (self.recomputePhaseOneReducedCosts(problem) != .optimal) return .numerical_failure;
    while (self.iterations < control.max_iterations) : (self.iterations += 1) {
        if (self.beginIteration(problem, control, .phase_one)) |status| return status;
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        if (self.prepareDegeneracyPolicy(basis.basic_value.len, basis.col_status.len) != .optimal)
            return .numerical_failure;
        // General bounded Phase I contains artificial columns and is not
        // the standard-form setting required by Bland's proof. Retain the
        // numerically stable Harris candidate set throughout this phase.
        var entering = self.choosePrimalEnteringWeightedTimed(
            basis.reduced_cost.len,
            self.numerical.dual_tolerance,
            true,
        );
        // A taboo cache is a search-order hint, never a feasibility
        // certificate. If it masks every improving column, expire the
        // cache and repeat the exact same pricing scan without exclusions
        // before allowing Phase I to terminate.
        if (entering == null and self.active_degeneracy_strategy == .perturbation_taboo and
            self.degeneracy.active)
        {
            self.degeneracy.invalidateTaboo();
            self.stats.taboo_retries += 1;
            entering = self.choosePrimalEnteringWeightedTimed(
                basis.reduced_cost.len,
                self.numerical.dual_tolerance,
                false,
            );
        }
        const entering_choice = entering orelse {
            if (self.factorization.update_count != 0) {
                self.factorization.recordReinversion(.solve_residual);
                if (self.refactorizeBasis(problem, .solve_residual) != .optimal) return .numerical_failure;
                if (self.recomputePhaseOneReducedCosts(problem) != .optimal) return .numerical_failure;
                continue;
            }
            break;
        };
        if (self.computeDirection(problem, entering_choice.column) != .optimal) return .numerical_failure;
        if (entering_choice.direction < 0) {
            for (basis.pivot_direction) |*value| value.* = -value.*;
        }
        const leaving = self.chooseLeavingWithPolicy(false);
        const own_step = if (entering_choice.direction > 0)
            basis.col_upper[entering_choice.column] - basis.primal[entering_choice.column]
        else
            basis.primal[entering_choice.column] - basis.col_lower[entering_choice.column];
        if (std.math.isFinite(own_step) and (leaving.status == .unbounded or own_step < leaving.step)) {
            for (basis.basic_value, basis.pivot_direction) |*value, direction| value.* -= own_step * direction;
            basis.primal[entering_choice.column] += entering_choice.direction * own_step;
            basis.col_status[entering_choice.column] = if (entering_choice.direction > 0) .at_upper else .at_lower;
            self.observeIterationStep(
                own_step,
                entering_choice.column,
                entering_choice.direction,
                null,
                0,
                null,
                @abs(basis.reduced_cost[entering_choice.column] * own_step),
                false,
                0,
            );
        } else {
            if (leaving.status != .optimal) return leaving.status;
            const leaving_column = basis.basic_index[@intCast(leaving.row.?)];
            const objective_change = @abs(basis.reduced_cost[entering_choice.column] * leaving.step);
            const ratio_tie_count = if ((self.statistics_io != null or self.active_degeneracy_trace.len != 0) and
                leaving.step <= self.numerical.primal_tolerance)
                self.countPrimalRatioTies(leaving.step)
            else
                0;
            if (self.updateReducedCostsAfterPrimalPivot(problem, entering_choice.column, entering_choice.direction, leaving.row.?) != .optimal)
                return .numerical_failure;
            if (self.performPivot(problem, entering_choice.column, entering_choice.direction, leaving.row.?, leaving.bound, leaving.step) != .optimal)
                return .numerical_failure;
            self.observeIterationStep(
                leaving.step,
                entering_choice.column,
                entering_choice.direction,
                leaving_column,
                ratio_tie_count,
                own_step,
                objective_change,
                false,
                0,
            );
        }
        const refresh_period = if (self.active_adaptive_reprice and self.degeneracy.active and
            self.numerical.consecutive_degenerate_pivots >= 32)
            @min(self.reduced_cost_refresh_period, 4)
        else
            self.reduced_cost_refresh_period;
        if (self.fresh_factorization_pivots_remaining != 0 or self.reduced_cost_update_count >= refresh_period) {
            const refresh_status = if (self.active_adaptive_reprice)
                self.recomputePhaseOneReducedCostsWithDrift(problem)
            else
                self.recomputePhaseOneReducedCosts(problem);
            if (refresh_status != .optimal) return .numerical_failure;
        }
    }
    if (self.iterations >= control.max_iterations) return .iteration_limit;
    if (self.recomputeBasicValues(problem) != .optimal) {
        // A bounded perturbation may leave the terminal basis just beyond
        // an original bound. That invalidates the perturbed path, not the
        // model. Restore the logical epoch so infeasibility is proved by
        // the unperturbed Phase-I path instead of escaping as a numerical
        // failure.
        if (self.degeneracy.ever_active) {
            self.stats.cold_restart_phase_one += 1;
            return self.restartPhaseOneWithoutPerturbation(problem, control);
        }
        return .numerical_failure;
    }
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    const artificial_begin = problem.num_cols + problem.num_rows;
    var infeasibility: f64 = 0.0;
    for (basis.primal[artificial_begin..]) |value| infeasibility += @max(value, 0.0);
    if (infeasibility > self.numerical.primal_tolerance) {
        // Perturbed or taboo-guided pricing may change the path, but an
        // infeasibility conclusion must come from the unperturbed logical
        // epoch. Restore it transactionally and let baseline Phase I
        // either prove infeasibility or continue searching.
        if (self.degeneracy.ever_active) {
            self.stats.cold_restart_phase_one += 1;
            return self.restartPhaseOneWithoutPerturbation(problem, control);
        }
        return .infeasible;
    }

    if (self.cleanupArtificialBasis(problem) != .optimal) {
        if (self.degeneracy.ever_active) {
            self.stats.cold_restart_phase_one += 1;
            return self.restartPhaseOneWithoutPerturbation(problem, control);
        }
        return .numerical_failure;
    }
    @memset(basis.col_upper[artificial_begin..], 0.0);
    for (basis.col_status[artificial_begin..]) |*status| {
        if (status.* != .basic) status.* = .fixed;
    }
    for (basis.basic_index, 0..) |column, row| {
        if (column >= artificial_begin) basis.basic_upper[row] = 0.0;
    }
    self.phase1_needed = false;
    self.numerical.clearAntiCyclingFallback();
    if (self.degeneracy.ever_active) {
        self.degeneracy.clearAfterProgress();
        self.stats.perturbation_cleanups += 1;
        if (self.refactorizeBasis(problem, .cleanup) != .optimal or
            self.recomputeBasicValues(problem) != .optimal or
            self.recomputePhaseOneReducedCosts(problem) != .optimal)
        {
            self.stats.cold_restart_phase_one += 1;
            return self.restartPhaseOneWithoutPerturbation(problem, control);
        }
        var cleanup_infeasibility: f64 = 0.0;
        for (basis.primal[artificial_begin..]) |value| cleanup_infeasibility += @max(value, 0.0);
        if (cleanup_infeasibility > self.numerical.primal_tolerance) {
            self.stats.cold_restart_phase_one += 1;
            return self.restartPhaseOneWithoutPerturbation(problem, control);
        }
    }
    return .optimal;
}

/// Transactional fallback for a failed perturbation cleanup. Restore the
/// validated logical epoch and rerun Phase I with the baseline policy;
/// no perturbed ordering state or reduced cost can escape to publication.
pub fn restartPhaseOneWithoutPerturbation(
    self: *SimplexEngine,
    problem: problem_module.ProblemView,
    control: SolveControl,
) SolveStatus {
    self.degeneracy.resetSolve();
    self.active_degeneracy_strategy = .baseline;
    self.numerical.clearAntiCyclingFallback();
    self.direction_requires_reinversion = false;
    self.fresh_factorization_pivots_remaining = 0;
    self.reduced_cost_update_count = 0;
    if (self.initializeProblemStorage(problem) != .optimal) return .numerical_failure;
    self.factorization.factorizeIdentity(problem.num_rows) catch return .numerical_failure;
    self.dual_edge_weights_valid = true;
    self.observeFactorizationStability();
    if (self.initializeLogicalBasicValues(problem) != .optimal) return .numerical_failure;
    if (self.installArtificialPhaseOneBasis(problem) != .optimal) return .numerical_failure;
    if (!self.phase1_needed) return .optimal;
    if (self.refactorizeBasis(problem, .phase_one_setup) != .optimal) return .numerical_failure;
    self.numerical.markRefactorized();
    return self.solvePhaseOne(problem, control);
}

/// Pivot zero-valued artificial basics out whenever a stable original or
/// logical column is available. An artificial that remains basic denotes
/// a rank-redundant row and is fixed at zero until presolve can remove that
/// row explicitly.
pub fn cleanupArtificialBasis(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
    const cleanup_started = self.statisticsTimestamp();
    defer {
        const elapsed = self.elapsedSince(cleanup_started);
        self.stats.cleanup_ns = std.math.add(u64, self.stats.cleanup_ns, elapsed) catch std.math.maxInt(u64);
    }
    self.cleanup_active = true;
    defer self.cleanup_active = false;
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    const artificial_begin = problem.num_cols + problem.num_rows;

    for (0..problem.num_rows) |leaving_row| {
        if (basis.basic_index[leaving_row] < artificial_begin) continue;

        var entering_col: ?usize = null;
        var best_pivot = self.numerical.pivot_tolerance;
        for (0..artificial_begin) |column| {
            if (basis.col_status[column] == .basic) continue;
            if (self.computeDirection(problem, column) != .optimal) return .numerical_failure;
            const candidate = @abs(basis.pivot_direction[leaving_row]);
            if (candidate > best_pivot) {
                best_pivot = candidate;
                entering_col = column;
            }
        }

        const column = entering_col orelse continue;
        if (self.computeDirection(problem, column) != .optimal) return .numerical_failure;
        if (self.performPivot(problem, column, 1.0, leaving_row, .fixed, 0.0) != .optimal)
            return .numerical_failure;
    }
    return self.recomputeBasicValues(problem);
}

/// Recompute exact reduced costs for the artificial Phase-I objective.
///
/// Basic costs are one only for artificial basics. Solving `B^-T c_B`
/// produces the dual vector, after which structural columns are priced through
/// the selected row- or column-oriented sparse representation. Logical and
/// artificial reduced costs are written from their known singleton columns.
/// A successful call resets the incremental reduced-cost update counter.
pub fn recomputePhaseOneReducedCosts(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    const artificial_begin = problem.num_cols + problem.num_rows;
    for (basis.basic_index, 0..) |global_col, row| basis.dual[row] = if (global_col >= artificial_begin) 1.0 else 0.0;
    self.factorization.solveTranspose(basis.dual) catch return .numerical_failure;
    const row_pricing = self.selectRowPricing(basis.dual);
    const pricing_started = self.statisticsTimestamp();
    defer if (row_pricing)
        self.recordRowPricingElapsed(pricing_started)
    else
        self.recordPricingElapsed(pricing_started);

    // Compute c - A^T*y directly from CSC. The previous generic column
    // loop cleared and dotted a rows-long dense vector for every internal
    // column, turning sparse Phase-I pricing into O(columns * rows).
    // Structural columns use one sparse pass; logical and artificial
    // columns have one known nonzero and are written directly.
    if (row_pricing) {
        @memset(basis.reduced_cost[0..problem.num_cols], 0.0);
        for (basis.dual, 0..) |dual_value, row_index| {
            if (@abs(dual_value) <= self.numerical.zero_tolerance) continue;
            const scaled_dual = dual_value * basis.row_scale[row_index];
            for (self.pricing_row_view.rowColumns(row_index), self.pricing_row_view.rowValues(row_index)) |column_u32, coefficient| {
                const column: usize = @intCast(column_u32);
                basis.reduced_cost[column] -= scaled_dual * coefficient * basis.column_scale[column];
            }
        }
    } else {
        for (0..problem.num_cols) |column| {
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            var reduced: f64 = 0.0;
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
                if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                const row_index = row.toUsize();
                reduced -= basis.row_scale[row_index] * coefficient * basis.column_scale[column] * basis.dual[row_index];
            }
            basis.reduced_cost[column] = reduced;
        }
    }
    for (0..problem.num_rows) |row| {
        basis.reduced_cost[problem.num_cols + row] = -basis.dual[row];
        basis.reduced_cost[artificial_begin + row] = 1.0 - basis.artificial_sign[row] * basis.dual[row];
    }
    self.reduced_cost_update_count = 0;
    return .optimal;
}

/// Reprice Phase I exactly and adapt the next refresh interval from drift.
///
/// The incrementally maintained reduced costs are snapshotted before the exact
/// recomputation. Their maximum relative discrepancy is accumulated for
/// diagnostics; large drift shortens the refresh period, while consistently
/// small drift may relax it only up to the validated eight-pivot ceiling.
pub fn recomputePhaseOneReducedCostsWithDrift(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    const count = basis.reduced_cost.len;
    if (self.degeneracy.reduced_cost_snapshot.len < count) return .numerical_failure;
    @memcpy(self.degeneracy.reduced_cost_snapshot[0..count], basis.reduced_cost);
    if (self.recomputePhaseOneReducedCosts(problem) != .optimal) return .numerical_failure;
    var maximum_drift: f64 = 0.0;
    for (basis.reduced_cost, self.degeneracy.reduced_cost_snapshot[0..count]) |exact, updated| {
        maximum_drift = @max(maximum_drift, @abs(exact - updated) / @max(1.0, @abs(exact)));
    }
    self.maximum_reduced_cost_drift = @max(self.maximum_reduced_cost_drift, maximum_drift);
    self.exact_reprices += 1;
    if (maximum_drift > self.numerical.dual_tolerance * 10.0)
        self.reduced_cost_refresh_period = @max(self.reduced_cost_refresh_period / 2, 1)
    else if (maximum_drift < self.numerical.dual_tolerance * 0.1)
        // Eight is the validated fixed-policy ceiling. Adaptive mode may
        // refresh more often after drift or a long degenerate chain, but
        // never weakens the existing numerical safeguard.
        self.reduced_cost_refresh_period = @min(self.reduced_cost_refresh_period * 2, 8);
    return .optimal;
}

test "engine solves a standard-form LP with primal revised simplex" {
    const rows = [_]foundation.RowId{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) };
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &[_]f64{ -1, -1 },
        .col_lower = &[_]f64{ 0, 0 },
        .col_upper = &[_]f64{ std.math.inf(f64), std.math.inf(f64) },
        .row_lower = &[_]f64{-std.math.inf(f64)},
        .row_upper = &[_]f64{4},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &[_]usize{ 0, 1, 2 }, &rows, &[_]f64{ 1, 2 }),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    const status = engine.solveProblem(problem, .{});
    try std.testing.expectEqual(SolveStatus.optimal, status);
    try std.testing.expectApproxEqAbs(@as(f64, 4), engine.basis.?.primal[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0), engine.basis.?.primal[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, -4), engine.objective_value, 1e-12);
    const solution = engine.solutionView(problem, status).?;
    try std.testing.expectEqual(@as(usize, 2), solution.primal.len);
    try std.testing.expectApproxEqAbs(@as(f64, -4), solution.objective_value, 1e-12);
    try std.testing.expectEqual(engine.iterations, engine.stats.phase_two_iterations);
    try std.testing.expectEqual(@as(usize, 0), engine.stats.phase_one_iterations);
    try std.testing.expect(engine.stats.pricing_calls > 0);
    try std.testing.expectEqual(@as(u64, 0), engine.stats.pricing_ns);
    try std.testing.expectEqual(engine.factorization.stats.ftran_calls, engine.factorization.stats.dense_ftran_dispatches);
    try std.testing.expectEqual(engine.stats.rebuild_calls, engine.stats.classifiedRebuilds());
    try std.testing.expect(engine.requestedBytes() > 0);
}

test "Bland fallback resolves the Beale degenerate cycling example" {
    const rows = [_]foundation.RowId{
        foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(2),
        foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(0),
        foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1),
    };
    const problem = problem_module.ProblemView{
        .num_rows = 3,
        .num_cols = 4,
        .col_cost = &[_]f64{ 10, -57, -9, -24 },
        .col_lower = &[_]f64{ 0, 0, 0, 0 },
        .col_upper = &[_]f64{ std.math.inf(f64), std.math.inf(f64), std.math.inf(f64), std.math.inf(f64) },
        .row_lower = &[_]f64{ -std.math.inf(f64), -std.math.inf(f64), -std.math.inf(f64) },
        .row_upper = &[_]f64{ 0, 0, 1 },
        .matrix = matrix.CscView.initAssumeValid(
            3,
            4,
            &[_]usize{ 0, 3, 5, 7, 9 },
            &rows,
            &[_]f64{ 0.5, 0.5, 1, -5.5, -1.5, -2.5, -0.5, 9, 1 },
        ),
        .objective_sense = .maximize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.pricing.rule = .dantzig;
    engine.numerical.degenerate_pivot_limit = 1;
    var degeneracy_trace: [100]DegeneracyTraceEvent = undefined;

    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{
        .max_iterations = 100,
        .degeneracy_trace = &degeneracy_trace,
    }));
    try std.testing.expect(engine.numerical.anti_cycling_activations > 0);
    try std.testing.expect(engine.numerical.degenerate_pivot_count > 0);
    try std.testing.expectEqual(engine.numerical.degenerate_pivot_count, engine.stats.classifiedDegeneratePivots());
    try std.testing.expectEqual(engine.numerical.degenerate_pivot_count, engine.degeneracy_trace_count);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), engine.objective_value, 1e-9);
}

test "Phase I restores feasibility for a greater-equal row" {
    const rows = [_]foundation.RowId{foundation.RowId.fromUsizeAssumeValid(0)};
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{1},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{2},
        .row_upper = &[_]f64{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &[_]usize{ 0, 1 }, &rows, &[_]f64{1}),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{}));
    try std.testing.expectApproxEqAbs(@as(f64, 2), engine.basis.?.primal[0], 1e-12);
}

test "Phase I removes zero artificial basics after redundant equalities" {
    const rows = [_]foundation.RowId{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1) };
    const problem = problem_module.ProblemView{
        .num_rows = 2,
        .num_cols = 1,
        .col_cost = &[_]f64{0},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{ 1, 1 },
        .row_upper = &[_]f64{ 1, 1 },
        .matrix = matrix.CscView.initAssumeValid(2, 1, &[_]usize{ 0, 2 }, &rows, &[_]f64{ 1, 1 }),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{}));
    try std.testing.expectApproxEqAbs(@as(f64, 1), engine.basis.?.primal[0], 1e-12);
    const artificial_begin = problem.num_cols + problem.num_rows;
    for (engine.basis.?.basic_index) |column| try std.testing.expect(column < artificial_begin);
}

test "Phase I detects an infeasible LP" {
    const rows = [_]foundation.RowId{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1) };
    const problem = problem_module.ProblemView{
        .num_rows = 2,
        .num_cols = 1,
        .col_cost = &[_]f64{0},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{ 2, -std.math.inf(f64) },
        .row_upper = &[_]f64{ std.math.inf(f64), 1 },
        .matrix = matrix.CscView.initAssumeValid(2, 1, &[_]usize{ 0, 2 }, &rows, &[_]f64{ 1, 1 }),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.infeasible, engine.solveProblem(problem, .{}));
}

test "Phase I respects a zero iteration limit" {
    const rows = [_]foundation.RowId{foundation.RowId.fromUsizeAssumeValid(0)};
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{0},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{2},
        .row_upper = &[_]f64{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &[_]usize{ 0, 1 }, &rows, &[_]f64{1}),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.iteration_limit, engine.solveProblem(problem, .{ .max_iterations = 0 }));
}

test "primal Devex framework updates every nonbasic tableau weight" {
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 2, 2);
    engine.basis.?.initializeSlackBasis();
    engine.initializePrimalDevexFramework();
    engine.basis.?.pivot_direction[0] = 2;
    engine.basis.?.pivot_direction[1] = 1;
    engine.basis.?.tableau[0] = 2;
    engine.basis.?.tableau[1] = 3;
    engine.basis.?.tableau[2] = 1;

    try std.testing.expectEqual(SolveStatus.optimal, engine.updatePrimalDevexFramework(0, 0, 2));
    // Initial logical basics are outside the frozen reference set, so the
    // pivotal reference norm is the entering coordinate alone: 1 / 2^2.
    try std.testing.expectApproxEqAbs(@as(f64, 3.25), engine.basis.?.col_edge_weight[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1), engine.basis.?.col_edge_weight[2], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1), engine.basis.?.col_edge_weight[0], 1e-12);
    try std.testing.expectEqual(@as(usize, 1), engine.stats.devex_framework_updates);
}
