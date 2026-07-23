//! Dual revised-simplex loop, dual Phase I, bound flips, and warm-basis
//! repair for `SimplexEngine`.
//!
//! ## Responsibility
//!
//! Owns the dual Phase-I/Phase-II iteration loops, dual leaving-row
//! selection, the tableau-row FTRAN, aggregate bound flips, and dual-based
//! repair of warm imported bases.

const std = @import("std");
const basis_module = @import("basis.zig");
const basis_snapshot_module = @import("basis_snapshot.zig");
const factorization_module = @import("factorization.zig");
const pricing_module = @import("pricing.zig");
const ratio_module = @import("ratio_test.zig");
const numerical_module = @import("numerical.zig");
const highs_random = @import("highs_random.zig");
const dual_phase_one_module = @import("dual_phase_one.zig");
const crash_module = @import("crash.zig");
const degeneracy_module = @import("degeneracy.zig");
const pricing_workspace_module = @import("pricing_workspace.zig");
const problem_module = @import("problem.zig");
const solution_module = @import("solution.zig");
const foundation = @import("foundation");
const matrix = @import("matrix");
const SimplexEngine = @import("engine.zig").SimplexEngine;
const Algorithm = @import("engine.zig").Algorithm;
const CallbackAction = @import("engine.zig").CallbackAction;
const ProgressEventView = @import("engine.zig").ProgressEventView;
const SolveStatus = @import("engine.zig").SolveStatus;
const DualPhaseOneCandidateReason = @import("engine.zig").DualPhaseOneCandidateReason;
const DualPhaseOneFailureDiagnostic = @import("engine.zig").DualPhaseOneFailureDiagnostic;
const SolveControl = @import("engine.zig").SolveControl;
const RebuildReason = @import("engine.zig").RebuildReason;

const enable_highs_rhs_infeasibility_list = false;

/// Enter the ordinary dual-simplex Phase II loop.
///
/// The caller must already have installed and factorized a valid basis.
/// This entry point binds the persistent dual workspace directly to the
/// basis's original bounds and costs, marks the phase transition in the
/// dual-control state machine, and then delegates all iteration decisions to
/// `solveDualLoop`. Unlike `solveDualWithCostShifts`, it never constructs a
/// perturbed auxiliary objective.
pub fn solveDual(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    _ = self.dual_work.bindPhaseTwo(basis, problem.num_cols + problem.num_rows) catch return .numerical_failure;
    self.dual_control.enterPhase(.phase_two);
    return solveDualLoop(self, problem, control, false);
}

/// Build a dual-feasible cold start using the same high-level mechanism as
/// HiGHS: flip avoidable boxed infeasibilities and shift one-sided costs.
/// Free-variable infeasibilities remain the responsibility of dual Phase I.
pub fn solveDualWithCostShifts(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    const count = problem.num_cols + problem.num_rows;
    self.dual_work.beginCostEpoch(basis, count) catch return .numerical_failure;
    self.buildDualPhaseOneCosts(problem);
    const correction = rebuildDualState(self, problem, true, null, true) orelse return .numerical_failure;
    if (correction.free_infeasibility_count != 0) {
        self.shifted_dual_exit = .setup_free_infeasibility;
        return .not_implemented;
    }
    if (!self.classifyFeasibility(problem).dual) {
        self.shifted_dual_exit = .setup_not_dual_feasible;
        return .not_implemented;
    }
    self.algorithm = .dual_revised;
    self.shifted_dual_exit = .running;
    const status = solveDualLoop(self, problem, control, true);
    // Every completed attempt must leave an actionable diagnostic. Individual
    // CHUZC/cleanup exits set a more specific reason inside solveDualLoop;
    // this wrapper closes the remaining numerical/external-stop paths.
    if (self.shifted_dual_exit == .running) {
        self.shifted_dual_failure_site = self.failure_site;
        self.shifted_dual_exit = switch (status) {
            .work_limit, .time_limit, .iteration_limit, .interrupted => .phase_two_stopped,
            else => .phase_two_numerical_failure,
        };
    }
    return status;
}

/// `HEkkDual::rebuild` in upstream order. Refactorization is optional because
/// HiGHS also rebuilds from a still-valid inverse for phase decisions.
fn rebuildDualState(
    self: *SimplexEngine,
    problem: problem_module.ProblemView,
    use_work_costs: bool,
    refactor_reason: ?RebuildReason,
    correct_dual: bool,
) ?dual_phase_one_module.DualCorrection {
    if (refactor_reason) |reason| {
        if (self.refactorizeBasis(problem, reason) != .optimal) return null;
    } else {
        _ = self.dual_control.beginRebuild();
    }
    self.dual_control.has_fresh_rebuild = false;
    const dual_status = if (use_work_costs)
        self.recomputeReducedCostsFromWork(problem)
    else
        self.recomputeReducedCosts(problem);
    if (dual_status != .optimal) return null;
    const basis = if (self.basis) |*value| value else return null;
    const count = problem.num_cols + problem.num_rows;
    const correction = if (correct_dual)
        self.dual_work.correctDualInfeasibilities(
            basis,
            count,
            self.numerical.dual_tolerance,
            self.dual_control.force_phase_two,
        )
    else
        dual_phase_one_module.DualCorrection{};
    if (correct_dual) {
        self.dual_control.force_phase_two = false;
        self.dual_control.costs_shifted = self.dual_control.costs_shifted or correction.shift_count != 0;
    }
    // Boxed flips change nonbasic values, so primal recomputation must follow
    // correction. Repricing again would consume a different state than HiGHS:
    // correctDual has already updated every shifted nonbasic workDual exactly.
    if (self.recomputeBasicValuesUnchecked(problem) != .optimal) return null;
    if (!self.dual_work.recordRebuildState(basis, count, self.numerical.primal_tolerance)) return null;
    if (enable_highs_rhs_infeasibility_list and self.active_dual_initialization_strategy == .highs and
        !self.dual_work.createPrimalInfeasibilityList(basis, problem.num_rows, self.numerical.primal_tolerance)) return null;
    self.dual_control.finishRebuild();
    return correction;
}

/// Measure the unperturbed starting basis before selecting dual Phase I or II.
///
/// The three dual and three primal statistics are retained for diagnostics
/// and for the HiGHS-compatible phase decision. `force_phase_two` accepts the
/// same tiny squared-dual-error case as an effectively feasible basis, while
/// `near_optimal` identifies starts for which cleanup should be inexpensive.
/// Returns `false` only when no basis is installed.
fn assessInitialDualState(self: *SimplexEngine, problem: problem_module.ProblemView) bool {
    const basis = if (self.basis) |*value| value else return false;
    const count = problem.num_cols + problem.num_rows;
    const dual_stats = self.dual_work.dualInfeasibilityStats(
        basis,
        count,
        self.numerical.dual_tolerance,
    );
    self.dual_control.unperturbed_dual_infeasibility_count = dual_stats.count;
    self.dual_control.unperturbed_dual_infeasibility_max = dual_stats.maximum;
    self.dual_control.unperturbed_dual_infeasibility_sum = dual_stats.sum;
    self.dual_control.force_phase_two = dual_stats.maximum * dual_stats.maximum < self.numerical.dual_tolerance;

    var primal_count: usize = 0;
    var primal_max: f64 = 0.0;
    var primal_sum: f64 = 0.0;
    for (basis.basic_value, basis.basic_lower, basis.basic_upper) |value, lower, upper| {
        const violation = @max(@max(lower - value, value - upper), 0.0);
        if (violation >= self.numerical.primal_tolerance) primal_count += 1;
        primal_max = @max(primal_max, violation);
        primal_sum += violation;
    }
    self.dual_control.initial_primal_infeasibility_count = primal_count;
    self.dual_control.initial_primal_infeasibility_max = primal_max;
    self.dual_control.initial_primal_infeasibility_sum = primal_sum;
    const no_dual_infeasibilities = dual_stats.count == 0 or self.dual_control.force_phase_two;
    self.dual_control.near_optimal = no_dual_infeasibilities and primal_count < 1000 and primal_max < 1e-3;
    return true;
}

/// Run the serial revised dual-simplex Phase II state machine.
///
/// Each iteration performs CHUZR (including exact DSE validation), BTRAN and
/// row pricing, CHUZC/BFRT, optional aggregate bound flips, and the accepted
/// basis pivot. Rebuilds and repricing are deliberately centralized here so
/// all numerical-recovery exits observe one consistent basis/work state.
///
/// When `initial_use_work_costs` is true, reduced costs come from the
/// perturbed/shifted workspace until primal feasibility is reached. The loop
/// then restores the original objective, performs a cleanup reinversion, and
/// only reports a terminal status after classifying that original problem.
fn solveDualLoop(
    self: *SimplexEngine,
    problem: problem_module.ProblemView,
    control: SolveControl,
    initial_use_work_costs: bool,
) SolveStatus {
    var use_work_costs = initial_use_work_costs;
    if (use_work_costs) self.dual_control.costs_shifted = true;
    defer self.dual_control.costs_shifted = false;
    const phase_started = self.statisticsTimestamp();
    const iteration_started = self.iterations;
    const saved_pricing_rule = self.beginDualEdgeWeightPhase();
    defer self.pricing.rule = saved_pricing_rule;
    defer self.recordPhaseElapsed(.phase_two, phase_started);
    defer self.recordPhaseIterations(.phase_two, iteration_started);
    const initial_reprice = if (use_work_costs)
        self.recomputeReducedCostsFromWork(problem)
    else
        self.recomputeReducedCosts(problem);
    if (initial_reprice != .optimal) {
        if (use_work_costs) self.shifted_dual_exit = .phase_two_numerical_failure;
        return .numerical_failure;
    }
    const shifted_iteration_started = self.iterations;
    var shifted_iterations_recorded = false;
    var no_entering_recovery_used = false;
    var fresh_no_entering_basis = false;
    defer if (use_work_costs and !shifted_iterations_recorded) {
        self.stats.shifted_dual_iterations += self.iterations -| shifted_iteration_started;
    };
    while (self.iterations < control.max_iterations) : (self.iterations += 1) {
        if (self.beginIteration(problem, control, .phase_two)) |status| return status;
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        const work = self.dual_work.view(basis, problem.num_cols + problem.num_rows) orelse return .numerical_failure;
        if (self.pricing.rule == .steepest_edge and self.ensureExactDualEdgeWeights() != .optimal) {
            self.failure_site = .pivot_edge_weights;
            return .numerical_failure;
        }
        var leaving_choice: ?pricing_module.DualLeavingChoice = null;
        if (chooseDualLeavingWithExactDseValidation(self, problem, &leaving_choice) != .optimal)
            return .numerical_failure;
        const leaving = leaving_choice orelse {
            if (!use_work_costs) return self.finishOptimal(problem);
            // Shifted costs are only a path to primal feasibility. Never
            // publish their objective: restore the original costs and perform
            // a fresh factorization, basic-value rebuild and classification
            // first. Reusing the updated factor chain here made brandy's
            // original-cost cleanup fail after an otherwise successful dual
            // pass and triggered a full cold fallback.
            if (rebuildDualState(self, problem, false, .cleanup, false) == null) {
                self.shifted_dual_exit = .phase_two_numerical_failure;
                return .numerical_failure;
            }
            self.stats.shifted_dual_iterations += self.iterations -| shifted_iteration_started;
            shifted_iterations_recorded = true;
            const original_feasibility = self.classifyFeasibility(problem);
            if (original_feasibility.dual) {
                if (original_feasibility.primal) {
                    self.shifted_dual_exit = .original_dual_feasible;
                    return self.finishOptimal(problem);
                }
                // Cost cleanup can preserve original-objective dual
                // feasibility while exposing a small primal infeasibility.
                // This is a valid original-cost dual Phase-II start, not an
                // optimum. Calling finishOptimal here caused brandy to fail
                // validation and cold-restart the entire solve.
                // Continue original-cost Phase II in this controller. Calling
                // solveDual recursively used to reset phase-local edge-weight,
                // accounting and rebuild state, unlike HEkkDual's single major
                // phase loop.
                self.stats.shifted_dual_iterations += self.iterations -| shifted_iteration_started;
                shifted_iterations_recorded = true;
                use_work_costs = false;
                self.dual_control.phase = .phase_two;
                self.dual_control.costs_perturbed = false;
                self.dual_control.costs_shifted = false;
                self.dual_control.has_fresh_rebuild = true;
                self.shifted_dual_exit = .original_dual_feasible;
                no_entering_recovery_used = false;
                fresh_no_entering_basis = true;
                continue;
            }
            if (original_feasibility.primal) {
                // Primal feasibility does not authorize an algorithm switch
                // on an explicitly forced-dual run. The original costs are
                // dual infeasible, so hand the restored basis to dual Phase I
                // just as HiGHS does instead of solving it with primal Phase
                // II and reporting that iteration count as a dual comparison.
                self.shifted_dual_exit = .original_dual_infeasible;
                return .not_implemented;
            }
            self.shifted_dual_exit = .cleanup_neither_feasible;
            return .numerical_failure;
        };

        if (self.active_dual_edge_weight_strategy == .steepest_devex and
            self.pricing.rule == .steepest_edge and !self.dual_edge_weights_valid)
            self.switchDualDseToDevex(.invalid);
        const original_cols = problem.num_cols + problem.num_rows;
        var active_ratio_test = self.ratio_test;
        if (self.numerical.anti_cycling_active) active_ratio_test.rule = .standard;
        const entering = active_ratio_test.chooseDualEntering(
            basis.tableau[0..original_cols],
            work.dual,
            basis.col_status[0..original_cols],
            work.move,
            work.lower,
            work.upper,
            work.value,
            leaving.bound,
            leaving.violation,
            basis.dual_ratio[0..original_cols],
            basis.dual_direction[0..original_cols],
            basis.flip_columns[0..original_cols],
            .{
                .highs_choose_possible = self.active_dual_initialization_strategy == .highs,
                .highs_large_step = self.active_dual_initialization_strategy == .highs,
                .highs_breakpoint_groups = enable_highs_rhs_infeasibility_list,
                .trace_chuzc = self.active_pivot_trace.len != 0,
                .permutation = if (self.active_dual_initialization_strategy == .highs and
                    self.dual_work.highs_permutation.len >= original_cols)
                    self.dual_work.highs_permutation[0..original_cols]
                else
                    &.{},
                .update_count = self.factorization.update_count,
                .dual_feasibility_tolerance = self.numerical.dual_tolerance,
            },
        );

        if (self.applyBoundFlips(problem, entering.flip_count) != .optimal) {
            self.failure_site = .pivot_update;
            return .numerical_failure;
        }
        const entering_col = entering.column orelse {
            if (use_work_costs) {
                if (!no_entering_recovery_used) {
                    no_entering_recovery_used = true;
                    if (rebuildDualState(self, problem, true, .solve_residual, true) != null and
                        self.classifyFeasibility(problem).dual)
                    {
                        fresh_no_entering_basis = true;
                        continue;
                    }
                }
                if (fresh_no_entering_basis and self.buildAndValidateInfeasibilityRay(problem, leaving)) {
                    self.shifted_dual_exit = .phase_two_fresh_infeasible;
                    return .infeasible;
                }
                self.shifted_dual_exit = .phase_two_no_entering;
                // Reuse the allocation-free candidate diagnostic for the
                // shifted Phase-II row. Capturing it before Phase I starts is
                // essential: otherwise the later fallback overwrites the
                // actual cold-start CHUZC failure site.
                self.recordDualPhaseOneNoEntering(leaving, entering.flip_count, original_cols);
            }
            return .infeasible;
        };
        const entering_index: usize = @intCast(entering_col);
        fresh_no_entering_basis = false;
        if (self.computeDirection(problem, entering_index) != .optimal) return .numerical_failure;
        if (entering.direction < 0.0) {
            for (basis.pivot_direction) |*value| value.* = -value.*;
        }
        const leaving_row: usize = @intCast(leaving.row);
        const target = if (leaving.bound == .at_lower) work.lower[basis.basic_index[leaving_row]] else work.upper[basis.basic_index[leaving_row]];
        const pivot = basis.pivot_direction[leaving_row];
        if (@abs(pivot) <= self.numerical.pivot_tolerance) {
            self.failure_site = .ratio_test;
            return .numerical_failure;
        }
        const step = (basis.basic_value[leaving_row] - target) / pivot;
        if (!std.math.isFinite(step) or step < -self.numerical.primal_tolerance) {
            self.failure_site = .primal_step;
            return .numerical_failure;
        }
        const leaving_column = basis.basic_index[leaving_row];
        if (self.updateReducedCostsAfterDualPivot(
            problem,
            entering_index,
            leaving.row,
            entering.zero_dual_step,
        ) != .optimal) {
            self.failure_site = .reduced_cost;
            return .numerical_failure;
        }
        if (self.performPivot(
            problem,
            entering_index,
            entering.direction,
            leaving_row,
            leaving.bound,
            @max(step, 0.0),
        ) != .optimal) return .numerical_failure;
        // `performPivot` may have reinverted after committing the basis. A
        // fresh factorization is only the first half of dual rebuild: refresh
        // all work solution state before another CHUZR decision.
        if (self.dual_control.has_fresh_rebuild) {
            const rebuild_dual = if (use_work_costs)
                self.recomputeReducedCostsFromWork(problem)
            else
                self.recomputeReducedCosts(problem);
            if (rebuild_dual != .optimal) return .numerical_failure;
        }
        self.observeIterationStep(
            @max(step, 0.0),
            entering_index,
            entering.direction,
            leaving_column,
            0,
            null,
            null,
            false,
            entering.flip_count,
        );
        if (self.fresh_factorization_pivots_remaining != 0 or
            self.reduced_cost_update_count >= self.reduced_cost_refresh_period)
        {
            const refresh = if (use_work_costs)
                self.recomputeReducedCostsFromWork(problem)
            else
                self.recomputeReducedCostsWithDrift(problem);
            if (refresh != .optimal) {
                self.failure_site = .reduced_cost;
                return .numerical_failure;
            }
            if (!self.classifyFeasibility(problem).dual) {
                self.failure_site = .dual_feasibility;
                return .numerical_failure;
            }
        }
    }
    return .iteration_limit;
}

/// Construct the HiGHS-style row-space proof for a fresh Phase-II
/// no-entering decision and validate it entirely in original coordinates.
/// Costs do not participate: this certifies primal infeasibility of the
/// original bounds, not of the shifted objective.
pub fn buildAndValidateInfeasibilityRay(
    self: *SimplexEngine,
    problem: problem_module.ProblemView,
    leaving: pricing_module.DualLeavingChoice,
) bool {
    self.infeasibility_certificate_failure = .none;
    self.infeasibility_certificate_infinite_row_mass = 0.0;
    self.infeasibility_certificate_infinite_column_mass = 0.0;
    const basis = if (self.basis) |*value| value else {
        self.infeasibility_certificate_failure = .invalid_workspace;
        return false;
    };
    if (basis.infeasibility_ray.len != problem.num_rows or basis.residual_work.len != problem.num_rows) {
        self.infeasibility_certificate_failure = .invalid_workspace;
        return false;
    }
    // `residual_work` is dead after the fresh BTRAN at this decision point.
    // Reuse it allocation-free to measure each ray component in the same
    // relative form as HiGHS: |y_i| * max_j |A_ij|.
    @memset(basis.residual_work, 0.0);
    for (0..problem.num_cols) |column| {
        const begin = problem.matrix.col_starts[column];
        const end = problem.matrix.col_starts[column + 1];
        for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, value| {
            const row_index = row.toUsize();
            basis.residual_work[row_index] = @max(basis.residual_work[row_index], @abs(value));
        }
    }
    const move_out: f64 = if (leaving.bound == .at_lower) -1.0 else 1.0;
    var proof_lower: f64 = 0.0;
    var ray_max: f64 = 0.0;
    for (basis.infeasibility_ray, basis.dual_row, basis.row_scale, basis.residual_work, problem.row_lower, problem.row_upper) |
        *ray_value,
        scaled_ep,
        row_scale,
        row_max,
        lower,
        upper,
    | {
        var value = move_out * scaled_ep * row_scale;
        if (!std.math.isFinite(value)) {
            self.infeasibility_certificate_failure = .nonfinite_ray;
            return false;
        }
        const relative_contribution = @abs(value) * row_max;
        // Match the model's canonical coefficient floor before validating the
        // proof. Any remaining ray is still re-evaluated against the original
        // CSC and original bounds below; this is not a feasibility-tolerance
        // relaxation.
        if (relative_contribution <= matrix.MatrixTargetPolicy.model_coefficient_tolerance) value = 0.0;
        if (value > 0.0) {
            if (!std.math.isFinite(lower)) {
                self.infeasibility_certificate_infinite_row_mass = @max(
                    self.infeasibility_certificate_infinite_row_mass,
                    relative_contribution,
                );
                value = 0.0;
            } else {
                proof_lower += value * lower;
            }
        } else if (value < 0.0) {
            if (!std.math.isFinite(upper)) {
                self.infeasibility_certificate_infinite_row_mass = @max(
                    self.infeasibility_certificate_infinite_row_mass,
                    relative_contribution,
                );
                value = 0.0;
            } else {
                proof_lower += value * upper;
            }
        }
        ray_value.* = value;
        ray_max = @max(ray_max, @abs(value));
    }
    if (ray_max <= self.numerical.zero_tolerance) {
        self.infeasibility_certificate_failure = .zero_ray;
        return false;
    }
    if (!std.math.isFinite(proof_lower)) {
        self.infeasibility_certificate_failure = .nonfinite_proof;
        return false;
    }

    var implied_upper: f64 = 0.0;
    var has_infinite_column_bound = false;
    for (0..problem.num_cols) |column| {
        const begin = problem.matrix.col_starts[column];
        const end = problem.matrix.col_starts[column + 1];
        var coefficient: f64 = 0.0;
        for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, value| {
            coefficient += basis.infeasibility_ray[row.toUsize()] * value;
        }
        if (!std.math.isFinite(coefficient)) {
            self.infeasibility_certificate_failure = .nonfinite_proof;
            return false;
        }
        if (coefficient > 0.0) {
            const upper = problem.col_upper[column];
            if (!std.math.isFinite(upper)) {
                has_infinite_column_bound = true;
                self.infeasibility_certificate_infinite_column_mass += @abs(coefficient);
                continue;
            }
            implied_upper += coefficient * upper;
        } else if (coefficient < 0.0) {
            const lower = problem.col_lower[column];
            if (!std.math.isFinite(lower)) {
                has_infinite_column_bound = true;
                self.infeasibility_certificate_infinite_column_mass += @abs(coefficient);
                continue;
            }
            implied_upper += coefficient * lower;
        }
    }
    if (has_infinite_column_bound and
        self.infeasibility_certificate_infinite_column_mass > matrix.MatrixTargetPolicy.model_coefficient_tolerance)
    {
        self.infeasibility_certificate_failure = .infinite_column_bound;
        return false;
    }
    if (!std.math.isFinite(implied_upper)) {
        self.infeasibility_certificate_failure = .nonfinite_proof;
        return false;
    }
    const gap = proof_lower - implied_upper;
    self.infeasibility_certificate_gap = gap;
    if (!std.math.isFinite(gap)) {
        self.infeasibility_certificate_failure = .nonfinite_proof;
        return false;
    }
    if (gap <= self.numerical.primal_tolerance) {
        self.infeasibility_certificate_failure = .nonpositive_gap;
        return false;
    }
    self.infeasibility_ray_valid = true;
    return true;
}

/// Dedicated dual Phase I. Special bounds encode the negated dual
/// infeasibility objective, while deterministic perturbed costs remove
/// ratio ties. The borrowed matrix is never copied and every iteration
/// uses engine-owned work arrays.
pub fn solveDualPhaseOne(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
    self.dual_control.enterPhase(.phase_one);
    const phase_started = self.statisticsTimestamp();
    const iteration_started = self.iterations;
    const saved_pricing_rule = self.beginDualEdgeWeightPhase();
    defer self.pricing.rule = saved_pricing_rule;
    defer self.recordPhaseElapsed(.dual_phase_one, phase_started);
    defer self.recordPhaseIterations(.dual_phase_one, iteration_started);
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    const original_count = problem.num_cols + problem.num_rows;
    self.dual_work.begin(basis, original_count) catch return .numerical_failure;
    defer if (self.dual_work.active)
        self.dual_work.restoreOriginalBounds(basis, original_count);

    initializeUnperturbedDualCosts(self, problem);
    if (self.recomputeReducedCostsFromWork(problem) != .optimal) return .numerical_failure;
    if (!assessInitialDualState(self, problem)) return .numerical_failure;
    self.buildDualPhaseOneCosts(problem);
    self.dual_control.costs_perturbed = true;
    if (self.recomputeReducedCostsFromWork(problem) != .optimal) return .not_implemented;
    // Compatibility path retained until the complete B2 phase-transition
    // lifecycle is enabled; the isolated source-equivalent no-flip/lower-side
    // experiment is recorded in todo.md.
    flipDualInfeasibleBoxedColumns(
        basis,
        self.dual_work.nonbasic_move[0..original_count],
        original_count,
        self.numerical.dual_tolerance,
    );
    // Infeasibility is encoded in nonbasic primal values (±1/0), not
    // shifted costs. The raw scaled LP cost + perturbation keeps the
    // reduced cost magnitude intact so the Phase-I objective correctly
    // measures remaining infeasibility.
    self.dual_work.installWorkingBounds(basis, original_count);
    if (self.active_dual_initialization_strategy == .highs) {
        self.dual_phase_one_initial_flips = self.dual_work.correctInitialDualInfeasibilities(
            basis,
            original_count,
            self.numerical.dual_tolerance,
        );
        self.dual_phase_one_initial_objective = self.dual_work.dualObjective(basis, original_count);
    }
    if (self.recomputeBasicValuesUnchecked(problem) != .optimal) return .not_implemented;
    if (enable_highs_rhs_infeasibility_list and self.active_dual_initialization_strategy == .highs and
        !self.dual_work.createPrimalInfeasibilityList(basis, problem.num_rows, self.numerical.primal_tolerance))
        return .not_implemented;
    for (basis.basic_value, basis.basic_lower, basis.basic_upper) |value, lower, upper| {
        self.dual_phase_one_initial_max_basic_violation = @max(
            self.dual_phase_one_initial_max_basic_violation,
            @max(lower - value, value - upper),
        );
    }
    self.algorithm = .dual_revised;
    self.dual_edge_weights_valid = false;

    while (self.iterations < control.max_iterations) : (self.iterations += 1) {
        if (self.beginIteration(problem, control, .dual_phase_one)) |status| return status;
        if (self.pricing.rule == .steepest_edge and self.ensureExactDualEdgeWeights() != .optimal)
            return .not_implemented;
        var leaving_choice: ?pricing_module.DualLeavingChoice = null;
        if (chooseDualLeavingWithExactDseValidation(self, problem, &leaving_choice) != .optimal)
            return .not_implemented;
        const leaving = leaving_choice orelse break;
        if (self.dual_phase_one_initial_leaving_row == null)
            self.dual_phase_one_initial_leaving_row = leaving.row;
        if (self.active_dual_edge_weight_strategy == .steepest_devex and
            self.pricing.rule == .steepest_edge and !self.dual_edge_weights_valid)
            self.switchDualDseToDevex(.invalid);
        var active_ratio_test = self.ratio_test;
        if (self.numerical.anti_cycling_active) active_ratio_test.rule = .standard;
        const work = self.dual_work.view(basis, original_count) orelse return .not_implemented;
        const entering = active_ratio_test.chooseDualEntering(
            basis.tableau[0..original_count],
            work.dual,
            basis.col_status[0..original_count],
            work.move,
            work.lower,
            work.upper,
            work.value,
            leaving.bound,
            leaving.violation,
            basis.dual_ratio[0..original_count],
            basis.dual_direction[0..original_count],
            basis.flip_columns[0..original_count],
            .{
                .highs_choose_possible = self.active_dual_initialization_strategy == .highs,
                .highs_large_step = self.active_dual_initialization_strategy == .highs,
                .highs_breakpoint_groups = enable_highs_rhs_infeasibility_list,
                .trace_chuzc = self.active_pivot_trace.len != 0,
                .permutation = if (self.dual_work.highs_permutation.len >= original_count)
                    self.dual_work.highs_permutation[0..original_count]
                else
                    &.{},
                .update_count = self.factorization.update_count,
                .dual_feasibility_tolerance = self.numerical.dual_tolerance,
            },
        );
        self.dual_work.recordRemainingViolation(
            @intCast(leaving.row),
            leaving.violation,
            basis.tableau[0..original_count],
            basis.flip_columns[0..entering.flip_count],
            basis.col_lower[0..original_count],
            basis.col_upper[0..original_count],
        );
        if (entering.column == null)
            self.recordDualPhaseOneNoEntering(leaving, entering.flip_count, original_count);
        if (self.applyBoundFlips(problem, entering.flip_count) != .optimal) return .not_implemented;
        // Flips were applied but no entering column survived the ratio test:
        // continue the loop so the updated bounds can produce a pivot candidate
        // on the next pass. This is the small-flip-capacity path (HiGHS-style
        // [0,1] bounds give ~1 capacity per flip vs dynamic radius ~720).
        if (entering.column == null) {
            if (entering.flip_count == 0) return .not_implemented;
            if (self.recomputeBasicValuesUnchecked(problem) != .optimal) return .not_implemented;
            continue;
        }
        const entering_col = entering.column orelse unreachable;
        const entering_index: usize = @intCast(entering_col);
        if (self.computeDirection(problem, entering_index) != .optimal) return .not_implemented;
        if (entering.direction < 0.0) {
            for (basis.pivot_direction) |*value| value.* = -value.*;
        }
        const leaving_row: usize = @intCast(leaving.row);
        const target = if (leaving.bound == .at_lower) basis.basic_lower[leaving_row] else basis.basic_upper[leaving_row];
        const pivot = basis.pivot_direction[leaving_row];
        if (@abs(pivot) <= self.numerical.pivot_tolerance) return .not_implemented;
        const step = (basis.basic_value[leaving_row] - target) / pivot;
        if (!std.math.isFinite(step) or step < -self.numerical.primal_tolerance) return .not_implemented;
        const leaving_column = basis.basic_index[leaving_row];
        if (self.updateReducedCostsAfterDualPivot(
            problem,
            entering_index,
            leaving.row,
            entering.zero_dual_step,
        ) != .optimal) return .not_implemented;
        if (self.performPivot(problem, entering_index, entering.direction, leaving_row, leaving.bound, @max(step, 0.0)) != .optimal)
            return .not_implemented;
        if (self.dual_control.has_fresh_rebuild and self.recomputeReducedCostsFromWork(problem) != .optimal)
            return .not_implemented;
        self.observeIterationStep(
            @max(step, 0.0),
            entering_index,
            entering.direction,
            leaving_column,
            0,
            null,
            null,
            false,
            entering.flip_count,
        );
    }
    if (self.iterations >= control.max_iterations) return .iteration_limit;

    // Assess the Phase-I objective while the special bounds and values are
    // still installed. HiGHS then reinstalls the original (Phase-II) bounds
    // before performing any Phase-II iteration.
    if (self.recomputeReducedCostsFromWork(problem) != .optimal) return .not_implemented;

    // Always restore the original model state before classifying feasibility
    // or entering either Phase-II algorithm. Boxed columns are placed on the
    // side implied by the freshly priced reduced cost during restoration.
    self.dual_work.restoreOriginalBounds(basis, original_count);
    if (self.refactorizeBasis(problem, .cleanup) != .optimal) return .not_implemented;
    if (self.recomputeBasicValuesUnchecked(problem) != .optimal) return .not_implemented;
    if (self.recomputeReducedCosts(problem) != .optimal) return .not_implemented;

    var dual_infeasible: usize = 0;
    for (basis.reduced_cost[0..original_count], basis.col_status[0..original_count]) |reduced, status| {
        const infeasible = switch (status) {
            .at_lower => reduced < -self.numerical.dual_tolerance,
            .at_upper => reduced > self.numerical.dual_tolerance,
            .free, .superbasic => @abs(reduced) > self.numerical.dual_tolerance,
            .basic, .fixed => false,
        };
        if (infeasible) dual_infeasible += 1;
    }
    if (dual_infeasible == 0) return self.solveDual(problem, control);

    var primal_feasible = true;
    for (basis.basic_value, basis.basic_lower, basis.basic_upper) |value, lower, upper| {
        if (value < lower - self.numerical.primal_tolerance or value > upper + self.numerical.primal_tolerance) {
            primal_feasible = false;
            break;
        }
    }
    if (primal_feasible) {
        self.pricing.rule = saved_pricing_rule;
        return self.solvePrimal(problem, control);
    }
    return .not_implemented;
}

/// Repair dual-infeasible finite boxed nonbasic columns by moving each one to
/// its opposite bound.
///
/// This is the allocation-free initialization step used before shifted-cost
/// dual Phase II. It updates the public basis status and primal value together
/// with the explicit HiGHS-style nonbasic move (`+1` from lower, `-1` from
/// upper). Free, fixed, basic, and one-sided columns are intentionally left
/// unchanged because a simple bound exchange cannot repair them.
fn flipDualInfeasibleBoxedColumns(
    basis: *basis_module.BasisState,
    nonbasic_move: []i8,
    column_count: usize,
    dual_tolerance: f64,
) void {
    std.debug.assert(nonbasic_move.len >= column_count);
    for (0..column_count) |column| {
        const lower = basis.col_lower[column];
        const upper = basis.col_upper[column];
        if (!std.math.isFinite(lower) or !std.math.isFinite(upper) or lower == upper) continue;
        const reduced = basis.reduced_cost[column];
        switch (basis.col_status[column]) {
            .at_lower => if (reduced < -dual_tolerance) {
                basis.col_status[column] = .at_upper;
                basis.primal[column] = upper;
                nonbasic_move[column] = -1;
            },
            .at_upper => if (reduced > dual_tolerance) {
                basis.col_status[column] = .at_lower;
                basis.primal[column] = lower;
                nonbasic_move[column] = 1;
            },
            else => {},
        }
    }
}

/// Capture the first dual Phase-I row for which the ratio test cannot
/// select an entering column. This post-decision scan is diagnostic only:
/// it reads the already materialized ep/tableau and caller-owned buffers.
pub fn recordDualPhaseOneNoEntering(
    self: *SimplexEngine,
    leaving: pricing_module.DualLeavingChoice,
    flip_count: usize,
    column_count: usize,
) void {
    if (self.dual_phase_one_failure != null) return;
    const basis = if (self.basis) |*value| value else return;
    var diagnostic = DualPhaseOneFailureDiagnostic{
        .iteration = self.iterations,
        .leaving_row = leaving.row,
        .leaving_column = basis.basic_index[leaving.row],
        .leaving_bound = leaving.bound,
        .violation = leaving.violation,
        .leaving_value = basis.basic_value[leaving.row],
        .working_lower = basis.basic_lower[leaving.row],
        .working_upper = basis.basic_upper[leaving.row],
        .original_lower = self.dual_work.saved_lower[basis.basic_index[leaving.row]],
        .original_upper = self.dual_work.saved_upper[basis.basic_index[leaving.row]],
        .ep_nonzeros = 0,
        .ep_max_abs = 0.0,
        .small_tableau = 0,
        .basic_or_fixed = 0,
        .wrong_pivot_sign = 0,
        .accepted_bound_flips = 0,
        .eligible_unselected = 0,
    };
    for (basis.dual_row, 0..) |value, row| {
        diagnostic.ep_max_abs = @max(diagnostic.ep_max_abs, @abs(value));
        if (@abs(value) <= self.numerical.zero_tolerance) continue;
        diagnostic.ep_nonzeros += 1;
        if (self.dual_phase_one_ep_trace_count < self.active_dual_phase_one_ep_trace.len) {
            self.active_dual_phase_one_ep_trace[self.dual_phase_one_ep_trace_count] = .{
                .row = @intCast(row),
                .value = value,
            };
            self.dual_phase_one_ep_trace_count += 1;
        }
    }

    for (0..column_count) |column| {
        const alpha = basis.tableau[column];
        const status = basis.col_status[column];
        var direction: f64 = 0.0;
        var signed_pivot: f64 = 0.0;
        const width = basis.col_upper[column] - basis.col_lower[column];
        const flip_capacity = if (std.math.isFinite(width)) @abs(alpha) * width else 0.0;
        const reason: DualPhaseOneCandidateReason = blk: {
            if (@abs(alpha) <= self.ratio_test.tolerance) {
                diagnostic.small_tableau += 1;
                break :blk .small_tableau;
            }
            direction = switch (status) {
                .at_lower => 1.0,
                .at_upper => -1.0,
                .free, .superbasic => if (leaving.bound == .at_lower)
                    (if (alpha < 0.0) 1.0 else -1.0)
                else
                    (if (alpha > 0.0) 1.0 else -1.0),
                .basic, .fixed => {
                    diagnostic.basic_or_fixed += 1;
                    break :blk .basic_or_fixed;
                },
            };
            signed_pivot = alpha * direction;
            const eligible = if (leaving.bound == .at_lower)
                signed_pivot < -self.ratio_test.tolerance
            else
                signed_pivot > self.ratio_test.tolerance;
            if (!eligible) {
                diagnostic.wrong_pivot_sign += 1;
                break :blk .wrong_pivot_sign;
            }
            var was_flipped = false;
            for (basis.flip_columns[0..flip_count]) |flipped| {
                if (@as(usize, @intCast(flipped)) == column) {
                    was_flipped = true;
                    break;
                }
            }
            if (was_flipped) {
                diagnostic.accepted_bound_flips += 1;
                break :blk .accepted_bound_flip;
            }
            diagnostic.eligible_unselected += 1;
            break :blk .eligible_unselected;
        };
        if (self.dual_phase_one_candidate_trace_count < self.active_dual_phase_one_candidate_trace.len) {
            self.active_dual_phase_one_candidate_trace[self.dual_phase_one_candidate_trace_count] = .{
                .column = @intCast(column),
                .status = status,
                .explicit_move = self.dual_work.nonbasic_move[column],
                .tableau = alpha,
                .direction = direction,
                .signed_pivot = signed_pivot,
                .reduced_cost = basis.reduced_cost[column],
                .lower = basis.col_lower[column],
                .upper = basis.col_upper[column],
                .primal = basis.primal[column],
                .flip_capacity = flip_capacity,
                .reason = reason,
            };
            self.dual_phase_one_candidate_trace_count += 1;
        }
    }
    self.dual_phase_one_failure = diagnostic;
}

/// Construct the persistent perturbed working objective used by dual startup.
///
/// Structural costs are first converted to the engine's minimization sign and
/// scaling convention. A bound-aware perturbation is then applied with either
/// the HiGHS random stream/permutation or the deterministic comparison stream.
/// Logical columns receive only a tiny perturbation. Original costs and bounds
/// remain available in `DualPhaseOneWorkspace`, so cleanup can restore the
/// exact LP objective without reallocating or reconstructing the basis.
pub fn buildDualPhaseOneCosts(self: *SimplexEngine, problem: problem_module.ProblemView) void {
    const basis = if (self.basis) |*value| value else return;
    const workspace = &self.dual_work;
    const original_count = problem.num_cols + problem.num_rows;
    const maximize = problem.objective_sense == .maximize;
    if (self.active_dual_initialization_strategy == .highs)
        highs_random.initializeVectors(
            workspace.perturbation[0..original_count],
            workspace.highs_permutation[0..original_count],
            problem.num_cols,
        );
    workspace.resetCorrectionRandom(problem.num_cols, original_count);
    var maximum: f64 = 0.0;
    var boxed: usize = 0;
    for (0..problem.num_cols) |column| {
        const cost = (if (maximize) -problem.col_cost[column] else problem.col_cost[column]) *
            self.objective_scale * basis.column_scale[column];
        workspace.work_cost[column] = cost;
        maximum = @max(maximum, @abs(cost));
        if (!enable_highs_rhs_infeasibility_list and
            std.math.isFinite(workspace.saved_lower[column]) and std.math.isFinite(workspace.saved_upper[column]) and
            workspace.saved_lower[column] != workspace.saved_upper[column]) boxed += 1;
    }
    if (enable_highs_rhs_infeasibility_list) for (workspace.work_range[0..original_count]) |range| {
        if (range < 1e30) boxed += 1;
    };
    if (maximum > 100.0) maximum = @sqrt(@sqrt(maximum));
    const boxed_denominator = if (enable_highs_rhs_infeasibility_list) original_count else problem.num_cols;
    if (boxed_denominator != 0 and boxed * 100 < boxed_denominator) maximum = @min(maximum, 1.0);
    const base = 5e-7 * if (enable_highs_rhs_infeasibility_list) maximum else @max(maximum, 1.0);
    for (0..problem.num_cols) |column| {
        const random_value = if (enable_highs_rhs_infeasibility_list)
            workspace.perturbation[column]
        else
            deterministicUnit(column, workspace.basis_epoch);
        const lower = workspace.saved_lower[column];
        const upper = workspace.saved_upper[column];
        var perturbation: f64 = 0.0;
        if (lower != upper) {
            const magnitude = if (enable_highs_rhs_infeasibility_list)
                (1.0 + random_value) * (@abs(workspace.work_cost[column]) + 1.0) * base
            else
                random_value * base;
            if (std.math.isFinite(lower) and !std.math.isFinite(upper)) {
                perturbation = magnitude;
            } else if (!std.math.isFinite(lower) and std.math.isFinite(upper)) {
                perturbation = -magnitude;
            } else if (std.math.isFinite(lower) and std.math.isFinite(upper)) {
                perturbation = if (workspace.work_cost[column] >= 0.0) magnitude else -magnitude;
            }
        }
        workspace.perturbation[column] = perturbation;
        workspace.work_cost[column] += perturbation;
    }
    for (problem.num_cols..original_count) |column| {
        const perturbation = if (self.active_dual_initialization_strategy == .highs)
            (0.5 - workspace.perturbation[column]) * 1e-12
        else
            deterministicUnit(column, workspace.basis_epoch) * 1e-12;
        workspace.perturbation[column] = perturbation;
        workspace.work_cost[column] = perturbation;
    }
}

/// First `initialiseCost(..., phase unknown)` pass from HEkkDual::solve.
/// This materializes the original scaled objective in the persistent work
/// arrays before the unperturbed dual/feasibility assessment. Logical costs
/// are zero in the original LP.
fn initializeUnperturbedDualCosts(self: *SimplexEngine, problem: problem_module.ProblemView) void {
    const basis = if (self.basis) |*value| value else return;
    const workspace = &self.dual_work;
    const original_count = problem.num_cols + problem.num_rows;
    const maximize = problem.objective_sense == .maximize;
    for (0..problem.num_cols) |column| {
        workspace.work_cost[column] =
            (if (maximize) -problem.col_cost[column] else problem.col_cost[column]) *
            self.objective_scale * basis.column_scale[column];
        workspace.perturbation[column] = 0.0;
    }
    @memset(workspace.work_cost[problem.num_cols..original_count], 0.0);
    @memset(workspace.perturbation[problem.num_cols..original_count], 0.0);
    @memset(workspace.work_shift[0..original_count], 0.0);
}

/// Map `(column, basis epoch)` to a reproducible sample in `(0, 1]`.
///
/// The SplitMix64-style integer avalanche avoids shared mutable RNG state and
/// is used only by the deterministic baseline initialization path. The upper
/// 53 bits are converted exactly to an IEEE-754-compatible unit sample.
fn deterministicUnit(column: usize, epoch: u64) f64 {
    var value = @as(u64, @intCast(column)) +% epoch *% 0x9e3779b97f4a7c15;
    value = (value ^ (value >> 30)) *% 0xbf58476d1ce4e5b9;
    value = (value ^ (value >> 27)) *% 0x94d049bb133111eb;
    value ^= value >> 31;
    return (@as(f64, @floatFromInt(value >> 11)) + 1.0) * 0x1.0p-53;
}

/// Repair a warm basis that is neither primal nor dual feasible by using
/// the zero auxiliary objective. Every reduced cost is then exactly zero,
/// so dual pivots can restore primal feasibility without discarding the
/// imported basis or factorization. The original objective is restored by
/// the caller before entering Phase II.
pub fn repairWarmBasisWithDual(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
    const phase_started = self.statisticsTimestamp();
    const iteration_started = self.iterations;
    defer self.recordPhaseElapsed(.dual_feasibility_repair, phase_started);
    defer self.recordPhaseIterations(.dual_feasibility_repair, iteration_started);
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    self.algorithm = .dual_revised;
    @memset(basis.reduced_cost, 0.0);
    while (self.iterations < control.max_iterations) : (self.iterations += 1) {
        if (self.beginIteration(problem, control, .dual_feasibility_repair)) |status| return status;
        const leaving = self.pricing.chooseDualLeavingWeighted(
            basis.basic_value,
            basis.basic_lower,
            basis.basic_upper,
            basis.row_edge_weight,
            self.numerical.primal_tolerance,
        ) orelse return .optimal;
        if (self.computeDualTableauRow(problem, leaving.row) != .optimal) return .numerical_failure;
        const original_cols = problem.num_cols + problem.num_rows;
        var active_ratio_test = self.ratio_test;
        if (self.numerical.anti_cycling_active) active_ratio_test.rule = .standard;
        const entering = active_ratio_test.chooseDualEntering(
            basis.tableau[0..original_cols],
            basis.reduced_cost[0..original_cols],
            basis.col_status[0..original_cols],
            &.{},
            basis.col_lower[0..original_cols],
            basis.col_upper[0..original_cols],
            basis.primal[0..original_cols],
            leaving.bound,
            leaving.violation,
            basis.dual_ratio[0..original_cols],
            basis.dual_direction[0..original_cols],
            basis.flip_columns[0..original_cols],
            .{},
        );
        if (self.applyBoundFlips(problem, entering.flip_count) != .optimal) return .numerical_failure;
        const entering_col = entering.column orelse return .infeasible;
        const entering_index: usize = @intCast(entering_col);
        if (self.computeDirection(problem, entering_index) != .optimal) return .numerical_failure;
        if (entering.direction < 0.0) {
            for (basis.pivot_direction) |*value| value.* = -value.*;
        }
        const leaving_row: usize = @intCast(leaving.row);
        const target = if (leaving.bound == .at_lower) basis.basic_lower[leaving_row] else basis.basic_upper[leaving_row];
        const pivot = basis.pivot_direction[leaving_row];
        if (@abs(pivot) <= self.numerical.pivot_tolerance) return .numerical_failure;
        const step = (basis.basic_value[leaving_row] - target) / pivot;
        if (!std.math.isFinite(step) or step < -self.numerical.primal_tolerance) return .numerical_failure;
        const leaving_column = basis.basic_index[leaving_row];
        if (self.performPivot(
            problem,
            entering_index,
            entering.direction,
            leaving_row,
            leaving.bound,
            @max(step, 0.0),
        ) != .optimal) return .numerical_failure;
        self.observeIterationStep(
            @max(step, 0.0),
            entering_index,
            entering.direction,
            leaving_column,
            0,
            null,
            null,
            false,
            entering.flip_count,
        );
        @memset(basis.reduced_cost, 0.0);
    }
    return .iteration_limit;
}

/// Materialize the pivotal dual tableau row for `leaving_row`.
///
/// The routine computes `e_p^T B^-1` by sparse BTRAN, validates/corrects the
/// selected row's exact steepest-edge weight, and prices that vector against
/// every structural, logical, and artificial column. Row-wise pricing is used
/// for sufficiently sparse BTRAN results; otherwise the column representation
/// is scanned. On success `basis.dual_row`, `basis.tableau`, and
/// `dual_row_index` describe the same leaving row.
pub fn computeDualTableauRow(self: *SimplexEngine, problem: problem_module.ProblemView, leaving_row: u32) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    const row: usize = @intCast(leaving_row);
    if (row >= problem.num_rows) return .numerical_failure;
    @memset(basis.dual_row, 0.0);
    basis.dual_row[row] = 1.0;
    // Unit-vector BTRAN for pivotal row: wire sparse-index dispatch.
    const ep_index = [_]u32{@intCast(row)};
    self.factorization.solveTransposeSparse(basis.dual_row, &ep_index) catch return .numerical_failure;
    self.observeEpDensity(basis.dual_row);
    if (self.pricing.rule == .steepest_edge and self.dual_edge_weights_valid) {
        var exact_weight: f64 = 0.0;
        for (basis.dual_row) |entry| exact_weight += entry * entry;
        if (!std.math.isFinite(exact_weight) or exact_weight <= self.numerical.zero_tolerance)
            return .numerical_failure;
        const updated_weight = basis.row_edge_weight[row];
        if (self.active_dual_initialization_strategy == .highs) {
            // HiGHS CHUZR always overwrites the selected row's recurrence
            // weight with the exact BTRAN norm. Candidate rejection is kept
            // in chooseDualLeavingWithExactDseValidation so the corrected
            // row participates in a fresh pricing scan.
            basis.row_edge_weight[row] = exact_weight;
            if (updated_weight != exact_weight)
                self.numerical.dual_edge_weight_corrections += 1;
        } else {
            const relative_error = @abs(updated_weight - exact_weight) / @max(exact_weight, self.numerical.zero_tolerance);
            if (relative_error > self.numerical.dual_edge_weight_error_tolerance) {
                basis.row_edge_weight[row] = exact_weight;
                self.numerical.dual_edge_weight_corrections += 1;
                // Legacy safeguard retained for the baseline comparison
                // track. The HiGHS-aligned track rejects and reprices before
                // accepting the leaving row instead of switching to Devex.
                if (updated_weight < 0.5 * exact_weight) self.dual_edge_weights_valid = false;
            }
        }
    }
    self.dual_row_index = leaving_row;
    @memset(basis.tableau, 0.0);
    const row_pricing = self.selectRowPricing(basis.dual_row);
    const pricing_started = self.statisticsTimestamp();
    defer if (row_pricing)
        self.recordRowPricingElapsed(pricing_started)
    else
        self.recordPricingElapsed(pricing_started);
    if (row_pricing) {
        for (basis.dual_row, 0..) |dual_value, row_index| {
            if (@abs(dual_value) <= self.numerical.zero_tolerance) continue;
            const scaled_dual = dual_value * basis.row_scale[row_index];
            for (self.pricing_row_view.rowColumns(row_index), self.pricing_row_view.rowValues(row_index)) |column_u32, coefficient| {
                const column: usize = @intCast(column_u32);
                basis.tableau[column] += scaled_dual * coefficient * basis.column_scale[column];
            }
        }
    } else {
        for (0..problem.num_cols) |column| {
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            var alpha: f64 = 0.0;
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |matrix_row, coefficient| {
                if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                const row_index = matrix_row.toUsize();
                alpha += basis.dual_row[row_index] * basis.row_scale[row_index] * coefficient * basis.column_scale[column];
            }
            basis.tableau[column] = alpha;
        }
    }
    for (0..problem.num_rows) |logical_row| {
        basis.tableau[problem.num_cols + logical_row] = basis.dual_row[logical_row];
    }
    var nonzero_count: usize = 0;
    for (basis.tableau[0 .. problem.num_cols + problem.num_rows]) |value| {
        if (@abs(value) > self.numerical.zero_tolerance) nonzero_count += 1;
    }
    const tableau_count = problem.num_cols + problem.num_rows;
    self.dual_hyper_sparse_active = tableau_count > 0 and nonzero_count * 10 < tableau_count;
    if (!self.dual_hyper_sparse_active) self.dual_candidate_count = 0;
    return .optimal;
}

/// HiGHS CHUZR safeguard for dual steepest-edge pricing. A recurrence weight
/// that is too small can make a row look artificially attractive, so BTRAN
/// the selected row, install its exact weight and repeat pricing when the old
/// value was below one quarter of the exact value. Overestimates are accepted.
fn chooseDualLeavingWithExactDseValidation(
    self: *SimplexEngine,
    problem: problem_module.ProblemView,
    leaving_choice: *?pricing_module.DualLeavingChoice,
) SolveStatus {
    leaving_choice.* = null;
    while (true) {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        const leaving = if (self.active_dual_initialization_strategy == .highs and
            self.dual_work.infeasibility_list_active)
        blk: {
            const row_u32 = self.dual_work.choosePrimalInfeasibleRow(basis, problem.num_rows) orelse return .optimal;
            const row: usize = @intCast(row_u32);
            const value = basis.basic_value[row];
            if (value < basis.basic_lower[row] - self.numerical.primal_tolerance) {
                break :blk pricing_module.DualLeavingChoice{
                    .row = row_u32,
                    .bound = .at_lower,
                    .violation = basis.basic_lower[row] - value,
                };
            }
            if (value > basis.basic_upper[row] + self.numerical.primal_tolerance) {
                break :blk pricing_module.DualLeavingChoice{
                    .row = row_u32,
                    .bound = .at_upper,
                    .violation = value - basis.basic_upper[row],
                };
            }
            return .optimal;
        } else self.pricing.chooseDualLeaving(self) orelse return .optimal;
        const row: usize = @intCast(leaving.row);
        if (row >= basis.row_edge_weight.len) return .numerical_failure;
        const updated_weight = basis.row_edge_weight[row];
        if (self.computeDualTableauRow(problem, leaving.row) != .optimal)
            return .numerical_failure;

        if (self.active_dual_initialization_strategy == .highs and
            self.pricing.rule == .steepest_edge and self.dual_edge_weights_valid)
        {
            const exact_weight = basis.row_edge_weight[row];
            if (updated_weight < 0.25 * exact_weight) {
                self.stats.dual_dse_weight_rejections += 1;
                continue;
            }
        }
        leaving_choice.* = leaving;
        return .optimal;
    }
}

test "HiGHS CHUZR rejects an underestimated DSE weight and reprices" {
    const problem = problem_module.ProblemView{
        .num_rows = 2,
        .num_cols = 0,
        .col_cost = &.{},
        .col_lower = &.{},
        .col_upper = &.{},
        .row_lower = &[_]f64{ 0, 0 },
        .row_upper = &[_]f64{ std.math.inf(f64), std.math.inf(f64) },
        .matrix = matrix.CscView.initAssumeValid(2, 0, &[_]usize{0}, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 2, 0);
    engine.basis.?.initializeSlackBasis();
    engine.basis.?.basic_value[0] = -10;
    engine.basis.?.basic_value[1] = -2;
    @memset(engine.basis.?.basic_lower, 0);
    @memset(engine.basis.?.basic_upper, std.math.inf(f64));
    @memset(engine.basis.?.row_edge_weight, 1);
    engine.pricing.rule = .steepest_edge;
    engine.dual_edge_weights_valid = true;
    engine.active_dual_initialization_strategy = .highs;
    // B^-T e_0 has squared norm 100, while B^-T e_1 has norm 1.
    // The stale weights first select row 0 (violation 10), but its corrected
    // score is 1, so the fresh CHUZR scan must select row 1 (score 2).
    try engine.factorization.factorize(2, &[_]f64{ 0.1, 0, 0, 1 });

    var leaving: ?pricing_module.DualLeavingChoice = null;
    try std.testing.expectEqual(
        SolveStatus.optimal,
        chooseDualLeavingWithExactDseValidation(&engine, problem, &leaving),
    );
    try std.testing.expectEqual(@as(u32, 1), leaving.?.row);
    try std.testing.expectApproxEqAbs(@as(f64, 100), engine.basis.?.row_edge_weight[0], 1e-12);
    try std.testing.expectEqual(@as(usize, 1), engine.stats.dual_dse_weight_rejections);
    try std.testing.expect(engine.dual_edge_weights_valid);
    try std.testing.expectEqual(pricing_module.PricingRule.steepest_edge, engine.pricing.rule);
}

/// Atomically apply the BFRT flip prefix selected by CHUZC.
///
/// All column displacements are accumulated in equation space and transformed
/// with one FTRAN. The solve is residual-checked and retried after reinversion
/// when necessary. No status, primal value, or infeasibility-list entry is
/// published until that aggregate solve succeeds, so numerical failure cannot
/// leave a partially committed flip batch. `flip_count` addresses the prefix
/// of `basis.flip_columns` produced by the ratio test.
pub fn applyBoundFlips(self: *SimplexEngine, problem: problem_module.ProblemView, flip_count: usize) SolveStatus {
    if (flip_count == 0) return .optimal;
    self.stats.bound_flips += flip_count;
    self.stats.bound_flip_batches += 1;
    self.stats.bound_flip_ftran_savings += flip_count - 1;
    const basis = if (self.basis) |*value| value else return .numerical_failure;

    // Accumulate every bound displacement in equation space, then apply
    // B^-1 once. A pivotal tableau row cannot update the whole basic
    // vector; batching the RHS is the exact allocation-free equivalent of
    // one FTRAN per flipped column.
    if (self.accumulateBoundFlipRhs(problem, flip_count, basis.residual_work) != .optimal)
        return .numerical_failure;
    @memcpy(basis.rhs_work, basis.residual_work);
    self.factorization.solve(basis.rhs_work) catch return .numerical_failure;
    for (basis.rhs_work) |value| if (!std.math.isFinite(value)) return .numerical_failure;
    if (!self.boundFlipResidualAcceptable(problem)) {
        self.factorization.recordReinversion(.solve_residual);
        if (self.refactorizeBasis(problem, .solve_residual) != .optimal) return .numerical_failure;
        // Reinversion and its validation reuse row workspaces, so rebuild
        // the immutable aggregate before retrying the exact same batch.
        if (self.accumulateBoundFlipRhs(problem, flip_count, basis.residual_work) != .optimal)
            return .numerical_failure;
        @memcpy(basis.rhs_work, basis.residual_work);
        self.factorization.solve(basis.rhs_work) catch return .numerical_failure;
        for (basis.rhs_work) |value| if (!std.math.isFinite(value)) return .numerical_failure;
        if (!self.boundFlipResidualAcceptable(problem)) return .numerical_failure;
    }
    self.dual_work.changed_row_count = self.factorization.gatherFtranResultIndices(
        basis.rhs_work,
        self.dual_work.changed_row_index,
    );

    // Commit only after the aggregate solve succeeds. Thus an invalid
    // width, column, or numerical solve never publishes a partial batch.
    for (basis.basic_value, basis.rhs_work) |*value, displacement| value.* -= displacement;
    for (basis.flip_columns[0..flip_count]) |column_u32| {
        const column: usize = @intCast(column_u32);
        const delta = switch (basis.col_status[column]) {
            .at_lower => basis.col_upper[column] - basis.col_lower[column],
            .at_upper => -(basis.col_upper[column] - basis.col_lower[column]),
            else => return .numerical_failure,
        };
        basis.primal[column] += delta;
        basis.col_status[column] = if (delta > 0.0) .at_upper else .at_lower;
        self.dual_work.noteBoundFlip(column);
    }
    for (basis.basic_index, basis.basic_value) |basic_column, value| basis.primal[basic_column] = value;
    self.dual_work.updatePrimalInfeasibilityList(
        basis,
        self.dual_work.changed_row_index[0..self.dual_work.changed_row_count],
        self.numerical.primal_tolerance,
    );
    return .optimal;
}

/// Check the aggregate solve against the current basis before any flip is
/// committed. The same normwise backward-error scale used by entering
/// FTRAN protects long FT chains without adding another solve.
pub fn boundFlipResidualAcceptable(self: *SimplexEngine, problem: problem_module.ProblemView) bool {
    const basis = if (self.basis) |*value| value else return false;
    @memcpy(basis.pivot_direction, basis.residual_work);
    for (basis.basic_margin, basis.residual_work) |*magnitude, rhs| magnitude.* = @abs(rhs);
    self.subtractBasisProductWithMagnitude(
        problem,
        basis.rhs_work,
        basis.pivot_direction,
        basis.basic_margin,
    ) catch return false;
    var residual_max: f64 = 0.0;
    var equation_scale: f64 = 0.0;
    for (basis.pivot_direction) |value| residual_max = @max(residual_max, @abs(value));
    for (basis.basic_margin) |value| equation_scale = @max(equation_scale, value);
    const relative = residual_max / @max(1.0, equation_scale);
    self.numerical.last_ftran_relative_residual = relative;
    self.numerical.max_ftran_relative_residual = @max(self.numerical.max_ftran_relative_residual, relative);
    return std.math.isFinite(relative) and relative <= self.numerical.residual_tolerance;
}

/// Form the combined scaled equation-space displacement for a batch of
/// boxed nonbasic bound flips. The caller owns `output`; no column is
/// materialized and no allocation occurs in this hot path.
pub fn accumulateBoundFlipRhs(
    self: *SimplexEngine,
    problem: problem_module.ProblemView,
    flip_count: usize,
    output: []f64,
) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    if (output.len != problem.num_rows or flip_count > basis.flip_columns.len) return .numerical_failure;
    @memset(output, 0.0);
    for (basis.flip_columns[0..flip_count]) |column_u32| {
        const column: usize = @intCast(column_u32);
        if (column >= problem.num_cols + problem.num_rows) return .numerical_failure;
        const delta = switch (basis.col_status[column]) {
            .at_lower => basis.col_upper[column] - basis.col_lower[column],
            .at_upper => -(basis.col_upper[column] - basis.col_lower[column]),
            else => return .numerical_failure,
        };
        if (!std.math.isFinite(delta)) return .numerical_failure;
        if (column < problem.num_cols) {
            const column_scale = basis.column_scale[column];
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
                if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                const row_index = row.toUsize();
                if (row_index >= output.len) return .numerical_failure;
                output[row_index] += delta * coefficient * basis.row_scale[row_index] * column_scale;
            }
        } else {
            output[column - problem.num_cols] += delta;
        }
    }
    for (output) |value| if (!std.math.isFinite(value)) return .numerical_failure;
    return .optimal;
}

/// Maintain `r' = r - (r_q / alpha_pq) alpha_p` from the tableau row of
/// the old basis. `pivot_direction` may have been negated for movement
/// from an upper bound, but basis replacement still installs the original
/// `A_q`, so the unsigned tableau pivot is the required denominator.
pub fn updateReducedCostsAfterDualPivot(
    self: *SimplexEngine,
    problem: problem_module.ProblemView,
    entering_col: usize,
    leaving_row: u32,
    zero_dual_step: bool,
) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    const row: usize = @intCast(leaving_row);
    const original_count = problem.num_cols + problem.num_rows;
    if (row >= problem.num_rows or entering_col >= original_count) return .numerical_failure;
    const leaving_col = basis.basic_index[row];
    if (leaving_col >= original_count) return .numerical_failure;
    const alpha_pq = basis.tableau[entering_col];
    if (!std.math.isFinite(alpha_pq) or @abs(alpha_pq) <= self.numerical.pivot_tolerance)
        return .numerical_failure;
    if (zero_dual_step) {
        const shift = -basis.reduced_cost[entering_col];
        if (!std.math.isFinite(shift) or !self.dual_work.shiftCost(entering_col, shift))
            return .numerical_failure;
    } else {
        const theta = basis.reduced_cost[entering_col] / alpha_pq;
        if (!std.math.isFinite(theta)) return .numerical_failure;
        for (basis.reduced_cost[0..original_count], basis.tableau[0..original_count]) |*reduced, alpha|
            reduced.* -= theta * alpha;
        const artificial_begin = original_count;
        for (0..problem.num_rows) |logical_row|
            basis.reduced_cost[artificial_begin + logical_row] -=
                theta * basis.dual_row[logical_row] * basis.artificial_sign[logical_row];
    }
    basis.reduced_cost[entering_col] = 0.0;
    if (!self.dual_work.shiftBack(basis.reduced_cost, leaving_col)) return .numerical_failure;
    for (basis.reduced_cost) |value| if (!std.math.isFinite(value)) return .numerical_failure;
    self.reduced_cost_update_count += 1;
    self.stats.dual_reduced_cost_updates += 1;
    return .optimal;
}

test "dual bound flips share one aggregate FTRAN" {
    const rows = [_]foundation.RowId{
        foundation.RowId.fromUsizeAssumeValid(0),
        foundation.RowId.fromUsizeAssumeValid(1),
        foundation.RowId.fromUsizeAssumeValid(0),
        foundation.RowId.fromUsizeAssumeValid(1),
    };
    const problem = problem_module.ProblemView{
        .num_rows = 2,
        .num_cols = 2,
        .col_cost = &[_]f64{ 0, 0 },
        .col_lower = &[_]f64{ 0, 0 },
        .col_upper = &[_]f64{ 2, 3 },
        .row_lower = &[_]f64{ -std.math.inf(f64), -std.math.inf(f64) },
        .row_upper = &[_]f64{ 10, 10 },
        .matrix = matrix.CscView.initAssumeValid(
            2,
            2,
            &[_]usize{ 0, 2, 4 },
            &rows,
            &[_]f64{ 1, 2, 3, 4 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 2, 2);
    engine.basis.?.initializeSlackBasis();
    engine.basis.?.col_lower[0] = 0;
    engine.basis.?.col_lower[1] = 0;
    engine.basis.?.col_upper[0] = 2;
    engine.basis.?.col_upper[1] = 3;
    engine.basis.?.basic_value[0] = 10;
    engine.basis.?.basic_value[1] = 10;
    engine.basis.?.primal[2] = 10;
    engine.basis.?.primal[3] = 10;
    engine.basis.?.flip_columns[0] = 0;
    engine.basis.?.flip_columns[1] = 1;
    try engine.factorization.factorizeIdentity(2);

    try std.testing.expectEqual(SolveStatus.optimal, engine.applyBoundFlips(problem, 2));
    // Combined equation displacement is [1*2 + 3*3, 2*2 + 4*3] = [11, 16].
    try std.testing.expectApproxEqAbs(@as(f64, -1), engine.basis.?.basic_value[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, -6), engine.basis.?.basic_value[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2), engine.basis.?.primal[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 3), engine.basis.?.primal[1], 1e-12);
    try std.testing.expectEqual(basis_module.BasisStatus.at_upper, engine.basis.?.col_status[0]);
    try std.testing.expectEqual(basis_module.BasisStatus.at_upper, engine.basis.?.col_status[1]);
    try std.testing.expectEqual(@as(usize, 1), engine.factorization.stats.ftran_calls);
    try std.testing.expectEqual(@as(usize, 1), engine.stats.bound_flip_batches);
    try std.testing.expectEqual(@as(usize, 1), engine.stats.bound_flip_ftran_savings);
}

test "neither-feasible imported basis is repaired without a crash restart" {
    const Context = struct {
        saw_repair: bool = false,
        /// Test callback detecting entry into dual-feasibility repair.
        fn callback(event: ProgressEventView, context_ptr: ?*anyopaque) CallbackAction {
            const self: *@This() = @ptrCast(@alignCast(context_ptr.?));
            self.saw_repair = self.saw_repair or event.phase == .dual_feasibility_repair;
            return .continue_solve;
        }
    };
    const rows = [_]foundation.RowId{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) };
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &[_]f64{ 0, -1 },
        .col_lower = &[_]f64{ 0, 0 },
        .col_upper = &[_]f64{ 1, std.math.inf(f64) },
        .row_lower = &[_]f64{2},
        .row_upper = &[_]f64{2},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &[_]usize{ 0, 1, 2 }, &rows, &[_]f64{ 1, 1 }),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    const initial_basis = basis_snapshot_module.BasisView{
        .structural_status = &[_]basis_module.BasisStatus{ .basic, .at_lower },
        .logical_status = &[_]basis_module.BasisStatus{.fixed},
        .basic_index = &[_]u32{0},
    };
    var context = Context{};
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{
        .initial_basis = initial_basis,
        .iteration_callback = Context.callback,
        .callback_user_data = &context,
    }));
    try std.testing.expect(context.saw_repair);
    try std.testing.expectEqual(Algorithm.primal_revised, engine.algorithm);
    try std.testing.expectApproxEqAbs(@as(f64, 0), engine.basis.?.primal[0], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 2), engine.basis.?.primal[1], 1e-10);
}

test "dual rank one reduced cost update uses the old tableau row" {
    const rows = [_]foundation.RowId{
        foundation.RowId.fromUsizeAssumeValid(0),
        foundation.RowId.fromUsizeAssumeValid(0),
    };
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &[_]f64{ 0, 0 },
        .col_lower = &[_]f64{ 0, 0 },
        .col_upper = &[_]f64{ 1, 1 },
        .row_lower = &[_]f64{0},
        .row_upper = &[_]f64{0},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &[_]usize{ 0, 1, 2 }, &rows, &[_]f64{ 4, 3 }),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 1, 2);
    engine.basis.?.reduced_cost[0] = 2;
    engine.basis.?.reduced_cost[1] = -1;
    engine.basis.?.tableau[0] = 4;
    engine.basis.?.tableau[1] = 3;
    engine.basis.?.tableau[2] = 1;
    engine.basis.?.dual_row[0] = 1;
    engine.basis.?.artificial_sign[0] = -1;

    try std.testing.expectEqual(SolveStatus.optimal, engine.updateReducedCostsAfterDualPivot(problem, 0, 0, false));
    try std.testing.expectApproxEqAbs(@as(f64, 0), engine.basis.?.reduced_cost[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, -2.5), engine.basis.?.reduced_cost[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, -0.5), engine.basis.?.reduced_cost[2], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), engine.basis.?.reduced_cost[3], 1e-12);
    try std.testing.expectEqual(@as(usize, 1), engine.reduced_cost_update_count);
    try std.testing.expectEqual(@as(usize, 1), engine.stats.dual_reduced_cost_updates);
}

test "fresh dual no-entering publishes a validated original-coordinate Farkas ray" {
    const rows = [_]foundation.RowId{
        foundation.RowId.fromUsizeAssumeValid(0),
        foundation.RowId.fromUsizeAssumeValid(1),
    };
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
    const status = engine.solveProblem(problem, .{ .phase_one_strategy = .dual });
    try std.testing.expectEqual(SolveStatus.infeasible, status);
    try std.testing.expect(engine.infeasibility_ray_valid);
    try std.testing.expect(engine.infeasibility_certificate_gap > engine.numerical.primal_tolerance);
    const solution = engine.solutionView(problem, status).?;
    try std.testing.expectEqual(@as(usize, 2), solution.infeasibility_ray.len);
    const coefficient = solution.infeasibility_ray[0] + solution.infeasibility_ray[1];
    try std.testing.expectApproxEqAbs(@as(f64, 0), coefficient, 1e-12);
}

test "Farkas validation ignores only sub-floor infinite-bound column mass" {
    const row = foundation.RowId.fromUsizeAssumeValid(0);
    const column_rows = [_]foundation.RowId{ row, row };
    const leaving = pricing_module.DualLeavingChoice{
        .row = 0,
        .bound = .at_upper,
        .violation = 1,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 1, 2);
    engine.basis.?.row_scale[0] = 1;
    engine.basis.?.dual_row[0] = 1;

    const accepted_problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &[_]f64{ 0, 0 },
        .col_lower = &[_]f64{ 0, 0 },
        .col_upper = &[_]f64{ 0, std.math.inf(f64) },
        .row_lower = &[_]f64{1},
        .row_upper = &[_]f64{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &[_]usize{ 0, 1, 2 }, &column_rows, &[_]f64{ 1, 5e-10 }),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    try std.testing.expect(engine.buildAndValidateInfeasibilityRay(accepted_problem, leaving));
    try std.testing.expectApproxEqAbs(@as(f64, 5e-10), engine.infeasibility_certificate_infinite_column_mass, 1e-20);
    try std.testing.expectApproxEqAbs(@as(f64, 1), engine.infeasibility_certificate_gap, 1e-12);

    engine.infeasibility_ray_valid = false;
    const rejected_problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &[_]f64{ 0, 0 },
        .col_lower = &[_]f64{ 0, 0 },
        .col_upper = &[_]f64{ 0, std.math.inf(f64) },
        .row_lower = &[_]f64{1},
        .row_upper = &[_]f64{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &[_]usize{ 0, 1, 2 }, &column_rows, &[_]f64{ 1, 2e-9 }),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    try std.testing.expect(!engine.buildAndValidateInfeasibilityRay(rejected_problem, leaving));
    try std.testing.expectEqual(
        @import("engine.zig").InfeasibilityCertificateFailure.infinite_column_bound,
        engine.infeasibility_certificate_failure,
    );
}
