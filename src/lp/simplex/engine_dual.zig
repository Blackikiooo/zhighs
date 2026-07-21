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

pub fn solveDual(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
    const phase_started = self.statisticsTimestamp();
    const iteration_started = self.iterations;
    const saved_pricing_rule = self.beginDualEdgeWeightPhase();
    defer self.pricing.rule = saved_pricing_rule;
    defer self.recordPhaseElapsed(.phase_two, phase_started);
    defer self.recordPhaseIterations(.phase_two, iteration_started);
    if (self.recomputeReducedCosts(problem) != .optimal) return .numerical_failure;
    while (self.iterations < control.max_iterations) : (self.iterations += 1) {
        if (self.beginIteration(problem, control, .phase_two)) |status| return status;
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        if (self.pricing.rule == .steepest_edge and self.ensureExactDualEdgeWeights() != .optimal)
            return .numerical_failure;
        const leaving = self.pricing.chooseDualLeaving(self) orelse return self.finishOptimal(problem);

        if (self.computeDualTableauRow(problem, leaving.row) != .optimal) return .numerical_failure;
        if (self.active_dual_edge_weight_strategy == .steepest_devex and
            self.pricing.rule == .steepest_edge and !self.dual_edge_weights_valid)
            self.switchDualDseToDevex(.invalid);
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
        if (self.updateReducedCostsAfterDualPivot(problem, entering_index, leaving.row) != .optimal)
            return .numerical_failure;
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
        if (self.fresh_factorization_pivots_remaining != 0 or
            self.reduced_cost_update_count >= self.reduced_cost_refresh_period)
        {
            if (self.recomputeReducedCostsWithDrift(problem) != .optimal) return .numerical_failure;
            if (!self.classifyFeasibility(problem).dual) return .numerical_failure;
        }
    }
    return .iteration_limit;
}

/// Dedicated dual Phase I. Special bounds encode the negated dual
/// infeasibility objective, while deterministic perturbed costs remove
/// ratio ties. The borrowed matrix is never copied and every iteration
/// uses engine-owned work arrays.
pub fn solveDualPhaseOne(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
    const phase_started = self.statisticsTimestamp();
    const iteration_started = self.iterations;
    const saved_pricing_rule = self.beginDualEdgeWeightPhase();
    defer self.pricing.rule = saved_pricing_rule;
    defer self.recordPhaseElapsed(.dual_phase_one, phase_started);
    defer self.recordPhaseIterations(.dual_phase_one, iteration_started);
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    const original_count = problem.num_cols + problem.num_rows;
    self.dual_phase_one.begin(basis, original_count) catch return .numerical_failure;
    defer if (self.dual_phase_one.active)
        self.dual_phase_one.restoreOriginalBounds(basis, original_count);

    self.buildDualPhaseOneCosts(problem);
    if (self.recomputeReducedCostsFromWork(problem) != .optimal) return .not_implemented;
    // Infeasibility is encoded in nonbasic primal values (±1/0), not
    // shifted costs. The raw scaled LP cost + perturbation keeps the
    // reduced cost magnitude intact so the Phase-I objective correctly
    // measures remaining infeasibility.
    self.dual_phase_one.installWorkingBounds(basis, original_count);
    if (self.recomputeBasicValuesUnchecked(problem) != .optimal) return .not_implemented;
    self.algorithm = .dual_revised;
    self.dual_edge_weights_valid = false;

    while (self.iterations < control.max_iterations) : (self.iterations += 1) {
        if (self.beginIteration(problem, control, .dual_phase_one)) |status| return status;
        if (self.pricing.rule == .steepest_edge and self.ensureExactDualEdgeWeights() != .optimal)
            return .not_implemented;
        const leaving = self.pricing.chooseDualLeaving(self) orelse break;
        if (self.computeDualTableauRow(problem, leaving.row) != .optimal) return .not_implemented;
        if (self.active_dual_edge_weight_strategy == .steepest_devex and
            self.pricing.rule == .steepest_edge and !self.dual_edge_weights_valid)
            self.switchDualDseToDevex(.invalid);
        var active_ratio_test = self.ratio_test;
        if (self.numerical.anti_cycling_active) active_ratio_test.rule = .standard;
        const entering = active_ratio_test.chooseDualEntering(
            basis.tableau[0..original_count],
            basis.reduced_cost[0..original_count],
            basis.col_status[0..original_count],
            self.dual_phase_one.nonbasic_move[0..original_count],
            basis.col_lower[0..original_count],
            basis.col_upper[0..original_count],
            basis.primal[0..original_count],
            leaving.bound,
            leaving.violation,
            basis.dual_ratio[0..original_count],
            basis.dual_direction[0..original_count],
            basis.flip_columns[0..original_count],
        );
        self.dual_phase_one.recordRemainingViolation(
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
        if (self.performPivot(problem, entering_index, entering.direction, leaving_row, leaving.bound, @max(step, 0.0)) != .optimal)
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
        if (self.recomputeReducedCostsFromWork(problem) != .optimal) return .not_implemented;
    }
    if (self.iterations >= control.max_iterations) return .iteration_limit;

    // HiGHS-style Phase-I → Phase-II: do NOT restore original bounds
    // yet. The working subproblem is already primal-feasible; Phase II
    // runs with working bounds + original LP costs. The deferred
    // restoreOriginalBounds (line 130) recovers original bounds after
    // Phase II, and finishOptimal validates with original bounds.
    if (self.recomputeReducedCostsFromWork(problem) != .optimal) return .not_implemented;

    // Compute dual objective in working coordinates: Σ primal_j · d_j
    var dual_obj: f64 = 0.0;
    for (0..original_count) |j| {
        if (basis.col_status[j] == .basic) continue;
        dual_obj += basis.primal[j] * basis.reduced_cost[j];
    }

    if (dual_obj == 0.0) {
        return self.solveDual(problem, control);
    }

    // Nonzero dual objective → restore and fall back
    self.dual_phase_one.restoreOriginalBounds(basis, original_count);
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
        .leaving_bound = leaving.bound,
        .violation = leaving.violation,
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

pub fn buildDualPhaseOneCosts(self: *SimplexEngine, problem: problem_module.ProblemView) void {
    const basis = if (self.basis) |*value| value else return;
    const workspace = &self.dual_phase_one;
    const original_count = problem.num_cols + problem.num_rows;
    const maximize = problem.objective_sense == .maximize;
    var maximum: f64 = 0.0;
    var boxed: usize = 0;
    for (0..problem.num_cols) |column| {
        const cost = (if (maximize) -problem.col_cost[column] else problem.col_cost[column]) *
            self.objective_scale * basis.column_scale[column];
        workspace.work_cost[column] = cost;
        maximum = @max(maximum, @abs(cost));
        if (std.math.isFinite(workspace.saved_lower[column]) and std.math.isFinite(workspace.saved_upper[column]) and
            workspace.saved_lower[column] != workspace.saved_upper[column])
            boxed += 1;
    }
    if (maximum > 100.0) maximum = @sqrt(@sqrt(maximum));
    if (problem.num_cols != 0 and boxed * 100 < problem.num_cols) maximum = @min(maximum, 1.0);
    const base = 5e-7 * @max(maximum, 1.0);
    for (0..problem.num_cols) |column| {
        const random = deterministicUnit(column, workspace.basis_epoch);
        const lower = workspace.saved_lower[column];
        const upper = workspace.saved_upper[column];
        var perturbation: f64 = 0.0;
        if (lower != upper) {
            if (std.math.isFinite(lower) and !std.math.isFinite(upper)) {
                perturbation = random * base;
            } else if (!std.math.isFinite(lower) and std.math.isFinite(upper)) {
                perturbation = -random * base;
            } else if (std.math.isFinite(lower) and std.math.isFinite(upper)) {
                perturbation = if (workspace.work_cost[column] >= 0.0) random * base else -random * base;
            }
        }
        workspace.perturbation[column] = perturbation;
        workspace.work_cost[column] += perturbation;
    }
    for (problem.num_cols..original_count) |column| {
        const perturbation = deterministicUnit(column, workspace.basis_epoch) * 1e-12;
        workspace.perturbation[column] = perturbation;
        workspace.work_cost[column] = perturbation;
    }
}

/// Shift only currently infeasible nonbasic costs to the feasibility
/// boundary. Unlike the former zero-objective repair, feasible reduced
/// costs and all basic costs remain intact, so Phase I retains a meaningful
/// dual objective and deterministic ratio ordering.
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
        const relative_error = @abs(updated_weight - exact_weight) / @max(exact_weight, self.numerical.zero_tolerance);
        if (relative_error > self.numerical.dual_edge_weight_error_tolerance) {
            basis.row_edge_weight[row] = exact_weight;
            self.numerical.dual_edge_weight_corrections += 1;
            // Huangfu's safeguard rejects a seriously underestimated
            // pivotal weight; force a full exact reset after this pivot.
            if (updated_weight < 0.5 * exact_weight) self.dual_edge_weights_valid = false;
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
        self.dual_phase_one.noteBoundFlip(column);
    }
    for (basis.basic_index, basis.basic_value) |basic_column, value| basis.primal[basic_column] = value;
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
) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    const row: usize = @intCast(leaving_row);
    const original_count = problem.num_cols + problem.num_rows;
    if (row >= problem.num_rows or entering_col >= original_count) return .numerical_failure;
    const alpha_pq = basis.tableau[entering_col];
    if (!std.math.isFinite(alpha_pq) or @abs(alpha_pq) <= self.numerical.pivot_tolerance)
        return .numerical_failure;
    const theta = basis.reduced_cost[entering_col] / alpha_pq;
    if (!std.math.isFinite(theta)) return .numerical_failure;

    for (basis.reduced_cost[0..original_count], basis.tableau[0..original_count]) |*reduced, alpha|
        reduced.* -= theta * alpha;
    const artificial_begin = original_count;
    for (0..problem.num_rows) |logical_row|
        basis.reduced_cost[artificial_begin + logical_row] -=
            theta * basis.dual_row[logical_row] * basis.artificial_sign[logical_row];
    basis.reduced_cost[entering_col] = 0.0;
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

    try std.testing.expectEqual(SolveStatus.optimal, engine.updateReducedCostsAfterDualPivot(problem, 0, 0));
    try std.testing.expectApproxEqAbs(@as(f64, 0), engine.basis.?.reduced_cost[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, -2.5), engine.basis.?.reduced_cost[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, -0.5), engine.basis.?.reduced_cost[2], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), engine.basis.?.reduced_cost[3], 1e-12);
    try std.testing.expectEqual(@as(usize, 1), engine.reduced_cost_update_count);
    try std.testing.expectEqual(@as(usize, 1), engine.stats.dual_reduced_cost_updates);
}
