//! Dual Devex/DSE edge-weight frameworks for `SimplexEngine`.
//!
//! ## Responsibility
//!
//! Owns exact dual steepest-edge initialization, the incremental DSE
//! recurrence with its transactional fallback to full dual Devex, and the
//! dual Devex weight update.

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

/// Full dual Devex recurrence used by HiGHS-style row pricing. The hot
/// FTRAN column is reused directly, so updating every row is allocation
/// free and requires no additional basis solve.
pub fn updateDualDevexWeights(self: *SimplexEngine, leaving_row: usize, pivot: f64) void {
    const basis = if (self.basis) |*value| value else return;
    const iteration = if (self.dual_iteration) |*value| value else return;
    if (leaving_row >= basis.row_edge_weight.len or @abs(pivot) <= self.numerical.pivot_tolerance) return;
    var pivotal_weight = basis.row_edge_weight[leaving_row] / (pivot * pivot);
    if (!std.math.isFinite(pivotal_weight)) pivotal_weight = 1.0;
    pivotal_weight = @max(1.0, pivotal_weight);
    for (basis.row_edge_weight, iteration.col_aq.values) |*weight, alpha| {
        const candidate = pivotal_weight * alpha * alpha;
        if (std.math.isFinite(candidate)) weight.* = @max(weight.*, candidate);
    }
    basis.row_edge_weight[leaving_row] = pivotal_weight;
    self.stats.dual_devex_updates += 1;
}

/// Initialize exact dual steepest-edge weights after a crash, imported
/// basis, reinversion, or detected recurrence drift.
pub fn ensureExactDualEdgeWeights(self: *SimplexEngine) SolveStatus {
    if (self.dual_edge_weights_valid) return .optimal;
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    const iteration = if (self.dual_iteration) |*value| value else return .numerical_failure;
    for (basis.row_edge_weight, 0..) |*weight, row| {
        iteration.row_ep.clear();
        iteration.row_ep.values[row] = 1.0;
        iteration.row_ep.appendIndex(@intCast(row)) catch return .numerical_failure;
        self.nla.btran(&iteration.row_ep) catch return .numerical_failure;
        iteration.row_ep.reIndex();
        const norm_squared = iteration.row_ep.norm2() catch return .numerical_failure;
        if (!std.math.isFinite(norm_squared)) return .numerical_failure;
        weight.* = @max(norm_squared, self.numerical.zero_tolerance);
    }
    self.dual_edge_weights_valid = true;
    self.dual_row_index = null;
    return .optimal;
}

/// Cause for abandoning exact DSE recurrence in favor of dual Devex.
const DualDseFallbackReason = enum { invalid, budget };

/// Start a fresh DSE framework at a dual-phase boundary. Returning the
/// previous pricing rule lets callers restore policy without retaining a
/// second engine-wide rule or introducing a hot-loop branch.
pub fn beginDualEdgeWeightPhase(self: *SimplexEngine) pricing_module.PricingRule {
    const saved = self.pricing.rule;
    if (self.active_dual_edge_weight_strategy == .steepest_devex) {
        if (self.dual_dse_certificate_fallback_active) {
            self.pricing.rule = .devex;
            return saved;
        }
        self.pricing.rule = .steepest_edge;
        // Preserve an exact framework across a phase boundary. In particular,
        // the logical basis starts with exact unit DSE weights; invalidating
        // them here made a no-pivot Phase-I probe and the following Phase II
        // each perform num_rows redundant BTRAN solves. Basis import, primal
        // pivots and Devex fallback explicitly clear this flag when the stored
        // weights no longer describe the current basis.
        self.dual_row_index = null;
        self.dual_dse_updates_since_start = 0;
        self.dual_dse_aq_density = 0.0;
        self.dual_dse_ep_density = 0.0;
        self.dual_dse_ap_density = 0.0;
        self.dual_dse_tau_density = 0.0;
        self.dual_dse_costly_iterations = 0;
        self.dual_dse_average_log_low_error = 0.0;
        self.dual_dse_average_log_high_error = 0.0;
    }
    return saved;
}

/// Replace a rejected or over-budget DSE framework with a deterministic
/// unit-weight dual Devex framework. Subsequent pivots use the full Devex
/// recurrence and never pay for another exact all-row DSE initialization
/// in the current phase.
pub fn switchDualDseToDevex(self: *SimplexEngine, reason: DualDseFallbackReason) void {
    const basis = if (self.basis) |*value| value else return;
    self.pricing.rule = .devex;
    @memset(basis.row_edge_weight, 1.0);
    self.dual_edge_weights_valid = false;
    self.dual_row_index = null;
    self.dual_candidate_count = 0;
    switch (reason) {
        .invalid => self.stats.dual_dse_invalid_fallbacks += 1,
        .budget => self.stats.dual_dse_budget_fallbacks += 1,
    }
}

/// Forrest--Goldfarb DSE recurrence. `dual_row` is the freshly computed
/// BTRAN result B^-T e_p. One additional FTRAN forms
/// tau = B^-1 B^-T e_p before the factor update changes B.
pub fn updateDualSteepestEdgeWeights(self: *SimplexEngine, leaving_row: usize, pivot: f64) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    const iteration = if (self.dual_iteration) |*value| value else return .numerical_failure;
    if (!self.dual_edge_weights_valid or leaving_row >= basis.row_edge_weight.len or
        @abs(pivot) <= self.numerical.pivot_tolerance)
        return .numerical_failure;
    var exact_pivot_weight = self.dual_row_exact_weight;
    var ep_nonzeros = self.dual_row_nonzero_count;
    if (self.dual_row_index != @as(u32, @intCast(leaving_row)) or
        !std.math.isFinite(exact_pivot_weight) or exact_pivot_weight <= self.numerical.zero_tolerance)
    {
        exact_pivot_weight = 0.0;
        ep_nonzeros = 0;
        for (iteration.row_ep.values) |entry| {
            exact_pivot_weight += entry * entry;
            if (@abs(entry) > self.numerical.zero_tolerance) ep_nonzeros += 1;
        }
    }
    if (!std.math.isFinite(exact_pivot_weight) or exact_pivot_weight <= self.numerical.zero_tolerance)
        return .numerical_failure;
    basis.row_edge_weight[leaving_row] = exact_pivot_weight;
    iteration.col_dse.copyFrom(&iteration.row_ep) catch return .numerical_failure;
    self.nla.ftran(&iteration.col_dse) catch return .numerical_failure;
    var aq_nonzeros: usize = 0;
    var tau_nonzeros: usize = 0;
    const new_pivotal_edge_weight = exact_pivot_weight / (pivot * pivot);
    const kai = -2.0 / pivot;
    for (basis.row_edge_weight, iteration.col_aq.values, iteration.col_dse.values, 0..) |*weight, alpha, tau, row| {
        if (@abs(alpha) > self.numerical.zero_tolerance) aq_nonzeros += 1;
        if (@abs(tau) > self.numerical.zero_tolerance) tau_nonzeros += 1;
        if (row == leaving_row) continue;
        // HiGHS traverses the indexed pivotal column and therefore never
        // evaluates the recurrence for exact structural zeros.
        if (alpha == 0.0) continue;
        // Preserve HiGHS' operation order: strict CHUZR ties are sensitive to
        // a one-ulp difference between this form and the algebraic ratio form.
        const updated = weight.* + alpha * (new_pivotal_edge_weight * alpha + kai * tau);
        if (!std.math.isFinite(updated)) return .numerical_failure;
        weight.* = @max(updated, 1e-4);
    }
    basis.row_edge_weight[leaving_row] = @max(new_pivotal_edge_weight, 1e-4);
    const row_count = basis.row_edge_weight.len;
    if (row_count != 0) {
        const inverse_rows = 1.0 / @as(f64, @floatFromInt(row_count));
        const multiplier = 0.05;
        self.dual_dse_aq_density = (1.0 - multiplier) * self.dual_dse_aq_density +
            multiplier * @as(f64, @floatFromInt(aq_nonzeros)) * inverse_rows;
        self.dual_dse_ep_density = (1.0 - multiplier) * self.dual_dse_ep_density +
            multiplier * @as(f64, @floatFromInt(ep_nonzeros)) * inverse_rows;
        self.dual_dse_tau_density = (1.0 - multiplier) * self.dual_dse_tau_density +
            multiplier * @as(f64, @floatFromInt(tau_nonzeros)) * inverse_rows;
        self.dual_dse_ap_density = (1.0 - multiplier) * self.dual_dse_ap_density +
            multiplier * self.current_dual_tableau_density;
        const denominator = @max(
            @max(self.dual_dse_ep_density, self.dual_dse_aq_density),
            self.dual_dse_ap_density,
        );
        const relative_density = if (denominator > 0.0) self.dual_dse_tau_density / denominator else 0.0;
        if (relative_density * relative_density > 1000.0 and self.dual_dse_tau_density > 0.01)
            self.dual_dse_costly_iterations += 1;
    }
    return .optimal;
}

/// HiGHS `switchToDevex` density/frequency criterion. The explicit budget is
/// a Phase-I safety cap: auxiliary costs can otherwise keep an inaccurate or
/// unproductive DSE framework alive until the controller abandons the dual
/// start. Original-objective Phase II follows HiGHS' cost/accuracy policy so
/// useful exact weights are not discarded merely because a fixed number of
/// pivots has elapsed. The unknown phase retains the hard-cap semantics for
/// direct framework users and unit tests.
pub fn shouldSwitchDualDseToDevex(self: *const SimplexEngine, total_columns: usize) bool {
    if (!self.active_adaptive_dual_dse_budget and
        self.active_dual_dse_update_budget != 0 and self.dual_control.phase != .phase_two)
        return self.dual_dse_updates_since_start >= self.active_dual_dse_update_budget;
    const costly = self.dual_dse_updates_since_start * 10 > total_columns and
        self.dual_dse_costly_iterations * 20 > self.dual_dse_updates_since_start;
    const inaccurate = self.dual_dse_average_log_low_error +
        self.dual_dse_average_log_high_error > 10.0;
    return costly or inaccurate;
}

test "dual DSE budget caps Phase I but not original-objective Phase II" {
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.active_dual_dse_update_budget = 64;
    engine.dual_dse_updates_since_start = 64;
    engine.active_adaptive_dual_dse_budget = false;

    // Direct framework users retain the explicit hard-cap behavior.
    try std.testing.expect(engine.shouldSwitchDualDseToDevex(1000));
    engine.dual_control.enterPhase(.phase_one);
    try std.testing.expect(engine.shouldSwitchDualDseToDevex(1000));

    // Automatic policy ignores the diagnostic cap and uses measured work.
    engine.active_adaptive_dual_dse_budget = true;
    try std.testing.expect(!engine.shouldSwitchDualDseToDevex(1000));

    // Phase II only switches for measured DSE cost or weight inaccuracy.
    engine.dual_control.enterPhase(.phase_two);
    try std.testing.expect(!engine.shouldSwitchDualDseToDevex(1000));
    engine.dual_dse_average_log_high_error = 10.01;
    try std.testing.expect(engine.shouldSwitchDualDseToDevex(1000));
}

test "exact dual steepest-edge weights use BTRAN row norms" {
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 2, 0);
    try engine.ensureDualIterationWorkspace(2, 0);
    engine.basis.?.basic_value[0] = -1.0;
    engine.basis.?.basic_value[1] = -2.0;
    @memset(engine.basis.?.basic_lower, 0.0);
    @memset(engine.basis.?.basic_upper, std.math.inf(f64));
    try engine.nla.factorization.factorize(2, &[_]f64{ 2, 0, 0, 0.5 });
    try std.testing.expectEqual(SolveStatus.optimal, engine.ensureExactDualEdgeWeights());
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), engine.basis.?.row_edge_weight[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), engine.basis.?.row_edge_weight[1], 1e-12);

    engine.pricing.rule = .steepest_edge;
    const choice = engine.pricing.chooseDualLeavingWeighted(
        engine.basis.?.basic_value,
        engine.basis.?.basic_lower,
        engine.basis.?.basic_upper,
        engine.basis.?.row_edge_weight,
        engine.numerical.primal_tolerance,
    ).?;
    try std.testing.expectEqual(@as(u32, 0), choice.row);
}

test "incremental dual steepest-edge recurrence matches the new inverse rows" {
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 2, 0);
    try engine.ensureDualIterationWorkspace(2, 0);
    try engine.nla.factorization.factorizeIdentity(2);
    engine.pricing.rule = .steepest_edge;
    engine.dual_edge_weights_valid = true;
    engine.dual_row_index = 0;
    engine.dual_iteration.?.row_ep.values[0] = 1.0;
    engine.dual_iteration.?.row_ep.values[1] = 0.0;
    try engine.dual_iteration.?.row_ep.appendIndex(0);
    engine.dual_iteration.?.col_aq.values[0] = 2.0;
    engine.dual_iteration.?.col_aq.values[1] = 1.0;
    try std.testing.expectEqual(SolveStatus.optimal, engine.updateDualSteepestEdgeWeights(0, 2.0));
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), engine.basis.?.row_edge_weight[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1.25), engine.basis.?.row_edge_weight[1], 1e-12);
}

test "dual steepest-edge recurrence matches exact weights for a nontrivial basis" {
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 2, 0);
    try engine.ensureDualIterationWorkspace(2, 0);
    try engine.nla.factorization.factorize(2, &[_]f64{ 2, 1, 1, 3 });
    engine.pricing.rule = .steepest_edge;
    engine.dual_edge_weights_valid = false;
    try std.testing.expectEqual(SolveStatus.optimal, engine.ensureExactDualEdgeWeights());

    engine.dual_iteration.?.row_ep.clear();
    engine.dual_iteration.?.row_ep.values[0] = 1;
    try engine.dual_iteration.?.row_ep.appendIndex(0);
    try engine.nla.btran(&engine.dual_iteration.?.row_ep);
    engine.dual_iteration.?.col_aq.values[0] = 4;
    engine.dual_iteration.?.col_aq.values[1] = -2;
    try engine.nla.factorization.solve(engine.dual_iteration.?.col_aq.values);
    const pivot = engine.dual_iteration.?.col_aq.values[0];
    try std.testing.expectEqual(SolveStatus.optimal, engine.updateDualSteepestEdgeWeights(0, pivot));
    const updated_weights = [_]f64{
        engine.basis.?.row_edge_weight[0],
        engine.basis.?.row_edge_weight[1],
    };
    try engine.nla.factorization.update(.{
        .leaving_row = 0,
        .entering_col = 0,
        .direction = engine.dual_iteration.?.col_aq.values,
        .column_scale = 1,
    });

    for (updated_weights, 0..) |updated, row| {
        @memset(engine.basis.?.rhs_work, 0);
        engine.basis.?.rhs_work[row] = 1;
        try engine.nla.factorization.solveTranspose(engine.basis.?.rhs_work);
        var exact: f64 = 0;
        for (engine.basis.?.rhs_work) |value| exact += value * value;
        try std.testing.expectApproxEqAbs(exact, updated, 1e-12);
    }
}

test "dual Devex updates every row from the hot FTRAN column" {
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 2, 0);
    try engine.ensureDualIterationWorkspace(2, 0);
    engine.basis.?.row_edge_weight[0] = 1;
    engine.basis.?.row_edge_weight[1] = 2;
    engine.dual_iteration.?.col_aq.values[0] = 2;
    engine.dual_iteration.?.col_aq.values[1] = 1;

    engine.updateDualDevexWeights(1, 1);

    try std.testing.expectApproxEqAbs(@as(f64, 8), engine.basis.?.row_edge_weight[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2), engine.basis.?.row_edge_weight[1], 1e-12);
    try std.testing.expectEqual(@as(usize, 1), engine.stats.dual_devex_updates);
}

test "rejected dual steepest-edge framework falls back transactionally" {
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 2, 0);
    engine.active_dual_edge_weight_strategy = .steepest_devex;
    engine.pricing.rule = .steepest_edge;
    engine.dual_edge_weights_valid = true;
    engine.dual_row_index = 1;
    engine.dual_candidate_count = 2;
    engine.basis.?.row_edge_weight[0] = 9;
    engine.basis.?.row_edge_weight[1] = 4;

    engine.switchDualDseToDevex(.invalid);

    try std.testing.expectEqual(pricing_module.PricingRule.devex, engine.pricing.rule);
    try std.testing.expectEqualSlices(f64, &[_]f64{ 1, 1 }, engine.basis.?.row_edge_weight);
    try std.testing.expect(!engine.dual_edge_weights_valid);
    try std.testing.expectEqual(@as(?u32, null), engine.dual_row_index);
    try std.testing.expectEqual(@as(usize, 0), engine.dual_candidate_count);
    try std.testing.expectEqual(@as(usize, 1), engine.stats.dual_dse_invalid_fallbacks);
}
