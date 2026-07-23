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
    if (leaving_row >= basis.row_edge_weight.len or @abs(pivot) <= self.numerical.pivot_tolerance) return;
    var pivotal_weight = basis.row_edge_weight[leaving_row] / (pivot * pivot);
    if (!std.math.isFinite(pivotal_weight)) pivotal_weight = 1.0;
    pivotal_weight = @max(1.0, pivotal_weight);
    for (basis.row_edge_weight, basis.pivot_direction) |*weight, alpha| {
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
    for (basis.row_edge_weight, 0..) |*weight, row| {
        @memset(basis.dual_row, 0.0);
        basis.dual_row[row] = 1.0;
        // Unit-vector BTRAN: pass the single nonzero index for
        // sparse-index adaptive dispatch (hyper-sparse L solve).
        const ep_index = [_]u32{@intCast(row)};
        self.factorization.solveTransposeSparse(basis.dual_row, &ep_index) catch return .numerical_failure;
        var norm_squared: f64 = 0.0;
        for (basis.dual_row) |entry| norm_squared += entry * entry;
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
        self.pricing.rule = .steepest_edge;
        self.dual_edge_weights_valid = false;
        self.dual_row_index = null;
        self.dual_dse_updates_since_start = 0;
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
    if (!self.dual_edge_weights_valid or leaving_row >= basis.row_edge_weight.len or
        @abs(pivot) <= self.numerical.pivot_tolerance)
        return .numerical_failure;
    var exact_pivot_weight: f64 = 0.0;
    for (basis.dual_row) |entry| exact_pivot_weight += entry * entry;
    if (!std.math.isFinite(exact_pivot_weight) or exact_pivot_weight <= self.numerical.zero_tolerance)
        return .numerical_failure;
    basis.row_edge_weight[leaving_row] = exact_pivot_weight;
    @memcpy(basis.residual_work, basis.dual_row);
    self.factorization.solve(basis.residual_work) catch return .numerical_failure;
    for (basis.row_edge_weight, basis.pivot_direction, basis.residual_work, 0..) |*weight, alpha, tau, row| {
        if (row == leaving_row) continue;
        const ratio = alpha / pivot;
        const updated = weight.* - 2.0 * ratio * tau + ratio * ratio * exact_pivot_weight;
        if (!std.math.isFinite(updated)) return .numerical_failure;
        weight.* = @max(updated, 1e-4);
    }
    basis.row_edge_weight[leaving_row] = @max(exact_pivot_weight / (pivot * pivot), 1e-4);
    return .optimal;
}

test "exact dual steepest-edge weights use BTRAN row norms" {
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 2, 0);
    engine.basis.?.basic_value[0] = -1.0;
    engine.basis.?.basic_value[1] = -2.0;
    @memset(engine.basis.?.basic_lower, 0.0);
    @memset(engine.basis.?.basic_upper, std.math.inf(f64));
    try engine.factorization.factorize(2, &[_]f64{ 2, 0, 0, 0.5 });
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
    try engine.factorization.factorizeIdentity(2);
    engine.pricing.rule = .steepest_edge;
    engine.dual_edge_weights_valid = true;
    engine.dual_row_index = 0;
    engine.basis.?.dual_row[0] = 1.0;
    engine.basis.?.dual_row[1] = 0.0;
    engine.basis.?.pivot_direction[0] = 2.0;
    engine.basis.?.pivot_direction[1] = 1.0;
    try std.testing.expectEqual(SolveStatus.optimal, engine.updateDualSteepestEdgeWeights(0, 2.0));
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), engine.basis.?.row_edge_weight[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1.25), engine.basis.?.row_edge_weight[1], 1e-12);
}

test "dual steepest-edge recurrence matches exact weights for a nontrivial basis" {
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 2, 0);
    try engine.factorization.factorize(2, &[_]f64{ 2, 1, 1, 3 });
    engine.pricing.rule = .steepest_edge;
    engine.dual_edge_weights_valid = false;
    try std.testing.expectEqual(SolveStatus.optimal, engine.ensureExactDualEdgeWeights());

    @memset(engine.basis.?.dual_row, 0);
    engine.basis.?.dual_row[0] = 1;
    try engine.factorization.solveTranspose(engine.basis.?.dual_row);
    engine.basis.?.pivot_direction[0] = 4;
    engine.basis.?.pivot_direction[1] = -2;
    try engine.factorization.solve(engine.basis.?.pivot_direction);
    const pivot = engine.basis.?.pivot_direction[0];
    try std.testing.expectEqual(SolveStatus.optimal, engine.updateDualSteepestEdgeWeights(0, pivot));
    const updated_weights = [_]f64{
        engine.basis.?.row_edge_weight[0],
        engine.basis.?.row_edge_weight[1],
    };
    try engine.factorization.update(.{
        .leaving_row = 0,
        .entering_col = 0,
        .direction = engine.basis.?.pivot_direction,
        .column_scale = 1,
    });

    for (updated_weights, 0..) |updated, row| {
        @memset(engine.basis.?.rhs_work, 0);
        engine.basis.?.rhs_work[row] = 1;
        try engine.factorization.solveTranspose(engine.basis.?.rhs_work);
        var exact: f64 = 0;
        for (engine.basis.?.rhs_work) |value| exact += value * value;
        try std.testing.expectApproxEqAbs(exact, updated, 1e-12);
    }
}

test "dual Devex updates every row from the hot FTRAN column" {
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 2, 0);
    engine.basis.?.row_edge_weight[0] = 1;
    engine.basis.?.row_edge_weight[1] = 2;
    engine.basis.?.pivot_direction[0] = 2;
    engine.basis.?.pivot_direction[1] = 1;

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
