//! Numerical tolerances and refactorization policy for simplex iterations.
//!
//! Keeping these values together makes numerical behavior reproducible and
//! allows later adaptive scaling/precision policies without changing the
//! pivot loop.

const std = @import("std");

pub const NumericalState = struct {
    primal_tolerance: f64 = 1e-7,
    dual_tolerance: f64 = 1e-7,
    pivot_tolerance: f64 = 1e-12,
    zero_tolerance: f64 = 1e-12,
    perturbation: f64 = 0.0,
    max_update_count: usize = 100,
    max_refinement_steps: usize = 2,
    residual_tolerance: f64 = 1e-10,
    update_count: usize = 0,
    numerical_warning: bool = false,
    last_relative_residual: f64 = 0.0,
    max_relative_residual: f64 = 0.0,
    refinement_count: usize = 0,
    last_ftran_relative_residual: f64 = 0.0,
    max_ftran_relative_residual: f64 = 0.0,
    /// A cheap pivot-spread warning indicator, not a formal condition number.
    pivot_condition_estimate: f64 = 1.0,
    dual_edge_weight_error_tolerance: f64 = 0.25,
    dual_edge_weight_corrections: usize = 0,
    degenerate_pivot_limit: usize = 8,
    consecutive_degenerate_pivots: usize = 0,
    degenerate_pivot_count: usize = 0,
    anti_cycling_activations: usize = 0,
    anti_cycling_active: bool = false,

    pub fn observePivot(self: *NumericalState, pivot: f64) void {
        self.update_count += 1;
        if (!std.math.isFinite(pivot) or @abs(pivot) <= self.pivot_tolerance) {
            self.numerical_warning = true;
        }
    }

    pub fn needsRefactor(self: NumericalState) bool {
        return self.numerical_warning or self.update_count >= self.max_update_count;
    }

    pub fn markRefactorized(self: *NumericalState) void {
        self.update_count = 0;
        self.numerical_warning = false;
    }

    pub fn observeResidual(self: *NumericalState, absolute: f64, rhs_scale: f64) void {
        self.last_relative_residual = absolute / @max(1.0, rhs_scale);
        self.max_relative_residual = @max(self.max_relative_residual, self.last_relative_residual);
        if (!std.math.isFinite(self.last_relative_residual) or self.last_relative_residual > self.residual_tolerance * 100.0)
            self.numerical_warning = true;
    }

    pub fn isPrimalFeasible(self: NumericalState, violation: f64) bool {
        return violation <= self.primal_tolerance;
    }

    pub fn isDualFeasible(self: NumericalState, violation: f64) bool {
        return violation <= self.dual_tolerance;
    }

    /// Observe primal movement after one simplex iteration. Repeated zero-step
    /// pivots activate deterministic lexicographic tie-breaking; the first
    /// material movement exits the fallback immediately.
    pub fn observeStep(self: *NumericalState, step: f64) void {
        if (std.math.isFinite(step) and step > self.primal_tolerance) {
            self.consecutive_degenerate_pivots = 0;
            self.anti_cycling_active = false;
            self.perturbation = 0.0;
            return;
        }
        self.degenerate_pivot_count += 1;
        self.consecutive_degenerate_pivots += 1;
        if (!self.anti_cycling_active and self.degenerate_pivot_limit != 0 and
            self.consecutive_degenerate_pivots >= self.degenerate_pivot_limit)
        {
            self.anti_cycling_active = true;
            self.anti_cycling_activations += 1;
            self.perturbation = @max(self.zero_tolerance, self.primal_tolerance * 0.1);
        }
    }

    pub fn resetAntiCycling(self: *NumericalState) void {
        self.consecutive_degenerate_pivots = 0;
        self.degenerate_pivot_count = 0;
        self.anti_cycling_activations = 0;
        self.anti_cycling_active = false;
        self.perturbation = 0.0;
        self.last_ftran_relative_residual = 0.0;
        self.max_ftran_relative_residual = 0.0;
    }
};

test "numerical state requests refactor after unstable pivot" {
    var state = NumericalState{ .max_update_count = 2 };
    state.observePivot(0.0);
    try std.testing.expect(state.needsRefactor());
    state.markRefactorized();
    try std.testing.expect(!state.needsRefactor());
}

test "repeated degenerate steps activate and progress clears anti cycling" {
    var state = NumericalState{ .degenerate_pivot_limit = 3 };
    state.observeStep(0.0);
    state.observeStep(1e-12);
    try std.testing.expect(!state.anti_cycling_active);
    state.observeStep(0.0);
    try std.testing.expect(state.anti_cycling_active);
    try std.testing.expectEqual(@as(usize, 1), state.anti_cycling_activations);
    try std.testing.expect(state.perturbation > 0.0);
    state.observeStep(1.0);
    try std.testing.expect(!state.anti_cycling_active);
    try std.testing.expectEqual(@as(usize, 0), state.consecutive_degenerate_pivots);
}
