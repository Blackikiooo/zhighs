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
    /// A cheap pivot-spread warning indicator, not a formal condition number.
    pivot_condition_estimate: f64 = 1.0,

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
};

test "numerical state requests refactor after unstable pivot" {
    var state = NumericalState{ .max_update_count = 2 };
    state.observePivot(0.0);
    try std.testing.expect(state.needsRefactor());
    state.markRefactorized();
    try std.testing.expect(!state.needsRefactor());
}
