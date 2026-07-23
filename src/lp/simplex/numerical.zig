//! Numerical tolerances and refactorization policy for simplex iterations.
//!
//! Keeping these values together makes numerical behavior reproducible and
//! allows later adaptive scaling/precision policies without changing the
//! pivot loop.

const std = @import("std");

/// Mutable numerical state shared across simplex iterations. Tracks pivot
/// stability, refactorization triggers, and anti-cycling fallbacks.
pub const NumericalState = struct {
    // --- Feasibility and pivot tolerances ---
    /// Maximum accepted violation of a primal bound.
    primal_tolerance: f64 = 1e-7,
    /// Maximum accepted sign violation of a nonbasic reduced cost.
    dual_tolerance: f64 = 1e-7,
    /// Pivot magnitudes at or below this threshold are treated as singular.
    pivot_tolerance: f64 = 1e-12,
    /// Generic structural-zero threshold for residual and scalar comparisons.
    zero_tolerance: f64 = 1e-12,
    /// Active anti-degeneracy perturbation; zero outside fallback mode.
    perturbation: f64 = 0.0,

    // --- Refactorization policy ---
    /// Maximum Forrest–Tomlin updates allowed before reinversion.
    max_update_count: usize = 100,
    /// Maximum correction solves attempted by iterative refinement.
    max_refinement_steps: usize = 2,
    /// Number of fresh-basis pivots required after a forward-accuracy warning
    /// before Forrest--Tomlin updates are tried again.
    fresh_factorization_recovery_pivots: usize = 32,
    /// Target relative residual used by solve validation.
    residual_tolerance: f64 = 1e-10,

    // --- Live diagnostics (reset by `markRefactorized`) ---
    /// Basis updates observed since the last successful factorization.
    update_count: usize = 0,
    /// Latched request for reinversion after an accuracy warning.
    numerical_warning: bool = false,
    /// Relative residual from the most recently observed generic solve.
    last_relative_residual: f64 = 0.0,
    /// Largest generic relative residual observed during the solve.
    max_relative_residual: f64 = 0.0,
    /// Total refinement correction steps performed during the solve.
    refinement_count: usize = 0,
    /// Relative residual from the most recent entering-column FTRAN.
    last_ftran_relative_residual: f64 = 0.0,
    /// Largest entering-column FTRAN residual observed during the solve.
    max_ftran_relative_residual: f64 = 0.0,
    /// A cheap pivot-spread warning indicator, not a formal condition number.
    pivot_condition_estimate: f64 = 1.0,
    /// Relative error threshold for replacing a recurrence DSE weight by its exact value.
    dual_edge_weight_error_tolerance: f64 = 0.25,
    /// Number of exact DSE weight corrections applied.
    dual_edge_weight_corrections: usize = 0,

    // --- Anti-cycling fallback ---
    /// Consecutive zero-length steps required to activate anti-cycling.
    degenerate_pivot_limit: usize = 8,
    /// Current run length of degenerate steps.
    consecutive_degenerate_pivots: usize = 0,
    /// Total degenerate pivots observed in the current solve.
    degenerate_pivot_count: usize = 0,
    /// Number of times anti-cycling mode was entered.
    anti_cycling_activations: usize = 0,
    /// Whether deterministic anti-cycling choices are currently required.
    anti_cycling_active: bool = false,

    /// Account for one pivot; flag a warning if it is non-finite or below tolerance.
    pub fn observePivot(self: *NumericalState, pivot: f64) void {
        self.update_count += 1;
        if (!std.math.isFinite(pivot) or @abs(pivot) <= self.pivot_tolerance) {
            self.numerical_warning = true;
        }
    }

    /// True when the factorization must be rebuilt before the next FTran/BTran.
    pub fn needsRefactor(self: NumericalState) bool {
        return self.numerical_warning or self.update_count >= self.max_update_count;
    }

    /// Reset all update/refactor counters after a fresh factorization.
    pub fn markRefactorized(self: *NumericalState) void {
        self.update_count = 0;
        self.numerical_warning = false;
    }

    /// Record a residual observation from a solve/refinement step.
    pub fn observeResidual(self: *NumericalState, absolute: f64, rhs_scale: f64) void {
        self.last_relative_residual = absolute / @max(1.0, rhs_scale);
        self.max_relative_residual = @max(self.max_relative_residual, self.last_relative_residual);
        if (!std.math.isFinite(self.last_relative_residual) or self.last_relative_residual > self.residual_tolerance * 100.0)
            self.numerical_warning = true;
    }

    /// Check primal feasibility of a violation.
    pub fn isPrimalFeasible(self: NumericalState, violation: f64) bool {
        return violation <= self.primal_tolerance;
    }

    /// Check dual feasibility of a violation.
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

    /// Full reset of anti-cycling counters (called between solves).
    pub fn resetAntiCycling(self: *NumericalState) void {
        self.consecutive_degenerate_pivots = 0;
        self.degenerate_pivot_count = 0;
        self.anti_cycling_activations = 0;
        self.anti_cycling_active = false;
        self.perturbation = 0.0;
        self.last_ftran_relative_residual = 0.0;
        self.max_ftran_relative_residual = 0.0;
    }

    /// Leave the current fallback after a phase/objective transition while
    /// preserving solve-level diagnostics accumulated so far.
    pub fn clearAntiCyclingFallback(self: *NumericalState) void {
        self.consecutive_degenerate_pivots = 0;
        self.anti_cycling_active = false;
        self.perturbation = 0.0;
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
