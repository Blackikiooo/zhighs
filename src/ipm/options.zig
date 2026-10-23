//! Tolerances and run controls for the interior-point engine.
//!
//! Every numerical constant the engine relies on lives here so call sites can
//! audit or override them as a group instead of hunting anonymous literals.

/// Convergence and numerical-safety tolerances. Defaults target the usual
/// 1e-8 LP optimality bar on well-scaled models.
pub const Tolerances = struct {
    /// Relative primal infeasibility ‖b − Ax‖∞ / (1 + ‖b‖∞) required for optimality.
    primal_feasibility: f64 = 1e-8,
    /// Relative dual infeasibility ‖c − Aᵀy − s‖∞ / (1 + ‖c‖∞) required for optimality.
    dual_feasibility: f64 = 1e-8,
    /// Relative duality gap |cᵀx − bᵀy| / (1 + |cᵀx|) required for optimality.
    relative_gap: f64 = 1e-8,
    /// Diagonal shift applied on the first Cholesky regularization retry.
    pivot_regularization: f64 = 1e-12,
    /// Multiplicative growth of the regularization shift between retries.
    regularization_growth: f64 = 100.0,
    /// Largest diagonal shift attempted before declaring numerical failure.
    max_regularization: f64 = 1e-6,
    /// Complementarity growth relative to the starting point that marks divergence.
    divergence_growth: f64 = 1e10,
    /// Variable magnitude beyond which the iterate is treated as escaping to infinity.
    escape_magnitude: f64 = 1e12,
    /// Residual slack factor granted to divergence heuristics relative to the
    /// tight convergence tolerances.
    heuristic_residual_relaxation: f64 = 100.0,
};

pub const RowOrdering = enum {
    automatic,
    natural,
    reverse_cuthill_mckee,
    minimum_degree,
    /// Approximate minimum degree: near-minimum fill at a fraction of the
    /// exact ordering's cost; the default for the automatic policy.
    amd,
};

/// Top-level engine controls.
pub const Options = struct {
    /// Keep the original dense kernel available for validation/benchmarking.
    dense_normal: bool = false,
    /// Apply one-pass row/column max-norm equilibration to standard form.
    equilibrate: bool = false,
    /// Retry a failed unscaled solve with equilibration while preserving the
    /// original total time budget.
    retry_equilibrated_on_failure: bool = true,
    /// Reorder equality rows before sparse symbolic factorization.
    row_ordering: RowOrdering = .automatic,
    time_limit_ns: ?u64 = null,
    interrupt_flag: ?*const @import("std").atomic.Value(bool) = null,
    /// Iteration work budget (one unit per predictor-corrector iteration).
    work_limit: ?u64 = null,
    /// Return true to interrupt before the next iteration.
    iteration_callback: ?*const fn (usize, ?*anyopaque) bool = null,
    callback_context: ?*anyopaque = null,
    /// Convergence and numerical-safety tolerances.
    tolerances: Tolerances = .{},
    /// Hard cap on Mehrotra predictor-corrector iterations.
    max_iterations: usize = 100,
    /// Iterative-refinement corrections after each normal-equation solve.
    refinement_steps: usize = 0,
    /// Gondzio centrality corrections per iteration; each reuses the current
    /// factorization to widen the step, so one is nearly free.
    centrality_corrections: usize = 2,
    /// Feasibility-correction steps when global residuals indicate a candidate
    /// optimum but an individually scaled row still misses tolerance.
    primal_correction_steps: usize = 3,
    /// Fraction-to-boundary backstep applied to primal and dual step lengths.
    step_fraction: f64 = 0.999,
    /// Minimal positive coordinate enforced when lifting the starting point.
    start_margin: f64 = 1e-4,

    pub fn valid(self: Options) bool {
        const std = @import("std");
        inline for (@typeInfo(Tolerances).@"struct".fields) |field| {
            const value = @field(self.tolerances, field.name);
            if (!std.math.isFinite(value) or value <= 0) return false;
        }
        return std.math.isFinite(self.step_fraction) and self.step_fraction > 0 and self.step_fraction < 1 and
            std.math.isFinite(self.start_margin) and self.start_margin > 0 and
            self.tolerances.regularization_growth > 1 and
            self.tolerances.max_regularization >= self.tolerances.pivot_regularization;
    }
};

test {
    const tolerances = Tolerances{};
    try @import("std").testing.expect(tolerances.relative_gap > 0.0);
    const options = Options{};
    try @import("std").testing.expect(options.max_iterations > 0);
}