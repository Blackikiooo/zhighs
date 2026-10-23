//! Model-facing coefficient tolerance shared by LP, presolve, and NLA kernels.
//!
//! This is a statement about the mathematical model, not about matrix storage
//! layout. Keeping it in `foundation` lets `matrix` stay a pure structural
//! layer while every consumer applies the same dropping rule.

const std = @import("std");

/// Canonical LP coefficient floor shared by model-facing sparse kernels.
/// Values at or below this magnitude are numerically zero before scaling;
/// applying the decision in original coordinates keeps assembly, FTRAN,
/// residual checks, and certificate validation structurally consistent.
pub const model_coefficient_tolerance: f64 = 1e-9;

/// Return whether a coefficient survives the shared model-dropping threshold.
pub inline fn retainsModelCoefficient(value: f64) bool {
    return @abs(value) > model_coefficient_tolerance;
}

test "model coefficient policy has a deterministic inclusive floor" {
    try std.testing.expect(!retainsModelCoefficient(model_coefficient_tolerance));
    try std.testing.expect(!retainsModelCoefficient(-model_coefficient_tolerance));
    try std.testing.expect(retainsModelCoefficient(2.0 * model_coefficient_tolerance));
}
