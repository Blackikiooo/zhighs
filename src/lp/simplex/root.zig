//! Revised-simplex algorithm family for sparse LPs.
//!
//! The package is intentionally split into policy/state modules so primal and
//! dual simplex, MIP re-optimization, and future parallel variants can share
//! the same basis and factorization infrastructure.

const std = @import("std");

/// Mutable basis membership, solution vectors and pivot workspaces.
pub const basis = @import("basis.zig");
/// Borrowed and owning basis representations used by warm starts.
pub const basis_snapshot = @import("basis_snapshot.zig");
/// Dense/sparse basis factorization and update policy.
pub const factorization = @import("factorization.zig");
/// Primal and dual entering/leaving selection policies.
pub const pricing = @import("pricing.zig");
/// Primal and dual ratio-test implementations.
pub const ratio_test = @import("ratio_test.zig");
/// Numerical tolerances, residual tracking and reinversion policy.
pub const numerical = @import("numerical.zig");
/// Persistent dual work arrays and Phase-I bound/cost transformations.
pub const dual_phase_one = @import("dual_phase_one.zig");
/// Preferred semantic alias for the persistent serial-dual work module.
pub const dual_work = dual_phase_one;
/// Unified dual phase, rebuild and cost-lifecycle controller.
pub const dual_state = @import("dual_state.zig");
/// Initial basis construction strategies.
pub const crash = @import("crash.zig");
/// Perturbation, taboo and anti-cycling state.
pub const degeneracy = @import("degeneracy.zig");
/// Reusable row-oriented matrix view for pricing.
pub const pricing_workspace = @import("pricing_workspace.zig");
/// Borrowed LP input representation.
pub const problem = @import("problem.zig");
/// Borrowed, engine-owned solve result representation.
pub const solution = @import("solution.zig");
/// Top-level revised-simplex controller.
pub const engine = @import("engine.zig");

/// Convenience alias for the borrowed LP input view.
pub const ProblemView = problem.ProblemView;
/// Convenience alias for the borrowed solve result.
pub const SolutionView = solution.SolutionView;
/// Convenience alias for a borrowed warm-start basis.
pub const BasisView = basis_snapshot.BasisView;
/// Convenience alias for an owning warm-start basis snapshot.
pub const BasisSnapshot = basis_snapshot.BasisSnapshot;

test {
    std.testing.refAllDecls(@This());
}
