//! Revised-simplex algorithm family for sparse LPs.
//!
//! The package is intentionally split into policy/state modules so primal and
//! dual simplex, MIP re-optimization, and future parallel variants can share
//! the same basis and factorization infrastructure.

const std = @import("std");

// Submodules. Each owns a distinct concern of the revised simplex machinery.
pub const basis = @import("basis.zig"); // Basis status tracking (basic/nonbasic, bounds)
pub const basis_snapshot = @import("basis_snapshot.zig"); // Immutable snapshots for MIP re-optimization
pub const factorization = @import("factorization.zig"); // LU factorization with Forrest-Tomlin updates
pub const pricing = @import("pricing.zig"); // Entering variable selection policies
pub const ratio_test = @import("ratio_test.zig"); // Leaving variable selection (ratio test)
pub const numerical = @import("numerical.zig"); // Tolerances and refactorization triggers
pub const dual_phase_one = @import("dual_phase_one.zig"); // Dual phase-one construction
pub const crash = @import("crash.zig"); // Advanced basis construction (crash start)
pub const degeneracy = @import("degeneracy.zig"); // Anti-cycling: perturbation + taboo search
pub const pricing_workspace = @import("pricing_workspace.zig"); // Row-oriented CSR view for pricing
pub const problem = @import("problem.zig"); // Borrowed LP input view
pub const solution = @import("solution.zig"); // Borrowed LP solution view
pub const engine = @import("engine.zig"); // Top-level simplex driver

// Convenience re-exports of the most commonly used view types.
pub const ProblemView = problem.ProblemView;
pub const SolutionView = solution.SolutionView;
pub const BasisView = basis_snapshot.BasisView;
pub const BasisSnapshot = basis_snapshot.BasisSnapshot;

test {
    std.testing.refAllDecls(@This());
}
