//! Linear-programming engines and shared LP solve state.

const std = @import("std");

pub const simplex = @import("simplex/root.zig");
pub const presolve = @import("presolve/root.zig");
pub const ProblemView = simplex.problem.ProblemView;
pub const ObjectiveSense = simplex.problem.ObjectiveSense;
pub const SolutionView = simplex.solution.SolutionView;
pub const BasisView = simplex.BasisView;
pub const BasisSnapshot = simplex.BasisSnapshot;

test {
    std.testing.refAllDecls(@This());
}
