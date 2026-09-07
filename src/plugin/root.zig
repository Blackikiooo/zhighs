//! Runtime component contracts.
//!
//! Contracts cover constraint handlers, presolvers, separators, propagators,
//! branching rules, primal heuristics, conflict handlers, node selectors, cut
//! selectors, variable pricers, Benders components, IIS finders, event
//! handlers, readers, and optional relaxation handlers. The types currently
//! materialised are the ones the MIP branch-and-bound engine dispatches
//! through; further contracts follow the same fat-pointer shape.

const std = @import("std");

pub const Kind = @import("kind.zig").Kind;
const contracts = @import("contracts.zig");

pub const Integrality = contracts.Integrality;
pub const ProblemView = contracts.ProblemView;
pub const BranchingKind = contracts.BranchingKind;
pub const BranchingCandidate = contracts.BranchingCandidate;
pub const BranchingContext = contracts.BranchingContext;
pub const BranchingRule = contracts.BranchingRule;
pub const bindBranching = contracts.bindBranching;
pub const HeuristicContext = contracts.HeuristicContext;
pub const PrimalHeuristic = contracts.PrimalHeuristic;
pub const bindHeuristic = contracts.bindHeuristic;
pub const OpenNode = contracts.OpenNode;
pub const NodeSelectionContext = contracts.NodeSelectionContext;
pub const NodeSelector = contracts.NodeSelector;
pub const bindNodeSelector = contracts.bindNodeSelector;

test {
    std.testing.refAllDecls(@This());
}
