//! Public facade for the modular mixed-integer solver.
//!
//! Implementation details are separated by ownership and hot-path role:
//! `search` orchestrates branch-and-bound, `domain` normalizes variable
//! domains, `propagation` tightens node bounds, `node` owns search nodes,
//! `incumbent` validates primal candidates, and `bounds` centralizes objective
//! and gap policy.

const types = @import("types.zig");
const search = @import("search.zig");

pub const Integrality = types.Integrality;
pub const ProblemView = types.ProblemView;
pub const SolveStatus = types.SolveStatus;
pub const SolveComponents = types.SolveComponents;
pub const SolveControl = types.SolveControl;
pub const SolveResult = types.SolveResult;

pub const solve = search.solve;

test {
    _ = search;
    _ = @import("bounds.zig");
    _ = @import("domain.zig");
    _ = @import("incumbent.zig");
    _ = @import("node.zig");
    _ = @import("node_queue.zig");
    _ = @import("propagation.zig");
    _ = @import("clique_table.zig");
    _ = @import("reduced_cost_fixing.zig");
    _ = @import("root_reduced_costs.zig");
    _ = @import("conflict_pool.zig");
    _ = @import("pseudocost.zig");
    _ = @import("strong_branching.zig");
    _ = @import("cut_pool.zig");
    _ = @import("cut_relaxation.zig");
    _ = @import("separation.zig");
    _ = @import("cut_generation.zig");
    _ = @import("relaxation.zig");
}
