//! Public control and result types for mixed-integer solves.

const std = @import("std");
const solver = @import("solver");
const plugin = @import("plugin");

pub const Integrality = plugin.Integrality;
pub const ProblemView = plugin.ProblemView;

pub const SolveStatus = enum {
    optimal,
    infeasible,
    unbounded,
    node_limit,
    solution_limit,
    iteration_limit,
    work_limit,
    time_limit,
    interrupted,
    unbounded_or_infeasible,
    invalid_model,
    numerical_failure,
};

/// Component overrides for a single solve. Supplied plugin instances must
/// outlive the call to `solve`. Unset entries fall back to built-in defaults.
pub const SolveComponents = struct {
    branching_rule: ?plugin.BranchingRule = null,
    node_selector: ?plugin.NodeSelector = null,
    heuristics: []const plugin.PrimalHeuristic = &.{},
};

pub const SolveControl = struct {
    node_limit: usize = std.math.maxInt(usize),
    solution_limit: usize = std.math.maxInt(usize),
    integer_tolerance: f64 = 1e-6,
    absolute_gap: f64 = 1e-10,
    relative_gap: f64 = 1e-4,
    initial_solution: ?[]const f64 = null,
    rounding_heuristic: bool = true,
    propagation: bool = true,
    max_propagation_rounds: usize = 8,
    clique_propagation: bool = true,
    max_clique_row_entries: usize = 256,
    max_clique_edges: usize = 1_000_000,
    reduced_cost_fixing: bool = true,
    conflict_learning: bool = true,
    max_conflict_entries: usize = 1_000_000,
    max_conflict_length: usize = 128,
    max_conflict_propagation_rounds: usize = 8,
    conflict_max_age: u16 = 32,
    pseudocost_min_reliable: usize = 4,
    reliability_branching: bool = true,
    strong_branching_max_candidates: usize = 2,
    strong_branching_max_iterations: usize = 20,
    strong_branching_max_work: u64 = 40,
    cut_separation: bool = true,
    max_cut_rounds: usize = 4,
    max_cuts_per_round: usize = 32,
    cut_max_age: u16 = 8,
    lp: solver.LpSolveControl = .{},
    components: SolveComponents = .{},
};

pub const SolveResult = struct {
    allocator: std.mem.Allocator,
    status: SolveStatus,
    primal: []f64 = &.{},
    objective_value: f64 = 0.0,
    objective_bound: f64 = 0.0,
    nodes: usize = 0,
    lp_iterations: usize = 0,
    lp_work: u64 = 0,
    lp_cold_solves: usize = 0,
    lp_reoptimizations: usize = 0,
    pseudocost_observations: usize = 0,
    strong_branching_probes: usize = 0,
    strong_branching_iterations: usize = 0,
    strong_branching_work: u64 = 0,
    cut_rounds: usize = 0,
    cuts_added: usize = 0,
    cuts_expired: usize = 0,
    clique_edges: usize = 0,
    clique_tightenings: usize = 0,
    reduced_cost_tightenings: usize = 0,
    root_reduced_cost_tightenings: usize = 0,
    conflicts_learned: usize = 0,
    conflict_tightenings: usize = 0,
    conflicts_expired: usize = 0,
    conflicts_compacted: usize = 0,
    conflict_checks: usize = 0,
    solution_count: usize = 0,

    pub fn deinit(self: *SolveResult) void {
        self.allocator.free(self.primal);
        self.* = undefined;
    }
};

pub fn validControl(control: SolveControl) bool {
    return std.math.isFinite(control.integer_tolerance) and control.integer_tolerance >= 0 and
        std.math.isFinite(control.absolute_gap) and control.absolute_gap >= 0 and
        std.math.isFinite(control.relative_gap) and control.relative_gap >= 0;
}
