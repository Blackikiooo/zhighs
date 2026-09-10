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
