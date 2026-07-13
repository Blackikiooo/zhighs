//! Borrowed simplex solution view.
//!
//! All slices reference engine-owned SoA storage. A `SolutionView` becomes
//! invalid when the engine is re-solved or deinitialized.

pub const SolveStatus = enum {
    optimal,
    infeasible,
    unbounded,
    iteration_limit,
    work_limit,
    time_limit,
    interrupted,
    numerical_failure,
    not_implemented,
};

pub const SolutionView = struct {
    status: SolveStatus,
    primal: []const f64,
    dual: []const f64,
    reduced_cost: []const f64,
    objective_value: f64,
    iterations: usize,
};
