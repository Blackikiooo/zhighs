//! Borrowed simplex solution view.
//!
//! All slices reference engine-owned SoA storage. A `SolutionView` becomes
//! invalid when the engine is re-solved or deinitialized.

/// Outcome of a simplex solve. Values are ordered roughly by severity:
/// success -> early-exit -> failure.
pub const SolveStatus = enum {
    optimal, // KKT conditions satisfied within tolerance
    infeasible, // No feasible point exists for the supplied bounds
    unbounded, // Objective improves without bound along a primal ray
    iteration_limit, // Stopped after exhausting the configured iteration budget
    work_limit, // Stopped after exhausting the configured floating-point work budget
    time_limit, // Stopped after exhausting the configured wall-clock budget
    interrupted, // Stopped in response to an external cancellation request
    numerical_failure, // Refactorization could not recover accuracy; aborting
    not_implemented, // Requested path is not yet supported by the engine
};

/// Read-only view of the solver result. All slices point at engine storage
/// and remain valid until the next `solve*` call or `deinit`.
pub const SolutionView = struct {
    status: SolveStatus,
    primal: []const f64, // Primal variable values (length = num_cols)
    dual: []const f64, // Row dual values (length = num_rows)
    reduced_cost: []const f64, // Column reduced costs (length = num_cols)
    unbounded_ray: []const f64, // Primal ray when `status == .unbounded`; empty otherwise
    infeasibility_ray: []const f64, // Row-space Farkas ray when `status == .infeasible`
    objective_value: f64, // Final objective value (including offset)
    iterations: usize, // Number of simplex pivots performed
};
