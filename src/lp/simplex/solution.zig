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
    /// Terminal or early-exit condition produced by the solve.
    status: SolveStatus,
    /// Structural primal values in original model coordinates.
    primal: []const f64,
    /// Row multipliers in original row coordinates and objective convention.
    dual: []const f64,
    /// Reduced costs for structural columns in original coordinates.
    reduced_cost: []const f64,
    /// Improving structural direction when `status == .unbounded`; empty otherwise.
    unbounded_ray: []const f64,
    /// Row-space Farkas certificate when `status == .infeasible`; empty otherwise.
    infeasibility_ray: []const f64,
    /// Final objective including the model offset; meaningful for an optimal solution.
    objective_value: f64,
    /// Number of committed simplex pivots across every phase and cleanup pass.
    iterations: usize,
};
