//! Revised-simplex orchestration boundary.
//!
//! `SimplexEngine` owns all mutable solve state and reusable workspaces but
//! borrows the problem matrix through `ProblemView`. Algorithm implementations
//! are split into `engine_*.zig` files; their methods share this struct so the
//! hot path does not allocate or copy large state objects.
const std = @import("std");
const basis_module = @import("basis.zig");
const basis_snapshot_module = @import("basis_snapshot.zig");
const factorization_module = @import("factorization.zig");
const pricing_module = @import("pricing.zig");
const ratio_module = @import("ratio_test.zig");
const numerical_module = @import("numerical.zig");
const dual_phase_one_module = @import("dual_phase_one.zig");
const dual_state_module = @import("dual_state.zig");
const crash_module = @import("crash.zig");
const degeneracy_module = @import("degeneracy.zig");
const pricing_workspace_module = @import("pricing_workspace.zig");
const problem_module = @import("problem.zig");
const solution_module = @import("solution.zig");
const foundation = @import("foundation");
const matrix = @import("matrix");

/// Revised-simplex direction used for the active phase.
pub const Algorithm = enum { primal_revised, dual_revised };
/// Cold-start method used when the logical basis is not immediately feasible.
pub const PhaseOneStrategy = enum { primal, dual, automatic };
/// Initial-basis construction policy.
pub const CrashStrategy = enum { logical, ltssf, bixby, automatic };
/// Degeneracy and cycling mitigation policy.
pub const DegeneracyStrategy = enum { baseline, perturbation, perturbation_taboo, automatic };
/// Pricing rule installed during primal Phase I.
pub const PhaseOnePricingStrategy = enum { inherit, dantzig, devex, steepest_edge };
/// Matrix traversal used to construct pricing results.
pub const PricingKernel = enum { column, row, automatic };
/// Primal Devex weight-maintenance implementation.
pub const DevexStrategy = enum { legacy, framework };
/// Optional segmented pricing policy for primal simplex.
pub const PrimalPricingStrategy = enum { inherit, partial };
/// Dual edge-weight lifecycle. `inherit` preserves the caller's pricing rule;
/// `steepest_devex` starts with exact DSE and deterministically falls back to
/// full dual Devex when the recurrence is rejected or exceeds its work budget.
pub const DualEdgeWeightStrategy = enum { inherit, steepest_devex };
/// Dual initialization formulas used before Phase I/II selection.
pub const DualInitializationStrategy = enum { baseline, highs };
/// Public progress phase reported to callbacks and traces.
pub const SolvePhase = enum { phase_one, dual_phase_one, dual_feasibility_repair, phase_two };
/// Return value by which a progress callback requests continuation or stop.
pub const CallbackAction = enum { continue_solve, stop };
/// Borrowed scalar progress snapshot. It owns no memory and is valid only for
/// the callback invocation.
pub const ProgressEventView = struct {
    /// Phase executing when the event was emitted.
    phase: SolvePhase,
    /// Primal or dual revised algorithm currently performing pivots.
    algorithm: Algorithm,
    /// Attempted-iteration counter visible to the public solve.
    iterations: usize,
    /// Deterministic work units consumed so far.
    work_used: u64,
    /// Current internally scaled objective estimate.
    objective_value: f64,
    /// Largest current primal bound violation.
    primal_infeasibility: f64,
    /// Largest current reduced-cost sign violation.
    dual_infeasibility: f64,
};
/// User callback invoked at configured work intervals; storage is caller-owned.
pub const IterationCallback = *const fn (event: ProgressEventView, user_data: ?*anyopaque) CallbackAction;
/// Informational callback that cannot stop the solve.
pub const IterationLogCallback = *const fn (event: ProgressEventView, user_data: ?*anyopaque) void;
/// Built-in iteration logging level.
pub const LogLevel = enum { off, iterations };
/// Public solve result status re-exported by the engine API.
pub const SolveStatus = solution_module.SolveStatus;
/// Internal stage at which an unrecovered numerical failure occurred.
pub const FailureSite = enum {
    none,
    direction_column,
    direction_solve,
    direction_refactor,
    direction_refinement,
    pivot_update,
    pivot_factorization,
    pivot_edge_weights,
    ratio_test,
    primal_step,
    dual_feasibility,
    reduced_cost,
    optimality_check,
};

/// Diagnostic exit route from the shifted-cost dual Phase-II path.
pub const ShiftedDualExit = enum {
    none,
    running,
    setup_free_infeasibility,
    setup_not_dual_feasible,
    phase_two_no_entering,
    phase_two_fresh_infeasible,
    phase_two_numerical_failure,
    phase_two_stopped,
    original_dual_feasible,
    original_dual_infeasible,
    original_dual_phase_two_optimal,
    original_dual_phase_two_failed,
    cleanup_primal_optimal,
    cleanup_primal_failed,
    cleanup_neither_feasible,
};
/// Reason validation of an infeasibility certificate was rejected.
pub const InfeasibilityCertificateFailure = enum {
    none,
    invalid_workspace,
    nonfinite_ray,
    zero_ray,
    infinite_row_bound,
    infinite_column_bound,
    nonfinite_proof,
    nonpositive_gap,
};
/// Deterministic record of one committed pivot for differential testing.
pub const PivotTraceEvent = struct {
    /// Phase in which the pivot was committed.
    phase: SolvePhase,
    /// Global attempted-iteration number.
    iteration: usize,
    /// Internal column entering the basis.
    entering_column: u32,
    /// Internal column leaving the basis.
    leaving_column: u32,
    /// Basis row replaced by the entering column.
    leaving_row: u32,
    /// Accepted tableau pivot after direction normalization.
    pivot: f64,
    /// Nonnegative primal movement applied during the pivot.
    step: f64,
    /// Factor updates present immediately before/after trace publication.
    update_count: usize,
    /// Relative residual measured for the entering-column FTRAN.
    ftran_relative_residual: f64,
    /// Current inexpensive pivot-spread stability estimate.
    condition_estimate: f64,
    /// Number of nonbasic bound flips associated with this iteration.
    bound_flip_count: usize,
};
/// Mutually exclusive primary cause assigned to a degenerate iteration.
pub const DegeneracyReason = enum {
    bound_tie,
    ratio_tie,
    zero_primal_step,
    phase_one_objective_stall,
    repeated_basis,
    small_pivot_retry,
    bound_flip,
};
/// Allocation-free diagnostic emitted only for a numerically degenerate
/// iteration. Storage is owned by the solve caller.
pub const DegeneracyTraceEvent = struct {
    /// Phase in which degeneracy was detected.
    phase: SolvePhase,
    /// Global attempted-iteration number.
    iteration: usize,
    /// Primary classification used for statistics.
    reason: DegeneracyReason,
    /// Internal entering column, or the sentinel chosen by the producer.
    entering_column: u32,
    /// Internal leaving column, or the sentinel chosen by the producer.
    leaving_column: u32,
    /// Primal movement associated with the event.
    step: f64,
    /// Objective movement associated with the event.
    objective_change: f64,
    /// Hash of the basis head used to detect repeated bases.
    basis_fingerprint: u64,
};
/// Classification assigned to one CHUZC candidate in failure diagnostics.
pub const DualPhaseOneCandidateReason = enum {
    small_tableau,
    basic_or_fixed,
    wrong_pivot_sign,
    accepted_bound_flip,
    eligible_unselected,
};
/// Snapshot of one dual Phase-I CHUZC candidate at a no-entering failure.
pub const DualPhaseOneCandidateTraceEvent = struct {
    /// Internal candidate column.
    column: u32,
    /// Basis/bound status at the time of CHUZC.
    status: basis_module.BasisStatus,
    /// Explicit persistent nonbasic movement in {-1,0,+1}.
    explicit_move: i8,
    /// Pivotal tableau coefficient.
    tableau: f64,
    /// Temporary CHUZC movement direction.
    direction: f64,
    /// Tableau coefficient after leaving/movement sign normalization.
    signed_pivot: f64,
    /// Working reduced cost.
    reduced_cost: f64,
    /// Working lower bound.
    lower: f64,
    /// Working upper bound.
    upper: f64,
    /// Working nonbasic primal value.
    primal: f64,
    /// Amount of leaving-row violation removable by flipping this column.
    flip_capacity: f64,
    /// Reason the candidate was rejected, flipped or left unselected.
    reason: DualPhaseOneCandidateReason,
};
/// One nonzero of the pivotal BTRAN row captured for diagnostics.
pub const DualPhaseOneEpTraceEvent = struct {
    /// Basis-row index of this retained BTRAN entry.
    row: u32,
    /// Numerical value of `e_p^T B^-1` at `row`.
    value: f64,
};
/// Aggregate context recorded when dual Phase I cannot choose an entering column.
pub const DualPhaseOneFailureDiagnostic = struct {
    /// Attempted iteration on which CHUZC failed.
    iteration: usize,
    /// Selected leaving basis row.
    leaving_row: u32,
    /// Internal basic column occupying `leaving_row`.
    leaving_column: u32,
    /// Working bound violated by the leaving variable.
    leaving_bound: basis_module.BasisStatus,
    /// Magnitude of the leaving variable's bound violation.
    violation: f64,
    /// Working value of the leaving variable.
    leaving_value: f64,
    /// Phase-I lower bound of the leaving variable.
    working_lower: f64,
    /// Phase-I upper bound of the leaving variable.
    working_upper: f64,
    /// Original LP lower bound of the leaving variable.
    original_lower: f64,
    /// Original LP upper bound of the leaving variable.
    original_upper: f64,
    /// Number of retained nonzeros in the pivotal BTRAN row.
    ep_nonzeros: usize,
    /// Largest absolute pivotal-row coefficient.
    ep_max_abs: f64,
    /// Candidates rejected for a coefficient below pivot tolerance.
    small_tableau: usize,
    /// Candidates rejected because they were basic or fixed.
    basic_or_fixed: usize,
    /// Candidates rejected by the required pivot sign.
    wrong_pivot_sign: usize,
    /// Candidates consumed as accepted BFRT flips.
    accepted_bound_flips: usize,
    /// Eligible candidates not selected before failure.
    eligible_unselected: usize,
};
/// Result of the primal ratio-test wrapper, including early failure status.
pub const PrimalLeavingResult = struct {
    /// `.optimal` means `row`/`step` contain a valid choice.
    status: SolveStatus,
    /// Leaving basis row, or null when the ray is unbounded.
    row: ?u32 = null,
    /// Maximum feasible primal movement.
    step: f64 = 0.0,
    /// Bound at which the leaving variable becomes nonbasic.
    bound: basis_module.BasisStatus = .at_lower,
};
/// Errors exposed while validating and installing a caller-provided basis.
pub const BasisImportError = basis_snapshot_module.BasisViewError || error{
    InvalidNonbasicStatus,
    SingularBasis,
    NumericalFailure,
};
/// Caller-owned policy, limits, callbacks and diagnostic buffers for one solve.
pub const SolveControl = struct {
    /// Maximum attempted simplex iterations across all internal phases.
    max_iterations: usize = 1_000_000,
    /// `null` disables the wall-clock limit; zero requests immediate stop.
    time_limit_ns: ?u64 = null,
    /// Optional caller-owned atomic flag. The flag must outlive `solveProblem`.
    interrupt_flag: ?*const std.atomic.Value(bool) = null,
    /// Optional borrowed warm-start basis.
    initial_basis: ?basis_snapshot_module.BasisView = null,
    /// Optional clock provider for deterministic tests or custom runtimes.
    /// The value is borrowed only for the duration of a solve call.
    clock_io: ?std.Io = null,
    /// Deterministic simplex work budget. One unit is charged for each
    /// attempted primal, dual, or Phase-I iteration.
    work_limit: ?u64 = null,
    /// Progress callback invoked according to `callback_interval_work`.
    iteration_callback: ?IterationCallback = null,
    /// Opaque pointer forwarded unchanged to `iteration_callback`.
    callback_user_data: ?*anyopaque = null,
    /// Deterministic work units between progress callbacks; zero is normalized.
    callback_interval_work: u64 = 1,
    /// Built-in logging mode.
    log_level: LogLevel = .off,
    /// Optional external iteration-log sink.
    log_callback: ?IterationLogCallback = null,
    /// Opaque pointer forwarded unchanged to `log_callback`.
    log_user_data: ?*anyopaque = null,
    /// Deterministic work units between log events.
    log_interval_work: u64 = 100,
    /// Optional caller-owned trace storage used for deterministic pivot-path
    /// differential tests. Events beyond the supplied capacity are dropped.
    pivot_trace: []PivotTraceEvent = &.{},
    /// Optional caller-owned storage for mutually exclusive degeneracy causes.
    /// Events beyond the supplied capacity are dropped while counters remain exact.
    degeneracy_trace: []DegeneracyTraceEvent = &.{},
    /// Optional caller-owned storage populated only when dual Phase I has no
    /// entering column. These buffers are diagnostic and never affect pricing.
    dual_phase_one_candidate_trace: []DualPhaseOneCandidateTraceEvent = &.{},
    /// Optional caller-owned storage for nonzeros of the failing BTRAN row.
    dual_phase_one_ep_trace: []DualPhaseOneEpTraceEvent = &.{},
    /// Collect phase and kernel timings. Disabled by default so production
    /// solves do not perform clock reads in simplex hot paths.
    collect_statistics: bool = false,
    /// Cold-start Phase-I policy. `primal` remains the compatibility default
    /// until corpus A/B validation allows deterministic automatic selection.
    phase_one_strategy: PhaseOneStrategy = .primal,
    /// Cold-start basis construction. Kept logical by default until the
    /// explicit LTSSF path passes its performance gate.
    crash_strategy: CrashStrategy = .logical,
    /// Optional structural-column cap for deterministic crash-prefix A/B.
    /// Null uses the longest numerically valid prefix within the logical-anchor limit.
    crash_max_columns: ?usize = null,
    /// Degenerate-pivot policy. Automatic enables bounded perturbation after
    /// the anti-cycling trigger; terminal validation cold-restarts baseline if
    /// the perturbed epoch cannot publish a valid status or certificate.
    degeneracy_strategy: DegeneracyStrategy = .automatic,
    /// Pricing rule override used only during primal Phase I.
    phase_one_pricing: PhaseOnePricingStrategy = .inherit,
    /// Enable drift-triggered exact reduced-cost repricing.
    adaptive_reprice: bool = false,
    /// Select once per pricing operation; no representation branch exists in
    /// either coefficient inner loop.
    pricing_kernel: PricingKernel = .column,
    /// Full reference-framework Devex is the default as of T4 (2026-07-20).
    /// 90/93 Stage 7 optimal, 88 common-completion models zero regression,
    /// and the framework→legacy cold restart provides a deterministic safety
    /// net. Explicit `.legacy` remains available for A/B isolation.
    devex_strategy: DevexStrategy = .framework,
    /// Explicit segmented primal pricing A/B policy. `inherit` preserves the
    /// engine pricing rule selected by the embedding application.
    primal_pricing_strategy: PrimalPricingStrategy = .inherit,
    /// Explicit dual DSE -> Devex A/B policy. Default is steepest-devex
    /// (HiGHS-equivalent: SteepestEdge with Devex fallback).
    dual_edge_weight_strategy: DualEdgeWeightStrategy = .steepest_devex,
    /// Maximum successful DSE recurrence updates per dual phase before
    /// switching to dual Devex. Zero disables the budget-triggered switch.
    dual_dse_update_budget: usize = 64,
    /// Source-alignment track for HiGHS Phase-I coordinates and initialization.
    /// Baseline remains the production default until the complete dual corpus
    /// clears the no-regression gate.
    dual_initialization_strategy: DualInitializationStrategy = .baseline,
};

/// Optional counters and timings accumulated during one solve.
pub const SimplexStats = struct {
    /// Attempted primal Phase-I iterations.
    phase_one_iterations: usize = 0,
    /// Attempted dual Phase-I iterations.
    dual_phase_one_iterations: usize = 0,
    /// Dual Phase-I paths that fell back to primal Phase I.
    dual_phase_one_fallbacks: usize = 0,
    /// Transactional basis restores followed by a dual Phase-I retry.
    dual_phase_one_snapshot_retries: usize = 0,
    /// Iterations executed with temporary shifted costs.
    shifted_dual_iterations: usize = 0,
    /// Cleanup iterations after removing shifted costs.
    shifted_cleanup_iterations: usize = 0,
    /// Sparse crash-basis construction attempts.
    crash_attempts: usize = 0,
    /// Crash attempts rejected in favor of the logical basis.
    crash_fallbacks: usize = 0,
    /// Structural columns considered by the crash planner.
    crash_planned_columns: usize = 0,
    /// Structural columns installed in the accepted crash basis.
    crash_structural_columns: usize = 0,
    /// Nonzeros in the accepted crash basis.
    crash_basis_nonzeros: usize = 0,
    /// Stability estimate of the accepted crash factorization.
    crash_condition_estimate: f64 = 1.0,
    /// Iterations spent repairing a warm basis.
    dual_repair_iterations: usize = 0,
    /// Attempted Phase-II iterations.
    phase_two_iterations: usize = 0,
    /// Nanoseconds spent in primal Phase I.
    phase_one_ns: u64 = 0,
    /// Nanoseconds spent in dual Phase I.
    dual_phase_one_ns: u64 = 0,
    /// Nanoseconds spent repairing warm-basis feasibility.
    dual_repair_ns: u64 = 0,
    /// Nanoseconds spent in Phase II.
    phase_two_ns: u64 = 0,
    /// Nanoseconds spent in terminal cleanup.
    cleanup_ns: u64 = 0,
    /// Complete basis rebuild pipelines executed.
    rebuild_calls: usize = 0,
    /// Nanoseconds spent in basis rebuild pipelines.
    rebuild_ns: u64 = 0,
    /// Pricing operations executed.
    pricing_calls: usize = 0,
    /// Nanoseconds spent in pricing operations.
    pricing_ns: u64 = 0,
    /// Vectors sampled for pricing-density statistics.
    pricing_samples: usize = 0,
    /// Retained nonzeros across sampled pricing vectors.
    pricing_nonzeros: usize = 0,
    /// Total entries across sampled pricing vectors.
    pricing_entries: usize = 0,
    /// Individual nonbasic columns flipped by BFRT.
    bound_flips: usize = 0,
    /// BFRT batches containing at least one flip.
    bound_flip_batches: usize = 0,
    /// FTRAN calls avoided by aggregating each flip batch.
    bound_flip_ftran_savings: usize = 0,
    /// Incremental dual reduced-cost updates.
    dual_reduced_cost_updates: usize = 0,
    /// Exact dual reduced-cost recomputations.
    dual_exact_reprices: usize = 0,
    /// Successful dual steepest-edge recurrence updates.
    dual_dse_updates: usize = 0,
    /// CHUZR candidates rejected after BTRAN exposed an updated DSE weight
    /// below one quarter of its exact value.
    dual_dse_weight_rejections: usize = 0,
    /// Successful dual Devex recurrence updates.
    dual_devex_updates: usize = 0,
    /// DSE-to-Devex switches caused by invalid recurrence state.
    dual_dse_invalid_fallbacks: usize = 0,
    /// DSE-to-Devex switches caused by the configured update budget.
    dual_dse_budget_fallbacks: usize = 0,
    /// Primal Devex reference frameworks initialized.
    devex_frameworks: usize = 0,
    /// Pivots updating a primal Devex framework.
    devex_framework_updates: usize = 0,
    /// Invalid primal Devex weights observed.
    devex_bad_weights: usize = 0,
    /// Degenerate iterations classified as bound ties.
    degeneracy_bound_ties: usize = 0,
    /// Degenerate iterations classified as ratio ties.
    degeneracy_ratio_ties: usize = 0,
    /// Degenerate iterations classified as zero primal steps.
    degeneracy_zero_primal_steps: usize = 0,
    /// Degenerate iterations classified as Phase-I objective stalls.
    degeneracy_phase_one_objective_stalls: usize = 0,
    /// Degenerate iterations classified as repeated bases.
    degeneracy_repeated_bases: usize = 0,
    /// Degenerate iterations requiring a small-pivot retry.
    degeneracy_small_pivot_retries: usize = 0,
    /// Degenerate iterations associated with a bound flip.
    degeneracy_bound_flips: usize = 0,
    /// Perturbation epochs activated.
    perturbation_activations: usize = 0,
    /// Perturbation epochs expired by policy.
    perturbation_expirations: usize = 0,
    /// Cleanup passes performed after perturbation.
    perturbation_cleanups: usize = 0,
    /// Internal cold restarts performed.
    cold_restart_solves: usize = 0,
    /// Cold restarts returning specifically to Phase I.
    cold_restart_phase_one: usize = 0,
    /// Basis/taboo entries recorded.
    taboo_records: usize = 0,
    /// Candidate choices retried due to taboo state.
    taboo_retries: usize = 0,
    /// Entering FTRAN vectors sampled for density.
    aq_samples: usize = 0,
    /// Nonzeros across sampled entering FTRAN vectors.
    aq_nonzeros: usize = 0,
    /// Pivotal BTRAN rows sampled for density.
    ep_samples: usize = 0,
    /// Nonzeros across sampled pivotal BTRAN rows.
    ep_nonzeros: usize = 0,
    /// Pricing operations dispatched to dense traversal.
    dense_pricing_dispatches: usize = 0,
    /// Pricing operations dispatched to hyper-sparse traversal.
    hyper_pricing_dispatches: usize = 0,
    /// Pricing operations using the row-oriented matrix.
    row_pricing_dispatches: usize = 0,
    /// Pricing operations using CSC column traversal.
    column_pricing_dispatches: usize = 0,
    /// Rebuilds requested to set up Phase I.
    rebuild_phase_one_setup: usize = 0,
    /// Rebuilds requested by solve residual validation.
    rebuild_solve_residual: usize = 0,
    /// Rebuilds requested after a small pivot.
    rebuild_small_pivot: usize = 0,
    /// Rebuilds requested by factor update-count limit.
    rebuild_update_limit: usize = 0,
    /// Rebuilds requested by factor growth.
    rebuild_update_growth: usize = 0,
    /// Rebuilds requested after entering-direction refinement.
    rebuild_direction_refinement: usize = 0,
    /// Rebuilds requested while operating in fresh-factor mode.
    rebuild_fresh_mode: usize = 0,
    /// Rebuilds requested because a factor update was rejected.
    rebuild_update_rejected: usize = 0,
    /// Rebuilds requested for original-objective cleanup.
    rebuild_cleanup: usize = 0,
    /// Rebuilds requested to initialize/reset edge weights.
    rebuild_edge_weight_reset: usize = 0,
    /// Rebuilds requested by the general numerical policy.
    rebuild_numerical_policy: usize = 0,
    /// Classification of the first factor-update failure.
    first_update_failure_kind: ?factorization_module.UpdateFailureKind = null,
    /// Iteration at which the first factor-update failure occurred.
    first_update_failure_iteration: usize = 0,
    /// Entering column involved in the first factor-update failure.
    first_update_failure_entering: u32 = 0,
    /// Leaving row involved in the first factor-update failure.
    first_update_failure_leaving_row: u32 = 0,

    /// Sum rebuild counters classified by reason.
    pub fn classifiedRebuilds(self: SimplexStats) usize {
        return self.rebuild_phase_one_setup + self.rebuild_solve_residual + self.rebuild_small_pivot +
            self.rebuild_update_limit + self.rebuild_update_growth + self.rebuild_direction_refinement +
            self.rebuild_fresh_mode + self.rebuild_update_rejected + self.rebuild_cleanup +
            self.rebuild_edge_weight_reset + self.rebuild_numerical_policy;
    }

    /// Sum degeneracy counters classified by mutually exclusive reason.
    pub fn classifiedDegeneratePivots(self: SimplexStats) usize {
        return self.degeneracy_bound_ties + self.degeneracy_ratio_ties +
            self.degeneracy_zero_primal_steps + self.degeneracy_phase_one_objective_stalls +
            self.degeneracy_repeated_bases + self.degeneracy_small_pivot_retries +
            self.degeneracy_bound_flips;
    }
};

/// Monotonic accounting for one public solve. Unlike `iterations`, these
/// counters survive every internal fallback, cleanup solve, and cold restart.
/// HiGHS' simplex iteration count corresponds to committed basis changes;
/// bound-only moves are therefore reported separately.
pub const IterationCounters = struct {
    /// Iteration bodies entered, including iterations that do not commit a basis change.
    attempted_iterations: usize = 0,
    /// Total committed basis changes across all phases.
    committed_pivots: usize = 0,
    /// Nonbasic bound changes that do not alter basis membership.
    bound_moves: usize = 0,
    /// Pivots performed while temporary shifted costs are installed.
    shifted_dual_pivots: usize = 0,
    /// Pivots performed by dual Phase I.
    dual_phase_one_pivots: usize = 0,
    /// Pivots performed by dual Phase II with original/working costs.
    dual_phase_two_pivots: usize = 0,
    /// Pivots performed by artificial primal Phase I.
    primal_phase_one_pivots: usize = 0,
    /// Pivots performed by primal Phase II.
    primal_phase_two_pivots: usize = 0,
    /// Pivots performed while repairing a warm basis for dual feasibility.
    dual_repair_pivots: usize = 0,
    /// Pivots performed during terminal cleanup.
    cleanup_pivots: usize = 0,

    /// Sum pivots assigned to a concrete phase; should equal `committed_pivots`.
    pub fn classifiedPivots(self: IterationCounters) usize {
        return self.shifted_dual_pivots + self.dual_phase_one_pivots +
            self.dual_phase_two_pivots + self.primal_phase_one_pivots +
            self.primal_phase_two_pivots + self.dual_repair_pivots +
            self.cleanup_pivots;
    }
};

/// Reason the numerical factorization layer requests a basis rebuild.
pub const RebuildReason = enum {
    phase_one_setup,
    solve_residual,
    small_pivot,
    update_limit,
    update_growth,
    direction_refinement,
    fresh_mode,
    update_rejected,
    cleanup,
    edge_weight_reset,
    numerical_policy,
};

/// Owning revised-simplex solver instance and reusable workspace.
pub const SimplexEngine = struct {
    /// Allocator owning the basis and every reusable sub-workspace.
    allocator: std.mem.Allocator,
    /// Revised-simplex direction currently executing.
    algorithm: Algorithm = .dual_revised,
    /// Model-shaped mutable basis, allocated at the start of a solve.
    basis: ?basis_module.BasisState = null,
    /// Current basis inverse/factor representation and update chain.
    factorization: factorization_module.Factorization,
    /// Persistent dual Phase-I/II bounds, moves, shifts and diagnostics.
    dual_work: dual_phase_one_module.DualWorkState,
    /// Unified dual phase and rebuild controller.
    dual_control: dual_state_module.ControlState = .{},
    /// Reusable sparse crash-basis planner.
    crash: crash_module.CrashWorkspace,
    /// Perturbation, taboo and repeated-basis workspace.
    degeneracy: degeneracy_module.Workspace,
    /// Optional CSR companion used by row-oriented pricing.
    pricing_row_view: pricing_workspace_module.RowView,
    /// Pricing rule and partial-pricing cursor state.
    pricing: pricing_module.Pricing = .{},
    /// Active primal/dual ratio-test policy.
    ratio_test: ratio_module.RatioTest = .{},
    /// Tolerances, residual observations and reinversion policy.
    numerical: numerical_module.NumericalState = .{},
    /// Attempted iterations in the current internal solve path.
    iterations: usize = 0,
    /// Monotonic public-solve pivot accounting across fallback paths.
    iteration_counters: IterationCounters = .{},
    /// Recursion depth used to distinguish public solves from internal fallbacks.
    solve_depth: usize = 0,
    /// Current scaled objective estimate or final published objective.
    objective_value: f64 = 0.0,
    /// Multiplier mapping working objective coefficients to original scale.
    objective_scale: f64 = 1.0,
    /// Whether `basis.unbounded_ray` contains a validated certificate.
    unbounded_ray_valid: bool = false,
    /// Whether `basis.infeasibility_ray` contains a validated certificate.
    infeasibility_ray_valid: bool = false,
    /// Positive separation established by the most recent Farkas validation.
    infeasibility_certificate_gap: f64 = 0.0,
    /// First reason the most recent Farkas validation failed.
    infeasibility_certificate_failure: InfeasibilityCertificateFailure = .none,
    /// Ray mass touching infinite row bounds during validation.
    infeasibility_certificate_infinite_row_mass: f64 = 0.0,
    /// Ray mass touching infinite column bounds during validation.
    infeasibility_certificate_infinite_column_mass: f64 = 0.0,
    /// Whether artificial primal Phase I is required after initialization.
    phase1_needed: bool = false,
    /// Monotonic start timestamp, absent when no time source is configured.
    solve_start_ns: ?i96 = null,
    /// Borrowed clock backend used for stop checks.
    solve_clock_io: ?std.Io = null,
    /// Deterministic work units consumed by the public solve.
    work_used: u64 = 0,
    /// Whether row dual edge weights match the current basis.
    dual_edge_weights_valid: bool = false,
    /// Basis row whose `dual_row`/tableau workspace is currently valid.
    dual_row_index: ?u32 = null,
    /// Active entries in the hyper-sparse dual leaving candidate list.
    dual_candidate_count: usize = 0,
    /// Score cutoff used to retain hyper-sparse dual candidates.
    dual_candidate_cutoff: f64 = 0.0,
    /// Whether CHUZR currently uses its sparse candidate representation.
    dual_hyper_sparse_active: bool = false,
    /// Structural basic columns replaced by logical columns while repairing a
    /// singular imported basis during the current solve.
    rank_repair_count: usize = 0,
    /// Number of pivot events written to caller trace storage.
    pivot_trace_count: usize = 0,
    /// Borrowed caller storage active for the current solve.
    active_pivot_trace: []PivotTraceEvent = &.{},
    /// Number of degeneracy events written to caller storage.
    degeneracy_trace_count: usize = 0,
    /// Borrowed caller degeneracy trace storage.
    active_degeneracy_trace: []DegeneracyTraceEvent = &.{},
    /// Fixed-size recent-basis hash ring for repeated-basis detection.
    degeneracy_basis_fingerprints: [32]u64 = @splat(0),
    /// Number of initialized entries in the basis hash ring.
    degeneracy_basis_fingerprint_count: usize = 0,
    /// Next hash-ring slot to overwrite.
    degeneracy_basis_fingerprint_cursor: usize = 0,
    /// Last captured dual Phase-I no-entering diagnostic.
    dual_phase_one_failure: ?DualPhaseOneFailureDiagnostic = null,
    /// Number of candidate diagnostic events written.
    dual_phase_one_candidate_trace_count: usize = 0,
    /// Borrowed candidate trace storage for the current solve.
    active_dual_phase_one_candidate_trace: []DualPhaseOneCandidateTraceEvent = &.{},
    /// Number of pivotal-row diagnostic entries written.
    dual_phase_one_ep_trace_count: usize = 0,
    /// Borrowed pivotal-row trace storage for the current solve.
    active_dual_phase_one_ep_trace: []DualPhaseOneEpTraceEvent = &.{},
    /// Boxed columns flipped during initial dual infeasibility correction.
    dual_phase_one_initial_flips: usize = 0,
    /// Working dual objective immediately after Phase-I initialization.
    dual_phase_one_initial_objective: f64 = 0.0,
    /// Largest basic-bound violation immediately after Phase-I initialization.
    dual_phase_one_initial_max_basic_violation: f64 = 0.0,
    /// First leaving row selected by dual Phase I, if any.
    dual_phase_one_initial_leaving_row: ?u32 = null,
    /// Public progress phase currently reported by callbacks.
    current_phase: SolvePhase = .phase_two,
    /// Entering-column refinement requested reinversion instead of FT update.
    direction_requires_reinversion: bool = false,
    /// Pivots still forced to use a fresh factorization after an accuracy warning.
    fresh_factorization_pivots_remaining: usize = 0,
    /// Incremental reduced-cost updates since the last exact reprice.
    reduced_cost_update_count: usize = 0,
    /// Fixed/adaptive maximum incremental updates between exact reprices.
    reduced_cost_refresh_period: usize = 8,
    /// Internal stage of the most recent unrecovered failure.
    failure_site: FailureSite = .none,
    /// Exit route taken by shifted-cost dual Phase II.
    shifted_dual_exit: ShiftedDualExit = .none,
    /// Failure stage saved before shifted-path fallback overwrites it.
    shifted_dual_failure_site: FailureSite = .none,
    /// Solve statistics; timing fields remain zero when collection is disabled.
    stats: SimplexStats = .{},
    /// Borrowed clock backend used only for optional statistics.
    statistics_io: ?std.Io = null,
    /// Whether the current primal/dual pass is original-objective cleanup.
    cleanup_active: bool = false,
    /// Resolved degeneracy policy for the current solve.
    active_degeneracy_strategy: DegeneracyStrategy = .baseline,
    /// Resolved reduced-cost drift repricing policy.
    active_adaptive_reprice: bool = false,
    /// Largest reduced-cost drift observed at an exact reprice.
    maximum_reduced_cost_drift: f64 = 0.0,
    /// Number of exact reprices performed in the current solve.
    exact_reprices: usize = 0,
    /// Resolved row/column pricing traversal.
    active_pricing_kernel: PricingKernel = .column,
    /// Resolved primal Devex implementation.
    active_devex_strategy: DevexStrategy = .legacy,
    /// Resolved segmented primal pricing policy.
    active_primal_pricing_strategy: PrimalPricingStrategy = .inherit,
    /// Resolved dual DSE-to-Devex lifecycle.
    active_dual_edge_weight_strategy: DualEdgeWeightStrategy = .inherit,
    /// Resolved DSE update budget before Devex fallback.
    active_dual_dse_update_budget: usize = 64,
    /// Resolved dual initialization formulas.
    active_dual_initialization_strategy: DualInitializationStrategy = .baseline,
    /// Successful DSE recurrences since the current dual phase began.
    dual_dse_updates_since_start: usize = 0,
    /// Pivots performed in the current primal Devex reference framework.
    devex_framework_iterations: usize = 0,
    /// Invalid/nonpositive Devex weights observed in the current framework.
    devex_bad_weight_count: usize = 0,
    /// Request to rebuild the Devex reference set after committing the pivot.
    devex_reset_after_pivot: bool = false,

    /// Construct an empty engine; model-shaped storage is allocated lazily.
    pub fn init(a: std.mem.Allocator) SimplexEngine {
        return .{
            .allocator = a,
            .factorization = factorization_module.Factorization.init(a),
            .dual_work = dual_phase_one_module.DualWorkState.init(a),
            .crash = crash_module.CrashWorkspace.init(a),
            .degeneracy = degeneracy_module.Workspace.init(a),
            .pricing_row_view = pricing_workspace_module.RowView.init(a),
        };
    }
    /// Release the active basis and every retained sub-workspace.
    pub fn deinit(self: *SimplexEngine) void {
        if (self.basis) |*b| b.deinit();
        self.factorization.deinit();
        self.dual_work.deinit();
        self.crash.deinit();
        self.degeneracy.deinit();
        self.pricing_row_view.deinit();
    }

    /// Return requested bytes held by retained engine workspaces.
    ///
    /// Allocator metadata and borrowed problem storage are excluded.
    pub fn requestedBytes(self: *const SimplexEngine) usize {
        const basis_bytes = if (self.basis) |*basis| basis.requestedBytes() else 0;
        return basis_bytes + self.factorization.requestedBytes() + self.dual_work.requestedBytes() +
            self.crash.requestedBytes() + self.degeneracy.requestedBytes() + self.pricing_row_view.requestedBytes();
    }
    /// Legacy dimension-only entry point retained for API compatibility.
    /// Use `solveProblem` for an actual LP.
    pub fn solve(_: *SimplexEngine, _: usize, _: usize, _: SolveControl) SolveStatus {
        return .not_implemented;
    }

    /// Solve entry point consuming a borrowed LP `ProblemView`.
    ///
    /// The view must outlive this call; the engine never takes ownership of
    /// model arrays. Basis storage is owned by the engine instance.
    pub fn solveProblem(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
        const root_solve = self.solve_depth == 0;
        self.solve_depth += 1;
        defer self.solve_depth -= 1;
        if (root_solve) self.iteration_counters = .{};
        if (root_solve) self.dual_control.reset();
        self.pricing.resetPartial();
        self.startSolveClock(control);
        self.resetStatistics(control);
        self.work_used = 0;
        self.rank_repair_count = 0;
        self.fresh_factorization_pivots_remaining = 0;
        self.reduced_cost_update_count = 0;
        self.failure_site = .none;
        self.shifted_dual_exit = .none;
        self.shifted_dual_failure_site = .none;
        self.unbounded_ray_valid = false;
        self.infeasibility_ray_valid = false;
        self.infeasibility_certificate_gap = 0.0;
        self.infeasibility_certificate_failure = .none;
        self.infeasibility_certificate_infinite_row_mass = 0.0;
        self.infeasibility_certificate_infinite_column_mass = 0.0;
        self.dual_phase_one_initial_flips = 0;
        self.dual_phase_one_initial_objective = 0.0;
        self.dual_phase_one_initial_max_basic_violation = 0.0;
        self.dual_phase_one_initial_leaving_row = null;
        self.numerical.resetAntiCycling();
        self.degeneracy.resetSolve();
        self.active_degeneracy_strategy = control.degeneracy_strategy;
        self.active_adaptive_reprice = control.adaptive_reprice;
        self.active_pricing_kernel = control.pricing_kernel;
        self.active_devex_strategy = control.devex_strategy;
        self.active_primal_pricing_strategy = control.primal_pricing_strategy;
        self.active_dual_edge_weight_strategy = control.dual_edge_weight_strategy;
        self.active_dual_dse_update_budget = control.dual_dse_update_budget;
        self.active_dual_initialization_strategy = control.dual_initialization_strategy;
        self.dual_dse_updates_since_start = 0;
        self.reduced_cost_refresh_period = 8;
        self.maximum_reduced_cost_drift = 0.0;
        self.exact_reprices = 0;
        if (self.controlledStop(control)) |status| return status;
        problem.validate() catch |err| return switch (err) {
            error.InvalidBounds => .infeasible,
            error.DimensionMismatch, error.InvalidMatrix => .numerical_failure,
        };
        if (self.active_pricing_kernel != .column)
            self.pricing_row_view.build(problem.matrix) catch return .numerical_failure;
        if (self.basis) |*old| old.deinit();
        self.basis = basis_module.BasisState.init(self.allocator, problem.num_rows, problem.num_cols) catch return .numerical_failure;
        self.degeneracy.ensureCapacity(problem.num_rows, self.basis.?.col_status.len) catch return .numerical_failure;
        if (self.initializeProblemStorage(problem) != .optimal) return .infeasible;
        self.factorization.factorizeIdentity(problem.num_rows) catch return .numerical_failure;
        self.dual_edge_weights_valid = true;
        self.observeFactorizationStability();
        self.iterations = 0;

        if (control.initial_basis) |initial_basis| {
            if (self.importBasis(problem, initial_basis)) |_| {
                if (self.recomputeReducedCosts(problem) != .optimal) return .numerical_failure;
                const feasibility = self.classifyFeasibility(problem);
                if (feasibility.primal) {
                    self.algorithm = .primal_revised;
                    return self.solvePrimal(problem, control);
                }
                if (feasibility.dual) {
                    self.algorithm = .dual_revised;
                    return self.solveDual(problem, control);
                }
                const repair_status = self.repairWarmBasisWithDual(problem, control);
                if (repair_status == .optimal) {
                    if (self.recomputeReducedCosts(problem) != .optimal) return .numerical_failure;
                    const repaired = self.classifyFeasibility(problem);
                    if (repaired.primal) {
                        self.algorithm = .primal_revised;
                        return self.solvePrimal(problem, control);
                    }
                    if (repaired.dual) {
                        self.algorithm = .dual_revised;
                        return self.solveDual(problem, control);
                    }
                } else if (repair_status == .work_limit or repair_status == .time_limit or
                    repair_status == .iteration_limit or repair_status == .interrupted)
                {
                    return repair_status;
                }
            } else |_| {}
            // Invalid, singular, or neither-feasible warm starts fall back to
            // the deterministic logical crash basis.
            if (self.initializeProblemStorage(problem) != .optimal) return .infeasible;
            self.factorization.factorizeIdentity(problem.num_rows) catch return .numerical_failure;
            self.dual_edge_weights_valid = true;
            self.observeFactorizationStability();
        }

        self.algorithm = .primal_revised;
        const logical_numerical_state = self.numerical;
        const logical_pricing_state = self.pricing;
        var crash_installed = false;
        const crash_strategy = switch (control.crash_strategy) {
            .logical => CrashStrategy.logical,
            .ltssf => CrashStrategy.ltssf,
            .bixby => CrashStrategy.bixby,
            // Corpus A/B rejects both sparse crash policies as a universal
            // default: Bixby benefits scsd1 but still regresses brandy. Keep
            // automatic pinned to logical until a broader policy clears it.
            .automatic => CrashStrategy.logical,
        };
        if (crash_strategy == .ltssf or crash_strategy == .bixby) {
            self.stats.crash_attempts += 1;
            crash_installed = self.installSparseCrashBasis(problem, control.crash_max_columns, if (crash_strategy == .bixby) .bixby else .ltssf);
            if (!crash_installed) {
                self.stats.crash_fallbacks += 1;
                if (self.initializeProblemStorage(problem) != .optimal) return .infeasible;
                self.factorization.factorizeIdentity(problem.num_rows) catch return .numerical_failure;
                self.numerical = logical_numerical_state;
                self.pricing = logical_pricing_state;
                self.dual_edge_weights_valid = true;
                self.observeFactorizationStability();
            }
        }
        if (!crash_installed) {
            if (self.initializeLogicalBasicValues(problem) != .optimal) return .numerical_failure;
        }
        if (self.recomputeReducedCosts(problem) != .optimal) return .numerical_failure;
        const crash_feasibility = self.classifyFeasibility(problem);
        // An explicitly requested dual path must not silently become primal
        // simplex merely because the logical basis is primal feasible. HiGHS
        // forced-dual still establishes dual feasibility and enters dual
        // Phase II; preserving the shortcut made the iteration comparison an
        // algorithm-vs-algorithm mismatch (notably blend).
        if (crash_feasibility.primal and control.phase_one_strategy != .dual)
            return self.solvePrimal(problem, control);
        // The artificial primal Phase I columns are tied to the logical
        // identity. A structural crash must use the general dual Phase I;
        // failure restores the logical basis before the primal fallback.
        const phase_one_strategy = if (crash_installed)
            PhaseOneStrategy.dual
        else switch (control.phase_one_strategy) {
            .primal => PhaseOneStrategy.primal,
            .dual => PhaseOneStrategy.dual,
            .automatic => self.chooseColdPhaseOneStrategy(problem),
        };
        if (phase_one_strategy == .dual) {
            const numerical_before_dual_phase_one = if (crash_installed) logical_numerical_state else self.numerical;
            const pricing_before_dual_phase_one = if (crash_installed) logical_pricing_state else self.pricing;
            const checkpoint_basis = if (self.basis) |*basis| basis else return .numerical_failure;
            self.dual_work.captureBasisCheckpoint(
                checkpoint_basis,
                problem.num_cols + problem.num_rows,
                problem.num_rows,
            ) catch return .numerical_failure;
            if (self.active_dual_initialization_strategy == .baseline) {
                const shifted_dual_status = self.solveDualWithCostShifts(problem, control);
                switch (shifted_dual_status) {
                    .optimal, .work_limit, .time_limit, .iteration_limit, .interrupted => return shifted_dual_status,
                    .infeasible => if (self.infeasibility_ray_valid) return .infeasible,
                    .unbounded, .not_implemented, .numerical_failure => {},
                }
            }
            const dual_phase_one_status = self.solveDualPhaseOne(problem, control);
            if (dual_phase_one_status != .not_implemented and dual_phase_one_status != .numerical_failure)
                return dual_phase_one_status;
            // The capacity-guarded shifted path may reach a different basis
            // from the pre-A1 path. If its Phase I fails, restore the exact
            // pre-shift basis transactionally and give the original dual
            // Phase I one bounded retry before the primal cold fallback.
            if (self.dual_work.checkpoint_valid) retry: {
                self.numerical = numerical_before_dual_phase_one;
                self.pricing = pricing_before_dual_phase_one;
                self.failure_site = .none;
                self.direction_requires_reinversion = false;
                self.reduced_cost_update_count = 0;
                self.dual_edge_weights_valid = true;
                const checkpoint_view = basis_snapshot_module.BasisView{
                    .structural_status = self.dual_work.checkpoint_status[0..problem.num_cols],
                    .logical_status = self.dual_work.checkpoint_status[problem.num_cols..][0..problem.num_rows],
                    .basic_index = self.dual_work.checkpoint_basic_index[0..problem.num_rows],
                };
                self.importBasis(problem, checkpoint_view) catch break :retry;
                if (self.recomputeReducedCosts(problem) != .optimal) break :retry;
                self.stats.dual_phase_one_snapshot_retries += 1;
                const retry_status = self.solveDualPhaseOne(problem, control);
                if (retry_status != .not_implemented and retry_status != .numerical_failure)
                    return retry_status;
            }
            if (crash_installed) self.stats.crash_fallbacks += 1;
            self.stats.dual_phase_one_fallbacks += 1;
            // A failed explicit experiment falls back to the frozen logical
            // crash and artificial Phase I, preserving the correctness path.
            if (self.initializeProblemStorage(problem) != .optimal) return .infeasible;
            self.factorization.factorizeIdentity(problem.num_rows) catch return .numerical_failure;
            self.numerical = numerical_before_dual_phase_one;
            self.pricing = pricing_before_dual_phase_one;
            self.failure_site = .none;
            self.direction_requires_reinversion = false;
            self.fresh_factorization_pivots_remaining = 0;
            self.reduced_cost_update_count = 0;
            self.dual_candidate_count = 0;
            self.dual_candidate_cutoff = 0.0;
            self.dual_hyper_sparse_active = false;
            self.dual_row_index = null;
            self.cleanup_active = false;
            self.dual_edge_weights_valid = true;
            self.observeFactorizationStability();
            if (self.initializeLogicalBasicValues(problem) != .optimal) return .numerical_failure;
            self.iterations = 0;
        }
        self.algorithm = .primal_revised;
        if (self.installArtificialPhaseOneBasis(problem) != .optimal) return .numerical_failure;
        if (self.phase1_needed and self.refactorizeBasis(problem, .phase_one_setup) != .optimal) return .numerical_failure;
        self.numerical.markRefactorized();
        self.iterations = 0;
        if (self.phase1_needed) {
            const phase1_status = self.solvePhaseOne(problem, control);
            if (phase1_status != .optimal) return phase1_status;
        }
        if (self.controlledStop(control)) |status| return status;
        return self.solvePrimal(problem, control);
    }

    // Method re-exports from the engine_*.zig implementation files.
    pub const initializeLogicalBasicValues = @import("engine_setup.zig").initializeLogicalBasicValues;
    pub const installArtificialPhaseOneBasis = @import("engine_setup.zig").installArtificialPhaseOneBasis;
    pub const installSparseCrashBasis = @import("engine_setup.zig").installSparseCrashBasis;
    pub const refreshProblemStorage = @import("engine_setup.zig").refreshProblemStorage;
    pub const initializeProblemStorage = @import("engine_setup.zig").initializeProblemStorage;

    pub const importBasis = @import("engine_basis.zig").importBasis;
    pub const exportBasisView = @import("engine_basis.zig").exportBasisView;
    pub const exportBasisSnapshot = @import("engine_basis.zig").exportBasisSnapshot;
    pub const classifyFeasibility = @import("engine_basis.zig").classifyFeasibility;

    pub const computeDirection = @import("engine_pivot.zig").computeDirection;
    pub const directionResidualAcceptable = @import("engine_pivot.zig").directionResidualAcceptable;
    pub const performPivot = @import("engine_pivot.zig").performPivot;
    pub const recomputeReducedCostsFromWork = @import("engine_pivot.zig").recomputeReducedCostsFromWork;

    pub const controlledStop = @import("engine_progress.zig").controlledStop;
    pub const beginIteration = @import("engine_progress.zig").beginIteration;
    pub const progressEvent = @import("engine_progress.zig").progressEvent;
    pub const startSolveClock = @import("engine_progress.zig").startSolveClock;
    pub const resetStatistics = @import("engine_progress.zig").resetStatistics;
    pub const statisticsTimestamp = @import("engine_progress.zig").statisticsTimestamp;
    pub const elapsedSince = @import("engine_progress.zig").elapsedSince;
    pub const recordPhaseElapsed = @import("engine_progress.zig").recordPhaseElapsed;
    pub const recordPhaseIterations = @import("engine_progress.zig").recordPhaseIterations;
    pub const recordPricingElapsed = @import("engine_progress.zig").recordPricingElapsed;
    pub const recordRowPricingElapsed = @import("engine_progress.zig").recordRowPricingElapsed;
    pub const selectRowPricing = @import("engine_progress.zig").selectRowPricing;
    pub const observeAqDensity = @import("engine_progress.zig").observeAqDensity;
    pub const observeEpDensity = @import("engine_progress.zig").observeEpDensity;
    pub const observePricingDensity = @import("engine_progress.zig").observePricingDensity;

    pub const recomputeReducedCostsWithDrift = @import("engine_pivot.zig").recomputeReducedCostsWithDrift;
    pub const updateReducedCostsAfterPrimalPivot = @import("engine_pivot.zig").updateReducedCostsAfterPrimalPivot;
    pub const recomputeReducedCosts = @import("engine_pivot.zig").recomputeReducedCosts;
    pub const refactorizeBasis = @import("engine_pivot.zig").refactorizeBasis;
    pub const factorizeCurrentBasis = @import("engine_pivot.zig").factorizeCurrentBasis;
    pub const finishRefactorization = @import("engine_pivot.zig").finishRefactorization;

    pub const repairRankDeficientBasis = @import("engine_basis.zig").repairRankDeficientBasis;
    pub const nonbasicStatusForBounds = @import("engine_basis.zig").nonbasicStatusForBounds;

    pub const observeFactorizationStability = @import("engine_pivot.zig").observeFactorizationStability;

    pub const fillInternalColumn = @import("engine_setup.zig").fillInternalColumn;

    pub const recomputeBasicValues = @import("engine_pivot.zig").recomputeBasicValues;
    pub const recomputeBasicValuesUnchecked = @import("engine_pivot.zig").recomputeBasicValuesUnchecked;
    pub const recomputeBasicValuesImpl = @import("engine_pivot.zig").recomputeBasicValuesImpl;
    pub const subtractBasisProduct = @import("engine_pivot.zig").subtractBasisProduct;
    pub const subtractBasisProductWithMagnitude = @import("engine_pivot.zig").subtractBasisProductWithMagnitude;

    pub const finishOptimal = @import("engine_finish.zig").finishOptimal;
    pub const finishUnbounded = @import("engine_finish.zig").finishUnbounded;
    pub const validateUnboundedRay = @import("engine_finish.zig").validateUnboundedRay;
    pub const validateOptimalSolution = @import("engine_finish.zig").validateOptimalSolution;
    pub const solutionView = @import("engine_finish.zig").solutionView;

    /// Reoptimize with the existing basis and factorization. The caller must
    /// guarantee unchanged matrix structure and values. No model view is
    /// retained after this call.
    pub fn reoptimizeProblem(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
        const root_solve = self.solve_depth == 0;
        self.solve_depth += 1;
        defer self.solve_depth -= 1;
        if (root_solve) self.iteration_counters = .{};
        self.pricing.resetPartial();
        self.startSolveClock(control);
        self.resetStatistics(control);
        self.work_used = 0;
        self.rank_repair_count = 0;
        self.fresh_factorization_pivots_remaining = 0;
        self.numerical.resetAntiCycling();
        self.degeneracy.resetSolve();
        self.active_degeneracy_strategy = switch (control.degeneracy_strategy) {
            .automatic => .automatic,
            .baseline => .baseline,
            .perturbation => .perturbation,
            .perturbation_taboo => .perturbation_taboo,
        };
        self.active_adaptive_reprice = control.adaptive_reprice;
        self.active_pricing_kernel = control.pricing_kernel;
        self.active_devex_strategy = control.devex_strategy;
        self.active_primal_pricing_strategy = control.primal_pricing_strategy;
        self.active_dual_edge_weight_strategy = control.dual_edge_weight_strategy;
        self.active_dual_dse_update_budget = control.dual_dse_update_budget;
        self.dual_dse_updates_since_start = 0;
        self.reduced_cost_refresh_period = 8;
        self.maximum_reduced_cost_drift = 0.0;
        self.exact_reprices = 0;
        if (self.controlledStop(control)) |status| return status;
        problem.validate() catch |err| return switch (err) {
            error.InvalidBounds => .infeasible,
            error.DimensionMismatch, error.InvalidMatrix => .numerical_failure,
        };
        if (self.active_pricing_kernel != .column)
            self.pricing_row_view.build(problem.matrix) catch return .numerical_failure;
        if (!self.refreshProblemStorage(problem)) {
            var cold_control = control;
            cold_control.initial_basis = null;
            return self.solveProblem(problem, cold_control);
        }
        self.iterations = 0;
        if (self.recomputeBasicValuesUnchecked(problem) != .optimal or self.recomputeReducedCosts(problem) != .optimal)
            return .numerical_failure;
        const feasibility = self.classifyFeasibility(problem);
        if (feasibility.primal) {
            self.algorithm = .primal_revised;
            return self.solvePrimal(problem, control);
        }
        if (feasibility.dual) {
            self.algorithm = .dual_revised;
            return self.solveDual(problem, control);
        }
        const repair_status = self.repairWarmBasisWithDual(problem, control);
        if (repair_status == .optimal) {
            if (self.recomputeReducedCosts(problem) != .optimal) return .numerical_failure;
            if (self.classifyFeasibility(problem).primal) {
                self.algorithm = .primal_revised;
                return self.solvePrimal(problem, control);
            }
        } else if (repair_status == .work_limit or repair_status == .time_limit or
            repair_status == .iteration_limit or repair_status == .interrupted)
        {
            return repair_status;
        }
        var cold_control = control;
        cold_control.initial_basis = null;
        return self.solveProblem(problem, cold_control);
    }

    // Method re-exports from the engine_*.zig implementation files.
    pub const chooseLeaving = @import("engine_primal.zig").chooseLeaving;
    pub const chooseLeavingWithPolicy = @import("engine_primal.zig").chooseLeavingWithPolicy;
    pub const chooseLeavingBland = @import("engine_primal.zig").chooseLeavingBland;
    pub const pivotNeedsFreshFactorization = @import("engine_primal.zig").pivotNeedsFreshFactorization;
    pub const updateLegacyDevexWeights = @import("engine_primal.zig").updateLegacyDevexWeights;
    pub const initializePrimalDevexFramework = @import("engine_primal.zig").initializePrimalDevexFramework;
    pub const updatePrimalDevexFramework = @import("engine_primal.zig").updatePrimalDevexFramework;
    pub const solvePrimal = @import("engine_primal.zig").solvePrimal;
    pub const scaledDualTolerance = @import("engine_primal.zig").scaledDualTolerance;
    pub const choosePrimalEnteringTimed = @import("engine_primal.zig").choosePrimalEnteringTimed;
    pub const choosePrimalEnteringWeightedTimed = @import("engine_primal.zig").choosePrimalEnteringWeightedTimed;
    pub const chooseColdPhaseOneStrategy = @import("engine_primal.zig").chooseColdPhaseOneStrategy;
    pub const solvePhaseOne = @import("engine_primal.zig").solvePhaseOne;
    pub const restartPhaseOneWithoutPerturbation = @import("engine_primal.zig").restartPhaseOneWithoutPerturbation;

    pub const restartSolveWithoutPerturbation = @import("engine_degeneracy.zig").restartSolveWithoutPerturbation;
    pub const countPrimalRatioTies = @import("engine_degeneracy.zig").countPrimalRatioTies;
    pub const prepareDegeneracyPolicy = @import("engine_degeneracy.zig").prepareDegeneracyPolicy;
    pub const observeIterationStep = @import("engine_degeneracy.zig").observeIterationStep;
    pub const basisFingerprint = @import("engine_degeneracy.zig").basisFingerprint;
    pub const rememberDegenerateBasis = @import("engine_degeneracy.zig").rememberDegenerateBasis;

    pub const cleanupArtificialBasis = @import("engine_primal.zig").cleanupArtificialBasis;
    pub const recomputePhaseOneReducedCosts = @import("engine_primal.zig").recomputePhaseOneReducedCosts;
    pub const recomputePhaseOneReducedCostsWithDrift = @import("engine_primal.zig").recomputePhaseOneReducedCostsWithDrift;

    // Method re-exports from the engine_*.zig implementation files.
    pub const updateDualDevexWeights = @import("engine_dual_edge_weight.zig").updateDualDevexWeights;

    pub const solveDual = @import("engine_dual.zig").solveDual;
    pub const solveDualWithCostShifts = @import("engine_dual.zig").solveDualWithCostShifts;
    pub const buildAndValidateInfeasibilityRay = @import("engine_dual.zig").buildAndValidateInfeasibilityRay;
    pub const solveDualPhaseOne = @import("engine_dual.zig").solveDualPhaseOne;
    pub const recordDualPhaseOneNoEntering = @import("engine_dual.zig").recordDualPhaseOneNoEntering;
    pub const buildDualPhaseOneCosts = @import("engine_dual.zig").buildDualPhaseOneCosts;
    pub const dualCandidate = @import("engine_dual_candidate.zig").dualCandidate;
    pub const dualCandidateScore = @import("engine_dual_candidate.zig").dualCandidateScore;
    pub const bestDualCandidate = @import("engine_dual_candidate.zig").bestDualCandidate;
    pub const rebuildDualCandidateList = @import("engine_dual_candidate.zig").rebuildDualCandidateList;

    pub const ensureExactDualEdgeWeights = @import("engine_dual_edge_weight.zig").ensureExactDualEdgeWeights;
    pub const beginDualEdgeWeightPhase = @import("engine_dual_edge_weight.zig").beginDualEdgeWeightPhase;
    pub const switchDualDseToDevex = @import("engine_dual_edge_weight.zig").switchDualDseToDevex;
    pub const updateDualSteepestEdgeWeights = @import("engine_dual_edge_weight.zig").updateDualSteepestEdgeWeights;

    pub const repairWarmBasisWithDual = @import("engine_dual.zig").repairWarmBasisWithDual;
    pub const computeDualTableauRow = @import("engine_dual.zig").computeDualTableauRow;
    pub const applyBoundFlips = @import("engine_dual.zig").applyBoundFlips;
    pub const boundFlipResidualAcceptable = @import("engine_dual.zig").boundFlipResidualAcceptable;
    pub const accumulateBoundFlipRhs = @import("engine_dual.zig").accumulateBoundFlipRhs;
    pub const updateReducedCostsAfterDualPivot = @import("engine_dual.zig").updateReducedCostsAfterDualPivot;
};

test {
    std.testing.refAllDecls(@This());
}

test "engine solves with a nonzero structural lower bound" {
    const rows = [_]foundation.RowId{foundation.RowId.fromUsizeAssumeValid(0)};
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{1},
        .col_lower = &[_]f64{2},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{-std.math.inf(f64)},
        .row_upper = &[_]f64{5},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &[_]usize{ 0, 1 }, &rows, &[_]f64{1}),
        .objective_sense = .maximize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.iteration_counters.attempted_iterations = 999;
    const status = engine.solveProblem(problem, .{});
    try std.testing.expectEqual(SolveStatus.optimal, status);
    try std.testing.expect(engine.iteration_counters.attempted_iterations < 999);
    try std.testing.expectEqual(
        engine.iteration_counters.committed_pivots,
        engine.iteration_counters.classifiedPivots(),
    );
    try std.testing.expectApproxEqAbs(@as(f64, 5), engine.basis.?.primal[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 5), engine.objective_value, 1e-12);
}

test "engine flips a nonbasic variable to its finite upper bound" {
    const rows = [_]foundation.RowId{foundation.RowId.fromUsizeAssumeValid(0)};
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{1},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{3},
        .row_lower = &[_]f64{-std.math.inf(f64)},
        .row_upper = &[_]f64{10},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &[_]usize{ 0, 1 }, &rows, &[_]f64{1}),
        .objective_sense = .maximize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{}));
    try std.testing.expectApproxEqAbs(@as(f64, 3), engine.basis.?.primal[0], 1e-12);
    try std.testing.expectEqual(basis_module.BasisStatus.at_upper, engine.basis.?.col_status[0]);
}

test "engine pivots a variable downward from its upper bound" {
    const rows = [_]foundation.RowId{foundation.RowId.fromUsizeAssumeValid(0)};
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{1},
        .col_lower = &[_]f64{-std.math.inf(f64)},
        .col_upper = &[_]f64{5},
        .row_lower = &[_]f64{2},
        .row_upper = &[_]f64{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &[_]usize{ 0, 1 }, &rows, &[_]f64{1}),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{}));
    try std.testing.expectApproxEqAbs(@as(f64, 2), engine.basis.?.primal[0], 1e-12);
    try std.testing.expectEqual(basis_module.BasisStatus.basic, engine.basis.?.col_status[0]);
}

test "engine normalizes and solves a ranged row" {
    const rows = [_]foundation.RowId{foundation.RowId.fromUsizeAssumeValid(0)};
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{1},
        .col_lower = &[_]f64{1},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{1},
        .row_upper = &[_]f64{5},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &[_]usize{ 0, 1 }, &rows, &[_]f64{1}),
        .objective_sense = .maximize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{}));
    try std.testing.expectApproxEqAbs(@as(f64, 5), engine.basis.?.primal[0], 1e-12);
}

test "engine rejects malformed borrowed problem dimensions" {
    const problem = problem_module.ProblemView{
        .num_rows = 0,
        .num_cols = 1,
        .col_cost = &.{},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{1},
        .row_lower = &.{},
        .row_upper = &.{},
        .matrix = matrix.CscView.initAssumeValid(0, 1, &[_]usize{ 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.numerical_failure, engine.solveProblem(problem, .{}));
}
