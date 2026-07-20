//! Revised-simplex orchestration boundary; policies remain replaceable.
const std = @import("std");
const basis_module = @import("basis.zig");
const basis_snapshot_module = @import("basis_snapshot.zig");
const factorization_module = @import("factorization.zig");
const pricing_module = @import("pricing.zig");
const ratio_module = @import("ratio_test.zig");
const numerical_module = @import("numerical.zig");
const dual_phase_one_module = @import("dual_phase_one.zig");
const crash_module = @import("crash.zig");
const degeneracy_module = @import("degeneracy.zig");
const pricing_workspace_module = @import("pricing_workspace.zig");
const problem_module = @import("problem.zig");
const solution_module = @import("solution.zig");
const foundation = @import("foundation");
const matrix = @import("matrix");

pub const Algorithm = enum { primal_revised, dual_revised };
pub const PhaseOneStrategy = enum { primal, dual, automatic };
pub const CrashStrategy = enum { logical, ltssf, bixby, automatic };
pub const DegeneracyStrategy = enum { baseline, perturbation, perturbation_taboo, automatic };
pub const PhaseOnePricingStrategy = enum { inherit, dantzig, devex, steepest_edge };
pub const PricingKernel = enum { column, row, automatic };
pub const DevexStrategy = enum { legacy, framework };
pub const PrimalPricingStrategy = enum { inherit, partial };
/// Dual edge-weight lifecycle. `inherit` preserves the caller's pricing rule;
/// `steepest_devex` starts with exact DSE and deterministically falls back to
/// full dual Devex when the recurrence is rejected or exceeds its work budget.
pub const DualEdgeWeightStrategy = enum { inherit, steepest_devex };
pub const SolvePhase = enum { phase_one, dual_phase_one, dual_feasibility_repair, phase_two };
pub const CallbackAction = enum { continue_solve, stop };
/// Borrowed scalar progress snapshot. It owns no memory and is valid only for
/// the callback invocation.
pub const ProgressEventView = struct {
    phase: SolvePhase,
    algorithm: Algorithm,
    iterations: usize,
    work_used: u64,
    objective_value: f64,
    primal_infeasibility: f64,
    dual_infeasibility: f64,
};
pub const IterationCallback = *const fn (event: ProgressEventView, user_data: ?*anyopaque) CallbackAction;
pub const IterationLogCallback = *const fn (event: ProgressEventView, user_data: ?*anyopaque) void;
pub const LogLevel = enum { off, iterations };
pub const SolveStatus = solution_module.SolveStatus;
pub const FailureSite = enum {
    none,
    direction_column,
    direction_solve,
    direction_refactor,
    direction_refinement,
    pivot_update,
    pivot_factorization,
    pivot_edge_weights,
    reduced_cost,
    optimality_check,
};
pub const PivotTraceEvent = struct {
    phase: SolvePhase,
    iteration: usize,
    entering_column: u32,
    leaving_column: u32,
    leaving_row: u32,
    pivot: f64,
    step: f64,
    update_count: usize,
    ftran_relative_residual: f64,
    condition_estimate: f64,
};
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
    phase: SolvePhase,
    iteration: usize,
    reason: DegeneracyReason,
    entering_column: u32,
    leaving_column: u32,
    step: f64,
    objective_change: f64,
    basis_fingerprint: u64,
};
pub const DualPhaseOneCandidateReason = enum {
    small_tableau,
    basic_or_fixed,
    wrong_pivot_sign,
    accepted_bound_flip,
    eligible_unselected,
};
pub const DualPhaseOneCandidateTraceEvent = struct {
    column: u32,
    status: basis_module.BasisStatus,
    tableau: f64,
    direction: f64,
    signed_pivot: f64,
    reduced_cost: f64,
    lower: f64,
    upper: f64,
    primal: f64,
    flip_capacity: f64,
    reason: DualPhaseOneCandidateReason,
};
pub const DualPhaseOneEpTraceEvent = struct { row: u32, value: f64 };
pub const DualPhaseOneFailureDiagnostic = struct {
    iteration: usize,
    leaving_row: u32,
    leaving_bound: basis_module.BasisStatus,
    violation: f64,
    ep_nonzeros: usize,
    ep_max_abs: f64,
    small_tableau: usize,
    basic_or_fixed: usize,
    wrong_pivot_sign: usize,
    accepted_bound_flips: usize,
    eligible_unselected: usize,
};
const PrimalLeavingResult = struct {
    status: SolveStatus,
    row: ?u32 = null,
    step: f64 = 0.0,
    bound: basis_module.BasisStatus = .at_lower,
};
pub const BasisImportError = basis_snapshot_module.BasisViewError || error{
    InvalidNonbasicStatus,
    SingularBasis,
    NumericalFailure,
};
pub const SolveControl = struct {
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
    iteration_callback: ?IterationCallback = null,
    callback_user_data: ?*anyopaque = null,
    callback_interval_work: u64 = 1,
    log_level: LogLevel = .off,
    log_callback: ?IterationLogCallback = null,
    log_user_data: ?*anyopaque = null,
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
    phase_one_pricing: PhaseOnePricingStrategy = .inherit,
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
    /// Explicit dual DSE -> Devex A/B policy. It is not enabled by default.
    dual_edge_weight_strategy: DualEdgeWeightStrategy = .inherit,
    /// Maximum successful DSE recurrence updates per dual phase before
    /// switching to dual Devex. Zero disables the budget-triggered switch.
    dual_dse_update_budget: usize = 64,
};

pub const SimplexStats = struct {
    phase_one_iterations: usize = 0,
    dual_phase_one_iterations: usize = 0,
    dual_phase_one_fallbacks: usize = 0,
    crash_attempts: usize = 0,
    crash_fallbacks: usize = 0,
    crash_planned_columns: usize = 0,
    crash_structural_columns: usize = 0,
    crash_basis_nonzeros: usize = 0,
    crash_condition_estimate: f64 = 1.0,
    dual_repair_iterations: usize = 0,
    phase_two_iterations: usize = 0,
    phase_one_ns: u64 = 0,
    dual_phase_one_ns: u64 = 0,
    dual_repair_ns: u64 = 0,
    phase_two_ns: u64 = 0,
    cleanup_ns: u64 = 0,
    rebuild_calls: usize = 0,
    rebuild_ns: u64 = 0,
    pricing_calls: usize = 0,
    pricing_ns: u64 = 0,
    pricing_samples: usize = 0,
    pricing_nonzeros: usize = 0,
    pricing_entries: usize = 0,
    bound_flips: usize = 0,
    bound_flip_batches: usize = 0,
    bound_flip_ftran_savings: usize = 0,
    dual_reduced_cost_updates: usize = 0,
    dual_exact_reprices: usize = 0,
    dual_dse_updates: usize = 0,
    dual_devex_updates: usize = 0,
    dual_dse_invalid_fallbacks: usize = 0,
    dual_dse_budget_fallbacks: usize = 0,
    devex_frameworks: usize = 0,
    devex_framework_updates: usize = 0,
    devex_bad_weights: usize = 0,
    degeneracy_bound_ties: usize = 0,
    degeneracy_ratio_ties: usize = 0,
    degeneracy_zero_primal_steps: usize = 0,
    degeneracy_phase_one_objective_stalls: usize = 0,
    degeneracy_repeated_bases: usize = 0,
    degeneracy_small_pivot_retries: usize = 0,
    degeneracy_bound_flips: usize = 0,
    perturbation_activations: usize = 0,
    perturbation_expirations: usize = 0,
    perturbation_cleanups: usize = 0,
    cold_restart_solves: usize = 0,
    cold_restart_phase_one: usize = 0,
    taboo_records: usize = 0,
    taboo_retries: usize = 0,
    aq_samples: usize = 0,
    aq_nonzeros: usize = 0,
    ep_samples: usize = 0,
    ep_nonzeros: usize = 0,
    dense_pricing_dispatches: usize = 0,
    hyper_pricing_dispatches: usize = 0,
    row_pricing_dispatches: usize = 0,
    column_pricing_dispatches: usize = 0,
    rebuild_phase_one_setup: usize = 0,
    rebuild_solve_residual: usize = 0,
    rebuild_small_pivot: usize = 0,
    rebuild_update_limit: usize = 0,
    rebuild_update_growth: usize = 0,
    rebuild_direction_refinement: usize = 0,
    rebuild_fresh_mode: usize = 0,
    rebuild_update_rejected: usize = 0,
    rebuild_cleanup: usize = 0,
    rebuild_edge_weight_reset: usize = 0,
    rebuild_numerical_policy: usize = 0,
    first_update_failure_kind: ?factorization_module.UpdateFailureKind = null,
    first_update_failure_iteration: usize = 0,
    first_update_failure_entering: u32 = 0,
    first_update_failure_leaving_row: u32 = 0,

    pub fn classifiedRebuilds(self: SimplexStats) usize {
        return self.rebuild_phase_one_setup + self.rebuild_solve_residual + self.rebuild_small_pivot +
            self.rebuild_update_limit + self.rebuild_update_growth + self.rebuild_direction_refinement +
            self.rebuild_fresh_mode + self.rebuild_update_rejected + self.rebuild_cleanup +
            self.rebuild_edge_weight_reset + self.rebuild_numerical_policy;
    }

    pub fn classifiedDegeneratePivots(self: SimplexStats) usize {
        return self.degeneracy_bound_ties + self.degeneracy_ratio_ties +
            self.degeneracy_zero_primal_steps + self.degeneracy_phase_one_objective_stalls +
            self.degeneracy_repeated_bases + self.degeneracy_small_pivot_retries +
            self.degeneracy_bound_flips;
    }
};

const RebuildReason = enum {
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

pub const SimplexEngine = struct {
    allocator: std.mem.Allocator,
    algorithm: Algorithm = .dual_revised,
    basis: ?basis_module.BasisState = null,
    factorization: factorization_module.Factorization,
    dual_phase_one: dual_phase_one_module.DualPhaseOneWorkspace,
    crash: crash_module.CrashWorkspace,
    degeneracy: degeneracy_module.Workspace,
    pricing_row_view: pricing_workspace_module.RowView,
    pricing: pricing_module.Pricing = .{},
    ratio_test: ratio_module.RatioTest = .{},
    numerical: numerical_module.NumericalState = .{},
    iterations: usize = 0,
    objective_value: f64 = 0.0,
    objective_scale: f64 = 1.0,
    unbounded_ray_valid: bool = false,
    phase1_needed: bool = false,
    solve_start_ns: ?i96 = null,
    solve_clock_io: ?std.Io = null,
    work_used: u64 = 0,
    dual_edge_weights_valid: bool = false,
    dual_row_index: ?u32 = null,
    dual_candidate_count: usize = 0,
    dual_candidate_cutoff: f64 = 0.0,
    dual_hyper_sparse_active: bool = false,
    /// Structural basic columns replaced by logical columns while repairing a
    /// singular imported basis during the current solve.
    rank_repair_count: usize = 0,
    pivot_trace_count: usize = 0,
    active_pivot_trace: []PivotTraceEvent = &.{},
    degeneracy_trace_count: usize = 0,
    active_degeneracy_trace: []DegeneracyTraceEvent = &.{},
    degeneracy_basis_fingerprints: [32]u64 = @splat(0),
    degeneracy_basis_fingerprint_count: usize = 0,
    degeneracy_basis_fingerprint_cursor: usize = 0,
    dual_phase_one_failure: ?DualPhaseOneFailureDiagnostic = null,
    dual_phase_one_candidate_trace_count: usize = 0,
    active_dual_phase_one_candidate_trace: []DualPhaseOneCandidateTraceEvent = &.{},
    dual_phase_one_ep_trace_count: usize = 0,
    active_dual_phase_one_ep_trace: []DualPhaseOneEpTraceEvent = &.{},
    current_phase: SolvePhase = .phase_two,
    direction_requires_reinversion: bool = false,
    fresh_factorization_pivots_remaining: usize = 0,
    reduced_cost_update_count: usize = 0,
    reduced_cost_refresh_period: usize = 8,
    failure_site: FailureSite = .none,
    stats: SimplexStats = .{},
    statistics_io: ?std.Io = null,
    cleanup_active: bool = false,
    active_degeneracy_strategy: DegeneracyStrategy = .baseline,
    active_adaptive_reprice: bool = false,
    maximum_reduced_cost_drift: f64 = 0.0,
    exact_reprices: usize = 0,
    active_pricing_kernel: PricingKernel = .column,
    active_devex_strategy: DevexStrategy = .legacy,
    active_primal_pricing_strategy: PrimalPricingStrategy = .inherit,
    active_dual_edge_weight_strategy: DualEdgeWeightStrategy = .inherit,
    active_dual_dse_update_budget: usize = 64,
    dual_dse_updates_since_start: usize = 0,
    devex_framework_iterations: usize = 0,
    devex_bad_weight_count: usize = 0,
    devex_reset_after_pivot: bool = false,

    pub fn init(a: std.mem.Allocator) SimplexEngine {
        return .{
            .allocator = a,
            .factorization = factorization_module.Factorization.init(a),
            .dual_phase_one = dual_phase_one_module.DualPhaseOneWorkspace.init(a),
            .crash = crash_module.CrashWorkspace.init(a),
            .degeneracy = degeneracy_module.Workspace.init(a),
            .pricing_row_view = pricing_workspace_module.RowView.init(a),
        };
    }
    pub fn deinit(self: *SimplexEngine) void {
        if (self.basis) |*b| b.deinit();
        self.factorization.deinit();
        self.dual_phase_one.deinit();
        self.crash.deinit();
        self.degeneracy.deinit();
        self.pricing_row_view.deinit();
    }

    pub fn requestedBytes(self: *const SimplexEngine) usize {
        const basis_bytes = if (self.basis) |*basis| basis.requestedBytes() else 0;
        return basis_bytes + self.factorization.requestedBytes() + self.dual_phase_one.requestedBytes() +
            self.crash.requestedBytes() + self.degeneracy.requestedBytes() + self.pricing_row_view.requestedBytes();
    }
    pub fn solve(_: *SimplexEngine, _: usize, _: usize, _: SolveControl) SolveStatus {
        return .not_implemented;
    }

    /// Solve entry point consuming a borrowed LP `ProblemView`.
    ///
    /// The view must outlive this call; the engine never takes ownership of
    /// model arrays. Basis storage is owned by the engine instance.
    pub fn solveProblem(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
        self.pricing.resetPartial();
        self.startSolveClock(control);
        self.resetStatistics(control);
        self.work_used = 0;
        self.rank_repair_count = 0;
        self.fresh_factorization_pivots_remaining = 0;
        self.reduced_cost_update_count = 0;
        self.failure_site = .none;
        self.unbounded_ray_valid = false;
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
        if (crash_feasibility.primal) return self.solvePrimal(problem, control);
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
            const dual_phase_one_status = self.solveDualPhaseOne(problem, control);
            if (dual_phase_one_status != .not_implemented and dual_phase_one_status != .numerical_failure)
                return dual_phase_one_status;
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

    /// Materialize the logical crash solution from current nonbasic structural
    /// values. The basis remains the identity; values may violate logical
    /// bounds so the caller can choose primal, dual, or artificial Phase I.
    fn initializeLogicalBasicValues(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        @memcpy(basis.basic_value, basis.row_rhs);
        for (0..problem.num_cols) |column| {
            const initial = basis.primal[column];
            if (initial == 0.0) continue;
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
                if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                const row_index = row.toUsize();
                basis.basic_value[row_index] -= basis.row_scale[row_index] * coefficient * basis.column_scale[column] * initial;
            }
        }
        for (basis.basic_index, basis.basic_value, 0..) |column, value, row| {
            if (!std.math.isFinite(value)) return .numerical_failure;
            basis.primal[column] = value;
            basis.basic_lower[row] = basis.col_lower[column];
            basis.basic_upper[row] = basis.col_upper[column];
        }
        return .optimal;
    }

    /// Replace each infeasible logical basic with its nonnegative artificial
    /// column. Initial Phase I and every transactional restart must use this
    /// exact state transition; rebuilding only the logical identity would
    /// leave Phase-I costs disconnected from the violated rows.
    fn installArtificialPhaseOneBasis(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        self.phase1_needed = false;
        for (basis.basic_value, basis.basic_lower, basis.basic_upper, 0..) |*value, lower, upper, row| {
            const logical_col = problem.num_cols + row;
            const artificial_col = problem.num_cols + problem.num_rows + row;
            if (!std.math.isFinite(value.*)) return .numerical_failure;
            if (value.* < lower - self.numerical.primal_tolerance) {
                basis.primal[logical_col] = lower;
                basis.col_status[logical_col] = .at_lower;
                basis.basic_pos[logical_col] = std.math.maxInt(u32);
                basis.artificial_sign[row] = -1.0;
                value.* = lower - value.*;
                basis.basic_index[row] = @intCast(artificial_col);
                basis.basic_pos[artificial_col] = @intCast(row);
                basis.col_status[artificial_col] = .basic;
                basis.col_upper[artificial_col] = std.math.inf(f64);
                basis.basic_lower[row] = 0.0;
                basis.basic_upper[row] = std.math.inf(f64);
                self.phase1_needed = true;
            } else if (value.* > upper + self.numerical.primal_tolerance) {
                basis.primal[logical_col] = upper;
                basis.col_status[logical_col] = .at_upper;
                basis.basic_pos[logical_col] = std.math.maxInt(u32);
                basis.artificial_sign[row] = 1.0;
                value.* -= upper;
                basis.basic_index[row] = @intCast(artificial_col);
                basis.basic_pos[artificial_col] = @intCast(row);
                basis.col_status[artificial_col] = .basic;
                basis.col_upper[artificial_col] = std.math.inf(f64);
                basis.basic_lower[row] = 0.0;
                basis.basic_upper[row] = std.math.inf(f64);
                self.phase1_needed = true;
            } else {
                basis.primal[logical_col] = value.*;
            }
        }
        return .optimal;
    }

    /// Install a sparse structural crash basis transactionally. The caller
    /// restores the logical identity on false, so no partially matched basis
    /// can escape numerical validation.
    fn installSparseCrashBasis(
        self: *SimplexEngine,
        problem: problem_module.ProblemView,
        maximum_columns: ?usize,
        scoring: crash_module.CrashScoring,
    ) bool {
        const basis = if (self.basis) |*value| value else return false;
        const average_column_degree = if (problem.num_cols == 0) 0 else (problem.matrix.values.len + problem.num_cols - 1) / problem.num_cols;
        const near_singleton_limit: u32 = @intCast(@min(@max(average_column_degree * 2, 4), 32));
        const plan_view = self.crash.plan(
            problem.matrix,
            problem.col_cost,
            basis.col_lower[0..problem.num_cols],
            basis.col_upper[0..problem.num_cols],
            basis.row_scale,
            basis.column_scale[0..problem.num_cols],
            near_singleton_limit,
            self.numerical.pivot_tolerance,
            scoring,
        ) catch return false;
        if (plan_view.rows.len == 0) return false;
        self.stats.crash_planned_columns = plan_view.rows.len;

        // A near-singleton symbolic matching can still be numerically rank
        // deficient. Validate the longest deterministic prefix first, then
        // halve it until a safe basis is found. No partial basis is published.
        // Retain at least one quarter of the logical columns as stable anchors.
        // A fully structural crash can be perfectly invertible yet leave the
        // current dual Phase-I update path without a safe reinversion route.
        const anchored_limit = if (problem.num_rows == 0) 0 else @max(@as(usize, 1), problem.num_rows * 3 / 4);
        var prefix_count = @min(plan_view.rows.len, maximum_columns orelse anchored_limit);
        while (prefix_count != 0) : (prefix_count /= 2) {
            if (self.initializeProblemStorage(problem) != .optimal) return false;
            var basis_nonzeros = problem.num_rows - prefix_count;
            for (plan_view.rows[0..prefix_count], plan_view.columns[0..prefix_count]) |row_u32, column_u32| {
                const row: usize = @intCast(row_u32);
                const column: usize = @intCast(column_u32);
                const logical = problem.num_cols + row;
                const leaving_status = nonbasicStatusForBounds(basis.col_lower[logical], basis.col_upper[logical]);
                basis.primal[logical] = switch (leaving_status) {
                    .at_lower, .fixed => basis.col_lower[logical],
                    .at_upper => basis.col_upper[logical],
                    .free, .superbasic => 0.0,
                    .basic => unreachable,
                };
                basis.applyPivot(row, column, leaving_status) catch return false;
                basis.basic_lower[row] = basis.col_lower[column];
                basis.basic_upper[row] = basis.col_upper[column];
                const begin = problem.matrix.col_starts[column];
                const end = problem.matrix.col_starts[column + 1];
                for (problem.matrix.values[begin..end]) |coefficient| {
                    if (matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) basis_nonzeros += 1;
                }
            }
            self.factorizeCurrentBasis(problem) catch continue;
            if (self.finishRefactorization() != .optimal) continue;
            const condition = self.factorization.pivotConditionEstimate();
            if (!std.math.isFinite(condition) or condition > 1e12) continue;
            if (self.recomputeBasicValuesUnchecked(problem) != .optimal) continue;
            if (!std.math.isFinite(self.numerical.last_relative_residual) or
                self.numerical.last_relative_residual > self.numerical.residual_tolerance * 100.0)
                continue;
            self.stats.crash_structural_columns = prefix_count;
            self.stats.crash_basis_nonzeros = basis_nonzeros;
            self.stats.crash_condition_estimate = condition;
            return true;
        }
        return false;
    }

    /// Reoptimize with the existing basis and factorization. The caller must
    /// guarantee unchanged matrix structure and values. No model view is
    /// retained after this call.
    pub fn reoptimizeProblem(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
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

    fn refreshProblemStorage(self: *SimplexEngine, problem: problem_module.ProblemView) bool {
        const basis = if (self.basis) |*value| value else return false;
        if (basis.num_rows != problem.num_rows or basis.num_structural_cols != problem.num_cols) return false;
        for (problem.col_lower, problem.col_upper, 0..) |lower, upper, column| {
            if (lower > upper) return false;
            const scale = basis.column_scale[column];
            basis.col_lower[column] = lower / scale;
            basis.col_upper[column] = upper / scale;
            switch (basis.col_status[column]) {
                .basic => {},
                .at_lower => {
                    if (!std.math.isFinite(lower)) return false;
                    basis.primal[column] = lower / scale;
                },
                .at_upper => {
                    if (!std.math.isFinite(upper)) return false;
                    basis.primal[column] = upper / scale;
                },
                .fixed => {
                    if (lower != upper) return false;
                    basis.primal[column] = lower / scale;
                },
                .free, .superbasic => basis.primal[column] = @min(@max(0.0, lower), upper) / scale,
            }
        }
        for (problem.row_lower, problem.row_upper, 0..) |lower, upper, row| {
            const sign: f64 = if (std.math.isFinite(upper)) 1.0 else if (std.math.isFinite(lower)) -1.0 else 0.0;
            if (std.math.sign(basis.row_scale[row]) != sign) return false;
            const magnitude = @abs(basis.row_scale[row]);
            basis.row_rhs[row] = if (sign > 0.0) upper * magnitude else if (sign < 0.0) -lower * magnitude else 0.0;
            const logical = problem.num_cols + row;
            basis.col_lower[logical] = 0.0;
            basis.col_upper[logical] = if (sign > 0.0 and std.math.isFinite(lower)) (upper - lower) * magnitude else std.math.inf(f64);
            switch (basis.col_status[logical]) {
                .basic => {},
                .at_lower, .fixed => basis.primal[logical] = basis.col_lower[logical],
                .at_upper => {
                    if (!std.math.isFinite(basis.col_upper[logical])) return false;
                    basis.primal[logical] = basis.col_upper[logical];
                },
                .free, .superbasic => basis.primal[logical] = 0.0,
            }
        }
        for (basis.basic_index, 0..) |column, row| {
            basis.basic_lower[row] = basis.col_lower[column];
            basis.basic_upper[row] = basis.col_upper[column];
        }
        self.objective_scale = objectiveScale(problem.col_cost);
        return true;
    }

    fn initializeProblemStorage(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        basis.initializeSlackBasis();
        @memset(basis.primal, 0.0);
        @memset(basis.artificial_sign, 0.0);
        for (problem.col_lower, problem.col_upper, 0..) |lower, upper, col| {
            if (lower > upper) return .infeasible;
            if (lower == upper) {
                basis.primal[col] = lower;
                basis.col_status[col] = .fixed;
            } else if (std.math.isFinite(lower)) {
                basis.primal[col] = lower;
                basis.col_status[col] = .at_lower;
            } else if (std.math.isFinite(upper)) {
                basis.primal[col] = upper;
                basis.col_status[col] = .at_upper;
            } else {
                basis.primal[col] = 0.0;
                basis.col_status[col] = .free;
            }
            basis.col_lower[col] = lower;
            basis.col_upper[col] = upper;
        }
        @memset(basis.residual_work, 0.0);
        for (problem.matrix.row_indices, problem.matrix.values) |row, value| {
            if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(value)) continue;
            const row_index = row.toUsize();
            basis.residual_work[row_index] = @max(basis.residual_work[row_index], @abs(value));
        }
        for (problem.row_lower, problem.row_upper, 0..) |lower, upper, row| {
            if (lower > upper) return .infeasible;
            const maximum = basis.residual_work[row];
            const magnitude = if (maximum == 0.0)
                1.0
            else
                @exp2(std.math.clamp(@round(-@log2(maximum)), -20.0, 20.0));
            if (std.math.isFinite(upper)) {
                basis.row_scale[row] = magnitude;
                basis.row_rhs[row] = upper * magnitude;
                basis.col_lower[problem.num_cols + row] = 0.0;
                basis.col_upper[problem.num_cols + row] = if (std.math.isFinite(lower)) (upper - lower) * magnitude else std.math.inf(f64);
            } else if (std.math.isFinite(lower)) {
                basis.row_scale[row] = -magnitude;
                basis.row_rhs[row] = -lower * magnitude;
                basis.col_lower[problem.num_cols + row] = 0.0;
                basis.col_upper[problem.num_cols + row] = std.math.inf(f64);
            } else {
                basis.row_scale[row] = 0.0;
                basis.row_rhs[row] = 0.0;
                basis.col_lower[problem.num_cols + row] = 0.0;
                basis.col_upper[problem.num_cols + row] = std.math.inf(f64);
            }
            basis.basic_lower[row] = basis.col_lower[problem.num_cols + row];
            basis.basic_upper[row] = basis.col_upper[problem.num_cols + row];
        }
        var scaled_matrix_min = std.math.inf(f64);
        var scaled_matrix_max: f64 = 0.0;
        for (problem.matrix.row_indices, problem.matrix.values) |row, coefficient| {
            if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
            const magnitude = @abs(basis.row_scale[row.toUsize()] * coefficient);
            if (magnitude != 0.0) scaled_matrix_min = @min(scaled_matrix_min, magnitude);
            scaled_matrix_max = @max(scaled_matrix_max, magnitude);
        }
        const use_column_scaling = scaled_matrix_min < std.math.inf(f64) and scaled_matrix_max / scaled_matrix_min > 1e6;
        for (0..problem.num_cols) |column| {
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            var maximum: f64 = 0.0;
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
                if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                maximum = @max(maximum, @abs(basis.row_scale[row.toUsize()] * coefficient));
            }
            const scale = if (use_column_scaling) powerOfTwoScale(maximum, 20) else 1.0;
            basis.column_scale[column] = scale;
            basis.col_lower[column] = problem.col_lower[column] / scale;
            basis.col_upper[column] = problem.col_upper[column] / scale;
            basis.primal[column] /= scale;
        }
        self.objective_scale = objectiveScale(problem.col_cost);
        return .optimal;
    }

    fn powerOfTwoScale(maximum: f64, exponent_limit: comptime_int) f64 {
        if (maximum == 0.0 or !std.math.isFinite(maximum)) return 1.0;
        return @exp2(std.math.clamp(@round(-@log2(maximum)), -@as(f64, exponent_limit), @as(f64, exponent_limit)));
    }

    fn objectiveScale(cost: []const f64) f64 {
        var maximum: f64 = 0.0;
        for (cost) |value| maximum = @max(maximum, @abs(value));
        return powerOfTwoScale(maximum, 15);
    }

    /// Restore a validated borrowed basis and rebuild factorization/basic
    /// values without retaining any caller memory.
    pub fn importBasis(self: *SimplexEngine, problem: problem_module.ProblemView, view: basis_snapshot_module.BasisView) BasisImportError!void {
        try view.validate(problem.num_cols, problem.num_rows);
        const basis = if (self.basis) |*value| value else return error.NumericalFailure;
        @memcpy(basis.col_status[0..problem.num_cols], view.structural_status);
        @memcpy(basis.col_status[problem.num_cols..][0..problem.num_rows], view.logical_status);
        @memset(basis.col_status[problem.num_cols + problem.num_rows ..], .fixed);

        for (basis.col_status[0 .. problem.num_cols + problem.num_rows], basis.col_lower[0 .. problem.num_cols + problem.num_rows], basis.col_upper[0 .. problem.num_cols + problem.num_rows]) |status, lower, upper| {
            const valid = switch (status) {
                .basic => true,
                .at_lower => std.math.isFinite(lower),
                .at_upper => std.math.isFinite(upper),
                .fixed => std.math.isFinite(lower) and lower == upper,
                .free => !std.math.isFinite(lower) and !std.math.isFinite(upper),
                .superbasic => 0.0 >= lower and 0.0 <= upper,
            };
            if (!valid) return error.InvalidNonbasicStatus;
        }
        @memset(basis.basic_pos, std.math.maxInt(u32));
        @memcpy(basis.basic_index, view.basic_index);
        for (basis.basic_index, 0..) |column, row| basis.basic_pos[column] = @intCast(row);

        for (basis.col_status[0 .. problem.num_cols + problem.num_rows], 0..) |status, column| {
            basis.primal[column] = switch (status) {
                .at_lower, .fixed => basis.col_lower[column],
                .at_upper => basis.col_upper[column],
                .free, .superbasic, .basic => 0.0,
            };
        }
        for (basis.basic_index, 0..) |column, row| {
            basis.basic_lower[row] = basis.col_lower[column];
            basis.basic_upper[row] = basis.col_upper[column];
        }
        self.factorizeCurrentBasis(problem) catch |err| switch (err) {
            error.Singular => try self.repairRankDeficientBasis(problem),
            error.OutOfMemory => return error.OutOfMemory,
            error.DimensionMismatch, error.NotImplemented, error.NumericalFailure => return error.NumericalFailure,
        };
        if (self.finishRefactorization() != .optimal) return error.NumericalFailure;
        if (self.recomputeBasicValues(problem) != .optimal) {
            // A valid dual warm start is allowed to be primal infeasible, so
            // rebuild values without enforcing bounds.
            if (self.recomputeBasicValuesUnchecked(problem) != .optimal) return error.NumericalFailure;
        }
    }

    pub fn exportBasisView(self: *const SimplexEngine, problem: problem_module.ProblemView) ?basis_snapshot_module.BasisView {
        const basis = if (self.basis) |*value| value else return null;
        const artificial_begin = problem.num_cols + problem.num_rows;
        for (basis.basic_index) |column| if (column >= artificial_begin) return null;
        return .{
            .structural_status = basis.col_status[0..problem.num_cols],
            .logical_status = basis.col_status[problem.num_cols..][0..problem.num_rows],
            .basic_index = basis.basic_index,
        };
    }

    pub fn exportBasisSnapshot(self: *const SimplexEngine, allocator: std.mem.Allocator, problem: problem_module.ProblemView) BasisImportError!basis_snapshot_module.BasisSnapshot {
        const basis_view = self.exportBasisView(problem) orelse return error.NumericalFailure;
        return basis_snapshot_module.BasisSnapshot.initFromView(allocator, basis_view);
    }

    pub fn classifyFeasibility(self: *const SimplexEngine, problem: problem_module.ProblemView) struct { primal: bool, dual: bool } {
        const basis = if (self.basis) |*value| value else return .{ .primal = false, .dual = false };
        var primal = true;
        for (basis.basic_value, basis.basic_lower, basis.basic_upper) |value, lower, upper| {
            if (value < lower - self.numerical.primal_tolerance or value > upper + self.numerical.primal_tolerance) {
                primal = false;
                break;
            }
        }
        var dual = true;
        for (basis.reduced_cost[0 .. problem.num_cols + problem.num_rows], basis.col_status[0 .. problem.num_cols + problem.num_rows]) |reduced, status| {
            const infeasible = switch (status) {
                .at_lower => reduced < -self.numerical.dual_tolerance,
                .at_upper => reduced > self.numerical.dual_tolerance,
                .free, .superbasic => @abs(reduced) > self.numerical.dual_tolerance,
                .basic, .fixed => false,
            };
            if (infeasible) {
                dual = false;
                break;
            }
        }
        return .{ .primal = primal, .dual = dual };
    }

    /// Compute the revised-simplex entering direction `B^-1 A_j` into the
    /// engine-owned SoA workspace without allocating or copying model data.
    pub fn computeDirection(self: *SimplexEngine, problem: problem_module.ProblemView, entering_col: usize) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        self.direction_requires_reinversion = false;
        if (self.fillInternalColumn(problem, entering_col, basis.pivot_direction) != .optimal) {
            self.failure_site = .direction_column;
            return .numerical_failure;
        }
        self.factorization.solveForUpdate(basis.pivot_direction) catch {
            self.failure_site = .direction_solve;
            return .numerical_failure;
        };
        if (!self.directionResidualAcceptable(problem, entering_col)) {
            self.factorization.recordReinversion(.solve_residual);
            if (self.refactorizeBasis(problem, .solve_residual) != .optimal) {
                self.failure_site = .direction_refactor;
                return .numerical_failure;
            }
            if (self.fillInternalColumn(problem, entering_col, basis.pivot_direction) != .optimal) return .numerical_failure;
            self.factorization.solveForUpdate(basis.pivot_direction) catch return .numerical_failure;
            var refinement_step: usize = 0;
            while (!self.directionResidualAcceptable(problem, entering_col)) : (refinement_step += 1) {
                if (refinement_step >= self.numerical.max_refinement_steps) {
                    self.failure_site = .direction_refinement;
                    return .numerical_failure;
                }
                self.factorization.solve(basis.residual_work) catch {
                    self.failure_site = .direction_refinement;
                    return .numerical_failure;
                };
                for (basis.pivot_direction, basis.residual_work) |*value, correction| value.* += correction;
                self.numerical.refinement_count += 1;
                // The captured partial aq predates refinement. Force a fresh
                // base factorization after this pivot instead of publishing an
                // FT update from inconsistent capture state.
                self.direction_requires_reinversion = true;
            }
        }
        self.observeAqDensity(basis.pivot_direction);
        return .optimal;
    }

    /// Validate an updated FTRAN before its direction reaches ratio testing.
    /// A drifting FT chain can otherwise select an invalid leaving row and
    /// permanently corrupt the basis before the periodic growth gate fires.
    fn directionResidualAcceptable(self: *SimplexEngine, problem: problem_module.ProblemView, entering_col: usize) bool {
        const basis = if (self.basis) |*value| value else return false;
        if (self.fillInternalColumn(problem, entering_col, basis.residual_work) != .optimal) return false;
        for (basis.rhs_work, basis.residual_work) |*magnitude, rhs| magnitude.* = @abs(rhs);
        self.subtractBasisProductWithMagnitude(
            problem,
            basis.pivot_direction,
            basis.residual_work,
            basis.rhs_work,
        ) catch return false;
        var residual_max: f64 = 0.0;
        var equation_scale: f64 = 0.0;
        for (basis.residual_work) |value| residual_max = @max(residual_max, @abs(value));
        for (basis.rhs_work) |value| equation_scale = @max(equation_scale, value);
        // Measure normwise backward error. Dividing only by |a_q| rejects a
        // backward-stable solve whenever a poorly scaled basis produces a
        // large x whose B*x terms cancel. |a_q| + |B|*|x| is the natural
        // scale of the equations and is also what iterative refinement can
        // meaningfully improve.
        const relative = residual_max / @max(1.0, equation_scale);
        self.numerical.last_ftran_relative_residual = relative;
        self.numerical.max_ftran_relative_residual = @max(self.numerical.max_ftran_relative_residual, relative);
        return std.math.isFinite(relative) and relative <= self.numerical.residual_tolerance;
    }

    /// Choose the leaving basic row for the currently materialized pivot
    /// direction. Returns `unbounded` when no positive pivot coefficient
    /// limits the entering variable.
    pub fn chooseLeaving(self: *SimplexEngine) PrimalLeavingResult {
        return self.chooseLeavingWithPolicy(true);
    }

    fn chooseLeavingWithPolicy(self: *SimplexEngine, allow_bland: bool) PrimalLeavingResult {
        const basis = if (self.basis) |*value| value else return .{ .status = .numerical_failure };
        for (basis.basic_margin, basis.ratio_direction, basis.basic_value, basis.basic_lower, basis.basic_upper, basis.pivot_direction) |*margin, *ratio_direction, value, lower, upper, direction| {
            if (direction > self.numerical.zero_tolerance) {
                margin.* = value - lower;
                ratio_direction.* = direction;
            } else if (direction < -self.numerical.zero_tolerance) {
                margin.* = upper - value;
                ratio_direction.* = -direction;
            } else {
                margin.* = 0.0;
                ratio_direction.* = 0.0;
            }
        }
        // Explicit baseline retains the proven Bland fallback (notably
        // feasible starts such as kb2). Automatic mode keeps Harris during its
        // 256-pivot hysteresis window, then uses bounded perturbation; this
        // avoids changing short local ties into weak Bland pivots.
        if (allow_bland and self.numerical.anti_cycling_active and
            self.active_degeneracy_strategy == .baseline)
            return self.chooseLeavingBland();
        const choice = if (self.degeneracy.active)
            self.ratio_test.chooseLeavingPerturbed(
                basis.ratio_direction,
                basis.basic_margin,
                self.degeneracy.row_rank[0..basis.basic_margin.len],
                self.numerical.primal_tolerance,
            )
        else
            self.ratio_test.chooseLeaving(basis.ratio_direction, basis.basic_margin);
        if (choice.row == null) return .{ .status = .unbounded };
        const row = choice.row.?;
        const bound: basis_module.BasisStatus = if (basis.pivot_direction[@intCast(row)] > 0) .at_lower else .at_upper;
        return .{ .status = .optimal, .row = row, .step = choice.step, .bound = bound };
    }

    fn chooseLeavingBland(self: *SimplexEngine) PrimalLeavingResult {
        const basis = if (self.basis) |*value| value else return .{ .status = .numerical_failure };
        var best_row: ?u32 = null;
        var best_column: u32 = std.math.maxInt(u32);
        var best_step = std.math.inf(f64);
        const tie_tolerance = self.numerical.perturbation;
        var maximum_direction: f64 = 0.0;
        for (basis.ratio_direction) |direction| maximum_direction = @max(maximum_direction, direction);
        const pivot_tolerance = @max(self.numerical.zero_tolerance, self.ratio_test.tolerance * maximum_direction);
        for (basis.ratio_direction, basis.basic_margin, basis.basic_index, 0..) |direction, margin, basic_column, row| {
            if (direction <= pivot_tolerance) continue;
            const step = @max(margin / direction, 0.0);
            if (!std.math.isFinite(step)) continue;
            if (step < best_step - tie_tolerance or
                (@abs(step - best_step) <= tie_tolerance and basic_column < best_column))
            {
                best_step = step;
                best_column = basic_column;
                best_row = @intCast(row);
            }
        }
        const row = best_row orelse return .{ .status = .unbounded };
        const bound: basis_module.BasisStatus = if (basis.pivot_direction[row] > 0.0) .at_lower else .at_upper;
        return .{ .status = .optimal, .row = row, .step = best_step, .bound = bound };
    }

    /// Forward-error guard for a pivotal entry extracted from an updated
    /// FTRAN. `sqrt(epsilon)` is the boundary where cancellation can leave
    /// fewer than half the significand bits. The unit floor also catches an
    /// entire entering column that is small before column equilibration.
    fn pivotNeedsFreshFactorization(self: *const SimplexEngine, leaving_row: u32) bool {
        const basis = if (self.basis) |*value| value else return true;
        const row: usize = @intCast(leaving_row);
        if (row >= basis.pivot_direction.len) return true;
        var maximum: f64 = 0.0;
        for (basis.pivot_direction) |value| maximum = @max(maximum, @abs(value));
        return @abs(basis.pivot_direction[row]) <= @sqrt(std.math.floatEps(f64)) * @max(1.0, maximum);
    }

    /// Apply one primal pivot and rebuild the dense basis factorization in
    /// existing storage. This is allocation-free; update factorizations can
    /// replace the reinversion later without changing the state transition.
    pub fn performPivot(self: *SimplexEngine, problem: problem_module.ProblemView, entering_col: usize, entering_direction: f64, leaving_row: usize, leaving_bound: basis_module.BasisStatus, step: f64) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        if (entering_col >= basis.col_status.len or leaving_row >= problem.num_rows) return .numerical_failure;
        const leaving_col = basis.basic_index[leaving_row];
        const pivot = basis.pivot_direction[leaving_row];
        if (self.pricing.rule == .steepest_edge and self.dual_row_index == @as(u32, @intCast(leaving_row))) {
            if (self.updateDualSteepestEdgeWeights(leaving_row, pivot) == .optimal) {
                self.stats.dual_dse_updates += 1;
                self.dual_dse_updates_since_start += 1;
                if (self.active_dual_edge_weight_strategy == .steepest_devex and
                    self.active_dual_dse_update_budget != 0 and
                    self.dual_dse_updates_since_start >= self.active_dual_dse_update_budget)
                    self.switchDualDseToDevex(.budget);
            } else if (self.active_dual_edge_weight_strategy == .steepest_devex) {
                self.switchDualDseToDevex(.invalid);
            } else {
                self.dual_edge_weights_valid = false;
            }
        } else {
            self.dual_edge_weights_valid = false;
            if (self.algorithm == .dual_revised and self.pricing.rule == .devex and
                self.active_dual_edge_weight_strategy == .steepest_devex)
                self.updateDualDevexWeights(leaving_row, pivot)
            else if (self.algorithm == .dual_revised or self.active_devex_strategy == .legacy)
                self.updateLegacyDevexWeights(entering_col, leaving_row, pivot);
        }
        self.dual_row_index = null;
        var forced_rebuild_reason: ?RebuildReason = null;
        const update_succeeded = blk: {
            if (self.direction_requires_reinversion) {
                forced_rebuild_reason = .direction_refinement;
                break :blk false;
            }
            if (self.fresh_factorization_pivots_remaining != 0) {
                forced_rebuild_reason = .fresh_mode;
                break :blk false;
            }
            self.factorization.update(.{
                .leaving_row = @intCast(leaving_row),
                .entering_col = @intCast(entering_col),
                .direction = basis.pivot_direction,
                .column_scale = entering_direction,
            }) catch |err| {
                const failure_kind = self.factorization.recordUpdateFailure(err);
                if (self.stats.first_update_failure_kind == null) {
                    self.stats.first_update_failure_kind = failure_kind;
                    self.stats.first_update_failure_iteration = self.iterations;
                    self.stats.first_update_failure_entering = @intCast(entering_col);
                    self.stats.first_update_failure_leaving_row = @intCast(leaving_row);
                }
                forced_rebuild_reason = if (self.cleanup_active) .cleanup else .update_rejected;
                break :blk false;
            };
            break :blk true;
        };
        self.direction_requires_reinversion = false;
        self.numerical.observePivot(pivot);
        for (basis.basic_value, basis.pivot_direction) |*value, direction| value.* -= step * direction;
        const entering_value = basis.primal[entering_col] + entering_direction * step;
        basis.basic_value[leaving_row] = entering_value;
        basis.basic_lower[leaving_row] = basis.col_lower[entering_col];
        basis.basic_upper[leaving_row] = basis.col_upper[entering_col];
        basis.primal[entering_col] = entering_value;
        basis.primal[leaving_col] = if (leaving_bound == .at_upper) basis.col_upper[leaving_col] else basis.col_lower[leaving_col];
        basis.applyPivot(leaving_row, entering_col, leaving_bound) catch return .numerical_failure;
        self.dual_phase_one.notePivot(entering_col, leaving_col, leaving_bound);
        if (self.devex_reset_after_pivot) {
            self.initializePrimalDevexFramework();
            self.devex_reset_after_pivot = false;
        }
        for (basis.basic_index, basis.basic_value) |global_col, value| basis.primal[global_col] = value;
        if (self.pivot_trace_count < self.active_pivot_trace.len) {
            self.active_pivot_trace[self.pivot_trace_count] = .{
                .phase = self.current_phase,
                .iteration = self.iterations,
                .entering_column = @intCast(entering_col),
                .leaving_column = leaving_col,
                .leaving_row = @intCast(leaving_row),
                .pivot = pivot,
                .step = step,
                .update_count = self.factorization.update_count,
                .ftran_relative_residual = self.numerical.last_ftran_relative_residual,
                .condition_estimate = self.numerical.pivot_condition_estimate,
            };
            self.pivot_trace_count += 1;
        }
        const reinversion_reason = self.factorization.reinversionReason(self.numerical.max_update_count);
        if (reinversion_reason) |reason| self.factorization.recordReinversion(reason);
        if (!update_succeeded or reinversion_reason != null or self.numerical.needsRefactor()) {
            const rebuild_reason = forced_rebuild_reason orelse if (reinversion_reason) |reason| switch (reason) {
                .update_limit => RebuildReason.update_limit,
                .update_growth => RebuildReason.update_growth,
                .solve_residual => RebuildReason.solve_residual,
                .small_pivot => RebuildReason.small_pivot,
            } else RebuildReason.numerical_policy;
            const status = self.refactorizeBasis(problem, rebuild_reason);
            if (status == .optimal and rebuild_reason == .fresh_mode)
                self.fresh_factorization_pivots_remaining -|= 1;
            if (status != .optimal and self.failure_site == .none) self.failure_site = .pivot_update;
            return status;
        }
        return .optimal;
    }

    /// Lightweight Devex reference update. It deliberately uses the already
    /// hot FTRAN direction and does not allocate. Exact steepest-edge weights
    /// can replace this policy without changing basis storage or pivot code.
    fn updateLegacyDevexWeights(self: *SimplexEngine, entering_col: usize, leaving_row: usize, pivot: f64) void {
        const basis = if (self.basis) |*value| value else return;
        if (entering_col >= basis.col_edge_weight.len or leaving_row >= basis.row_edge_weight.len) return;
        var norm_squared: f64 = 1.0;
        for (basis.pivot_direction) |value| norm_squared += value * value;
        if (!std.math.isFinite(norm_squared)) norm_squared = 1.0;
        const entering_weight = @max(1.0, norm_squared);
        basis.col_edge_weight[entering_col] = entering_weight;
        basis.row_edge_weight[leaving_row] = if (@abs(pivot) > self.numerical.pivot_tolerance)
            @max(1.0, entering_weight / (pivot * pivot))
        else
            entering_weight;

        if (self.pricing.devex_reset_period != 0 and self.pricing.iterations % self.pricing.devex_reset_period == 0) {
            @memset(basis.col_edge_weight, 1.0);
            @memset(basis.row_edge_weight, 1.0);
        }
    }

    /// Full dual Devex recurrence used by HiGHS-style row pricing. The hot
    /// FTRAN column is reused directly, so updating every row is allocation
    /// free and requires no additional basis solve.
    fn updateDualDevexWeights(self: *SimplexEngine, leaving_row: usize, pivot: f64) void {
        const basis = if (self.basis) |*value| value else return;
        if (leaving_row >= basis.row_edge_weight.len or @abs(pivot) <= self.numerical.pivot_tolerance) return;
        var pivotal_weight = basis.row_edge_weight[leaving_row] / (pivot * pivot);
        if (!std.math.isFinite(pivotal_weight)) pivotal_weight = 1.0;
        pivotal_weight = @max(1.0, pivotal_weight);
        for (basis.row_edge_weight, basis.pivot_direction) |*weight, alpha| {
            const candidate = pivotal_weight * alpha * alpha;
            if (std.math.isFinite(candidate)) weight.* = @max(weight.*, candidate);
        }
        basis.row_edge_weight[leaving_row] = pivotal_weight;
        self.stats.dual_devex_updates += 1;
    }

    /// Freeze the current nonbasic set as a primal Devex reference framework.
    /// Subsequent pivots retain these bits until accumulated weight error
    /// triggers a deterministic reset after the combinatorial pivot commits.
    fn initializePrimalDevexFramework(self: *SimplexEngine) void {
        const basis = if (self.basis) |*value| value else return;
        @memset(basis.col_edge_weight, 1.0);
        for (basis.devex_reference, basis.col_status) |*reference, status|
            reference.* = @intFromBool(status != .basic);
        self.devex_framework_iterations = 0;
        self.devex_bad_weight_count = 0;
        self.stats.devex_frameworks += 1;
    }

    /// Full primal Devex reference recurrence. The pivotal norm is evaluated
    /// from the frozen reference set and the hot FTRAN direction, then every
    /// nonbasic weight touched by the complete old-basis tableau row is raised
    /// to its framework lower bound. No candidate list or allocation is used.
    fn updatePrimalDevexFramework(
        self: *SimplexEngine,
        entering_col: usize,
        leaving_row: usize,
        pivot: f64,
    ) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        if (entering_col >= basis.col_edge_weight.len or leaving_row >= basis.basic_index.len or
            @abs(pivot) <= self.numerical.pivot_tolerance)
            return .numerical_failure;
        var pivotal_norm: f64 = @floatFromInt(basis.devex_reference[entering_col]);
        for (basis.basic_index, basis.pivot_direction) |basic_column, direction| {
            if (basis.devex_reference[basic_column] != 0) pivotal_norm += direction * direction;
        }
        if (!std.math.isFinite(pivotal_norm) or pivotal_norm <= 0.0) return .numerical_failure;
        if (basis.col_edge_weight[entering_col] > 3.0 * pivotal_norm) {
            self.devex_bad_weight_count += 1;
            self.stats.devex_bad_weights += 1;
        }
        const new_pivotal_weight = pivotal_norm / (pivot * pivot);
        if (!std.math.isFinite(new_pivotal_weight)) return .numerical_failure;
        for (basis.col_edge_weight, basis.devex_reference, basis.tableau, basis.col_status) |*weight, reference, alpha, status| {
            if (status == .basic) continue;
            const candidate = new_pivotal_weight * alpha * alpha + @as(f64, @floatFromInt(reference));
            if (!std.math.isFinite(candidate)) return .numerical_failure;
            weight.* = @max(weight.*, candidate);
        }
        const leaving_col: usize = @intCast(basis.basic_index[leaving_row]);
        basis.col_edge_weight[leaving_col] = @max(1.0, new_pivotal_weight);
        basis.col_edge_weight[entering_col] = 1.0;
        self.devex_framework_iterations += 1;
        self.stats.devex_framework_updates += 1;
        if (self.devex_bad_weight_count > 3) self.devex_reset_after_pivot = true;
        return .optimal;
    }

    fn solvePrimal(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
        const phase_started = self.statisticsTimestamp();
        const iteration_started = self.iterations;
        const saved_pricing_rule = self.pricing.rule;
        if (self.active_primal_pricing_strategy == .partial) self.pricing.rule = .partial;
        defer self.pricing.rule = saved_pricing_rule;
        defer self.recordPhaseElapsed(.phase_two, phase_started);
        defer self.recordPhaseIterations(.phase_two, iteration_started);
        if (self.recomputeReducedCosts(problem) != .optimal) {
            self.failure_site = .reduced_cost;
            return .numerical_failure;
        }
        if (self.pricing.rule == .devex and self.active_devex_strategy == .framework)
            self.initializePrimalDevexFramework();
        while (self.iterations < control.max_iterations) : (self.iterations += 1) {
            if (self.beginIteration(problem, control, .phase_two)) |status| return status;
            const basis = if (self.basis) |*value| value else return .numerical_failure;
            const original_cols = problem.num_cols + problem.num_rows;
            const pricing_tolerance = self.scaledDualTolerance(problem);
            if (self.prepareDegeneracyPolicy(basis.basic_value.len, original_cols) != .optimal)
                return .numerical_failure;
            var entering_candidate = if (self.active_degeneracy_strategy == .baseline)
                self.choosePrimalEnteringTimed(original_cols, pricing_tolerance)
            else
                self.choosePrimalEnteringWeightedTimed(original_cols, pricing_tolerance, true);
            if (entering_candidate == null and self.active_degeneracy_strategy == .perturbation_taboo and
                self.degeneracy.active)
            {
                self.degeneracy.invalidateTaboo();
                self.stats.taboo_retries += 1;
                entering_candidate = self.choosePrimalEnteringWeightedTimed(original_cols, pricing_tolerance, false);
            }
            const entering = entering_candidate orelse {
                // Never certify dual feasibility from an updated BTRAN chain.
                // Reinvert the unchanged basis once and price again; a stale
                // transpose solve can otherwise publish a false optimum.
                if (self.factorization.update_count != 0) {
                    self.factorization.recordReinversion(.solve_residual);
                    if (self.refactorizeBasis(problem, .solve_residual) != .optimal) return .numerical_failure;
                    if (self.recomputeReducedCosts(problem) != .optimal) {
                        self.failure_site = .reduced_cost;
                        return .numerical_failure;
                    }
                    continue;
                }
                const status = self.finishOptimal(problem);
                if (status == .numerical_failure) {
                    if (self.active_degeneracy_strategy != .baseline) {
                        self.stats.cold_restart_solves += 1;
                        return self.restartSolveWithoutPerturbation(problem, control);
                    }
                    // finishOptimal still fails at baseline degeneracy.
                    // When the Devex framework changes the pivot path into a
                    // basis that fails the residual check, fall back to the
                    // well-tested legacy pricing with a clean cold start.
                    if (self.active_devex_strategy == .framework) {
                        self.stats.cold_restart_solves += 1;
                        var legacy_control = control;
                        legacy_control.initial_basis = null;
                        legacy_control.degeneracy_strategy = .baseline;
                        legacy_control.devex_strategy = .legacy;
                        return self.solveProblem(problem, legacy_control);
                    }
                    self.failure_site = .optimality_check;
                }
                return status;
            };
            if (self.computeDirection(problem, entering.column) != .optimal) return .numerical_failure;
            if (entering.direction < 0) {
                for (basis.pivot_direction) |*value| value.* = -value.*;
            }
            var leaving = self.chooseLeaving();
            var small_pivot_retry = false;
            // Backward error bounds the full solve but not a tiny pivotal
            // component's forward error. Recompute suspicious pivots from a
            // fresh basis before mutating basis membership.
            if (leaving.status == .optimal and self.pivotNeedsFreshFactorization(leaving.row.?)) {
                small_pivot_retry = true;
                if (self.factorization.update_count != 0) {
                    // Once a basis exposes a component at the forward-accuracy
                    // boundary, occasional FT updates can steer later ratio
                    // tests onto a different numerical path. Recompute from a
                    // fresh basis before deciding whether the pivot is real.
                    self.fresh_factorization_pivots_remaining = self.numerical.fresh_factorization_recovery_pivots;
                    self.factorization.recordReinversion(.small_pivot);
                    if (self.refactorizeBasis(problem, .small_pivot) != .optimal) return .numerical_failure;
                    if (self.computeDirection(problem, entering.column) != .optimal) return .numerical_failure;
                    if (entering.direction < 0) {
                        for (basis.pivot_direction) |*value| value.* = -value.*;
                    }
                    leaving = self.chooseLeaving();
                }
                // A warning that survives fresh FTRAN is not update drift.
                // Keep Bland's entering choice, but let Harris replace the
                // numerically weak leaving row before basis membership is
                // mutated. This prevents a known near-rank pivot from being
                // handed to reinversion as an already committed basis.
                if (leaving.status == .optimal and self.pivotNeedsFreshFactorization(leaving.row.?))
                    leaving = self.chooseLeavingWithPolicy(false);
            }
            const own_step = if (entering.direction > 0)
                basis.col_upper[entering.column] - basis.primal[entering.column]
            else
                basis.primal[entering.column] - basis.col_lower[entering.column];
            if (std.math.isFinite(own_step) and (leaving.status == .unbounded or own_step < leaving.step)) {
                for (basis.basic_value, basis.pivot_direction) |*value, direction| value.* -= own_step * direction;
                basis.primal[entering.column] = if (entering.direction > 0) basis.col_upper[entering.column] else basis.col_lower[entering.column];
                basis.col_status[entering.column] = if (entering.direction > 0) .at_upper else .at_lower;
                self.observeIterationStep(
                    own_step,
                    entering.column,
                    entering.direction,
                    null,
                    0,
                    null,
                    @abs(basis.reduced_cost[entering.column] * own_step),
                    small_pivot_retry,
                    0,
                );
            } else {
                if (leaving.status == .unbounded) {
                    const status = self.finishUnbounded(problem, entering.column, entering.direction);
                    if (status == .numerical_failure and self.active_degeneracy_strategy != .baseline) {
                        self.stats.cold_restart_solves += 1;
                        return self.restartSolveWithoutPerturbation(problem, control);
                    }
                    return status;
                }
                if (leaving.status != .optimal) return leaving.status;
                const leaving_column = basis.basic_index[@intCast(leaving.row.?)];
                const objective_change = @abs(basis.reduced_cost[entering.column] * leaving.step);
                const ratio_tie_count = if ((self.statistics_io != null or self.active_degeneracy_trace.len != 0) and
                    leaving.step <= self.numerical.primal_tolerance)
                    self.countPrimalRatioTies(leaving.step)
                else
                    0;
                if (self.updateReducedCostsAfterPrimalPivot(problem, entering.column, entering.direction, leaving.row.?) != .optimal)
                    return .numerical_failure;
                if (self.performPivot(problem, entering.column, entering.direction, leaving.row.?, leaving.bound, leaving.step) != .optimal) {
                    if (self.failure_site == .none) self.failure_site = .pivot_update;
                    return .numerical_failure;
                }
                self.observeIterationStep(
                    leaving.step,
                    entering.column,
                    entering.direction,
                    leaving_column,
                    ratio_tie_count,
                    own_step,
                    objective_change,
                    small_pivot_retry,
                    0,
                );
            }
            if ((self.fresh_factorization_pivots_remaining != 0 or self.reduced_cost_update_count >= self.reduced_cost_refresh_period) and
                self.recomputeReducedCosts(problem) != .optimal)
            {
                self.failure_site = .reduced_cost;
                return .numerical_failure;
            }
        }
        return .iteration_limit;
    }

    fn scaledDualTolerance(self: *const SimplexEngine, problem: problem_module.ProblemView) f64 {
        const basis = if (self.basis) |*value| value else return self.numerical.dual_tolerance;
        var minimum_column_scale: f64 = 1.0;
        for (basis.column_scale[0..problem.num_cols]) |scale| minimum_column_scale = @min(minimum_column_scale, scale);
        return @max(std.math.floatEps(f64), self.numerical.dual_tolerance * self.objective_scale * minimum_column_scale);
    }

    fn choosePrimalEnteringTimed(self: *SimplexEngine, column_count: usize, tolerance: f64) ?pricing_module.EnteringChoice {
        const started = self.statisticsTimestamp();
        defer self.recordPricingElapsed(started);
        const basis = if (self.basis) |*value| value else return null;
        self.observePricingDensity(basis.reduced_cost[0..column_count]);
        self.stats.dense_pricing_dispatches += 1;
        return if (self.numerical.anti_cycling_active)
            self.pricing.choosePrimalEnteringBland(
                basis.reduced_cost[0..column_count],
                basis.col_status[0..column_count],
                tolerance,
            )
        else if (self.pricing.rule == .partial)
            self.pricing.choosePrimalEnteringMultiple(
                basis.reduced_cost[0..column_count],
                basis.col_status[0..column_count],
                basis.col_edge_weight[0..column_count],
                basis.flip_columns[0..column_count],
                tolerance,
            )
        else
            self.pricing.choosePrimalEnteringWeighted(
                basis.reduced_cost[0..column_count],
                basis.col_status[0..column_count],
                basis.col_edge_weight[0..column_count],
                tolerance,
            );
    }

    fn choosePrimalEnteringWeightedTimed(
        self: *SimplexEngine,
        column_count: usize,
        tolerance: f64,
        respect_taboo: bool,
    ) ?pricing_module.EnteringChoice {
        const started = self.statisticsTimestamp();
        defer self.recordPricingElapsed(started);
        const basis = if (self.basis) |*value| value else return null;
        self.observePricingDensity(basis.reduced_cost[0..column_count]);
        self.stats.dense_pricing_dispatches += 1;
        return if (self.degeneracy.active)
            self.pricing.choosePrimalEnteringPerturbed(
                basis.reduced_cost[0..column_count],
                basis.col_status[0..column_count],
                basis.col_edge_weight[0..column_count],
                self.degeneracy.column_rank[0..column_count],
                if (respect_taboo and self.active_degeneracy_strategy == .perturbation_taboo)
                    self.degeneracy.taboo_until[0..column_count]
                else
                    &.{},
                self.iterations,
                tolerance,
            )
        else if (self.pricing.rule == .partial)
            self.pricing.choosePrimalEnteringMultiple(
                basis.reduced_cost[0..column_count],
                basis.col_status[0..column_count],
                basis.col_edge_weight[0..column_count],
                basis.flip_columns[0..column_count],
                tolerance,
            )
        else
            self.pricing.choosePrimalEnteringWeighted(
                basis.reduced_cost[0..column_count],
                basis.col_status[0..column_count],
                basis.col_edge_weight[0..column_count],
                tolerance,
            );
    }

    fn solveDual(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
        const phase_started = self.statisticsTimestamp();
        const iteration_started = self.iterations;
        const saved_pricing_rule = self.beginDualEdgeWeightPhase();
        defer self.pricing.rule = saved_pricing_rule;
        defer self.recordPhaseElapsed(.phase_two, phase_started);
        defer self.recordPhaseIterations(.phase_two, iteration_started);
        if (self.recomputeReducedCosts(problem) != .optimal) return .numerical_failure;
        while (self.iterations < control.max_iterations) : (self.iterations += 1) {
            if (self.beginIteration(problem, control, .phase_two)) |status| return status;
            const basis = if (self.basis) |*value| value else return .numerical_failure;
            if (self.pricing.rule == .steepest_edge and self.ensureExactDualEdgeWeights() != .optimal)
                return .numerical_failure;
            const leaving = self.chooseDualLeavingRow() orelse return self.finishOptimal(problem);

            if (self.computeDualTableauRow(problem, leaving.row) != .optimal) return .numerical_failure;
            if (self.active_dual_edge_weight_strategy == .steepest_devex and
                self.pricing.rule == .steepest_edge and !self.dual_edge_weights_valid)
                self.switchDualDseToDevex(.invalid);
            const original_cols = problem.num_cols + problem.num_rows;
            var active_ratio_test = self.ratio_test;
            if (self.numerical.anti_cycling_active) active_ratio_test.rule = .standard;
            const entering = active_ratio_test.chooseDualEntering(
                basis.tableau[0..original_cols],
                basis.reduced_cost[0..original_cols],
                basis.col_status[0..original_cols],
                &.{},
                basis.col_lower[0..original_cols],
                basis.col_upper[0..original_cols],
                basis.primal[0..original_cols],
                leaving.bound,
                leaving.violation,
                basis.dual_ratio[0..original_cols],
                basis.dual_direction[0..original_cols],
                basis.flip_columns[0..original_cols],
            );

            if (self.applyBoundFlips(problem, entering.flip_count) != .optimal) return .numerical_failure;
            const entering_col = entering.column orelse return .infeasible;
            const entering_index: usize = @intCast(entering_col);
            if (self.computeDirection(problem, entering_index) != .optimal) return .numerical_failure;
            if (entering.direction < 0.0) {
                for (basis.pivot_direction) |*value| value.* = -value.*;
            }
            const leaving_row: usize = @intCast(leaving.row);
            const target = if (leaving.bound == .at_lower) basis.basic_lower[leaving_row] else basis.basic_upper[leaving_row];
            const pivot = basis.pivot_direction[leaving_row];
            if (@abs(pivot) <= self.numerical.pivot_tolerance) return .numerical_failure;
            const step = (basis.basic_value[leaving_row] - target) / pivot;
            if (!std.math.isFinite(step) or step < -self.numerical.primal_tolerance) return .numerical_failure;
            const leaving_column = basis.basic_index[leaving_row];
            if (self.updateReducedCostsAfterDualPivot(problem, entering_index, leaving.row) != .optimal)
                return .numerical_failure;
            if (self.performPivot(
                problem,
                entering_index,
                entering.direction,
                leaving_row,
                leaving.bound,
                @max(step, 0.0),
            ) != .optimal) return .numerical_failure;
            self.observeIterationStep(
                @max(step, 0.0),
                entering_index,
                entering.direction,
                leaving_column,
                0,
                null,
                null,
                false,
                entering.flip_count,
            );
            if (self.fresh_factorization_pivots_remaining != 0 or
                self.reduced_cost_update_count >= self.reduced_cost_refresh_period)
            {
                if (self.recomputeReducedCostsWithDrift(problem) != .optimal) return .numerical_failure;
                if (!self.classifyFeasibility(problem).dual) return .numerical_failure;
            }
        }
        return .iteration_limit;
    }

    /// Dedicated dual Phase I. Special bounds encode the negated dual
    /// infeasibility objective, while deterministic perturbed costs remove
    /// ratio ties. The borrowed matrix is never copied and every iteration
    /// uses engine-owned work arrays.
    fn solveDualPhaseOne(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
        const phase_started = self.statisticsTimestamp();
        const iteration_started = self.iterations;
        const saved_pricing_rule = self.beginDualEdgeWeightPhase();
        defer self.pricing.rule = saved_pricing_rule;
        defer self.recordPhaseElapsed(.dual_phase_one, phase_started);
        defer self.recordPhaseIterations(.dual_phase_one, iteration_started);
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        const original_count = problem.num_cols + problem.num_rows;
        self.dual_phase_one.begin(basis, original_count) catch return .numerical_failure;
        defer if (self.dual_phase_one.active)
            self.dual_phase_one.restoreOriginalBounds(basis, original_count);

        self.buildDualPhaseOneCosts(problem);
        if (self.recomputeReducedCostsFromWork(problem) != .optimal) return .not_implemented;
        self.shiftDualPhaseOneCosts(problem);
        if (self.recomputeReducedCostsFromWork(problem) != .optimal) return .not_implemented;
        self.dual_phase_one.installWorkingBounds(basis, original_count);
        if (self.recomputeBasicValuesUnchecked(problem) != .optimal) return .not_implemented;
        self.algorithm = .dual_revised;
        self.dual_edge_weights_valid = false;

        while (self.iterations < control.max_iterations) : (self.iterations += 1) {
            if (self.beginIteration(problem, control, .dual_phase_one)) |status| return status;
            if (self.pricing.rule == .steepest_edge and self.ensureExactDualEdgeWeights() != .optimal)
                return .not_implemented;
            const leaving = self.chooseDualLeavingRow() orelse break;
            if (self.computeDualTableauRow(problem, leaving.row) != .optimal) return .not_implemented;
            if (self.active_dual_edge_weight_strategy == .steepest_devex and
                self.pricing.rule == .steepest_edge and !self.dual_edge_weights_valid)
                self.switchDualDseToDevex(.invalid);
            var active_ratio_test = self.ratio_test;
            if (self.numerical.anti_cycling_active) active_ratio_test.rule = .standard;
            const entering = active_ratio_test.chooseDualEntering(
                basis.tableau[0..original_count],
                basis.reduced_cost[0..original_count],
                basis.col_status[0..original_count],
                self.dual_phase_one.nonbasic_move[0..original_count],
                basis.col_lower[0..original_count],
                basis.col_upper[0..original_count],
                basis.primal[0..original_count],
                leaving.bound,
                leaving.violation,
                basis.dual_ratio[0..original_count],
                basis.dual_direction[0..original_count],
                basis.flip_columns[0..original_count],
            );
            self.dual_phase_one.recordRemainingViolation(
                @intCast(leaving.row),
                leaving.violation,
                basis.tableau[0..original_count],
                basis.flip_columns[0..entering.flip_count],
                basis.col_lower[0..original_count],
                basis.col_upper[0..original_count],
            );
            if (entering.column == null)
                self.recordDualPhaseOneNoEntering(leaving, entering.flip_count, original_count);
            if (self.applyBoundFlips(problem, entering.flip_count) != .optimal) return .not_implemented;
            const entering_col = entering.column orelse return .not_implemented;
            const entering_index: usize = @intCast(entering_col);
            if (self.computeDirection(problem, entering_index) != .optimal) return .not_implemented;
            if (entering.direction < 0.0) {
                for (basis.pivot_direction) |*value| value.* = -value.*;
            }
            const leaving_row: usize = @intCast(leaving.row);
            const target = if (leaving.bound == .at_lower) basis.basic_lower[leaving_row] else basis.basic_upper[leaving_row];
            const pivot = basis.pivot_direction[leaving_row];
            if (@abs(pivot) <= self.numerical.pivot_tolerance) return .not_implemented;
            const step = (basis.basic_value[leaving_row] - target) / pivot;
            if (!std.math.isFinite(step) or step < -self.numerical.primal_tolerance) return .not_implemented;
            const leaving_column = basis.basic_index[leaving_row];
            if (self.performPivot(problem, entering_index, entering.direction, leaving_row, leaving.bound, @max(step, 0.0)) != .optimal)
                return .not_implemented;
            self.observeIterationStep(
                @max(step, 0.0),
                entering_index,
                entering.direction,
                leaving_column,
                0,
                null,
                null,
                false,
                entering.flip_count,
            );
            if (self.recomputeReducedCostsFromWork(problem) != .optimal) return .not_implemented;
        }
        if (self.iterations >= control.max_iterations) return .iteration_limit;

        // Remove perturbation and Phase-I bounds before making any claim
        // about the original LP. A fresh factorization and fresh BTRAN price
        // the restored problem independently of the update chain.
        self.dual_phase_one.restoreOriginalBounds(basis, original_count);
        if (self.refactorizeBasis(problem, .cleanup) != .optimal) return .not_implemented;
        if (self.recomputeBasicValuesUnchecked(problem) != .optimal) return .not_implemented;
        if (self.recomputeReducedCosts(problem) != .optimal) return .not_implemented;
        const feasibility = self.classifyFeasibility(problem);
        if (feasibility.dual) return self.solveDual(problem, control);
        if (feasibility.primal) {
            // The dual-phase edge-weight forcing policy must not leak into a
            // primal Phase-II selected after restoring the original model.
            self.pricing.rule = saved_pricing_rule;
            return self.solvePrimal(problem, control);
        }
        return .not_implemented;
    }

    /// Capture the first dual Phase-I row for which the ratio test cannot
    /// select an entering column. This post-decision scan is diagnostic only:
    /// it reads the already materialized ep/tableau and caller-owned buffers.
    fn recordDualPhaseOneNoEntering(
        self: *SimplexEngine,
        leaving: pricing_module.DualLeavingChoice,
        flip_count: usize,
        column_count: usize,
    ) void {
        if (self.dual_phase_one_failure != null) return;
        const basis = if (self.basis) |*value| value else return;
        var diagnostic = DualPhaseOneFailureDiagnostic{
            .iteration = self.iterations,
            .leaving_row = leaving.row,
            .leaving_bound = leaving.bound,
            .violation = leaving.violation,
            .ep_nonzeros = 0,
            .ep_max_abs = 0.0,
            .small_tableau = 0,
            .basic_or_fixed = 0,
            .wrong_pivot_sign = 0,
            .accepted_bound_flips = 0,
            .eligible_unselected = 0,
        };
        for (basis.dual_row, 0..) |value, row| {
            diagnostic.ep_max_abs = @max(diagnostic.ep_max_abs, @abs(value));
            if (@abs(value) <= self.numerical.zero_tolerance) continue;
            diagnostic.ep_nonzeros += 1;
            if (self.dual_phase_one_ep_trace_count < self.active_dual_phase_one_ep_trace.len) {
                self.active_dual_phase_one_ep_trace[self.dual_phase_one_ep_trace_count] = .{
                    .row = @intCast(row),
                    .value = value,
                };
                self.dual_phase_one_ep_trace_count += 1;
            }
        }

        for (0..column_count) |column| {
            const alpha = basis.tableau[column];
            const status = basis.col_status[column];
            var direction: f64 = 0.0;
            var signed_pivot: f64 = 0.0;
            const width = basis.col_upper[column] - basis.col_lower[column];
            const flip_capacity = if (std.math.isFinite(width)) @abs(alpha) * width else 0.0;
            const reason: DualPhaseOneCandidateReason = blk: {
                if (@abs(alpha) <= self.ratio_test.tolerance) {
                    diagnostic.small_tableau += 1;
                    break :blk .small_tableau;
                }
                direction = switch (status) {
                    .at_lower => 1.0,
                    .at_upper => -1.0,
                    .free, .superbasic => if (leaving.bound == .at_lower)
                        (if (alpha < 0.0) 1.0 else -1.0)
                    else
                        (if (alpha > 0.0) 1.0 else -1.0),
                    .basic, .fixed => {
                        diagnostic.basic_or_fixed += 1;
                        break :blk .basic_or_fixed;
                    },
                };
                signed_pivot = alpha * direction;
                const eligible = if (leaving.bound == .at_lower)
                    signed_pivot < -self.ratio_test.tolerance
                else
                    signed_pivot > self.ratio_test.tolerance;
                if (!eligible) {
                    diagnostic.wrong_pivot_sign += 1;
                    break :blk .wrong_pivot_sign;
                }
                var was_flipped = false;
                for (basis.flip_columns[0..flip_count]) |flipped| {
                    if (@as(usize, @intCast(flipped)) == column) {
                        was_flipped = true;
                        break;
                    }
                }
                if (was_flipped) {
                    diagnostic.accepted_bound_flips += 1;
                    break :blk .accepted_bound_flip;
                }
                diagnostic.eligible_unselected += 1;
                break :blk .eligible_unselected;
            };
            if (self.dual_phase_one_candidate_trace_count < self.active_dual_phase_one_candidate_trace.len) {
                self.active_dual_phase_one_candidate_trace[self.dual_phase_one_candidate_trace_count] = .{
                    .column = @intCast(column),
                    .status = status,
                    .tableau = alpha,
                    .direction = direction,
                    .signed_pivot = signed_pivot,
                    .reduced_cost = basis.reduced_cost[column],
                    .lower = basis.col_lower[column],
                    .upper = basis.col_upper[column],
                    .primal = basis.primal[column],
                    .flip_capacity = flip_capacity,
                    .reason = reason,
                };
                self.dual_phase_one_candidate_trace_count += 1;
            }
        }
        self.dual_phase_one_failure = diagnostic;
    }

    /// Deterministic cold-start selector. Automatic mode remains opt-in until
    /// corpus timing establishes a benefit over the frozen primal baseline.
    fn chooseColdPhaseOneStrategy(self: *const SimplexEngine, problem: problem_module.ProblemView) PhaseOneStrategy {
        const basis = if (self.basis) |*value| value else return .primal;
        var primal_count: usize = 0;
        var primal_sum: f64 = 0.0;
        for (basis.basic_value, basis.basic_lower, basis.basic_upper) |value, lower, upper| {
            const amount = @max(@max(lower - value, value - upper), 0.0);
            if (amount > self.numerical.primal_tolerance) {
                primal_count += 1;
                primal_sum += amount;
            }
        }
        var dual_count: usize = 0;
        var free_count: usize = 0;
        var boxed_count: usize = 0;
        for (0..problem.num_cols + problem.num_rows) |column| {
            const lower = basis.col_lower[column];
            const upper = basis.col_upper[column];
            if (!std.math.isFinite(lower) and !std.math.isFinite(upper)) free_count += 1;
            if (std.math.isFinite(lower) and std.math.isFinite(upper) and lower != upper) boxed_count += 1;
            const reduced = basis.reduced_cost[column];
            const infeasible = switch (basis.col_status[column]) {
                .at_lower => reduced < -self.numerical.dual_tolerance,
                .at_upper => reduced > self.numerical.dual_tolerance,
                .free, .superbasic => @abs(reduced) > self.numerical.dual_tolerance,
                else => false,
            };
            if (infeasible) dual_count += 1;
        }
        var coefficient_min = std.math.inf(f64);
        var coefficient_max: f64 = 0.0;
        for (0..problem.num_cols) |column| {
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
                if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                const magnitude = @abs(coefficient * basis.row_scale[row.toUsize()] * basis.column_scale[column]);
                coefficient_min = @min(coefficient_min, magnitude);
                coefficient_max = @max(coefficient_max, magnitude);
            }
        }
        var cost_min = std.math.inf(f64);
        var cost_max: f64 = 0.0;
        for (problem.col_cost) |cost| {
            const magnitude = @abs(cost);
            if (magnitude == 0.0) continue;
            cost_min = @min(cost_min, magnitude);
            cost_max = @max(cost_max, magnitude);
        }
        const coefficient_range = if (coefficient_min < std.math.inf(f64)) coefficient_max / coefficient_min else 1.0;
        const cost_range = if (cost_min < std.math.inf(f64)) cost_max / cost_min else 1.0;
        // Dual Phase I is favoured for many violated logical rows, provided
        // free columns do not dominate its amplified infeasibility objective.
        // Wide sparse models amortize BTRAN row pricing particularly well;
        // highly scaled models retain the more defensive primal path.
        const width_score = if (problem.num_rows == 0) problem.num_cols else problem.num_cols / problem.num_rows;
        const scale_risk = (if (coefficient_range > 1e12) problem.num_rows else 0) +
            (if (cost_range > 1e12) problem.num_rows else 0);
        const risk_score = free_count * 8 + boxed_count + scale_risk;
        const benefit_score = primal_count * 4 + @min(dual_count, problem.num_cols) + width_score;
        if (primal_count > 0 and primal_sum > self.numerical.primal_tolerance and
            problem.num_cols >= problem.num_rows * 2 and benefit_score > risk_score)
            return .dual;
        return .primal;
    }

    fn buildDualPhaseOneCosts(self: *SimplexEngine, problem: problem_module.ProblemView) void {
        const basis = if (self.basis) |*value| value else return;
        const workspace = &self.dual_phase_one;
        const original_count = problem.num_cols + problem.num_rows;
        const maximize = problem.objective_sense == .maximize;
        var maximum: f64 = 0.0;
        var boxed: usize = 0;
        for (0..problem.num_cols) |column| {
            const cost = (if (maximize) -problem.col_cost[column] else problem.col_cost[column]) *
                self.objective_scale * basis.column_scale[column];
            workspace.work_cost[column] = cost;
            maximum = @max(maximum, @abs(cost));
            if (std.math.isFinite(workspace.saved_lower[column]) and std.math.isFinite(workspace.saved_upper[column]) and
                workspace.saved_lower[column] != workspace.saved_upper[column])
                boxed += 1;
        }
        if (maximum > 100.0) maximum = @sqrt(@sqrt(maximum));
        if (problem.num_cols != 0 and boxed * 100 < problem.num_cols) maximum = @min(maximum, 1.0);
        const base = 5e-7 * @max(maximum, 1.0);
        for (0..problem.num_cols) |column| {
            const random = deterministicUnit(column, workspace.basis_epoch);
            const lower = workspace.saved_lower[column];
            const upper = workspace.saved_upper[column];
            var perturbation: f64 = 0.0;
            if (lower != upper) {
                if (std.math.isFinite(lower) and !std.math.isFinite(upper)) {
                    perturbation = random * base;
                } else if (!std.math.isFinite(lower) and std.math.isFinite(upper)) {
                    perturbation = -random * base;
                } else if (std.math.isFinite(lower) and std.math.isFinite(upper)) {
                    perturbation = if (workspace.work_cost[column] >= 0.0) random * base else -random * base;
                }
            }
            workspace.perturbation[column] = perturbation;
            workspace.work_cost[column] += perturbation;
        }
        for (problem.num_cols..original_count) |column| {
            const perturbation = deterministicUnit(column, workspace.basis_epoch) * 1e-12;
            workspace.perturbation[column] = perturbation;
            workspace.work_cost[column] = perturbation;
        }
    }

    /// Shift only currently infeasible nonbasic costs to the feasibility
    /// boundary. Unlike the former zero-objective repair, feasible reduced
    /// costs and all basic costs remain intact, so Phase I retains a meaningful
    /// dual objective and deterministic ratio ordering.
    fn shiftDualPhaseOneCosts(self: *SimplexEngine, problem: problem_module.ProblemView) void {
        const basis = if (self.basis) |*value| value else return;
        const workspace = &self.dual_phase_one;
        const original_count = problem.num_cols + problem.num_rows;
        for (0..original_count) |column| {
            if (basis.col_status[column] == .basic or basis.col_status[column] == .fixed) continue;
            const reduced = basis.reduced_cost[column];
            const shift = switch (basis.col_status[column]) {
                .at_lower => if (reduced < self.numerical.dual_tolerance)
                    self.numerical.dual_tolerance - reduced
                else
                    0.0,
                .at_upper => if (reduced > -self.numerical.dual_tolerance)
                    -self.numerical.dual_tolerance - reduced
                else
                    0.0,
                .free, .superbasic => -reduced,
                else => 0.0,
            };
            workspace.dual_infeasibility[column] = @abs(shift);
            workspace.work_cost[column] += shift;
            workspace.perturbation[column] += shift;
        }
    }

    fn deterministicUnit(column: usize, epoch: u64) f64 {
        var value = @as(u64, @intCast(column)) +% epoch *% 0x9e3779b97f4a7c15;
        value = (value ^ (value >> 30)) *% 0xbf58476d1ce4e5b9;
        value = (value ^ (value >> 27)) *% 0x94d049bb133111eb;
        value ^= value >> 31;
        return (@as(f64, @floatFromInt(value >> 11)) + 1.0) * 0x1.0p-53;
    }

    fn recomputeReducedCostsFromWork(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        const costs = self.dual_phase_one.work_cost;
        const original_count = problem.num_cols + problem.num_rows;
        if (costs.len < original_count) return .numerical_failure;
        for (basis.basic_index, 0..) |column, row| basis.dual[row] = if (column < original_count) costs[column] else 0.0;
        self.factorization.solveTranspose(basis.dual) catch return .numerical_failure;
        const pricing_started = self.statisticsTimestamp();
        defer self.recordPricingElapsed(pricing_started);
        for (0..problem.num_cols) |column| {
            var reduced = costs[column];
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
                if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                const row_index = row.toUsize();
                reduced -= basis.row_scale[row_index] * coefficient * basis.column_scale[column] * basis.dual[row_index];
            }
            basis.reduced_cost[column] = reduced;
        }
        for (0..problem.num_rows) |row| basis.reduced_cost[problem.num_cols + row] = costs[problem.num_cols + row] - basis.dual[row];
        @memset(basis.reduced_cost[original_count..], 0.0);
        self.reduced_cost_update_count = 0;
        return .optimal;
    }

    fn chooseDualLeavingRow(self: *SimplexEngine) ?pricing_module.DualLeavingChoice {
        const pricing_started = self.statisticsTimestamp();
        defer self.recordRowPricingElapsed(pricing_started);
        const basis = if (self.basis) |*value| value else return null;
        if (self.pricing.rule == .hyper_sparse and self.dual_hyper_sparse_active) {
            self.stats.hyper_pricing_dispatches += 1;
        } else {
            self.stats.dense_pricing_dispatches += 1;
        }
        if (self.numerical.anti_cycling_active) return self.chooseDualLeavingBland();
        if (self.pricing.rule != .hyper_sparse or !self.dual_hyper_sparse_active)
            return self.pricing.chooseDualLeavingWeighted(
                basis.basic_value,
                basis.basic_lower,
                basis.basic_upper,
                basis.row_edge_weight,
                self.numerical.primal_tolerance,
            );

        var best = self.bestDualCandidate();
        if (best == null or self.dualCandidateScore(best.?.row) + self.numerical.primal_tolerance < self.dual_candidate_cutoff) {
            self.rebuildDualCandidateList();
            best = self.bestDualCandidate();
        }
        return best;
    }

    fn chooseDualLeavingBland(self: *const SimplexEngine) ?pricing_module.DualLeavingChoice {
        const basis = if (self.basis) |*value| value else return null;
        var best: ?pricing_module.DualLeavingChoice = null;
        var best_basic_column: u32 = std.math.maxInt(u32);
        for (0..basis.num_rows) |row| {
            const candidate = self.dualCandidate(row) orelse continue;
            const basic_column = basis.basic_index[row];
            if (basic_column < best_basic_column) {
                best_basic_column = basic_column;
                best = candidate;
            }
        }
        return best;
    }

    fn dualCandidate(self: *const SimplexEngine, row: usize) ?pricing_module.DualLeavingChoice {
        const basis = if (self.basis) |*value| value else return null;
        const value = basis.basic_value[row];
        if (value < basis.basic_lower[row] - self.numerical.primal_tolerance)
            return .{ .row = @intCast(row), .bound = .at_lower, .violation = basis.basic_lower[row] - value };
        if (value > basis.basic_upper[row] + self.numerical.primal_tolerance)
            return .{ .row = @intCast(row), .bound = .at_upper, .violation = value - basis.basic_upper[row] };
        return null;
    }

    fn dualCandidateScore(self: *const SimplexEngine, row_u32: u32) f64 {
        const row: usize = @intCast(row_u32);
        const basis = if (self.basis) |*value| value else return 0.0;
        const candidate = self.dualCandidate(row) orelse return 0.0;
        return candidate.violation / @sqrt(@max(basis.row_edge_weight[row], 1.0));
    }

    fn bestDualCandidate(self: *SimplexEngine) ?pricing_module.DualLeavingChoice {
        const basis = if (self.basis) |*value| value else return null;
        var best: ?pricing_module.DualLeavingChoice = null;
        var best_score = self.numerical.primal_tolerance;
        for (basis.dual_candidate_rows[0..self.dual_candidate_count], basis.dual_candidate_score[0..self.dual_candidate_count]) |row, *stored_score| {
            const score = self.dualCandidateScore(row);
            stored_score.* = score;
            if (score > best_score) {
                best_score = score;
                best = self.dualCandidate(@intCast(row));
            }
        }
        return best;
    }

    fn rebuildDualCandidateList(self: *SimplexEngine) void {
        const basis = if (self.basis) |*value| value else return;
        const capacity = @min(basis.num_rows, 32);
        self.dual_candidate_count = 0;
        self.dual_candidate_cutoff = 0.0;
        if (capacity == 0) return;
        for (0..basis.num_rows) |row| {
            const score = self.dualCandidateScore(@intCast(row));
            if (score <= self.numerical.primal_tolerance) continue;
            if (self.dual_candidate_count < capacity) {
                const slot = self.dual_candidate_count;
                basis.dual_candidate_rows[slot] = @intCast(row);
                basis.dual_candidate_score[slot] = score;
                self.dual_candidate_count += 1;
                continue;
            }
            var weakest: usize = 0;
            for (basis.dual_candidate_score[1..capacity], 1..) |candidate_score, slot| {
                if (candidate_score < basis.dual_candidate_score[weakest]) weakest = slot;
            }
            if (score > basis.dual_candidate_score[weakest]) {
                basis.dual_candidate_rows[weakest] = @intCast(row);
                basis.dual_candidate_score[weakest] = score;
            }
        }
        if (self.dual_candidate_count > 0) {
            self.dual_candidate_cutoff = basis.dual_candidate_score[0];
            for (basis.dual_candidate_score[1..self.dual_candidate_count]) |score|
                self.dual_candidate_cutoff = @min(self.dual_candidate_cutoff, score);
        }
    }

    /// Initialize exact dual steepest-edge weights after a crash, imported
    /// basis, reinversion, or detected recurrence drift.
    fn ensureExactDualEdgeWeights(self: *SimplexEngine) SolveStatus {
        if (self.dual_edge_weights_valid) return .optimal;
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        for (basis.row_edge_weight, 0..) |*weight, row| {
            @memset(basis.dual_row, 0.0);
            basis.dual_row[row] = 1.0;
            // Unit-vector BTRAN: pass the single nonzero index for
            // sparse-index adaptive dispatch (hyper-sparse L solve).
            const ep_index = [_]u32{@intCast(row)};
            self.factorization.solveTransposeSparse(basis.dual_row, &ep_index) catch return .numerical_failure;
            var norm_squared: f64 = 0.0;
            for (basis.dual_row) |entry| norm_squared += entry * entry;
            if (!std.math.isFinite(norm_squared)) return .numerical_failure;
            weight.* = @max(norm_squared, self.numerical.zero_tolerance);
        }
        self.dual_edge_weights_valid = true;
        self.dual_row_index = null;
        return .optimal;
    }

    const DualDseFallbackReason = enum { invalid, budget };

    /// Start a fresh DSE framework at a dual-phase boundary. Returning the
    /// previous pricing rule lets callers restore policy without retaining a
    /// second engine-wide rule or introducing a hot-loop branch.
    fn beginDualEdgeWeightPhase(self: *SimplexEngine) pricing_module.PricingRule {
        const saved = self.pricing.rule;
        if (self.active_dual_edge_weight_strategy == .steepest_devex) {
            self.pricing.rule = .steepest_edge;
            self.dual_edge_weights_valid = false;
            self.dual_row_index = null;
            self.dual_dse_updates_since_start = 0;
        }
        return saved;
    }

    /// Replace a rejected or over-budget DSE framework with a deterministic
    /// unit-weight dual Devex framework. Subsequent pivots use the full Devex
    /// recurrence and never pay for another exact all-row DSE initialization
    /// in the current phase.
    fn switchDualDseToDevex(self: *SimplexEngine, reason: DualDseFallbackReason) void {
        const basis = if (self.basis) |*value| value else return;
        self.pricing.rule = .devex;
        @memset(basis.row_edge_weight, 1.0);
        self.dual_edge_weights_valid = false;
        self.dual_row_index = null;
        self.dual_candidate_count = 0;
        switch (reason) {
            .invalid => self.stats.dual_dse_invalid_fallbacks += 1,
            .budget => self.stats.dual_dse_budget_fallbacks += 1,
        }
    }

    /// Forrest--Goldfarb DSE recurrence. `dual_row` is the freshly computed
    /// BTRAN result B^-T e_p. One additional FTRAN forms
    /// tau = B^-1 B^-T e_p before the factor update changes B.
    fn updateDualSteepestEdgeWeights(self: *SimplexEngine, leaving_row: usize, pivot: f64) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        if (!self.dual_edge_weights_valid or leaving_row >= basis.row_edge_weight.len or
            @abs(pivot) <= self.numerical.pivot_tolerance)
            return .numerical_failure;
        var exact_pivot_weight: f64 = 0.0;
        for (basis.dual_row) |entry| exact_pivot_weight += entry * entry;
        if (!std.math.isFinite(exact_pivot_weight) or exact_pivot_weight <= self.numerical.zero_tolerance)
            return .numerical_failure;
        basis.row_edge_weight[leaving_row] = exact_pivot_weight;
        @memcpy(basis.residual_work, basis.dual_row);
        self.factorization.solve(basis.residual_work) catch return .numerical_failure;

        for (basis.row_edge_weight, basis.pivot_direction, basis.residual_work, 0..) |*weight, alpha, tau, row| {
            if (row == leaving_row) continue;
            const ratio = alpha / pivot;
            const updated = weight.* - 2.0 * ratio * tau + ratio * ratio * exact_pivot_weight;
            if (!std.math.isFinite(updated)) return .numerical_failure;
            weight.* = @max(updated, 1e-4);
        }
        basis.row_edge_weight[leaving_row] = @max(exact_pivot_weight / (pivot * pivot), 1e-4);
        return .optimal;
    }

    /// Repair a warm basis that is neither primal nor dual feasible by using
    /// the zero auxiliary objective. Every reduced cost is then exactly zero,
    /// so dual pivots can restore primal feasibility without discarding the
    /// imported basis or factorization. The original objective is restored by
    /// the caller before entering Phase II.
    fn repairWarmBasisWithDual(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
        const phase_started = self.statisticsTimestamp();
        const iteration_started = self.iterations;
        defer self.recordPhaseElapsed(.dual_feasibility_repair, phase_started);
        defer self.recordPhaseIterations(.dual_feasibility_repair, iteration_started);
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        self.algorithm = .dual_revised;
        @memset(basis.reduced_cost, 0.0);
        while (self.iterations < control.max_iterations) : (self.iterations += 1) {
            if (self.beginIteration(problem, control, .dual_feasibility_repair)) |status| return status;
            const leaving = self.pricing.chooseDualLeavingWeighted(
                basis.basic_value,
                basis.basic_lower,
                basis.basic_upper,
                basis.row_edge_weight,
                self.numerical.primal_tolerance,
            ) orelse return .optimal;
            if (self.computeDualTableauRow(problem, leaving.row) != .optimal) return .numerical_failure;
            const original_cols = problem.num_cols + problem.num_rows;
            var active_ratio_test = self.ratio_test;
            if (self.numerical.anti_cycling_active) active_ratio_test.rule = .standard;
            const entering = active_ratio_test.chooseDualEntering(
                basis.tableau[0..original_cols],
                basis.reduced_cost[0..original_cols],
                basis.col_status[0..original_cols],
                &.{},
                basis.col_lower[0..original_cols],
                basis.col_upper[0..original_cols],
                basis.primal[0..original_cols],
                leaving.bound,
                leaving.violation,
                basis.dual_ratio[0..original_cols],
                basis.dual_direction[0..original_cols],
                basis.flip_columns[0..original_cols],
            );
            if (self.applyBoundFlips(problem, entering.flip_count) != .optimal) return .numerical_failure;
            const entering_col = entering.column orelse return .infeasible;
            const entering_index: usize = @intCast(entering_col);
            if (self.computeDirection(problem, entering_index) != .optimal) return .numerical_failure;
            if (entering.direction < 0.0) {
                for (basis.pivot_direction) |*value| value.* = -value.*;
            }
            const leaving_row: usize = @intCast(leaving.row);
            const target = if (leaving.bound == .at_lower) basis.basic_lower[leaving_row] else basis.basic_upper[leaving_row];
            const pivot = basis.pivot_direction[leaving_row];
            if (@abs(pivot) <= self.numerical.pivot_tolerance) return .numerical_failure;
            const step = (basis.basic_value[leaving_row] - target) / pivot;
            if (!std.math.isFinite(step) or step < -self.numerical.primal_tolerance) return .numerical_failure;
            const leaving_column = basis.basic_index[leaving_row];
            if (self.performPivot(
                problem,
                entering_index,
                entering.direction,
                leaving_row,
                leaving.bound,
                @max(step, 0.0),
            ) != .optimal) return .numerical_failure;
            self.observeIterationStep(
                @max(step, 0.0),
                entering_index,
                entering.direction,
                leaving_column,
                0,
                null,
                null,
                false,
                entering.flip_count,
            );
            @memset(basis.reduced_cost, 0.0);
        }
        return .iteration_limit;
    }

    fn computeDualTableauRow(self: *SimplexEngine, problem: problem_module.ProblemView, leaving_row: u32) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        const row: usize = @intCast(leaving_row);
        if (row >= problem.num_rows) return .numerical_failure;
        @memset(basis.dual_row, 0.0);
        basis.dual_row[row] = 1.0;
        // Unit-vector BTRAN for pivotal row: wire sparse-index dispatch.
        const ep_index = [_]u32{@intCast(row)};
        self.factorization.solveTransposeSparse(basis.dual_row, &ep_index) catch return .numerical_failure;
        self.observeEpDensity(basis.dual_row);
        if (self.pricing.rule == .steepest_edge and self.dual_edge_weights_valid) {
            var exact_weight: f64 = 0.0;
            for (basis.dual_row) |entry| exact_weight += entry * entry;
            if (!std.math.isFinite(exact_weight) or exact_weight <= self.numerical.zero_tolerance)
                return .numerical_failure;
            const updated_weight = basis.row_edge_weight[row];
            const relative_error = @abs(updated_weight - exact_weight) / @max(exact_weight, self.numerical.zero_tolerance);
            if (relative_error > self.numerical.dual_edge_weight_error_tolerance) {
                basis.row_edge_weight[row] = exact_weight;
                self.numerical.dual_edge_weight_corrections += 1;
                // Huangfu's safeguard rejects a seriously underestimated
                // pivotal weight; force a full exact reset after this pivot.
                if (updated_weight < 0.5 * exact_weight) self.dual_edge_weights_valid = false;
            }
        }
        self.dual_row_index = leaving_row;
        @memset(basis.tableau, 0.0);
        const row_pricing = self.selectRowPricing(basis.dual_row);
        const pricing_started = self.statisticsTimestamp();
        defer if (row_pricing)
            self.recordRowPricingElapsed(pricing_started)
        else
            self.recordPricingElapsed(pricing_started);
        if (row_pricing) {
            for (basis.dual_row, 0..) |dual_value, row_index| {
                if (@abs(dual_value) <= self.numerical.zero_tolerance) continue;
                const scaled_dual = dual_value * basis.row_scale[row_index];
                for (self.pricing_row_view.rowColumns(row_index), self.pricing_row_view.rowValues(row_index)) |column_u32, coefficient| {
                    const column: usize = @intCast(column_u32);
                    basis.tableau[column] += scaled_dual * coefficient * basis.column_scale[column];
                }
            }
        } else {
            for (0..problem.num_cols) |column| {
                const begin = problem.matrix.col_starts[column];
                const end = problem.matrix.col_starts[column + 1];
                var alpha: f64 = 0.0;
                for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |matrix_row, coefficient| {
                    if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                    const row_index = matrix_row.toUsize();
                    alpha += basis.dual_row[row_index] * basis.row_scale[row_index] * coefficient * basis.column_scale[column];
                }
                basis.tableau[column] = alpha;
            }
        }
        for (0..problem.num_rows) |logical_row| {
            basis.tableau[problem.num_cols + logical_row] = basis.dual_row[logical_row];
        }
        var nonzero_count: usize = 0;
        for (basis.tableau[0 .. problem.num_cols + problem.num_rows]) |value| {
            if (@abs(value) > self.numerical.zero_tolerance) nonzero_count += 1;
        }
        const tableau_count = problem.num_cols + problem.num_rows;
        self.dual_hyper_sparse_active = tableau_count > 0 and nonzero_count * 10 < tableau_count;
        if (!self.dual_hyper_sparse_active) self.dual_candidate_count = 0;
        return .optimal;
    }

    fn applyBoundFlips(self: *SimplexEngine, problem: problem_module.ProblemView, flip_count: usize) SolveStatus {
        if (flip_count == 0) return .optimal;
        self.stats.bound_flips += flip_count;
        self.stats.bound_flip_batches += 1;
        self.stats.bound_flip_ftran_savings += flip_count - 1;
        const basis = if (self.basis) |*value| value else return .numerical_failure;

        // Accumulate every bound displacement in equation space, then apply
        // B^-1 once. A pivotal tableau row cannot update the whole basic
        // vector; batching the RHS is the exact allocation-free equivalent of
        // one FTRAN per flipped column.
        if (self.accumulateBoundFlipRhs(problem, flip_count, basis.residual_work) != .optimal)
            return .numerical_failure;
        @memcpy(basis.rhs_work, basis.residual_work);
        self.factorization.solve(basis.rhs_work) catch return .numerical_failure;
        for (basis.rhs_work) |value| if (!std.math.isFinite(value)) return .numerical_failure;
        if (!self.boundFlipResidualAcceptable(problem)) {
            self.factorization.recordReinversion(.solve_residual);
            if (self.refactorizeBasis(problem, .solve_residual) != .optimal) return .numerical_failure;
            // Reinversion and its validation reuse row workspaces, so rebuild
            // the immutable aggregate before retrying the exact same batch.
            if (self.accumulateBoundFlipRhs(problem, flip_count, basis.residual_work) != .optimal)
                return .numerical_failure;
            @memcpy(basis.rhs_work, basis.residual_work);
            self.factorization.solve(basis.rhs_work) catch return .numerical_failure;
            for (basis.rhs_work) |value| if (!std.math.isFinite(value)) return .numerical_failure;
            if (!self.boundFlipResidualAcceptable(problem)) return .numerical_failure;
        }

        // Commit only after the aggregate solve succeeds. Thus an invalid
        // width, column, or numerical solve never publishes a partial batch.
        for (basis.basic_value, basis.rhs_work) |*value, displacement| value.* -= displacement;
        for (basis.flip_columns[0..flip_count]) |column_u32| {
            const column: usize = @intCast(column_u32);
            const delta = switch (basis.col_status[column]) {
                .at_lower => basis.col_upper[column] - basis.col_lower[column],
                .at_upper => -(basis.col_upper[column] - basis.col_lower[column]),
                else => return .numerical_failure,
            };
            basis.primal[column] += delta;
            basis.col_status[column] = if (delta > 0.0) .at_upper else .at_lower;
            self.dual_phase_one.noteBoundFlip(column);
        }
        for (basis.basic_index, basis.basic_value) |basic_column, value| basis.primal[basic_column] = value;
        return .optimal;
    }

    /// Check the aggregate solve against the current basis before any flip is
    /// committed. The same normwise backward-error scale used by entering
    /// FTRAN protects long FT chains without adding another solve.
    fn boundFlipResidualAcceptable(self: *SimplexEngine, problem: problem_module.ProblemView) bool {
        const basis = if (self.basis) |*value| value else return false;
        @memcpy(basis.pivot_direction, basis.residual_work);
        for (basis.basic_margin, basis.residual_work) |*magnitude, rhs| magnitude.* = @abs(rhs);
        self.subtractBasisProductWithMagnitude(
            problem,
            basis.rhs_work,
            basis.pivot_direction,
            basis.basic_margin,
        ) catch return false;
        var residual_max: f64 = 0.0;
        var equation_scale: f64 = 0.0;
        for (basis.pivot_direction) |value| residual_max = @max(residual_max, @abs(value));
        for (basis.basic_margin) |value| equation_scale = @max(equation_scale, value);
        const relative = residual_max / @max(1.0, equation_scale);
        self.numerical.last_ftran_relative_residual = relative;
        self.numerical.max_ftran_relative_residual = @max(self.numerical.max_ftran_relative_residual, relative);
        return std.math.isFinite(relative) and relative <= self.numerical.residual_tolerance;
    }

    /// Form the combined scaled equation-space displacement for a batch of
    /// boxed nonbasic bound flips. The caller owns `output`; no column is
    /// materialized and no allocation occurs in this hot path.
    fn accumulateBoundFlipRhs(
        self: *SimplexEngine,
        problem: problem_module.ProblemView,
        flip_count: usize,
        output: []f64,
    ) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        if (output.len != problem.num_rows or flip_count > basis.flip_columns.len) return .numerical_failure;
        @memset(output, 0.0);
        for (basis.flip_columns[0..flip_count]) |column_u32| {
            const column: usize = @intCast(column_u32);
            if (column >= problem.num_cols + problem.num_rows) return .numerical_failure;
            const delta = switch (basis.col_status[column]) {
                .at_lower => basis.col_upper[column] - basis.col_lower[column],
                .at_upper => -(basis.col_upper[column] - basis.col_lower[column]),
                else => return .numerical_failure,
            };
            if (!std.math.isFinite(delta)) return .numerical_failure;
            if (column < problem.num_cols) {
                const column_scale = basis.column_scale[column];
                const begin = problem.matrix.col_starts[column];
                const end = problem.matrix.col_starts[column + 1];
                for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
                    if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                    const row_index = row.toUsize();
                    if (row_index >= output.len) return .numerical_failure;
                    output[row_index] += delta * coefficient * basis.row_scale[row_index] * column_scale;
                }
            } else {
                output[column - problem.num_cols] += delta;
            }
        }
        for (output) |value| if (!std.math.isFinite(value)) return .numerical_failure;
        return .optimal;
    }

    fn solvePhaseOne(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
        const phase_started = self.statisticsTimestamp();
        const iteration_started = self.iterations;
        const saved_pricing_rule = self.pricing.rule;
        const phase_two_refresh_period = self.reduced_cost_refresh_period;
        self.pricing.rule = switch (control.phase_one_pricing) {
            .inherit => if (self.active_primal_pricing_strategy == .partial) .partial else saved_pricing_rule,
            .dantzig => .dantzig,
            .devex => .devex,
            .steepest_edge => .steepest_edge,
        };
        defer self.pricing.rule = saved_pricing_rule;
        // Phase-I drift observations do not justify relaxing the original
        // objective's refresh cadence across the phase boundary.
        defer self.reduced_cost_refresh_period = phase_two_refresh_period;
        defer self.recordPhaseElapsed(.phase_one, phase_started);
        defer self.recordPhaseIterations(.phase_one, iteration_started);
        if (self.pricing.rule == .devex and self.active_devex_strategy == .framework)
            self.initializePrimalDevexFramework();
        if (self.recomputePhaseOneReducedCosts(problem) != .optimal) return .numerical_failure;
        while (self.iterations < control.max_iterations) : (self.iterations += 1) {
            if (self.beginIteration(problem, control, .phase_one)) |status| return status;
            const basis = if (self.basis) |*value| value else return .numerical_failure;
            if (self.prepareDegeneracyPolicy(basis.basic_value.len, basis.col_status.len) != .optimal)
                return .numerical_failure;
            // General bounded Phase I contains artificial columns and is not
            // the standard-form setting required by Bland's proof. Retain the
            // numerically stable Harris candidate set throughout this phase.
            var entering = self.choosePrimalEnteringWeightedTimed(
                basis.reduced_cost.len,
                self.numerical.dual_tolerance,
                true,
            );
            // A taboo cache is a search-order hint, never a feasibility
            // certificate. If it masks every improving column, expire the
            // cache and repeat the exact same pricing scan without exclusions
            // before allowing Phase I to terminate.
            if (entering == null and self.active_degeneracy_strategy == .perturbation_taboo and
                self.degeneracy.active)
            {
                self.degeneracy.invalidateTaboo();
                self.stats.taboo_retries += 1;
                entering = self.choosePrimalEnteringWeightedTimed(
                    basis.reduced_cost.len,
                    self.numerical.dual_tolerance,
                    false,
                );
            }
            const entering_choice = entering orelse {
                if (self.factorization.update_count != 0) {
                    self.factorization.recordReinversion(.solve_residual);
                    if (self.refactorizeBasis(problem, .solve_residual) != .optimal) return .numerical_failure;
                    if (self.recomputePhaseOneReducedCosts(problem) != .optimal) return .numerical_failure;
                    continue;
                }
                break;
            };
            if (self.computeDirection(problem, entering_choice.column) != .optimal) return .numerical_failure;
            if (entering_choice.direction < 0) {
                for (basis.pivot_direction) |*value| value.* = -value.*;
            }
            const leaving = self.chooseLeavingWithPolicy(false);
            const own_step = if (entering_choice.direction > 0)
                basis.col_upper[entering_choice.column] - basis.primal[entering_choice.column]
            else
                basis.primal[entering_choice.column] - basis.col_lower[entering_choice.column];
            if (std.math.isFinite(own_step) and (leaving.status == .unbounded or own_step < leaving.step)) {
                for (basis.basic_value, basis.pivot_direction) |*value, direction| value.* -= own_step * direction;
                basis.primal[entering_choice.column] += entering_choice.direction * own_step;
                basis.col_status[entering_choice.column] = if (entering_choice.direction > 0) .at_upper else .at_lower;
                self.observeIterationStep(
                    own_step,
                    entering_choice.column,
                    entering_choice.direction,
                    null,
                    0,
                    null,
                    @abs(basis.reduced_cost[entering_choice.column] * own_step),
                    false,
                    0,
                );
            } else {
                if (leaving.status != .optimal) return leaving.status;
                const leaving_column = basis.basic_index[@intCast(leaving.row.?)];
                const objective_change = @abs(basis.reduced_cost[entering_choice.column] * leaving.step);
                const ratio_tie_count = if ((self.statistics_io != null or self.active_degeneracy_trace.len != 0) and
                    leaving.step <= self.numerical.primal_tolerance)
                    self.countPrimalRatioTies(leaving.step)
                else
                    0;
                if (self.updateReducedCostsAfterPrimalPivot(problem, entering_choice.column, entering_choice.direction, leaving.row.?) != .optimal)
                    return .numerical_failure;
                if (self.performPivot(problem, entering_choice.column, entering_choice.direction, leaving.row.?, leaving.bound, leaving.step) != .optimal)
                    return .numerical_failure;
                self.observeIterationStep(
                    leaving.step,
                    entering_choice.column,
                    entering_choice.direction,
                    leaving_column,
                    ratio_tie_count,
                    own_step,
                    objective_change,
                    false,
                    0,
                );
            }
            const refresh_period = if (self.active_adaptive_reprice and self.degeneracy.active and
                self.numerical.consecutive_degenerate_pivots >= 32)
                @min(self.reduced_cost_refresh_period, 4)
            else
                self.reduced_cost_refresh_period;
            if (self.fresh_factorization_pivots_remaining != 0 or self.reduced_cost_update_count >= refresh_period) {
                const refresh_status = if (self.active_adaptive_reprice)
                    self.recomputePhaseOneReducedCostsWithDrift(problem)
                else
                    self.recomputePhaseOneReducedCosts(problem);
                if (refresh_status != .optimal) return .numerical_failure;
            }
        }
        if (self.iterations >= control.max_iterations) return .iteration_limit;
        if (self.recomputeBasicValues(problem) != .optimal) {
            // A bounded perturbation may leave the terminal basis just beyond
            // an original bound. That invalidates the perturbed path, not the
            // model. Restore the logical epoch so infeasibility is proved by
            // the unperturbed Phase-I path instead of escaping as a numerical
            // failure.
            if (self.degeneracy.ever_active) {
                self.stats.cold_restart_phase_one += 1;
                return self.restartPhaseOneWithoutPerturbation(problem, control);
            }
            return .numerical_failure;
        }
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        const artificial_begin = problem.num_cols + problem.num_rows;
        var infeasibility: f64 = 0.0;
        for (basis.primal[artificial_begin..]) |value| infeasibility += @max(value, 0.0);
        if (infeasibility > self.numerical.primal_tolerance) {
            // Perturbed or taboo-guided pricing may change the path, but an
            // infeasibility conclusion must come from the unperturbed logical
            // epoch. Restore it transactionally and let baseline Phase I
            // either prove infeasibility or continue searching.
            if (self.degeneracy.ever_active) {
                self.stats.cold_restart_phase_one += 1;
                return self.restartPhaseOneWithoutPerturbation(problem, control);
            }
            return .infeasible;
        }

        if (self.cleanupArtificialBasis(problem) != .optimal) {
            if (self.degeneracy.ever_active) {
                self.stats.cold_restart_phase_one += 1;
                return self.restartPhaseOneWithoutPerturbation(problem, control);
            }
            return .numerical_failure;
        }
        @memset(basis.col_upper[artificial_begin..], 0.0);
        for (basis.col_status[artificial_begin..]) |*status| {
            if (status.* != .basic) status.* = .fixed;
        }
        for (basis.basic_index, 0..) |column, row| {
            if (column >= artificial_begin) basis.basic_upper[row] = 0.0;
        }
        self.phase1_needed = false;
        self.numerical.clearAntiCyclingFallback();
        if (self.degeneracy.ever_active) {
            self.degeneracy.clearAfterProgress();
            self.stats.perturbation_cleanups += 1;
            if (self.refactorizeBasis(problem, .cleanup) != .optimal or
                self.recomputeBasicValues(problem) != .optimal or
                self.recomputePhaseOneReducedCosts(problem) != .optimal)
            {
                self.stats.cold_restart_phase_one += 1;
                return self.restartPhaseOneWithoutPerturbation(problem, control);
            }
            var cleanup_infeasibility: f64 = 0.0;
            for (basis.primal[artificial_begin..]) |value| cleanup_infeasibility += @max(value, 0.0);
            if (cleanup_infeasibility > self.numerical.primal_tolerance) {
                self.stats.cold_restart_phase_one += 1;
                return self.restartPhaseOneWithoutPerturbation(problem, control);
            }
        }
        return .optimal;
    }

    /// Transactional fallback for a failed perturbation cleanup. Restore the
    /// validated logical epoch and rerun Phase I with the baseline policy;
    /// no perturbed ordering state or reduced cost can escape to publication.
    fn restartPhaseOneWithoutPerturbation(
        self: *SimplexEngine,
        problem: problem_module.ProblemView,
        control: SolveControl,
    ) SolveStatus {
        self.degeneracy.resetSolve();
        self.active_degeneracy_strategy = .baseline;
        self.numerical.clearAntiCyclingFallback();
        self.direction_requires_reinversion = false;
        self.fresh_factorization_pivots_remaining = 0;
        self.reduced_cost_update_count = 0;
        if (self.initializeProblemStorage(problem) != .optimal) return .numerical_failure;
        self.factorization.factorizeIdentity(problem.num_rows) catch return .numerical_failure;
        self.dual_edge_weights_valid = true;
        self.observeFactorizationStability();
        if (self.initializeLogicalBasicValues(problem) != .optimal) return .numerical_failure;
        if (self.installArtificialPhaseOneBasis(problem) != .optimal) return .numerical_failure;
        if (!self.phase1_needed) return .optimal;
        if (self.refactorizeBasis(problem, .phase_one_setup) != .optimal) return .numerical_failure;
        self.numerical.markRefactorized();
        return self.solvePhaseOne(problem, control);
    }

    /// Terminal certificates are always checked in original coordinates. If
    /// a perturbed path reaches a basis that cannot pass that check, discard
    /// the entire experimental epoch and cold-start the proven baseline
    /// policy. This is deliberately broader than Phase-I cleanup because a
    /// Phase-II perturbation can change basis membership after artificials are
    /// gone, leaving no smaller transactional snapshot to restore.
    fn restartSolveWithoutPerturbation(
        self: *SimplexEngine,
        problem: problem_module.ProblemView,
        control: SolveControl,
    ) SolveStatus {
        var baseline_control = control;
        baseline_control.initial_basis = null;
        baseline_control.degeneracy_strategy = .baseline;
        return self.solveProblem(problem, baseline_control);
    }

    fn controlledStop(self: *SimplexEngine, control: SolveControl) ?SolveStatus {
        if (control.interrupt_flag) |flag| {
            if (flag.load(.acquire)) return .interrupted;
        }
        if (control.work_limit) |limit| {
            if (self.work_used >= limit) return .work_limit;
        }
        if (control.time_limit_ns) |limit| {
            if (limit == 0) return .time_limit;
            if (self.solve_start_ns) |start| {
                const io = self.solve_clock_io orelse return null;
                const now = std.Io.Clock.awake.now(io).nanoseconds;
                if (now >= start and @as(u128, @intCast(now - start)) >= limit) return .time_limit;
            }
        }
        return null;
    }

    fn beginIteration(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl, phase: SolvePhase) ?SolveStatus {
        self.current_phase = phase;
        if (self.controlledStop(control)) |status| return status;
        const callback_due = control.iteration_callback != null and
            self.work_used % @max(control.callback_interval_work, 1) == 0;
        const log_due = control.log_level == .iterations and control.log_callback != null and
            self.work_used % @max(control.log_interval_work, 1) == 0;
        if (callback_due or log_due) {
            const event = self.progressEvent(problem, phase);
            if (callback_due) {
                if (control.iteration_callback.?(event, control.callback_user_data) == .stop) return .interrupted;
            }
            if (log_due) control.log_callback.?(event, control.log_user_data);
        }
        self.work_used = std.math.add(u64, self.work_used, 1) catch std.math.maxInt(u64);
        return null;
    }

    fn progressEvent(self: *const SimplexEngine, problem: problem_module.ProblemView, phase: SolvePhase) ProgressEventView {
        var primal_infeasibility: f64 = 0.0;
        var dual_infeasibility: f64 = 0.0;
        var current_objective = problem.objective_offset;
        if (self.basis) |*basis| {
            for (problem.col_cost, basis.primal[0..problem.num_cols], basis.column_scale[0..problem.num_cols]) |cost, value, scale|
                current_objective += cost * value * scale;
            for (basis.basic_value, basis.basic_lower, basis.basic_upper) |value, lower, upper| {
                primal_infeasibility = @max(primal_infeasibility, @max(lower - value, value - upper));
            }
            for (basis.reduced_cost, basis.col_status) |reduced, status| {
                const violation: f64 = switch (status) {
                    .at_lower => -reduced,
                    .at_upper => reduced,
                    .free, .superbasic => @abs(reduced),
                    .basic, .fixed => 0.0,
                };
                dual_infeasibility = @max(dual_infeasibility, violation);
            }
        }
        return .{
            .phase = phase,
            .algorithm = self.algorithm,
            .iterations = self.iterations,
            .work_used = self.work_used,
            .objective_value = current_objective,
            .primal_infeasibility = @max(primal_infeasibility, 0.0),
            .dual_infeasibility = @max(dual_infeasibility, 0.0),
        };
    }

    /// Count exact-ratio ties after the ratio policy has already selected its
    /// leaving row. Keeping this read-only diagnostic outside RatioTest makes
    /// it impossible for instrumentation to perturb pivot selection.
    fn countPrimalRatioTies(self: *const SimplexEngine, selected_step: f64) u32 {
        const basis = if (self.basis) |*value| value else return 0;
        var count: u32 = 0;
        for (basis.ratio_direction, basis.basic_margin) |direction, margin| {
            if (direction <= self.ratio_test.tolerance) continue;
            const step = @max(margin / direction, 0.0);
            const scale = @max(1.0, @max(@abs(step), @abs(selected_step)));
            if (@abs(step - selected_step) <= self.ratio_test.tolerance * scale) count += 1;
        }
        return count;
    }

    fn prepareDegeneracyPolicy(self: *SimplexEngine, rows: usize, columns: usize) SolveStatus {
        if (self.active_degeneracy_strategy == .baseline or !self.numerical.anti_cycling_active or
            self.degeneracy.active)
            return .optimal;
        // Automatic mode distinguishes a truly long degenerate face from
        // frequent short local ties. Explicit perturbation remains available
        // at the normal anti-cycling threshold for controlled A/B runs.
        if (self.active_degeneracy_strategy == .automatic and
            self.numerical.consecutive_degenerate_pivots < 256)
            return .optimal;
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        self.degeneracy.activate(
            basis.row_scale[0..rows],
            basis.column_scale[0..columns],
            self.numerical.primal_tolerance,
            self.numerical.dual_tolerance,
            self.objective_scale,
        ) catch return .numerical_failure;
        self.stats.perturbation_activations += 1;
        return .optimal;
    }

    fn observeIterationStep(
        self: *SimplexEngine,
        step: f64,
        entering_column: usize,
        entering_direction: f64,
        leaving_column: ?u32,
        ratio_tie_count: u32,
        own_bound_step: ?f64,
        objective_change: ?f64,
        small_pivot_retry: bool,
        bound_flip_count: usize,
    ) void {
        self.numerical.observeStep(step);
        if (std.math.isFinite(step) and step > self.numerical.primal_tolerance) {
            self.degeneracy.clearAfterProgress();
            return;
        }
        if (self.degeneracy.advance(256)) self.stats.perturbation_expirations += 1;

        var fingerprint: u64 = 0;
        var repeated_basis = false;
        // Fingerprinting scans basis membership and is intentionally opt-in.
        // Scalar statistics keep the normal benchmark path O(1) per pivot.
        if (self.active_degeneracy_trace.len != 0 or self.active_degeneracy_strategy == .perturbation_taboo) {
            fingerprint = self.basisFingerprint();
            repeated_basis = self.rememberDegenerateBasis(fingerprint);
            if (repeated_basis and self.active_degeneracy_strategy == .perturbation_taboo)
                self.degeneracy.escalateAfterRepeat();
        }
        const bound_tie = if (own_bound_step) |own_step|
            leaving_column != null and std.math.isFinite(own_step) and
                @abs(own_step - step) <= self.ratio_test.tolerance * @max(1.0, @max(@abs(own_step), @abs(step)))
        else
            false;
        const change = objective_change orelse std.math.inf(f64);
        const reason: DegeneracyReason = if (bound_flip_count != 0)
            .bound_flip
        else if (small_pivot_retry)
            .small_pivot_retry
        else if (repeated_basis)
            .repeated_basis
        else if (bound_tie)
            .bound_tie
        else if (ratio_tie_count > 1)
            .ratio_tie
        else if (self.current_phase == .phase_one and change <= self.numerical.dual_tolerance)
            .phase_one_objective_stall
        else
            .zero_primal_step;
        if (self.active_degeneracy_strategy == .perturbation_taboo) {
            if (leaving_column) |leaving| {
                self.degeneracy.recordTaboo(
                    @intCast(entering_column),
                    leaving,
                    entering_direction,
                    self.iterations,
                    16,
                );
                self.stats.taboo_records += 1;
            }
        }
        switch (reason) {
            .bound_tie => self.stats.degeneracy_bound_ties += 1,
            .ratio_tie => self.stats.degeneracy_ratio_ties += 1,
            .zero_primal_step => self.stats.degeneracy_zero_primal_steps += 1,
            .phase_one_objective_stall => self.stats.degeneracy_phase_one_objective_stalls += 1,
            .repeated_basis => self.stats.degeneracy_repeated_bases += 1,
            .small_pivot_retry => self.stats.degeneracy_small_pivot_retries += 1,
            .bound_flip => self.stats.degeneracy_bound_flips += 1,
        }
        if (self.degeneracy_trace_count < self.active_degeneracy_trace.len) {
            self.active_degeneracy_trace[self.degeneracy_trace_count] = .{
                .phase = self.current_phase,
                .iteration = self.iterations,
                .reason = reason,
                .entering_column = @intCast(entering_column),
                .leaving_column = leaving_column orelse std.math.maxInt(u32),
                .step = step,
                .objective_change = change,
                .basis_fingerprint = fingerprint,
            };
            self.degeneracy_trace_count += 1;
        }
    }

    fn basisFingerprint(self: *const SimplexEngine) u64 {
        const basis = if (self.basis) |*value| value else return 0;
        var fingerprint: u64 = 0xcbf29ce484222325;
        for (basis.basic_index) |column| {
            fingerprint = (fingerprint ^ @as(u64, column)) *% 0x100000001b3;
        }
        for (basis.col_status) |status| {
            const status_bits: u8 = @bitCast(@as(i8, @intFromEnum(status)));
            fingerprint = (fingerprint ^ @as(u64, status_bits)) *% 0x100000001b3;
        }
        return fingerprint;
    }

    fn rememberDegenerateBasis(self: *SimplexEngine, fingerprint: u64) bool {
        for (self.degeneracy_basis_fingerprints[0..self.degeneracy_basis_fingerprint_count]) |previous| {
            if (previous == fingerprint) return true;
        }
        self.degeneracy_basis_fingerprints[self.degeneracy_basis_fingerprint_cursor] = fingerprint;
        self.degeneracy_basis_fingerprint_cursor = (self.degeneracy_basis_fingerprint_cursor + 1) % self.degeneracy_basis_fingerprints.len;
        self.degeneracy_basis_fingerprint_count = @min(
            self.degeneracy_basis_fingerprint_count + 1,
            self.degeneracy_basis_fingerprints.len,
        );
        return false;
    }

    fn startSolveClock(self: *SimplexEngine, control: SolveControl) void {
        self.active_pivot_trace = control.pivot_trace;
        self.pivot_trace_count = 0;
        self.active_degeneracy_trace = control.degeneracy_trace;
        self.degeneracy_trace_count = 0;
        self.degeneracy_basis_fingerprints = @splat(0);
        self.degeneracy_basis_fingerprint_count = 0;
        self.degeneracy_basis_fingerprint_cursor = 0;
        self.dual_phase_one_failure = null;
        self.active_dual_phase_one_candidate_trace = control.dual_phase_one_candidate_trace;
        self.dual_phase_one_candidate_trace_count = 0;
        self.active_dual_phase_one_ep_trace = control.dual_phase_one_ep_trace;
        self.dual_phase_one_ep_trace_count = 0;
        if (control.time_limit_ns == null) {
            self.solve_start_ns = null;
            self.solve_clock_io = null;
            return;
        }
        const io = control.clock_io orelse std.Io.Threaded.global_single_threaded.io();
        self.solve_clock_io = io;
        self.solve_start_ns = std.Io.Clock.awake.now(io).nanoseconds;
    }

    fn resetStatistics(self: *SimplexEngine, control: SolveControl) void {
        const saved_cold_restarts = .{
            self.stats.cold_restart_solves,
            self.stats.cold_restart_phase_one,
        };
        self.stats = .{};
        // Preserve cold-restart counters across recursive solveProblem calls
        // so that restartSolveWithoutPerturbation and restartPhaseOneWithout-
        // Perturbation counts survive the nested resetStatistics.
        self.stats.cold_restart_solves = saved_cold_restarts[0];
        self.stats.cold_restart_phase_one = saved_cold_restarts[1];
        self.statistics_io = if (control.collect_statistics)
            control.clock_io orelse std.Io.Threaded.global_single_threaded.io()
        else
            null;
        self.factorization.resetStatistics(self.statistics_io);
    }

    fn statisticsTimestamp(self: *const SimplexEngine) ?i96 {
        const io = self.statistics_io orelse return null;
        return std.Io.Clock.awake.now(io).nanoseconds;
    }

    fn elapsedSince(self: *const SimplexEngine, started: ?i96) u64 {
        const begin = started orelse return 0;
        const io = self.statistics_io orelse return 0;
        const end = std.Io.Clock.awake.now(io).nanoseconds;
        if (end <= begin) return 0;
        return @intCast(end - begin);
    }

    fn recordPhaseElapsed(self: *SimplexEngine, phase: SolvePhase, started: ?i96) void {
        const elapsed = self.elapsedSince(started);
        const target = switch (phase) {
            .phase_one => &self.stats.phase_one_ns,
            .dual_phase_one => &self.stats.dual_phase_one_ns,
            .dual_feasibility_repair => &self.stats.dual_repair_ns,
            .phase_two => &self.stats.phase_two_ns,
        };
        target.* = std.math.add(u64, target.*, elapsed) catch std.math.maxInt(u64);
    }

    fn recordPhaseIterations(self: *SimplexEngine, phase: SolvePhase, started: usize) void {
        const count = self.iterations -| started;
        const target = switch (phase) {
            .phase_one => &self.stats.phase_one_iterations,
            .dual_phase_one => &self.stats.dual_phase_one_iterations,
            .dual_feasibility_repair => &self.stats.dual_repair_iterations,
            .phase_two => &self.stats.phase_two_iterations,
        };
        target.* = std.math.add(usize, target.*, count) catch std.math.maxInt(usize);
    }

    fn recordPricingElapsed(self: *SimplexEngine, started: ?i96) void {
        self.stats.pricing_calls += 1;
        self.stats.column_pricing_dispatches += 1;
        self.stats.pricing_ns = std.math.add(u64, self.stats.pricing_ns, self.elapsedSince(started)) catch std.math.maxInt(u64);
    }

    fn recordRowPricingElapsed(self: *SimplexEngine, started: ?i96) void {
        self.stats.pricing_calls += 1;
        self.stats.row_pricing_dispatches += 1;
        self.stats.pricing_ns = std.math.add(u64, self.stats.pricing_ns, self.elapsedSince(started)) catch std.math.maxInt(u64);
    }

    /// Choose the matrix orientation once for the complete pricing operation.
    /// Automatic mode estimates actual CSR work from the active row support;
    /// coefficient loops remain branch-free with respect to representation.
    fn selectRowPricing(self: *const SimplexEngine, vector: []const f64) bool {
        switch (self.active_pricing_kernel) {
            .column => return false,
            .row => return self.pricing_row_view.num_rows == vector.len,
            .automatic => {},
        }
        if (self.pricing_row_view.num_rows != vector.len or self.pricing_row_view.nonzeros == 0) return false;
        var touched_entries: usize = 0;
        for (vector, 0..) |value, row| {
            if (@abs(value) <= self.numerical.zero_tolerance) continue;
            touched_entries = std.math.add(usize, touched_entries, self.pricing_row_view.rowDegree(row)) catch return false;
        }
        // Row traversal wins only with a useful work margin, absorbing scatter
        // stores and the support scan that CSC dot products do not pay.
        return touched_entries * 4 <= self.pricing_row_view.nonzeros * 3;
    }

    fn observeAqDensity(self: *SimplexEngine, vector: []const f64) void {
        if (self.statistics_io == null) return;
        self.stats.aq_samples += 1;
        for (vector) |value| if (@abs(value) > self.numerical.zero_tolerance) {
            self.stats.aq_nonzeros += 1;
        };
    }

    fn observeEpDensity(self: *SimplexEngine, vector: []const f64) void {
        if (self.statistics_io == null) return;
        self.stats.ep_samples += 1;
        for (vector) |value| if (@abs(value) > self.numerical.zero_tolerance) {
            self.stats.ep_nonzeros += 1;
        };
    }

    fn observePricingDensity(self: *SimplexEngine, vector: []const f64) void {
        if (self.statistics_io == null) return;
        self.stats.pricing_samples += 1;
        self.stats.pricing_entries += vector.len;
        for (vector) |value| if (@abs(value) > self.numerical.dual_tolerance) {
            self.stats.pricing_nonzeros += 1;
        };
    }

    /// Pivot zero-valued artificial basics out whenever a stable original or
    /// logical column is available. An artificial that remains basic denotes
    /// a rank-redundant row and is fixed at zero until presolve can remove that
    /// row explicitly.
    fn cleanupArtificialBasis(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
        const cleanup_started = self.statisticsTimestamp();
        defer {
            const elapsed = self.elapsedSince(cleanup_started);
            self.stats.cleanup_ns = std.math.add(u64, self.stats.cleanup_ns, elapsed) catch std.math.maxInt(u64);
        }
        self.cleanup_active = true;
        defer self.cleanup_active = false;
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        const artificial_begin = problem.num_cols + problem.num_rows;

        for (0..problem.num_rows) |leaving_row| {
            if (basis.basic_index[leaving_row] < artificial_begin) continue;

            var entering_col: ?usize = null;
            var best_pivot = self.numerical.pivot_tolerance;
            for (0..artificial_begin) |column| {
                if (basis.col_status[column] == .basic) continue;
                if (self.computeDirection(problem, column) != .optimal) return .numerical_failure;
                const candidate = @abs(basis.pivot_direction[leaving_row]);
                if (candidate > best_pivot) {
                    best_pivot = candidate;
                    entering_col = column;
                }
            }

            const column = entering_col orelse continue;
            if (self.computeDirection(problem, column) != .optimal) return .numerical_failure;
            if (self.performPivot(problem, column, 1.0, leaving_row, .fixed, 0.0) != .optimal)
                return .numerical_failure;
        }
        return self.recomputeBasicValues(problem);
    }

    fn recomputePhaseOneReducedCosts(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        const artificial_begin = problem.num_cols + problem.num_rows;
        for (basis.basic_index, 0..) |global_col, row| basis.dual[row] = if (global_col >= artificial_begin) 1.0 else 0.0;
        self.factorization.solveTranspose(basis.dual) catch return .numerical_failure;
        const row_pricing = self.selectRowPricing(basis.dual);
        const pricing_started = self.statisticsTimestamp();
        defer if (row_pricing)
            self.recordRowPricingElapsed(pricing_started)
        else
            self.recordPricingElapsed(pricing_started);

        // Compute c - A^T*y directly from CSC. The previous generic column
        // loop cleared and dotted a rows-long dense vector for every internal
        // column, turning sparse Phase-I pricing into O(columns * rows).
        // Structural columns use one sparse pass; logical and artificial
        // columns have one known nonzero and are written directly.
        if (row_pricing) {
            @memset(basis.reduced_cost[0..problem.num_cols], 0.0);
            for (basis.dual, 0..) |dual_value, row_index| {
                if (@abs(dual_value) <= self.numerical.zero_tolerance) continue;
                const scaled_dual = dual_value * basis.row_scale[row_index];
                for (self.pricing_row_view.rowColumns(row_index), self.pricing_row_view.rowValues(row_index)) |column_u32, coefficient| {
                    const column: usize = @intCast(column_u32);
                    basis.reduced_cost[column] -= scaled_dual * coefficient * basis.column_scale[column];
                }
            }
        } else {
            for (0..problem.num_cols) |column| {
                const begin = problem.matrix.col_starts[column];
                const end = problem.matrix.col_starts[column + 1];
                var reduced: f64 = 0.0;
                for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
                    if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                    const row_index = row.toUsize();
                    reduced -= basis.row_scale[row_index] * coefficient * basis.column_scale[column] * basis.dual[row_index];
                }
                basis.reduced_cost[column] = reduced;
            }
        }
        for (0..problem.num_rows) |row| {
            basis.reduced_cost[problem.num_cols + row] = -basis.dual[row];
            basis.reduced_cost[artificial_begin + row] = 1.0 - basis.artificial_sign[row] * basis.dual[row];
        }
        self.reduced_cost_update_count = 0;
        return .optimal;
    }

    fn recomputePhaseOneReducedCostsWithDrift(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        const count = basis.reduced_cost.len;
        if (self.degeneracy.reduced_cost_snapshot.len < count) return .numerical_failure;
        @memcpy(self.degeneracy.reduced_cost_snapshot[0..count], basis.reduced_cost);
        if (self.recomputePhaseOneReducedCosts(problem) != .optimal) return .numerical_failure;
        var maximum_drift: f64 = 0.0;
        for (basis.reduced_cost, self.degeneracy.reduced_cost_snapshot[0..count]) |exact, updated| {
            maximum_drift = @max(maximum_drift, @abs(exact - updated) / @max(1.0, @abs(exact)));
        }
        self.maximum_reduced_cost_drift = @max(self.maximum_reduced_cost_drift, maximum_drift);
        self.exact_reprices += 1;
        if (maximum_drift > self.numerical.dual_tolerance * 10.0)
            self.reduced_cost_refresh_period = @max(self.reduced_cost_refresh_period / 2, 1)
        else if (maximum_drift < self.numerical.dual_tolerance * 0.1)
            // Eight is the validated fixed-policy ceiling. Adaptive mode may
            // refresh more often after drift or a long degenerate chain, but
            // never weakens the existing numerical safeguard.
            self.reduced_cost_refresh_period = @min(self.reduced_cost_refresh_period * 2, 8);
        return .optimal;
    }

    /// Exact original-objective reprice with drift observation. This is the
    /// safety boundary for dual rank-1 updates: exact values replace the
    /// incrementally maintained vector before feasibility is revalidated.
    fn recomputeReducedCostsWithDrift(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        const count = basis.reduced_cost.len;
        if (self.degeneracy.reduced_cost_snapshot.len < count) return .numerical_failure;
        @memcpy(self.degeneracy.reduced_cost_snapshot[0..count], basis.reduced_cost);
        if (self.recomputeReducedCosts(problem) != .optimal) return .numerical_failure;
        var maximum_drift: f64 = 0.0;
        for (basis.reduced_cost, self.degeneracy.reduced_cost_snapshot[0..count]) |exact, updated| {
            maximum_drift = @max(maximum_drift, @abs(exact - updated) / @max(1.0, @abs(exact)));
        }
        self.maximum_reduced_cost_drift = @max(self.maximum_reduced_cost_drift, maximum_drift);
        self.exact_reprices += 1;
        self.stats.dual_exact_reprices += 1;
        if (self.active_adaptive_reprice) {
            if (maximum_drift > self.numerical.dual_tolerance * 10.0)
                self.reduced_cost_refresh_period = @max(self.reduced_cost_refresh_period / 2, 1)
            else if (maximum_drift < self.numerical.dual_tolerance * 0.1)
                self.reduced_cost_refresh_period = @min(self.reduced_cost_refresh_period * 2, 8);
        }
        return .optimal;
    }

    /// Maintain `r' = r - (r_q / alpha_pq) alpha_p` from the tableau row of
    /// the old basis. `pivot_direction` may have been negated for movement
    /// from an upper bound, but basis replacement still installs the original
    /// `A_q`, so the unsigned tableau pivot is the required denominator.
    fn updateReducedCostsAfterDualPivot(
        self: *SimplexEngine,
        problem: problem_module.ProblemView,
        entering_col: usize,
        leaving_row: u32,
    ) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        const row: usize = @intCast(leaving_row);
        const original_count = problem.num_cols + problem.num_rows;
        if (row >= problem.num_rows or entering_col >= original_count) return .numerical_failure;
        const alpha_pq = basis.tableau[entering_col];
        if (!std.math.isFinite(alpha_pq) or @abs(alpha_pq) <= self.numerical.pivot_tolerance)
            return .numerical_failure;
        const theta = basis.reduced_cost[entering_col] / alpha_pq;
        if (!std.math.isFinite(theta)) return .numerical_failure;

        for (basis.reduced_cost[0..original_count], basis.tableau[0..original_count]) |*reduced, alpha|
            reduced.* -= theta * alpha;
        const artificial_begin = original_count;
        for (0..problem.num_rows) |logical_row|
            basis.reduced_cost[artificial_begin + logical_row] -=
                theta * basis.dual_row[logical_row] * basis.artificial_sign[logical_row];
        basis.reduced_cost[entering_col] = 0.0;
        for (basis.reduced_cost) |value| if (!std.math.isFinite(value)) return .numerical_failure;
        self.reduced_cost_update_count += 1;
        self.stats.dual_reduced_cost_updates += 1;
        return .optimal;
    }

    /// Update all reduced costs after a primal basis replacement using the
    /// pivotal row of the old basis: `r' = r - (r_q / alpha_pq) * alpha_p`.
    /// One BTRAN and one CSC scan replace rebuilding `c_B`, solving for the
    /// full dual vector, and repricing from scratch. A periodic exact refresh
    /// bounds accumulated roundoff while keeping the hot path allocation-free.
    fn updateReducedCostsAfterPrimalPivot(
        self: *SimplexEngine,
        problem: problem_module.ProblemView,
        entering_col: usize,
        entering_direction: f64,
        leaving_row: u32,
    ) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        const row: usize = @intCast(leaving_row);
        if (row >= problem.num_rows or entering_col >= basis.reduced_cost.len) return .numerical_failure;
        const actual_pivot = basis.pivot_direction[row] * entering_direction;
        if (!std.math.isFinite(actual_pivot) or @abs(actual_pivot) <= self.numerical.pivot_tolerance)
            return .numerical_failure;
        const theta = basis.reduced_cost[entering_col] / actual_pivot;
        if (!std.math.isFinite(theta)) return .numerical_failure;

        @memset(basis.residual_work, 0.0);
        basis.residual_work[row] = 1.0;
        self.factorization.solveTranspose(basis.residual_work) catch return .numerical_failure;
        const row_pricing = self.selectRowPricing(basis.residual_work);
        const pricing_started = self.statisticsTimestamp();
        defer if (row_pricing)
            self.recordRowPricingElapsed(pricing_started)
        else
            self.recordPricingElapsed(pricing_started);
        if (row_pricing) {
            @memset(basis.tableau[0..problem.num_cols], 0.0);
            for (basis.residual_work, 0..) |dual_value, row_index| {
                if (@abs(dual_value) <= self.numerical.zero_tolerance) continue;
                const scaled_dual = dual_value * basis.row_scale[row_index];
                for (self.pricing_row_view.rowColumns(row_index), self.pricing_row_view.rowValues(row_index)) |column_u32, coefficient| {
                    const column: usize = @intCast(column_u32);
                    basis.tableau[column] += scaled_dual * coefficient * basis.column_scale[column];
                }
            }
            for (basis.reduced_cost[0..problem.num_cols], basis.tableau[0..problem.num_cols]) |*reduced, tableau_entry|
                reduced.* -= theta * tableau_entry;
        } else {
            for (0..problem.num_cols) |column| {
                const begin = problem.matrix.col_starts[column];
                const end = problem.matrix.col_starts[column + 1];
                var tableau_entry: f64 = 0.0;
                for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |matrix_row, coefficient| {
                    if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                    const row_index = matrix_row.toUsize();
                    tableau_entry += basis.residual_work[row_index] * basis.row_scale[row_index] * coefficient * basis.column_scale[column];
                }
                basis.tableau[column] = tableau_entry;
                basis.reduced_cost[column] -= theta * tableau_entry;
            }
        }
        const artificial_begin = problem.num_cols + problem.num_rows;
        for (0..problem.num_rows) |logical_row| {
            const multiplier = theta * basis.residual_work[logical_row];
            basis.tableau[problem.num_cols + logical_row] = basis.residual_work[logical_row];
            basis.tableau[artificial_begin + logical_row] = basis.residual_work[logical_row] * basis.artificial_sign[logical_row];
            basis.reduced_cost[problem.num_cols + logical_row] -= multiplier;
            basis.reduced_cost[artificial_begin + logical_row] -= multiplier * basis.artificial_sign[logical_row];
        }
        basis.reduced_cost[entering_col] = 0.0;
        for (basis.reduced_cost) |value| if (!std.math.isFinite(value)) return .numerical_failure;
        if (self.pricing.rule == .devex and self.active_devex_strategy == .framework and
            self.updatePrimalDevexFramework(entering_col, row, actual_pivot) != .optimal)
            return .numerical_failure;
        self.reduced_cost_update_count += 1;
        return .optimal;
    }

    fn recomputeReducedCosts(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        const maximize = problem.objective_sense == .maximize;
        for (basis.basic_index, 0..) |global_col, row| {
            basis.dual[row] = if (global_col < problem.num_cols)
                (if (maximize) -problem.col_cost[global_col] else problem.col_cost[global_col]) * self.objective_scale * basis.column_scale[global_col]
            else
                0.0;
        }
        self.factorization.solveTranspose(basis.dual) catch return .numerical_failure;
        const row_pricing = self.selectRowPricing(basis.dual);
        const pricing_started = self.statisticsTimestamp();
        defer if (row_pricing)
            self.recordRowPricingElapsed(pricing_started)
        else
            self.recordPricingElapsed(pricing_started);
        for (0..problem.num_cols) |col|
            basis.reduced_cost[col] = (if (maximize) -problem.col_cost[col] else problem.col_cost[col]) * self.objective_scale * basis.column_scale[col];
        if (row_pricing) {
            for (basis.dual, 0..) |dual_value, row_index| {
                if (@abs(dual_value) <= self.numerical.zero_tolerance) continue;
                const scaled_dual = dual_value * basis.row_scale[row_index];
                for (self.pricing_row_view.rowColumns(row_index), self.pricing_row_view.rowValues(row_index)) |column_u32, coefficient| {
                    const column: usize = @intCast(column_u32);
                    basis.reduced_cost[column] -= scaled_dual * coefficient * basis.column_scale[column];
                }
            }
        } else {
            for (0..problem.num_cols) |col| {
                const begin = problem.matrix.col_starts[col];
                const end = problem.matrix.col_starts[col + 1];
                for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, value| {
                    if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(value)) continue;
                    const row_index = row.toUsize();
                    basis.reduced_cost[col] -= basis.row_scale[row_index] * value * basis.column_scale[col] * basis.dual[row_index];
                }
            }
        }
        for (0..problem.num_rows) |row| basis.reduced_cost[problem.num_cols + row] = -basis.dual[row];
        @memset(basis.reduced_cost[problem.num_cols + problem.num_rows ..], 0.0);
        self.reduced_cost_update_count = 0;
        return .optimal;
    }

    fn refactorizeBasis(self: *SimplexEngine, problem: problem_module.ProblemView, reason: RebuildReason) SolveStatus {
        const rebuild_started = self.statisticsTimestamp();
        defer {
            self.stats.rebuild_calls += 1;
            self.stats.rebuild_ns = std.math.add(u64, self.stats.rebuild_ns, self.elapsedSince(rebuild_started)) catch std.math.maxInt(u64);
        }
        switch (reason) {
            .phase_one_setup => self.stats.rebuild_phase_one_setup += 1,
            .solve_residual => self.stats.rebuild_solve_residual += 1,
            .small_pivot => self.stats.rebuild_small_pivot += 1,
            .update_limit => self.stats.rebuild_update_limit += 1,
            .update_growth => self.stats.rebuild_update_growth += 1,
            .direction_refinement => self.stats.rebuild_direction_refinement += 1,
            .fresh_mode => self.stats.rebuild_fresh_mode += 1,
            .update_rejected => self.stats.rebuild_update_rejected += 1,
            .cleanup => self.stats.rebuild_cleanup += 1,
            .edge_weight_reset => self.stats.rebuild_edge_weight_reset += 1,
            .numerical_policy => self.stats.rebuild_numerical_policy += 1,
        }
        self.factorizeCurrentBasis(problem) catch {
            self.failure_site = .pivot_factorization;
            return .numerical_failure;
        };
        return self.finishRefactorization();
    }

    fn factorizeCurrentBasis(self: *SimplexEngine, problem: problem_module.ProblemView) factorization_module.FactorizationError!void {
        const basis = if (self.basis) |*value| value else return error.DimensionMismatch;
        return self.factorization.factorizeBasis(
            problem.matrix,
            basis.basic_index,
            basis.row_scale,
            basis.column_scale[0..problem.num_cols],
            basis.artificial_sign,
        );
    }

    fn finishRefactorization(self: *SimplexEngine) SolveStatus {
        self.numerical.markRefactorized();
        self.observeFactorizationStability();
        self.dual_edge_weights_valid = false;
        self.dual_row_index = null;
        if (self.pricing.rule == .steepest_edge and self.ensureExactDualEdgeWeights() != .optimal) {
            self.failure_site = .pivot_edge_weights;
            return .numerical_failure;
        }
        return .optimal;
    }

    /// Repair a singular imported basis by cumulatively replacing structural
    /// basics with currently nonbasic logical columns. Unique logical columns
    /// form an identity basis, so this deterministic process must recover a
    /// nonsingular basis unless factorization fails for a non-rank reason.
    fn repairRankDeficientBasis(self: *SimplexEngine, problem: problem_module.ProblemView) BasisImportError!void {
        const basis = if (self.basis) |*value| value else return error.NumericalFailure;
        const logical_begin = problem.num_cols;
        const logical_end = logical_begin + problem.num_rows;
        const maximum_incremental_trials = 8;
        var incremental_trials: usize = 0;
        for (0..problem.num_rows) |leaving_row| {
            const leaving_column: usize = basis.basic_index[leaving_row];
            if (leaving_column >= problem.num_cols) continue;

            var entering_column = logical_begin;
            while (entering_column < logical_end and basis.col_status[entering_column] == .basic) : (entering_column += 1) {}
            if (entering_column == logical_end) return error.SingularBasis;

            const leaving_status = nonbasicStatusForBounds(basis.col_lower[leaving_column], basis.col_upper[leaving_column]);
            basis.primal[leaving_column] = switch (leaving_status) {
                .at_upper => basis.col_upper[leaving_column],
                .at_lower, .fixed => basis.col_lower[leaving_column],
                .free, .superbasic => 0.0,
                .basic => unreachable,
            };
            basis.primal[entering_column] = 0.0;
            basis.applyPivot(leaving_row, entering_column, leaving_status) catch return error.NumericalFailure;
            basis.basic_lower[leaving_row] = basis.col_lower[entering_column];
            basis.basic_upper[leaving_row] = basis.col_upper[entering_column];
            self.rank_repair_count += 1;

            // Small deficiencies usually recover after one replacement. Cap
            // repeated INVERT trials; severe deficiency then falls through to
            // one final logical-basis factorization instead of O(n) retries.
            if (incremental_trials < maximum_incremental_trials) {
                incremental_trials += 1;
                self.factorizeCurrentBasis(problem) catch |err| switch (err) {
                    error.Singular => continue,
                    error.OutOfMemory => return error.OutOfMemory,
                    error.DimensionMismatch, error.NotImplemented, error.NumericalFailure => return error.NumericalFailure,
                };
                return;
            }
        }
        self.factorizeCurrentBasis(problem) catch |err| return switch (err) {
            error.Singular => error.SingularBasis,
            error.OutOfMemory => error.OutOfMemory,
            error.DimensionMismatch, error.NotImplemented, error.NumericalFailure => error.NumericalFailure,
        };
    }

    fn nonbasicStatusForBounds(lower: f64, upper: f64) basis_module.BasisStatus {
        if (std.math.isFinite(lower) and std.math.isFinite(upper) and lower == upper) return .fixed;
        if (std.math.isFinite(lower)) return .at_lower;
        if (std.math.isFinite(upper)) return .at_upper;
        return .free;
    }

    fn observeFactorizationStability(self: *SimplexEngine) void {
        self.numerical.pivot_condition_estimate = self.factorization.pivotConditionEstimate();
        if (!std.math.isFinite(self.numerical.pivot_condition_estimate) or self.numerical.pivot_condition_estimate > 1e12)
            self.numerical.numerical_warning = true;
    }

    fn fillInternalColumn(self: *SimplexEngine, problem: problem_module.ProblemView, column: usize, output: []f64) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        @memset(output, 0.0);
        if (column < problem.num_cols) {
            const column_scale = basis.column_scale[column];
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
                if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                const row_index = row.toUsize();
                output[row_index] = coefficient * basis.row_scale[row_index] * column_scale;
            }
        } else if (column < problem.num_cols + problem.num_rows) {
            output[column - problem.num_cols] = 1.0;
        } else if (column < problem.num_cols + 2 * problem.num_rows) {
            const row = column - problem.num_cols - problem.num_rows;
            output[row] = basis.artificial_sign[row];
        } else return .numerical_failure;
        return .optimal;
    }

    /// Recompute the primal basic solution from the immutable problem and
    /// current nonbasic values. Uses engine-owned workspace and performs no
    /// allocation; this limits drift from long Eta update chains.
    fn recomputeBasicValues(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
        return self.recomputeBasicValuesImpl(problem, true);
    }

    fn recomputeBasicValuesUnchecked(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
        return self.recomputeBasicValuesImpl(problem, false);
    }

    fn recomputeBasicValuesImpl(self: *SimplexEngine, problem: problem_module.ProblemView, enforce_bounds: bool) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        @memcpy(basis.residual_work, basis.row_rhs);
        for (basis.col_status, basis.primal, 0..) |status, value, column| {
            if (status == .basic or value == 0.0) continue;
            if (column < problem.num_cols) {
                const begin = problem.matrix.col_starts[column];
                const end = problem.matrix.col_starts[column + 1];
                for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
                    if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                    const row_index = row.toUsize();
                    basis.residual_work[row_index] -= basis.row_scale[row_index] * coefficient * basis.column_scale[column] * value;
                }
            } else if (column < problem.num_cols + problem.num_rows) {
                basis.residual_work[column - problem.num_cols] -= value;
            } else {
                const row = column - problem.num_cols - problem.num_rows;
                basis.residual_work[row] -= basis.artificial_sign[row] * value;
            }
        }
        @memcpy(basis.rhs_work, basis.residual_work);
        self.factorization.solve(basis.rhs_work) catch return .numerical_failure;

        var refinement_step: usize = 0;
        while (true) {
            @memcpy(basis.pivot_direction, basis.residual_work);
            self.subtractBasisProduct(problem, basis.rhs_work, basis.pivot_direction) catch return .numerical_failure;
            var residual_max: f64 = 0.0;
            var rhs_max: f64 = 0.0;
            for (basis.pivot_direction, basis.residual_work) |residual, rhs_value| {
                residual_max = @max(residual_max, @abs(residual));
                rhs_max = @max(rhs_max, @abs(rhs_value));
            }
            self.numerical.observeResidual(residual_max, rhs_max);
            if (!std.math.isFinite(residual_max)) return .numerical_failure;
            if (self.numerical.last_relative_residual <= self.numerical.residual_tolerance or
                refinement_step >= self.numerical.max_refinement_steps)
                break;
            self.factorization.solve(basis.pivot_direction) catch return .numerical_failure;
            for (basis.rhs_work, basis.pivot_direction) |*value, correction| value.* += correction;
            refinement_step += 1;
            self.numerical.refinement_count += 1;
        }
        for (basis.basic_index, basis.basic_value, basis.rhs_work, 0..) |column, *value, recomputed, row| {
            value.* = recomputed;
            basis.primal[column] = recomputed;
            if (enforce_bounds and (recomputed < basis.basic_lower[row] - self.numerical.primal_tolerance or
                recomputed > basis.basic_upper[row] + self.numerical.primal_tolerance))
                return .numerical_failure;
        }
        return .optimal;
    }

    /// Compute `residual -= B * x` directly from the borrowed CSC model and
    /// basis membership. No dense basis copy or temporary column is required.
    fn subtractBasisProduct(self: *SimplexEngine, problem: problem_module.ProblemView, x: []const f64, residual: []f64) !void {
        const basis = if (self.basis) |*value| value else return error.NumericalFailure;
        if (x.len != problem.num_rows or residual.len != problem.num_rows) return error.NumericalFailure;
        for (basis.basic_index, x) |global_col, value| {
            if (value == 0.0) continue;
            const column: usize = @intCast(global_col);
            if (column < problem.num_cols) {
                const begin = problem.matrix.col_starts[column];
                const end = problem.matrix.col_starts[column + 1];
                for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
                    if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                    const row_index = row.toUsize();
                    if (row_index >= residual.len) return error.NumericalFailure;
                    residual[row_index] -= basis.row_scale[row_index] * coefficient * basis.column_scale[column] * value;
                }
            } else if (column < problem.num_cols + problem.num_rows) {
                residual[column - problem.num_cols] -= value;
            } else if (column < problem.num_cols + 2 * problem.num_rows) {
                const row = column - problem.num_cols - problem.num_rows;
                residual[row] -= basis.artificial_sign[row] * value;
            } else return error.NumericalFailure;
        }
    }

    /// Subtract `B*x` and accumulate `|B|*|x|` in the same sparse traversal.
    /// The magnitude output lets FTRAN use a scale-aware backward-error test
    /// without constructing a dense basis or making another matrix pass.
    fn subtractBasisProductWithMagnitude(
        self: *SimplexEngine,
        problem: problem_module.ProblemView,
        x: []const f64,
        residual: []f64,
        magnitude: []f64,
    ) !void {
        const basis = if (self.basis) |*value| value else return error.NumericalFailure;
        if (x.len != problem.num_rows or residual.len != problem.num_rows or magnitude.len != problem.num_rows)
            return error.NumericalFailure;
        for (basis.basic_index, x) |global_col, value| {
            if (!std.math.isFinite(value)) return error.NumericalFailure;
            const column: usize = @intCast(global_col);
            if (column < problem.num_cols) {
                const start = problem.matrix.col_starts[column];
                const end = problem.matrix.col_starts[column + 1];
                for (problem.matrix.row_indices[start..end], problem.matrix.values[start..end]) |row, coefficient| {
                    if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                    const row_index = row.toUsize();
                    if (row_index >= residual.len) return error.NumericalFailure;
                    const product = basis.row_scale[row_index] * coefficient * basis.column_scale[column] * value;
                    residual[row_index] -= product;
                    magnitude[row_index] += @abs(product);
                }
            } else if (column < problem.num_cols + problem.num_rows) {
                const row = column - problem.num_cols;
                residual[row] -= value;
                magnitude[row] += @abs(value);
            } else {
                const row = column - problem.num_cols - problem.num_rows;
                if (row >= residual.len) return error.NumericalFailure;
                const product = basis.artificial_sign[row] * value;
                residual[row] -= product;
                magnitude[row] += @abs(product);
            }
        }
    }

    fn finishOptimal(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
        if (self.recomputeBasicValues(problem) != .optimal) return .numerical_failure;
        if (self.recomputeReducedCosts(problem) != .optimal) return .numerical_failure;
        if (self.validateOptimalSolution(problem) != .optimal) return .numerical_failure;
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        self.objective_value = problem.objective_offset;
        for (problem.col_cost, basis.primal[0..problem.num_cols], basis.column_scale[0..problem.num_cols]) |cost, value, scale|
            self.objective_value += cost * value * scale;
        if (!std.math.isFinite(self.objective_value)) return .numerical_failure;
        return .optimal;
    }

    /// Construct and validate an original-coordinate primal ray from the
    /// signed entering direction and current FTRAN column. No finite iterate
    /// is published as evidence: unboundedness is accepted only when variable
    /// bounds, row recession directions, and objective improvement all hold.
    fn finishUnbounded(self: *SimplexEngine, problem: problem_module.ProblemView, entering_col: usize, entering_direction: f64) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        @memset(basis.unbounded_ray, 0.0);
        if (entering_col < problem.num_cols)
            basis.unbounded_ray[entering_col] = entering_direction * basis.column_scale[entering_col];
        for (basis.basic_index, basis.pivot_direction) |column_u32, direction| {
            const column: usize = @intCast(column_u32);
            if (column < problem.num_cols)
                basis.unbounded_ray[column] = -direction * basis.column_scale[column];
        }
        if (!self.validateUnboundedRay(problem, basis.unbounded_ray)) return .numerical_failure;
        self.unbounded_ray_valid = true;
        return .unbounded;
    }

    fn validateUnboundedRay(self: *SimplexEngine, problem: problem_module.ProblemView, ray: []const f64) bool {
        const basis = if (self.basis) |*value| value else return false;
        if (ray.len != problem.num_cols) return false;
        var ray_max: f64 = 0.0;
        for (ray, problem.col_lower, problem.col_upper) |direction, lower, upper| {
            if (!std.math.isFinite(direction)) return false;
            ray_max = @max(ray_max, @abs(direction));
            if (direction > self.numerical.zero_tolerance and std.math.isFinite(upper)) return false;
            if (direction < -self.numerical.zero_tolerance and std.math.isFinite(lower)) return false;
        }
        if (ray_max <= self.numerical.zero_tolerance) return false;

        @memset(basis.rhs_work, 0.0);
        @memset(basis.residual_work, 0.0);
        for (0..problem.num_cols) |column| {
            const direction = ray[column];
            if (direction == 0.0) continue;
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
                if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                const row_index = row.toUsize();
                const product = coefficient * direction;
                basis.rhs_work[row_index] += product;
                basis.residual_work[row_index] += @abs(product);
            }
        }
        for (basis.rhs_work, basis.residual_work, problem.row_lower, problem.row_upper) |direction, magnitude, lower, upper| {
            const tolerance = self.numerical.residual_tolerance * @max(1.0, magnitude);
            if (std.math.isFinite(upper) and direction > tolerance) return false;
            if (std.math.isFinite(lower) and direction < -tolerance) return false;
        }

        var objective_direction: f64 = 0.0;
        var objective_magnitude: f64 = 0.0;
        for (problem.col_cost, ray) |cost, direction| {
            objective_direction += cost * direction;
            objective_magnitude += @abs(cost * direction);
        }
        const objective_tolerance = self.numerical.dual_tolerance * @max(1.0, objective_magnitude);
        return if (problem.objective_sense == .minimize)
            objective_direction < -objective_tolerance
        else
            objective_direction > objective_tolerance;
    }

    /// Allocation-free KKT feasibility check used before publishing an
    /// optimal result. This catches accumulated update drift and incorrect
    /// bound/status transitions at the solver boundary.
    fn validateOptimalSolution(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        const primal_tolerance = self.numerical.primal_tolerance;
        const dual_tolerance = self.numerical.dual_tolerance;

        for (basis.primal[0..problem.num_cols], basis.column_scale[0..problem.num_cols], problem.col_lower, problem.col_upper) |internal, scale, lower, upper| {
            const value = internal * scale;
            if (!std.math.isFinite(value) or value < lower - primal_tolerance or value > upper + primal_tolerance)
                return .numerical_failure;
        }

        @memset(basis.rhs_work, 0.0);
        for (0..problem.num_cols) |column| {
            const value = basis.primal[column] * basis.column_scale[column];
            if (value == 0.0) continue;
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
                if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                basis.rhs_work[row.toUsize()] += coefficient * value;
            }
        }
        for (basis.rhs_work, problem.row_lower, problem.row_upper) |activity, lower, upper| {
            if (!std.math.isFinite(activity) or activity < lower - primal_tolerance or activity > upper + primal_tolerance)
                return .numerical_failure;
        }

        for (basis.reduced_cost[0 .. problem.num_cols + problem.num_rows], basis.col_status[0 .. problem.num_cols + problem.num_rows], 0..) |internal_reduced, status, column| {
            const reduced = if (column < problem.num_cols)
                internal_reduced / (basis.column_scale[column] * self.objective_scale)
            else
                internal_reduced / self.objective_scale;
            if (!std.math.isFinite(reduced)) return .numerical_failure;
            const infeasible = switch (status) {
                .at_lower => reduced < -dual_tolerance,
                .at_upper => reduced > dual_tolerance,
                .free, .superbasic, .basic => @abs(reduced) > dual_tolerance,
                .fixed => false,
            };
            if (infeasible) return .numerical_failure;
        }

        const artificial_begin = problem.num_cols + problem.num_rows;
        for (basis.primal[artificial_begin..]) |value| {
            if (!std.math.isFinite(value) or @abs(value) > primal_tolerance) return .numerical_failure;
        }
        return .optimal;
    }

    /// Borrow the current engine-owned solution arrays without copying.
    pub fn solutionView(self: *const SimplexEngine, problem: problem_module.ProblemView, status: SolveStatus) ?solution_module.SolutionView {
        const basis = if (self.basis) |*value| value else return null;
        for (basis.published_primal, basis.primal[0..problem.num_cols], basis.column_scale[0..problem.num_cols]) |*published, internal, scale|
            published.* = internal * scale;
        for (basis.published_dual, basis.dual) |*published, internal|
            published.* = internal / self.objective_scale;
        for (basis.published_reduced_cost, basis.reduced_cost[0..problem.num_cols], basis.column_scale[0..problem.num_cols]) |*published, internal, scale|
            published.* = internal / (scale * self.objective_scale);
        return .{
            .status = status,
            .primal = basis.published_primal,
            .dual = basis.published_dual,
            .reduced_cost = basis.published_reduced_cost,
            .unbounded_ray = if (self.unbounded_ray_valid) basis.unbounded_ray else &.{},
            .objective_value = self.objective_value,
            .iterations = self.iterations,
        };
    }
};

test {
    std.testing.refAllDecls(@This());
}

test "anti cycling uses basic-column order for degenerate primal and dual ties" {
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 2, 4);
    engine.basis.?.initializeSlackBasis();
    engine.basis.?.basic_index[0] = 5;
    engine.basis.?.basic_index[1] = 4;
    @memset(engine.basis.?.basic_value, 0.0);
    @memset(engine.basis.?.basic_lower, 0.0);
    @memset(engine.basis.?.basic_upper, std.math.inf(f64));
    @memset(engine.basis.?.pivot_direction, 1.0);
    engine.numerical.anti_cycling_active = true;
    engine.numerical.perturbation = 1e-8;

    const primal = engine.chooseLeaving();
    try std.testing.expectEqual(SolveStatus.optimal, primal.status);
    try std.testing.expectEqual(@as(?u32, 1), primal.row);

    // Bland ordering must not admit a numerically unsafe pivot merely because
    // its basic-column ID is lower than every stable candidate.
    engine.basis.?.basic_index[0] = 3;
    engine.basis.?.pivot_direction[0] = engine.ratio_test.tolerance * 50.0;
    engine.basis.?.pivot_direction[1] = 100.0;
    const stable_primal = engine.chooseLeaving();
    try std.testing.expectEqual(@as(?u32, 1), stable_primal.row);

    engine.basis.?.pivot_direction[0] = 2.0;
    engine.basis.?.pivot_direction[1] = 100.0;
    try std.testing.expectEqual(@as(?u32, 0), engine.chooseLeaving().row);
    try std.testing.expectEqual(@as(?u32, 1), engine.chooseLeavingWithPolicy(false).row);

    engine.basis.?.basic_value[0] = -1.0;
    engine.basis.?.basic_value[1] = -1.0;
    const dual = engine.chooseDualLeavingRow().?;
    try std.testing.expectEqual(@as(u32, 1), dual.row);
}

test "engine solves a standard-form LP with primal revised simplex" {
    const rows = [_]foundation.RowId{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) };
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &[_]f64{ -1, -1 },
        .col_lower = &[_]f64{ 0, 0 },
        .col_upper = &[_]f64{ std.math.inf(f64), std.math.inf(f64) },
        .row_lower = &[_]f64{-std.math.inf(f64)},
        .row_upper = &[_]f64{4},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &[_]usize{ 0, 1, 2 }, &rows, &[_]f64{ 1, 2 }),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    const status = engine.solveProblem(problem, .{});
    try std.testing.expectEqual(SolveStatus.optimal, status);
    try std.testing.expectApproxEqAbs(@as(f64, 4), engine.basis.?.primal[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0), engine.basis.?.primal[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, -4), engine.objective_value, 1e-12);
    const solution = engine.solutionView(problem, status).?;
    try std.testing.expectEqual(@as(usize, 2), solution.primal.len);
    try std.testing.expectApproxEqAbs(@as(f64, -4), solution.objective_value, 1e-12);
    try std.testing.expectEqual(engine.iterations, engine.stats.phase_two_iterations);
    try std.testing.expectEqual(@as(usize, 0), engine.stats.phase_one_iterations);
    try std.testing.expect(engine.stats.pricing_calls > 0);
    try std.testing.expectEqual(@as(u64, 0), engine.stats.pricing_ns);
    try std.testing.expectEqual(engine.factorization.stats.ftran_calls, engine.factorization.stats.dense_ftran_dispatches);
    try std.testing.expectEqual(engine.stats.rebuild_calls, engine.stats.classifiedRebuilds());
    try std.testing.expect(engine.requestedBytes() > 0);
}

test "Bland fallback resolves the Beale degenerate cycling example" {
    const rows = [_]foundation.RowId{
        foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(2),
        foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(0),
        foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1),
    };
    const problem = problem_module.ProblemView{
        .num_rows = 3,
        .num_cols = 4,
        .col_cost = &[_]f64{ 10, -57, -9, -24 },
        .col_lower = &[_]f64{ 0, 0, 0, 0 },
        .col_upper = &[_]f64{ std.math.inf(f64), std.math.inf(f64), std.math.inf(f64), std.math.inf(f64) },
        .row_lower = &[_]f64{ -std.math.inf(f64), -std.math.inf(f64), -std.math.inf(f64) },
        .row_upper = &[_]f64{ 0, 0, 1 },
        .matrix = matrix.CscView.initAssumeValid(
            3,
            4,
            &[_]usize{ 0, 3, 5, 7, 9 },
            &rows,
            &[_]f64{ 0.5, 0.5, 1, -5.5, -1.5, -2.5, -0.5, 9, 1 },
        ),
        .objective_sense = .maximize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.pricing.rule = .dantzig;
    engine.numerical.degenerate_pivot_limit = 1;
    var degeneracy_trace: [100]DegeneracyTraceEvent = undefined;

    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{
        .max_iterations = 100,
        .degeneracy_trace = &degeneracy_trace,
    }));
    try std.testing.expect(engine.numerical.anti_cycling_activations > 0);
    try std.testing.expect(engine.numerical.degenerate_pivot_count > 0);
    try std.testing.expectEqual(engine.numerical.degenerate_pivot_count, engine.stats.classifiedDegeneratePivots());
    try std.testing.expectEqual(engine.numerical.degenerate_pivot_count, engine.degeneracy_trace_count);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), engine.objective_value, 1e-9);
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
    const status = engine.solveProblem(problem, .{});
    try std.testing.expectEqual(SolveStatus.optimal, status);
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

test "dual bound flips share one aggregate FTRAN" {
    const rows = [_]foundation.RowId{
        foundation.RowId.fromUsizeAssumeValid(0),
        foundation.RowId.fromUsizeAssumeValid(1),
        foundation.RowId.fromUsizeAssumeValid(0),
        foundation.RowId.fromUsizeAssumeValid(1),
    };
    const problem = problem_module.ProblemView{
        .num_rows = 2,
        .num_cols = 2,
        .col_cost = &[_]f64{ 0, 0 },
        .col_lower = &[_]f64{ 0, 0 },
        .col_upper = &[_]f64{ 2, 3 },
        .row_lower = &[_]f64{ -std.math.inf(f64), -std.math.inf(f64) },
        .row_upper = &[_]f64{ 10, 10 },
        .matrix = matrix.CscView.initAssumeValid(
            2,
            2,
            &[_]usize{ 0, 2, 4 },
            &rows,
            &[_]f64{ 1, 2, 3, 4 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 2, 2);
    engine.basis.?.initializeSlackBasis();
    engine.basis.?.col_lower[0] = 0;
    engine.basis.?.col_lower[1] = 0;
    engine.basis.?.col_upper[0] = 2;
    engine.basis.?.col_upper[1] = 3;
    engine.basis.?.basic_value[0] = 10;
    engine.basis.?.basic_value[1] = 10;
    engine.basis.?.primal[2] = 10;
    engine.basis.?.primal[3] = 10;
    engine.basis.?.flip_columns[0] = 0;
    engine.basis.?.flip_columns[1] = 1;
    try engine.factorization.factorizeIdentity(2);

    try std.testing.expectEqual(SolveStatus.optimal, engine.applyBoundFlips(problem, 2));
    // Combined equation displacement is [1*2 + 3*3, 2*2 + 4*3] = [11, 16].
    try std.testing.expectApproxEqAbs(@as(f64, -1), engine.basis.?.basic_value[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, -6), engine.basis.?.basic_value[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2), engine.basis.?.primal[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 3), engine.basis.?.primal[1], 1e-12);
    try std.testing.expectEqual(basis_module.BasisStatus.at_upper, engine.basis.?.col_status[0]);
    try std.testing.expectEqual(basis_module.BasisStatus.at_upper, engine.basis.?.col_status[1]);
    try std.testing.expectEqual(@as(usize, 1), engine.factorization.stats.ftran_calls);
    try std.testing.expectEqual(@as(usize, 1), engine.stats.bound_flip_batches);
    try std.testing.expectEqual(@as(usize, 1), engine.stats.bound_flip_ftran_savings);
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

test "Phase I restores feasibility for a greater-equal row" {
    const rows = [_]foundation.RowId{foundation.RowId.fromUsizeAssumeValid(0)};
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{1},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{std.math.inf(f64)},
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
}

test "Phase I removes zero artificial basics after redundant equalities" {
    const rows = [_]foundation.RowId{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1) };
    const problem = problem_module.ProblemView{
        .num_rows = 2,
        .num_cols = 1,
        .col_cost = &[_]f64{0},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{ 1, 1 },
        .row_upper = &[_]f64{ 1, 1 },
        .matrix = matrix.CscView.initAssumeValid(2, 1, &[_]usize{ 0, 2 }, &rows, &[_]f64{ 1, 1 }),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{}));
    try std.testing.expectApproxEqAbs(@as(f64, 1), engine.basis.?.primal[0], 1e-12);
    const artificial_begin = problem.num_cols + problem.num_rows;
    for (engine.basis.?.basic_index) |column| try std.testing.expect(column < artificial_begin);
}

test "Phase I detects an infeasible LP" {
    const rows = [_]foundation.RowId{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1) };
    const problem = problem_module.ProblemView{
        .num_rows = 2,
        .num_cols = 1,
        .col_cost = &[_]f64{0},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{ 2, -std.math.inf(f64) },
        .row_upper = &[_]f64{ std.math.inf(f64), 1 },
        .matrix = matrix.CscView.initAssumeValid(2, 1, &[_]usize{ 0, 2 }, &rows, &[_]f64{ 1, 1 }),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.infeasible, engine.solveProblem(problem, .{}));
}

test "engine detects an unbounded improving structural column" {
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{-1},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{-std.math.inf(f64)},
        .row_upper = &[_]f64{1},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &[_]usize{ 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    const status = engine.solveProblem(problem, .{});
    try std.testing.expectEqual(SolveStatus.unbounded, status);
    const solution = engine.solutionView(problem, status).?;
    try std.testing.expectEqual(@as(usize, 1), solution.unbounded_ray.len);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), solution.unbounded_ray[0], 1e-12);
}

test "Phase I respects a zero iteration limit" {
    const rows = [_]foundation.RowId{foundation.RowId.fromUsizeAssumeValid(0)};
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{0},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{2},
        .row_upper = &[_]f64{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &[_]usize{ 0, 1 }, &rows, &[_]f64{1}),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.iteration_limit, engine.solveProblem(problem, .{ .max_iterations = 0 }));
}

test "deterministic work limit spans Phase I and Phase II" {
    const rows = [_]foundation.RowId{foundation.RowId.fromUsizeAssumeValid(0)};
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{1},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{2},
        .row_upper = &[_]f64{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &[_]usize{ 0, 1 }, &rows, &[_]f64{1}),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.work_limit, engine.solveProblem(problem, .{ .work_limit = 1 }));
    try std.testing.expectEqual(@as(u64, 1), engine.work_used);
}

test "iteration callback can stop without allocating callback state" {
    const Context = struct {
        calls: usize = 0,
        last_work: u64 = 0,

        fn callback(event: ProgressEventView, context_ptr: ?*anyopaque) CallbackAction {
            const self: *@This() = @ptrCast(@alignCast(context_ptr.?));
            self.calls += 1;
            self.last_work = event.work_used;
            return .stop;
        }
    };
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{-1},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{1},
        .row_lower = &[_]f64{-std.math.inf(f64)},
        .row_upper = &[_]f64{1},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &[_]usize{ 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var context = Context{};
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.interrupted, engine.solveProblem(problem, .{
        .iteration_callback = Context.callback,
        .callback_user_data = &context,
    }));
    try std.testing.expectEqual(@as(usize, 1), context.calls);
    try std.testing.expectEqual(@as(u64, 0), context.last_work);
    try std.testing.expectEqual(@as(u64, 0), engine.work_used);
}

test "structured iteration logging obeys its deterministic interval" {
    const Context = struct {
        calls: usize = 0,
        fn log(_: ProgressEventView, context_ptr: ?*anyopaque) void {
            const self: *@This() = @ptrCast(@alignCast(context_ptr.?));
            self.calls += 1;
        }
    };
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{-1},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{1},
        .row_lower = &[_]f64{-std.math.inf(f64)},
        .row_upper = &[_]f64{1},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &[_]usize{ 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var context = Context{};
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{
        .log_level = .iterations,
        .log_callback = Context.log,
        .log_user_data = &context,
        .log_interval_work = 1,
    }));
    try std.testing.expect(context.calls >= 1);
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

test "engine honors an immediate time limit" {
    const problem = problem_module.ProblemView{
        .num_rows = 0,
        .num_cols = 1,
        .col_cost = &[_]f64{0},
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
    try std.testing.expectEqual(SolveStatus.time_limit, engine.solveProblem(problem, .{ .time_limit_ns = 0 }));
}

test "engine honors a caller-owned atomic interrupt flag" {
    const problem = problem_module.ProblemView{
        .num_rows = 0,
        .num_cols = 1,
        .col_cost = &[_]f64{0},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{1},
        .row_lower = &.{},
        .row_upper = &.{},
        .matrix = matrix.CscView.initAssumeValid(0, 1, &[_]usize{ 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var interrupted = std.atomic.Value(bool).init(true);
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.interrupted, engine.solveProblem(problem, .{ .interrupt_flag = &interrupted }));
}

test "imported dual-feasible basis reoptimizes a changed RHS" {
    const rows = [_]foundation.RowId{foundation.RowId.fromUsizeAssumeValid(0)};
    const original = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{1},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{2},
        .row_upper = &[_]f64{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &[_]usize{ 0, 1 }, &rows, &[_]f64{1}),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var first_engine = SimplexEngine.init(std.testing.allocator);
    defer first_engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, first_engine.solveProblem(original, .{}));
    const exported = first_engine.exportBasisView(original).?;
    var snapshot = try basis_snapshot_module.BasisSnapshot.initFromView(std.testing.allocator, exported);
    defer snapshot.deinit();

    const modified = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = original.col_cost,
        .col_lower = original.col_lower,
        .col_upper = original.col_upper,
        .row_lower = &[_]f64{-1},
        .row_upper = original.row_upper,
        .matrix = original.matrix,
        .objective_sense = original.objective_sense,
        .objective_offset = original.objective_offset,
    };
    var second_engine = SimplexEngine.init(std.testing.allocator);
    defer second_engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, second_engine.solveProblem(modified, .{ .initial_basis = snapshot.view() }));
    try std.testing.expectEqual(Algorithm.dual_revised, second_engine.algorithm);
    try std.testing.expectApproxEqAbs(@as(f64, 0), second_engine.basis.?.primal[0], 1e-12);
    try std.testing.expect(second_engine.stats.dual_reduced_cost_updates > 0);

    // Exact refresh replaces a deliberately drifted incremental value and
    // records the normalized discrepancy before terminal publication.
    second_engine.basis.?.reduced_cost[0] += 1e-4;
    try std.testing.expectEqual(SolveStatus.optimal, second_engine.recomputeReducedCostsWithDrift(modified));
    try std.testing.expect(second_engine.maximum_reduced_cost_drift >= 1e-4);
    try std.testing.expectEqual(@as(usize, 1), second_engine.stats.dual_exact_reprices);

    // The explicit steepest-devex lifecycle starts with exact DSE and uses a
    // one-update budget here to force a deterministic dual Devex fallback.
    var fallback_engine = SimplexEngine.init(std.testing.allocator);
    defer fallback_engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, fallback_engine.solveProblem(modified, .{
        .initial_basis = snapshot.view(),
        .dual_edge_weight_strategy = .steepest_devex,
        .dual_dse_update_budget = 1,
    }));
    try std.testing.expectEqual(Algorithm.dual_revised, fallback_engine.algorithm);
    try std.testing.expect(fallback_engine.stats.dual_dse_updates > 0);
    try std.testing.expectEqual(@as(usize, 1), fallback_engine.stats.dual_dse_budget_fallbacks);
}

test "singular imported basis is repaired without discarding independent structural basics" {
    const rows = [_]foundation.RowId{
        foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1),
        foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1),
    };
    const problem = problem_module.ProblemView{
        .num_rows = 2,
        .num_cols = 2,
        .col_cost = &[_]f64{ 0, 0 },
        .col_lower = &[_]f64{ 0, 0 },
        .col_upper = &[_]f64{ std.math.inf(f64), std.math.inf(f64) },
        .row_lower = &[_]f64{ 1, 1 },
        .row_upper = &[_]f64{ 1, 1 },
        .matrix = matrix.CscView.initAssumeValid(
            2,
            2,
            &[_]usize{ 0, 2, 4 },
            &rows,
            &[_]f64{ 1, 1, 1, 1 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    const initial_basis = basis_snapshot_module.BasisView{
        .structural_status = &[_]basis_module.BasisStatus{ .basic, .basic },
        .logical_status = &[_]basis_module.BasisStatus{ .fixed, .fixed },
        .basic_index = &[_]u32{ 0, 1 },
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();

    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{ .initial_basis = initial_basis }));
    try std.testing.expectEqual(@as(usize, 1), engine.rank_repair_count);
    var structural_basics: usize = 0;
    for (engine.basis.?.basic_index) |column| structural_basics += @intFromBool(column < problem.num_cols);
    try std.testing.expectEqual(@as(usize, 1), structural_basics);
    try std.testing.expect(engine.factorization.pivotConditionEstimate() < 1e12);
}

test "rank repair cumulatively replaces multiple dependent structural basics" {
    const rows = [_]foundation.RowId{
        foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(2),
        foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(2),
        foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(2),
    };
    const problem = problem_module.ProblemView{
        .num_rows = 3,
        .num_cols = 3,
        .col_cost = &[_]f64{ 0, 0, 0 },
        .col_lower = &[_]f64{ 0, 0, 0 },
        .col_upper = &[_]f64{ std.math.inf(f64), std.math.inf(f64), std.math.inf(f64) },
        .row_lower = &[_]f64{ 1, 1, 1 },
        .row_upper = &[_]f64{ 1, 1, 1 },
        .matrix = matrix.CscView.initAssumeValid(
            3,
            3,
            &[_]usize{ 0, 3, 6, 9 },
            &rows,
            &[_]f64{ 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    const initial_basis = basis_snapshot_module.BasisView{
        .structural_status = &[_]basis_module.BasisStatus{ .basic, .basic, .basic },
        .logical_status = &[_]basis_module.BasisStatus{ .fixed, .fixed, .fixed },
        .basic_index = &[_]u32{ 0, 1, 2 },
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();

    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{ .initial_basis = initial_basis }));
    try std.testing.expectEqual(@as(usize, 2), engine.rank_repair_count);
    var structural_basics: usize = 0;
    for (engine.basis.?.basic_index) |column| structural_basics += @intFromBool(column < problem.num_cols);
    try std.testing.expectEqual(@as(usize, 1), structural_basics);
}

test "rank repair restores a singular sparse backend basis" {
    const n = 64;
    const allocator = std.testing.allocator;
    const starts = try allocator.alloc(usize, n + 1);
    defer allocator.free(starts);
    const rows = try allocator.alloc(foundation.RowId, n);
    defer allocator.free(rows);
    const values = try allocator.alloc(f64, n);
    defer allocator.free(values);
    const zeros = try allocator.alloc(f64, n);
    defer allocator.free(zeros);
    const infinities = try allocator.alloc(f64, n);
    defer allocator.free(infinities);
    const structural_status = try allocator.alloc(basis_module.BasisStatus, n);
    defer allocator.free(structural_status);
    const logical_status = try allocator.alloc(basis_module.BasisStatus, n);
    defer allocator.free(logical_status);
    const basic_index = try allocator.alloc(u32, n);
    defer allocator.free(basic_index);
    for (0..n) |column| {
        starts[column] = column;
        rows[column] = foundation.RowId.fromUsizeAssumeValid(if (column == 1) 0 else column);
        values[column] = 1.0;
        basic_index[column] = @intCast(column);
    }
    starts[n] = n;
    @memset(zeros, 0.0);
    @memset(infinities, std.math.inf(f64));
    @memset(structural_status, .basic);
    @memset(logical_status, .fixed);

    const problem = problem_module.ProblemView{
        .num_rows = n,
        .num_cols = n,
        .col_cost = zeros,
        .col_lower = zeros,
        .col_upper = infinities,
        .row_lower = zeros,
        .row_upper = zeros,
        .matrix = matrix.CscView.initAssumeValid(n, n, starts, rows, values),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    const initial_basis = basis_snapshot_module.BasisView{
        .structural_status = structural_status,
        .logical_status = logical_status,
        .basic_index = basic_index,
    };
    var engine = SimplexEngine.init(allocator);
    defer engine.deinit();

    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{ .initial_basis = initial_basis }));
    try std.testing.expectEqual(factorization_module.BackendKind.sparse_lu, engine.factorization.backend_kind);
    try std.testing.expectEqual(@as(usize, 2), engine.rank_repair_count);
}

test "neither-feasible imported basis is repaired without a crash restart" {
    const Context = struct {
        saw_repair: bool = false,
        fn callback(event: ProgressEventView, context_ptr: ?*anyopaque) CallbackAction {
            const self: *@This() = @ptrCast(@alignCast(context_ptr.?));
            self.saw_repair = self.saw_repair or event.phase == .dual_feasibility_repair;
            return .continue_solve;
        }
    };
    const rows = [_]foundation.RowId{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) };
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &[_]f64{ 0, -1 },
        .col_lower = &[_]f64{ 0, 0 },
        .col_upper = &[_]f64{ 1, std.math.inf(f64) },
        .row_lower = &[_]f64{2},
        .row_upper = &[_]f64{2},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &[_]usize{ 0, 1, 2 }, &rows, &[_]f64{ 1, 1 }),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    const initial_basis = basis_snapshot_module.BasisView{
        .structural_status = &[_]basis_module.BasisStatus{ .basic, .at_lower },
        .logical_status = &[_]basis_module.BasisStatus{.fixed},
        .basic_index = &[_]u32{0},
    };
    var context = Context{};
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{
        .initial_basis = initial_basis,
        .iteration_callback = Context.callback,
        .callback_user_data = &context,
    }));
    try std.testing.expect(context.saw_repair);
    try std.testing.expectEqual(Algorithm.primal_revised, engine.algorithm);
    try std.testing.expectApproxEqAbs(@as(f64, 0), engine.basis.?.primal[0], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 2), engine.basis.?.primal[1], 1e-10);
}

test "exact dual steepest-edge weights use BTRAN row norms" {
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 2, 0);
    engine.basis.?.basic_value[0] = -1.0;
    engine.basis.?.basic_value[1] = -2.0;
    @memset(engine.basis.?.basic_lower, 0.0);
    @memset(engine.basis.?.basic_upper, std.math.inf(f64));
    try engine.factorization.factorize(2, &[_]f64{ 2, 0, 0, 0.5 });
    try std.testing.expectEqual(SolveStatus.optimal, engine.ensureExactDualEdgeWeights());
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), engine.basis.?.row_edge_weight[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), engine.basis.?.row_edge_weight[1], 1e-12);

    engine.pricing.rule = .steepest_edge;
    const choice = engine.pricing.chooseDualLeavingWeighted(
        engine.basis.?.basic_value,
        engine.basis.?.basic_lower,
        engine.basis.?.basic_upper,
        engine.basis.?.row_edge_weight,
        engine.numerical.primal_tolerance,
    ).?;
    try std.testing.expectEqual(@as(u32, 0), choice.row);
}

test "incremental dual steepest-edge recurrence matches the new inverse rows" {
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 2, 0);
    try engine.factorization.factorizeIdentity(2);
    engine.pricing.rule = .steepest_edge;
    engine.dual_edge_weights_valid = true;
    engine.dual_row_index = 0;
    engine.basis.?.dual_row[0] = 1.0;
    engine.basis.?.dual_row[1] = 0.0;
    engine.basis.?.pivot_direction[0] = 2.0;
    engine.basis.?.pivot_direction[1] = 1.0;
    try std.testing.expectEqual(SolveStatus.optimal, engine.updateDualSteepestEdgeWeights(0, 2.0));
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), engine.basis.?.row_edge_weight[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1.25), engine.basis.?.row_edge_weight[1], 1e-12);
}

test "dual Devex updates every row from the hot FTRAN column" {
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 2, 0);
    engine.basis.?.row_edge_weight[0] = 1;
    engine.basis.?.row_edge_weight[1] = 2;
    engine.basis.?.pivot_direction[0] = 2;
    engine.basis.?.pivot_direction[1] = 1;

    engine.updateDualDevexWeights(1, 1);

    try std.testing.expectApproxEqAbs(@as(f64, 8), engine.basis.?.row_edge_weight[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2), engine.basis.?.row_edge_weight[1], 1e-12);
    try std.testing.expectEqual(@as(usize, 1), engine.stats.dual_devex_updates);
}

test "rejected dual steepest-edge framework falls back transactionally" {
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 2, 0);
    engine.active_dual_edge_weight_strategy = .steepest_devex;
    engine.pricing.rule = .steepest_edge;
    engine.dual_edge_weights_valid = true;
    engine.dual_row_index = 1;
    engine.dual_candidate_count = 2;
    engine.basis.?.row_edge_weight[0] = 9;
    engine.basis.?.row_edge_weight[1] = 4;

    engine.switchDualDseToDevex(.invalid);

    try std.testing.expectEqual(pricing_module.PricingRule.devex, engine.pricing.rule);
    try std.testing.expectEqualSlices(f64, &[_]f64{ 1, 1 }, engine.basis.?.row_edge_weight);
    try std.testing.expect(!engine.dual_edge_weights_valid);
    try std.testing.expectEqual(@as(?u32, null), engine.dual_row_index);
    try std.testing.expectEqual(@as(usize, 0), engine.dual_candidate_count);
    try std.testing.expectEqual(@as(usize, 1), engine.stats.dual_dse_invalid_fallbacks);
}

test "dual rank one reduced cost update uses the old tableau row" {
    const rows = [_]foundation.RowId{
        foundation.RowId.fromUsizeAssumeValid(0),
        foundation.RowId.fromUsizeAssumeValid(0),
    };
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &[_]f64{ 0, 0 },
        .col_lower = &[_]f64{ 0, 0 },
        .col_upper = &[_]f64{ 1, 1 },
        .row_lower = &[_]f64{0},
        .row_upper = &[_]f64{0},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &[_]usize{ 0, 1, 2 }, &rows, &[_]f64{ 4, 3 }),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 1, 2);
    engine.basis.?.reduced_cost[0] = 2;
    engine.basis.?.reduced_cost[1] = -1;
    engine.basis.?.tableau[0] = 4;
    engine.basis.?.tableau[1] = 3;
    engine.basis.?.tableau[2] = 1;
    engine.basis.?.dual_row[0] = 1;
    engine.basis.?.artificial_sign[0] = -1;

    try std.testing.expectEqual(SolveStatus.optimal, engine.updateReducedCostsAfterDualPivot(problem, 0, 0));
    try std.testing.expectApproxEqAbs(@as(f64, 0), engine.basis.?.reduced_cost[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, -2.5), engine.basis.?.reduced_cost[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, -0.5), engine.basis.?.reduced_cost[2], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), engine.basis.?.reduced_cost[3], 1e-12);
    try std.testing.expectEqual(@as(usize, 1), engine.reduced_cost_update_count);
    try std.testing.expectEqual(@as(usize, 1), engine.stats.dual_reduced_cost_updates);
}

test "primal Devex framework updates every nonbasic tableau weight" {
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 2, 2);
    engine.basis.?.initializeSlackBasis();
    engine.initializePrimalDevexFramework();
    engine.basis.?.pivot_direction[0] = 2;
    engine.basis.?.pivot_direction[1] = 1;
    engine.basis.?.tableau[0] = 2;
    engine.basis.?.tableau[1] = 3;
    engine.basis.?.tableau[2] = 1;

    try std.testing.expectEqual(SolveStatus.optimal, engine.updatePrimalDevexFramework(0, 0, 2));
    // Initial logical basics are outside the frozen reference set, so the
    // pivotal reference norm is the entering coordinate alone: 1 / 2^2.
    try std.testing.expectApproxEqAbs(@as(f64, 3.25), engine.basis.?.col_edge_weight[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1), engine.basis.?.col_edge_weight[2], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1), engine.basis.?.col_edge_weight[0], 1e-12);
    try std.testing.expectEqual(@as(usize, 1), engine.stats.devex_framework_updates);
}

test "hyper-sparse dual candidate list retains the most attractive rows" {
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 64, 0);
    @memset(engine.basis.?.basic_lower, 0.0);
    @memset(engine.basis.?.basic_upper, std.math.inf(f64));
    @memset(engine.basis.?.basic_value, 0.0);
    for (24..64) |row| engine.basis.?.basic_value[row] = -@as(f64, @floatFromInt(row - 23));
    engine.pricing.rule = .hyper_sparse;
    engine.dual_hyper_sparse_active = true;
    const choice = engine.chooseDualLeavingRow().?;
    try std.testing.expectEqual(@as(u32, 63), choice.row);
    try std.testing.expectEqual(@as(usize, 32), engine.dual_candidate_count);
    try std.testing.expect(engine.dual_candidate_cutoff > 0.0);
}
