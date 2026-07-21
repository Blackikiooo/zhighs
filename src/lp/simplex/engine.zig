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
    ratio_test,
    primal_step,
    dual_feasibility,
    reduced_cost,
    optimality_check,
};

pub const ShiftedDualExit = enum {
    none,
    running,
    setup_free_infeasibility,
    setup_not_dual_feasible,
    phase_two_no_entering,
    phase_two_numerical_failure,
    phase_two_stopped,
    original_dual_feasible,
    original_dual_phase_two_optimal,
    original_dual_phase_two_failed,
    cleanup_primal_optimal,
    cleanup_primal_failed,
    cleanup_neither_feasible,
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
pub const PrimalLeavingResult = struct {
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
    /// Explicit dual DSE -> Devex A/B policy. Default is steepest-devex
    /// (HiGHS-equivalent: SteepestEdge with Devex fallback).
    dual_edge_weight_strategy: DualEdgeWeightStrategy = .steepest_devex,
    /// Maximum successful DSE recurrence updates per dual phase before
    /// switching to dual Devex. Zero disables the budget-triggered switch.
    dual_dse_update_budget: usize = 64,
};

pub const SimplexStats = struct {
    phase_one_iterations: usize = 0,
    dual_phase_one_iterations: usize = 0,
    dual_phase_one_fallbacks: usize = 0,
    shifted_dual_iterations: usize = 0,
    shifted_cleanup_iterations: usize = 0,
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
    shifted_dual_exit: ShiftedDualExit = .none,
    shifted_dual_failure_site: FailureSite = .none,
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
        self.shifted_dual_exit = .none;
        self.shifted_dual_failure_site = .none;
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
            const shifted_dual_status = self.solveDualWithCostShifts(problem, control);
            // A shifted-cost CHUZC failure is not an original-model primal
            // infeasibility certificate. Only proven optimality and external
            // stop limits may escape this experimental path; all algorithmic
            // failures continue through Phase I / the frozen primal fallback.
            switch (shifted_dual_status) {
                .optimal, .work_limit, .time_limit, .iteration_limit, .interrupted => return shifted_dual_status,
                .infeasible, .unbounded, .not_implemented, .numerical_failure => {},
            }
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
