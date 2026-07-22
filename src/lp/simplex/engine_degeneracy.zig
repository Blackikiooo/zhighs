//! Degeneracy and anti-cycling policy driver for `SimplexEngine`.
//!
//! ## Responsibility
//!
//! Owns degeneracy-strategy preparation, per-iteration tie observation,
//! degenerate-basis fingerprinting, and transactional restarts without
//! perturbation.

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
const SimplexEngine = @import("engine.zig").SimplexEngine;
const SolveStatus = @import("engine.zig").SolveStatus;
const DegeneracyReason = @import("engine.zig").DegeneracyReason;
const SolveControl = @import("engine.zig").SolveControl;

/// Terminal certificates are always checked in original coordinates. If
/// a perturbed path reaches a basis that cannot pass that check, discard
/// the entire experimental epoch and cold-start the proven baseline
/// policy. This is deliberately broader than Phase-I cleanup because a
/// Phase-II perturbation can change basis membership after artificials are
/// gone, leaving no smaller transactional snapshot to restore.
pub fn restartSolveWithoutPerturbation(
    self: *SimplexEngine,
    problem: problem_module.ProblemView,
    control: SolveControl,
) SolveStatus {
    var baseline_control = control;
    baseline_control.initial_basis = null;
    baseline_control.degeneracy_strategy = .baseline;
    return self.solveProblem(problem, baseline_control);
}

/// Count exact-ratio ties after the ratio policy has already selected its
/// leaving row. Keeping this read-only diagnostic outside RatioTest makes
/// it impossible for instrumentation to perturb pivot selection.
pub fn countPrimalRatioTies(self: *const SimplexEngine, selected_step: f64) u32 {
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

pub fn prepareDegeneracyPolicy(self: *SimplexEngine, rows: usize, columns: usize) SolveStatus {
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

pub fn observeIterationStep(
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
    if (leaving_column != null and self.pivot_trace_count != 0)
        self.active_pivot_trace[self.pivot_trace_count - 1].bound_flip_count = bound_flip_count;
    if (leaving_column != null) {
        self.iteration_counters.committed_pivots += 1;
        if (self.cleanup_active) {
            self.iteration_counters.cleanup_pivots += 1;
        } else switch (self.current_phase) {
            .phase_one => self.iteration_counters.primal_phase_one_pivots += 1,
            .dual_phase_one => self.iteration_counters.dual_phase_one_pivots += 1,
            .dual_feasibility_repair => self.iteration_counters.dual_repair_pivots += 1,
            .phase_two => {
                if (self.shifted_dual_accounting_active)
                    self.iteration_counters.shifted_dual_pivots += 1
                else if (self.algorithm == .dual_revised)
                    self.iteration_counters.dual_phase_two_pivots += 1
                else
                    self.iteration_counters.primal_phase_two_pivots += 1;
            },
        }
    } else {
        self.iteration_counters.bound_moves += 1;
    }
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

pub fn basisFingerprint(self: *const SimplexEngine) u64 {
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

pub fn rememberDegenerateBasis(self: *SimplexEngine, fingerprint: u64) bool {
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
    const dual = engine.pricing.chooseDualLeaving(engine).?;
    try std.testing.expectEqual(@as(u32, 1), dual.row);
}
