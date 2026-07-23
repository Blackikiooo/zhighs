//! HiGHS-compatible control state for the serial dual simplex.
//!
//! This state deliberately contains no matrix or basis storage. It is the
//! single owner of phase, rebuild and cost-lifecycle decisions that otherwise
//! become scattered boolean combinations across Phase I, Phase II and cleanup.

const std = @import("std");

/// Major state of the serial dual-simplex controller.
pub const SolvePhase = enum {
    unknown,
    phase_one,
    phase_two,
    optimal,
    optimal_cleanup,
    primal_infeasible_cleanup,
    exit,
    taboo_basis,
    failed,
};

/// First event that requires the current iteration loop to stop and rebuild.
pub const RebuildReason = enum {
    none,
    possibly_optimal,
    possibly_dual_unbounded,
    choose_column_failed,
    excessive_primal_value,
    possibly_singular_basis,
    update_limit,
    cleanup,
    phase_transition,
};

/// Phase and cost-lifecycle state shared by dual Phase I, Phase II and cleanup.
pub const ControlState = struct {
    /// Active major phase; `.unknown` is the reset state before initialization.
    phase: SolvePhase = .unknown,
    /// Pending rebuild request. The first request is preserved until consumed.
    rebuild_reason: RebuildReason = .none,
    /// True only after factorization and all rebuild-dependent vectors are fresh.
    has_fresh_rebuild: bool = false,
    /// Skip Phase I when the unperturbed assessment safely permits Phase II.
    force_phase_two: bool = false,
    /// Whether `workCost` currently contains deterministic perturbations.
    costs_perturbed: bool = false,
    /// Whether one or more temporary dual-feasibility shifts remain installed.
    costs_shifted: bool = false,
    /// Policy switch allowing perturbation during this solve.
    allow_cost_perturbation: bool = true,
    /// Policy switch allowing temporary cost shifts during this solve.
    allow_cost_shifting: bool = true,
    /// Number of increasingly strict cleanup attempts made while leaving Phase I.
    phase_one_cleanup_level: u16 = 0,
    /// Number of original-objective cleanup attempts made after Phase II.
    cleanup_level: u16 = 0,
    /// Count of unperturbed dual-infeasible nonbasic variables at initialization.
    unperturbed_dual_infeasibility_count: usize = 0,
    /// Largest unperturbed dual infeasibility.
    unperturbed_dual_infeasibility_max: f64 = 0.0,
    /// Sum of unperturbed dual infeasibilities.
    unperturbed_dual_infeasibility_sum: f64 = 0.0,
    /// Count of primal-infeasible basic variables at initialization.
    initial_primal_infeasibility_count: usize = 0,
    /// Largest initial basic-bound violation.
    initial_primal_infeasibility_max: f64 = 0.0,
    /// Sum of initial basic-bound violations.
    initial_primal_infeasibility_sum: f64 = 0.0,
    /// Cached initialization decision used by perturbation and phase selection.
    near_optimal: bool = false,

    /// Restore all controller fields to their cold-solve defaults.
    pub fn reset(self: *ControlState) void {
        self.* = .{};
    }

    /// Enter dual Phase I or II and require a transition rebuild.
    pub fn enterPhase(self: *ControlState, phase: SolvePhase) void {
        std.debug.assert(phase == .phase_one or phase == .phase_two);
        self.phase = phase;
        self.rebuild_reason = .phase_transition;
        self.has_fresh_rebuild = false;
    }

    /// Preserve the first reason, matching HEkkDual's early-return pipeline:
    /// later iteration stages must not overwrite the event that stopped it.
    pub fn requestRebuild(self: *ControlState, reason: RebuildReason) void {
        std.debug.assert(reason != .none);
        if (self.rebuild_reason == .none) self.rebuild_reason = reason;
        self.has_fresh_rebuild = false;
    }

    /// Consume the pending reason and mark rebuild-dependent state stale.
    pub fn beginRebuild(self: *ControlState) RebuildReason {
        const reason = self.rebuild_reason;
        self.rebuild_reason = .none;
        self.has_fresh_rebuild = false;
        return reason;
    }

    /// Publish successful completion of the complete rebuild pipeline.
    pub fn finishRebuild(self: *ControlState) void {
        std.debug.assert(self.rebuild_reason == .none);
        self.has_fresh_rebuild = true;
    }
};

test "dual control state preserves the first rebuild reason" {
    var state = ControlState{};
    state.enterPhase(.phase_one);
    try std.testing.expectEqual(RebuildReason.phase_transition, state.beginRebuild());
    state.finishRebuild();
    try std.testing.expect(state.has_fresh_rebuild);

    state.requestRebuild(.possibly_optimal);
    state.requestRebuild(.update_limit);
    try std.testing.expectEqual(RebuildReason.possibly_optimal, state.beginRebuild());
    state.finishRebuild();
    try std.testing.expectEqual(SolvePhase.phase_one, state.phase);
}
