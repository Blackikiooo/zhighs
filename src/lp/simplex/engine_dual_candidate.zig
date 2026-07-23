//! Hyper-sparse dual candidate-list maintenance for `SimplexEngine`.
//!
//! ## Responsibility
//!
//! Owns scoring, selection, and periodic rebuild of the small attractive-row
//! candidate set used by hyper-sparse dual pricing.

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

/// Classify one row as a dual leaving candidate under current basic bounds.
pub fn dualCandidate(self: *const SimplexEngine, row: usize) ?pricing_module.DualLeavingChoice {
    const basis = if (self.basis) |*value| value else return null;
    const value = basis.basic_value[row];
    if (value < basis.basic_lower[row] - self.numerical.primal_tolerance)
        return .{ .row = @intCast(row), .bound = .at_lower, .violation = basis.basic_lower[row] - value };
    if (value > basis.basic_upper[row] + self.numerical.primal_tolerance)
        return .{ .row = @intCast(row), .bound = .at_upper, .violation = value - basis.basic_upper[row] };
    return null;
}

/// Return the weighted CHUZR merit for a candidate row, or negative infinity.
pub fn dualCandidateScore(self: *const SimplexEngine, row_u32: u32) f64 {
    const row: usize = @intCast(row_u32);
    const basis = if (self.basis) |*value| value else return 0.0;
    const candidate = self.dualCandidate(row) orelse return 0.0;
    return candidate.violation / @sqrt(@max(basis.row_edge_weight[row], 1.0));
}

/// Select the highest-merit row from the active sparse candidate list.
pub fn bestDualCandidate(self: *SimplexEngine) ?pricing_module.DualLeavingChoice {
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

/// Rebuild the sparse CHUZR candidate set and density cutoff from all rows.
pub fn rebuildDualCandidateList(self: *SimplexEngine) void {
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
    const choice = engine.pricing.chooseDualLeaving(engine).?;
    try std.testing.expectEqual(@as(u32, 63), choice.row);
    try std.testing.expectEqual(@as(usize, 32), engine.dual_candidate_count);
    try std.testing.expect(engine.dual_candidate_cutoff > 0.0);
}
