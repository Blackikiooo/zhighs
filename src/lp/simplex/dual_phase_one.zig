//! Persistent SoA workspace for dual-simplex Phase I.
//!
//! Phase I changes bounds and costs, but never the borrowed model matrix. The
//! saved model bounds make entry/exit transactional, while the other arrays
//! are reused by every basis epoch and keep the iteration path allocation-free.

const std = @import("std");
const basis_module = @import("basis.zig");

pub const DualPhaseOneWorkspace = struct {
    allocator: std.mem.Allocator,
    saved_lower: []f64 = &.{},
    saved_upper: []f64 = &.{},
    work_cost: []f64 = &.{},
    dual_infeasibility: []f64 = &.{},
    perturbation: []f64 = &.{},
    /// Explicit Phase-I nonbasic direction: +1 moves up from the lower
    /// bound, -1 moves down from the upper bound, and zero is basic/fixed/free.
    nonbasic_move: []i8 = &.{},
    /// Residual basic-bound violation after the current ratio-test flip set.
    remaining_violation: []f64 = &.{},
    basis_epoch: u64 = 0,
    /// Basic-coordinate envelope in the engine's scaled coordinates when the
    /// current Phase-I epoch was installed.
    working_radius: f64 = 1.0,
    active: bool = false,

    pub fn init(allocator: std.mem.Allocator) DualPhaseOneWorkspace {
        return .{ .allocator = allocator };
    }

    pub fn deinit(self: *DualPhaseOneWorkspace) void {
        self.allocator.free(self.saved_lower);
        self.allocator.free(self.saved_upper);
        self.allocator.free(self.work_cost);
        self.allocator.free(self.dual_infeasibility);
        self.allocator.free(self.perturbation);
        self.allocator.free(self.nonbasic_move);
        self.allocator.free(self.remaining_violation);
        self.* = .{ .allocator = self.allocator };
    }

    /// Grow all SoA arrays together. Existing capacity is retained across
    /// solves, so phase transitions allocate only when the model grows.
    pub fn ensureCapacity(self: *DualPhaseOneWorkspace, count: usize) !void {
        if (self.saved_lower.len >= count) return;
        const lower = try self.allocator.alloc(f64, count);
        errdefer self.allocator.free(lower);
        const upper = try self.allocator.alloc(f64, count);
        errdefer self.allocator.free(upper);
        const cost = try self.allocator.alloc(f64, count);
        errdefer self.allocator.free(cost);
        const infeasibility = try self.allocator.alloc(f64, count);
        errdefer self.allocator.free(infeasibility);
        const perturbation = try self.allocator.alloc(f64, count);
        errdefer self.allocator.free(perturbation);
        const move = try self.allocator.alloc(i8, count);
        errdefer self.allocator.free(move);
        const remaining = try self.allocator.alloc(f64, count);
        errdefer self.allocator.free(remaining);

        self.allocator.free(self.saved_lower);
        self.allocator.free(self.saved_upper);
        self.allocator.free(self.work_cost);
        self.allocator.free(self.dual_infeasibility);
        self.allocator.free(self.perturbation);
        self.allocator.free(self.nonbasic_move);
        self.allocator.free(self.remaining_violation);
        self.saved_lower = lower;
        self.saved_upper = upper;
        self.work_cost = cost;
        self.dual_infeasibility = infeasibility;
        self.perturbation = perturbation;
        self.nonbasic_move = move;
        self.remaining_violation = remaining;
    }

    pub fn begin(self: *DualPhaseOneWorkspace, basis: *basis_module.BasisState, original_count: usize) !void {
        try self.ensureCapacity(original_count);
        @memcpy(self.saved_lower[0..original_count], basis.col_lower[0..original_count]);
        @memcpy(self.saved_upper[0..original_count], basis.col_upper[0..original_count]);
        @memset(self.dual_infeasibility[0..original_count], 0.0);
        @memset(self.perturbation[0..original_count], 0.0);
        @memset(self.nonbasic_move[0..original_count], 0);
        @memset(self.remaining_violation[0..original_count], 0.0);
        self.basis_epoch +%= 1;
        self.active = true;
    }

    /// Start a work-cost-only epoch. Bounds are snapshotted for perturbation
    /// policy, but no Phase-I bounds are installed and no deferred restore is
    /// required.
    pub fn beginCostEpoch(self: *DualPhaseOneWorkspace, basis: *basis_module.BasisState, original_count: usize) !void {
        try self.ensureCapacity(original_count);
        @memcpy(self.saved_lower[0..original_count], basis.col_lower[0..original_count]);
        @memcpy(self.saved_upper[0..original_count], basis.col_upper[0..original_count]);
        @memset(self.dual_infeasibility[0..original_count], 0.0);
        @memset(self.perturbation[0..original_count], 0.0);
        @memset(self.nonbasic_move[0..original_count], 0);
        @memset(self.remaining_violation[0..original_count], 0.0);
        self.basis_epoch +%= 1;
        self.active = false;
    }

    /// Map original bounds to the dual Phase-I subproblem (HiGHS-style).
    ///
    /// Working bounds:
    ///   FREE           → [-1000, 1000]
    ///   lower-unbounded → [-1, 0]
    ///   upper-unbounded → [0, 1]
    ///   boxed / fixed  → [0, 0]    (collapsed)
    ///
    /// Nonbasic primal values encode infeasibility (NOT the cost):
    ///   dual-infeasible → ±1  (sign from required direction)
    ///   dual-feasible   →  0
    /// The dual objective is then Σ value_j · d_j = −Σ infeasibilities.
    ///
    /// Boxed columns collapse to [0, 0] and have zero Phase-I move, exactly
    /// like HiGHS' initialiseNonbasicValueAndMove. Their original side is
    /// selected again only after the Phase-II bounds have been restored.
    pub fn installWorkingBounds(
        self: *DualPhaseOneWorkspace,
        basis: *basis_module.BasisState,
        original_count: usize,
    ) void {
        self.working_radius = 1000.0;

        for (0..original_count) |column| {
            const lower = self.saved_lower[column];
            const upper = self.saved_upper[column];
            const lower_finite = std.math.isFinite(lower);
            const upper_finite = std.math.isFinite(upper);

            // ── Set working bounds ──
            if (!lower_finite and !upper_finite) {
                basis.col_lower[column] = -1000.0;
                basis.col_upper[column] = 1000.0;
            } else if (!lower_finite) {
                basis.col_lower[column] = -1.0;
                basis.col_upper[column] = 0.0;
            } else if (!upper_finite) {
                basis.col_lower[column] = 0.0;
                basis.col_upper[column] = 1.0;
            } else {
                // Both finite: collapse to [0, 0] regardless of equality.
                // nonbasic_move records the original side.
                basis.col_lower[column] = 0.0;
                basis.col_upper[column] = 0.0;
            }

            if (basis.col_status[column] == .basic) continue;

            // ── Infeasibility from original reduced cost ──
            const reduced = basis.reduced_cost[column];
            const infeasibility = if (!lower_finite and !upper_finite)
                @abs(reduced)
            else if (!lower_finite)
                @max(reduced, 0.0)
            else if (!upper_finite)
                @max(-reduced, 0.0)
            else switch (basis.col_status[column]) {
                .at_lower => @max(-reduced, 0.0),
                .at_upper => @max(reduced, 0.0),
                else => 0.0,
            };
            self.dual_infeasibility[column] = infeasibility;

            // ── Nonbasic primal: ±1 for infeasible, 0 for feasible ──
            if (basis.col_lower[column] == basis.col_upper[column]) {
                // Fixed or collapsed boxed variables do not participate in
                // the Phase-I ratio test or bound flips.
                self.nonbasic_move[column] = 0;
                basis.col_status[column] = .fixed;
                basis.primal[column] = 0.0;
            } else if (!lower_finite and !upper_finite) {
                // Free: value = sign of reduced cost (or 0 if none)
                basis.col_status[column] = .free;
                if (reduced < -0.0) {
                    basis.primal[column] = -1.0;
                    self.nonbasic_move[column] = -1;
                } else if (reduced > 0.0) {
                    basis.primal[column] = 1.0;
                    self.nonbasic_move[column] = 1;
                } else {
                    basis.primal[column] = 0.0;
                    self.nonbasic_move[column] = 0;
                }
            } else if (!lower_finite) {
                // Lower-unbounded: working [-1, 0], value = sign(reduced)
                basis.col_status[column] = .at_upper;
                if (reduced > 0.0) {
                    basis.primal[column] = -1.0;
                    self.nonbasic_move[column] = -1;
                } else {
                    basis.primal[column] = 0.0;
                    self.nonbasic_move[column] = 0;
                }
            } else {
                // Upper-unbounded: working [0, 1]
                basis.col_status[column] = .at_lower;
                if (reduced < -0.0) {
                    basis.primal[column] = 1.0;
                    self.nonbasic_move[column] = 1;
                } else {
                    basis.primal[column] = 0.0;
                    self.nonbasic_move[column] = 0;
                }
            }
        }
        for (basis.basic_index, 0..) |column, row| {
            basis.basic_lower[row] = basis.col_lower[column];
            basis.basic_upper[row] = basis.col_upper[column];
        }
    }

    /// Restore the original model bounds. nonbasic_move (tracked through
    /// the Phase-I pivot path) determines the side for boxed columns that
    /// were collapsed to [0, 0]; for other columns the fresh reduced-cost
    /// sign is used. Basic membership is unchanged.
    pub fn restoreOriginalBounds(self: *DualPhaseOneWorkspace, basis: *basis_module.BasisState, original_count: usize) void {
        @memcpy(basis.col_lower[0..original_count], self.saved_lower[0..original_count]);
        @memcpy(basis.col_upper[0..original_count], self.saved_upper[0..original_count]);
        for (0..original_count) |column| {
            if (basis.col_status[column] == .basic) continue;
            const lower = basis.col_lower[column];
            const upper = basis.col_upper[column];
            const lower_finite = std.math.isFinite(lower);
            const upper_finite = std.math.isFinite(upper);
            const move = self.nonbasic_move[column];
            if (lower == upper) {
                basis.col_status[column] = .fixed;
                basis.primal[column] = lower;
            } else if (!lower_finite and !upper_finite) {
                basis.col_status[column] = .free;
                basis.primal[column] = 0.0;
            } else if (!lower_finite) {
                basis.col_status[column] = .at_upper;
                basis.primal[column] = upper;
            } else if (!upper_finite) {
                basis.col_status[column] = .at_lower;
                basis.primal[column] = lower;
            } else if (move != 0) {
                // Boxed column that was collapsed: use the Phase-I-tracked
                // direction. move > 0 → at_lower, move < 0 → at_upper.
                if (move > 0) {
                    basis.col_status[column] = .at_lower;
                    basis.primal[column] = lower;
                } else {
                    basis.col_status[column] = .at_upper;
                    basis.primal[column] = upper;
                }
            } else if (basis.reduced_cost[column] >= 0.0) {
                basis.col_status[column] = .at_lower;
                basis.primal[column] = lower;
            } else {
                basis.col_status[column] = .at_upper;
                basis.primal[column] = upper;
            }
        }
        for (basis.basic_index, 0..) |column, row| {
            basis.basic_lower[row] = basis.col_lower[column];
            basis.basic_upper[row] = basis.col_upper[column];
        }
        self.active = false;
    }

    pub fn noteBoundFlip(self: *DualPhaseOneWorkspace, column: usize) void {
        if (!self.active or column >= self.nonbasic_move.len) return;
        self.nonbasic_move[column] = -self.nonbasic_move[column];
    }

    pub fn notePivot(
        self: *DualPhaseOneWorkspace,
        entering_column: usize,
        leaving_column: usize,
        leaving_bound: basis_module.BasisStatus,
    ) void {
        if (!self.active or entering_column >= self.nonbasic_move.len or leaving_column >= self.nonbasic_move.len) return;
        self.nonbasic_move[entering_column] = 0;
        self.nonbasic_move[leaving_column] = switch (leaving_bound) {
            .at_lower => 1,
            .at_upper => -1,
            .fixed => 0,
            else => 0,
        };
    }

    pub fn recordRemainingViolation(
        self: *DualPhaseOneWorkspace,
        row: usize,
        violation: f64,
        tableau: []const f64,
        flip_columns: []const u32,
        lower: []const f64,
        upper: []const f64,
    ) void {
        if (!self.active or row >= self.remaining_violation.len) return;
        var corrected: f64 = 0.0;
        for (flip_columns) |column_u32| {
            const column: usize = @intCast(column_u32);
            if (column >= tableau.len or column >= lower.len or column >= upper.len) continue;
            corrected += @abs(tableau[column]) * (upper[column] - lower[column]);
        }
        self.remaining_violation[row] = @max(violation - corrected, 0.0);
    }

    pub fn requestedBytes(self: *const DualPhaseOneWorkspace) usize {
        return std.mem.sliceAsBytes(self.saved_lower).len +
            std.mem.sliceAsBytes(self.saved_upper).len +
            std.mem.sliceAsBytes(self.work_cost).len +
            std.mem.sliceAsBytes(self.dual_infeasibility).len +
            std.mem.sliceAsBytes(self.perturbation).len +
            std.mem.sliceAsBytes(self.nonbasic_move).len +
            std.mem.sliceAsBytes(self.remaining_violation).len;
    }
};

test "dual Phase-I HiGHS-style bounds: boxed collapsed, values encode infeasibility" {
    var basis = try basis_module.BasisState.init(std.testing.allocator, 1, 4);
    defer basis.deinit();
    var workspace = DualPhaseOneWorkspace.init(std.testing.allocator);
    defer workspace.deinit();
    // col0: LOWER  [0, inf]   rc=-2 at_lower   → infeasible → value= 1
    // col1: UPPER  [-inf, 3]  rc= 3 at_upper   → infeasible → value=-1
    // col2: FREE   [-inf, inf] rc=-4 free       → infeasible → value=-1
    // col3: BOXED  [-2, 2]    rc=-5 at_lower   → collapsed [0,0], move=1
    basis.col_lower[0..4].* = .{ 0.0, -std.math.inf(f64), -std.math.inf(f64), -2.0 };
    basis.col_upper[0..4].* = .{ std.math.inf(f64), 3.0, std.math.inf(f64), 2.0 };
    basis.col_status[0..4].* = .{ .at_lower, .at_upper, .free, .at_lower };
    basis.reduced_cost[0..4].* = .{ -2.0, 3.0, -4.0, -5.0 };
    try workspace.begin(&basis, 4);
    workspace.installWorkingBounds(&basis, 4);
    // Working bounds: LOWER [0,1]  UPPER [-1,0]  FREE [-1000,1000]  BOXED [0,0]
    try std.testing.expectEqualSlices(f64, &.{ 0.0, -1.0, -1000.0, 0.0 }, basis.col_lower[0..4]);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 0.0, 1000.0, 0.0 }, basis.col_upper[0..4]);
    // Primal values encode infeasibility: ±1 infeasible, 0 feasible
    try std.testing.expectEqualSlices(f64, &.{ 1.0, -1.0, -1.0, 0.0 }, basis.primal[0..4]);
    // nonbasic_move: LOWER=+1, UPPER=-1, FREE=-1, BOXED=+1 (original at_lower)
    try std.testing.expectEqualSlices(i8, &.{ 1, -1, -1, 1 }, workspace.nonbasic_move[0..4]);
    workspace.noteBoundFlip(0);
    try std.testing.expectEqual(@as(i8, -1), workspace.nonbasic_move[0]);
    workspace.notePivot(1, 3, .at_upper);
    try std.testing.expectEqual(@as(i8, 0), workspace.nonbasic_move[1]);
    try std.testing.expectEqual(@as(i8, -1), workspace.nonbasic_move[3]);
    workspace.restoreOriginalBounds(&basis, 4);
    try std.testing.expectEqualSlices(f64, &.{ 0.0, -std.math.inf(f64), -std.math.inf(f64), -2.0 }, basis.col_lower[0..4]);
}

test "dual Phase-I records violation left after bound flips" {
    var workspace = DualPhaseOneWorkspace.init(std.testing.allocator);
    defer workspace.deinit();
    var basis = try basis_module.BasisState.init(std.testing.allocator, 1, 1);
    defer basis.deinit();
    try workspace.begin(&basis, 2);

    workspace.recordRemainingViolation(
        0,
        10,
        &.{ 2, -1 },
        &.{0},
        &.{ 0, 0 },
        &.{ 3, 4 },
    );
    try std.testing.expectApproxEqAbs(@as(f64, 4), workspace.remaining_violation[0], 1e-12);
}

test "dual Phase-I HiGHS bounds: boxed collapsed, basic bounds follow" {
    var basis = try basis_module.BasisState.init(std.testing.allocator, 1, 0);
    defer basis.deinit();
    var workspace = DualPhaseOneWorkspace.init(std.testing.allocator);
    defer workspace.deinit();
    basis.basic_index[0] = 0;
    basis.col_status[0] = .basic;
    basis.col_lower[0] = 0;
    basis.col_upper[0] = 1;
    basis.basic_value[0] = 12; // irrelevant: bounds are hardcoded

    try workspace.begin(&basis, 1);
    workspace.installWorkingBounds(&basis, 1);

    // BOXED column → collapsed [0, 0]; basic bounds follow column bounds
    try std.testing.expectApproxEqAbs(@as(f64, 1000), workspace.working_radius, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0), basis.col_lower[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0), basis.col_upper[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0), basis.basic_lower[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0), basis.basic_upper[0], 1e-12);
}
