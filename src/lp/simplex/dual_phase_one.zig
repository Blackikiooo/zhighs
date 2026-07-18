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
    basis_epoch: u64 = 0,
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

        self.allocator.free(self.saved_lower);
        self.allocator.free(self.saved_upper);
        self.allocator.free(self.work_cost);
        self.allocator.free(self.dual_infeasibility);
        self.allocator.free(self.perturbation);
        self.saved_lower = lower;
        self.saved_upper = upper;
        self.work_cost = cost;
        self.dual_infeasibility = infeasibility;
        self.perturbation = perturbation;
    }

    pub fn begin(self: *DualPhaseOneWorkspace, basis: *basis_module.BasisState, original_count: usize) !void {
        try self.ensureCapacity(original_count);
        @memcpy(self.saved_lower[0..original_count], basis.col_lower[0..original_count]);
        @memcpy(self.saved_upper[0..original_count], basis.col_upper[0..original_count]);
        @memset(self.dual_infeasibility[0..original_count], 0.0);
        @memset(self.perturbation[0..original_count], 0.0);
        self.basis_epoch +%= 1;
        self.active = true;
    }

    /// Map original bounds to the bounded dual Phase-I problem. Nonbasic
    /// values are selected from the reduced-cost sign, making the working
    /// basis dual feasible while its objective is minus dual infeasibility.
    pub fn installWorkingBounds(
        self: *DualPhaseOneWorkspace,
        basis: *basis_module.BasisState,
        original_count: usize,
    ) void {
        for (0..original_count) |column| {
            const lower = self.saved_lower[column];
            const upper = self.saved_upper[column];
            const lower_finite = std.math.isFinite(lower);
            const upper_finite = std.math.isFinite(upper);
            if (!lower_finite and !upper_finite) {
                basis.col_lower[column] = -1000.0;
                basis.col_upper[column] = 1000.0;
            } else if (!lower_finite) {
                basis.col_lower[column] = -1.0;
                basis.col_upper[column] = 0.0;
            } else if (!upper_finite) {
                basis.col_lower[column] = 0.0;
                basis.col_upper[column] = 1.0;
            } else if (lower == upper) {
                basis.col_lower[column] = 0.0;
                basis.col_upper[column] = 0.0;
            } else {
                // HiGHS can collapse boxed columns because its separate
                // nonbasicMove machinery flips their original bound before
                // Phase I. zhighs stores move in BasisStatus, so retain a
                // symmetric unit interval to preserve the same entering
                // choices without adding a second membership representation.
                basis.col_lower[column] = -1.0;
                basis.col_upper[column] = 1.0;
            }
            if (basis.col_status[column] == .basic) continue;

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

            if (basis.col_lower[column] == basis.col_upper[column]) {
                basis.col_status[column] = .fixed;
                basis.primal[column] = basis.col_lower[column];
            } else if (!lower_finite and !upper_finite) {
                basis.col_status[column] = .free;
                basis.primal[column] = 0.0;
            } else if (!lower_finite) {
                basis.col_status[column] = .at_upper;
                basis.primal[column] = basis.col_upper[column];
            } else if (!upper_finite) {
                basis.col_status[column] = .at_lower;
                basis.primal[column] = basis.col_lower[column];
            } else {
                // Preserve the original boxed move. Cost correction performed
                // before this transition makes either move dual feasible.
                basis.primal[column] = if (basis.col_status[column] == .at_upper)
                    basis.col_upper[column]
                else
                    basis.col_lower[column];
            }
        }
        for (basis.basic_index, 0..) |column, row| {
            basis.basic_lower[row] = basis.col_lower[column];
            basis.basic_upper[row] = basis.col_upper[column];
        }
    }

    /// Restore the model bounds and choose a deterministic original-bound
    /// status from the fresh reduced cost. Basic membership is unchanged.
    pub fn restoreOriginalBounds(self: *DualPhaseOneWorkspace, basis: *basis_module.BasisState, original_count: usize) void {
        @memcpy(basis.col_lower[0..original_count], self.saved_lower[0..original_count]);
        @memcpy(basis.col_upper[0..original_count], self.saved_upper[0..original_count]);
        for (0..original_count) |column| {
            if (basis.col_status[column] == .basic) continue;
            const lower = basis.col_lower[column];
            const upper = basis.col_upper[column];
            if (lower == upper) {
                basis.col_status[column] = .fixed;
                basis.primal[column] = lower;
            } else if (!std.math.isFinite(lower) and !std.math.isFinite(upper)) {
                basis.col_status[column] = .free;
                basis.primal[column] = 0.0;
            } else if (!std.math.isFinite(lower)) {
                basis.col_status[column] = .at_upper;
                basis.primal[column] = upper;
            } else if (!std.math.isFinite(upper) or basis.reduced_cost[column] >= 0.0) {
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

    pub fn requestedBytes(self: *const DualPhaseOneWorkspace) usize {
        return std.mem.sliceAsBytes(self.saved_lower).len +
            std.mem.sliceAsBytes(self.saved_upper).len +
            std.mem.sliceAsBytes(self.work_cost).len +
            std.mem.sliceAsBytes(self.dual_infeasibility).len +
            std.mem.sliceAsBytes(self.perturbation).len;
    }
};

test "dual Phase-I workspace maps lower, upper, free and boxed bounds" {
    var basis = try basis_module.BasisState.init(std.testing.allocator, 1, 4);
    defer basis.deinit();
    var workspace = DualPhaseOneWorkspace.init(std.testing.allocator);
    defer workspace.deinit();
    basis.col_lower[0..4].* = .{ 0.0, -std.math.inf(f64), -std.math.inf(f64), -2.0 };
    basis.col_upper[0..4].* = .{ std.math.inf(f64), 3.0, std.math.inf(f64), 2.0 };
    basis.col_status[0..4].* = .{ .at_lower, .at_upper, .free, .at_lower };
    basis.reduced_cost[0..4].* = .{ -2.0, 3.0, -4.0, -5.0 };
    try workspace.begin(&basis, 4);
    workspace.installWorkingBounds(&basis, 4);
    try std.testing.expectEqualSlices(f64, &.{ 0.0, -1.0, -1000.0, -1.0 }, basis.col_lower[0..4]);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 0.0, 1000.0, 1.0 }, basis.col_upper[0..4]);
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 0.0, 0.0, -1.0 }, basis.primal[0..4]);
    workspace.restoreOriginalBounds(&basis, 4);
    try std.testing.expectEqualSlices(f64, &.{ 0.0, -std.math.inf(f64), -std.math.inf(f64), -2.0 }, basis.col_lower[0..4]);
}
