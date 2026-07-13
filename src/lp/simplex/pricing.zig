//! Entering-variable pricing policies.

const std = @import("std");
const basis = @import("basis.zig");

pub const PricingRule = enum { dantzig, devex, steepest_edge, partial, hyper_sparse };
pub const EnteringChoice = struct { column: u32, direction: f64 };
pub const DualLeavingChoice = struct {
    row: u32,
    bound: basis.BasisStatus,
    violation: f64,
};
pub const Pricing = struct {
    rule: PricingRule = .devex,
    devex_reset_period: usize = 100,
    iterations: usize = 0,

    fn normalizedScore(self: Pricing, violation: f64, weight: f64) f64 {
        return switch (self.rule) {
            .dantzig, .partial => violation,
            .devex, .steepest_edge, .hyper_sparse => violation / @sqrt(@max(weight, 1.0)),
        };
    }

    pub fn chooseEntering(self: *Pricing, reduced_cost: []const f64, tolerance: f64) ?u32 {
        self.iterations += 1;
        var best: ?u32 = null;
        var best_value = tolerance;
        for (reduced_cost, 0..) |value, i| {
            const magnitude = @abs(value);
            if (magnitude > best_value) {
                best_value = magnitude;
                best = @intCast(i);
            }
        }
        return best;
    }

    /// Select an entering nonbasic column without allocating. Basic columns
    /// are excluded to preserve basis invariants during revised iterations.
    pub fn chooseEnteringStatus(self: *Pricing, reduced_cost: []const f64, status: []const basis.BasisStatus, tolerance: f64) ?u32 {
        if (status.len != reduced_cost.len) return null;
        self.iterations += 1;
        var best: ?u32 = null;
        var best_value = tolerance;
        for (reduced_cost, status, 0..) |value, column_status, i| {
            if (column_status == .basic) continue;
            const magnitude = @abs(value);
            if (magnitude > best_value) {
                best_value = magnitude;
                best = @intCast(i);
            }
        }
        return best;
    }

    /// Primal minimization pricing with bound-aware reduced-cost signs.
    pub fn choosePrimalEntering(self: *Pricing, reduced_cost: []const f64, status: []const basis.BasisStatus, tolerance: f64) ?u32 {
        const choice = self.choosePrimalEnteringDirection(reduced_cost, status, tolerance) orelse return null;
        return choice.column;
    }

    pub fn choosePrimalEnteringDirection(self: *Pricing, reduced_cost: []const f64, status: []const basis.BasisStatus, tolerance: f64) ?EnteringChoice {
        return self.choosePrimalEnteringWeighted(reduced_cost, status, &.{}, tolerance);
    }

    /// Bound-aware primal pricing with optional caller-owned edge weights.
    /// Empty weights select Dantzig-compatible unit weights. No candidates are
    /// materialized, which keeps the scan allocation-free and cache-linear.
    pub fn choosePrimalEnteringWeighted(self: *Pricing, reduced_cost: []const f64, status: []const basis.BasisStatus, edge_weight: []const f64, tolerance: f64) ?EnteringChoice {
        if (status.len < reduced_cost.len) return null;
        if (edge_weight.len != 0 and edge_weight.len < reduced_cost.len) return null;
        self.iterations += 1;
        var best: ?EnteringChoice = null;
        var best_score = tolerance;
        for (reduced_cost, status[0..reduced_cost.len], 0..) |value, column_status, i| {
            const violation = switch (column_status) {
                .at_lower => -value,
                .at_upper => value,
                .free, .superbasic => @abs(value),
                .basic, .fixed => continue,
            };
            const weight = if (edge_weight.len == 0) 1.0 else edge_weight[i];
            const score = self.normalizedScore(violation, weight);
            if (score > best_score) {
                best_score = score;
                const direction: f64 = switch (column_status) {
                    .at_upper => -1.0,
                    .free, .superbasic => if (value < 0) 1.0 else -1.0,
                    else => 1.0,
                };
                best = .{ .column = @intCast(i), .direction = direction };
            }
        }
        return best;
    }

    /// Dantzig-style dual pricing over the dense basic-value SoA. The most
    /// primal-infeasible row leaves; no temporary candidate objects are built.
    pub fn chooseDualLeaving(self: *Pricing, value: []const f64, lower: []const f64, upper: []const f64, tolerance: f64) ?DualLeavingChoice {
        return self.chooseDualLeavingWeighted(value, lower, upper, &.{}, tolerance);
    }

    pub fn chooseDualLeavingWeighted(self: *Pricing, value: []const f64, lower: []const f64, upper: []const f64, edge_weight: []const f64, tolerance: f64) ?DualLeavingChoice {
        if (value.len != lower.len or value.len != upper.len) return null;
        if (edge_weight.len != 0 and edge_weight.len != value.len) return null;
        self.iterations += 1;
        var best: ?DualLeavingChoice = null;
        var best_score = tolerance;
        for (value, lower, upper, 0..) |basic_value, lb, ub, row| {
            const choice: ?DualLeavingChoice = if (basic_value < lb - tolerance)
                .{ .row = @intCast(row), .bound = .at_lower, .violation = lb - basic_value }
            else if (basic_value > ub + tolerance)
                .{ .row = @intCast(row), .bound = .at_upper, .violation = basic_value - ub }
            else
                null;
            if (choice) |candidate| {
                const weight = if (edge_weight.len == 0) 1.0 else edge_weight[row];
                const score = self.normalizedScore(candidate.violation, weight);
                if (score > best_score) {
                    best_score = score;
                    best = candidate;
                }
            }
        }
        return best;
    }
};

test {
    std.testing.refAllDecls(@This());
}

test "pricing excludes basic columns" {
    var pricing = Pricing{};
    const costs = [_]f64{ 100.0, 3.0, 8.0 };
    const statuses = [_]basis.BasisStatus{ .basic, .at_lower, .at_lower };
    try std.testing.expectEqual(@as(?u32, 2), pricing.chooseEnteringStatus(&costs, &statuses, 1e-9));
}

test "dual pricing chooses the largest bound violation" {
    var pricing = Pricing{};
    const choice = pricing.chooseDualLeaving(
        &[_]f64{ -2, 5 },
        &[_]f64{ 0, 0 },
        &[_]f64{ 10, 3 },
        1e-9,
    ).?;
    try std.testing.expectEqual(@as(u32, 0), choice.row);
    try std.testing.expectEqual(basis.BasisStatus.at_lower, choice.bound);
}

test "devex weights can change the chosen entering column" {
    var pricing = Pricing{ .rule = .devex };
    const choice = pricing.choosePrimalEnteringWeighted(
        &[_]f64{ -4, -3 },
        &[_]basis.BasisStatus{ .at_lower, .at_lower },
        &[_]f64{ 100, 1 },
        1e-9,
    ).?;
    try std.testing.expectEqual(@as(u32, 1), choice.column);
}
