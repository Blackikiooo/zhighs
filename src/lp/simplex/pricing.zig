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
    /// Persistent segmented-pricing cursor. A segment is selected once per
    /// pricing operation; only an all-segment miss can certify optimality.
    partial_candidate_count: usize = 0,
    partial_cached_searches: usize = 0,
    partial_refill_interval: usize = 1,
    partial_searches: usize = 0,
    partial_scanned_entries: usize = 0,
    partial_full_scans: usize = 0,

    pub fn resetPartial(self: *Pricing) void {
        self.partial_candidate_count = 0;
        self.partial_cached_searches = 0;
        self.partial_searches = 0;
        self.partial_scanned_entries = 0;
        self.partial_full_scans = 0;
    }

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

    /// Bland fallback used only after repeated degenerate pivots. The first
    /// improving nonbasic column is selected, making the entering decision
    /// independent of weights and previous candidate ordering.
    pub fn choosePrimalEnteringBland(self: *Pricing, reduced_cost: []const f64, status: []const basis.BasisStatus, tolerance: f64) ?EnteringChoice {
        if (status.len < reduced_cost.len) return null;
        self.iterations += 1;
        for (reduced_cost, status[0..reduced_cost.len], 0..) |value, column_status, column| {
            const direction: f64 = switch (column_status) {
                .at_lower => if (-value > tolerance) 1.0 else continue,
                .at_upper => if (value > tolerance) -1.0 else continue,
                .free, .superbasic => if (@abs(value) > tolerance) (if (value < 0.0) 1.0 else -1.0) else continue,
                .basic, .fixed => continue,
            };
            return .{ .column = @intCast(column), .direction = direction };
        }
        return null;
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

    /// Multiple pricing over a caller-owned persistent candidate array. A
    /// global refill records every currently improving column. Later pivots
    /// rescore only that pool until it is exhausted, at which point another
    /// full scan is required before optimality can be certified.
    pub fn choosePrimalEnteringMultiple(
        self: *Pricing,
        reduced_cost: []const f64,
        status: []const basis.BasisStatus,
        edge_weight: []const f64,
        candidates: []u32,
        tolerance: f64,
    ) ?EnteringChoice {
        if (status.len < reduced_cost.len or candidates.len < reduced_cost.len) return null;
        if (edge_weight.len != 0 and edge_weight.len < reduced_cost.len) return null;
        self.iterations += 1;
        self.partial_searches += 1;

        if (self.partial_candidate_count != 0 and
            (self.partial_refill_interval == 0 or self.partial_cached_searches < self.partial_refill_interval))
        {
            var retained: usize = 0;
            var best: ?EnteringChoice = null;
            var best_score = tolerance;
            for (candidates[0..self.partial_candidate_count]) |column_u32| {
                const column: usize = @intCast(column_u32);
                const choice = primalCandidate(reduced_cost[column], status[column], tolerance) orelse continue;
                candidates[retained] = column_u32;
                retained += 1;
                const weight = if (edge_weight.len == 0) 1.0 else edge_weight[column];
                const score = self.normalizedScore(choice.violation, weight);
                if (score > best_score) {
                    best_score = score;
                    best = .{ .column = column_u32, .direction = choice.direction };
                }
            }
            self.partial_scanned_entries += self.partial_candidate_count;
            self.partial_candidate_count = retained;
            if (best) |choice| {
                self.partial_cached_searches += 1;
                return choice;
            }
        }

        var best: ?EnteringChoice = null;
        var best_score = tolerance;
        var count: usize = 0;
        for (reduced_cost, status[0..reduced_cost.len], 0..) |value, column_status, column| {
            const choice = primalCandidate(value, column_status, tolerance) orelse continue;
            candidates[count] = @intCast(column);
            count += 1;
            const weight = if (edge_weight.len == 0) 1.0 else edge_weight[column];
            const score = self.normalizedScore(choice.violation, weight);
            if (score > best_score) {
                best_score = score;
                best = .{ .column = @intCast(column), .direction = choice.direction };
            }
        }
        self.partial_scanned_entries += reduced_cost.len;
        self.partial_candidate_count = count;
        self.partial_cached_searches = 0;
        self.partial_full_scans += 1;
        return best;
    }

    const PrimalCandidate = struct { violation: f64, direction: f64 };

    fn primalCandidate(value: f64, column_status: basis.BasisStatus, tolerance: f64) ?PrimalCandidate {
        return switch (column_status) {
            .at_lower => if (-value > tolerance) .{ .violation = -value, .direction = 1 } else null,
            .at_upper => if (value > tolerance) .{ .violation = value, .direction = -1 } else null,
            .free, .superbasic => if (@abs(value) > tolerance) .{ .violation = @abs(value), .direction = if (value < 0) 1 else -1 } else null,
            .basic, .fixed => null,
        };
    }

    /// Weighted pricing with a deterministic virtual cost perturbation used
    /// solely as a tie-break. Candidate eligibility and reduced costs remain
    /// in the original coordinates.
    pub fn choosePrimalEnteringPerturbed(
        self: *Pricing,
        reduced_cost: []const f64,
        status: []const basis.BasisStatus,
        edge_weight: []const f64,
        column_rank: []const f64,
        taboo_until: []const usize,
        current_iteration: usize,
        tolerance: f64,
    ) ?EnteringChoice {
        if (status.len < reduced_cost.len or column_rank.len < reduced_cost.len) return null;
        if (edge_weight.len != 0 and edge_weight.len < reduced_cost.len) return null;
        if (taboo_until.len != 0 and taboo_until.len < reduced_cost.len) return null;
        // Degeneracy perturbation requires a globally stable rank tie-break;
        // keep its validated full scan until segmented candidate caching also
        // stores rank and taboo generations.
        self.iterations += 1;
        var best: ?EnteringChoice = null;
        var best_score = tolerance;
        var best_rank = std.math.inf(f64);
        for (reduced_cost, status[0..reduced_cost.len], column_rank[0..reduced_cost.len], 0..) |value, column_status, rank, column| {
            if (taboo_until.len != 0 and taboo_until[column] > current_iteration) continue;
            const violation = switch (column_status) {
                .at_lower => -value,
                .at_upper => value,
                .free, .superbasic => @abs(value),
                .basic, .fixed => continue,
            };
            const weight = if (edge_weight.len == 0) 1.0 else edge_weight[column];
            const score = self.normalizedScore(violation, weight);
            if (score > best_score + tolerance or
                (score > tolerance and @abs(score - best_score) <= tolerance and rank < best_rank))
            {
                best_score = score;
                best_rank = rank;
                const direction: f64 = switch (column_status) {
                    .at_upper => -1.0,
                    .free, .superbasic => if (value < 0) 1.0 else -1.0,
                    else => 1.0,
                };
                best = .{ .column = @intCast(column), .direction = direction };
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

test "Bland pricing selects the first improving nonbasic column" {
    var pricing = Pricing{ .rule = .devex };
    const choice = pricing.choosePrimalEnteringBland(
        &[_]f64{ -2, -100, -3 },
        &[_]basis.BasisStatus{ .basic, .at_lower, .at_lower },
        1e-9,
    ).?;
    try std.testing.expectEqual(@as(u32, 1), choice.column);
    try std.testing.expectEqual(@as(f64, 1), choice.direction);
}

test "multiple pricing reuses one validated candidate pool then refills" {
    var pricing = Pricing{ .rule = .partial, .partial_refill_interval = 1 };
    const costs = [_]f64{ -2, -5, -9, 1 };
    var statuses: [4]basis.BasisStatus = @splat(.at_lower);
    var candidates: [4]u32 = undefined;

    const first = pricing.choosePrimalEnteringMultiple(&costs, &statuses, &.{}, &candidates, 1e-9).?;
    statuses[first.column] = .basic;
    const second = pricing.choosePrimalEnteringMultiple(&costs, &statuses, &.{}, &candidates, 1e-9).?;
    statuses[second.column] = .basic;
    const third = pricing.choosePrimalEnteringMultiple(&costs, &statuses, &.{}, &candidates, 1e-9).?;

    try std.testing.expectEqual(@as(u32, 2), first.column);
    try std.testing.expectEqual(@as(u32, 1), second.column);
    try std.testing.expectEqual(@as(u32, 0), third.column);
    try std.testing.expectEqual(@as(usize, 2), pricing.partial_full_scans);
    try std.testing.expectEqual(@as(usize, 11), pricing.partial_scanned_entries);
}
