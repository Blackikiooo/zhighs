//! Entering-variable pricing policies.

const std = @import("std");
const basis = @import("basis.zig");

/// Available pricing rules. `dantzig` is the textbook rule; `devex` and
/// `steepest_edge` use weighted norms; `partial`/`hyper_sparse` use segmented
/// candidate pools to amortize the scan cost.
pub const PricingRule = enum { dantzig, devex, steepest_edge, partial, hyper_sparse };

/// Result of a primal pricing scan: which column enters and the sign of the
/// primal step that improves the objective.
pub const EnteringChoice = struct { column: u32, direction: f64 };

/// Result of a dual pricing scan: which basic row leaves, which bound it
/// violates, and by how much.
pub const DualLeavingChoice = struct {
    row: u32,
    bound: basis.BasisStatus,
    violation: f64,
};

pub const Pricing = struct {
    rule: PricingRule = .devex,
    devex_reset_period: usize = 100, // Devex weights are reset after this many iterations
    iterations: usize = 0,
    /// Persistent segmented-pricing cursor. A segment is selected once per
    /// pricing operation; only an all-segment miss can certify optimality.
    partial_candidate_count: usize = 0, // Live entries in the candidate pool
    partial_cached_searches: usize = 0, // Scans served from the pool since last refill
    partial_refill_interval: usize = 1, // Max cached searches before a forced refill
    partial_searches: usize = 0, // Total pricing calls (cached + full)
    partial_scanned_entries: usize = 0, // Total nonzeros inspected across all scans
    partial_full_scans: usize = 0, // Number of full refill scans performed

    /// Clear all partial-pricing statistics (called between solves).
    pub fn resetPartial(self: *Pricing) void {
        self.partial_candidate_count = 0;
        self.partial_cached_searches = 0;
        self.partial_searches = 0;
        self.partial_scanned_entries = 0;
        self.partial_full_scans = 0;
    }

    /// Apply the rule-specific normalization: Dantzig uses the raw violation,
    /// weighted rules divide by sqrt(weight) to approximate the steepest edge.
    fn normalizedScore(self: Pricing, violation: f64, weight: f64) f64 {
        return switch (self.rule) {
            .dantzig, .partial => violation,
            .devex, .steepest_edge, .hyper_sparse => violation / @sqrt(@max(weight, 1.0)),
        };
    }

    /// Simplest Dantzig scan: pick the column with the largest |reduced cost|.
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

    /// Same as `choosePrimalEntering` but also returns the primal direction.
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

        // Try to serve the request from the cached pool first.
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

        // Pool exhausted (or first call): perform a full scan and refill.
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

    /// Intermediate result for one primal candidate: its violation magnitude
    /// and the direction that improves the objective.
    const PrimalCandidate = struct { violation: f64, direction: f64 };

    /// Compute the violation and direction for a single (value, status) pair.
    /// Returns null when the column is not eligible to enter.
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
        column_rank: []const f64, // Deterministic tie-break rank per column
        taboo_until: []const usize, // Iteration until which a column is forbidden
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
            // Accept a new candidate when it strictly improves the score, or
            // when it ties within tolerance and has a lexicographically smaller rank.
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

    /// Bland dual pricing fallback: selects the eligible basic variable with
    /// the smallest basic index. Guarantees deterministic termination during
    /// anti-cycling.
    fn chooseDualLeavingBland(
        self: *Pricing,
        basic_value: []const f64,
        basic_lower: []const f64,
        basic_upper: []const f64,
        basic_index: []const u32,
        primal_tolerance: f64,
    ) ?DualLeavingChoice {
        _ = self;
        var best: ?DualLeavingChoice = null;
        var best_basic_column: u32 = std.math.maxInt(u32);
        for (basic_value, basic_lower, basic_upper, basic_index, 0..) |value, lower, upper, basic_col, row| {
            const choice: ?DualLeavingChoice = if (value < lower - primal_tolerance)
                .{ .row = @intCast(row), .bound = .at_lower, .violation = lower - value }
            else if (value > upper + primal_tolerance)
                .{ .row = @intCast(row), .bound = .at_upper, .violation = value - upper }
            else
                null;
            if (choice) |candidate| {
                if (basic_col < best_basic_column) {
                    best_basic_column = basic_col;
                    best = candidate;
                }
            }
        }
        return best;
    }

    /// Unified dual leaving-row dispatch. Routes to Bland, weighted, or
    /// hyper-sparse candidate selection based on engine state and the
    /// active pricing rule.
    pub fn chooseDualLeaving(self: *Pricing, engine: anytype) ?DualLeavingChoice {
        const pricing_started = engine.statisticsTimestamp();
        defer engine.recordRowPricingElapsed(pricing_started);
        const bs = if (engine.basis) |*value| value else return null;
        if (self.rule == .hyper_sparse and engine.dual_hyper_sparse_active) {
            engine.stats.hyper_pricing_dispatches += 1;
        } else {
            engine.stats.dense_pricing_dispatches += 1;
        }
        if (engine.numerical.anti_cycling_active) {
            return self.chooseDualLeavingBland(
                bs.basic_value,
                bs.basic_lower,
                bs.basic_upper,
                bs.basic_index,
                engine.numerical.primal_tolerance,
            );
        }
        if (self.rule != .hyper_sparse or !engine.dual_hyper_sparse_active) {
            return self.chooseDualLeavingWeighted(
                bs.basic_value,
                bs.basic_lower,
                bs.basic_upper,
                bs.row_edge_weight,
                engine.numerical.primal_tolerance,
            );
        }
        var best = engine.bestDualCandidate();
        if (best == null or engine.dualCandidateScore(best.?.row) + engine.numerical.primal_tolerance < engine.dual_candidate_cutoff) {
            engine.rebuildDualCandidateList();
            best = engine.bestDualCandidate();
        }
        return best;
    }

    /// Dantzig-style dual pricing over the dense basic-value SoA. The most
    /// primal-infeasible row leaves; no temporary candidate objects are built.
    pub fn chooseDualLeavingDantzig(self: *Pricing, value: []const f64, lower: []const f64, upper: []const f64, tolerance: f64) ?DualLeavingChoice {
        return self.chooseDualLeavingWeighted(value, lower, upper, &.{}, tolerance);
    }

    /// Weighted dual pricing. Score is the bound violation normalized by the
    /// row's edge weight; the largest score wins.
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
