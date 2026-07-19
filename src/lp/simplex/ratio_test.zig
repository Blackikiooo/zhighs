//! Primal/dual ratio-test policies.

const std = @import("std");
const basis = @import("basis.zig");

pub const RatioRule = enum { standard, harris_two_pass, bound_flipping };
pub const LeavingChoice = struct { row: ?u32 = null, step: f64 = 0.0 };
pub const DualEnteringChoice = struct {
    column: ?u32 = null,
    direction: f64 = 0.0,
    theta: f64 = 0.0,
    flip_count: usize = 0,
};

pub const RatioTest = struct {
    rule: RatioRule = .bound_flipping,
    tolerance: f64 = 1e-9,

    pub fn chooseLeaving(self: *const RatioTest, direction: []const f64, rhs: []const f64) LeavingChoice {
        if (direction.len != rhs.len) return .{};
        return switch (self.rule) {
            .standard => self.chooseStandard(direction, rhs),
            .harris_two_pass, .bound_flipping => self.chooseHarris(direction, rhs),
        };
    }

    /// Harris selection with deterministic bounded primal perturbations.
    /// Positive ranks provide a tiny virtual margin to zero-valued basics, so
    /// a degenerate face produces real lexicographic progress. The caller
    /// limits ranks below primal feasibility tolerance and performs a fresh
    /// unperturbed cleanup before publishing a status.
    pub fn chooseLeavingPerturbed(
        self: *const RatioTest,
        direction: []const f64,
        rhs: []const f64,
        row_rank: []const f64,
        degeneracy_tolerance: f64,
    ) LeavingChoice {
        if (direction.len != rhs.len or row_rank.len < direction.len) return .{};
        var relaxed_step = std.math.inf(f64);
        var minimum_step = std.math.inf(f64);
        for (direction, rhs, row_rank[0..direction.len]) |coefficient, value, rank| {
            if (coefficient <= self.tolerance) continue;
            const virtual_margin = @max(value + rank, 0.0);
            minimum_step = @min(minimum_step, virtual_margin / coefficient);
            relaxed_step = @min(relaxed_step, (virtual_margin + self.tolerance) / coefficient);
        }
        if (!std.math.isFinite(relaxed_step)) return .{ .step = relaxed_step };
        if (minimum_step > degeneracy_tolerance) return self.chooseHarris(direction, rhs);

        var maximum_pivot: f64 = 0.0;
        for (direction, rhs, row_rank[0..direction.len]) |coefficient, value, rank| {
            if (coefficient <= self.tolerance) continue;
            const virtual_step = @max(value + rank, 0.0) / coefficient;
            if (virtual_step <= relaxed_step and virtual_step <= minimum_step + degeneracy_tolerance)
                maximum_pivot = @max(maximum_pivot, coefficient);
        }
        var choice = LeavingChoice{ .step = std.math.inf(f64) };
        var best_rank = std.math.inf(f64);
        var best_pivot: f64 = 0.0;
        for (direction, rhs, row_rank[0..direction.len], 0..) |coefficient, value, rank, row| {
            if (coefficient <= self.tolerance) continue;
            const virtual_step = @max(value + rank, 0.0) / coefficient;
            if (virtual_step > relaxed_step or virtual_step > minimum_step + degeneracy_tolerance) continue;
            // Threshold stability: perturb only candidates retaining at least
            // half the strongest pivot in the degenerate Harris set.
            if (coefficient < maximum_pivot * 0.5) continue;
            if (rank < best_rank or (rank == best_rank and coefficient > best_pivot)) {
                best_rank = rank;
                best_pivot = coefficient;
                choice = .{ .row = @intCast(row), .step = virtual_step };
            }
        }
        return choice;
    }

    fn chooseStandard(self: *const RatioTest, direction: []const f64, rhs: []const f64) LeavingChoice {
        var choice = LeavingChoice{ .step = std.math.inf(f64) };
        for (direction, rhs, 0..) |coefficient, value, i| {
            if (coefficient > self.tolerance) {
                const step = value / coefficient;
                if (step >= -self.tolerance and step < choice.step) {
                    choice = .{ .row = @intCast(i), .step = @max(step, 0.0) };
                }
            }
        }
        return choice;
    }

    fn chooseHarris(self: *const RatioTest, direction: []const f64, rhs: []const f64) LeavingChoice {
        var relaxed_step = std.math.inf(f64);
        for (direction, rhs) |coefficient, value| {
            if (coefficient > self.tolerance) {
                relaxed_step = @min(relaxed_step, (@max(value, 0.0) + self.tolerance) / coefficient);
            }
        }
        if (!std.math.isFinite(relaxed_step)) return .{ .step = relaxed_step };

        var choice = LeavingChoice{ .step = std.math.inf(f64) };
        var best_pivot: f64 = 0.0;
        for (direction, rhs, 0..) |coefficient, value, i| {
            if (coefficient <= self.tolerance) continue;
            const exact_step = @max(value / coefficient, 0.0);
            if (exact_step <= relaxed_step and coefficient > best_pivot) {
                best_pivot = coefficient;
                choice = .{ .row = @intCast(i), .step = exact_step };
            }
        }
        return choice;
    }

    /// Dual ratio test with boxed-variable bound flipping. Candidate ratios
    /// and directions are written into caller-owned SoA workspaces; candidate
    /// indices are sorted in place without allocation.
    pub fn chooseDualEntering(
        self: *const RatioTest,
        tableau: []const f64,
        reduced_cost: []const f64,
        status: []const basis.BasisStatus,
        lower: []const f64,
        upper: []const f64,
        primal: []const f64,
        leaving_bound: basis.BasisStatus,
        primal_infeasibility: f64,
        ratio_work: []f64,
        direction_work: []f64,
        candidate_work: []u32,
    ) DualEnteringChoice {
        const count = tableau.len;
        if (reduced_cost.len < count or status.len < count or lower.len < count or upper.len < count or
            primal.len < count or ratio_work.len < count or direction_work.len < count or candidate_work.len < count)
            return .{};

        var candidate_count: usize = 0;
        for (0..count) |column| {
            ratio_work[column] = std.math.inf(f64);
            direction_work[column] = 0.0;
            const alpha = tableau[column];
            if (@abs(alpha) <= self.tolerance) continue;

            const direction: f64 = switch (status[column]) {
                .at_lower => 1.0,
                .at_upper => -1.0,
                .free, .superbasic => if (leaving_bound == .at_lower)
                    (if (alpha < 0.0) 1.0 else -1.0)
                else
                    (if (alpha > 0.0) 1.0 else -1.0),
                .basic, .fixed => continue,
            };
            const signed_alpha = alpha * direction;
            const eligible = if (leaving_bound == .at_lower)
                signed_alpha < -self.tolerance
            else
                signed_alpha > self.tolerance;
            if (!eligible) continue;

            const rc = reduced_cost[column];
            const ratio = switch (status[column]) {
                .at_lower => @max(rc, 0.0) / @abs(alpha),
                .at_upper => @max(-rc, 0.0) / @abs(alpha),
                .free, .superbasic => @abs(rc) / @abs(alpha),
                else => unreachable,
            };
            ratio_work[column] = ratio;
            direction_work[column] = direction;
            candidate_work[candidate_count] = @intCast(column);
            candidate_count += 1;
        }
        if (candidate_count == 0) return .{};

        std.sort.pdq(u32, candidate_work[0..candidate_count], ratio_work, struct {
            fn lessThan(ratios: []f64, lhs: u32, rhs: u32) bool {
                return ratios[lhs] < ratios[rhs] or (ratios[lhs] == ratios[rhs] and lhs < rhs);
            }
        }.lessThan);

        var corrected: f64 = 0.0;
        var flip_count: usize = 0;
        for (candidate_work[0..candidate_count]) |column_u32| {
            const column: usize = @intCast(column_u32);
            const width = upper[column] - lower[column];
            const boxed = std.math.isFinite(width) and width > self.tolerance and
                (status[column] == .at_lower or status[column] == .at_upper);
            const capacity = if (boxed) @abs(tableau[column]) * width else 0.0;
            if (self.rule == .bound_flipping and boxed and corrected + capacity < primal_infeasibility - self.tolerance) {
                candidate_work[flip_count] = column_u32;
                flip_count += 1;
                corrected += capacity;
                continue;
            }
            const ratio = ratio_work[column];
            return .{
                .column = column_u32,
                .direction = direction_work[column],
                .theta = if (leaving_bound == .at_lower) -ratio else ratio,
                .flip_count = flip_count,
            };
        }
        return .{ .flip_count = flip_count };
    }
};

test {
    std.testing.refAllDecls(@This());
}

test "Harris ratio test prefers stable pivot within relaxed bound" {
    const test_rule = RatioTest{ .rule = .harris_two_pass, .tolerance = 1e-6 };
    const choice = test_rule.chooseLeaving(&[_]f64{ 1e-4, 1.0 }, &[_]f64{ 1e-4, 1.000001 });
    try std.testing.expectEqual(@as(?u32, 1), choice.row);
    try std.testing.expectApproxEqAbs(@as(f64, 1.000001), choice.step, 1e-12);
}

test "perturbed Harris gives a bounded positive step on a degenerate tie" {
    const test_rule = RatioTest{ .rule = .harris_two_pass, .tolerance = 1e-9 };
    const choice = test_rule.chooseLeavingPerturbed(&.{ 1, 2 }, &.{ 0, 0 }, &.{ 0.2, 0.1 }, 1e-7);
    try std.testing.expectEqual(@as(?u32, 1), choice.row);
    try std.testing.expect(choice.step > 0.0 and choice.step <= 1e-7);
}

test "dual bound-flipping ratio test records boxed breakpoints" {
    const test_rule = RatioTest{ .rule = .bound_flipping, .tolerance = 1e-9 };
    var ratios: [2]f64 = undefined;
    var directions: [2]f64 = undefined;
    var candidates: [2]u32 = undefined;
    const choice = test_rule.chooseDualEntering(
        &[_]f64{ -1, -1 },
        &[_]f64{ 0, 2 },
        &[_]basis.BasisStatus{ .at_lower, .at_lower },
        &[_]f64{ 0, 0 },
        &[_]f64{ 1, std.math.inf(f64) },
        &[_]f64{ 0, 0 },
        .at_lower,
        2,
        &ratios,
        &directions,
        &candidates,
    );
    try std.testing.expectEqual(@as(usize, 1), choice.flip_count);
    try std.testing.expectEqual(@as(?u32, 1), choice.column);
}

test "dual ratio ties use the lowest column index" {
    const test_rule = RatioTest{ .rule = .standard, .tolerance = 1e-9 };
    var ratios: [2]f64 = undefined;
    var directions: [2]f64 = undefined;
    var candidates: [2]u32 = undefined;
    const choice = test_rule.chooseDualEntering(
        &[_]f64{ -1, -2 },
        &[_]f64{ 1, 2 },
        &[_]basis.BasisStatus{ .at_lower, .at_lower },
        &[_]f64{ 0, 0 },
        &[_]f64{ std.math.inf(f64), std.math.inf(f64) },
        &[_]f64{ 0, 0 },
        .at_lower,
        1,
        &ratios,
        &directions,
        &candidates,
    );
    try std.testing.expectEqual(@as(?u32, 0), choice.column);
    try std.testing.expectEqual(@as(usize, 0), choice.flip_count);
}
