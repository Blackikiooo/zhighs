//! Incumbent-cutoff domain tightening from an optimal LP dual solution.
//!
//! Reduced costs use the solver's normalized sign convention for both objective
//! senses. Tightenings are node-local and appended to the persistent domain
//! arena so descendants inherit them without copying complete bound arrays.

const std = @import("std");
const solver = @import("solver");
const types = @import("types.zig");
const domain = @import("domain.zig");

pub const Result = struct {
    feasible: bool,
    change_head: ?domain.ChangeId,
    tightenings: usize,
};

pub const ImpliedBounds = struct {
    lower: ?f64 = null,
    upper: ?f64 = null,
};

pub fn cutoffGap(sense: solver.LpObjectiveSense, lp_objective: f64, incumbent_objective: f64) f64 {
    return switch (sense) {
        .minimize => incumbent_objective - lp_objective,
        .maximize => lp_objective - incumbent_objective,
    };
}

pub fn impliedBounds(kind: types.Integrality, lower: f64, upper: f64, cost: f64, gap: f64, tolerance: f64) ImpliedBounds {
    if ((kind != .binary and kind != .integer) or lower >= upper - tolerance or
        !std.math.isFinite(cost) or @abs(cost) <= 10.0 * tolerance or gap < 0.0)
        return .{};
    const maximum_degradation = cost * (upper - lower);
    if (maximum_degradation > gap and std.math.isFinite(lower)) {
        const candidate = @floor(gap / cost + lower + tolerance);
        if (std.math.isFinite(candidate) and candidate < upper - tolerance) return .{ .upper = @max(lower, candidate) };
    } else if (maximum_degradation < -gap and std.math.isFinite(upper)) {
        const candidate = @ceil(gap / cost + upper - tolerance);
        if (std.math.isFinite(candidate) and candidate > lower + tolerance) return .{ .lower = @min(upper, candidate) };
    }
    return .{};
}

pub fn tighten(
    allocator: std.mem.Allocator,
    sense: solver.LpObjectiveSense,
    integrality: []const types.Integrality,
    lower: []f64,
    upper: []f64,
    reduced_cost: []const f64,
    lp_objective: f64,
    incumbent_objective: f64,
    arena: *domain.DomainChangeArena,
    initial_head: ?domain.ChangeId,
    tolerance: f64,
) std.mem.Allocator.Error!Result {
    std.debug.assert(integrality.len == lower.len);
    std.debug.assert(upper.len == lower.len);
    std.debug.assert(reduced_cost.len == lower.len);
    const gap = cutoffGap(sense, lp_objective, incumbent_objective);
    if (!std.math.isFinite(gap) or gap < -tolerance)
        return .{ .feasible = gap >= -tolerance, .change_head = initial_head, .tightenings = 0 };

    const nonnegative_gap = @max(0.0, gap);
    var head = initial_head;
    var tightenings: usize = 0;
    for (integrality, lower, upper, reduced_cost, 0..) |kind, *lo, *up, cost, column| {
        const implied = impliedBounds(kind, lo.*, up.*, cost, nonnegative_gap, tolerance);
        if (implied.upper) |candidate| {
            up.* = candidate;
            head = try arena.appendWithReason(allocator, head, .{ .column = column, .upper = up.* }, .reduced_cost);
            tightenings += 1;
        } else if (implied.lower) |candidate| {
            lo.* = candidate;
            head = try arena.appendWithReason(allocator, head, .{ .column = column, .lower = lo.* }, .reduced_cost);
            tightenings += 1;
        }
    }
    return .{ .feasible = true, .change_head = head, .tightenings = tightenings };
}

test "minimization reduced cost tightens an integer upper bound" {
    var lower = [_]f64{0.0};
    var upper = [_]f64{10.0};
    var arena = domain.DomainChangeArena{};
    defer arena.deinit(std.testing.allocator);
    const result = try tighten(std.testing.allocator, .minimize, &.{.integer}, &lower, &upper, &.{2.0}, 0.0, 5.0, &arena, null, 1e-6);
    try std.testing.expect(result.feasible);
    try std.testing.expectEqual(@as(usize, 1), result.tightenings);
    try std.testing.expectEqual(@as(f64, 2.0), upper[0]);
}

test "maximization reduced cost tightens an integer lower bound" {
    var lower = [_]f64{0.0};
    var upper = [_]f64{10.0};
    var arena = domain.DomainChangeArena{};
    defer arena.deinit(std.testing.allocator);
    const result = try tighten(std.testing.allocator, .maximize, &.{.integer}, &lower, &upper, &.{-3.0}, 30.0, 23.0, &arena, null, 1e-6);
    try std.testing.expect(result.feasible);
    try std.testing.expectEqual(@as(usize, 1), result.tightenings);
    try std.testing.expectEqual(@as(f64, 8.0), lower[0]);
}

test "semi variables are skipped until disjunctive reduced-cost fixing exists" {
    var lower = [_]f64{0.0};
    var upper = [_]f64{10.0};
    var arena = domain.DomainChangeArena{};
    defer arena.deinit(std.testing.allocator);
    const result = try tighten(std.testing.allocator, .minimize, &.{.semi_integer}, &lower, &upper, &.{2.0}, 0.0, 1.0, &arena, null, 1e-6);
    try std.testing.expectEqual(@as(usize, 0), result.tightenings);
    try std.testing.expectEqual(@as(f64, 10.0), upper[0]);
}
