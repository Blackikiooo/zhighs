//! Built-in primal heuristics.
//!
//! Algorithm implementations live here; registration, ordering, timing, and
//! lifecycle remain the responsibility of `framework`.

const std = @import("std");
const plugin = @import("plugin");

/// Simple LP-solution rounding. Snaps each integer-like column to its nearest
/// integral value inside the active column bounds and resolves semi-continuous
/// / semi-integer columns to either their zero branch or their active value,
/// whichever is closer to the LP point. The engine validates the candidate, so
/// an infeasible rounding is simply discarded.
pub const Rounding = struct {
    pub const name = "rounding";

    pub fn trySolution(_: *Rounding, ctx: plugin.HeuristicContext) std.mem.Allocator.Error!?[]f64 {
        const rounded = try ctx.allocator.dupe(f64, ctx.lp_primal);
        errdefer ctx.allocator.free(rounded);
        for (ctx.problem.integrality, rounded, 0..) |kind, *value, column| switch (kind) {
            .continuous => {},
            .integer => value.* = std.math.clamp(@round(value.*), ctx.active_lower[column], ctx.active_upper[column]),
            .binary => value.* = std.math.clamp(@round(value.*), @max(0.0, ctx.problem.relaxation.col_lower[column]), @min(1.0, ctx.problem.relaxation.col_upper[column])),
            .semi_continuous, .semi_integer => {
                const lower = ctx.active_lower[column];
                const upper = ctx.active_upper[column];
                if (lower > upper) {
                    value.* = 0;
                    continue;
                }
                const active = if (kind == .semi_integer)
                    std.math.clamp(@round(value.*), lower, upper)
                else
                    std.math.clamp(value.*, lower, upper);
                value.* = if (@abs(value.*) <= @abs(value.* - active)) 0 else active;
            },
        };
        return rounded;
    }
};

const matrix = @import("matrix");

fn contextFor(
    comptime integrality: []const plugin.Integrality,
    primal: []const f64,
    col_lower: []const f64,
    col_upper: []const f64,
    active_lower: []const f64,
    active_upper: []const f64,
) plugin.HeuristicContext {
    const num_cols = integrality.len;
    return .{
        .allocator = std.testing.allocator,
        .problem = .{
            .relaxation = .{
                .num_rows = 0,
                .num_cols = num_cols,
                .col_cost = &.{},
                .col_lower = col_lower,
                .col_upper = col_upper,
                .row_lower = &.{},
                .row_upper = &.{},
                .matrix = matrix.CscView.initAssumeValid(0, num_cols, &([_]usize{0} ** (num_cols + 1)), &.{}, &.{}),
                .objective_sense = .minimize,
                .objective_offset = 0.0,
            },
            .integrality = integrality,
        },
        .node_lower = &.{},
        .node_upper = &.{},
        .lp_primal = primal,
        .active_lower = active_lower,
        .active_upper = active_upper,
        .tolerance = 1e-6,
    };
}

test "rounding snaps integer columns to the nearest active integral value" {
    var heuristic = Rounding{};
    const ctx = contextFor(
        &.{ .integer, .integer },
        &.{ 0.4, 1.6 },
        &.{ 0.0, 0.0 },
        &.{ 5.0, 5.0 },
        &.{ 0.0, 0.0 },
        &.{ 5.0, 5.0 },
    );
    const candidate = try heuristic.trySolution(ctx);
    defer std.testing.allocator.free(candidate.?);
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 2.0 }, candidate.?);
}

test "rounding clamps to the active bounds when the rounded value escapes" {
    var heuristic = Rounding{};
    const ctx = contextFor(
        &.{.integer},
        &.{4.7},
        &.{0.0},
        &.{3.0},
        &.{0.0},
        &.{3.0},
    );
    const candidate = try heuristic.trySolution(ctx);
    defer std.testing.allocator.free(candidate.?);
    try std.testing.expectEqualSlices(f64, &.{3.0}, candidate.?);
}

test "rounding resolves semi-continuous columns to zero or the active value" {
    var heuristic = Rounding{};
    const ctx = contextFor(
        &.{ .semi_continuous, .semi_integer },
        &.{ 0.4, 4.6 },
        &.{ -1.0, -1.0 },
        &.{ 6.0, 6.0 },
        &.{ 2.0, 2.0 },
        &.{ 6.0, 6.0 },
    );
    const candidate = try heuristic.trySolution(ctx);
    defer std.testing.allocator.free(candidate.?);
    // 0.4 is closer to zero; 4.6 rounds to 5 inside the active range.
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 5.0 }, candidate.?);
}

test "rounding forces a semi column to zero when its active range is empty" {
    var heuristic = Rounding{};
    const ctx = contextFor(
        &.{.semi_continuous},
        &.{3.0},
        &.{ -1.0, },
        &.{ 6.0, },
        &.{ 6.0, },
        &.{ 2.0, },
    );
    const candidate = try heuristic.trySolution(ctx);
    defer std.testing.allocator.free(candidate.?);
    try std.testing.expectEqualSlices(f64, &.{0.0}, candidate.?);
}
