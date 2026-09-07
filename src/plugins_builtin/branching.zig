//! Built-in branching rules.
//!
//! Algorithm implementations live here; registration, ordering, timing, and
//! lifecycle remain the responsibility of `framework`.

const std = @import("std");
const plugin = @import("plugin");

/// Classic most-fractional-variable branching. Among integer-like columns whose
/// LP value is fractional, selects the one with maximum fractionality.
/// Semi-continuous / semi-integer columns that neither sit at an allowed zero
/// nor inside their active range are reported as immediate semi-activation
/// violations, mirroring the pre-plugin engine's behaviour.
pub const MostFractional = struct {
    pub const name = "most-fractional";

    pub fn select(_: *MostFractional, ctx: plugin.BranchingContext) ?plugin.BranchingCandidate {
        var best: ?plugin.BranchingCandidate = null;
        var score: f64 = 0;
        for (ctx.problem.integrality, ctx.lp_primal, 0..) |kind, value, column| {
            if (kind == .continuous) continue;
            if (kind == .semi_continuous or kind == .semi_integer) {
                const zero_allowed = ctx.node_lower[column] <= ctx.tolerance and ctx.node_upper[column] >= -ctx.tolerance;
                const active = value >= ctx.active_lower[column] - ctx.tolerance and
                    value <= ctx.active_upper[column] + ctx.tolerance;
                if (zero_allowed and @abs(value) <= ctx.tolerance) continue;
                if (!active) return .{ .column = column, .kind = .semi_activation };
                if (kind == .semi_continuous) continue;
            }
            const fractionality = @abs(value - @round(value));
            if (fractionality > ctx.tolerance and fractionality > score) {
                best = .{ .column = column, .kind = .integer };
                score = fractionality;
            }
        }
        return best;
    }
};

const matrix = @import("matrix");

fn contextFor(
    comptime integrality: []const plugin.Integrality,
    primal: []const f64,
    node_lower: []const f64,
    node_upper: []const f64,
    active_lower: []const f64,
    active_upper: []const f64,
) plugin.BranchingContext {
    const num_cols = integrality.len;
    return .{
        .problem = .{
            .relaxation = .{
                .num_rows = 0,
                .num_cols = num_cols,
                .col_cost = &.{},
                .col_lower = &.{},
                .col_upper = &.{},
                .row_lower = &.{},
                .row_upper = &.{},
                .matrix = matrix.CscView.initAssumeValid(0, num_cols, &([_]usize{0} ** (num_cols + 1)), &.{}, &.{}),
                .objective_sense = .minimize,
                .objective_offset = 0.0,
            },
            .integrality = integrality,
        },
        .node_lower = node_lower,
        .node_upper = node_upper,
        .lp_primal = primal,
        .active_lower = active_lower,
        .active_upper = active_upper,
        .tolerance = 1e-6,
    };
}

test "most-fractional selects the column with maximum fractionality" {
    var rule = MostFractional{};
    const ctx = contextFor(
        &.{ .integer, .integer, .integer },
        &.{ 0.5, 0.25, 0.9 },
        &.{ 0, 0, 0 },
        &.{ 1, 1, 1 },
        &.{ 0, 0, 0 },
        &.{ 1, 1, 1 },
    );
    const candidate = rule.select(ctx).?;
    try std.testing.expectEqual(@as(usize, 0), candidate.column);
    try std.testing.expectEqual(plugin.BranchingKind.integer, candidate.kind);
}

test "most-fractional ignores integral and continuous columns" {
    var rule = MostFractional{};
    const ctx = contextFor(
        &.{ .continuous, .integer, .binary },
        &.{ 0.3, 2.0, 1.0 },
        &.{ 0, 0, 0 },
        &.{ 1, 10, 1 },
        &.{ 0, 0, 0 },
        &.{ 1, 10, 1 },
    );
    try std.testing.expectEqual(@as(?plugin.BranchingCandidate, null), rule.select(ctx));
}

test "most-fractional reports semi-activation violations before fractionality" {
    var rule = MostFractional{};
    const ctx = contextFor(
        &.{ .semi_continuous, .integer },
        &.{ 0.5, 0.5 },
        &.{ -1, 0 },
        &.{ 5, 1 },
        &.{ 1, 0 },
        &.{ 5, 1 },
    );
    const candidate = rule.select(ctx).?;
    try std.testing.expectEqual(@as(usize, 0), candidate.column);
    try std.testing.expectEqual(plugin.BranchingKind.semi_activation, candidate.kind);
}

test "most-fractional accepts a semi column resting at an allowed zero" {
    var rule = MostFractional{};
    const ctx = contextFor(
        &.{ .semi_continuous, .integer },
        &.{ 0.0, 0.5 },
        &.{ -1, 0 },
        &.{ 5, 1 },
        &.{ 1, 0 },
        &.{ 5, 1 },
    );
    const candidate = rule.select(ctx).?;
    try std.testing.expectEqual(@as(usize, 1), candidate.column);
    try std.testing.expectEqual(plugin.BranchingKind.integer, candidate.kind);
}
