//! Exact bound transfer for rows containing one structural coefficient.

const std = @import("std");

pub const RowSide = enum {
    lower,
    upper,
};

pub const Reduction = struct {
    lower: f64,
    upper: f64,
    lower_source: ?RowSide,
    upper_source: ?RowSide,
};

pub const Action = union(enum) {
    keep,
    infeasible,
    reduce: Reduction,
};

/// Transfer a singleton row `row_lower <= coefficient * x <= row_upper`
/// into the bounds of `x`. A feasible singleton row is then redundant and can
/// be removed; source sides are recorded only for bounds made strictly tighter.
pub fn classify(
    coefficient: f64,
    row_lower: f64,
    row_upper: f64,
    col_lower: f64,
    col_upper: f64,
) Action {
    if (coefficient == 0.0 or !std.math.isFinite(coefficient)) return .keep;

    const implied_lower = if (coefficient > 0.0) row_lower / coefficient else row_upper / coefficient;
    const implied_upper = if (coefficient > 0.0) row_upper / coefficient else row_lower / coefficient;
    const lower = @max(col_lower, implied_lower);
    const upper = @min(col_upper, implied_upper);
    if (lower > upper) return .infeasible;
    return .{ .reduce = .{
        .lower = lower,
        .upper = upper,
        .lower_source = if (implied_lower > col_lower)
            (if (coefficient > 0.0) .lower else .upper)
        else
            null,
        .upper_source = if (implied_upper < col_upper)
            (if (coefficient > 0.0) .upper else .lower)
        else
            null,
    } };
}

test "singleton row tightens positive and negative coefficients" {
    const positive = classify(2.0, 4.0, 10.0, 0.0, 10.0).reduce;
    try std.testing.expectEqual(@as(f64, 2.0), positive.lower);
    try std.testing.expectEqual(@as(f64, 5.0), positive.upper);
    try std.testing.expectEqual(RowSide.lower, positive.lower_source.?);
    try std.testing.expectEqual(RowSide.upper, positive.upper_source.?);

    const negative = classify(-2.0, -10.0, -4.0, 0.0, 10.0).reduce;
    try std.testing.expectEqual(@as(f64, 2.0), negative.lower);
    try std.testing.expectEqual(@as(f64, 5.0), negative.upper);
    try std.testing.expectEqual(RowSide.upper, negative.lower_source.?);
    try std.testing.expectEqual(RowSide.lower, negative.upper_source.?);
}

test "singleton row distinguishes redundancy and infeasibility" {
    const redundant = classify(1.0, 0.0, 10.0, 2.0, 8.0).reduce;
    try std.testing.expectEqual(@as(?RowSide, null), redundant.lower_source);
    try std.testing.expectEqual(@as(?RowSide, null), redundant.upper_source);
    try std.testing.expect(classify(1.0, 4.0, 5.0, 0.0, 3.0) == .infeasible);
}
