//! Row-activity implied bounds and forced variable fixing.
//!
//! For a row `l <= a_j x_j + S <= u` with `S in [S_l, S_u]`, the column
//! bounds imply
//!
//!     a_j > 0:  (l - S_u) / a_j <= x_j <= (u - S_l) / a_j
//!     a_j < 0:  (u - S_l) / a_j <= x_j <= (l - S_u) / a_j
//!
//! All arithmetic is outward rounded: division is stepped by one ulp in the
//! requested direction, which is conservative for any denominator sign. The
//! caller only commits a fixing when the implied interval collapses, so an
//! implied bound never has to be published as a looser/tighter column bound.
//! That keeps dual postsolve unchanged: a fixed column is already handled by
//! the exact column-substitution record.

const std = @import("std");
const rules = @import("redundant_row.zig");

pub const Bounds = struct {
    lower: f64 = -std.math.inf(f64),
    upper: f64 = std.math.inf(f64),
};

pub const Kind = enum {
    continuous,
    integer,
    binary,
};

/// Relative width below which a continuous implied interval is treated as a
/// single point. The value is kept well inside simplex feasibility tolerances.
pub const FIX_TOLERANCE: f64 = 1e-12;

/// Conservative interval of one row term `coefficient * x`.
pub fn termBounds(coefficient: f64, column_lower: f64, column_upper: f64) Bounds {
    if (coefficient == 0.0 or !std.math.isFinite(coefficient)) return .{};
    if (std.math.isNan(column_lower) or std.math.isNan(column_upper)) return .{ .lower = std.math.nan(f64), .upper = std.math.nan(f64) };
    const lower = if (coefficient > 0.0) column_lower else column_upper;
    const upper = if (coefficient > 0.0) column_upper else column_lower;
    return .{
        .lower = rules.outwardProduct(coefficient, lower, -std.math.inf(f64)),
        .upper = rules.outwardProduct(coefficient, upper, std.math.inf(f64)),
    };
}

/// Upper bound on the term's minimum contribution. Subtracting this from a
/// lower row-activity bound gives a rigorous lower bound on all other terms
/// without the looseness of subtracting the term's maximum.
pub fn minimumUpper(coefficient: f64, column_lower: f64, column_upper: f64) f64 {
    if (coefficient == 0.0 or !std.math.isFinite(coefficient)) return 0.0;
    if (std.math.isNan(column_lower) or std.math.isNan(column_upper)) return std.math.nan(f64);
    const value = if (coefficient > 0.0) column_lower else column_upper;
    return rules.outwardProduct(coefficient, value, std.math.inf(f64));
}

/// Lower bound on the term's maximum contribution, used symmetrically to
/// exclude one term from an upper row-activity bound.
pub fn maximumLower(coefficient: f64, column_lower: f64, column_upper: f64) f64 {
    if (coefficient == 0.0 or !std.math.isFinite(coefficient)) return 0.0;
    if (std.math.isNan(column_lower) or std.math.isNan(column_upper)) return std.math.nan(f64);
    const value = if (coefficient > 0.0) column_upper else column_lower;
    return rules.outwardProduct(coefficient, value, -std.math.inf(f64));
}

/// Bounds implied on one column by `l <= coefficient*x + others <= u`.
/// `others` must bound the sum of every other term in the row.
pub fn impliedBounds(coefficient: f64, row_lower: f64, row_upper: f64, others: Bounds) Bounds {
    if (coefficient == 0.0 or !std.math.isFinite(coefficient)) return .{};
    if (std.math.isNan(row_lower) or std.math.isNan(row_upper) or
        std.math.isNan(others.lower) or std.math.isNan(others.upper))
        return .{ .lower = std.math.nan(f64), .upper = std.math.nan(f64) };

    var result = Bounds{};
    if (coefficient > 0.0) {
        if (std.math.isFinite(row_lower) and std.math.isFinite(others.upper)) {
            const numerator = rules.outwardAdd(row_lower, -others.upper, -std.math.inf(f64));
            result.lower = @max(result.lower, divideOutward(numerator, coefficient, -std.math.inf(f64)));
        }
        if (std.math.isFinite(row_upper) and std.math.isFinite(others.lower)) {
            const numerator = rules.outwardAdd(row_upper, -others.lower, std.math.inf(f64));
            result.upper = @min(result.upper, divideOutward(numerator, coefficient, std.math.inf(f64)));
        }
    } else {
        if (std.math.isFinite(row_upper) and std.math.isFinite(others.lower)) {
            // a < 0: the lower bound comes from the row upper side. The
            // numerator must be rounded up before division by a negative.
            const numerator = rules.outwardAdd(row_upper, -others.lower, std.math.inf(f64));
            result.lower = @max(result.lower, divideOutward(numerator, coefficient, -std.math.inf(f64)));
        }
        if (std.math.isFinite(row_lower) and std.math.isFinite(others.upper)) {
            const numerator = rules.outwardAdd(row_lower, -others.upper, -std.math.inf(f64));
            result.upper = @min(result.upper, divideOutward(numerator, coefficient, std.math.inf(f64)));
        }
    }
    return result;
}

/// Return the forced value when the interval collapses to one feasible point
/// for the given domain. Integer and binary domains use `ceil`/`floor`, which
/// is valid because `Bounds` is outward rounded.
pub fn fixValue(bounds: Bounds, kind: Kind) ?f64 {
    if (std.math.isNan(bounds.lower) or std.math.isNan(bounds.upper)) return null;
    var lower = bounds.lower;
    var upper = bounds.upper;
    if (kind == .binary) {
        lower = @max(lower, 0.0);
        upper = @min(upper, 1.0);
    }
    if (lower > upper) return null;
    // An infinite side can never define a fixing; in particular the relative
    // tolerance test below would otherwise compare inf <= inf.
    if (!std.math.isFinite(lower) or !std.math.isFinite(upper)) return null;

    switch (kind) {
        .continuous => {
            if (lower == upper) return lower;
            const scale = @max(1.0, @max(@abs(lower), @abs(upper)));
            if (upper - lower <= FIX_TOLERANCE * scale) {
                const midpoint = lower + 0.5 * (upper - lower);
                return if (std.math.isFinite(midpoint)) midpoint else null;
            }
            return null;
        },
        .integer, .binary => {
            const integral_lower = @ceil(lower);
            const integral_upper = @floor(upper);
            if (!(integral_lower <= integral_upper) or integral_lower != integral_upper) return null;
            return integral_lower;
        },
    }
}

fn divideOutward(numerator: f64, denominator: f64, direction: f64) f64 {
    const quotient = numerator / denominator;
    if (!std.math.isFinite(quotient)) return quotient;
    return std.math.nextAfter(f64, quotient, direction);
}

test "forcing equality fixes every column at its extreme" {
    // x + y = 20 with y fixed at 10 fixes x at 10.
    const x = impliedBounds(1.0, 20.0, 20.0, .{ .lower = 10.0, .upper = 10.0 });
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), x.lower, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), x.upper, 1e-12);
    try std.testing.expectEqual(@as(?f64, 10.0), fixValue(x, .continuous));
    try std.testing.expectEqual(@as(?f64, 10.0), fixValue(x, .integer));
}

test "negative coefficient swaps the row sides" {
    // -2x + y = -4 with y in [0,2] gives 2 <= x <= 3.
    const x = impliedBounds(-2.0, -4.0, -4.0, .{ .lower = 0.0, .upper = 2.0 });
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), x.lower, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), x.upper, 1e-12);
    try std.testing.expectEqual(@as(?f64, null), fixValue(x, .continuous));
}

test "negative coefficient can force a fixing" {
    // -2x + y = -4 with y fixed at 2 gives x = 3.
    const x = impliedBounds(-2.0, -4.0, -4.0, .{ .lower = 2.0, .upper = 2.0 });
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), x.lower, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), x.upper, 1e-12);
    const fixed = fixValue(x, .continuous).?;
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), fixed, 1e-12);
}

test "integer rounding fixes the unique integer in an implied interval" {
    const bounds = Bounds{ .lower = 2.2, .upper = 3.1 };
    try std.testing.expectEqual(@as(?f64, 3.0), fixValue(bounds, .integer));
    try std.testing.expectEqual(@as(?f64, null), fixValue(bounds, .continuous));
}

test "binary domains clamp and never fix outside zero one" {
    try std.testing.expectEqual(@as(?f64, 1.0), fixValue(.{ .lower = 0.7, .upper = 1.4 }, .binary));
    try std.testing.expectEqual(@as(?f64, null), fixValue(.{ .lower = 1.2, .upper = 2.0 }, .binary));
}

test "infinite other activity cannot imply a bound" {
    const bounds = impliedBounds(1.0, 0.0, 10.0, .{ .lower = -std.math.inf(f64), .upper = std.math.inf(f64) });
    try std.testing.expectEqual(@as(f64, -std.math.inf(f64)), bounds.lower);
    try std.testing.expectEqual(@as(f64, std.math.inf(f64)), bounds.upper);
    try std.testing.expectEqual(@as(?f64, null), fixValue(bounds, .continuous));
}
