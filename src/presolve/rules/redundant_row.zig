//! Conservative row-activity bounds for redundant-row presolve.

const std = @import("std");

pub const Action = enum {
    keep,
    remove,
    infeasible,
};

/// An outward-rounded interval containing every value of one row activity.
pub const ActivityBounds = struct {
    lower: f64 = 0.0,
    upper: f64 = 0.0,

    pub fn add(self: *ActivityBounds, coefficient: f64, column_lower: f64, column_upper: f64) void {
        if (coefficient == 0.0) return;
        if (!std.math.isFinite(coefficient) or std.math.isNan(column_lower) or std.math.isNan(column_upper)) {
            self.lower = std.math.nan(f64);
            self.upper = std.math.nan(f64);
            return;
        }
        const lower_bound = if (coefficient > 0.0) column_lower else column_upper;
        const upper_bound = if (coefficient > 0.0) column_upper else column_lower;
        self.lower = outwardAdd(self.lower, outwardProduct(coefficient, lower_bound, -std.math.inf(f64)), -std.math.inf(f64));
        self.upper = outwardAdd(self.upper, outwardProduct(coefficient, upper_bound, std.math.inf(f64)), std.math.inf(f64));
    }
};

/// Classify `row_lower <= activity <= row_upper` from a conservative activity
/// interval. `remove` is returned only when every point in the column box
/// satisfies the row. `infeasible` is informative; the current pipeline keeps
/// such rows so simplex can publish a validated Farkas certificate.
pub fn classify(activity: ActivityBounds, row_lower: f64, row_upper: f64) Action {
    if (std.math.isNan(activity.lower) or std.math.isNan(activity.upper) or
        std.math.isNan(row_lower) or std.math.isNan(row_upper)) return .keep;
    if (activity.upper < row_lower or activity.lower > row_upper) return .infeasible;
    if (row_lower <= activity.lower and activity.upper <= row_upper) return .remove;
    return .keep;
}

fn outwardProduct(left: f64, right: f64, direction: f64) f64 {
    const product = left * right;
    if (!std.math.isFinite(product)) return product;
    return std.math.nextAfter(f64, product, direction);
}

fn outwardAdd(left: f64, right: f64, direction: f64) f64 {
    const sum = left + right;
    if (!std.math.isFinite(sum)) return sum;
    return std.math.nextAfter(f64, sum, direction);
}

test "activity bounds classify redundant infeasible and active rows" {
    var activity = ActivityBounds{};
    activity.add(2.0, 0.0, 1.0);
    activity.add(1.0, -2.0, 3.0);
    try std.testing.expect(classify(activity, -3.0, 6.0) == .remove);
    try std.testing.expect(classify(activity, 6.0, 8.0) == .infeasible);
    try std.testing.expect(classify(activity, -1.0, 4.0) == .keep);
}

test "activity bounds remain conservative with negative and unbounded columns" {
    var negative = ActivityBounds{};
    negative.add(-2.0, 1.0, 4.0);
    try std.testing.expect(negative.lower < -8.0);
    try std.testing.expect(negative.upper > -2.0);
    try std.testing.expect(classify(negative, -9.0, -1.0) == .remove);

    var unbounded = ActivityBounds{};
    unbounded.add(1.0, 0.0, std.math.inf(f64));
    try std.testing.expect(classify(unbounded, -std.math.inf(f64), 10.0) == .keep);
}
