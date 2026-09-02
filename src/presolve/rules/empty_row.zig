//! Classification for rows with no structural coefficients.

const std = @import("std");

pub const Action = enum {
    keep,
    remove,
    infeasible,
};

/// Classify a row after earlier reductions have updated its bounds and
/// structural nonzero count.
pub fn classify(nonzero_count: usize, lower: f64, upper: f64) Action {
    if (nonzero_count != 0) return .keep;
    if (lower > 0.0 or upper < 0.0) return .infeasible;
    return .remove;
}

test "empty row classification distinguishes redundancy and infeasibility" {
    try std.testing.expectEqual(Action.keep, classify(1, -1.0, 1.0));
    try std.testing.expectEqual(Action.remove, classify(0, -1.0, 1.0));
    try std.testing.expectEqual(Action.remove, classify(0, 0.0, 0.0));
    try std.testing.expectEqual(Action.infeasible, classify(0, 1.0, std.math.inf(f64)));
    try std.testing.expectEqual(Action.infeasible, classify(0, -std.math.inf(f64), -1.0));
}
