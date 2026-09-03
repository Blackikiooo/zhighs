//! Exact substitutions for fixed and structurally empty columns.

const std = @import("std");

pub const Bound = enum {
    fixed,
    lower,
    upper,
    free,
};

pub const Substitution = struct {
    value: f64,
    bound: Bound,
};

pub const Action = union(enum) {
    keep,
    substitute: Substitution,
};

/// A zero-cost, fully free structural column with one matrix entry can serve
/// as the basic variable for that row. Presolve may delete the column and its
/// only row, then reconstruct the column value from the remaining row activity.
pub fn isFreeSingleton(nonzero_count: usize, cost: f64, lower: f64, upper: f64) bool {
    return nonzero_count == 1 and cost == 0.0 and
        lower == -std.math.inf(f64) and upper == std.math.inf(f64);
}

/// Classify a column for exact substitution. `objective_sign` is `1` for
/// minimization and `-1` for maximization, yielding a normalized minimization
/// cost without coupling this rule to the LP engine types.
pub fn classify(
    nonzero_count: usize,
    cost: f64,
    lower: f64,
    upper: f64,
    objective_sign: f64,
) Action {
    if (lower == upper and std.math.isFinite(lower))
        return .{ .substitute = .{ .value = lower, .bound = .fixed } };
    if (nonzero_count != 0) return .keep;

    const normalized_cost = objective_sign * cost;
    if (normalized_cost > 0.0 and std.math.isFinite(lower))
        return .{ .substitute = .{ .value = lower, .bound = .lower } };
    if (normalized_cost < 0.0 and std.math.isFinite(upper))
        return .{ .substitute = .{ .value = upper, .bound = .upper } };
    if (normalized_cost == 0.0) {
        // Objective-neutral empty columns admit any feasible value. Prefer a
        // finite bound for a stable nonbasic basis status; only a truly free
        // column needs an interior canonical zero.
        if (std.math.isFinite(lower))
            return .{ .substitute = .{ .value = lower, .bound = .lower } };
        if (std.math.isFinite(upper))
            return .{ .substitute = .{ .value = upper, .bound = .upper } };
        if (lower == -std.math.inf(f64) and upper == std.math.inf(f64))
            return .{ .substitute = .{ .value = 0.0, .bound = .free } };
    }
    return .keep;
}

test "column substitution selects the finite objective-optimal bound" {
    const lower = classify(0, 2.0, -3.0, 7.0, 1.0);
    try std.testing.expectEqual(@as(f64, -3.0), lower.substitute.value);
    try std.testing.expectEqual(Bound.lower, lower.substitute.bound);

    const upper = classify(0, 2.0, -3.0, 7.0, -1.0);
    try std.testing.expectEqual(@as(f64, 7.0), upper.substitute.value);
    try std.testing.expectEqual(Bound.upper, upper.substitute.bound);
}

test "column substitution canonicalizes every feasible zero-cost empty column" {
    try std.testing.expect(classify(1, 2.0, -3.0, 7.0, 1.0) == .keep);
    const boxed = classify(0, 0.0, -3.0, 7.0, 1.0).substitute;
    try std.testing.expectEqual(@as(f64, -3.0), boxed.value);
    try std.testing.expectEqual(Bound.lower, boxed.bound);
    const upper = classify(0, 0.0, -std.math.inf(f64), 7.0, 1.0).substitute;
    try std.testing.expectEqual(@as(f64, 7.0), upper.value);
    try std.testing.expectEqual(Bound.upper, upper.bound);
    const free = classify(0, -0.0, -std.math.inf(f64), std.math.inf(f64), -1.0).substitute;
    try std.testing.expectEqual(@as(f64, 0.0), free.value);
    try std.testing.expectEqual(Bound.free, free.bound);
}

test "column substitution keeps coupled and objective-improving unbounded columns" {
    try std.testing.expect(classify(1, 0.0, -3.0, 7.0, 1.0) == .keep);
    try std.testing.expect(classify(0, -2.0, 0.0, std.math.inf(f64), 1.0) == .keep);
}

test "fixed substitution takes precedence over structural coupling" {
    const fixed = classify(3, -4.0, 2.0, 2.0, 1.0);
    try std.testing.expectEqual(@as(f64, 2.0), fixed.substitute.value);
    try std.testing.expectEqual(Bound.fixed, fixed.substitute.bound);
}

test "free singleton classification requires exact structural freedom and zero cost" {
    try std.testing.expect(isFreeSingleton(1, -0.0, -std.math.inf(f64), std.math.inf(f64)));
    try std.testing.expect(!isFreeSingleton(0, 0.0, -std.math.inf(f64), std.math.inf(f64)));
    try std.testing.expect(!isFreeSingleton(1, 1.0, -std.math.inf(f64), std.math.inf(f64)));
    try std.testing.expect(!isFreeSingleton(1, 0.0, 0.0, std.math.inf(f64)));
}
