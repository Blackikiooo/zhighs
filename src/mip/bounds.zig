//! Objective-bound comparison and MIP gap policy.

const std = @import("std");
const solver = @import("solver");
const Node = @import("node.zig").Node;
const SolveControl = @import("types.zig").SolveControl;

pub fn betterObjective(sense: solver.LpObjectiveSense, candidate: f64, incumbent: f64) bool {
    return switch (sense) {
        .minimize => candidate < incumbent,
        .maximize => candidate > incumbent,
    };
}

pub fn unknown(sense: solver.LpObjectiveSense) f64 {
    return switch (sense) {
        .minimize => -std.math.inf(f64),
        .maximize => std.math.inf(f64),
    };
}

pub fn cannotImprove(sense: solver.LpObjectiveSense, bound: f64, incumbent: f64) bool {
    return switch (sense) {
        .minimize => bound >= incumbent,
        .maximize => bound <= incumbent,
    };
}

pub fn gapSatisfied(bound: f64, incumbent: f64, control: SolveControl) bool {
    const gap = @abs(incumbent - bound);
    return gap <= control.absolute_gap or gap / @max(1.0, @abs(incumbent)) <= control.relative_gap;
}

pub fn bestOpen(sense: solver.LpObjectiveSense, nodes: []const Node, incumbent: f64, have_incumbent: bool) f64 {
    var result: ?f64 = null;
    for (nodes) |node| {
        if (node.bound) |bound| {
            result = if (result) |old| switch (sense) {
                .minimize => @min(old, bound),
                .maximize => @max(old, bound),
            } else bound;
        }
    }
    return result orelse if (have_incumbent) incumbent else 0.0;
}

pub fn combine(sense: solver.LpObjectiveSense, left: f64, right: f64) f64 {
    return switch (sense) {
        .minimize => @min(left, right),
        .maximize => @max(left, right),
    };
}

test "gap policy supports absolute and relative termination" {
    try std.testing.expect(gapSatisfied(99.95, 100.0, .{ .absolute_gap = 0.1, .relative_gap = 0.0 }));
    try std.testing.expect(gapSatisfied(99.0, 100.0, .{ .absolute_gap = 0.0, .relative_gap = 0.02 }));
    try std.testing.expect(!gapSatisfied(90.0, 100.0, .{ .absolute_gap = 0.1, .relative_gap = 0.02 }));
}
