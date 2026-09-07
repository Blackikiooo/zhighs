//! Built-in node selectors.
//!
//! Algorithm implementations live here; registration, ordering, timing, and
//! lifecycle remain the responsibility of `framework`.

const std = @import("std");
const solver = @import("solver");
const plugin = @import("plugin");

/// Depth-first search: always processes the most recently created node, which
/// drives the open set towards a leaf and finds incumbents early.
pub const Dfs = struct {
    pub const name = "dfs";

    pub fn select(_: *Dfs, ctx: plugin.NodeSelectionContext) ?usize {
        if (ctx.open.len == 0) return null;
        return ctx.open.len - 1;
    }
};

/// Best-bound search: always processes the open node with the strongest dual
/// bound, so the global dual bound improves as early as possible. Nodes whose
/// bound is still unknown (the root, before any relaxation has been solved)
/// are treated optimistically as the best possible bound.
pub const BestBound = struct {
    pub const name = "best-bound";

    pub fn select(_: *BestBound, ctx: plugin.NodeSelectionContext) ?usize {
        if (ctx.open.len == 0) return null;
        var best_index: usize = 0;
        var best_bound = optimisticBound(ctx.objective_sense, ctx.open[0].bound);
        for (ctx.open[1..], 1..) |node, index| {
            const bound = optimisticBound(ctx.objective_sense, node.bound);
            const better = switch (ctx.objective_sense) {
                .minimize => bound < best_bound,
                .maximize => bound > best_bound,
            };
            if (better) {
                best_index = index;
                best_bound = bound;
            }
        }
        return best_index;
    }

    fn optimisticBound(sense: solver.LpObjectiveSense, bound: ?f64) f64 {
        return bound orelse switch (sense) {
            .minimize => -std.math.inf(f64),
            .maximize => std.math.inf(f64),
        };
    }
};

test "dfs selects the newest node" {
    var selector = Dfs{};
    const open = [_]plugin.OpenNode{
        .{ .bound = 1.0, .depth = 0 },
        .{ .bound = 2.0, .depth = 1 },
        .{ .bound = 3.0, .depth = 2 },
    };
    try std.testing.expectEqual(@as(?usize, 2), selector.select(.{ .objective_sense = .minimize, .open = &open }));
}

test "dfs returns null on an empty open set" {
    var selector = Dfs{};
    try std.testing.expectEqual(@as(?usize, null), selector.select(.{ .objective_sense = .minimize, .open = &.{} }));
}

test "best-bound prefers the strongest bound when minimizing" {
    var selector = BestBound{};
    const open = [_]plugin.OpenNode{
        .{ .bound = 5.0, .depth = 0 },
        .{ .bound = 3.0, .depth = 1 },
        .{ .bound = 4.0, .depth = 1 },
    };
    try std.testing.expectEqual(@as(?usize, 1), selector.select(.{ .objective_sense = .minimize, .open = &open }));
}

test "best-bound prefers the strongest bound when maximizing" {
    var selector = BestBound{};
    const open = [_]plugin.OpenNode{
        .{ .bound = 5.0, .depth = 0 },
        .{ .bound = 3.0, .depth = 1 },
        .{ .bound = 4.0, .depth = 1 },
    };
    try std.testing.expectEqual(@as(?usize, 0), selector.select(.{ .objective_sense = .maximize, .open = &open }));
}

test "best-bound treats an unknown bound as the best possible" {
    var selector = BestBound{};
    const open = [_]plugin.OpenNode{
        .{ .bound = 5.0, .depth = 0 },
        .{ .bound = null, .depth = 1 },
    };
    try std.testing.expectEqual(@as(?usize, 1), selector.select(.{ .objective_sense = .minimize, .open = &open }));
}
