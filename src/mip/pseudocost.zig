//! Solve-local directional pseudocost observations and score preparation.
//!
//! The search owns this mutable state. Branching plugins receive only a
//! read-only score slice, preserving the plugin boundary and avoiding runtime
//! allocation during candidate selection.

const std = @import("std");
const solver = @import("solver");
const plugin = @import("plugin");
const BranchOrigin = @import("node.zig").BranchOrigin;

const DirectionStats = struct {
    unit_cost: f64 = 0.0,
    observations: usize = 0,

    fn add(self: *DirectionStats, value: f64) void {
        self.observations += 1;
        self.unit_cost += (value - self.unit_cost) / @as(f64, @floatFromInt(self.observations));
    }
};

const ColumnStats = struct {
    down: DirectionStats = .{},
    up: DirectionStats = .{},
};

pub const Pseudocost = struct {
    columns: []ColumnStats,
    scores: []plugin.PseudocostScore,
    average_unit_cost: f64 = 1.0,
    observations: usize = 0,
    min_reliable: usize,

    pub fn init(allocator: std.mem.Allocator, num_columns: usize, min_reliable: usize) std.mem.Allocator.Error!Pseudocost {
        const columns = try allocator.alloc(ColumnStats, num_columns);
        errdefer allocator.free(columns);
        @memset(columns, .{});
        const scores = try allocator.alloc(plugin.PseudocostScore, num_columns);
        return .{ .columns = columns, .scores = scores, .min_reliable = min_reliable };
    }

    pub fn deinit(self: *Pseudocost, allocator: std.mem.Allocator) void {
        allocator.free(self.columns);
        allocator.free(self.scores);
        self.* = undefined;
    }

    pub fn observe(self: *Pseudocost, sense: solver.LpObjectiveSense, origin: BranchOrigin, child_objective: f64) void {
        if (origin.column >= self.columns.len or !std.math.isFinite(origin.distance) or origin.distance <= 0.0 or
            !std.math.isFinite(origin.parent_objective) or !std.math.isFinite(child_objective)) return;
        const raw_degradation = switch (sense) {
            .minimize => child_objective - origin.parent_objective,
            .maximize => origin.parent_objective - child_objective,
        };
        const degradation = @max(0.0, raw_degradation);
        const unit_cost = degradation / origin.distance;
        if (!std.math.isFinite(unit_cost)) return;

        switch (origin.direction) {
            .down => self.columns[origin.column].down.add(unit_cost),
            .up => self.columns[origin.column].up.add(unit_cost),
        }
        self.observations += 1;
        self.average_unit_cost += (unit_cost - self.average_unit_cost) / @as(f64, @floatFromInt(self.observations));
    }

    pub fn prepareScores(self: *Pseudocost, primal: []const f64) []const plugin.PseudocostScore {
        std.debug.assert(primal.len == self.columns.len);
        for (self.columns, primal, self.scores) |stats, value, *score| {
            const down_distance = value - @floor(value);
            const up_distance = @ceil(value) - value;
            score.* = .{
                .down_gain = down_distance * self.estimate(stats.down),
                .up_gain = up_distance * self.estimate(stats.up),
                .down_reliable = stats.down.observations >= self.min_reliable,
                .up_reliable = stats.up.observations >= self.min_reliable,
            };
        }
        return self.scores;
    }

    fn estimate(self: *const Pseudocost, stats: DirectionStats) f64 {
        if (stats.observations == 0) return self.average_unit_cost;
        if (self.min_reliable == 0 or stats.observations >= self.min_reliable) return stats.unit_cost;
        const fraction = @as(f64, @floatFromInt(stats.observations)) / @as(f64, @floatFromInt(self.min_reliable));
        const local_weight = 0.9 + 0.1 * fraction;
        return local_weight * stats.unit_cost + (1.0 - local_weight) * self.average_unit_cost;
    }
};

test "directional observations use objective degradation per branch distance" {
    var model = try Pseudocost.init(std.testing.allocator, 2, 2);
    defer model.deinit(std.testing.allocator);
    model.observe(.minimize, .{ .column = 0, .direction = .down, .distance = 0.25, .parent_objective = 10.0 }, 11.0);
    model.observe(.maximize, .{ .column = 0, .direction = .up, .distance = 0.5, .parent_objective = 10.0 }, 8.0);
    const scores = model.prepareScores(&.{ 0.25, 0.5 });
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), scores[0].down_gain, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), scores[0].up_gain, 1e-12);
    try std.testing.expect(!scores[0].down_reliable);
    try std.testing.expectEqual(@as(usize, 2), model.observations);
}

test "zero degradation remains a valid directional observation" {
    var model = try Pseudocost.init(std.testing.allocator, 1, 1);
    defer model.deinit(std.testing.allocator);
    model.observe(.minimize, .{ .column = 0, .direction = .up, .distance = 0.5, .parent_objective = 3.0 }, 2.0);
    const scores = model.prepareScores(&.{0.5});
    try std.testing.expectEqual(@as(f64, 0.0), scores[0].up_gain);
    try std.testing.expect(scores[0].up_reliable);
}
