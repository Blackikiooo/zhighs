//! Persistent root-LP reduced-cost certificate and global lurking bounds.

const std = @import("std");
const solver = @import("solver");
const types = @import("types.zig");
const domain = @import("domain.zig");
const fixing = @import("reduced_cost_fixing.zig");

pub const ApplyResult = struct { feasible: bool = true, tightenings: usize = 0 };

pub const RootReducedCosts = struct {
    allocator: std.mem.Allocator,
    lower: []f64,
    upper: []f64,
    costs: []f64,
    objective: f64 = 0.0,
    valid: bool = false,

    pub fn init(allocator: std.mem.Allocator, num_columns: usize) std.mem.Allocator.Error!RootReducedCosts {
        const lower = try allocator.alloc(f64, num_columns);
        errdefer allocator.free(lower);
        const upper = try allocator.alloc(f64, num_columns);
        errdefer allocator.free(upper);
        const costs = try allocator.alloc(f64, num_columns);
        return .{ .allocator = allocator, .lower = lower, .upper = upper, .costs = costs };
    }

    pub fn deinit(self: *RootReducedCosts) void {
        self.allocator.free(self.lower);
        self.allocator.free(self.upper);
        self.allocator.free(self.costs);
        self.* = undefined;
    }

    pub fn capture(self: *RootReducedCosts, lower: []const f64, upper: []const f64, costs: []const f64, objective: f64) void {
        std.debug.assert(lower.len == self.lower.len and upper.len == self.upper.len and costs.len == self.costs.len);
        @memcpy(self.lower, lower);
        @memcpy(self.upper, upper);
        @memcpy(self.costs, costs);
        self.objective = objective;
        self.valid = true;
    }

    pub fn apply(self: *const RootReducedCosts, sense: solver.LpObjectiveSense, integrality: []const types.Integrality, target: *domain.Domain, incumbent: f64, tolerance: f64) ApplyResult {
        if (!self.valid) return .{};
        const gap = fixing.cutoffGap(sense, self.objective, incumbent);
        if (!std.math.isFinite(gap) or gap < -tolerance) return .{};
        var count: usize = 0;
        for (integrality, self.lower, self.upper, self.costs, 0..) |kind, lower, upper, cost, column| {
            const implied = fixing.impliedBounds(kind, lower, upper, cost, @max(0.0, gap), tolerance);
            if (implied.lower) |value| if (value > target.root_lower[column] + tolerance) {
                target.root_lower[column] = value;
                target.lower[column] = @max(target.lower[column], value);
                count += 1;
            };
            if (implied.upper) |value| if (value < target.root_upper[column] - tolerance) {
                target.root_upper[column] = value;
                target.upper[column] = @min(target.upper[column], value);
                count += 1;
            };
        }
        for (target.lower, target.upper) |lower, upper| if (lower > upper + tolerance) return .{ .feasible = false, .tightenings = count };
        if (count != 0) target.root_revision +%= 1;
        return .{ .tightenings = count };
    }
};

test "root certificate applies stronger bounds as the cutoff improves" {
    var root = try RootReducedCosts.init(std.testing.allocator, 1);
    defer root.deinit();
    root.capture(&.{0.0}, &.{10.0}, &.{2.0}, 0.0);
    var target = (try domain.Domain.init(std.testing.allocator, &.{.integer}, &.{0.0}, &.{10.0}, 1e-6)).?;
    defer target.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 1), root.apply(.minimize, &.{.integer}, &target, 7.0, 1e-6).tightenings);
    try std.testing.expectEqual(@as(f64, 3.0), target.root_upper[0]);
    try std.testing.expectEqual(@as(usize, 1), root.apply(.minimize, &.{.integer}, &target, 3.0, 1e-6).tightenings);
    try std.testing.expectEqual(@as(f64, 1.0), target.root_upper[0]);
}
