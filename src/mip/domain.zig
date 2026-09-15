//! Root-domain normalization and discrete-domain predicates.

const std = @import("std");
const types = @import("types.zig");

pub const BoundChange = struct {
    column: usize,
    lower: ?f64 = null,
    upper: ?f64 = null,
};

pub const ChangeReason = union(enum) {
    branch,
    row: usize,
    clique,
    conflict: usize,
    reduced_cost,
    unspecified,
};

/// Strongly typed index into a solve-local persistent change arena.
pub const ChangeId = enum(usize) {
    _,

    fn fromIndex(index: usize) ChangeId {
        return @enumFromInt(index);
    }

    fn toIndex(self: ChangeId) usize {
        return @intFromEnum(self);
    }
};

pub const DomainChangeArena = struct {
    const Entry = struct {
        change: BoundChange,
        parent: ?ChangeId,
        reason: ChangeReason,
    };

    entries: std.ArrayListUnmanaged(Entry) = .empty,

    pub fn deinit(self: *DomainChangeArena, allocator: std.mem.Allocator) void {
        self.entries.deinit(allocator);
        self.* = undefined;
    }

    pub fn append(
        self: *DomainChangeArena,
        allocator: std.mem.Allocator,
        parent: ?ChangeId,
        change: BoundChange,
    ) std.mem.Allocator.Error!ChangeId {
        return self.appendWithReason(allocator, parent, change, .unspecified);
    }

    pub fn appendWithReason(
        self: *DomainChangeArena,
        allocator: std.mem.Allocator,
        parent: ?ChangeId,
        change: BoundChange,
        reason: ChangeReason,
    ) std.mem.Allocator.Error!ChangeId {
        const id = ChangeId.fromIndex(self.entries.items.len);
        try self.entries.append(allocator, .{ .change = change, .parent = parent, .reason = reason });
        return id;
    }

    pub fn changeAt(self: *const DomainChangeArena, id: ChangeId) BoundChange {
        return self.entries.items[id.toIndex()].change;
    }

    pub fn parentOf(self: *const DomainChangeArena, id: ChangeId) ?ChangeId {
        return self.entries.items[id.toIndex()].parent;
    }

    pub fn reasonAt(self: *const DomainChangeArena, id: ChangeId) ChangeReason {
        return self.entries.items[id.toIndex()].reason;
    }
};

/// Reusable full-domain workspace reconstructed from sparse node changes.
/// Root and active bounds are immutable after initialization; `load` mutates
/// only the cache passed to the LP relaxation and plugin contexts.
pub const Domain = struct {
    root_lower: []f64,
    root_upper: []f64,
    active_lower: []f64,
    active_upper: []f64,
    lower: []f64,
    upper: []f64,
    root_revision: usize = 0,

    pub fn init(
        allocator: std.mem.Allocator,
        integrality: []const types.Integrality,
        original_lower: []const f64,
        original_upper: []const f64,
        tolerance: f64,
    ) std.mem.Allocator.Error!?Domain {
        const root_lower = try allocator.dupe(f64, original_lower);
        errdefer allocator.free(root_lower);
        const root_upper = try allocator.dupe(f64, original_upper);
        errdefer allocator.free(root_upper);
        const active_lower = try allocator.alloc(f64, original_lower.len);
        errdefer allocator.free(active_lower);
        const active_upper = try allocator.alloc(f64, original_upper.len);
        errdefer allocator.free(active_upper);
        const lower = try allocator.alloc(f64, original_lower.len);
        errdefer allocator.free(lower);
        const upper = try allocator.alloc(f64, original_upper.len);
        errdefer allocator.free(upper);

        if (!normalizeRootBounds(integrality, root_lower, root_upper, tolerance)) {
            allocator.free(root_lower);
            allocator.free(root_upper);
            allocator.free(active_lower);
            allocator.free(active_upper);
            allocator.free(lower);
            allocator.free(upper);
            return null;
        }
        computeActiveBounds(integrality, original_lower, original_upper, active_lower, active_upper, tolerance);
        @memcpy(lower, root_lower);
        @memcpy(upper, root_upper);
        return .{
            .root_lower = root_lower,
            .root_upper = root_upper,
            .active_lower = active_lower,
            .active_upper = active_upper,
            .lower = lower,
            .upper = upper,
        };
    }

    pub fn deinit(self: *Domain, allocator: std.mem.Allocator) void {
        allocator.free(self.root_lower);
        allocator.free(self.root_upper);
        allocator.free(self.active_lower);
        allocator.free(self.active_upper);
        allocator.free(self.lower);
        allocator.free(self.upper);
        self.* = undefined;
    }

    pub fn load(self: *Domain, arena: *const DomainChangeArena, change_head: ?ChangeId) bool {
        @memcpy(self.lower, self.root_lower);
        @memcpy(self.upper, self.root_upper);
        var cursor = change_head;
        while (cursor) |id| {
            const entry = arena.entries.items[id.toIndex()];
            const change = entry.change;
            if (change.column >= self.lower.len) return false;
            if (change.lower) |value| self.lower[change.column] = @max(self.lower[change.column], value);
            if (change.upper) |value| self.upper[change.column] = @min(self.upper[change.column], value);
            if (self.lower[change.column] > self.upper[change.column]) return false;
            cursor = entry.parent;
        }
        return true;
    }
};

pub fn computeActiveBounds(
    integrality: []const types.Integrality,
    lower: []const f64,
    upper: []const f64,
    active_lower: []f64,
    active_upper: []f64,
    tolerance: f64,
) void {
    for (integrality, lower, upper, active_lower, active_upper) |kind, lo, up, *active_lo, *active_up| {
        active_lo.* = if (kind == .semi_integer and std.math.isFinite(lo)) @ceil(lo - tolerance) else lo;
        active_up.* = if (kind == .semi_integer and std.math.isFinite(up)) @floor(up + tolerance) else up;
    }
}

pub fn normalizeRootBounds(integrality: []const types.Integrality, lower: []f64, upper: []f64, tolerance: f64) bool {
    for (integrality, lower, upper) |kind, *lo, *up| {
        if (kind == .continuous) continue;
        switch (kind) {
            .continuous => unreachable,
            .binary => {
                lo.* = @max(lo.*, 0.0);
                up.* = @min(up.*, 1.0);
            },
            .integer => {},
            .semi_continuous, .semi_integer => {
                const active_lo = if (kind == .semi_integer and std.math.isFinite(lo.*)) @ceil(lo.* - tolerance) else lo.*;
                const active_up = if (kind == .semi_integer and std.math.isFinite(up.*)) @floor(up.* + tolerance) else up.*;
                lo.* = @min(0.0, active_lo);
                up.* = @max(0.0, active_up);
                continue;
            },
        }
        if (std.math.isFinite(lo.*)) lo.* = @ceil(lo.* - tolerance);
        if (std.math.isFinite(up.*)) up.* = @floor(up.* + tolerance);
        if (lo.* > up.*) return false;
    }
    return true;
}

pub fn hasDiscreteColumns(integrality: []const types.Integrality) bool {
    for (integrality) |kind| if (kind != .continuous) return true;
    return false;
}

test "root domains normalize binary and integer bounds" {
    var lower = [_]f64{ -5.0, 0.2, 2.2 };
    var upper = [_]f64{ 10.0, 4.8, 4.8 };
    try std.testing.expect(normalizeRootBounds(&.{ .binary, .integer, .semi_integer }, &lower, &upper, 1e-6));
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 1.0, 0.0 }, &lower);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 4.0, 4.0 }, &upper);
}

test "empty integer lattice is rejected" {
    var lower = [_]f64{1.2};
    var upper = [_]f64{1.8};
    try std.testing.expect(!normalizeRootBounds(&.{.integer}, &lower, &upper, 1e-6));
}

test "domain replays sparse changes without mutating root bounds" {
    var value = (try Domain.init(
        std.testing.allocator,
        &.{ .integer, .continuous, .binary },
        &.{ 0.0, -10.0, -2.0 },
        &.{ 10.0, 10.0, 3.0 },
        1e-6,
    )).?;
    defer value.deinit(std.testing.allocator);
    var arena = DomainChangeArena{};
    defer arena.deinit(std.testing.allocator);
    const first = try arena.append(std.testing.allocator, null, .{ .column = 0, .lower = 4.0 });
    const second = try arena.append(std.testing.allocator, first, .{ .column = 1, .upper = 2.0 });
    try std.testing.expect(value.load(&arena, second));
    try std.testing.expectEqualSlices(f64, &.{ 4.0, -10.0, 0.0 }, value.lower);
    try std.testing.expectEqualSlices(f64, &.{ 10.0, 2.0, 1.0 }, value.upper);
    try std.testing.expect(value.load(&arena, null));
    try std.testing.expectEqualSlices(f64, &.{ 0.0, -10.0, 0.0 }, value.lower);
    try std.testing.expectEqualSlices(f64, &.{ 10.0, 10.0, 1.0 }, value.upper);
}
