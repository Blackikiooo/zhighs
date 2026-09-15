//! Solve-local learned integer no-goods with sparse storage and aging.

const std = @import("std");
const types = @import("types.zig");
const domain = @import("domain.zig");

pub const BoundKind = enum { lower, upper };
pub const Condition = struct { column: usize, kind: BoundKind, value: f64 };
pub const ConflictId = enum(usize) { _ };
const Metadata = struct { start: usize, len: usize, hash: u64, age: u16, active: bool };
const MetadataList = std.MultiArrayList(Metadata);

pub const AddResult = union(enum) { added: ConflictId, existing: ConflictId, rejected };
pub const PropagationResult = struct { feasible: bool, change_head: ?domain.ChangeId, tightenings: usize, checks: usize = 0 };

pub const ConflictPool = struct {
    allocator: std.mem.Allocator,
    entries: std.ArrayListUnmanaged(Condition) = .empty,
    metadata: MetadataList = .empty,
    scratch: std.ArrayListUnmanaged(Condition) = .empty,
    occurrences: []std.ArrayListUnmanaged(ConflictId),
    work: std.ArrayListUnmanaged(ConflictId) = .empty,
    next_work: std.ArrayListUnmanaged(ConflictId) = .empty,
    scheduled_generation: std.ArrayListUnmanaged(u32) = .empty,
    unary_conflicts: std.ArrayListUnmanaged(ConflictId) = .empty,
    generation: u32 = 0,
    last_root_revision: ?usize = null,
    max_entries: usize,
    max_conflict_length: usize,

    pub fn init(allocator: std.mem.Allocator, num_columns: usize, max_entries: usize, max_conflict_length: usize) std.mem.Allocator.Error!ConflictPool {
        const occurrences = try allocator.alloc(std.ArrayListUnmanaged(ConflictId), num_columns);
        @memset(occurrences, .empty);
        return .{ .allocator = allocator, .occurrences = occurrences, .max_entries = max_entries, .max_conflict_length = max_conflict_length };
    }

    pub fn deinit(self: *ConflictPool) void {
        self.entries.deinit(self.allocator);
        self.metadata.deinit(self.allocator);
        self.scratch.deinit(self.allocator);
        for (self.occurrences) |*list| list.deinit(self.allocator);
        self.allocator.free(self.occurrences);
        self.work.deinit(self.allocator);
        self.next_work.deinit(self.allocator);
        self.scheduled_generation.deinit(self.allocator);
        self.unary_conflicts.deinit(self.allocator);
        self.* = undefined;
    }

    pub fn activeCount(self: *const ConflictPool) usize {
        var count: usize = 0;
        for (self.metadata.items(.active)) |active| count += @intFromBool(active);
        return count;
    }

    pub fn totalCount(self: *const ConflictPool) usize {
        return self.metadata.len;
    }

    pub fn shouldCompact(self: *const ConflictPool) bool {
        return self.metadata.len != 0 and self.activeCount() <= self.metadata.len / 2;
    }

    pub fn learnPath(
        self: *ConflictPool,
        arena: *const domain.DomainChangeArena,
        initial_head: ?domain.ChangeId,
        integrality: []const types.Integrality,
        root_lower: []const f64,
        root_upper: []const f64,
        tolerance: f64,
    ) std.mem.Allocator.Error!AddResult {
        std.debug.assert(root_lower.len == integrality.len and root_upper.len == integrality.len);
        self.scratch.clearRetainingCapacity();
        var cursor = initial_head;
        while (cursor) |id| {
            const change = arena.changeAt(id);
            const reason = arena.reasonAt(id);
            cursor = arena.parentOf(id);
            switch (reason) {
                .branch => {},
                else => continue,
            }
            if (change.column >= integrality.len or (integrality[change.column] != .binary and integrality[change.column] != .integer)) return .rejected;
            if (change.lower) |value| {
                if (!std.math.isFinite(value)) return .rejected;
                try self.scratch.append(self.allocator, .{ .column = change.column, .kind = .lower, .value = value });
            }
            if (change.upper) |value| {
                if (!std.math.isFinite(value)) return .rejected;
                try self.scratch.append(self.allocator, .{ .column = change.column, .kind = .upper, .value = value });
            }
        }
        if (self.scratch.items.len == 0) return .rejected;
        std.mem.sort(Condition, self.scratch.items, {}, conditionBefore);
        var canonical_len: usize = 0;
        for (self.scratch.items) |condition| {
            if (canonical_len != 0 and sameKey(self.scratch.items[canonical_len - 1], condition)) {
                const previous = &self.scratch.items[canonical_len - 1];
                previous.value = switch (condition.kind) {
                    .lower => @max(previous.value, condition.value),
                    .upper => @min(previous.value, condition.value),
                };
            } else {
                self.scratch.items[canonical_len] = condition;
                canonical_len += 1;
            }
        }
        const first_canonical = self.scratch.items[0];
        var minimized_len: usize = 0;
        for (self.scratch.items[0..canonical_len]) |condition| {
            if (conditionTrue(condition, root_lower, root_upper, tolerance)) continue;
            self.scratch.items[minimized_len] = condition;
            minimized_len += 1;
        }
        // An empty no-good means the shared root domain itself is infeasible.
        // Keep one globally true condition until search exposes an explicit
        // global-infeasibility result path rather than silently dropping it.
        if (minimized_len == 0) {
            self.scratch.items[0] = first_canonical;
            minimized_len = 1;
        }
        return self.addCanonical(self.scratch.items[0..minimized_len]);
    }

    pub fn addCanonical(self: *ConflictPool, conditions: []const Condition) std.mem.Allocator.Error!AddResult {
        if (conditions.len == 0 or conditions.len > self.max_conflict_length or conditions.len > self.max_entries -| self.entries.items.len) return .rejected;
        for (conditions) |condition| if (condition.column >= self.occurrences.len) return .rejected;
        const hash = hashConditions(conditions);
        for (0..self.metadata.len) |index| {
            if (!self.metadata.items(.active)[index] or self.metadata.items(.hash)[index] != hash or self.metadata.items(.len)[index] != conditions.len) continue;
            const start = self.metadata.items(.start)[index];
            if (conditionsEqual(self.entries.items[start..][0..conditions.len], conditions)) {
                self.metadata.items(.age)[index] = 0;
                return .{ .existing = @enumFromInt(index) };
            }
        }
        try self.entries.ensureUnusedCapacity(self.allocator, conditions.len);
        try self.metadata.ensureUnusedCapacity(self.allocator, 1);
        try self.scheduled_generation.ensureUnusedCapacity(self.allocator, 1);
        if (conditions.len == 1) try self.unary_conflicts.ensureUnusedCapacity(self.allocator, 1);
        for (conditions) |condition| try self.occurrences[condition.column].ensureUnusedCapacity(self.allocator, 1);
        const start = self.entries.items.len;
        const id: ConflictId = @enumFromInt(self.metadata.len);
        self.entries.appendSliceAssumeCapacity(conditions);
        self.metadata.appendAssumeCapacity(.{ .start = start, .len = conditions.len, .hash = hash, .age = 0, .active = true });
        self.scheduled_generation.appendAssumeCapacity(0);
        if (conditions.len == 1) self.unary_conflicts.appendAssumeCapacity(id);
        for (conditions) |condition| self.occurrences[condition.column].appendAssumeCapacity(id);
        return .{ .added = id };
    }

    pub fn propagate(self: *ConflictPool, allocator: std.mem.Allocator, lower: []f64, upper: []f64, arena: *domain.DomainChangeArena, initial_head: ?domain.ChangeId, root_revision: usize, tolerance: f64, max_rounds: usize) std.mem.Allocator.Error!PropagationResult {
        if (max_rounds == 0 or self.metadata.len == 0)
            return .{ .feasible = true, .change_head = initial_head, .tightenings = 0 };
        self.work.clearRetainingCapacity();
        self.next_work.clearRetainingCapacity();
        self.advanceGeneration();
        if (self.last_root_revision == null or self.last_root_revision.? != root_revision) {
            try self.work.ensureUnusedCapacity(self.allocator, self.activeCount());
            for (self.metadata.items(.active), 0..) |active, index| if (active) {
                self.scheduled_generation.items[index] = self.generation;
                self.work.appendAssumeCapacity(@enumFromInt(index));
            };
            self.last_root_revision = root_revision;
        } else {
            var cursor = initial_head;
            while (cursor) |id| {
                try self.scheduleColumnInto(arena.changeAt(id).column, &self.work);
                cursor = arena.parentOf(id);
            }
            for (self.unary_conflicts.items) |id| try self.scheduleConflictInto(id, &self.work);
        }
        var head = initial_head;
        var tightenings: usize = 0;
        var checks: usize = 0;
        for (0..max_rounds) |_| {
            self.advanceGeneration();
            for (self.work.items) |id| {
                checks += 1;
                const index = @intFromEnum(id);
                if (!self.metadata.items(.active)[index]) continue;
                const start = self.metadata.items(.start)[index];
                const conditions = self.entries.items[start..][0..self.metadata.items(.len)[index]];
                var undecided: ?Condition = null;
                var satisfied = false;
                for (conditions) |condition| {
                    if (conditionFalse(condition, lower, upper, tolerance)) {
                        satisfied = true;
                        break;
                    }
                    if (!conditionTrue(condition, lower, upper, tolerance)) {
                        if (undecided != null) {
                            satisfied = true;
                            break;
                        }
                        undecided = condition;
                    }
                }
                if (satisfied) continue;
                self.metadata.items(.age)[index] = 0;
                const condition = undecided orelse return .{ .feasible = false, .change_head = head, .tightenings = tightenings };
                const column = condition.column;
                switch (condition.kind) {
                    .lower => {
                        const value = condition.value - 1.0;
                        if (value < lower[column] - tolerance) return .{ .feasible = false, .change_head = head, .tightenings = tightenings };
                        if (value < upper[column] - tolerance) {
                            upper[column] = value;
                            head = try arena.appendWithReason(allocator, head, .{ .column = column, .upper = value }, .{ .conflict = index });
                            tightenings += 1;
                            try self.scheduleColumn(column);
                        }
                    },
                    .upper => {
                        const value = condition.value + 1.0;
                        if (value > upper[column] + tolerance) return .{ .feasible = false, .change_head = head, .tightenings = tightenings };
                        if (value > lower[column] + tolerance) {
                            lower[column] = value;
                            head = try arena.appendWithReason(allocator, head, .{ .column = column, .lower = value }, .{ .conflict = index });
                            tightenings += 1;
                            try self.scheduleColumn(column);
                        }
                    },
                }
            }
            if (self.next_work.items.len == 0) break;
            const previous = self.work;
            self.work = self.next_work;
            self.next_work = previous;
            self.next_work.clearRetainingCapacity();
        }
        return .{ .feasible = true, .change_head = head, .tightenings = tightenings, .checks = checks };
    }

    pub fn age(self: *ConflictPool, max_age: u16) usize {
        var expired: usize = 0;
        for (self.metadata.items(.active), self.metadata.items(.age)) |*active, *value| {
            if (!active.*) continue;
            if (value.* >= max_age) {
                active.* = false;
                expired += 1;
            } else value.* += 1;
        }
        return expired;
    }

    /// Rebuilds flat storage and column occurrences at an explicit safe point.
    /// Conflict identifiers are solve-internal and no identifier survives this call.
    pub fn compact(self: *ConflictPool) std.mem.Allocator.Error!usize {
        const active_count = self.activeCount();
        if (active_count == self.metadata.len) return 0;
        const removed = self.metadata.len - active_count;
        var replacement = try ConflictPool.init(self.allocator, self.occurrences.len, self.max_entries, self.max_conflict_length);
        errdefer replacement.deinit();
        for (0..self.metadata.len) |index| {
            if (!self.metadata.items(.active)[index]) continue;
            const start = self.metadata.items(.start)[index];
            const conditions = self.entries.items[start..][0..self.metadata.items(.len)[index]];
            const result = try replacement.addCanonical(conditions);
            const new_id = switch (result) {
                .added => |id| id,
                else => unreachable,
            };
            replacement.metadata.items(.age)[@intFromEnum(new_id)] = self.metadata.items(.age)[index];
        }
        self.deinit();
        self.* = replacement;
        return removed;
    }

    fn advanceGeneration(self: *ConflictPool) void {
        self.generation +%= 1;
        if (self.generation == 0) {
            @memset(self.scheduled_generation.items, 0);
            self.generation = 1;
        }
    }

    fn scheduleColumn(self: *ConflictPool, column: usize) std.mem.Allocator.Error!void {
        try self.scheduleColumnInto(column, &self.next_work);
    }

    fn scheduleColumnInto(self: *ConflictPool, column: usize, target: *std.ArrayListUnmanaged(ConflictId)) std.mem.Allocator.Error!void {
        for (self.occurrences[column].items) |id| try self.scheduleConflictInto(id, target);
    }

    fn scheduleConflictInto(self: *ConflictPool, id: ConflictId, target: *std.ArrayListUnmanaged(ConflictId)) std.mem.Allocator.Error!void {
        const index = @intFromEnum(id);
        if (!self.metadata.items(.active)[index] or self.scheduled_generation.items[index] == self.generation) return;
        self.scheduled_generation.items[index] = self.generation;
        try target.append(self.allocator, id);
    }
};

fn conditionTrue(c: Condition, lower: []const f64, upper: []const f64, tolerance: f64) bool {
    return switch (c.kind) {
        .lower => lower[c.column] >= c.value - tolerance,
        .upper => upper[c.column] <= c.value + tolerance,
    };
}
fn conditionFalse(c: Condition, lower: []const f64, upper: []const f64, tolerance: f64) bool {
    return switch (c.kind) {
        .lower => upper[c.column] <= c.value - 1.0 + tolerance,
        .upper => lower[c.column] >= c.value + 1.0 - tolerance,
    };
}
fn sameKey(a: Condition, b: Condition) bool {
    return a.column == b.column and a.kind == b.kind;
}
fn conditionsEqual(a: []const Condition, b: []const Condition) bool {
    if (a.len != b.len) return false;
    for (a, b) |left, right| if (left.column != right.column or left.kind != right.kind or left.value != right.value) return false;
    return true;
}
fn conditionBefore(_: void, a: Condition, b: Condition) bool {
    return a.column < b.column or (a.column == b.column and @intFromEnum(a.kind) < @intFromEnum(b.kind));
}
fn hashConditions(conditions: []const Condition) u64 {
    var h = std.hash.Wyhash.init(0x636f_6e66_6c69_6374);
    for (conditions) |condition| {
        h.update(std.mem.asBytes(&condition.column));
        const kind: u8 = @intFromEnum(condition.kind);
        h.update(std.mem.asBytes(&kind));
        h.update(std.mem.asBytes(&condition.value));
    }
    return h.final();
}

test "learned integer no-good propagates its last undecided condition" {
    var pool = try ConflictPool.init(std.testing.allocator, 2, 100, 10);
    defer pool.deinit();
    _ = try pool.addCanonical(&.{ .{ .column = 0, .kind = .lower, .value = 1.0 }, .{ .column = 1, .kind = .upper, .value = 0.0 } });
    var lower = [_]f64{ 1.0, 0.0 };
    var upper = [_]f64{ 1.0, 1.0 };
    var arena = domain.DomainChangeArena{};
    defer arena.deinit(std.testing.allocator);
    const result = try pool.propagate(std.testing.allocator, &lower, &upper, &arena, null, 0, 1e-6, 4);
    try std.testing.expect(result.feasible);
    try std.testing.expectEqual(@as(f64, 1.0), lower[1]);
    try std.testing.expectEqual(@as(usize, 1), result.tightenings);
}

test "column occurrences schedule transitive no-good propagation" {
    var pool = try ConflictPool.init(std.testing.allocator, 3, 100, 10);
    defer pool.deinit();
    _ = try pool.addCanonical(&.{ .{ .column = 0, .kind = .lower, .value = 1.0 }, .{ .column = 1, .kind = .upper, .value = 0.0 } });
    _ = try pool.addCanonical(&.{ .{ .column = 1, .kind = .lower, .value = 1.0 }, .{ .column = 2, .kind = .lower, .value = 1.0 } });
    var lower = [_]f64{ 1.0, 0.0, 0.0 };
    var upper = [_]f64{ 1.0, 1.0, 1.0 };
    var arena = domain.DomainChangeArena{};
    defer arena.deinit(std.testing.allocator);
    const result = try pool.propagate(std.testing.allocator, &lower, &upper, &arena, null, 0, 1e-6, 4);
    try std.testing.expect(result.feasible);
    try std.testing.expectEqual(@as(usize, 2), result.tightenings);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0, 0.0 }, &lower);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0, 0.0 }, &upper);
}

test "stable root revisions schedule only sparse path occurrences" {
    var pool = try ConflictPool.init(std.testing.allocator, 6, 100, 10);
    defer pool.deinit();
    _ = try pool.addCanonical(&.{ .{ .column = 0, .kind = .lower, .value = 1.0 }, .{ .column = 1, .kind = .upper, .value = 0.0 } });
    _ = try pool.addCanonical(&.{ .{ .column = 2, .kind = .lower, .value = 1.0 }, .{ .column = 3, .kind = .upper, .value = 0.0 } });
    _ = try pool.addCanonical(&.{ .{ .column = 4, .kind = .lower, .value = 1.0 }, .{ .column = 5, .kind = .upper, .value = 0.0 } });
    var arena = domain.DomainChangeArena{};
    defer arena.deinit(std.testing.allocator);
    var lower = [_]f64{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    var upper = [_]f64{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
    _ = try pool.propagate(std.testing.allocator, &lower, &upper, &arena, null, 7, 1e-6, 4);
    const head = try arena.append(std.testing.allocator, null, .{ .column = 0, .lower = 1.0 });
    lower[0] = 1.0;
    const result = try pool.propagate(std.testing.allocator, &lower, &upper, &arena, head, 7, 1e-6, 4);
    try std.testing.expectEqual(@as(usize, 1), result.tightenings);
    try std.testing.expectEqual(@as(f64, 1.0), lower[1]);
    try std.testing.expect(result.checks < pool.activeCount());
}

test "root revision changes refresh every active conflict" {
    var pool = try ConflictPool.init(std.testing.allocator, 2, 100, 10);
    defer pool.deinit();
    _ = try pool.addCanonical(&.{ .{ .column = 0, .kind = .lower, .value = 1.0 }, .{ .column = 1, .kind = .upper, .value = 0.0 } });
    var arena = domain.DomainChangeArena{};
    defer arena.deinit(std.testing.allocator);
    var lower = [_]f64{ 0.0, 0.0 };
    var upper = [_]f64{ 1.0, 1.0 };
    _ = try pool.propagate(std.testing.allocator, &lower, &upper, &arena, null, 0, 1e-6, 4);
    lower[0] = 1.0;
    const result = try pool.propagate(std.testing.allocator, &lower, &upper, &arena, null, 1, 1e-6, 4);
    try std.testing.expectEqual(@as(usize, 1), result.tightenings);
    try std.testing.expectEqual(@as(f64, 1.0), lower[1]);
}

test "path learning canonicalizes repeated bound changes" {
    var arena = domain.DomainChangeArena{};
    defer arena.deinit(std.testing.allocator);
    const first = try arena.appendWithReason(std.testing.allocator, null, .{ .column = 0, .lower = 1.0 }, .branch);
    const second = try arena.appendWithReason(std.testing.allocator, first, .{ .column = 0, .lower = 2.0 }, .branch);
    const third = try arena.appendWithReason(std.testing.allocator, second, .{ .column = 1, .upper = 3.0 }, .branch);
    var pool = try ConflictPool.init(std.testing.allocator, 2, 100, 10);
    defer pool.deinit();
    try std.testing.expectEqual(std.meta.Tag(AddResult).added, std.meta.activeTag(try pool.learnPath(&arena, third, &.{ .integer, .integer }, &.{ 0.0, 0.0 }, &.{ 10.0, 10.0 }, 1e-6)));
    try std.testing.expectEqual(@as(usize, 2), pool.entries.items.len);
    try std.testing.expectEqual(@as(f64, 2.0), pool.entries.items[0].value);
}

test "path learning removes conditions implied by the shared root domain" {
    var arena = domain.DomainChangeArena{};
    defer arena.deinit(std.testing.allocator);
    const root_implied = try arena.appendWithReason(std.testing.allocator, null, .{ .column = 0, .lower = 1.0 }, .branch);
    const local = try arena.appendWithReason(std.testing.allocator, root_implied, .{ .column = 1, .upper = 2.0 }, .branch);
    var pool = try ConflictPool.init(std.testing.allocator, 2, 100, 10);
    defer pool.deinit();
    _ = try pool.learnPath(&arena, local, &.{ .integer, .integer }, &.{ 1.0, 0.0 }, &.{ 10.0, 10.0 }, 1e-6);
    try std.testing.expectEqual(@as(usize, 1), pool.entries.items.len);
    try std.testing.expectEqual(@as(usize, 1), pool.entries.items[0].column);
}

test "compaction reclaims expired ranges and rebuilds occurrences" {
    var pool = try ConflictPool.init(std.testing.allocator, 3, 100, 10);
    defer pool.deinit();
    _ = try pool.addCanonical(&.{ .{ .column = 0, .kind = .lower, .value = 1.0 }, .{ .column = 1, .kind = .upper, .value = 0.0 } });
    _ = try pool.addCanonical(&.{ .{ .column = 0, .kind = .upper, .value = 0.0 }, .{ .column = 2, .kind = .upper, .value = 0.0 } });
    try std.testing.expectEqual(@as(usize, 0), pool.age(1));
    var lower = [_]f64{ 1.0, 0.0, 0.0 };
    var upper = [_]f64{ 1.0, 1.0, 1.0 };
    var arena = domain.DomainChangeArena{};
    defer arena.deinit(std.testing.allocator);
    _ = try pool.propagate(std.testing.allocator, &lower, &upper, &arena, null, 0, 1e-6, 4);
    try std.testing.expectEqual(@as(usize, 1), pool.age(1));
    try std.testing.expect(pool.shouldCompact());
    try std.testing.expectEqual(@as(usize, 1), try pool.compact());
    try std.testing.expectEqual(@as(usize, 1), pool.totalCount());
    try std.testing.expectEqual(@as(usize, 1), pool.activeCount());
    lower = .{ 1.0, 0.0, 0.0 };
    upper = .{ 1.0, 1.0, 1.0 };
    const result = try pool.propagate(std.testing.allocator, &lower, &upper, &arena, null, 0, 1e-6, 4);
    try std.testing.expect(result.tightenings > 0);
}
