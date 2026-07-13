//! Stable entity handles and dense-index mapping infrastructure.
//!
//! ## Responsibility
//!
//! Provides the DOD-friendly identity layer used by mutable model entities.
//! Handles remain stable while the model compacts its dense arrays; solver
//! kernels should consume dense indices directly and never resolve handles.
//! Slot metadata is stored in separate arrays and released slots are recycled
//! through a free list, so long-running models do not grow monotonically.

const std = @import("std");

pub const EntityKind = enum { generic, variable, constraint, quadratic_constraint, sos, general_constraint };

/// Compile-time-distinct stable identity of one model entity.
///
/// `slot` addresses the metadata arrays; `generation` prevents a released
/// slot from making an old handle silently refer to a newly-created entity.
pub fn TypedEntityId(comptime kind: EntityKind) type {
    return struct {
        slot: u32,
        generation: u64,

        pub const Kind = kind;

        pub fn eql(a: @This(), b: @This()) bool {
            return a.slot == b.slot and a.generation == b.generation;
        }
    };
}

pub const EntityId = TypedEntityId(.generic);
pub const VarId = TypedEntityId(.variable);
pub const ConstrId = TypedEntityId(.constraint);
pub const QConstrId = TypedEntityId(.quadratic_constraint);
pub const SosId = TypedEntityId(.sos);
pub const GenConstrId = TypedEntityId(.general_constraint);

pub const HandleError = error{
    InvalidHandle,
    HandleExhausted,
    DenseIndexOutOfRange,
};

/// Reusable stable-id table with a dense-index reverse map.
///
/// The table deliberately stores each field separately. Resolution touches
/// only `generations` and `dense_indices`; allocation additionally touches the
/// free-list arrays. No per-entity heap allocation is performed.
pub fn HandleTableFor(comptime Id: type) type {
    return struct {
        generations: std.ArrayListUnmanaged(u64) = .empty,
        dense_indices: std.ArrayListUnmanaged(?u32) = .empty,
        next_free: std.ArrayListUnmanaged(?u32) = .empty,
        dense_slots: std.ArrayListUnmanaged(u32) = .empty,
        free_head: ?u32 = null,

        const Self = @This();

        pub fn len(self: Self) usize {
            return self.generations.items.len;
        }

        pub fn liveLen(self: Self) usize {
            return self.dense_slots.items.len;
        }

        pub fn idAtDense(self: Self, dense_index: usize) HandleError!Id {
            if (dense_index >= self.dense_slots.items.len) return error.InvalidHandle;
            const slot = self.dense_slots.items[dense_index];
            return .{ .slot = slot, .generation = self.generations.items[slot] };
        }

        pub fn allocate(self: *Self, allocator: std.mem.Allocator) HandleError!Id {
            if (self.free_head) |slot| {
                const next = self.next_free.items[slot];
                self.free_head = next;
                self.next_free.items[slot] = null;
                self.dense_indices.items[slot] = null;
                return .{ .slot = slot, .generation = self.generations.items[slot] };
            }

            if (self.generations.items.len > std.math.maxInt(u32)) return error.HandleExhausted;
            const slot: u32 = @intCast(self.generations.items.len);
            self.generations.append(allocator, 0) catch return error.HandleExhausted;
            errdefer self.generations.items.len -= 1;
            self.dense_indices.append(allocator, null) catch return error.HandleExhausted;
            errdefer self.dense_indices.items.len -= 1;
            self.next_free.append(allocator, null) catch return error.HandleExhausted;
            return .{ .slot = slot, .generation = 0 };
        }

        /// Bind using the caller's allocator; preferred during model updates.
        pub fn bindDenseWithAllocator(self: *Self, allocator: std.mem.Allocator, id: Id, dense_index: u32) HandleError!void {
            const slot = try self.resolveSlot(id);
            if (self.dense_indices.items[slot] == null) {
                if (dense_index != self.dense_slots.items.len) return error.DenseIndexOutOfRange;
                self.dense_slots.append(allocator, id.slot) catch return error.HandleExhausted;
            }
            self.dense_indices.items[slot] = dense_index;
        }

        pub fn resolve(self: Self, id: Id) HandleError!u32 {
            const slot = try self.resolveSlot(id);
            return self.dense_indices.items[slot] orelse error.InvalidHandle;
        }

        pub fn updateDense(self: *Self, id: Id, dense_index: u32) HandleError!void {
            const slot = try self.resolveSlot(id);
            if (self.dense_indices.items[slot] == null) return error.InvalidHandle;
            self.dense_indices.items[slot] = dense_index;
        }

        pub fn release(self: *Self, id: Id) HandleError!void {
            const slot = try self.resolveSlot(id);
            const dense_index = self.dense_indices.items[slot] orelse return error.InvalidHandle;
            const dense: usize = dense_index;
            const last = self.dense_slots.items.len - 1;
            if (dense >= self.dense_slots.items.len) return error.InvalidHandle;
            if (dense != last) {
                const moved_slot = self.dense_slots.items[last];
                self.dense_slots.items[dense] = moved_slot;
                self.dense_indices.items[moved_slot] = dense_index;
            }
            _ = self.dense_slots.pop();
            self.dense_indices.items[slot] = null;
            if (self.generations.items[slot] == std.math.maxInt(u64)) {
                self.next_free.items[slot] = null;
                return;
            }
            self.generations.items[slot] += 1;
            self.next_free.items[slot] = self.free_head;
            self.free_head = id.slot;
        }

        /// Compact dense entities according to a deletion bitmap.
        ///
        /// Returns an `old_dense -> new_dense` map. Deleted entries contain
        /// `std.math.maxInt(u32)`; surviving entries are updated in-place in the
        /// reverse map. The caller owns the returned map.
        pub fn compact(self: *Self, allocator: std.mem.Allocator, deleted: []const bool) HandleError![]u32 {
            if (deleted.len != self.dense_slots.items.len) return error.DenseIndexOutOfRange;
            const remap = allocator.alloc(u32, deleted.len) catch return error.HandleExhausted;
            errdefer allocator.free(remap);

            var next_dense: usize = 0;
            for (deleted, 0..) |is_deleted, old_dense| {
                const slot = self.dense_slots.items[old_dense];
                if (is_deleted) {
                    remap[old_dense] = std.math.maxInt(u32);
                    self.dense_indices.items[slot] = null;
                    if (self.generations.items[slot] != std.math.maxInt(u64)) {
                        self.generations.items[slot] += 1;
                        self.next_free.items[slot] = self.free_head;
                        self.free_head = slot;
                    }
                } else {
                    remap[old_dense] = @intCast(next_dense);
                    self.dense_slots.items[next_dense] = slot;
                    self.dense_indices.items[slot] = @intCast(next_dense);
                    next_dense += 1;
                }
            }
            self.dense_slots.items.len = next_dense;
            return remap;
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            self.generations.deinit(allocator);
            self.dense_indices.deinit(allocator);
            self.next_free.deinit(allocator);
            self.dense_slots.deinit(allocator);
            self.* = undefined;
        }

        fn resolveSlot(self: Self, id: Id) HandleError!usize {
            if (id.slot >= self.generations.items.len) return error.InvalidHandle;
            const slot: usize = id.slot;
            if (self.generations.items[slot] != id.generation) return error.InvalidHandle;
            return slot;
        }
    };
}

pub const HandleTable = HandleTableFor(EntityId);

test "HandleTable recycles slots without reviving old handles" {
    var table: HandleTable = .{};
    defer table.deinit(std.testing.allocator);

    const first = try table.allocate(std.testing.allocator);
    try table.bindDenseWithAllocator(std.testing.allocator, first, 0);
    try std.testing.expectEqual(@as(u32, 0), try table.resolve(first));

    try table.release(first);
    try std.testing.expectError(error.InvalidHandle, table.resolve(first));

    const second = try table.allocate(std.testing.allocator);
    try std.testing.expectEqual(first.slot, second.slot);
    try std.testing.expect(second.generation != first.generation);
    try table.bindDenseWithAllocator(std.testing.allocator, second, 0);
    try std.testing.expectEqual(@as(u32, 0), try table.resolve(second));
}

test "HandleTable updates the moved dense entity on swap-remove" {
    var table: HandleTable = .{};
    defer table.deinit(std.testing.allocator);

    const first = try table.allocate(std.testing.allocator);
    const second = try table.allocate(std.testing.allocator);
    try table.bindDenseWithAllocator(std.testing.allocator, first, 0);
    try table.bindDenseWithAllocator(std.testing.allocator, second, 1);

    try table.release(first);
    try std.testing.expectError(error.InvalidHandle, table.resolve(first));
    try std.testing.expectEqual(@as(u32, 0), try table.resolve(second));
    try std.testing.expectEqual(@as(usize, 1), table.liveLen());
}

test "HandleTable compact invalidates deleted ids and returns remap" {
    var table: HandleTable = .{};
    defer table.deinit(std.testing.allocator);

    const first = try table.allocate(std.testing.allocator);
    const second = try table.allocate(std.testing.allocator);
    const third = try table.allocate(std.testing.allocator);
    try table.bindDenseWithAllocator(std.testing.allocator, first, 0);
    try table.bindDenseWithAllocator(std.testing.allocator, second, 1);
    try table.bindDenseWithAllocator(std.testing.allocator, third, 2);

    const remap = try table.compact(std.testing.allocator, &[_]bool{ false, true, false });
    defer std.testing.allocator.free(remap);
    try std.testing.expectEqual(@as(u32, 0), remap[0]);
    try std.testing.expectEqual(std.math.maxInt(u32), remap[1]);
    try std.testing.expectEqual(@as(u32, 1), remap[2]);
    try std.testing.expectError(error.InvalidHandle, table.resolve(second));
    try std.testing.expectEqual(@as(u32, 1), try table.resolve(third));
}
