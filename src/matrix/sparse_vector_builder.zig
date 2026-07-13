//! Mutable construction and canonicalization of sparse vectors.

const std = @import("std");
const Foundation = @import("foundation");
const sparse = @import("sparse_vector.zig");

pub fn SparseVectorBuilder(comptime Id: type) type {
    // Instantiating the output type performs the RowId/ColId validation.
    const Vector = sparse.SparseVector(Id);
    const Entry = struct {
        id: Id,
        value: f64,
    };
    const EntryList = std.MultiArrayList(Entry);

    return struct {
        dimension: usize,
        entries: EntryList = .empty,

        const Self = @This();

        const SortContext = struct {
            ids: []const Id,

            fn lessThan(self: @This(), lhs: usize, rhs: usize) bool {
                return self.ids[lhs].raw() < self.ids[rhs].raw();
            }
        };

        pub fn init(dimension: usize) Self {
            return .{ .dimension = dimension };
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            self.entries.deinit(allocator);
            self.* = undefined;
        }

        pub fn clearRetainingCapacity(self: *Self) void {
            self.entries.clearRetainingCapacity();
        }

        /// Reserves the single SoA entry allocation for a hot append loop.
        pub fn reserve(self: *Self, allocator: std.mem.Allocator, additional: usize) std.mem.Allocator.Error!void {
            try self.entries.ensureUnusedCapacity(allocator, additional);
        }

        /// Appends an unordered entry. Explicit zeros and duplicates are
        /// accepted and removed or merged by `freeze`.
        pub fn append(
            self: *Self,
            allocator: std.mem.Allocator,
            id: Id,
            value: f64,
        ) (std.mem.Allocator.Error || sparse.SparseVectorError)!void {
            if (id.toUsize() >= self.dimension) {
                return sparse.SparseVectorError.IndexOutOfBounds;
            }
            if (!std.math.isFinite(value)) {
                return sparse.SparseVectorError.NonFiniteValue;
            }
            if (self.entries.len >= self.entries.capacity)
                try self.entries.ensureUnusedCapacity(allocator, 1);
            self.appendPreReserved(id, value);
        }

        /// Trusted append for callers that hoist validation and reservation.
        /// Caller guarantees capacity, an in-range ID, and a finite value.
        pub fn appendPreReserved(self: *Self, id: Id, value: f64) void {
            std.debug.assert(self.entries.len < self.entries.capacity);
            std.debug.assert(id.toUsize() < self.dimension);
            std.debug.assert(std.math.isFinite(value));
            self.entries.appendAssumeCapacity(.{ .id = id, .value = value });
        }

        /// Sorts, merges duplicates, removes values within `drop_tolerance`,
        /// and returns an allocator-owned canonical vector.
        pub fn freeze(
            self: *Self,
            allocator: std.mem.Allocator,
            drop_tolerance: f64,
        ) (std.mem.Allocator.Error || sparse.SparseVectorError)!Vector {
            if (!std.math.isFinite(drop_tolerance) or drop_tolerance < 0.0) {
                return sparse.SparseVectorError.InvalidTolerance;
            }
            if (self.entries.len == 0) {
                return Vector.empty(self.dimension);
            }

            // Stable ordering preserves append order among duplicate IDs, so
            // floating-point merging is deterministic with respect to input.
            var fields = self.entries.slice();
            self.entries.sort(SortContext{ .ids = fields.items(.id) });
            fields = self.entries.slice();
            const ids = fields.items(.id);
            const entry_values = fields.items(.value);

            var read: usize = 0;
            var write: usize = 0;
            while (read < self.entries.len) {
                const id = ids[read];
                var sum: f64 = 0.0;

                while (read < self.entries.len and
                    ids[read].raw() == id.raw()) : (read += 1)
                {
                    sum += entry_values[read];
                }

                if (!std.math.isFinite(sum)) {
                    return sparse.SparseVectorError.NonFiniteValue;
                }
                if (@abs(sum) <= drop_tolerance) continue;

                fields.set(write, .{ .id = id, .value = sum });
                write += 1;
            }
            self.entries.shrinkRetainingCapacity(write);

            if (write == 0) return Vector.empty(self.dimension);

            var result = try Vector.initPackedUninitialized(allocator, self.dimension, write);
            errdefer result.deinit(allocator);
            fields = self.entries.slice();
            @memcpy(result.indices, fields.items(.id));
            @memcpy(result.values, fields.items(.value));
            return result;
        }
    };
}

test "builder sorts merges duplicates and removes exact zeros" {
    const RowId = Foundation.RowId;
    var builder = SparseVectorBuilder(RowId).init(5);
    defer builder.deinit(std.testing.allocator);

    try builder.append(std.testing.allocator, RowId.fromUsizeAssumeValid(3), 2.0);
    try builder.append(std.testing.allocator, RowId.fromUsizeAssumeValid(1), 4.0);
    try builder.append(std.testing.allocator, RowId.fromUsizeAssumeValid(3), -2.0);
    try builder.append(std.testing.allocator, RowId.fromUsizeAssumeValid(2), 5.0);
    try builder.append(std.testing.allocator, RowId.fromUsizeAssumeValid(4), 0.0);

    var vector = try builder.freeze(std.testing.allocator, 0.0);
    defer vector.deinit(std.testing.allocator);

    try vector.validate();
    try std.testing.expect(vector.storage != null);
    try std.testing.expectEqual(@as(usize, 0), @intFromPtr(vector.storage.?.ptr) % 64);
    try std.testing.expectEqual(@as(usize, 2), vector.nnz());
    try std.testing.expectEqual(@as(usize, 1), vector.indices[0].toUsize());
    try std.testing.expectEqual(@as(f64, 4.0), vector.values[0]);
    try std.testing.expectEqual(@as(usize, 2), vector.indices[1].toUsize());
    try std.testing.expectEqual(@as(f64, 5.0), vector.values[1]);
}

test "builder applies an explicit drop tolerance" {
    const ColId = Foundation.ColId;
    var builder = SparseVectorBuilder(ColId).init(3);
    defer builder.deinit(std.testing.allocator);

    try builder.append(std.testing.allocator, ColId.fromUsizeAssumeValid(0), 1e-10);
    try builder.append(std.testing.allocator, ColId.fromUsizeAssumeValid(1), 2.0);

    var vector = try builder.freeze(std.testing.allocator, 1e-9);
    defer vector.deinit(std.testing.allocator);

    try vector.validate();
    try std.testing.expectEqual(@as(usize, 1), vector.nnz());
    try std.testing.expectEqual(@as(usize, 1), vector.indices[0].toUsize());
}

test "builder pre-reserved SoA append preserves duplicate order" {
    const RowId = Foundation.RowId;
    var builder = SparseVectorBuilder(RowId).init(2);
    defer builder.deinit(std.testing.allocator);

    try builder.reserve(std.testing.allocator, 4);
    builder.appendPreReserved(RowId.fromUsizeAssumeValid(1), 7.0);
    builder.appendPreReserved(RowId.fromUsizeAssumeValid(0), 1e16);
    builder.appendPreReserved(RowId.fromUsizeAssumeValid(0), -1e16);
    builder.appendPreReserved(RowId.fromUsizeAssumeValid(0), 1.0);

    var vector = try builder.freeze(std.testing.allocator, 0.0);
    defer vector.deinit(std.testing.allocator);
    try vector.validate();
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 7.0 }, vector.values);
}

test "builder validates append input and freeze tolerance" {
    const RowId = Foundation.RowId;
    var builder = SparseVectorBuilder(RowId).init(2);
    defer builder.deinit(std.testing.allocator);

    try std.testing.expectError(
        sparse.SparseVectorError.IndexOutOfBounds,
        builder.append(
            std.testing.allocator,
            RowId.fromUsizeAssumeValid(2),
            1.0,
        ),
    );
    try std.testing.expectError(
        sparse.SparseVectorError.NonFiniteValue,
        builder.append(
            std.testing.allocator,
            RowId.fromUsizeAssumeValid(0),
            std.math.inf(f64),
        ),
    );
    try std.testing.expectError(
        sparse.SparseVectorError.InvalidTolerance,
        builder.freeze(std.testing.allocator, -1.0),
    );
}
