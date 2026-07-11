//! Mutable construction and canonicalization of sparse vectors.

const std = @import("std");
const Foundation = @import("foundation");
const sparse = @import("sparse_vector.zig");

pub fn SparseVectorBuilder(comptime Id: type) type {
    // Instantiating the output type performs the RowId/ColId validation.
    const Vector = sparse.SparseVector(Id);

    return struct {
        dimension: usize,
        entries: std.ArrayList(Entry) = .empty,

        const Self = @This();

        const Entry = struct {
            id: Id,
            value: f64,

            fn lessThan(_: void, lhs: Entry, rhs: Entry) bool {
                return lhs.id.raw() < rhs.id.raw();
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
            try self.entries.append(allocator, .{ .id = id, .value = value });
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
            if (self.entries.items.len == 0) {
                return Vector.empty(self.dimension);
            }

            // Stable ordering preserves append order among duplicate IDs, so
            // floating-point merging is deterministic with respect to input.
            std.sort.block(Entry, self.entries.items, {}, Entry.lessThan);

            var read: usize = 0;
            var write: usize = 0;
            while (read < self.entries.items.len) {
                const id = self.entries.items[read].id;
                var sum: f64 = 0.0;

                while (read < self.entries.items.len and
                    self.entries.items[read].id.raw() == id.raw()) : (read += 1)
                {
                    sum += self.entries.items[read].value;
                }

                if (!std.math.isFinite(sum)) {
                    return sparse.SparseVectorError.NonFiniteValue;
                }
                if (@abs(sum) <= drop_tolerance) continue;

                self.entries.items[write] = .{ .id = id, .value = sum };
                write += 1;
            }
            self.entries.items.len = write;

            if (write == 0) return Vector.empty(self.dimension);

            const indices = try allocator.alloc(Id, write);
            errdefer allocator.free(indices);
            const values = try allocator.alloc(f64, write);
            errdefer allocator.free(values);

            for (self.entries.items, 0..) |entry, position| {
                indices[position] = entry.id;
                values[position] = entry.value;
            }

            return .{
                .dimension = self.dimension,
                .indices = indices,
                .values = values,
            };
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
