//! Reusable sparse accumulator for row/column aggregation.
//!
//! Dense values provide O(1) updates; generation marks avoid clearing the full
//! dense arrays between sums; active IDs keep extraction proportional to the
//! touched set. This is intended for presolve aggregation and cut construction,
//! not canonical persistent storage.

const std = @import("std");
const foundation = @import("foundation");
const sparse_vector = @import("sparse_vector.zig");
const csc = @import("csc.zig");

pub fn SparseAccumulator(comptime Id: type) type {
    // Instantiation validates that Id is RowId or ColId.
    const OwnedVector = sparse_vector.SparseVector(Id);
    const View = sparse_vector.SparseVectorView(Id);

    return struct {
        dimension: usize,
        dense_values: []f64,
        marks: []u32,
        active: std.ArrayList(Id) = .empty,
        generation: u32 = 1,

        const Self = @This();

        pub fn init(allocator: std.mem.Allocator, dimension: usize) (std.mem.Allocator.Error || csc.MatrixError)!Self {
            if (dimension != 0) _ = Id.fromUsize(dimension - 1) catch return error.DimensionTooLarge;
            const dense_values = try allocator.alloc(f64, dimension);
            errdefer allocator.free(dense_values);
            const marks = try allocator.alloc(u32, dimension);
            errdefer allocator.free(marks);
            @memset(marks, 0);
            return .{ .dimension = dimension, .dense_values = dense_values, .marks = marks };
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.dense_values);
            allocator.free(self.marks);
            self.active.deinit(allocator);
            self.* = undefined;
        }

        pub inline fn touchedCount(self: Self) usize {
            return self.active.items.len;
        }

        pub fn reserve(self: *Self, allocator: std.mem.Allocator, touched_capacity: usize) std.mem.Allocator.Error!void {
            try self.active.ensureTotalCapacity(allocator, touched_capacity);
        }

        /// Logical clear is O(1) except once per 2^32 generations.
        pub fn clear(self: *Self) void {
            self.active.clearRetainingCapacity();
            self.generation +%= 1;
            if (self.generation == 0) {
                @memset(self.marks, 0);
                self.generation = 1;
            }
        }

        pub fn add(self: *Self, allocator: std.mem.Allocator, id: Id, value: f64) (std.mem.Allocator.Error || csc.MatrixError)!void {
            const index = id.toUsize();
            if (index >= self.dimension) return error.IndexOutOfBounds;
            if (!std.math.isFinite(value)) return error.NonFiniteValue;
            if (value == 0.0) return;

            if (self.marks[index] != self.generation) {
                // Reserve before changing marks/value so OOM leaves no hidden
                // touched entry missing from the active list.
                try self.active.ensureUnusedCapacity(allocator, 1);
                self.marks[index] = self.generation;
                self.dense_values[index] = value;
                self.active.appendAssumeCapacity(id);
            } else {
                const sum = self.dense_values[index] + value;
                if (!std.math.isFinite(sum)) return error.NonFiniteValue;
                self.dense_values[index] = sum;
            }
        }

        /// Allocation-free update after reserve. Caller guarantees a valid ID,
        /// finite value/sum, and enough capacity for a newly touched index.
        pub fn addAssumeValid(self: *Self, id: Id, value: f64) void {
            if (value == 0.0) return;
            const index = id.toUsize();
            if (self.marks[index] != self.generation) {
                self.marks[index] = self.generation;
                self.dense_values[index] = value;
                self.active.appendAssumeCapacity(id);
            } else {
                self.dense_values[index] += value;
            }
        }

        pub fn addVector(self: *Self, allocator: std.mem.Allocator, alpha: f64, vector: View) (std.mem.Allocator.Error || csc.MatrixError)!void {
            if (vector.dimension != self.dimension) return error.DimensionMismatch;
            if (!std.math.isFinite(alpha)) return error.NonFiniteValue;
            try vector.validate();
            for (vector.indices, vector.values) |id, value| try self.add(allocator, id, alpha * value);
        }

        pub inline fn get(self: Self, id: Id) f64 {
            const index = id.toUsize();
            std.debug.assert(index < self.dimension);
            return if (self.marks[index] == self.generation) self.dense_values[index] else 0.0;
        }

        /// Extracts sorted, unique canonical storage and removes small results.
        pub fn freeze(self: *Self, allocator: std.mem.Allocator, zero_tolerance: f64) (std.mem.Allocator.Error || csc.MatrixError)!OwnedVector {
            if (!std.math.isFinite(zero_tolerance) or zero_tolerance < 0.0) return error.InvalidTolerance;
            std.sort.pdq(Id, self.active.items, {}, lessThanId);
            var count: usize = 0;
            for (self.active.items) |id| if (@abs(self.dense_values[id.toUsize()]) > zero_tolerance) {
                count += 1;
            };
            const indices = try allocator.alloc(Id, count);
            errdefer allocator.free(indices);
            const values = try allocator.alloc(f64, count);
            errdefer allocator.free(values);
            var destination: usize = 0;
            for (self.active.items) |id| {
                const value = self.dense_values[id.toUsize()];
                if (@abs(value) <= zero_tolerance) continue;
                indices[destination] = id;
                values[destination] = value;
                destination += 1;
            }
            return .{ .dimension = self.dimension, .indices = indices, .values = values };
        }

        fn lessThanId(_: void, lhs: Id, rhs: Id) bool {
            return lhs.toUsize() < rhs.toUsize();
        }
    };
}

test "sparse accumulator merges cancels sorts and clears" {
    const Accumulator = SparseAccumulator(foundation.ColId);
    var sum = try Accumulator.init(std.testing.allocator, 5);
    defer sum.deinit(std.testing.allocator);
    try sum.add(std.testing.allocator, try foundation.ColId.init(3), 4.0);
    try sum.add(std.testing.allocator, try foundation.ColId.init(1), 2.0);
    try sum.add(std.testing.allocator, try foundation.ColId.init(3), -4.0);
    try std.testing.expectEqual(@as(usize, 2), sum.touchedCount());
    var vector = try sum.freeze(std.testing.allocator, 0.0);
    defer vector.deinit(std.testing.allocator);
    try vector.validate();
    try std.testing.expectEqual(@as(usize, 1), vector.nnz());
    try std.testing.expectEqual(@as(usize, 1), vector.indices[0].toUsize());
    sum.clear();
    try std.testing.expectEqual(@as(f64, 0.0), sum.get(try foundation.ColId.init(1)));
}

test "sparse accumulator generation wrap resets marks" {
    const Accumulator = SparseAccumulator(foundation.RowId);
    var sum = try Accumulator.init(std.testing.allocator, 2);
    defer sum.deinit(std.testing.allocator);
    try sum.add(std.testing.allocator, try foundation.RowId.init(0), 7.0);
    sum.generation = std.math.maxInt(u32);
    sum.clear();
    try std.testing.expectEqual(@as(u32, 1), sum.generation);
    try std.testing.expectEqual(@as(f64, 0.0), sum.get(try foundation.RowId.init(0)));
}

test "sparse accumulator adds canonical vectors" {
    const Accumulator = SparseAccumulator(foundation.RowId);
    var sum = try Accumulator.init(std.testing.allocator, 3);
    defer sum.deinit(std.testing.allocator);
    var ids = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(2) };
    var values = [_]f64{ 2.0, -3.0 };
    try sum.addVector(std.testing.allocator, -2.0, .{ .dimension = 3, .indices = &ids, .values = &values });
    try std.testing.expectEqual(@as(f64, -4.0), sum.get(ids[0]));
    try std.testing.expectEqual(@as(f64, 6.0), sum.get(ids[1]));
}
