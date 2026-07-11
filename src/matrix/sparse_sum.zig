//! Reusable sparse accumulator for row/column aggregation.
//!
//! Dense values provide O(1) updates; a soft sentinel tracks entries that
//! cancel to exactly zero; active IDs keep extraction proportional to the
//! touched set. This design matches HiGHS HighsSparseVectorSum: instead of a
//! separate marks array with generation counters, we check the value itself --
//! untouched entries are zero, touched entries are non-zero. This eliminates
//! one memory load and one store per add() compared to a generation-based
//! scheme, at the cost of needing to zero touched entries in clear().
//!
//! The sentinel (floatMin) handles the cancellation edge case: when a
//! sequence of adds produces exactly 0.0, we store floatMin so the active
//! list entry remains valid; freeze() and get() both recognize the sentinel
//! as effective zero.
//!
//! Intended for presolve aggregation and cut construction, not canonical
//! persistent storage.

const std = @import("std");
const foundation = @import("foundation");
const sparse_vector = @import("sparse_vector.zig");
const csc = @import("csc.zig");

/// Sentinel value stored when an accumulated entry cancels to exactly zero.
/// floatMin (~2.2e-308) is far below any practical LP value and serves as a
/// reliable "was-touched-but-now-zero" marker that preserves the active-list
/// entry without requiring a separate marks array.
const sentinel: f64 = std.math.floatMin(f64);

pub fn SparseAccumulator(comptime Id: type) type {
    // Instantiation validates that Id is RowId or ColId.
    const OwnedVector = sparse_vector.SparseVector(Id);
    const View = sparse_vector.SparseVectorView(Id);

    return struct {
        dimension: usize,
        dense_values: []f64,
        active: std.ArrayList(Id) = .empty,

        const Self = @This();

        pub fn init(allocator: std.mem.Allocator, dimension: usize) (std.mem.Allocator.Error || csc.MatrixError)!Self {
            if (dimension != 0) _ = Id.fromUsize(dimension - 1) catch return error.DimensionTooLarge;
            const dense_values = try allocator.alloc(f64, dimension);
            errdefer allocator.free(dense_values);
            @memset(dense_values, 0);
            return .{ .dimension = dimension, .dense_values = dense_values };
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.dense_values);
            self.active.deinit(allocator);
            self.* = undefined;
        }

        pub inline fn touchedCount(self: Self) usize {
            return self.active.items.len;
        }

        pub fn reserve(self: *Self, allocator: std.mem.Allocator, touched_capacity: usize) std.mem.Allocator.Error!void {
            try self.active.ensureTotalCapacity(allocator, touched_capacity);
        }

        /// Clears all touched entries using a sparse-vs-dense heuristic.
        ///
        /// If fewer than 30% of entries were touched, iterate only the active
        /// list; otherwise memset the entire dense array. This avoids the
        /// per-add overhead of a generation/marks scheme while keeping clear()
        /// proportional to the smaller of {touched, dimension}.
        pub fn clear(self: *Self) void {
            if (10 * self.active.items.len < 3 * self.dimension) {
                for (self.active.items) |id| self.dense_values[id.toUsize()] = 0.0;
            } else {
                @memset(self.dense_values, 0);
            }
            self.active.clearRetainingCapacity();
        }

        pub fn add(self: *Self, allocator: std.mem.Allocator, id: Id, value: f64) (std.mem.Allocator.Error || csc.MatrixError)!void {
            const index = id.toUsize();
            if (index >= self.dimension) return error.IndexOutOfBounds;
            if (!std.math.isFinite(value)) return error.NonFiniteValue;
            if (value == 0.0) return;

            if (self.dense_values[index] != 0.0) {
                const sum = self.dense_values[index] + value;
                if (!std.math.isFinite(sum)) return error.NonFiniteValue;
                self.dense_values[index] = if (sum == 0.0) sentinel else sum;
            } else {
                // Reserve before changing value so OOM leaves no hidden
                // touched entry missing from the active list.
                try self.active.ensureUnusedCapacity(allocator, 1);
                self.dense_values[index] = value;
                self.active.appendAssumeCapacity(id);
            }
        }

        /// Allocation-free update after reserve. Caller guarantees a valid ID,
        /// finite value/sum, and enough capacity for a newly touched index.
        pub fn addAssumeValid(self: *Self, id: Id, value: f64) void {
            if (value == 0.0) return;
            const index = id.toUsize();
            const dv = self.dense_values;
            if (dv[index] != 0.0) {
                dv[index] += value;
                if (dv[index] == 0.0) dv[index] = sentinel;
            } else {
                dv[index] = value;
                const active_items = &self.active.items;
                active_items.ptr[active_items.len] = id;
                active_items.len += 1;
            }
        }

        pub fn addVector(self: *Self, allocator: std.mem.Allocator, alpha: f64, vector: View) (std.mem.Allocator.Error || csc.MatrixError)!void {
            if (vector.dimension != self.dimension) return error.DimensionMismatch;
            if (!std.math.isFinite(alpha)) return error.NonFiniteValue;
            try vector.validate();
            for (vector.indices, vector.values) |id, value| try self.add(allocator, id, alpha * value);
        }

        /// Returns the current accumulated value. Untouched and cancelled-to-zero
        /// entries both return 0.0.
        pub inline fn get(self: Self, id: Id) f64 {
            const index = id.toUsize();
            std.debug.assert(index < self.dimension);
            const value = self.dense_values[index];
            return if (value == 0.0 or value == sentinel) 0.0 else value;
        }

        /// Extracts sorted, unique canonical storage and removes small results.
        /// Sentinel entries (cancelled to zero) are always discarded, even when
        /// zero_tolerance is 0.
        pub fn freeze(self: *Self, allocator: std.mem.Allocator, zero_tolerance: f64) (std.mem.Allocator.Error || csc.MatrixError)!OwnedVector {
            if (!std.math.isFinite(zero_tolerance) or zero_tolerance < 0.0) return error.InvalidTolerance;
            std.sort.pdq(Id, self.active.items, {}, lessThanId);
            var count: usize = 0;
            for (self.active.items) |id| {
                const value = self.dense_values[id.toUsize()];
                if (value != sentinel and @abs(value) > zero_tolerance) count += 1;
            }
            const indices = try allocator.alloc(Id, count);
            errdefer allocator.free(indices);
            const values = try allocator.alloc(f64, count);
            errdefer allocator.free(values);
            var destination: usize = 0;
            for (self.active.items) |id| {
                const value = self.dense_values[id.toUsize()];
                if (value == sentinel or @abs(value) <= zero_tolerance) continue;
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

test "sparse accumulator recovers after sentinel cancellation" {
    const Accumulator = SparseAccumulator(foundation.RowId);
    var sum = try Accumulator.init(std.testing.allocator, 2);
    defer sum.deinit(std.testing.allocator);

    // Add and cancel to zero → triggers sentinel
    try sum.add(std.testing.allocator, try foundation.RowId.init(0), 1e10);
    try sum.add(std.testing.allocator, try foundation.RowId.init(0), -1e10);
    try std.testing.expectEqual(@as(f64, 0.0), sum.get(try foundation.RowId.init(0)));

    // Clear and re-add
    sum.clear();
    try sum.add(std.testing.allocator, try foundation.RowId.init(0), 7.0);
    try std.testing.expectEqual(@as(f64, 7.0), sum.get(try foundation.RowId.init(0)));

    // Dense clear path (2 of 2 entries touched = 100%)
    sum.clear();
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
