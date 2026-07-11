//! Reusable sparse accumulator for row/column aggregation.
//!
//! Dense values provide O(1) updates; a soft sentinel tracks entries that
//! cancel to exactly zero; active IDs keep extraction proportional to the
//! touched set. This design matches HiGHS HighsSparseVectorSum.
//!
//! Active list uses raw ptr/cap/len + inline hot-path for minimal overhead.
//! The `addAssumeValid` function is marked `inline` to eliminate function
//! call overhead (~40 instructions of prologue/epilogue/stack-spill that
//! LLVM otherwise inserts) and to enable cross-call optimization across
//! paired "add first touch / add accumulate" call sites.

const std = @import("std");
const foundation = @import("foundation");
const sparse_vector = @import("sparse_vector.zig");
const csc = @import("csc.zig");
const memory = @import("memory.zig");

const sentinel: f64 = std.math.floatMin(f64);

pub fn SparseAccumulator(comptime Id: type) type {
    const OwnedVector = sparse_vector.SparseVector(Id);
    const View = sparse_vector.SparseVectorView(Id);

    return struct {
        dimension: usize,
        dense_values: []f64,
        active_ptr: [*]Id = undefined,
        active_len: usize = 0,
        active_cap: usize = 0,
        alloc: std.mem.Allocator = undefined,

        const Self = @This();

        pub fn init(allocator: std.mem.Allocator, dimension: usize) (std.mem.Allocator.Error || csc.MatrixError)!Self {
            if (dimension != 0) _ = Id.fromUsize(dimension - 1) catch return error.DimensionTooLarge;
            const dense_values = try allocator.alloc(f64, dimension);
            @memset(dense_values, 0);
            return .{ .dimension = dimension, .dense_values = dense_values, .alloc = allocator };
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.dense_values);
            if (self.active_cap > 0) allocator.free(self.active_ptr[0..self.active_cap]);
            self.* = undefined;
        }

        pub inline fn touchedCount(self: Self) usize { return self.active_len; }

        pub fn reserve(self: *Self, allocator: std.mem.Allocator, cap: usize) std.mem.Allocator.Error!void {
            if (cap > self.active_cap) {
                const buf = try allocator.alloc(Id, cap);
                if (self.active_cap > 0) {
                    @memcpy(buf[0..self.active_len], self.active_ptr[0..self.active_len]);
                    allocator.free(self.active_ptr[0..self.active_cap]);
                }
                self.active_ptr = buf.ptr;
                self.active_cap = cap;
            }
        }

        pub fn clear(self: *Self) void {
            if (10 * self.active_len < 3 * self.dimension) {
                var i: usize = 0;
                while (i < self.active_len) : (i += 1)
                    self.dense_values[@intFromEnum(self.active_ptr[i])] = 0.0;
            } else {
                memory.clearF64(self.dense_values);
            }
            self.active_len = 0;
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
                if (self.active_len == self.active_cap) try growActive(allocator, self);
                self.dense_values[index] = value;
                self.active_ptr[self.active_len] = id;
                self.active_len += 1;
            }
        }

        /// Allocation-free update after reserve. Inline eliminates function
        /// call overhead (~40 stack spill instructions) and enables LLVM to
        /// keep dense_values ptr, active_ptr, active_len in registers across
        /// paired touch/accumulate call sites.
        pub inline fn addAssumeValid(self: *Self, id: Id, value: f64) void {
            if (value == 0.0) return;
            const index = id.toUsize();
            const dv = self.dense_values;
            if (dv[index] != 0.0) {
                dv[index] += value;
                if (dv[index] == 0.0) dv[index] = sentinel;
            } else {
                dv[index] = value;
                self.active_ptr[self.active_len] = id;
                self.active_len += 1;
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
            const value = self.dense_values[index];
            return if (value == 0.0 or value == sentinel) 0.0 else value;
        }

        pub fn freeze(self: *Self, allocator: std.mem.Allocator, zero_tolerance: f64) (std.mem.Allocator.Error || csc.MatrixError)!OwnedVector {
            if (!std.math.isFinite(zero_tolerance) or zero_tolerance < 0.0) return error.InvalidTolerance;
            const slice = self.active_ptr[0..self.active_len];
            std.sort.pdq(Id, slice, {}, lessThanId);
            var count: usize = 0;
            for (slice) |id| {
                if (self.dense_values[@intFromEnum(id)] != sentinel and
                    @abs(self.dense_values[@intFromEnum(id)]) > zero_tolerance) count += 1;
            }
            const indices = try allocator.alloc(Id, count);
            errdefer allocator.free(indices);
            const vs = try allocator.alloc(f64, count);
            errdefer allocator.free(vs);
            var dst: usize = 0;
            for (slice) |id| {
                const v = self.dense_values[@intFromEnum(id)];
                if (v == sentinel or @abs(v) <= zero_tolerance) continue;
                indices[dst] = id;
                vs[dst] = v;
                dst += 1;
            }
            return .{ .dimension = self.dimension, .indices = indices, .values = vs };
        }

        fn lessThanId(_: void, lhs: Id, rhs: Id) bool { return lhs.toUsize() < rhs.toUsize(); }

        fn growActive(allocator: std.mem.Allocator, self: *Self) std.mem.Allocator.Error!void {
            const new_cap = @max(self.active_cap + 16, self.active_cap * 2);
            const buf = try allocator.alloc(Id, new_cap);
            if (self.active_cap > 0) {
                @memcpy(buf[0..self.active_len], self.active_ptr[0..self.active_len]);
                allocator.free(self.active_ptr[0..self.active_cap]);
            }
            self.active_ptr = buf.ptr;
            self.active_cap = new_cap;
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

test "sparse accumulator handles sentinel and dense clear" {
    const Accumulator = SparseAccumulator(foundation.RowId);
    var sum = try Accumulator.init(std.testing.allocator, 2);
    defer sum.deinit(std.testing.allocator);
    try sum.add(std.testing.allocator, try foundation.RowId.init(0), 1e10);
    try sum.add(std.testing.allocator, try foundation.RowId.init(0), -1e10);
    try std.testing.expectEqual(@as(f64, 0.0), sum.get(try foundation.RowId.init(0)));
    sum.clear();
    try sum.add(std.testing.allocator, try foundation.RowId.init(0), 7.0);
    try std.testing.expectEqual(@as(f64, 7.0), sum.get(try foundation.RowId.init(0)));
    sum.clear();
    try std.testing.expectEqual(@as(f64, 0.0), sum.get(try foundation.RowId.init(0)));
}

test "sparse accumulator adds canonical vectors" {
    const Accumulator = SparseAccumulator(foundation.RowId);
    var sum = try Accumulator.init(std.testing.allocator, 3);
    defer sum.deinit(std.testing.allocator);
    var ids = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(2) };
    var vals = [_]f64{ 2.0, -3.0 };
    try sum.addVector(std.testing.allocator, -2.0, .{ .dimension = 3, .indices = &ids, .values = &vals });
    try std.testing.expectEqual(@as(f64, -4.0), sum.get(ids[0]));
    try std.testing.expectEqual(@as(f64, 6.0), sum.get(ids[1]));
}
