//! Reusable sparse accumulator for row/column aggregation.
//!
//! Dense values provide O(1) updates; negative zero tracks entries that
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

// Positive zero means untouched. Negative zero is numerically zero but has a
// distinct bit pattern, so it can mark an entry that was touched and then
// cancelled without stealing any finite non-zero f64 value from callers.
const sentinel: f64 = -0.0;
const sentinel_bits: u64 = @bitCast(sentinel);

/// Test the bit-level untouched-positive-zero convention.
inline fn isTouchedValue(value: f64) bool {
    return @as(u64, @bitCast(value)) != 0;
}

/// Test the negative-zero touched-but-cancelled sentinel.
inline fn isSentinel(value: f64) bool {
    return @as(u64, @bitCast(value)) == sentinel_bits;
}

/// Return a dense-plus-active sparse accumulator specialized for matrix IDs.
pub fn SparseAccumulator(comptime Id: type) type {
    const OwnedVector = sparse_vector.SparseVector(Id);
    const View = sparse_vector.SparseVectorView(Id);

    return struct {
        /// Logical dense dimension addressed by `Id`.
        dimension: usize,
        /// Dense O(1)-lookup values; positive zero denotes untouched.
        dense_values: []f64,
        /// Allocation base for separately allocated padded active IDs.
        active_base_ptr: [*]Id = undefined,
        /// First usable active-ID slot after front padding.
        active_ptr: [*]Id = undefined,
        /// Number of currently touched IDs.
        active_len: usize = 0,
        /// Usable active-ID capacity excluding padding.
        active_cap: usize = 0,
        /// Allocator retained for allocation-free hot-path ownership.
        alloc: std.mem.Allocator = undefined,
        /// Optional packed allocation containing dense and active streams.
        combined_storage: ?[]align(64) u8 = null,
        /// Whether the active-ID stream must be freed independently.
        active_separate: bool = true,

        const Self = @This();
        const active_padding = @max(1, 64 / @sizeOf(Id));

        /// Allocate a zeroed dense accumulator; active IDs grow lazily.
        pub fn init(allocator: std.mem.Allocator, dimension: usize) (std.mem.Allocator.Error || csc.MatrixError)!Self {
            if (dimension != 0) _ = Id.fromUsize(dimension - 1) catch return error.DimensionTooLarge;
            const dense_values = try allocator.alloc(f64, dimension);
            @memset(dense_values, 0);
            return .{ .dimension = dimension, .dense_values = dense_values, .alloc = allocator };
        }

        /// Constructs the dense values and active-ID capacity in one
        /// cache-colored allocation. Prefer this when the touched-set upper
        /// bound is known before entering a hot loop.
        pub fn initWithCapacity(allocator: std.mem.Allocator, dimension: usize, touched_capacity: usize) (std.mem.Allocator.Error || csc.MatrixError)!Self {
            if (dimension != 0) _ = Id.fromUsize(dimension - 1) catch return error.DimensionTooLarge;
            const dense_bytes = std.math.mul(usize, dimension, @sizeOf(f64)) catch return error.DimensionTooLarge;
            const active_bytes = std.math.mul(usize, touched_capacity, @sizeOf(Id)) catch return error.DimensionTooLarge;
            const active_offset = std.mem.alignForward(usize, dense_bytes, 4096) + 64;
            const storage_len = std.math.add(usize, active_offset, active_bytes) catch return error.DimensionTooLarge;
            const storage = try allocator.alignedAlloc(u8, .@"64", storage_len);
            const dense_ptr: [*]f64 = @ptrCast(@alignCast(storage.ptr));
            const active_ptr: [*]Id = @ptrCast(@alignCast(storage.ptr + active_offset));
            const dense_values = dense_ptr[0..dimension];
            memory.clearF64(dense_values);
            return .{
                .dimension = dimension,
                .dense_values = dense_values,
                .active_base_ptr = active_ptr,
                .active_ptr = active_ptr,
                .active_cap = touched_capacity,
                .alloc = allocator,
                .combined_storage = storage,
                .active_separate = false,
            };
        }

        /// Release packed or separate dense/active storage.
        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            if (self.combined_storage) |storage| allocator.free(storage) else allocator.free(self.dense_values);
            if (self.active_separate and self.active_cap > 0)
                allocator.free(self.active_base_ptr[0 .. self.active_cap + active_padding]);
            self.* = undefined;
        }

        /// Number of unique IDs touched since the last clear.
        pub inline fn touchedCount(self: Self) usize {
            return self.active_len;
        }

        /// Ensure the active-ID list can hold at least `cap` unique touches.
        pub fn reserve(self: *Self, allocator: std.mem.Allocator, cap: usize) std.mem.Allocator.Error!void {
            if (cap > self.active_cap) {
                const buf = try allocator.alloc(Id, cap + active_padding);
                if (self.active_cap > 0) {
                    @memcpy(buf[active_padding..][0..self.active_len], self.active_ptr[0..self.active_len]);
                    if (self.active_separate)
                        allocator.free(self.active_base_ptr[0 .. self.active_cap + active_padding]);
                }
                self.active_base_ptr = buf.ptr;
                self.active_ptr = buf.ptr + active_padding;
                self.active_cap = cap;
                self.active_separate = true;
            }
        }

        /// Reset touched values, choosing sparse or dense clearing by density.
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

        /// Checked addition to one coordinate, growing the active list if needed.
        pub fn add(self: *Self, allocator: std.mem.Allocator, id: Id, value: f64) (std.mem.Allocator.Error || csc.MatrixError)!void {
            const index = id.toUsize();
            if (index >= self.dimension) return error.IndexOutOfBounds;
            if (!std.math.isFinite(value)) return error.NonFiniteValue;
            if (value == 0.0) return;

            if (isTouchedValue(self.dense_values[index])) {
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
            if (isTouchedValue(dv[index])) {
                dv[index] += value;
                if (dv[index] == 0.0) dv[index] = sentinel;
            } else {
                dv[index] = value;
                self.active_ptr[self.active_len] = id;
                self.active_len += 1;
            }
        }

        /// Accumulate `alpha * vector` after validating dimensions and values.
        pub fn addVector(self: *Self, allocator: std.mem.Allocator, alpha: f64, vector: View) (std.mem.Allocator.Error || csc.MatrixError)!void {
            if (vector.dimension != self.dimension) return error.DimensionMismatch;
            if (!std.math.isFinite(alpha)) return error.NonFiniteValue;
            try vector.validate();
            for (vector.indices, vector.values) |id, value| try self.add(allocator, id, alpha * value);
        }

        /// Read one coordinate, translating the touched-zero sentinel to `0.0`.
        pub inline fn get(self: Self, id: Id) f64 {
            const index = id.toUsize();
            std.debug.assert(index < self.dimension);
            const value = self.dense_values[index];
            return if (isSentinel(value)) 0.0 else value;
        }

        /// Sort touched IDs and copy retained values into an owned canonical vector.
        pub fn freeze(self: *Self, allocator: std.mem.Allocator, zero_tolerance: f64) (std.mem.Allocator.Error || csc.MatrixError)!OwnedVector {
            if (!std.math.isFinite(zero_tolerance) or zero_tolerance < 0.0) return error.InvalidTolerance;
            const slice = self.active_ptr[0..self.active_len];
            std.sort.pdq(Id, slice, {}, lessThanId);
            var count: usize = 0;
            for (slice) |id| {
                if (!isSentinel(self.dense_values[@intFromEnum(id)]) and
                    @abs(self.dense_values[@intFromEnum(id)]) > zero_tolerance) count += 1;
            }
            const indices = try allocator.alloc(Id, count);
            errdefer allocator.free(indices);
            const vs = try allocator.alloc(f64, count);
            errdefer allocator.free(vs);
            var dst: usize = 0;
            for (slice) |id| {
                const v = self.dense_values[@intFromEnum(id)];
                if (isSentinel(v) or @abs(v) <= zero_tolerance) continue;
                indices[dst] = id;
                vs[dst] = v;
                dst += 1;
            }
            return OwnedVector.initOwnedSlicesAssumeValid(self.dimension, indices, vs);
        }

        /// Sort active IDs into canonical ascending order.
        fn lessThanId(_: void, lhs: Id, rhs: Id) bool {
            return lhs.toUsize() < rhs.toUsize();
        }

        /// Geometrically grow the padded active-ID list.
        fn growActive(allocator: std.mem.Allocator, self: *Self) std.mem.Allocator.Error!void {
            const new_cap = @max(self.active_cap + 16, self.active_cap * 2);
            const buf = try allocator.alloc(Id, new_cap + active_padding);
            if (self.active_cap > 0) {
                @memcpy(buf[active_padding..][0..self.active_len], self.active_ptr[0..self.active_len]);
                if (self.active_separate)
                    allocator.free(self.active_base_ptr[0 .. self.active_cap + active_padding]);
            }
            self.active_base_ptr = buf.ptr;
            self.active_ptr = buf.ptr + active_padding;
            self.active_cap = new_cap;
            self.active_separate = true;
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

test "sparse accumulator preserves the smallest positive normal value" {
    const Accumulator = SparseAccumulator(foundation.RowId);
    var sum = try Accumulator.init(std.testing.allocator, 1);
    defer sum.deinit(std.testing.allocator);
    const id = try foundation.RowId.init(0);
    try sum.add(std.testing.allocator, id, std.math.floatMin(f64));
    try std.testing.expectEqual(std.math.floatMin(f64), sum.get(id));
    var vector = try sum.freeze(std.testing.allocator, 0.0);
    defer vector.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 1), vector.nnz());
    try std.testing.expectEqual(std.math.floatMin(f64), vector.values[0]);
}

test "sparse accumulator adds canonical vectors" {
    const Accumulator = SparseAccumulator(foundation.RowId);
    var sum = try Accumulator.init(std.testing.allocator, 3);
    defer sum.deinit(std.testing.allocator);
    var ids = [_]foundation.RowId{ try foundation.RowId.init(0), try foundation.RowId.init(2) };
    var vals = [_]f64{ 2.0, -3.0 };
    try sum.addVector(std.testing.allocator, -2.0, sparse_vector.SparseVectorView(foundation.RowId).initAssumeValid(3, &ids, &vals));
    try std.testing.expectEqual(@as(f64, -4.0), sum.get(ids[0]));
    try std.testing.expectEqual(@as(f64, 6.0), sum.get(ids[1]));
}
