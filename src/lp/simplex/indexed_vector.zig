//! Dense values with an optional producer-maintained index stream.
//!
//! This is the serial Zig counterpart of HiGHS `HVectorBase<Real>`.  The
//! numerical array is always addressable densely, while `state == .indexed`
//! makes `indices[0..index_count]` the ordered set of possible nonzeros.
//! Cancellation is deliberately not tracked until `tight`: factor kernels can
//! append cheaply without maintaining a second membership structure.

const std = @import("std");

pub const IndexState = enum {
    /// `indices[0..index_count]` describes every possible nonzero.
    indexed,
    /// No valid index stream exists; consumers must traverse all values.
    dense,
};

pub const IndexedVectorError = error{ DimensionMismatch, IndexOutOfBounds, IndexCapacityExceeded, DenseIndexState, Overflow };

/// Allocator-owned HiGHS-style numerical vector specialized for `Real`.
pub fn IndexedVector(comptime Real: type) type {
    return struct {
        allocator: std.mem.Allocator,
        values: []Real = &.{},
        indices: []u32 = &.{},
        index_count: usize = 0,
        state: IndexState = .indexed,

        /// UPDATE scratch retained with the vector, matching HVector ownership.
        byte_work: []u8 = &.{},
        index_work: []u32 = &.{},

        /// Packed snapshot requested by PRICE/UPDATE consumers.
        pack_requested: bool = false,
        packed_indices: []u32 = &.{},
        packed_values: []Real = &.{},
        packed_count: usize = 0,

        /// Abstract operation estimate used by HiGHS reinversion policy.
        synthetic_tick: f64 = 0.0,

        const Self = @This();
        pub const Scalar = Real;
        pub const highs_tiny: Real = 1e-14;
        pub const highs_zero: Real = 1e-50;
        pub const update_byte_padding: usize = 6400;

        /// Allocate every retained stream once for a fixed vector dimension.
        pub fn init(allocator: std.mem.Allocator, size: usize) !Self {
            if (size > std.math.maxInt(u32)) return error.Overflow;
            var self = Self{ .allocator = allocator };
            errdefer self.deinit();
            const index_work_count = std.math.mul(usize, size, 4) catch return error.Overflow;
            const byte_work_count = std.math.add(usize, size, update_byte_padding) catch return error.Overflow;
            self.values = try allocator.alloc(Real, size);
            self.indices = try allocator.alloc(u32, size);
            self.byte_work = try allocator.alloc(u8, byte_work_count);
            self.index_work = try allocator.alloc(u32, index_work_count);
            self.packed_indices = try allocator.alloc(u32, size);
            self.packed_values = try allocator.alloc(Real, size);
            @memset(self.values, 0);
            @memset(self.indices, 0);
            @memset(self.byte_work, 0);
            @memset(self.index_work, 0);
            @memset(self.packed_indices, 0);
            @memset(self.packed_values, 0);
            return self;
        }

        /// Release all vector and factor-update work storage.
        pub fn deinit(self: *Self) void {
            self.allocator.free(self.values);
            self.allocator.free(self.indices);
            self.allocator.free(self.byte_work);
            self.allocator.free(self.index_work);
            self.allocator.free(self.packed_indices);
            self.allocator.free(self.packed_values);
            self.* = .{ .allocator = self.allocator };
        }

        pub inline fn dimension(self: *const Self) usize {
            return self.values.len;
        }

        /// Borrow the valid index prefix, or null when the vector is dense.
        pub inline fn activeIndices(self: *const Self) ?[]const u32 {
            return if (self.state == .indexed) self.indices[0..self.index_count] else null;
        }

        /// Append one possible nonzero in producer-defined order.
        pub inline fn appendIndex(self: *Self, index: u32) IndexedVectorError!void {
            if (index >= self.values.len) return error.IndexOutOfBounds;
            if (self.state != .indexed) return error.DenseIndexState;
            if (self.index_count == self.indices.len) return error.IndexCapacityExceeded;
            self.indices[self.index_count] = index;
            self.index_count += 1;
        }

        /// Clear values by their valid support below 30% density, otherwise
        /// clear the full dense array, exactly matching `HVectorBase::clear`.
        pub fn clear(self: *Self) void {
            const dense_clear = self.state == .dense or self.index_count * 10 > self.values.len * 3;
            if (dense_clear) {
                @memset(self.values, 0);
            } else {
                for (self.indices[0..self.index_count]) |index| self.values[index] = 0;
            }
            self.clearScalars();
        }

        /// Reset lifecycle scalars without touching numerical storage.
        pub fn clearScalars(self: *Self) void {
            self.state = .indexed;
            self.index_count = 0;
            self.pack_requested = false;
            self.synthetic_tick = 0.0;
        }

        /// Remove entries below HiGHS' tiny threshold while preserving order.
        pub fn tight(self: *Self) void {
            if (self.state == .dense) {
                for (self.values) |*value| if (@abs(value.*) < highs_tiny) {
                    value.* = 0;
                };
                return;
            }
            var output: usize = 0;
            for (self.indices[0..self.index_count]) |index| {
                if (@abs(self.values[index]) >= highs_tiny) {
                    self.indices[output] = index;
                    output += 1;
                } else {
                    self.values[index] = 0;
                }
            }
            self.index_count = output;
        }

        /// Rebuild natural-order indices only when support is unknown or over
        /// 10% dense, preserving HiGHS' deterministic sparse ordering rule.
        pub fn reIndex(self: *Self) void {
            if (self.state == .indexed and self.index_count * 10 <= self.values.len) return;
            self.state = .indexed;
            self.index_count = 0;
            for (self.values, 0..) |value, index| if (value != 0) {
                self.indices[self.index_count] = @intCast(index);
                self.index_count += 1;
            };
        }

        /// Materialize a value/index snapshot when requested by a consumer.
        pub fn pack(self: *Self) IndexedVectorError!void {
            if (!self.pack_requested) return;
            if (self.state != .indexed) return error.DenseIndexState;
            self.pack_requested = false;
            self.packed_count = self.index_count;
            for (self.indices[0..self.index_count], 0..) |index, output| {
                self.packed_indices[output] = index;
                self.packed_values[output] = self.values[index];
            }
        }

        /// Zero the stale complement before a dense consumer reads `values`.
        /// Packed-value storage is dead until the next explicit pack request,
        /// so it serves as allocation-free preservation for active entries.
        pub fn materializeZeroComplement(self: *Self) void {
            if (self.state == .dense) return;
            for (self.indices[0..self.index_count], 0..) |index, output|
                self.packed_values[output] = self.values[index];
            @memset(self.values, 0);
            for (self.indices[0..self.index_count], 0..) |index, output|
                self.values[index] = self.packed_values[output];
            self.packed_count = 0;
        }

        /// Copy an indexed vector, converting scalar precision at compile time.
        pub fn copyFrom(self: *Self, from: anytype) IndexedVectorError!void {
            if (self.values.len != from.values.len) return error.DimensionMismatch;
            if (from.state != .indexed) return error.DenseIndexState;
            self.clear();
            self.synthetic_tick = from.synthetic_tick;
            for (from.indices[0..from.index_count]) |index| {
                try self.appendIndex(index);
                self.values[index] = @floatCast(from.values[index]);
            }
        }

        /// Squared Euclidean norm over the producer-maintained support.
        pub fn norm2(self: *const Self) IndexedVectorError!Real {
            if (self.state != .indexed) return error.DenseIndexState;
            var result: Real = 0;
            for (self.indices[0..self.index_count]) |index| {
                const value = self.values[index];
                result += value * value;
            }
            return result;
        }

        /// Add `pivot_x * pivot` without tracking cancellation, matching the
        /// exact HVector index and tiny-value semantics.
        pub fn saxpy(self: *Self, pivot_x: anytype, pivot: anytype) IndexedVectorError!void {
            if (self.values.len != pivot.values.len) return error.DimensionMismatch;
            if (self.state != .indexed or pivot.state != .indexed) return error.DenseIndexState;
            for (pivot.indices[0..pivot.index_count]) |index| {
                const old = self.values[index];
                const Product = @TypeOf(pivot_x * pivot.values[index]);
                const updated: Real = @floatCast(@as(Product, @floatCast(old)) + pivot_x * pivot.values[index]);
                if (old == 0) try self.appendIndex(index);
                self.values[index] = if (@abs(updated) < highs_tiny) highs_zero else updated;
            }
        }

        /// Bytes represented by retained slices, excluding allocator metadata.
        pub fn retainedBytes(self: *const Self) usize {
            return std.mem.sliceAsBytes(self.values).len +
                std.mem.sliceAsBytes(self.indices).len +
                std.mem.sliceAsBytes(self.byte_work).len +
                std.mem.sliceAsBytes(self.index_work).len +
                std.mem.sliceAsBytes(self.packed_indices).len +
                std.mem.sliceAsBytes(self.packed_values).len;
        }
    };
}

test "indexed vector sparse and dense clear preserve HVector lifecycle" {
    const Vector = IndexedVector(f64);
    var vector = try Vector.init(std.testing.allocator, 10);
    defer vector.deinit();

    vector.values[2] = 3;
    vector.values[7] = -4;
    try vector.appendIndex(2);
    try vector.appendIndex(7);
    vector.values[5] = 99; // Outside valid support is deliberately stale.
    vector.clear();
    try std.testing.expectEqual(@as(f64, 0), vector.values[2]);
    try std.testing.expectEqual(@as(f64, 0), vector.values[7]);
    try std.testing.expectEqual(@as(f64, 99), vector.values[5]);

    vector.state = .dense;
    vector.clear();
    for (vector.values) |value| try std.testing.expectEqual(@as(f64, 0), value);
}

test "indexed vector materializes a zero complement for dense consumers" {
    const Vector = IndexedVector(f64);
    var vector = try Vector.init(std.testing.allocator, 5);
    defer vector.deinit();
    @memset(vector.values, 99);
    vector.values[3] = -7;
    try vector.appendIndex(3);
    vector.materializeZeroComplement();
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 0, -7, 0 }, vector.values);
    try std.testing.expectEqualSlices(u32, &.{3}, vector.activeIndices().?);
}

test "indexed vector reindex pack norm and tight match HVector ordering" {
    const Vector = IndexedVector(f64);
    var vector = try Vector.init(std.testing.allocator, 20);
    defer vector.deinit();
    vector.state = .dense;
    vector.values[11] = 4;
    vector.values[3] = -2;
    vector.values[17] = 1e-16;
    vector.reIndex();
    try std.testing.expectEqualSlices(u32, &.{ 3, 11, 17 }, vector.activeIndices().?);
    try std.testing.expectApproxEqAbs(@as(f64, 20), try vector.norm2(), 1e-12);
    vector.tight();
    try std.testing.expectEqualSlices(u32, &.{ 3, 11 }, vector.activeIndices().?);
    vector.pack_requested = true;
    try vector.pack();
    try std.testing.expectEqualSlices(u32, &.{ 3, 11 }, vector.packed_indices[0..vector.packed_count]);
    try std.testing.expectEqualSlices(f64, &.{ -2, 4 }, vector.packed_values[0..vector.packed_count]);
}

test "indexed vector generic copy and saxpy retain possible nonzeros" {
    const Wide = IndexedVector(f128);
    const Narrow = IndexedVector(f64);
    var pivot = try Wide.init(std.testing.allocator, 4);
    defer pivot.deinit();
    pivot.values[1] = 2;
    pivot.values[3] = -1;
    try pivot.appendIndex(1);
    try pivot.appendIndex(3);

    var vector = try Narrow.init(std.testing.allocator, 4);
    defer vector.deinit();
    try vector.copyFrom(&pivot);
    try vector.saxpy(@as(f128, -1), &pivot);
    try std.testing.expectEqualSlices(u32, &.{ 1, 3 }, vector.activeIndices().?);
    try std.testing.expectEqual(@as(f64, Narrow.highs_zero), vector.values[1]);
    vector.tight();
    try std.testing.expectEqual(@as(usize, 0), vector.index_count);
}
