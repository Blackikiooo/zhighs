//! Canonical sparse-vector storage and zero-allocation borrowed views.
//!
//! This module represents matrix rows and columns as sorted `(index, value)`
//! pairs. It is intentionally different from the dense-plus-active workspace
//! used later by Simplex FTRAN/BTRAN.

const std = @import("std");
const Foundation = @import("foundation");

pub const SparseVectorError = error{
    /// The vector has a nonzero value at an index outside its dimension.
    DimensionMismatch,
    /// The vector contains an explicit zero value, which is not allowed in canonical form.
    ExplicitZero,
    /// The vector contains an index outside its dimension.
    IndexOutOfBounds,
    /// The vector contains indices that are not strictly increasing.
    IndicesNotStrictlyIncreasing,
    /// The tolerance value is invalid.
    InvalidTolerance,
    /// The lengths of the index and value slices do not match.
    LengthMismatch,
    /// The vector contains a non-finite value.
    NonFiniteValue,
};

/// compile-time requirement that the index type is either `RowId` or `ColId`.
fn requireMatrixId(comptime Id: type) void {
    if (Id != Foundation.RowId and Id != Foundation.ColId) {
        @compileError("SparseVector index type must be RowId or ColId");
    }
}

/// A borrowed canonical sparse vector.
///
/// Valid views have equally sized index/value slices, strictly increasing and
/// unique indices, finite nonzero values, and indices inside `dimension`.
pub fn SparseVectorView(comptime Id: type) type {
    requireMatrixId(Id);

    return struct {
        dimension: usize,
        indices: []const Id,
        values: []const f64,

        const Self = @This();

        pub inline fn nnz(self: Self) usize {
            return self.indices.len;
        }

        pub inline fn isEmpty(self: Self) bool {
            return self.indices.len == 0;
        }

        pub fn validate(self: Self) SparseVectorError!void {
            if (self.indices.len != self.values.len) {
                return SparseVectorError.LengthMismatch;
            }

            var previous_raw: ?Foundation.HUInt = null;

            for (self.indices, self.values) |id, value| {
                const position = id.toUsize();
                if (position >= self.dimension) {
                    return SparseVectorError.IndexOutOfBounds;
                }
                if (!std.math.isFinite(value)) {
                    return SparseVectorError.NonFiniteValue;
                }
                if (value == 0.0) {
                    return SparseVectorError.ExplicitZero;
                }

                const raw_id = id.raw();
                if (previous_raw) |previous| {
                    if (raw_id <= previous) {
                        return SparseVectorError.IndicesNotStrictlyIncreasing;
                    }
                }
                previous_raw = raw_id;
            }
        }

        /// Computes a sparse/dense dot product after checking dimensions.
        pub fn dotDense(self: Self, dense: []const f64) SparseVectorError!f64 {
            if (dense.len != self.dimension) {
                return SparseVectorError.DimensionMismatch;
            }
            return self.dotDenseAssumeValid(dense);
        }

        /// Computes a sparse/dense dot product for an already validated view.
        pub fn dotDenseAssumeValid(self: Self, dense: []const f64) f64 {
            std.debug.assert(self.indices.len == self.values.len);
            std.debug.assert(dense.len == self.dimension);

            var result: f64 = 0.0;
            for (self.indices, self.values) |id, value| {
                const position = id.toUsize();
                std.debug.assert(position < dense.len);
                result += value * dense[position];
            }
            return result;
        }

        /// Adds `alpha * self` to a dense vector after checking dimensions.
        pub fn addToDense(
            self: Self,
            alpha: f64,
            dense: []f64,
        ) SparseVectorError!void {
            if (dense.len != self.dimension) {
                return SparseVectorError.DimensionMismatch;
            }
            self.addToDenseAssumeValid(alpha, dense);
        }

        /// Adds `alpha * self` to an already dimension-compatible dense vector.
        pub fn addToDenseAssumeValid(self: Self, alpha: f64, dense: []f64) void {
            std.debug.assert(self.indices.len == self.values.len);
            std.debug.assert(dense.len == self.dimension);

            for (self.indices, self.values) |id, value| {
                const position = id.toUsize();
                std.debug.assert(position < dense.len);
                dense[position] += alpha * value;
            }
        }
    };
}

/// An allocator-owned canonical sparse vector.
pub fn SparseVector(comptime Id: type) type {
    requireMatrixId(Id);

    return struct {
        dimension: usize,
        indices: []Id,
        values: []f64,

        const Self = @This();
        pub const View = SparseVectorView(Id);

        pub fn empty(dimension: usize) Self {
            return .{
                .dimension = dimension,
                .indices = &.{},
                .values = &.{},
            };
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.indices);
            allocator.free(self.values);
            self.* = empty(0);
        }

        pub inline fn nnz(self: *const Self) usize {
            return self.indices.len;
        }

        pub inline fn view(self: *const Self) View {
            return .{
                .dimension = self.dimension,
                .indices = self.indices,
                .values = self.values,
            };
        }

        pub fn validate(self: *const Self) SparseVectorError!void {
            return self.view().validate();
        }
    };
}

test "canonical RowId sparse-vector view validates and computes" {
    const RowId = Foundation.RowId;
    const View = SparseVectorView(RowId);
    const indices = [_]RowId{
        RowId.fromUsizeAssumeValid(1),
        RowId.fromUsizeAssumeValid(3),
    };
    const values = [_]f64{ 2.0, -4.0 };
    const vector: View = .{
        .dimension = 5,
        .indices = &indices,
        .values = &values,
    };

    try vector.validate();
    try std.testing.expectEqual(@as(usize, 2), vector.nnz());
    try std.testing.expect(!vector.isEmpty());

    const dense = [_]f64{ 10.0, 3.0, 20.0, 5.0, 30.0 };
    try std.testing.expectEqual(@as(f64, -14.0), try vector.dotDense(&dense));

    var target = [_]f64{ 0.0, 1.0, 0.0, 2.0, 0.0 };
    try vector.addToDense(0.5, &target);
    try std.testing.expectEqualSlices(
        f64,
        &[_]f64{ 0.0, 2.0, 0.0, 0.0, 0.0 },
        &target,
    );
}

test "canonical ColId sparse-vector view is supported" {
    const ColId = Foundation.ColId;
    const indices = [_]ColId{ColId.fromUsizeAssumeValid(2)};
    const values = [_]f64{7.0};
    const vector: SparseVectorView(ColId) = .{
        .dimension = 4,
        .indices = &indices,
        .values = &values,
    };

    try vector.validate();
    try std.testing.expectEqual(@as(f64, 21.0), try vector.dotDense(&.{ 0.0, 0.0, 3.0, 0.0 }));
}

test "sparse-vector validation rejects broken invariants" {
    const RowId = Foundation.RowId;
    const id0 = RowId.fromUsizeAssumeValid(0);
    const id1 = RowId.fromUsizeAssumeValid(1);
    const id2 = RowId.fromUsizeAssumeValid(2);

    try std.testing.expectError(
        SparseVectorError.LengthMismatch,
        (SparseVectorView(RowId){
            .dimension = 3,
            .indices = &.{ id0, id1 },
            .values = &.{1.0},
        }).validate(),
    );
    try std.testing.expectError(
        SparseVectorError.IndicesNotStrictlyIncreasing,
        (SparseVectorView(RowId){
            .dimension = 3,
            .indices = &.{ id1, id0 },
            .values = &.{ 1.0, 2.0 },
        }).validate(),
    );
    try std.testing.expectError(
        SparseVectorError.IndicesNotStrictlyIncreasing,
        (SparseVectorView(RowId){
            .dimension = 3,
            .indices = &.{ id1, id1 },
            .values = &.{ 1.0, 2.0 },
        }).validate(),
    );
    try std.testing.expectError(
        SparseVectorError.IndexOutOfBounds,
        (SparseVectorView(RowId){
            .dimension = 2,
            .indices = &.{id2},
            .values = &.{1.0},
        }).validate(),
    );
    try std.testing.expectError(
        SparseVectorError.ExplicitZero,
        (SparseVectorView(RowId){
            .dimension = 3,
            .indices = &.{id0},
            .values = &.{0.0},
        }).validate(),
    );
    try std.testing.expectError(
        SparseVectorError.NonFiniteValue,
        (SparseVectorView(RowId){
            .dimension = 3,
            .indices = &.{id0},
            .values = &.{std.math.nan(f64)},
        }).validate(),
    );
}

test "checked dense operations reject dimension mismatch" {
    const RowId = Foundation.RowId;
    const vector: SparseVectorView(RowId) = .{
        .dimension = 2,
        .indices = &.{RowId.fromUsizeAssumeValid(0)},
        .values = &.{1.0},
    };

    try std.testing.expectError(
        SparseVectorError.DimensionMismatch,
        vector.dotDense(&.{1.0}),
    );
    var target = [_]f64{1.0};
    try std.testing.expectError(
        SparseVectorError.DimensionMismatch,
        vector.addToDense(1.0, &target),
    );
}
