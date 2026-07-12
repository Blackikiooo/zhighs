//! Hessian (quadratic objective matrix) storage for the solver-internal layer.
//!
//! Represents the quadratic term `½ xᵀ Q x` in the objective function.
//! Only one triangle of the symmetric matrix `Q` is stored — by convention the
//! lower triangle in canonical CSC order.
//!
//! ## Storage conventions
//!
//! - `format == .triangular` : CSC over the lower triangle.  Column `j` holds
//!   entries for rows `i ≥ j`.  Row indices are strictly increasing per column.
//!   Values are the raw matrix entries (not doubled).
//!
//! - `format == .diagonal` : Column `j` holds exactly one entry at row `j`.
//!   Has exactly `dimension` non‑zero values.  Ensures O(dimension) storage
//!   for separable quadratic objectives.
//!
//! ## Mathematical semantics
//!
//! ```text
//! objective += ½ · Σᵢ Σⱼ Qᵢⱼ · xᵢ · xⱼ
//! ```
//!
//! Experimental API.  The triangular CSC layout is stable; the format enum
//! may gain a `.full` variant for non‑convex heuristics.

const std = @import("std");
const foundation = @import("foundation");

const ColId = foundation.ColId;

// ── HessianFormat ──────────────────────────────────────────────────────────

/// Storage format for the triangular Hessian matrix.
pub const HessianFormat = enum(u8) {
    /// CSC over the lower triangle (including diagonal).
    triangular,
    /// Only the diagonal values are stored (dimension entries).
    diagonal,
};

// ── Curvature ──────────────────────────────────────────────────────────────

/// Declared or detected curvature of a quadratic function.
///
/// Defaults to `unknown`.  Solvers may update this field after an eigenvalue
/// analysis; it is **not** set automatically by model construction.
pub const Curvature = enum(u8) {
    /// Not yet analysed.
    unknown,
    /// Positive semi‑definite.
    convex,
    /// Negative semi‑definite.
    concave,
    /// Neither convex nor concave.
    indefinite,
};

// ── Hessian ────────────────────────────────────────────────────────────────

/// Experimental API: owning triangular Hessian matrix.
pub const Hessian = struct {
    allocator: std.mem.Allocator,
    /// Square matrix dimension.
    dimension: usize,
    /// CSC column starts (length = dimension + 1).
    starts: []usize,
    /// Sorted, unique row indices (lower‑triangle convention).
    indices: []ColId,
    /// Non‑zero values.
    values: []f64,
    /// Storage layout.
    format: HessianFormat,

    const Self = @This();

    /// Creates an empty (zero‑dimension, zero‑entry) Hessian.
    pub fn initEmpty(allocator: std.mem.Allocator) !Self {
        const dim: usize = 0;
        const starts = try allocator.alloc(usize, 1);
        errdefer allocator.free(starts);
        starts[0] = 0;
        const indices = try allocator.alloc(ColId, 0);
        errdefer allocator.free(indices);
        const values = try allocator.alloc(f64, 0);
        errdefer allocator.free(values);
        return Self{
            .allocator = allocator,
            .dimension = dim,
            .starts = starts,
            .indices = indices,
            .values = values,
            .format = .triangular,
        };
    }

    /// Releases all owned arrays.
    pub fn deinit(self: *Self) void {
        const allocator = self.allocator;
        if (self.starts.len > 0) allocator.free(self.starts);
        if (self.indices.len > 0) allocator.free(self.indices);
        if (self.values.len > 0) allocator.free(self.values);
        self.* = undefined;
    }

    /// Number of stored non‑zero entries.
    pub fn nnz(self: Self) usize {
        return self.values.len;
    }

    /// Structurally validates the Hessian.
    pub fn validate(self: *const Self) !void {
        if (self.starts.len != self.dimension + 1)
            return error.DimensionMismatch;
        if (self.indices.len != self.values.len)
            return error.DimensionMismatch;

        if (self.dimension == 0) {
            if (self.nnz() != 0) return error.DimensionMismatch;
            return;
        }

        if (self.starts[0] != 0) return error.InvalidHessian;
        if (self.starts[self.dimension] != self.nnz()) return error.InvalidHessian;

        for (0..self.dimension) |col| {
            const begin = self.starts[col];
            const end = self.starts[col + 1];
            if (begin > end or end > self.nnz()) return error.InvalidHessian;

            var previous_row: ?usize = null;
            for (begin..end) |pos| {
                const row = self.indices[pos].toUsize();
                const value = self.values[pos];
                if (!std.math.isFinite(value)) return error.NonFiniteValue;

                if (self.format == .triangular) {
                    if (row < col) return error.InvalidHessian;
                } else if (self.format == .diagonal) {
                    if (row != col) return error.InvalidHessian;
                    if (end - begin != 1) return error.InvalidHessian;
                }

                if (previous_row) |prev| {
                    if (row <= prev) return error.IndicesNotIncreasing;
                }
                previous_row = row;
            }
        }
    }
};

// ── Error set ──────────────────────────────────────────────────────────────

pub const HessianError = error{
    DimensionMismatch,
    InvalidHessian,
    NonFiniteValue,
    IndicesNotIncreasing,
};

// ── Tests ──────────────────────────────────────────────────────────────────

test "Hessian.initEmpty / deinit does not leak" {
    var h = try Hessian.initEmpty(std.testing.allocator);
    defer h.deinit();
    try std.testing.expectEqual(@as(usize, 0), h.dimension);
    try std.testing.expectEqual(@as(usize, 0), h.nnz());
}

test "Hessian.empty passes validation" {
    var h = try Hessian.initEmpty(std.testing.allocator);
    defer h.deinit();
    try h.validate();
}

test "Hessian.diagonal storage" {
    const allocator = std.testing.allocator;
    const starts = try allocator.dupe(usize, &.{ 0, 1, 2, 3 });
    const indices = try allocator.alloc(ColId, 3);
    const values = try allocator.alloc(f64, 3);
    for ([_]usize{ 0, 1, 2 }, 0..) |row, i| {
        indices[i] = ColId.fromUsizeAssumeValid(row);
        values[i] = @floatFromInt(i + 1);
    }
    var h = Hessian{
        .allocator = allocator,
        .dimension = 3,
        .starts = starts,
        .indices = indices,
        .values = values,
        .format = .diagonal,
    };
    defer h.deinit();
    try h.validate();
    try std.testing.expectEqual(@as(usize, 3), h.nnz());
}

test "Hessian.triangular validates lower-triangle constraint" {
    // Manually build a triangular Hessian where a column has row < col (invalid).
    const allocator = std.testing.allocator;
    const starts = try allocator.dupe(usize, &.{ 0, 1, 2 });
    const indices = try allocator.alloc(ColId, 2);
    indices[0] = ColId.fromUsizeAssumeValid(1); // col 0, row 1 (valid: 1 >= 0)
    indices[1] = ColId.fromUsizeAssumeValid(0); // col 1, row 0 (invalid: 0 < 1)
    const values = try allocator.alloc(f64, 2);
    values[0] = 2.0;
    values[1] = 3.0;
    var h = Hessian{
        .allocator = allocator,
        .dimension = 2,
        .starts = starts,
        .indices = indices,
        .values = values,
        .format = .triangular,
    };
    defer h.deinit();
    try std.testing.expectError(error.InvalidHessian, h.validate());
}

test "Hessian.rejects non-finite value" {
    const allocator = std.testing.allocator;
    const starts = try allocator.dupe(usize, &.{ 0, 1, 2 });
    const indices = try allocator.alloc(ColId, 2);
    indices[0] = ColId.fromUsizeAssumeValid(0);
    indices[1] = ColId.fromUsizeAssumeValid(1);
    const values = try allocator.alloc(f64, 2);
    values[0] = 1.0;
    values[1] = std.math.nan(f64);
    var h = Hessian{
        .allocator = allocator,
        .dimension = 2,
        .starts = starts,
        .indices = indices,
        .values = values,
        .format = .triangular,
    };
    defer h.deinit();
    try std.testing.expectError(error.NonFiniteValue, h.validate());
}

test "Hessian.rejects duplicate row indices" {
    const allocator = std.testing.allocator;
    const starts = try allocator.dupe(usize, &.{ 0, 2, 2 });
    const indices = try allocator.alloc(ColId, 2);
    indices[0] = ColId.fromUsizeAssumeValid(0);
    indices[1] = ColId.fromUsizeAssumeValid(0); // duplicate in same column
    const values = try allocator.alloc(f64, 2);
    values[0] = 1.0;
    values[1] = 2.0;
    var h = Hessian{
        .allocator = allocator,
        .dimension = 2,
        .starts = starts,
        .indices = indices,
        .values = values,
        .format = .triangular,
    };
    defer h.deinit();
    try std.testing.expectError(error.IndicesNotIncreasing, h.validate());
}

test "Curvature defaults to unknown" {
    const c: Curvature = .unknown;
    try std.testing.expectEqual(Curvature.unknown, c);
}
