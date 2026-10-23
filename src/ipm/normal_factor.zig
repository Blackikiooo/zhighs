const std = @import("std");
const Dense = @import("cholesky.zig").DenseCholesky;
const Sparse = @import("sparse_cholesky.zig").SparseCholesky;
const StandardForm = @import("standard_form.zig").StandardForm;
const graph_mod = @import("graph.zig");

pub const NormalFactor = union(enum) {
    dense: Dense,
    sparse: Sparse,

    /// The sparse path consumes the shared normal-equations graph and the
    /// selected elimination order; both are ignored (and may be null) for
    /// the dense reference kernel.
    pub fn init(a: std.mem.Allocator, sf: *const StandardForm, dense: bool, g: ?*const graph_mod.Graph, order: ?*const graph_mod.Order) !NormalFactor {
        return if (dense)
            .{ .dense = try Dense.init(a, sf.num_rows) }
        else
            .{ .sparse = try Sparse.init(a, sf, g.?, order.?) };
    }
    pub fn deinit(self: *NormalFactor) void {
        switch (self.*) {
            inline else => |*value| value.deinit(),
        }
    }
    pub fn storageLen(self: *const NormalFactor) usize {
        return switch (self.*) {
            inline else => |value| value.data.len,
        };
    }
    pub fn productCount(self: *const NormalFactor) usize {
        return switch (self.*) {
            .dense => 0,
            .sparse => |value| value.productCount(),
        };
    }
    pub fn indirectProductCount(self: *const NormalFactor) usize {
        return switch (self.*) {
            .dense => 0,
            .sparse => |value| value.indirect_count,
        };
    }
    pub fn index(self: *const NormalFactor, row: usize, col: usize) usize {
        return switch (self.*) {
            .dense => |value| row * value.n + col,
            .sparse => |*value| value.index(row, col),
        };
    }
    /// Factor the matrix, lifting degenerate pivots to `pivot_floor`
    /// (per-pivot regularization); pass 0 to reject them instead.
    pub fn factor(self: *NormalFactor, values: []const f64, pivot_floor: f64) error{NotPositiveDefinite}!void {
        switch (self.*) {
            inline else => |*value| try value.factor(values, pivot_floor),
        }
    }
    pub fn solve(self: *const NormalFactor, rhs: []f64) void {
        switch (self.*) {
            inline else => |*value| value.solve(rhs),
        }
    }

    /// out := M*x for a symmetric matrix stored by its lower triangle.
    pub fn multiply(self: *const NormalFactor, matrix: []const f64, x: []const f64, out: []f64) void {
        @memset(out, 0.0);
        switch (self.*) {
            .dense => |value| {
                for (0..value.n) |i| {
                    for (0..i + 1) |j| {
                        const entry = matrix[i * value.n + j];
                        out[i] += entry * x[j];
                        if (i != j) out[j] += entry * x[i];
                    }
                }
            },
            .sparse => |value| {
                for (0..value.row_starts.len - 1) |i| {
                    for (value.row_starts[i]..value.row_starts[i + 1]) |position| {
                        const j = value.col_indices[position];
                        const entry = matrix[position];
                        out[i] += entry * x[j];
                        if (i != j) out[j] += entry * x[i];
                    }
                }
            },
        }
    }
};