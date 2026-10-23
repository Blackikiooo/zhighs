//! Dense Cholesky factorization for the SPD normal-equations system.
//!
//! The normal matrix A·Θ·Aᵀ is assembled densely (m is the row count), so the
//! factorization here is the classical right-looking dense variant over a
//! full n×n row-major buffer of which only the lower triangle is referenced.

const std = @import("std");
const Allocator = std.mem.Allocator;

pub const FactorError = error{NotPositiveDefinite};

/// Owning dense lower-triangular factor storage with an in-place solve.
pub const DenseCholesky = struct {
    /// Logical order of the factored system.
    n: usize,
    /// Row-major n×n buffer; the lower triangle holds L after `factor`.
    data: []f64,
    allocator: Allocator,

    pub fn init(allocator: Allocator, n: usize) !DenseCholesky {
        return .{
            .n = n,
            .data = try allocator.alloc(f64, n * n),
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *DenseCholesky) void {
        self.allocator.free(self.data);
    }

    /// Factor the SPD matrix whose lower triangle is given in row-major
    /// `matrix` (n×n). Pivots below `pivot_floor` are lifted to the floor
    /// (per-pivot regularization); pass 0 to reject non-positive pivots
    /// instead. Non-finite pivots always fail. Fails without modifying the
    /// caller's buffer.
    pub fn factor(self: *DenseCholesky, matrix: []const f64, pivot_floor: f64) FactorError!void {
        const n = self.n;
        std.debug.assert(matrix.len == n * n);
        @memcpy(self.data, matrix);
        const l = self.data;
        for (0..n) |j| {
            var pivot = l[j * n + j];
            for (0..j) |k| {
                pivot -= l[j * n + k] * l[j * n + k];
            }
            if (!std.math.isFinite(pivot) or (pivot_floor <= 0.0 and pivot <= 0.0)) return error.NotPositiveDefinite;
            const diag = @sqrt(@max(pivot, pivot_floor));
            l[j * n + j] = diag;
            for (j + 1..n) |i| {
                var acc = l[i * n + j];
                for (0..j) |k| {
                    acc -= l[i * n + k] * l[j * n + k];
                }
                l[i * n + j] = acc / diag;
            }
        }
    }

    /// Solve L·Lᵀ·x = b in place, overwriting `rhs`.
    pub fn solve(self: *const DenseCholesky, rhs: []f64) void {
        const n = self.n;
        std.debug.assert(rhs.len == n);
        const l = self.data;
        for (0..n) |i| {
            var acc = rhs[i];
            for (0..i) |k| {
                acc -= l[i * n + k] * rhs[k];
            }
            rhs[i] = acc / l[i * n + i];
        }
        var i = n;
        while (i > 0) {
            i -= 1;
            var acc = rhs[i];
            for (i + 1..n) |k| {
                acc -= l[k * n + i] * rhs[k];
            }
            rhs[i] = acc / l[i * n + i];
        }
    }
};

test "DenseCholesky factors and solves a small SPD system" {
    const allocator = std.testing.allocator;
    // [[4, 2], [2, 3]]; solve against rhs [6, 5] -> x = [8/8, ...] = [1, 1].
    const matrix = [_]f64{ 4.0, 2.0, 2.0, 3.0 };
    var chol = try DenseCholesky.init(allocator, 2);
    defer chol.deinit();
    try chol.factor(&matrix, 0.0);
    var rhs = [_]f64{ 6.0, 5.0 };
    chol.solve(&rhs);
    try std.testing.expectApproxEqAbs(1.0, rhs[0], 1e-12);
    try std.testing.expectApproxEqAbs(1.0, rhs[1], 1e-12);
}

test "DenseCholesky rejects indefinite input" {
    const allocator = std.testing.allocator;
    const matrix = [_]f64{ 1.0, 2.0, 2.0, 1.0 };
    var chol = try DenseCholesky.init(allocator, 2);
    defer chol.deinit();
    try std.testing.expectError(error.NotPositiveDefinite, chol.factor(&matrix, 0.0));
}