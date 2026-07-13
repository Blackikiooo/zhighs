//! Dense LU factorization with partial pivoting.
//!
//! This is a correctness baseline for simplex and other linear-system users.
//! The storage is row-major and owned by the factorization; sparse LU and
//! update-based factorizations can implement the same solver boundary later.

const std = @import("std");

pub const DenseLuError = error{ DimensionMismatch, Singular, NumericalFailure, OutOfMemory };

pub const DenseLU = struct {
    allocator: std.mem.Allocator,
    n: usize = 0,
    lu: []f64 = &.{},
    pivots: []usize = &.{},
    work: []f64 = &.{},
    pivot_tolerance: f64 = 1e-12,
    min_abs_pivot: f64 = std.math.inf(f64),
    max_abs_pivot: f64 = 0.0,

    pub fn init(allocator: std.mem.Allocator) DenseLU {
        return .{ .allocator = allocator };
    }

    pub fn deinit(self: *DenseLU) void {
        self.allocator.free(self.lu);
        self.allocator.free(self.pivots);
        self.allocator.free(self.work);
        self.* = .{ .allocator = self.allocator };
    }

    /// Factorize an n-by-n row-major matrix without modifying the input.
    pub fn factorize(self: *DenseLU, n: usize, matrix: []const f64) DenseLuError!void {
        if (matrix.len != n * n) return error.DimensionMismatch;
        const next_lu = self.allocator.alloc(f64, matrix.len) catch return error.OutOfMemory;
        errdefer self.allocator.free(next_lu);
        const next_pivots = self.allocator.alloc(usize, n) catch return error.OutOfMemory;
        errdefer self.allocator.free(next_pivots);
        const next_work = self.allocator.alloc(f64, n) catch return error.OutOfMemory;
        errdefer self.allocator.free(next_work);
        @memcpy(next_lu, matrix);
        for (next_pivots, 0..) |*pivot, i| pivot.* = i;

        var min_abs_pivot = std.math.inf(f64);
        var max_abs_pivot: f64 = 0.0;
        for (0..n) |column| {
            var pivot_row = column;
            var pivot_abs: f64 = 0.0;
            for (column..n) |row| {
                const value = @abs(next_lu[row * n + column]);
                if (value > pivot_abs) {
                    pivot_abs = value;
                    pivot_row = row;
                }
            }
            if (!std.math.isFinite(pivot_abs) or pivot_abs <= self.pivot_tolerance) return error.Singular;
            if (pivot_row != column) {
                for (0..n) |j| std.mem.swap(f64, &next_lu[column * n + j], &next_lu[pivot_row * n + j]);
                std.mem.swap(usize, &next_pivots[column], &next_pivots[pivot_row]);
            }
            const diagonal = next_lu[column * n + column];
            min_abs_pivot = @min(min_abs_pivot, @abs(diagonal));
            max_abs_pivot = @max(max_abs_pivot, @abs(diagonal));
            for (column + 1..n) |row| {
                next_lu[row * n + column] /= diagonal;
                const multiplier = next_lu[row * n + column];
                for (column + 1..n) |j| next_lu[row * n + j] -= multiplier * next_lu[column * n + j];
            }
        }
        self.allocator.free(self.lu);
        self.allocator.free(self.pivots);
        self.allocator.free(self.work);
        self.lu = next_lu;
        self.pivots = next_pivots;
        self.work = next_work;
        self.n = n;
        self.min_abs_pivot = if (n == 0) 1.0 else min_abs_pivot;
        self.max_abs_pivot = if (n == 0) 1.0 else max_abs_pivot;
    }

    /// Zero-copy variant. Ownership of `matrix` is transferred to this LU
    /// object; the caller must not free or reuse it after success. The buffer
    /// is overwritten with LU data and released on failure.
    pub fn factorizeOwned(self: *DenseLU, n: usize, matrix: []f64) DenseLuError!void {
        if (matrix.len != n * n) return error.DimensionMismatch;
        errdefer self.allocator.free(matrix);
        const next_pivots = self.allocator.alloc(usize, n) catch return error.OutOfMemory;
        errdefer self.allocator.free(next_pivots);
        const next_work = self.allocator.alloc(f64, n) catch return error.OutOfMemory;
        errdefer self.allocator.free(next_work);
        for (next_pivots, 0..) |*pivot, i| pivot.* = i;
        var min_abs_pivot = std.math.inf(f64);
        var max_abs_pivot: f64 = 0.0;
        for (0..n) |column| {
            var pivot_row = column;
            var pivot_abs: f64 = 0.0;
            for (column..n) |row| {
                const value = @abs(matrix[row * n + column]);
                if (value > pivot_abs) {
                    pivot_abs = value;
                    pivot_row = row;
                }
            }
            if (!std.math.isFinite(pivot_abs) or pivot_abs <= self.pivot_tolerance) return error.Singular;
            if (pivot_row != column) {
                for (0..n) |j| std.mem.swap(f64, &matrix[column * n + j], &matrix[pivot_row * n + j]);
                std.mem.swap(usize, &next_pivots[column], &next_pivots[pivot_row]);
            }
            const diagonal = matrix[column * n + column];
            min_abs_pivot = @min(min_abs_pivot, @abs(diagonal));
            max_abs_pivot = @max(max_abs_pivot, @abs(diagonal));
            for (column + 1..n) |row| {
                matrix[row * n + column] /= diagonal;
                const multiplier = matrix[row * n + column];
                for (column + 1..n) |j| matrix[row * n + j] -= multiplier * matrix[column * n + j];
            }
        }
        self.allocator.free(self.lu);
        self.allocator.free(self.pivots);
        self.allocator.free(self.work);
        self.lu = matrix;
        self.pivots = next_pivots;
        self.work = next_work;
        self.n = n;
        self.min_abs_pivot = if (n == 0) 1.0 else min_abs_pivot;
        self.max_abs_pivot = if (n == 0) 1.0 else max_abs_pivot;
    }

    /// Re-factorize the existing `lu` buffer after the caller overwrites it
    /// with a new row-major matrix. No allocation is performed.
    pub fn refactorizeInPlace(self: *DenseLU) DenseLuError!void {
        const n = self.n;
        if (self.lu.len != n * n or self.pivots.len != n or self.work.len != n) return error.DimensionMismatch;
        for (self.pivots, 0..) |*pivot, i| pivot.* = i;
        var min_abs_pivot = std.math.inf(f64);
        var max_abs_pivot: f64 = 0.0;
        for (0..n) |column| {
            var pivot_row = column;
            var pivot_abs: f64 = 0.0;
            for (column..n) |row| {
                const value = @abs(self.lu[row * n + column]);
                if (value > pivot_abs) {
                    pivot_abs = value;
                    pivot_row = row;
                }
            }
            if (!std.math.isFinite(pivot_abs) or pivot_abs <= self.pivot_tolerance) return error.Singular;
            if (pivot_row != column) {
                for (0..n) |j| std.mem.swap(f64, &self.lu[column * n + j], &self.lu[pivot_row * n + j]);
                std.mem.swap(usize, &self.pivots[column], &self.pivots[pivot_row]);
            }
            const diagonal = self.lu[column * n + column];
            min_abs_pivot = @min(min_abs_pivot, @abs(diagonal));
            max_abs_pivot = @max(max_abs_pivot, @abs(diagonal));
            for (column + 1..n) |row| {
                self.lu[row * n + column] /= diagonal;
                const multiplier = self.lu[row * n + column];
                for (column + 1..n) |j| self.lu[row * n + j] -= multiplier * self.lu[column * n + j];
            }
        }
        self.min_abs_pivot = if (n == 0) 1.0 else min_abs_pivot;
        self.max_abs_pivot = if (n == 0) 1.0 else max_abs_pivot;
    }

    /// Cheap stability signal used for refactorization policy. This is pivot
    /// spread, deliberately not presented as a rigorous matrix condition number.
    pub fn pivotConditionEstimate(self: DenseLU) f64 {
        if (self.min_abs_pivot <= 0.0) return std.math.inf(f64);
        return self.max_abs_pivot / self.min_abs_pivot;
    }

    pub fn solve(self: *DenseLU, rhs: []f64) DenseLuError!void {
        if (rhs.len != self.n) return error.DimensionMismatch;
        for (0..self.n) |i| {
            self.work[i] = rhs[self.pivots[i]];
            for (0..i) |j| self.work[i] -= self.lu[i * self.n + j] * self.work[j];
        }
        var row = self.n;
        while (row > 0) {
            row -= 1;
            for (row + 1..self.n) |j| self.work[row] -= self.lu[row * self.n + j] * self.work[j];
            self.work[row] /= self.lu[row * self.n + row];
        }
        @memcpy(rhs, self.work);
    }

    pub fn solveTranspose(self: *DenseLU, rhs: []f64) DenseLuError!void {
        if (rhs.len != self.n) return error.DimensionMismatch;
        @memcpy(self.work, rhs);

        // Solve U^T z = rhs by scattering each completed component through
        // one contiguous row of U.  The equivalent gather formulation reads
        // lu[j * n + i] down a column of row-major storage and becomes
        // cache/TLB bound once the factorization no longer fits in cache.
        for (0..self.n) |i| {
            const row = self.lu[i * self.n ..][0..self.n];
            const solved = self.work[i] / row[i];
            self.work[i] = solved;
            for (self.work[i + 1 ..], row[i + 1 ..]) |*pending, coefficient|
                pending.* -= coefficient * solved;
        }

        // Solve L^T y = z in reverse order.  L has an implicit unit diagonal;
        // scattering y[row] into earlier entries again traverses one
        // contiguous row instead of gathering a strided column.
        var i = self.n;
        while (i > 0) {
            i -= 1;
            const solved = self.work[i];
            const row = self.lu[i * self.n ..][0..i];
            for (self.work[0..i], row) |*pending, coefficient|
                pending.* -= coefficient * solved;
        }
        for (self.pivots, 0..) |original, pivot_index| rhs[original] = self.work[pivot_index];
    }
};

test "DenseLU solves and transpose-solves" {
    var lu = DenseLU.init(std.testing.allocator);
    defer lu.deinit();
    try lu.factorize(2, &[_]f64{ 4, 3, 6, 3 });
    var rhs = [_]f64{ 10, 12 };
    try lu.solve(&rhs);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), rhs[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), rhs[1], 1e-12);
    var trhs = [_]f64{ 10, 12 };
    try lu.solveTranspose(&trhs);
    try std.testing.expectApproxEqAbs(@as(f64, 7.0), trhs[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, -3.0), trhs[1], 1e-12);
}

test "DenseLU owned factorization avoids matrix copy" {
    var lu = DenseLU.init(std.testing.allocator);
    defer lu.deinit();
    const owned = try std.testing.allocator.alloc(f64, 4);
    @memcpy(owned, &[_]f64{ 2, 0, 0, 3 });
    try lu.factorizeOwned(2, owned);
    try std.testing.expectEqual(@as(usize, 2), lu.n);
    var rhs = [_]f64{ 4, 6 };
    try lu.solve(&rhs);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), rhs[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), rhs[1], 1e-12);
}

test "DenseLU transpose solve handles a larger dense factorization" {
    const dimension = 17;
    var matrix: [dimension * dimension]f64 = undefined;
    for (0..dimension) |row| {
        for (0..dimension) |col| {
            const pattern: i32 = @intCast((row * 17 + col * 13) % 11);
            matrix[row * dimension + col] = if (row == col)
                20.0 + @as(f64, @floatFromInt(row))
            else
                @as(f64, @floatFromInt(pattern - 5)) * 0.1;
        }
    }

    var lu = DenseLU.init(std.testing.allocator);
    defer lu.deinit();
    try lu.factorize(dimension, &matrix);

    var rhs: [dimension]f64 = undefined;
    for (&rhs, 0..) |*value, i| value.* = @as(f64, @floatFromInt(i + 1)) * 0.25 - 1.0;
    const original_rhs = rhs;
    try lu.solveTranspose(&rhs);

    for (0..dimension) |col| {
        var product: f64 = 0.0;
        for (0..dimension) |row| product += matrix[row * dimension + col] * rhs[row];
        try std.testing.expectApproxEqAbs(original_rhs[col], product, 1e-10);
    }
}

/// The original column-gather BTRAN formulation is intentionally retained as
/// a test oracle. The production implementation traverses rows contiguously,
/// so comparing both paths catches indexing and final-pivot-permutation errors.
fn solveTransposeGatherReference(lu: *const DenseLU, rhs: []f64, work: []f64) void {
    const n = lu.n;
    for (0..n) |i| {
        var value = rhs[i];
        for (0..i) |j| value -= lu.lu[j * n + i] * work[j];
        work[i] = value / lu.lu[i * n + i];
    }
    var i = n;
    while (i > 0) {
        i -= 1;
        var value = work[i];
        for (i + 1..n) |j| value -= lu.lu[j * n + i] * work[j];
        work[i] = value;
    }
    for (lu.pivots, 0..) |original, pivot_index| rhs[original] = work[pivot_index];
}

fn expectRelative(expected: f64, actual: f64, scale: f64) !void {
    const tolerance = scale * @max(1.0, @abs(expected));
    try std.testing.expectApproxEqAbs(expected, actual, tolerance);
}

test "DenseLU randomized FTRAN and BTRAN differential validates pivot permutation" {
    var prng = std.Random.DefaultPrng.init(0x6c75_7069_766f_7473);
    const random = prng.random();
    const dimensions = [_]usize{ 2, 3, 5, 8, 16, 31 };

    for (dimensions) |n| {
        for (0..4) |_| {
            const base = try std.testing.allocator.alloc(f64, n * n);
            defer std.testing.allocator.free(base);
            const matrix = try std.testing.allocator.alloc(f64, n * n);
            defer std.testing.allocator.free(matrix);
            const expected = try std.testing.allocator.alloc(f64, n);
            defer std.testing.allocator.free(expected);
            const ftran_rhs = try std.testing.allocator.alloc(f64, n);
            defer std.testing.allocator.free(ftran_rhs);
            const btran_rhs = try std.testing.allocator.alloc(f64, n);
            defer std.testing.allocator.free(btran_rhs);
            const btran_reference = try std.testing.allocator.alloc(f64, n);
            defer std.testing.allocator.free(btran_reference);
            const reference_work = try std.testing.allocator.alloc(f64, n);
            defer std.testing.allocator.free(reference_work);

            // Build a well-conditioned matrix, then cyclically permute its
            // rows. The largest first-column entry is no longer in row zero,
            // which guarantees that the test exercises partial pivoting.
            for (0..n) |row| {
                var off_diagonal_sum: f64 = 0.0;
                for (0..n) |col| {
                    if (row == col) continue;
                    const value = random.float(f64) * 0.5 - 0.25;
                    base[row * n + col] = value;
                    off_diagonal_sum += @abs(value);
                }
                base[row * n + row] = off_diagonal_sum + 1.0;
            }
            for (0..n) |physical_row| {
                const source_row = (physical_row + 1) % n;
                @memcpy(matrix[physical_row * n ..][0..n], base[source_row * n ..][0..n]);
            }
            for (expected) |*value| value.* = random.float(f64) * 2.0 - 1.0;

            for (0..n) |row| {
                var product: f64 = 0.0;
                for (0..n) |col| product += matrix[row * n + col] * expected[col];
                ftran_rhs[row] = product;
            }
            for (0..n) |col| {
                var product: f64 = 0.0;
                for (0..n) |row| product += matrix[row * n + col] * expected[row];
                btran_rhs[col] = product;
            }
            @memcpy(btran_reference, btran_rhs);

            var lu = DenseLU.init(std.testing.allocator);
            defer lu.deinit();
            try lu.factorize(n, matrix);
            var has_row_pivot = false;
            for (lu.pivots, 0..) |original, pivot_index|
                has_row_pivot = has_row_pivot or original != pivot_index;
            try std.testing.expect(has_row_pivot);

            try lu.solve(ftran_rhs);
            try lu.solveTranspose(btran_rhs);
            solveTransposeGatherReference(&lu, btran_reference, reference_work);

            const tolerance = 5e-12 * @as(f64, @floatFromInt(n));
            for (expected, ftran_rhs, btran_rhs, btran_reference) |want, ftran, btran, reference| {
                try expectRelative(want, ftran, tolerance);
                try expectRelative(want, btran, tolerance);
                try expectRelative(reference, btran, tolerance);
            }
        }
    }
}

test "DenseLU pivot tolerance controls near-singular acceptance" {
    const tiny = [_]f64{1e-13};

    var default_lu = DenseLU.init(std.testing.allocator);
    defer default_lu.deinit();
    try std.testing.expectError(error.Singular, default_lu.factorize(1, &tiny));

    var permissive_lu = DenseLU.init(std.testing.allocator);
    defer permissive_lu.deinit();
    permissive_lu.pivot_tolerance = 1e-15;
    try permissive_lu.factorize(1, &tiny);
    var rhs = [_]f64{2e-13};
    try permissive_lu.solve(&rhs);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), rhs[0], 1e-12);
}
