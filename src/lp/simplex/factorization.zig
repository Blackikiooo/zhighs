//! Basis factorization policy boundary.
//!
//! The first implementation may use dense LU for correctness tests; the
//! interface deliberately permits sparse LU, eta/Forrest--Tomlin updates, and
//! periodic reinversion without changing the simplex engine.

const std = @import("std");
const matrix = @import("matrix");

pub const FactorizationError = error{ DimensionMismatch, NotImplemented, Singular, NumericalFailure, OutOfMemory };
/// Selected numeric kernel. Sparse LU is intentionally a separate future
/// backend rather than being embedded in simplex orchestration.
pub const BackendKind = enum { dense_lu };
pub const FactorizationStats = struct {
    factorizations: usize = 0,
    ftran_calls: usize = 0,
    btran_calls: usize = 0,
    eta_updates: usize = 0,
};
pub const PivotUpdateView = struct {
    leaving_row: u32,
    entering_col: u32,
    direction: []const f64,
    /// Converts a signed movement direction back to the actual entering
    /// basis column. Primal simplex uses -1 when entering from an upper bound.
    column_scale: f64 = 1.0,
};

pub const Factorization = struct {
    allocator: std.mem.Allocator,
    update_count: usize = 0,
    dense_lu: matrix.DenseLU,
    eta_values: []f64 = &.{},
    eta_rows: []u32 = &.{},
    eta_count: usize = 0,
    eta_capacity: usize = 64,
    dimension: usize = 0,
    backend_kind: BackendKind = .dense_lu,
    stats: FactorizationStats = .{},

    pub fn init(allocator: std.mem.Allocator) Factorization {
        return .{ .allocator = allocator, .dense_lu = matrix.DenseLU.init(allocator) };
    }
    pub fn deinit(self: *Factorization) void {
        self.dense_lu.deinit();
        self.allocator.free(self.eta_values);
        self.allocator.free(self.eta_rows);
    }
    pub fn factorize(self: *Factorization, n: usize, matrix_data: []const f64) FactorizationError!void {
        self.dense_lu.factorize(n, matrix_data) catch |err| return switch (err) {
            error.DimensionMismatch => error.DimensionMismatch,
            error.Singular => error.Singular,
            error.NumericalFailure => error.NumericalFailure,
            error.OutOfMemory => error.OutOfMemory,
        };
        try self.prepareEtaStorage(n);
        self.stats.factorizations += 1;
    }

    /// Build and factorize an identity basis without an intermediate matrix
    /// copy. Used by the canonical slack-basis crash initializer.
    pub fn factorizeIdentity(self: *Factorization, n: usize) FactorizationError!void {
        const data = self.allocator.alloc(f64, n * n) catch return error.OutOfMemory;
        @memset(data, 0.0);
        for (0..n) |i| data[i * n + i] = 1.0;
        self.dense_lu.factorizeOwned(n, data) catch |err| return switch (err) {
            error.DimensionMismatch => error.DimensionMismatch,
            error.Singular => error.Singular,
            error.NumericalFailure => error.NumericalFailure,
            error.OutOfMemory => error.OutOfMemory,
        };
        try self.prepareEtaStorage(n);
        self.stats.factorizations += 1;
    }

    pub fn mutableBasisBuffer(self: *Factorization, n: usize) FactorizationError![]f64 {
        if (self.dense_lu.n != n or self.dense_lu.lu.len != n * n) return error.DimensionMismatch;
        return self.dense_lu.lu;
    }

    pub fn refactorizeInPlace(self: *Factorization) FactorizationError!void {
        self.dense_lu.refactorizeInPlace() catch |err| return switch (err) {
            error.DimensionMismatch => error.DimensionMismatch,
            error.Singular => error.Singular,
            error.NumericalFailure => error.NumericalFailure,
            error.OutOfMemory => error.OutOfMemory,
        };
        self.update_count = 0;
        self.eta_count = 0;
        self.stats.factorizations += 1;
    }
    pub fn solve(self: *Factorization, rhs: []f64) FactorizationError!void {
        self.stats.ftran_calls += 1;
        self.dense_lu.solve(rhs) catch |err| return switch (err) {
            error.DimensionMismatch => error.DimensionMismatch,
            error.Singular => error.Singular,
            error.NumericalFailure => error.NumericalFailure,
            error.OutOfMemory => error.OutOfMemory,
        };
        for (0..self.eta_count) |update_index| self.applyEtaInverse(update_index, rhs);
    }
    pub fn solveTranspose(self: *Factorization, rhs: []f64) FactorizationError!void {
        self.stats.btran_calls += 1;
        var update_index = self.eta_count;
        while (update_index > 0) {
            update_index -= 1;
            self.applyEtaInverseTranspose(update_index, rhs);
        }
        self.dense_lu.solveTranspose(rhs) catch |err| return switch (err) {
            error.DimensionMismatch => error.DimensionMismatch,
            error.Singular => error.Singular,
            error.NumericalFailure => error.NumericalFailure,
            error.OutOfMemory => error.OutOfMemory,
        };
    }
    pub fn update(self: *Factorization, update_view: PivotUpdateView) FactorizationError!void {
        if (update_view.direction.len != self.dimension or update_view.leaving_row >= self.dimension) return error.DimensionMismatch;
        if (self.eta_count >= self.eta_capacity) return error.NotImplemented;
        const pivot_row: usize = @intCast(update_view.leaving_row);
        if (!std.math.isFinite(update_view.column_scale) or
            @abs(update_view.direction[pivot_row] * update_view.column_scale) <= self.dense_lu.pivot_tolerance)
            return error.Singular;
        const eta = self.eta_values[self.eta_count * self.dimension ..][0..self.dimension];
        for (eta, update_view.direction) |*value, direction| value.* = direction * update_view.column_scale;
        self.eta_rows[self.eta_count] = update_view.leaving_row;
        self.eta_count += 1;
        self.update_count += 1;
        self.stats.eta_updates += 1;
    }
    pub fn needsRefactor(self: Factorization, limit: usize) bool {
        return self.update_count >= limit;
    }

    fn prepareEtaStorage(self: *Factorization, dimension: usize) FactorizationError!void {
        self.eta_values = self.allocator.realloc(self.eta_values, dimension * self.eta_capacity) catch return error.OutOfMemory;
        self.eta_rows = self.allocator.realloc(self.eta_rows, self.eta_capacity) catch return error.OutOfMemory;
        self.dimension = dimension;
        self.eta_count = 0;
        self.update_count = 0;
    }

    fn applyEtaInverse(self: *const Factorization, update_index: usize, vector: []f64) void {
        const row: usize = @intCast(self.eta_rows[update_index]);
        const eta = self.eta_values[update_index * self.dimension ..][0..self.dimension];
        const pivot_value = vector[row] / eta[row];
        for (vector, eta, 0..) |*value, coefficient, i| {
            if (i != row) value.* -= coefficient * pivot_value;
        }
        vector[row] = pivot_value;
    }

    fn applyEtaInverseTranspose(self: *const Factorization, update_index: usize, vector: []f64) void {
        const row: usize = @intCast(self.eta_rows[update_index]);
        const eta = self.eta_values[update_index * self.dimension ..][0..self.dimension];
        var value = vector[row];
        for (vector, eta, 0..) |entry, coefficient, i| {
            if (i != row) value -= coefficient * entry;
        }
        vector[row] = value / eta[row];
    }
};

test {
    std.testing.refAllDecls(@This());
}

test "eta update applies FTRAN and BTRAN without refactorization" {
    var factorization = Factorization.init(std.testing.allocator);
    defer factorization.deinit();
    try factorization.factorizeIdentity(2);
    try factorization.update(.{ .leaving_row = 0, .entering_col = 0, .direction = &[_]f64{ 2, 1 } });
    var rhs = [_]f64{ 4, 3 };
    try factorization.solve(&rhs);
    try std.testing.expectApproxEqAbs(@as(f64, 2), rhs[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1), rhs[1], 1e-12);
    var transpose_rhs = [_]f64{ 5, 1 };
    try factorization.solveTranspose(&transpose_rhs);
    try std.testing.expectApproxEqAbs(@as(f64, 2), transpose_rhs[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1), transpose_rhs[1], 1e-12);
}

test "eta update restores the actual column from a signed movement direction" {
    var factorization = Factorization.init(std.testing.allocator);
    defer factorization.deinit();
    try factorization.factorizeIdentity(2);
    try factorization.update(.{
        .leaving_row = 0,
        .entering_col = 0,
        .direction = &[_]f64{ 2, 1 },
        .column_scale = -1,
    });
    var rhs = [_]f64{ -4, 3 };
    try factorization.solve(&rhs);
    try std.testing.expectApproxEqAbs(@as(f64, 2), rhs[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 5), rhs[1], 1e-12);
}
