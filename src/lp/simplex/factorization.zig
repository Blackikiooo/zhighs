//! Basis factorization policy boundary.
//!
//! The first implementation may use dense LU for correctness tests; the
//! interface deliberately permits sparse LU, eta/Forrest--Tomlin updates, and
//! periodic reinversion without changing the simplex engine.

const std = @import("std");
const matrix = @import("matrix");

pub const FactorizationError = error{ DimensionMismatch, NotImplemented, Singular, NumericalFailure, OutOfMemory };
/// Selected base factorization. Updates are applied above either backend and
/// are cleared whenever a new base factorization is installed.
pub const BackendKind = enum { dense_lu, sparse_lu };
pub const ReinversionReason = enum { update_limit, update_growth, solve_residual };
pub const FactorizationStats = struct {
    factorizations: usize = 0,
    ftran_calls: usize = 0,
    btran_calls: usize = 0,
    eta_updates: usize = 0,
    ft_updates: usize = 0,
    maximum_update_growth: f64 = 1.0,
    update_limit_reinversions: usize = 0,
    update_growth_reinversions: usize = 0,
    solve_residual_reinversions: usize = 0,
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
    sparse_lu: matrix.SparseLU,
    sparse_basis: matrix.SparseBasisBuffers,
    eta_values: []f64 = &.{},
    eta_rows: []u32 = &.{},
    identity_basic: []u32 = &.{},
    identity_scale: []f64 = &.{},
    identity_sign: []f64 = &.{},
    eta_count: usize = 0,
    eta_capacity: usize = 64,
    dimension: usize = 0,
    backend_kind: BackendKind = .dense_lu,
    stats: FactorizationStats = .{},
    /// Largest max(|d|)/|d[p]| observed since reinversion. This inexpensive
    /// signal estimates how strongly an update can amplify solve error.
    maximum_update_growth: f64 = 1.0,

    pub fn init(allocator: std.mem.Allocator) Factorization {
        return .{
            .allocator = allocator,
            .dense_lu = matrix.DenseLU.init(allocator),
            .sparse_lu = matrix.SparseLU.init(allocator),
            .sparse_basis = matrix.SparseBasisBuffers.init(allocator),
        };
    }
    pub fn deinit(self: *Factorization) void {
        self.dense_lu.deinit();
        self.sparse_lu.deinit();
        self.sparse_basis.deinit();
        self.allocator.free(self.eta_values);
        self.allocator.free(self.eta_rows);
        self.allocator.free(self.identity_basic);
        self.allocator.free(self.identity_scale);
        self.allocator.free(self.identity_sign);
    }
    pub fn factorize(self: *Factorization, n: usize, matrix_data: []const f64) FactorizationError!void {
        self.dense_lu.factorize(n, matrix_data) catch |err| return switch (err) {
            error.DimensionMismatch => error.DimensionMismatch,
            error.Singular => error.Singular,
            error.NumericalFailure => error.NumericalFailure,
            error.OutOfMemory => error.OutOfMemory,
        };
        try self.prepareEtaStorage(n);
        self.backend_kind = .dense_lu;
        self.stats.factorizations += 1;
    }

    /// Assemble and factorize the current simplex basis directly as CSC. The
    /// retained buffers make subsequent reinversions allocation-free whenever
    /// their previous capacities are sufficient.
    pub fn factorizeSparseBasis(
        self: *Factorization,
        problem_matrix: matrix.CscView,
        basic_index: []const u32,
        row_scale: []const f64,
        artificial_sign: []const f64,
    ) FactorizationError!void {
        const basis = self.sparse_basis.assemble(problem_matrix, basic_index, row_scale, artificial_sign) catch |err| return switch (err) {
            error.DimensionMismatch, error.InvalidBasisColumn => error.DimensionMismatch,
            error.CapacityOverflow, error.OutOfMemory => error.OutOfMemory,
            error.NonFiniteValue => error.NumericalFailure,
        };
        self.sparse_lu.factorizeAssumeValid(basis) catch |err| return switch (err) {
            error.DimensionMismatch, error.DimensionTooLarge => error.DimensionMismatch,
            error.Singular => error.Singular,
            error.NumericalFailure, error.InvalidBasis => error.NumericalFailure,
            error.OutOfMemory, error.CapacityOverflow => error.OutOfMemory,
        };
        self.prepareSparseUpdateStorage(problem_matrix.num_rows);
        self.backend_kind = .sparse_lu;
        self.stats.factorizations += 1;
    }

    /// Production basis reinversion policy. Small bases retain the dense
    /// oracle because its contiguous kernel is cheaper than sparse metadata;
    /// larger bases assemble CSC directly for SparseLU.
    pub fn factorizeBasis(
        self: *Factorization,
        problem_matrix: matrix.CscView,
        basic_index: []const u32,
        row_scale: []const f64,
        artificial_sign: []const f64,
    ) FactorizationError!void {
        const n = problem_matrix.num_rows;
        if (n >= 64)
            return self.factorizeSparseBasis(problem_matrix, basic_index, row_scale, artificial_sign);
        if (self.dense_lu.n != n or self.dense_lu.lu.len != n * n) return error.DimensionMismatch;
        const buffer = self.dense_lu.lu;
        @memset(buffer, 0.0);
        if (basic_index.len != n or row_scale.len != n or artificial_sign.len != n) return error.DimensionMismatch;
        for (basic_index, 0..) |global_column_u32, basis_column| {
            const global_column: usize = global_column_u32;
            if (global_column < problem_matrix.num_cols) {
                const begin = problem_matrix.col_starts[global_column];
                const end = problem_matrix.col_starts[global_column + 1];
                for (problem_matrix.row_indices[begin..end], problem_matrix.values[begin..end]) |row, coefficient| {
                    const row_index = row.toUsize();
                    const scaled = coefficient * row_scale[row_index];
                    if (!std.math.isFinite(scaled)) return error.NumericalFailure;
                    buffer[row_index * n + basis_column] = scaled;
                }
            } else {
                const internal = global_column - problem_matrix.num_cols;
                if (internal >= 2 * n) return error.DimensionMismatch;
                const row = if (internal < n) internal else internal - n;
                const value = if (internal < n) 1.0 else artificial_sign[row];
                if (!std.math.isFinite(value) or value == 0.0) return error.NumericalFailure;
                buffer[row * n + basis_column] = value;
            }
        }
        self.backend_kind = .dense_lu;
        return self.refactorizeInPlace();
    }

    /// Build and factorize an identity basis without an intermediate matrix
    /// copy. Used by the canonical slack-basis crash initializer.
    pub fn factorizeIdentity(self: *Factorization, n: usize) FactorizationError!void {
        if (n >= 64) {
            self.identity_basic = self.allocator.realloc(self.identity_basic, n) catch return error.OutOfMemory;
            self.identity_scale = self.allocator.realloc(self.identity_scale, n) catch return error.OutOfMemory;
            self.identity_sign = self.allocator.realloc(self.identity_sign, n) catch return error.OutOfMemory;
            for (self.identity_basic, 0..) |*column, index| column.* = @intCast(index);
            @memset(self.identity_scale, 1.0);
            @memset(self.identity_sign, 1.0);
            const empty_matrix = matrix.CscView.initAssumeValid(n, 0, &[_]usize{0}, &.{}, &.{});
            return self.factorizeSparseBasis(empty_matrix, self.identity_basic, self.identity_scale, self.identity_sign);
        }
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
        self.backend_kind = .dense_lu;
        self.stats.factorizations += 1;
    }

    pub fn mutableBasisBuffer(self: *Factorization, n: usize) FactorizationError![]f64 {
        if (self.dense_lu.n != n or self.dense_lu.lu.len != n * n) return error.DimensionMismatch;
        return self.dense_lu.lu;
    }

    pub fn refactorizeInPlace(self: *Factorization) FactorizationError!void {
        if (self.backend_kind != .dense_lu) return error.NotImplemented;
        self.dense_lu.refactorizeInPlace() catch |err| return switch (err) {
            error.DimensionMismatch => error.DimensionMismatch,
            error.Singular => error.Singular,
            error.NumericalFailure => error.NumericalFailure,
            error.OutOfMemory => error.OutOfMemory,
        };
        self.update_count = 0;
        self.eta_count = 0;
        self.maximum_update_growth = 1.0;
        self.stats.factorizations += 1;
    }
    pub fn solve(self: *Factorization, rhs: []f64) FactorizationError!void {
        self.stats.ftran_calls += 1;
        switch (self.backend_kind) {
            .dense_lu => self.dense_lu.solve(rhs) catch |err| return switch (err) {
                error.DimensionMismatch => error.DimensionMismatch,
                error.Singular => error.Singular,
                error.NumericalFailure => error.NumericalFailure,
                error.OutOfMemory => error.OutOfMemory,
            },
            .sparse_lu => self.sparse_lu.solve(rhs) catch |err| return switch (err) {
                error.DimensionMismatch, error.DimensionTooLarge => error.DimensionMismatch,
                error.Singular => error.Singular,
                error.NumericalFailure, error.InvalidBasis => error.NumericalFailure,
                error.OutOfMemory, error.CapacityOverflow => error.OutOfMemory,
            },
        }
        for (0..self.eta_count) |update_index| self.applyEtaInverse(update_index, rhs);
    }

    /// FTRAN used for an entering column. SparseLU retains the partial `aq`
    /// spike required by its next Forrest--Tomlin update; dense LU continues
    /// to use the product-form Eta path.
    pub fn solveForUpdate(self: *Factorization, rhs: []f64) FactorizationError!void {
        if (self.backend_kind == .dense_lu) return self.solve(rhs);
        self.stats.ftran_calls += 1;
        self.sparse_lu.solveForUpdate(rhs) catch |err| return switch (err) {
            error.DimensionMismatch, error.DimensionTooLarge => error.DimensionMismatch,
            error.Singular => error.Singular,
            error.NumericalFailure, error.InvalidBasis => error.NumericalFailure,
            error.OutOfMemory, error.CapacityOverflow => error.OutOfMemory,
        };
    }
    pub fn solveTranspose(self: *Factorization, rhs: []f64) FactorizationError!void {
        self.stats.btran_calls += 1;
        var update_index = self.eta_count;
        while (update_index > 0) {
            update_index -= 1;
            self.applyEtaInverseTranspose(update_index, rhs);
        }
        switch (self.backend_kind) {
            .dense_lu => self.dense_lu.solveTranspose(rhs) catch |err| return switch (err) {
                error.DimensionMismatch => error.DimensionMismatch,
                error.Singular => error.Singular,
                error.NumericalFailure => error.NumericalFailure,
                error.OutOfMemory => error.OutOfMemory,
            },
            .sparse_lu => self.sparse_lu.solveTranspose(rhs) catch |err| return switch (err) {
                error.DimensionMismatch, error.DimensionTooLarge => error.DimensionMismatch,
                error.Singular => error.Singular,
                error.NumericalFailure, error.InvalidBasis => error.NumericalFailure,
                error.OutOfMemory, error.CapacityOverflow => error.OutOfMemory,
            },
        }
    }
    pub fn update(self: *Factorization, update_view: PivotUpdateView) FactorizationError!void {
        if (update_view.direction.len != self.dimension or update_view.leaving_row >= self.dimension) return error.DimensionMismatch;
        if (self.backend_kind == .dense_lu and self.eta_count >= self.eta_capacity) return error.NotImplemented;
        const pivot_row: usize = @intCast(update_view.leaving_row);
        if (!std.math.isFinite(update_view.column_scale)) return error.Singular;
        const pivot_magnitude = @abs(update_view.direction[pivot_row] * update_view.column_scale);
        var maximum_magnitude: f64 = 0.0;
        for (update_view.direction) |value| {
            if (!std.math.isFinite(value)) return error.NumericalFailure;
            maximum_magnitude = @max(maximum_magnitude, @abs(value * update_view.column_scale));
        }
        if (pivot_magnitude <= self.dense_lu.pivot_tolerance)
            return error.Singular;
        const update_growth = maximum_magnitude / pivot_magnitude;
        if (!std.math.isFinite(update_growth)) return error.NumericalFailure;
        switch (self.backend_kind) {
            .dense_lu => {
                const eta = self.eta_values[self.eta_count * self.dimension ..][0..self.dimension];
                for (eta, update_view.direction) |*value, direction| value.* = direction * update_view.column_scale;
                self.eta_rows[self.eta_count] = update_view.leaving_row;
                self.eta_count += 1;
                self.stats.eta_updates += 1;
            },
            .sparse_lu => {
                self.sparse_lu.applyForrestTomlinUpdate(
                    update_view.leaving_row,
                    update_view.direction,
                    update_view.column_scale,
                ) catch |err| return switch (err) {
                    error.DimensionMismatch, error.DimensionTooLarge => error.DimensionMismatch,
                    error.Singular => error.Singular,
                    error.NumericalFailure, error.InvalidBasis => error.NumericalFailure,
                    error.OutOfMemory, error.CapacityOverflow => error.OutOfMemory,
                };
                self.stats.ft_updates += 1;
            },
        }
        self.update_count += 1;
        self.maximum_update_growth = @max(self.maximum_update_growth, update_growth);
        self.stats.maximum_update_growth = @max(self.stats.maximum_update_growth, update_growth);
    }

    /// Adapt the update-chain length to observed numerical growth. Benign
    /// updates retain the caller's hard limit; increasingly oblique entering
    /// columns request earlier reinversion before residual drift accumulates.
    pub fn recommendedUpdateLimit(self: Factorization, hard_limit: usize) usize {
        const storage_limit = if (self.backend_kind == .dense_lu) @min(hard_limit, self.eta_capacity) else hard_limit;
        if (self.maximum_update_growth >= 1e10) return @min(storage_limit, 8);
        if (self.maximum_update_growth >= 1e7) return @min(storage_limit, 16);
        if (self.maximum_update_growth >= 1e4) return @min(storage_limit, 32);
        if (self.maximum_update_growth >= 1e2) return @min(storage_limit, 64);
        return storage_limit;
    }
    pub fn reinversionReason(self: Factorization, hard_limit: usize) ?ReinversionReason {
        const recommended = self.recommendedUpdateLimit(hard_limit);
        if (self.update_count < recommended) return null;
        const storage_limit = if (self.backend_kind == .dense_lu) @min(hard_limit, self.eta_capacity) else hard_limit;
        return if (recommended < storage_limit) .update_growth else .update_limit;
    }
    pub fn needsRefactor(self: Factorization, limit: usize) bool {
        return self.reinversionReason(limit) != null;
    }
    pub fn recordReinversion(self: *Factorization, reason: ReinversionReason) void {
        switch (reason) {
            .update_limit => self.stats.update_limit_reinversions += 1,
            .update_growth => self.stats.update_growth_reinversions += 1,
            .solve_residual => self.stats.solve_residual_reinversions += 1,
        }
    }

    pub fn pivotConditionEstimate(self: *const Factorization) f64 {
        return switch (self.backend_kind) {
            .dense_lu => self.dense_lu.pivotConditionEstimate(),
            .sparse_lu => self.sparse_lu.pivotConditionEstimate(),
        };
    }

    fn prepareEtaStorage(self: *Factorization, dimension: usize) FactorizationError!void {
        self.eta_values = self.allocator.realloc(self.eta_values, dimension * self.eta_capacity) catch return error.OutOfMemory;
        self.eta_rows = self.allocator.realloc(self.eta_rows, self.eta_capacity) catch return error.OutOfMemory;
        self.dimension = dimension;
        self.eta_count = 0;
        self.update_count = 0;
        self.maximum_update_growth = 1.0;
    }

    fn prepareSparseUpdateStorage(self: *Factorization, dimension: usize) void {
        self.dimension = dimension;
        self.eta_count = 0;
        self.update_count = 0;
        self.maximum_update_growth = 1.0;
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

test "relative update growth tightens the reinversion interval" {
    var factorization = Factorization.init(std.testing.allocator);
    defer factorization.deinit();
    try factorization.factorizeIdentity(2);
    try factorization.update(.{
        .leaving_row = 0,
        .entering_col = 0,
        .direction = &[_]f64{ 1e-5, 1 },
    });
    try std.testing.expectApproxEqRel(@as(f64, 1e5), factorization.maximum_update_growth, 1e-12);
    try std.testing.expectEqual(@as(usize, 32), factorization.recommendedUpdateLimit(100));
    factorization.update_count = 31;
    try std.testing.expect(!factorization.needsRefactor(100));
    factorization.update_count = 32;
    try std.testing.expect(factorization.needsRefactor(100));
    try std.testing.expectEqual(ReinversionReason.update_growth, factorization.reinversionReason(100).?);
    factorization.recordReinversion(.update_growth);
    try std.testing.expectEqual(@as(usize, 1), factorization.stats.update_growth_reinversions);
    try factorization.refactorizeInPlace();
    try std.testing.expectEqual(@as(f64, 1.0), factorization.maximum_update_growth);
}

test "reinversion is requested before retained update storage is exhausted" {
    var factorization = Factorization.init(std.testing.allocator);
    defer factorization.deinit();
    try factorization.factorizeIdentity(1);
    factorization.update_count = factorization.eta_capacity;
    try std.testing.expectEqual(factorization.eta_capacity, factorization.recommendedUpdateLimit(100));
    try std.testing.expectEqual(ReinversionReason.update_limit, factorization.reinversionReason(100).?);
}

test "factorization update rejects non-finite direction entries" {
    var factorization = Factorization.init(std.testing.allocator);
    defer factorization.deinit();
    try factorization.factorizeIdentity(2);
    try std.testing.expectError(error.NumericalFailure, factorization.update(.{
        .leaving_row = 0,
        .entering_col = 0,
        .direction = &[_]f64{ 1, std.math.nan(f64) },
    }));
}

test "sparse backend factors an assembled simplex basis" {
    const foundation = @import("foundation");
    const rows = [_]foundation.RowId{
        foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1),
        foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1),
    };
    const problem_matrix = matrix.CscView.initAssumeValid(
        2,
        2,
        &[_]usize{ 0, 2, 4 },
        &rows,
        &[_]f64{ 4, 1, 2, 3 },
    );
    var factorization = Factorization.init(std.testing.allocator);
    defer factorization.deinit();
    try factorization.factorizeSparseBasis(problem_matrix, &[_]u32{ 0, 1 }, &[_]f64{ 1, 1 }, &[_]f64{ 1, 1 });
    try std.testing.expectEqual(BackendKind.sparse_lu, factorization.backend_kind);
    var rhs = [_]f64{ 6, 7 };
    try factorization.solve(&rhs);
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), rhs[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2.2), rhs[1], 1e-12);
    var transpose_rhs = [_]f64{ 5, 8 };
    try factorization.solveTranspose(&transpose_rhs);
    try std.testing.expectApproxEqAbs(@as(f64, 0.7), transpose_rhs[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2.2), transpose_rhs[1], 1e-12);
    try std.testing.expect(factorization.pivotConditionEstimate() >= 1.0);

    var aq = [_]f64{ 2, 1 };
    try factorization.solveForUpdate(&aq);
    try factorization.update(.{ .leaving_row = 0, .entering_col = 2, .direction = &aq });
    try std.testing.expectEqual(@as(usize, 1), factorization.stats.ft_updates);
    try std.testing.expectEqual(@as(usize, 0), factorization.eta_count);
    rhs = .{ 4, 7 };
    try factorization.solve(&rhs);
    // Updated basis [[2, 2], [1, 3]].
    try std.testing.expectApproxEqAbs(@as(f64, -0.5), rhs[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2.5), rhs[1], 1e-12);
}

test "large crash identity uses sparse backend without Eta slab" {
    var factorization = Factorization.init(std.testing.allocator);
    defer factorization.deinit();
    try factorization.factorizeIdentity(64);
    try std.testing.expectEqual(BackendKind.sparse_lu, factorization.backend_kind);
    try std.testing.expectEqual(@as(usize, 0), factorization.eta_values.len);
    var rhs = [_]f64{1.0} ** 64;
    try factorization.solve(&rhs);
    try std.testing.expectEqualSlices(f64, &[_]f64{1.0} ** 64, &rhs);
}
