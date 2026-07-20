//! Basis factorization policy boundary.
//!
//! The first implementation may use dense LU for correctness tests; the
//! interface deliberately permits sparse LU, eta/Forrest--Tomlin updates, and
//! periodic reinversion without changing the simplex engine.

const std = @import("std");
const matrix = @import("matrix");

/// Errors raised by factorization and update operations.
pub const FactorizationError = error{ DimensionMismatch, NotImplemented, Singular, NumericalFailure, OutOfMemory };

/// Selected base factorization. Updates are applied above either backend and
/// are cleared whenever a new base factorization is installed.
pub const BackendKind = enum { dense_lu, sparse_lu };

/// Reason a reinversion was triggered. Reported via `recordReinversion`.
pub const ReinversionReason = enum { update_limit, update_growth, solve_residual, small_pivot };

/// Classification of an update failure for diagnostic counters.
pub const UpdateFailureKind = enum { dimension_mismatch, unsupported, singular, numerical, out_of_memory };

/// Aggregate statistics across one solve. Counts are always populated; the
/// `*_ns` timing fields require an explicit `statistics_io` clock to be set.
pub const FactorizationStats = struct {
    factorizations: usize = 0, // Total base factorizations built
    ftran_calls: usize = 0, // Forward solves (B^-1 * x)
    btran_calls: usize = 0, // Transpose solves (B^-T * x)
    eta_updates: usize = 0, // Product-form Eta updates applied (dense backend)
    ft_updates: usize = 0, // Forrest-Tomlin updates applied (sparse backend)
    maximum_update_growth: f64 = 1.0, // Worst max/diagonal ratio observed
    update_limit_reinversions: usize = 0,
    update_growth_reinversions: usize = 0,
    solve_residual_reinversions: usize = 0,
    small_pivot_reinversions: usize = 0,
    maximum_update_count: usize = 0, // Peak update chain length seen
    invert_ns: u64 = 0, // Time spent inside factorize*
    ftran_ns: u64 = 0, // Time spent inside solve/solveForUpdate
    btran_ns: u64 = 0, // Time spent inside solveTranspose
    update_ns: u64 = 0, // Time spent inside update
    ftran_rhs_samples: usize = 0,
    ftran_rhs_nonzeros: usize = 0, // Total nonzeros across sampled FTRAN RHS vectors
    btran_rhs_samples: usize = 0,
    btran_rhs_nonzeros: usize = 0,
    dense_ftran_dispatches: usize = 0,
    dense_btran_dispatches: usize = 0,
    hyper_ftran_dispatches: usize = 0, // Reserved for hyper-sparse dispatches
    hyper_btran_dispatches: usize = 0,
    /// Sparse-index dispatches where the caller provided nonzero positions
    /// and the sparse LU backend used adaptive (hyper-sparse vs dense) solve.
    sparse_ftran_dispatches: usize = 0,
    sparse_btran_dispatches: usize = 0,
    update_dimension_failures: usize = 0,
    update_unsupported_failures: usize = 0,
    update_singular_failures: usize = 0,
    update_numerical_failures: usize = 0,
    update_out_of_memory_failures: usize = 0,
    dense_update_failures: usize = 0,
    sparse_update_failures: usize = 0,
};

/// View of one pivot update supplied by the simplex engine.
pub const PivotUpdateView = struct {
    leaving_row: u32,
    entering_col: u32,
    direction: []const f64, // Updated tableau column (B^-1 * a_entering)
    /// Converts a signed movement direction back to the actual entering
    /// basis column. Primal simplex uses -1 when entering from an upper bound.
    column_scale: f64 = 1.0,
};

pub const Factorization = struct {
    allocator: std.mem.Allocator,
    update_count: usize = 0, // Updates applied since last factorize
    dense_lu: matrix.DenseLU, // Dense LU backend (small bases)
    sparse_lu: matrix.SparseLU, // Sparse LU backend (larger bases)
    sparse_basis: matrix.SparseBasisBuffers, // Reusable CSC assembly buffers
    eta_values: []f64 = &.{}, // Slab of Eta vectors (dense backend only)
    eta_rows: []u32 = &.{}, // Leaving-row index for each Eta
    identity_basic: []u32 = &.{}, // Identity basis head for `factorizeIdentity`
    identity_scale: []f64 = &.{}, // Column scales for the identity basis
    identity_sign: []f64 = &.{}, // Sign convention for artificial identity columns
    eta_count: usize = 0, // Live Eta vectors in the slab
    eta_capacity: usize = 64, // Max Eta vectors before forced refactor (dense backend)
    dimension: usize = 0, // Current basis dimension (n = num_rows)
    backend_kind: BackendKind = .dense_lu, // Active backend for the current basis
    sparse_dimension_threshold: usize = 64, // Bases at or above this use the sparse backend
    stats: FactorizationStats = .{},
    /// Largest max(|d|)/|d[p]| observed since reinversion. This inexpensive
    /// signal estimates how strongly an update can amplify solve error.
    maximum_update_growth: f64 = 1.0,
    /// Optional solve-owned clock. It remains null in production hot paths
    /// unless the caller explicitly requests benchmark statistics.
    statistics_io: ?std.Io = null,

    /// Construct an empty factorization. Backends are lazily sized on first use.
    pub fn init(allocator: std.mem.Allocator) Factorization {
        return .{
            .allocator = allocator,
            .dense_lu = matrix.DenseLU.init(allocator),
            .sparse_lu = matrix.SparseLU.init(allocator),
            .sparse_basis = matrix.SparseBasisBuffers.init(allocator),
        };
    }

    /// Release all backend and Eta storage.
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

    /// Factorize a dense `n x n` matrix in row-major order via the dense backend.
    pub fn factorize(self: *Factorization, n: usize, matrix_data: []const f64) FactorizationError!void {
        const started = self.statisticsTimestamp();
        defer self.recordElapsed(&self.stats.invert_ns, started);
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
        column_scale: []const f64,
        artificial_sign: []const f64,
    ) FactorizationError!void {
        const started = self.statisticsTimestamp();
        defer self.recordElapsed(&self.stats.invert_ns, started);
        const basis = self.sparse_basis.assemble(problem_matrix, basic_index, row_scale, column_scale, artificial_sign) catch |err| return switch (err) {
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
        column_scale: []const f64,
        artificial_sign: []const f64,
    ) FactorizationError!void {
        const n = problem_matrix.num_rows;
        if (n >= self.sparse_dimension_threshold)
            return self.factorizeSparseBasis(problem_matrix, basic_index, row_scale, column_scale, artificial_sign);
        if (self.dense_lu.n != n or self.dense_lu.lu.len != n * n) return error.DimensionMismatch;
        const buffer = self.dense_lu.lu;
        @memset(buffer, 0.0);
        if (basic_index.len != n or row_scale.len != n or column_scale.len != problem_matrix.num_cols or artificial_sign.len != n) return error.DimensionMismatch;
        // Scatter the basis columns into the dense buffer. Structural columns
        // come from the problem matrix; logical/artificial columns are unit
        // vectors with an optional sign convention.
        for (basic_index, 0..) |global_column_u32, basis_column| {
            const global_column: usize = global_column_u32;
            if (global_column < problem_matrix.num_cols) {
                const begin = problem_matrix.col_starts[global_column];
                const end = problem_matrix.col_starts[global_column + 1];
                for (problem_matrix.row_indices[begin..end], problem_matrix.values[begin..end]) |row, coefficient| {
                    if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                    const row_index = row.toUsize();
                    const scaled = coefficient * row_scale[row_index] * column_scale[global_column];
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
        if (n >= self.sparse_dimension_threshold) {
            self.identity_basic = self.allocator.realloc(self.identity_basic, n) catch return error.OutOfMemory;
            self.identity_scale = self.allocator.realloc(self.identity_scale, n) catch return error.OutOfMemory;
            self.identity_sign = self.allocator.realloc(self.identity_sign, n) catch return error.OutOfMemory;
            for (self.identity_basic, 0..) |*column, index| column.* = @intCast(index);
            @memset(self.identity_scale, 1.0);
            @memset(self.identity_sign, 1.0);
            const empty_matrix = matrix.CscView.initAssumeValid(n, 0, &[_]usize{0}, &.{}, &.{});
            return self.factorizeSparseBasis(empty_matrix, self.identity_basic, self.identity_scale, &.{}, self.identity_sign);
        }
        const started = self.statisticsTimestamp();
        defer self.recordElapsed(&self.stats.invert_ns, started);
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

    /// Expose the dense LU backing buffer for in-place fills by the caller.
    pub fn mutableBasisBuffer(self: *Factorization, n: usize) FactorizationError![]f64 {
        if (self.dense_lu.n != n or self.dense_lu.lu.len != n * n) return error.DimensionMismatch;
        return self.dense_lu.lu;
    }

    /// Refactorize the dense LU in place, using the existing buffer contents.
    /// Clears all Eta updates.
    pub fn refactorizeInPlace(self: *Factorization) FactorizationError!void {
        const started = self.statisticsTimestamp();
        defer self.recordElapsed(&self.stats.invert_ns, started);
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

    /// Forward solve `rhs = B^-1 * rhs` in place, applying any pending Eta updates.
    pub fn solve(self: *Factorization, rhs: []f64) FactorizationError!void {
        const started = self.statisticsTimestamp();
        defer self.recordElapsed(&self.stats.ftran_ns, started);
        self.stats.ftran_calls += 1;
        self.observeFtranRhs(rhs);
        self.stats.dense_ftran_dispatches += 1;
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
        // Apply pending Eta updates (dense backend only).
        for (0..self.eta_count) |update_index| self.applyEtaInverse(update_index, rhs);
    }

    /// FTRAN with caller-provided nonzero positions for the sparse-index
    /// adaptive dispatch. Pass `null` input_indices for the standard dense
    /// solve (e.g. iterative refinement residuals). The returned RHS is
    /// always dense.
    pub fn solveSparse(self: *Factorization, rhs: []f64, input_indices: []const u32) FactorizationError!void {
        const started = self.statisticsTimestamp();
        defer self.recordElapsed(&self.stats.ftran_ns, started);
        self.stats.ftran_calls += 1;
        self.observeFtranRhs(rhs);
        if (self.backend_kind == .sparse_lu and input_indices.len * 8 < self.dimension) {
            self.stats.hyper_ftran_dispatches += 1;
            self.stats.sparse_ftran_dispatches += 1;
            self.sparse_lu.solveAdaptive(rhs, input_indices, false) catch |err| return switch (err) {
                error.DimensionMismatch, error.DimensionTooLarge => error.DimensionMismatch,
                error.Singular => error.Singular,
                error.NumericalFailure, error.InvalidBasis => error.NumericalFailure,
                error.OutOfMemory, error.CapacityOverflow => error.OutOfMemory,
            };
        } else {
            self.stats.dense_ftran_dispatches += 1;
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
        }
        for (0..self.eta_count) |update_index| self.applyEtaInverse(update_index, rhs);
    }

    pub fn solveTransposeSparse(self: *Factorization, rhs: []f64, input_indices: []const u32) FactorizationError!void {
        const started = self.statisticsTimestamp();
        defer self.recordElapsed(&self.stats.btran_ns, started);
        self.stats.btran_calls += 1;
        self.observeBtranRhs(rhs);
        var update_index = self.eta_count;
        while (update_index > 0) {
            update_index -= 1;
            self.applyEtaInverseTranspose(update_index, rhs);
        }
        if (self.backend_kind == .sparse_lu and input_indices.len * 8 < self.dimension) {
            self.stats.hyper_btran_dispatches += 1;
            self.stats.sparse_btran_dispatches += 1;
            self.sparse_lu.solveAdaptive(rhs, input_indices, true) catch |err| return switch (err) {
                error.DimensionMismatch, error.DimensionTooLarge => error.DimensionMismatch,
                error.Singular => error.Singular,
                error.NumericalFailure, error.InvalidBasis => error.NumericalFailure,
                error.OutOfMemory, error.CapacityOverflow => error.OutOfMemory,
            };
        } else {
            self.stats.dense_btran_dispatches += 1;
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
    }

    /// FTRAN used for an entering column. SparseLU retains the partial `aq`
    /// spike required by its next Forrest--Tomlin update; dense LU continues
    /// to use the product-form Eta path.
    pub fn solveForUpdate(self: *Factorization, rhs: []f64) FactorizationError!void {
        if (self.backend_kind == .dense_lu) return self.solve(rhs);
        const started = self.statisticsTimestamp();
        defer self.recordElapsed(&self.stats.ftran_ns, started);
        self.stats.ftran_calls += 1;
        self.observeFtranRhs(rhs);
        self.stats.dense_ftran_dispatches += 1;
        self.sparse_lu.solveForUpdate(rhs) catch |err| return switch (err) {
            error.DimensionMismatch, error.DimensionTooLarge => error.DimensionMismatch,
            error.Singular => error.Singular,
            error.NumericalFailure, error.InvalidBasis => error.NumericalFailure,
            error.OutOfMemory, error.CapacityOverflow => error.OutOfMemory,
        };
    }

    /// Transpose solve `rhs = B^-T * rhs` in place. Eta updates are applied
    /// in reverse order before the backend transpose solve.
    pub fn solveTranspose(self: *Factorization, rhs: []f64) FactorizationError!void {
        const started = self.statisticsTimestamp();
        defer self.recordElapsed(&self.stats.btran_ns, started);
        self.stats.btran_calls += 1;
        self.observeBtranRhs(rhs);
        self.stats.dense_btran_dispatches += 1;
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

    /// Apply one pivot update to the factorization. Records growth statistics
    /// and dispatches to either the Eta (dense) or Forrest-Tomlin (sparse) path.
    pub fn update(self: *Factorization, update_view: PivotUpdateView) FactorizationError!void {
        const started = self.statisticsTimestamp();
        defer self.recordElapsed(&self.stats.update_ns, started);
        if (update_view.direction.len != self.dimension or update_view.leaving_row >= self.dimension) return error.DimensionMismatch;
        if (self.backend_kind == .dense_lu and self.eta_count >= self.eta_capacity) return error.NotImplemented;
        const pivot_row: usize = @intCast(update_view.leaving_row);
        if (!std.math.isFinite(update_view.column_scale)) return error.Singular;
        const pivot_magnitude = @abs(update_view.direction[pivot_row] * update_view.column_scale);
        // Scan the direction vector for non-finite entries and the largest magnitude.
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
                // Append a new Eta vector to the slab.
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
        self.stats.maximum_update_count = @max(self.stats.maximum_update_count, self.update_count);
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

    /// Return the reason a reinversion is currently required, if any.
    pub fn reinversionReason(self: Factorization, hard_limit: usize) ?ReinversionReason {
        const recommended = self.recommendedUpdateLimit(hard_limit);
        if (self.update_count < recommended) return null;
        const storage_limit = if (self.backend_kind == .dense_lu) @min(hard_limit, self.eta_capacity) else hard_limit;
        return if (recommended < storage_limit) .update_growth else .update_limit;
    }

    /// Convenience wrapper: true if a reinversion is required.
    pub fn needsRefactor(self: Factorization, limit: usize) bool {
        return self.reinversionReason(limit) != null;
    }

    /// Record that a reinversion was performed for `reason`. Called by the
    /// engine after a successful refactor.
    pub fn recordReinversion(self: *Factorization, reason: ReinversionReason) void {
        switch (reason) {
            .update_limit => self.stats.update_limit_reinversions += 1,
            .update_growth => self.stats.update_growth_reinversions += 1,
            .solve_residual => self.stats.solve_residual_reinversions += 1,
            .small_pivot => self.stats.small_pivot_reinversions += 1,
        }
    }

    /// Classify and record a failed update for diagnostic counters.
    pub fn recordUpdateFailure(self: *Factorization, err: FactorizationError) UpdateFailureKind {
        switch (self.backend_kind) {
            .dense_lu => self.stats.dense_update_failures += 1,
            .sparse_lu => self.stats.sparse_update_failures += 1,
        }
        return switch (err) {
            error.DimensionMismatch => blk: {
                self.stats.update_dimension_failures += 1;
                break :blk .dimension_mismatch;
            },
            error.NotImplemented => blk: {
                self.stats.update_unsupported_failures += 1;
                break :blk .unsupported;
            },
            error.Singular => blk: {
                self.stats.update_singular_failures += 1;
                break :blk .singular;
            },
            error.NumericalFailure => blk: {
                self.stats.update_numerical_failures += 1;
                break :blk .numerical;
            },
            error.OutOfMemory => blk: {
                self.stats.update_out_of_memory_failures += 1;
                break :blk .out_of_memory;
            },
        };
    }

    /// Backend-specific pivot condition estimate (cheap heuristic, not a true
    /// condition number).
    pub fn pivotConditionEstimate(self: *const Factorization) f64 {
        return switch (self.backend_kind) {
            .dense_lu => self.dense_lu.pivotConditionEstimate(),
            .sparse_lu => self.sparse_lu.pivotConditionEstimate(),
        };
    }

    /// Bytes held by all retained factorization buffers.
    pub fn requestedBytes(self: *const Factorization) usize {
        var total = self.sparse_lu.requestedBytes();
        total += std.mem.sliceAsBytes(self.dense_lu.lu).len;
        total += std.mem.sliceAsBytes(self.dense_lu.pivots).len;
        total += std.mem.sliceAsBytes(self.dense_lu.work).len;
        total += std.mem.sliceAsBytes(self.sparse_basis.starts).len;
        total += std.mem.sliceAsBytes(self.sparse_basis.rows).len;
        total += std.mem.sliceAsBytes(self.sparse_basis.values).len;
        total += std.mem.sliceAsBytes(self.eta_values).len;
        total += std.mem.sliceAsBytes(self.eta_rows).len;
        total += std.mem.sliceAsBytes(self.identity_basic).len;
        total += std.mem.sliceAsBytes(self.identity_scale).len;
        total += std.mem.sliceAsBytes(self.identity_sign).len;
        return total;
    }

    /// Reset per-solve counters and enable timing only when an I/O clock is
    /// supplied. Call counts remain available even when timing is disabled.
    pub fn resetStatistics(self: *Factorization, io: ?std.Io) void {
        self.stats = .{};
        self.statistics_io = io;
    }

    /// Read the current timer, or return null when timing is disabled.
    fn statisticsTimestamp(self: *const Factorization) ?i96 {
        const io = self.statistics_io orelse return null;
        return std.Io.Clock.awake.now(io).nanoseconds;
    }

    /// Accumulate elapsed time into `target` when timing is enabled.
    fn recordElapsed(self: *const Factorization, target: *u64, started: ?i96) void {
        const begin = started orelse return;
        const io = self.statistics_io orelse return;
        const end = std.Io.Clock.awake.now(io).nanoseconds;
        if (end <= begin) return;
        target.* = std.math.add(u64, target.*, @intCast(end - begin)) catch std.math.maxInt(u64);
    }

    /// Sample FTRAN RHS density for hyper-sparsity heuristics.
    fn observeFtranRhs(self: *Factorization, rhs: []const f64) void {
        if (self.statistics_io == null) return;
        self.stats.ftran_rhs_samples += 1;
        for (rhs) |value| if (value != 0.0) {
            self.stats.ftran_rhs_nonzeros += 1;
        };
    }

    /// Sample BTRAN RHS density for hyper-sparsity heuristics.
    fn observeBtranRhs(self: *Factorization, rhs: []const f64) void {
        if (self.statistics_io == null) return;
        self.stats.btran_rhs_samples += 1;
        for (rhs) |value| if (value != 0.0) {
            self.stats.btran_rhs_nonzeros += 1;
        };
    }

    /// Allocate (or reallocate) the dense-backend Eta slab for `dimension` rows.
    fn prepareEtaStorage(self: *Factorization, dimension: usize) FactorizationError!void {
        self.eta_values = self.allocator.realloc(self.eta_values, dimension * self.eta_capacity) catch return error.OutOfMemory;
        self.eta_rows = self.allocator.realloc(self.eta_rows, self.eta_capacity) catch return error.OutOfMemory;
        self.dimension = dimension;
        self.eta_count = 0;
        self.update_count = 0;
        self.maximum_update_growth = 1.0;
    }

    /// Reset update counters for the sparse backend (no Eta slab is needed).
    fn prepareSparseUpdateStorage(self: *Factorization, dimension: usize) void {
        self.dimension = dimension;
        self.eta_count = 0;
        self.update_count = 0;
        self.maximum_update_growth = 1.0;
    }

    /// Apply the `update_index`-th Eta to a forward-solved vector: solve
    /// `(I + e_p * eta^T / eta_p) * x = rhs` in place.
    fn applyEtaInverse(self: *const Factorization, update_index: usize, vector: []f64) void {
        const row: usize = @intCast(self.eta_rows[update_index]);
        const eta = self.eta_values[update_index * self.dimension ..][0..self.dimension];
        const pivot_value = vector[row] / eta[row];
        for (vector, eta, 0..) |*value, coefficient, i| {
            if (i != row) value.* -= coefficient * pivot_value;
        }
        vector[row] = pivot_value;
    }

    /// Apply the transpose of `applyEtaInverse` for BTran.
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
    try factorization.factorizeSparseBasis(problem_matrix, &[_]u32{ 0, 1 }, &[_]f64{ 1, 1 }, &[_]f64{ 1, 1 }, &[_]f64{ 1, 1 });
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
