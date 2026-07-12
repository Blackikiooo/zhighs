//! Incremental builder for [`LinearModel`](@import("linear_model.zig")).
//!
//! Collects costs, bounds, integrality, and sparse matrix coefficients via
//! triplet / coordinate form, then `freeze()` validates, canonicalises, and
//! returns an owning `LinearModel`.
//!
//! ## Usage
//! ```zig
//! var bld = try LinearModelBuilder.init(allocator, num_rows, num_cols);
//! defer bld.deinit();
//! try bld.setColCost(0, 1.0);
//! try bld.setColBounds(0, 0.0, 10.0);
//! try bld.appendCoefficient(0, 0, 2.0);  // row 0, col 0 = 2.0
//! var model = try bld.freeze();
//! defer model.deinit();
//! ```
//!
//! Experimental API.  The build-and-freeze lifecycle is stable; individual
//! setter signatures may be extended as presolve integration reveals needs.

const std = @import("std");
const types = @import("types.zig");
const problem_class = @import("problem_class.zig");
const linear_model = @import("linear_model.zig");
const matrix = @import("matrix");
const foundation = @import("foundation");

const INFINITY = types.INFINITY;
const ObjectiveSense = types.ObjectiveSense;
const RowId = foundation.RowId;
const ColId = foundation.ColId;
const MatrixBuilder = matrix.MatrixBuilder;
const CscMatrix = matrix.CscMatrix;
const MatrixStore = matrix.MatrixStore;
const Integrality = linear_model.Integrality;
const LinearModel = linear_model.LinearModel;
const ProblemClass = problem_class.ProblemClass;
const DomainClass = problem_class.DomainClass;

/// Zero-tolerance for the matrix builder freeze.
const ZERO_TOLERANCE: f64 = 1e-12;

/// Experimental API: incremental builder for LinearModel.
pub const LinearModelBuilder = struct {
    allocator: std.mem.Allocator,

    num_rows: usize,
    num_cols: usize,

    objective_sense: ObjectiveSense,
    objective_offset: f64,

    col_cost: []f64,
    col_lower: []f64,
    col_upper: []f64,

    row_lower: []f64,
    row_upper: []f64,

    /// `null` means all continuous.
    integrality: ?[]Integrality,

    /// Accumulates (row, col, value) triplets.
    matrix_builder: MatrixBuilder,

    /// After `freeze()` the model arrays are moved out and this flag prevents
    /// the destructor from freeing them (they are now owned by the model).
    frozen: bool,

    const Self = @This();

    // ── Construction ───────────────────────────────────────────────────────

    /// Creates a builder with `num_rows` rows and `num_cols` columns.
    ///
    /// Defaults:
    /// - col_cost = 0
    /// - col_lower = 0, col_upper = +∞
    /// - row_lower = −∞, row_upper = +∞
    /// - integrality = null (all continuous)
    /// - objective sense = minimise
    /// - objective offset = 0
    pub fn init(allocator: std.mem.Allocator, num_rows: usize, num_cols: usize) !Self {
        const col_cost = try allocator.alloc(f64, num_cols);
        errdefer allocator.free(col_cost);
        @memset(col_cost, 0.0);

        const col_lower = try allocator.alloc(f64, num_cols);
        errdefer allocator.free(col_lower);
        for (col_lower) |*v| v.* = 0.0;

        const col_upper = try allocator.alloc(f64, num_cols);
        errdefer allocator.free(col_upper);
        for (col_upper) |*v| v.* = INFINITY;

        const row_lower = try allocator.alloc(f64, num_rows);
        errdefer allocator.free(row_lower);
        for (row_lower) |*v| v.* = -INFINITY;

        const row_upper = try allocator.alloc(f64, num_rows);
        errdefer allocator.free(row_upper);
        for (row_upper) |*v| v.* = INFINITY;

        const mb = try MatrixBuilder.init(num_rows, num_cols);
        errdefer mb.deinit(allocator);

        return Self{
            .allocator = allocator,
            .num_rows = num_rows,
            .num_cols = num_cols,
            .objective_sense = .minimize,
            .objective_offset = 0.0,
            .col_cost = col_cost,
            .col_lower = col_lower,
            .col_upper = col_upper,
            .row_lower = row_lower,
            .row_upper = row_upper,
            .integrality = null,
            .matrix_builder = mb,
            .frozen = false,
        };
    }

    // ── Destruction ────────────────────────────────────────────────────────

    /// Releases builder-owned resources.
    ///
    /// After `freeze()` the model arrays are owned by the returned
    /// `LinearModel`; this destructor skips them but still cleans up
    /// the internal `MatrixBuilder`.
    pub fn deinit(self: *Self) void {
        if (!self.frozen) {
            self.allocator.free(self.col_cost);
            self.allocator.free(self.col_lower);
            self.allocator.free(self.col_upper);
            self.allocator.free(self.row_lower);
            self.allocator.free(self.row_upper);
            if (self.integrality) |int| self.allocator.free(int);
        }
        self.matrix_builder.deinit(self.allocator);
        self.* = undefined;
    }

    // ── Setters ────────────────────────────────────────────────────────────

    /// Experimental API.
    pub fn setObjectiveSense(self: *Self, sense: ObjectiveSense) void {
        self.objective_sense = sense;
    }

    /// Experimental API.
    pub fn setObjectiveOffset(self: *Self, offset: f64) void {
        self.objective_offset = offset;
    }

    /// Experimental API: set the cost coefficient for column `col`.
    pub fn setColCost(self: *Self, col: usize, cost: f64) !void {
        if (col >= self.num_cols) return error.IndexOutOfRange;
        if (!std.math.isFinite(cost)) return error.NonFiniteValue;
        self.col_cost[col] = cost;
    }

    /// Experimental API: set both bounds for column `col`.
    /// `lower` and `upper` must satisfy `lower ≤ upper` and be finite
    /// or ±∞ respectively.
    pub fn setColBounds(self: *Self, col: usize, lower: f64, upper: f64) !void {
        if (col >= self.num_cols) return error.IndexOutOfRange;
        if (lower > upper) return error.InvalidBounds;
        if (lower != -INFINITY and !std.math.isFinite(lower)) return error.NonFiniteValue;
        if (upper != INFINITY and !std.math.isFinite(upper)) return error.NonFiniteValue;
        self.col_lower[col] = lower;
        self.col_upper[col] = upper;
    }

    /// Experimental API: set the lower and upper bounds for row `row`.
    /// `lower` and `upper` must satisfy `lower ≤ upper`.
    pub fn setRowBounds(self: *Self, row: usize, lower: f64, upper: f64) !void {
        if (row >= self.num_rows) return error.IndexOutOfRange;
        if (lower > upper) return error.InvalidBounds;
        if (lower != -INFINITY and !std.math.isFinite(lower)) return error.NonFiniteValue;
        if (upper != INFINITY and !std.math.isFinite(upper)) return error.NonFiniteValue;
        self.row_lower[row] = lower;
        self.row_upper[row] = upper;
    }

    /// Experimental API: set integrality for column `col`.
    /// Calling this ensures `integrality` is allocated (initialised to
    /// `.continuous` for every column).
    pub fn setColIntegrality(self: *Self, col: usize, int: Integrality) !void {
        if (col >= self.num_cols) return error.IndexOutOfRange;
        // Lazily allocate integrality array on first call.
        if (self.integrality == null) {
            self.integrality = try self.allocator.alloc(Integrality, self.num_cols);
            @memset(self.integrality.?, .continuous);
        }
        self.integrality.?[col] = int;
    }

    /// Experimental API: batch-set all integrality at once.
    /// The slice must have length `num_cols`.
    pub fn setAllIntegrality(self: *Self, integrality: []const Integrality) !void {
        if (integrality.len != self.num_cols) return error.DimensionMismatch;
        if (self.integrality) |existing| {
            self.allocator.free(existing);
            self.integrality = null;
        }
        self.integrality = try self.allocator.dupe(Integrality, integrality);
    }

    /// Experimental API: append a matrix coefficient.
    /// Duplicate coordinates are merged during freeze; explicit zeros are
    /// removed if their absolute value falls below `ZERO_TOLERANCE`.
    pub fn appendCoefficient(self: *Self, row: usize, col: usize, value: f64) !void {
        const r = try RowId.fromUsize(row);
        const c = try ColId.fromUsize(col);
        try self.matrix_builder.append(self.allocator, r, c, value);
    }

    // ── Freeze ─────────────────────────────────────────────────────────────

    /// Validate, canonicalise, and return an owning `LinearModel`.
    ///
    /// After a successful call the builder is marked `frozen` and its
    /// destructor will not free the arrays that were moved into the model.
    ///
    /// On failure the builder state is unchanged (no leak).
    pub fn freeze(self: *Self) !LinearModel {
        // ── Pre-freeze validation ──────────────────────────────────────────
        try self.validateBuiltData();

        // ── Build the canonical CSC matrix ───────────────────────────────
        // MatrixBuilder.freeze sorts, merges duplicates, and removes
        // explicit zeros — the returned CscMatrix is structurally valid.
        var csc: CscMatrix = try self.matrix_builder.freeze(self.allocator, ZERO_TOLERANCE);
        errdefer csc.deinit(self.allocator);

        // Wrap in MatrixStore (no re-validation — builder guarantees it).
        const store = MatrixStore.initAssumeValid(csc);

        // ── Transfer ownership to LinearModel ──────────────────────────────
        const model = LinearModel{
            .allocator = self.allocator,
            .objective_sense = self.objective_sense,
            .objective_offset = self.objective_offset,
            .num_rows = self.num_rows,
            .num_cols = self.num_cols,
            .col_cost = self.col_cost,
            .col_lower = self.col_lower,
            .col_upper = self.col_upper,
            .row_lower = self.row_lower,
            .row_upper = self.row_upper,
            .integrality = self.integrality,
            .matrix = store,
            .revision = 1,
        };

        // Prevent deinit from freeing the moved arrays.
        self.frozen = true;

        return model;
    }

    // ── Internal validation ────────────────────────────────────────────────

    /// Checks array dimensions and value validity before building the matrix.
    fn validateBuiltData(self: *const Self) !void {
        if (self.col_cost.len != self.num_cols) return error.DimensionMismatch;
        if (self.col_lower.len != self.num_cols) return error.DimensionMismatch;
        if (self.col_upper.len != self.num_cols) return error.DimensionMismatch;
        if (self.row_lower.len != self.num_rows) return error.DimensionMismatch;
        if (self.row_upper.len != self.num_rows) return error.DimensionMismatch;

        // Validate all costs are finite.
        for (self.col_cost, 0..) |cost, col| {
            if (!std.math.isFinite(cost)) return error.NonFiniteObjective;
            _ = col;
        }

        // Validate bounds: lower ≤ upper, and any non-infinite bound is finite.
        for (0..self.num_cols) |col| {
            const lb = self.col_lower[col];
            const ub = self.col_upper[col];
            if (lb > ub) return error.InvalidBounds;
            if (lb != -INFINITY and !std.math.isFinite(lb)) return error.NonFiniteValue;
            if (ub != INFINITY and !std.math.isFinite(ub)) return error.NonFiniteValue;
        }
        for (0..self.num_rows) |row| {
            const lb = self.row_lower[row];
            const ub = self.row_upper[row];
            if (lb > ub) return error.InvalidBounds;
            if (lb != -INFINITY and !std.math.isFinite(lb)) return error.NonFiniteValue;
            if (ub != INFINITY and !std.math.isFinite(ub)) return error.NonFiniteValue;
        }

        // Validate integrality length if present.
        if (self.integrality) |int| {
            if (int.len != self.num_cols) return error.DimensionMismatch;
        }
    }
};

// ── Tests ──────────────────────────────────────────────────────────────────

test "LinearModelBuilder.init and deinit does not leak" {
    var bld = try LinearModelBuilder.init(std.testing.allocator, 2, 3);
    defer bld.deinit();
    try std.testing.expectEqual(@as(usize, 3), bld.num_cols);
    try std.testing.expectEqual(@as(usize, 2), bld.num_rows);
}

test "LinearModelBuilder.builds empty model" {
    var bld = try LinearModelBuilder.init(std.testing.allocator, 0, 0);
    defer bld.deinit();
    var lm = try bld.freeze();
    defer lm.deinit();
    try std.testing.expectEqual(@as(usize, 0), lm.num_rows);
    try std.testing.expectEqual(@as(usize, 0), lm.num_cols);
    try std.testing.expect(!lm.isMixedInteger());
    try std.testing.expectEqual(ProblemClass.lp, lm.problemClass());
}

test "LinearModelBuilder.buildSimpleLP" {
    var bld = try LinearModelBuilder.init(std.testing.allocator, 2, 3);
    defer bld.deinit();
    try bld.setColCost(0, 1.0);
    try bld.setColCost(1, 2.0);
    try bld.setColCost(2, 3.0);
    try bld.setColBounds(0, 0.0, 1.0);
    try bld.setColBounds(1, -5.0, 5.0);
    try bld.setColBounds(2, 0.0, INFINITY);
    try bld.setRowBounds(0, -INFINITY, 10.0);
    try bld.setRowBounds(1, 5.0, INFINITY);
    try bld.appendCoefficient(0, 0, 1.0);
    try bld.appendCoefficient(0, 1, 1.0);
    try bld.appendCoefficient(1, 1, -1.0);
    try bld.appendCoefficient(1, 2, 2.0);

    var lm = try bld.freeze();
    defer lm.deinit();

    try std.testing.expectEqual(@as(usize, 2), lm.num_rows);
    try std.testing.expectEqual(@as(usize, 3), lm.num_cols);
    try std.testing.expectEqual(@as(u64, 1), lm.revision);
    try std.testing.expect(!lm.isMixedInteger());
    try std.testing.expectEqual(ProblemClass.lp, lm.problemClass());
}

test "LinearModelBuilder.buildMILP" {
    var bld = try LinearModelBuilder.init(std.testing.allocator, 1, 2);
    defer bld.deinit();
    try bld.setColIntegrality(1, .integer);
    try bld.appendCoefficient(0, 0, 1.0);
    try bld.appendCoefficient(0, 1, 1.0);

    var lm = try bld.freeze();
    defer lm.deinit();

    try std.testing.expect(lm.isMixedInteger());
    try std.testing.expectEqual(DomainClass.mixed_integer, lm.domainClass());
    try std.testing.expectEqual(ProblemClass.milp, lm.problemClass());
    try std.testing.expectEqual(Integrality.continuous, lm.integralityAt(0));
    try std.testing.expectEqual(Integrality.integer, lm.integralityAt(1));
}

test "LinearModelBuilder.duplicate coefficients merged" {
    var bld = try LinearModelBuilder.init(std.testing.allocator, 1, 1);
    defer bld.deinit();
    try bld.appendCoefficient(0, 0, 2.0);
    try bld.appendCoefficient(0, 0, 3.0); // sums to 5.0

    var lm = try bld.freeze();
    defer lm.deinit();

    const csc = lm.matrix.csc();
    try std.testing.expectEqual(@as(usize, 1), csc.nnz());
    try std.testing.expectEqual(@as(f64, 5.0), csc.values[0]);
}

test "LinearModelBuilder.explicit zeros removed" {
    var bld = try LinearModelBuilder.init(std.testing.allocator, 1, 1);
    defer bld.deinit();
    try bld.appendCoefficient(0, 0, 1.0);
    try bld.appendCoefficient(0, 0, -1.0); // cancels to 0

    var lm = try bld.freeze();
    defer lm.deinit();

    const csc = lm.matrix.csc();
    try std.testing.expectEqual(@as(usize, 0), csc.nnz());
}

test "LinearModelBuilder.rejects invalid bounds" {
    var bld = try LinearModelBuilder.init(std.testing.allocator, 1, 1);
    defer bld.deinit();
    try std.testing.expectError(error.InvalidBounds, bld.setColBounds(0, 5.0, 3.0));
    try std.testing.expectError(error.InvalidBounds, bld.setRowBounds(0, 5.0, 3.0));
}

test "LinearModelBuilder.rejects non-finite cost" {
    var bld = try LinearModelBuilder.init(std.testing.allocator, 1, 1);
    defer bld.deinit();
    try std.testing.expectError(error.NonFiniteValue, bld.setColCost(0, std.math.nan(f64)));
    try std.testing.expectError(error.NonFiniteValue, bld.setColCost(0, std.math.inf(f64)));
}

test "LinearModelBuilder.rejects out-of-range column" {
    var bld = try LinearModelBuilder.init(std.testing.allocator, 1, 1);
    defer bld.deinit();
    try std.testing.expectError(error.IndexOutOfRange, bld.setColCost(1, 1.0));
    try std.testing.expectError(error.IndexOutOfRange, bld.setColBounds(1, 0.0, 1.0));
    try std.testing.expectError(error.IndexOutOfRange, bld.setColIntegrality(1, .integer));
}

test "LinearModelBuilder.freeze fails on dimension mismatch" {
    // Manually create an inconsistent builder state.
    var bld = try LinearModelBuilder.init(std.testing.allocator, 1, 1);
    defer bld.deinit();
    // Corrupt the length to trigger validation failure.
    bld.num_cols = 2;
    try std.testing.expectError(error.DimensionMismatch, bld.freeze());
}

test "LinearModelBuilder.setAllIntegrality" {
    var bld = try LinearModelBuilder.init(std.testing.allocator, 1, 3);
    defer bld.deinit();
    try bld.setAllIntegrality(&.{ .integer, .continuous, .semi_integer });
    var lm = try bld.freeze();
    defer lm.deinit();
    try std.testing.expect(lm.isMixedInteger());
    try std.testing.expectEqual(Integrality.integer, lm.integralityAt(0));
    try std.testing.expectEqual(Integrality.continuous, lm.integralityAt(1));
    try std.testing.expectEqual(Integrality.semi_integer, lm.integralityAt(2));
}

test "LinearModelBuilder.rejects setAllIntegrality length mismatch" {
    var bld = try LinearModelBuilder.init(std.testing.allocator, 1, 3);
    defer bld.deinit();
    try std.testing.expectError(error.DimensionMismatch, bld.setAllIntegrality(&.{.integer}));
}

test "LinearModelBuilder.builder can be deinitialized after freeze without double-free" {
    var bld = try LinearModelBuilder.init(std.testing.allocator, 1, 1);
    try bld.appendCoefficient(0, 0, 1.0);
    var lm = try bld.freeze();
    defer lm.deinit();
    // bld.deinit() is now safe because frozen=true.
    bld.deinit();
}
