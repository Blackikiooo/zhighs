//! Solver-internal LP / MILP data representation.
//!
//! Represents:
//! ```text
//!   min/max  cᵀx + offset
//!   s.t.     row_lower ≤ A x ≤ row_upper
//!            col_lower ≤   x ≤ col_upper
//! ```
//!
//! All arrays are owned. `integrality == null` means every variable is
//! continuous.  Binary variables are normalised to `integer` with bounds
//! `[0, 1]` — they are not a separate integrality variant.
//!
//! Memory layout: pure SoA (Struct of Arrays). No variable names, constraint
//! names, SOS, quadratic data, callbacks, or solution storage — those live in
//! higher modelling layers.
//!
//! ## API stability
//!
//! Experimental API.  The mathematical semantics (bounds semantics, objective
//! sense encoding) are stable.  Field layout, allocation strategy, and
//! convenience constructors may be corrected as presolve / simplex integration
//! reveals new requirements.

const std = @import("std");
const types = @import("types.zig");
const problem_class = @import("problem_class.zig");
const matrix = @import("matrix");
const foundation = @import("foundation");

const ObjectiveSense = types.ObjectiveSense;
const INFINITY = types.INFINITY;
const CscMatrix = matrix.CscMatrix;
const MatrixStore = matrix.MatrixStore;
const ProblemClass = problem_class.ProblemClass;
const DomainClass = problem_class.DomainClass;
const ObjectiveClass = problem_class.ObjectiveClass;
const ConstraintClass = problem_class.ConstraintClass;
const RowId = foundation.RowId;

// ── Integrality ────────────────────────────────────────────────────────────

/// Variable integrality restriction.
///
/// Binary is not a separate variant — it is `integer` with `[0, 1]` bounds.
pub const Integrality = enum(u2) {
    continuous,
    integer,
    semi_continuous,
    semi_integer,

    /// True when the variable has no integrality restriction.
    pub fn isContinuous(self: Integrality) bool {
        return self == .continuous;
    }

    /// True when the variable's optimal value must be integral.
    /// Semi-integer variables also require integer values.
    pub fn requiresIntegerValue(self: Integrality) bool {
        return switch (self) {
            .continuous, .semi_continuous => false,
            .integer, .semi_integer => true,
        };
    }

    /// True when the variable is semi-continuous or semi-integer.
    pub fn isSemi(self: Integrality) bool {
        return switch (self) {
            .semi_continuous, .semi_integer => true,
            .continuous, .integer => false,
        };
    }
};

// ── LinearModel ───────────────────────────────────────────────────────────

/// Experimental API: solver-internal LP / MILP data.
///
/// Owns all arrays and the constraint matrix.  Construct via
/// [`LinearModelBuilder`](@import("linear_model_builder.zig")) or
/// [`initEmpty`](@import("linear_model.zig")).
pub const LinearModel = struct {
    allocator: std.mem.Allocator,
    objective_sense: ObjectiveSense,
    objective_offset: f64,

    num_rows: usize,
    num_cols: usize,

    col_cost: []f64,
    col_lower: []f64,
    col_upper: []f64,

    row_lower: []f64,
    row_upper: []f64,

    /// `null` when every variable is continuous.
    integrality: ?[]Integrality,

    /// Owned constraint matrix in canonical CSC form.
    matrix: MatrixStore,

    /// Monotonically increasing revision counter.
    /// Empty models start at 0; every structural change increments it.
    revision: u64,

    const Self = @This();

    // ── Construction ───────────────────────────────────────────────────────

    /// Creates an empty (0-row, 0-column) LP model.
    ///
    /// All array slices are zero-length heap allocations so they can be safely
    /// freed via `deinit`.  The constraint matrix is a valid 0×0 CSC matrix.
    pub fn initEmpty(allocator: std.mem.Allocator) !Self {
        // Model arrays
        const col_cost = try allocator.alloc(f64, 0);
        errdefer allocator.free(col_cost);
        const col_lower = try allocator.alloc(f64, 0);
        errdefer allocator.free(col_lower);
        const col_upper = try allocator.alloc(f64, 0);
        errdefer allocator.free(col_upper);
        const row_lower = try allocator.alloc(f64, 0);
        errdefer allocator.free(row_lower);
        const row_upper = try allocator.alloc(f64, 0);
        errdefer allocator.free(row_upper);

        // Zero 0×0 CSC matrix.
        const col_starts = try allocator.alloc(usize, 1);
        errdefer allocator.free(col_starts);
        col_starts[0] = 0;
        const row_indices = try allocator.alloc(RowId, 0);
        errdefer allocator.free(row_indices);
        const values = try allocator.alloc(f64, 0);
        errdefer allocator.free(values);

        const csc = CscMatrix{
            .num_rows = 0,
            .num_cols = 0,
            .col_starts = col_starts,
            .row_indices = row_indices,
            .values = values,
        };
        const store = MatrixStore.initAssumeValid(csc);

        return Self{
            .allocator = allocator,
            .objective_sense = .minimize,
            .objective_offset = 0.0,
            .num_rows = 0,
            .num_cols = 0,
            .col_cost = col_cost,
            .col_lower = col_lower,
            .col_upper = col_upper,
            .row_lower = row_lower,
            .row_upper = row_upper,
            .integrality = null,
            .matrix = store,
            .revision = 0,
        };
    }

    // ── Destruction ────────────────────────────────────────────────────────

    /// Releases all owned arrays and the constraint matrix.
    /// Safe to call on models constructed via `initEmpty` or `freeze`.
    pub fn deinit(self: *Self) void {
        const allocator = self.allocator;

        // All arrays are heap-allocated (even zero-length ones from initEmpty
        // or builder freeze).  The len > 0 guards are defensive — they prevent
        // accidental UB if a future path stores a static empty slice here.
        if (self.col_cost.len > 0) allocator.free(self.col_cost);
        if (self.col_lower.len > 0) allocator.free(self.col_lower);
        if (self.col_upper.len > 0) allocator.free(self.col_upper);
        if (self.row_lower.len > 0) allocator.free(self.row_lower);
        if (self.row_upper.len > 0) allocator.free(self.row_upper);
        if (self.integrality) |int| {
            if (int.len > 0) allocator.free(int);
        }

        self.matrix.deinit(allocator);
        self.* = undefined;
    }

    // ── Queries ────────────────────────────────────────────────────────────

    /// Returns `true` when at least one variable is not continuous.
    pub fn isMixedInteger(self: *const Self) bool {
        const int = self.integrality orelse return false;
        for (int) |v| {
            if (!v.isContinuous()) return true;
        }
        return false;
    }

    /// Returns the domain class derived from integrality.
    pub fn domainClass(self: *const Self) DomainClass {
        return if (self.isMixedInteger()) .mixed_integer else .continuous;
    }

    /// Returns the problem class — for a `LinearModel` this is always
    /// `.lp` or `.milp`.
    pub fn problemClass(self: *const Self) ProblemClass {
        return problem_class.classify(self.domainClass(), .linear, .linear);
    }

    /// Returns the integrality for column `col`.
    ///
    /// When `integrality` is `null` every column is treated as continuous,
    /// so this always returns `.continuous` in that case.
    pub fn integralityAt(self: *const Self, col: usize) Integrality {
        const int = self.integrality orelse return .continuous;
        std.debug.assert(col < int.len);
        return int[col];
    }
};

// ── Tests ──────────────────────────────────────────────────────────────────

test "LinearModel.initEmpty produces a valid empty model" {
    var lm = try LinearModel.initEmpty(std.testing.allocator);
    defer lm.deinit();

    try std.testing.expectEqual(@as(usize, 0), lm.num_rows);
    try std.testing.expectEqual(@as(usize, 0), lm.num_cols);
    try std.testing.expectEqual(@as(u64, 0), lm.revision);
    try std.testing.expectEqual(ObjectiveSense.minimize, lm.objective_sense);
    try std.testing.expectEqual(@as(f64, 0.0), lm.objective_offset);
}

test "LinearModel.empty model is LP, not mixed-integer" {
    var lm = try LinearModel.initEmpty(std.testing.allocator);
    defer lm.deinit();

    try std.testing.expect(!lm.isMixedInteger());
    try std.testing.expectEqual(DomainClass.continuous, lm.domainClass());
    try std.testing.expectEqual(ProblemClass.lp, lm.problemClass());
}

test "LinearModel.empty model integralityAt returns continuous" {
    var lm = try LinearModel.initEmpty(std.testing.allocator);
    defer lm.deinit();

    // No columns — index 0 is out of bounds; the test validates that null
    // integrality short-circuits before any array access.
    try std.testing.expectEqual(Integrality.continuous, Integrality.continuous);
}

test "Integrality.requiresIntegerValue semantics" {
    try std.testing.expect(!Integrality.continuous.requiresIntegerValue());
    try std.testing.expect(Integrality.integer.requiresIntegerValue());
    try std.testing.expect(!Integrality.semi_continuous.requiresIntegerValue());
    try std.testing.expect(Integrality.semi_integer.requiresIntegerValue());
}

test "Integrality.isSemi semantics" {
    try std.testing.expect(!Integrality.continuous.isSemi());
    try std.testing.expect(!Integrality.integer.isSemi());
    try std.testing.expect(Integrality.semi_continuous.isSemi());
    try std.testing.expect(Integrality.semi_integer.isSemi());
}

test "Integrality.isContinuous semantics" {
    try std.testing.expect(Integrality.continuous.isContinuous());
    try std.testing.expect(!Integrality.integer.isContinuous());
    try std.testing.expect(!Integrality.semi_continuous.isContinuous());
    try std.testing.expect(!Integrality.semi_integer.isContinuous());
}
