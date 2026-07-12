//! Quadratic programming model (QP / MIQP / QCP / MIQCP).
//!
//! Combines a [`LinearModel`] with an optional objective Hessian and an
//! optional set of quadratic constraints.
//!
//! Classification is derived automatically:
//! - No quadratic constraints, no objective Hessian          → delegates to LinearModel
//! - Objective Hessian present, no quadratic constraints      → QP / MIQP
//! - At least one quadratic constraint                       → QCP / MIQCP
//!
//! Experimental API.  The composition model is stable; constraint batch
//! storage may be optimised from slice‑of‑structs to a SoA layout.

const std = @import("std");
const hessian_module = @import("hessian.zig");
const linear_model_module = @import("linear_model.zig");
const problem_class_module = @import("problem_class.zig");

const Hessian = hessian_module.Hessian;
const LinearModel = linear_model_module.LinearModel;
const ProblemClass = problem_class_module.ProblemClass;
const DomainClass = problem_class_module.DomainClass;
const ObjectiveClass = problem_class_module.ObjectiveClass;
const ConstraintClass = problem_class_module.ConstraintClass;
const classify = problem_class_module.classify;

// ── QuadraticConstraint ────────────────────────────────────────────────────

/// Experimental API: a single quadratic constraint.
///
/// Represents:
/// ```text
/// lower ≤ ½ xᵀ Q x + aᵀ x ≤ upper
/// ```
///
/// `linear` gives the coefficients `a` as (index, value) pairs.
/// `hessian` stores the lower triangle of `Q`.
///
/// All column indices in both `linear` and `hessian` are relative to the
/// parent model's column space and must be in `[0, model.num_cols)`.
pub const QuadraticConstraint = struct {
    allocator: std.mem.Allocator,
    /// Linear term indices (sorted, unique).
    linear_indices: []foundation.ColId,
    /// Linear term values (same length as `linear_indices`).
    linear_values: []f64,
    /// Lower triangle of the quadratic term.
    hessian: Hessian,
    /// Constraint lower bound.
    lower: f64,
    /// Constraint upper bound.
    upper: f64,

    const Self = @This();

    /// Releases all owned memory.
    pub fn deinit(self: *Self) void {
        const allocator = self.allocator;
        if (self.linear_indices.len > 0) allocator.free(self.linear_indices);
        if (self.linear_values.len > 0) allocator.free(self.linear_values);
        self.hessian.deinit();
        self.* = undefined;
    }

    /// Number of non‑zero linear coefficients.
    pub fn linearNnz(self: Self) usize {
        return self.linear_indices.len;
    }

    /// Structural and numerical validation.
    pub fn validate(self: *const Self, num_cols: usize) !void {
        if (self.linear_indices.len != self.linear_values.len)
            return error.DimensionMismatch;
        try self.hessian.validate();

        // Bounds: lower ≤ upper.
        if (self.lower > self.upper) return error.InvalidBounds;
        if (self.lower != -std.math.inf(f64) and !std.math.isFinite(self.lower))
            return error.NonFiniteValue;
        if (self.upper != std.math.inf(f64) and !std.math.isFinite(self.upper))
            return error.NonFiniteValue;

        // Linear indices: sorted, unique, in range.
        var prev: ?usize = null;
        for (self.linear_indices) |idx| {
            const col = idx.toUsize();
            if (col >= num_cols) return error.ColumnOutOfRange;
            if (prev) |p| {
                if (col <= p) return error.IndicesNotIncreasing;
            }
            prev = col;
        }

        // Linear values: all finite.
        for (self.linear_values) |v| {
            if (!std.math.isFinite(v)) return error.NonFiniteValue;
        }
    }
};

// ── QuadraticModel ─────────────────────────────────────────────────────────

/// Experimental API: owning QP / MIQP / QCP / MIQCP model.
///
/// Owns the linear data, one objective Hessian, and an array of quadratic
/// constraints.  The problem class is derived automatically — do not cache
/// a separate `is_mip` flag.
pub const QuadraticModel = struct {
    allocator: std.mem.Allocator,
    /// Underlying LP / MILP data (costs, bounds, constraint matrix).
    linear: LinearModel,
    /// Objective Hessian.  `null` means purely linear objective.
    objective_hessian: ?Hessian,
    /// Quadratic constraints.  Empty slice means no quadratic constraints.
    quadratic_constraints: []QuadraticConstraint,

    const Self = @This();

    /// Releases all owned resources in dependency order.
    pub fn deinit(self: *Self) void {
        // Quadratic constraints first (they may reference the linear model).
        for (self.quadratic_constraints) |*qc| qc.deinit();
        if (self.quadratic_constraints.len > 0)
            self.allocator.free(self.quadratic_constraints);
        // Objective Hessian.
        if (self.objective_hessian) |*h| h.deinit();
        // Linear model last.
        self.linear.deinit();
        self.* = undefined;
    }

    /// Number of rows (delegated to the linear model).
    pub fn numRows(self: *const Self) usize {
        return self.linear.num_rows;
    }

    /// Number of columns (delegated to the linear model).
    pub fn numCols(self: *const Self) usize {
        return self.linear.num_cols;
    }

    /// Domain class derived from linear model integrality.
    pub fn domainClass(self: *const Self) DomainClass {
        return self.linear.domainClass();
    }

    /// Objective class derived from Hessian presence.
    pub fn objectiveClass(self: *const Self) ObjectiveClass {
        return if (self.objective_hessian != null) .quadratic else .linear;
    }

    /// Constraint class derived from quadratic constraint presence.
    pub fn constraintClass(self: *const Self) ConstraintClass {
        return if (self.quadratic_constraints.len > 0) .quadratic else .linear;
    }

    /// Derived problem class (never cached — always computed).
    pub fn problemClass(self: *const Self) ProblemClass {
        return classify(
            self.domainClass(),
            self.objectiveClass(),
            self.constraintClass(),
        );
    }
};

// ── Error set ──────────────────────────────────────────────────────────────

pub const QuadraticError = error{
    DimensionMismatch,
    InvalidBounds,
    NonFiniteValue,
    ColumnOutOfRange,
    IndicesNotIncreasing,
};

// ── Tests ──────────────────────────────────────────────────────────────────

const foundation = @import("foundation");
const ColId = foundation.ColId;

test "QuadraticModel.empty linear model is LP" {
    var lm = try LinearModel.initEmpty(std.testing.allocator);
    errdefer lm.deinit();

    var qm = QuadraticModel{
        .allocator = std.testing.allocator,
        .linear = lm,
        .objective_hessian = null,
        .quadratic_constraints = &.{},
    };
    defer qm.deinit();

    try std.testing.expectEqual(@as(usize, 0), qm.numRows());
    try std.testing.expectEqual(@as(usize, 0), qm.numCols());
    try std.testing.expectEqual(ProblemClass.lp, qm.problemClass());
}

test "QuadraticModel.with objective Hessian is QP" {
    var lm = try LinearModel.initEmpty(std.testing.allocator);
    errdefer lm.deinit();

    var qm = QuadraticModel{
        .allocator = std.testing.allocator,
        .linear = lm,
        .objective_hessian = hessian_module.Hessian.initEmpty(std.testing.allocator) catch null,
        .quadratic_constraints = &.{},
    };
    defer qm.deinit();

    try std.testing.expectEqual(ProblemClass.qp, qm.problemClass());
}

test "QuadraticModel.with empty Hessian still classifies as QP" {
    var lm = try LinearModel.initEmpty(std.testing.allocator);
    errdefer lm.deinit();

    var qm = QuadraticModel{
        .allocator = std.testing.allocator,
        .linear = lm,
        .objective_hessian = Hessian.initEmpty(std.testing.allocator) catch null,
        .quadratic_constraints = &.{},
    };
    defer qm.deinit();

    // An empty Hessian is still "not null", so classification is QP even if
    // the quadratic contribution is empty. This is a minor edge case noted
    // in the issue tracker.
    try std.testing.expectEqual(ObjectiveClass.quadratic, qm.objectiveClass());
    try std.testing.expectEqual(ProblemClass.qp, qm.problemClass());
}

test "QuadraticConstraint.validate linear indices" {
    var lm = try LinearModel.initEmpty(std.testing.allocator);
    defer lm.deinit();

    const qc_alloc = std.testing.allocator;
    var qc = QuadraticConstraint{
        .allocator = qc_alloc,
        .linear_indices = try qc_alloc.dupe(ColId, &.{
            ColId.fromUsizeAssumeValid(0),
            ColId.fromUsizeAssumeValid(2),
        }),
        .linear_values = try qc_alloc.dupe(f64, &.{ 1.0, -1.5 }),
        .hessian = try Hessian.initEmpty(qc_alloc),
        .lower = -std.math.inf(f64),
        .upper = 10.0,
    };
    defer qc.deinit();
    try qc.validate(3);
}

test "QuadraticConstraint.rejects out-of-range column" {
    const qc_alloc = std.testing.allocator;
    var qc = QuadraticConstraint{
        .allocator = qc_alloc,
        .linear_indices = try qc_alloc.dupe(ColId, &.{
            ColId.fromUsizeAssumeValid(5),
        }),
        .linear_values = try qc_alloc.dupe(f64, &.{1.0}),
        .hessian = try Hessian.initEmpty(qc_alloc),
        .lower = -std.math.inf(f64),
        .upper = 10.0,
    };
    defer qc.deinit();
    try std.testing.expectError(error.ColumnOutOfRange, qc.validate(3));
}

test "QuadraticConstraint.default deinit after initEmpty" {
    const qc_alloc = std.testing.allocator;
    var qc = QuadraticConstraint{
        .allocator = qc_alloc,
        .linear_indices = &.{},
        .linear_values = &.{},
        .hessian = try Hessian.initEmpty(qc_alloc),
        .lower = -std.math.inf(f64),
        .upper = std.math.inf(f64),
    };
    defer qc.deinit();
    try qc.validate(0);
}
