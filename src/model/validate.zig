//! Layered validation for the solver-internal model types.
//!
//! Error sets grow more specific at each layer so callers can distinguish
//! structural defects (dimensions, ownership) from numerical defects (NaN,
//! infinity) from semantic defects (duplicate, cycle, bound inconsistency).
//!
//! Each function here is a thin wrapper around the corresponding type's
//! intrinsic `validate` method, composing checks when a model aggregates
//! multiple types.
//!
//! Experimental API.  Validator signatures are stable; error sets may acquire
//! new members as additional model types are validated.

const std = @import("std");
const types = @import("types.zig");
const linear_model = @import("linear_model.zig");
const hessian_module = @import("hessian.zig");
const quadratic_model_module = @import("quadratic_model.zig");
const expression_graph_module = @import("expression_graph.zig");
const nonlinear_model_module = @import("nonlinear_model.zig");
const compiled_model_module = @import("compiled_model.zig");

const INFINITY = types.INFINITY;
const LinearModel = linear_model.LinearModel;
const Integrality = linear_model.Integrality;
const Hessian = hessian_module.Hessian;
const QuadraticConstraint = quadratic_model_module.QuadraticConstraint;
const QuadraticModel = quadratic_model_module.QuadraticModel;
const ExpressionGraph = expression_graph_module.ExpressionGraph;
const NonlinearModel = nonlinear_model_module.NonlinearModel;
const NonlinearConstraint = nonlinear_model_module.NonlinearConstraint;
const CompiledModel = compiled_model_module.CompiledModel;

// ── Error set ──────────────────────────────────────────────────────────────

/// Union of all validation errors across model layers.
pub const ValidationError = error{
    DimensionMismatch,
    InvalidBounds,
    InvalidMatrix,
    NonFiniteObjective,
    NonFiniteCoefficient,
    InvalidIntegrality,
    InvalidHessian,
    InvalidQuadraticConstraint,
    InvalidExpressionNode,
    InvalidExpressionRoot,
    CyclicExpression,
    VariableIndexOutOfRange,
    EmptyNonlinearModel,
    NonFiniteValue,
    IndicesNotIncreasing,
    ColumnOutOfRange,
    IndexOutOfBounds,
};

// ── LinearModel validation ─────────────────────────────────────────────────

/// Validates the numerical array contents of a `LinearModel`.
///
/// Checks:
/// - Matrix dimensions match model dimensions.
/// - All costs are finite.
/// - All bounds satisfy `lower ≤ upper` and are finite (or ±∞).
/// - Integrality length matches `num_cols`.
/// - Constraint matrix is structurally valid.
///
/// This function does **not** check convexity, degeneracy, or solver
/// compatibility — those are higher-level concerns.
pub fn validateLinearModel(model: *const LinearModel) ValidationError!void {
    const allocator = model.allocator;
    _ = allocator;

    // ── Matrix structural validity ─────────────────────────────────────────
    model.matrix.csc().validate() catch return error.InvalidMatrix;

    // ── Matrix dimensions ──────────────────────────────────────────────────
    const csc = model.matrix.csc();
    if (csc.num_rows != model.num_rows or csc.num_cols != model.num_cols)
        return error.DimensionMismatch;

    // ── Objective ──────────────────────────────────────────────────────────
    if (model.col_cost.len != model.num_cols) return error.DimensionMismatch;
    for (model.col_cost) |cost| {
        if (!std.math.isFinite(cost)) return error.NonFiniteObjective;
    }

    // ── Column bounds ──────────────────────────────────────────────────────
    if (model.col_lower.len != model.num_cols) return error.DimensionMismatch;
    if (model.col_upper.len != model.num_cols) return error.DimensionMismatch;
    for (0..model.num_cols) |col| {
        const lb = model.col_lower[col];
        const ub = model.col_upper[col];
        if (lb > ub) return error.InvalidBounds;
        if (lb != -INFINITY and !std.math.isFinite(lb)) return error.NonFiniteCoefficient;
        if (ub != INFINITY and !std.math.isFinite(ub)) return error.NonFiniteCoefficient;
    }

    // ── Row bounds ─────────────────────────────────────────────────────────
    if (model.row_lower.len != model.num_rows) return error.DimensionMismatch;
    if (model.row_upper.len != model.num_rows) return error.DimensionMismatch;
    for (0..model.num_rows) |row| {
        const lb = model.row_lower[row];
        const ub = model.row_upper[row];
        if (lb > ub) return error.InvalidBounds;
        if (lb != -INFINITY and !std.math.isFinite(lb)) return error.NonFiniteCoefficient;
        if (ub != INFINITY and !std.math.isFinite(ub)) return error.NonFiniteCoefficient;
    }

    // ── Integrality ────────────────────────────────────────────────────────
    if (model.integrality) |int| {
        if (int.len != model.num_cols) return error.InvalidIntegrality;
    }
}

// ── Hessian validation ─────────────────────────────────────────────────────

/// Validates a Hessian matrix.
pub fn validateHessian(hessian: *const Hessian) ValidationError!void {
    hessian.validate() catch |err| return switch (err) {
        error.DimensionMismatch => error.DimensionMismatch,
        error.InvalidHessian => error.InvalidHessian,
        error.NonFiniteValue => error.NonFiniteValue,
        error.IndicesNotIncreasing => error.IndicesNotIncreasing,
    };
}

// ── QuadraticConstraint validation ─────────────────────────────────────────

/// Validates a single quadratic constraint.
pub fn validateQuadraticConstraint(qc: *const QuadraticConstraint, num_cols: usize) ValidationError!void {
    qc.validate(num_cols) catch |err| return switch (err) {
        error.DimensionMismatch => error.DimensionMismatch,
        error.InvalidBounds => error.InvalidBounds,
        error.InvalidHessian => error.InvalidHessian,
        error.NonFiniteValue => error.NonFiniteValue,
        error.ColumnOutOfRange => error.ColumnOutOfRange,
        error.IndicesNotIncreasing => error.IndicesNotIncreasing,
    };
}

// ── QuadraticModel validation ──────────────────────────────────────────────

/// Validates a quadratic model, delegating to sub-validators.
pub fn validateQuadraticModel(model: *const QuadraticModel) ValidationError!void {
    try validateLinearModel(&model.linear);
    if (model.objective_hessian) |*h| try validateHessian(h);
    for (model.quadratic_constraints) |*qc| {
        try validateQuadraticConstraint(qc, model.linear.num_cols);
    }
}

// ── ExpressionGraph validation ─────────────────────────────────────────────

/// Validates an expression graph.
pub fn validateExpressionGraph(graph: *const ExpressionGraph, num_cols: ?usize) (ValidationError || std.mem.Allocator.Error)!void {
    try graph.validate(num_cols);
}

// ── NonlinearModel validation ──────────────────────────────────────────────

/// Validates a nonlinear model.
pub fn validateNonlinearModel(model: *const NonlinearModel) (ValidationError || std.mem.Allocator.Error)!void {
    try validateQuadraticModel(&model.quadratic);
    try validateExpressionGraph(&model.graph, model.quadratic.linear.num_cols);

    // Check nonlinear roots belong to the graph.
    if (model.nonlinear_objective) |root| {
        const idx = @intFromEnum(root);
        if (idx >= model.graph.num_nodes) return error.InvalidExpressionRoot;
    }
    for (model.nonlinear_constraints) |nc| {
        const idx = @intFromEnum(nc.root);
        if (idx >= model.graph.num_nodes) return error.InvalidExpressionRoot;
    }
    if (!model.hasNonlinearRoot()) return error.EmptyNonlinearModel;
}

// ── CompiledModel validation ───────────────────────────────────────────────

/// Validates a compiled model by dispatching to the active variant's validator.
pub fn validateCompiledModel(model: *const CompiledModel) (ValidationError || std.mem.Allocator.Error)!void {
    switch (model.*) {
        .linear => |*m| try validateLinearModel(m),
        .quadratic => |*m| try validateQuadraticModel(m),
        .nonlinear => |*m| try validateNonlinearModel(m),
    }
}

// ── Tests ──────────────────────────────────────────────────────────────────

test "validateLinearModel on empty model" {
    var lm = try LinearModel.initEmpty(std.testing.allocator);
    defer lm.deinit();
    try validateLinearModel(&lm);
}

test "validateLinearModel on typical LP" {
    const Builder = @import("linear_model_builder.zig").LinearModelBuilder;
    var builder = try Builder.init(std.testing.allocator, 2, 3);
    defer builder.deinit();
    try builder.appendCoefficient(0, 0, 1.0);
    try builder.appendCoefficient(0, 1, 1.0);
    try builder.appendCoefficient(1, 1, -1.0);
    try builder.appendCoefficient(1, 2, 2.0);
    var lm = try builder.freeze();
    defer lm.deinit();
    try validateLinearModel(&lm);
}

test "validateLinearModel rejects invalid bounds" {
    const Builder = @import("linear_model_builder.zig").LinearModelBuilder;
    var builder = try Builder.init(std.testing.allocator, 1, 1);
    defer builder.deinit();
    try builder.setColBounds(0, 0.0, 1.0);

    var lm = try builder.freeze();
    defer lm.deinit();
    lm.col_lower[0] = 5.0;
    lm.col_upper[0] = 3.0;
    try std.testing.expectError(error.InvalidBounds, validateLinearModel(&lm));
}

test "validateLinearModel rejects non-finite cost" {
    const Builder = @import("linear_model_builder.zig").LinearModelBuilder;
    var builder = try Builder.init(std.testing.allocator, 1, 1);
    defer builder.deinit();
    var lm = try builder.freeze();
    defer lm.deinit();
    lm.col_cost[0] = std.math.nan(f64);
    try std.testing.expectError(error.NonFiniteObjective, validateLinearModel(&lm));
}

test "validateHessian on empty Hessian" {
    var h = try Hessian.initEmpty(std.testing.allocator);
    defer h.deinit();
    try validateHessian(&h);
}

test "validateQuadraticModel delegates to linear" {
    const lm = try LinearModel.initEmpty(std.testing.allocator);
    const qm = QuadraticModel{
        .allocator = std.testing.allocator,
        .linear = lm,
        .objective_hessian = null,
        .quadratic_constraints = &.{},
    };
    var cm = CompiledModel{ .quadratic = qm };
    defer cm.deinit();
    try validateCompiledModel(&cm);
}

test "validateCompiledModel nonlinear variant" {
    var bld = expression_graph_module.ExpressionGraphBuilder.init(std.testing.allocator);
    defer bld.deinit();
    const c1 = try bld.addConstant(1.0);
    const graph = try bld.freeze();
    const lm = try LinearModel.initEmpty(std.testing.allocator);
    const nm = NonlinearModel{
        .allocator = std.testing.allocator,
        .quadratic = QuadraticModel{
            .allocator = std.testing.allocator,
            .linear = lm,
            .objective_hessian = null,
            .quadratic_constraints = &.{},
        },
        .graph = graph,
        .nonlinear_objective = c1,
        .nonlinear_constraints = &.{},
    };
    var cm = CompiledModel{ .nonlinear = nm };
    defer cm.deinit();
    try validateCompiledModel(&cm);
}
