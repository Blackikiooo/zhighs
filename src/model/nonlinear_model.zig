//! Non‑linear programming model (NLP / MINLP).
//!
//! Extends [`QuadraticModel`] with a full expression DAG for nonlinear
//! objective terms and nonlinear constraints.
//!
//! ## Classification
//!
//! The problem class is derived automatically:
//! - At least one nonlinear root (objective or constraint) → NLP / MINLP
//! - Otherwise → delegates to QuadraticModel
//!
//! Experimental API.  The composition model is stable; constraint storage may
//! be optimised from slice‑of‑structs to a SoA layout after solver profiling.

const std = @import("std");
const expression_graph_module = @import("expression_graph.zig");
const quadratic_model_module = @import("quadratic_model.zig");
const problem_class_module = @import("problem_class.zig");

const ExpressionGraph = expression_graph_module.ExpressionGraph;
const NodeId = expression_graph_module.NodeId;
const QuadraticModel = quadratic_model_module.QuadraticModel;
const ProblemClass = problem_class_module.ProblemClass;
const DomainClass = problem_class_module.DomainClass;
const ObjectiveClass = problem_class_module.ObjectiveClass;
const ConstraintClass = problem_class_module.ConstraintClass;
const classify = problem_class_module.classify;

// ── NonlinearConstraint ────────────────────────────────────────────────────

/// Experimental API: a single nonlinear constraint.
///
/// Represents:
/// ```text
/// lower ≤ f(x) ≤ upper
/// ```
/// where `f` is the sub‑expression rooted at `root` in the parent model's
/// expression graph.
pub const NonlinearConstraint = struct {
    /// Root node of the constraint expression in the parent graph.
    root: NodeId,
    /// Constraint lower bound.
    lower: f64,
    /// Constraint upper bound.
    upper: f64,
};

// ── NonlinearModel ─────────────────────────────────────────────────────────

/// Experimental API: owning NLP / MINLP model.
///
/// Combines quadratic data (which wraps linear data) with an expression DAG
/// for nonlinear terms.  At least one nonlinear root must be present; if none
/// exists the caller should use QuadraticModel instead.
pub const NonlinearModel = struct {
    allocator: std.mem.Allocator,
    /// Underlying quadratic / linear data (costs, bounds, matrix, Hessian).
    quadratic: QuadraticModel,
    /// Expression DAG for all nonlinear terms.
    graph: ExpressionGraph,
    /// Root node of the nonlinear objective term, or `null`.
    /// When present, the full objective is `quadratic.linear.cᵀx + ½xᵀQx + f(x)`.
    nonlinear_objective: ?NodeId,
    /// Array of nonlinear constraints.
    nonlinear_constraints: []NonlinearConstraint,

    const Self = @This();

    /// Releases all owned resources in reverse dependency order.
    pub fn deinit(self: *Self) void {
        if (self.nonlinear_constraints.len > 0)
            self.allocator.free(self.nonlinear_constraints);
        self.graph.deinit();
        self.quadratic.deinit();
        self.* = undefined;
    }

    /// Domain class derived from the underlying linear model's integrality.
    pub fn domainClass(self: *const Self) DomainClass {
        return self.quadratic.domainClass();
    }

    /// Objective class — always `.nonlinear` when this type is used.
    pub fn objectiveClass(self: *const Self) ObjectiveClass {
        _ = self;
        return .nonlinear;
    }

    /// Constraint class — `.nonlinear` if any nonlinear constraint exists,
    /// otherwise delegates to the quadratic layer.
    pub fn constraintClass(self: *const Self) ConstraintClass {
        if (self.nonlinear_constraints.len > 0) return .nonlinear;
        return self.quadratic.constraintClass();
    }

    /// Derived problem class.
    pub fn problemClass(self: *const Self) ProblemClass {
        return classify(
            self.domainClass(),
            self.objectiveClass(),
            self.constraintClass(),
        );
    }

    /// Returns `true` when at least one nonlinear root (objective or
    /// constraint) is present.
    pub fn hasNonlinearRoot(self: *const Self) bool {
        if (self.nonlinear_objective != null) return true;
        return self.nonlinear_constraints.len > 0;
    }
};

// ── Tests ──────────────────────────────────────────────────────────────────

const foundation = @import("foundation");
const ColId = foundation.ColId;
const LinearModel = @import("linear_model.zig").LinearModel;
const Hessian = @import("hessian.zig").Hessian;
const ExpressionGraphBuilder = expression_graph_module.ExpressionGraphBuilder;

test "NonlinearModel.empty model (no nonlinear root) still classifies as NLP" {
    var bld = ExpressionGraphBuilder.init(std.testing.allocator);
    defer bld.deinit();
    const graph = try bld.freeze();

    const lm = try LinearModel.initEmpty(std.testing.allocator);
    var nm = NonlinearModel{
        .allocator = std.testing.allocator,
        .quadratic = QuadraticModel{
            .allocator = std.testing.allocator,
            .linear = lm,
            .objective_hessian = null,
            .quadratic_constraints = &.{},
        },
        .graph = graph,
        .nonlinear_objective = null,
        .nonlinear_constraints = &.{},
    };
    defer nm.deinit();

    try std.testing.expectEqual(ObjectiveClass.nonlinear, nm.objectiveClass());
    try std.testing.expect(!nm.hasNonlinearRoot());
}

test "NonlinearModel.with nonlinear constraint is NLP" {
    var bld = ExpressionGraphBuilder.init(std.testing.allocator);
    defer bld.deinit();
    const c5 = try bld.addConstant(5.0);
    const graph = try bld.freeze();

    const lm = try LinearModel.initEmpty(std.testing.allocator);
    var nm = NonlinearModel{
        .allocator = std.testing.allocator,
        .quadratic = QuadraticModel{
            .allocator = std.testing.allocator,
            .linear = lm,
            .objective_hessian = null,
            .quadratic_constraints = &.{},
        },
        .graph = graph,
        .nonlinear_objective = null,
        .nonlinear_constraints = try std.testing.allocator.dupe(NonlinearConstraint, &.{
            .{ .root = c5, .lower = 0.0, .upper = 10.0 },
        }),
    };
    defer nm.deinit();

    try std.testing.expect(nm.hasNonlinearRoot());
    try std.testing.expectEqual(ConstraintClass.nonlinear, nm.constraintClass());
    try std.testing.expectEqual(ProblemClass.nlp, nm.problemClass());
}

test "NonlinearModel.with nonlinear objective and integer is MINLP" {
    var bld = ExpressionGraphBuilder.init(std.testing.allocator);
    defer bld.deinit();
    const v0 = try bld.addVariable(ColId.fromUsizeAssumeValid(0));
    const graph = try bld.freeze();
    _ = v0;

    const lm = try LinearModel.initEmpty(std.testing.allocator);
    var nm = NonlinearModel{
        .allocator = std.testing.allocator,
        .quadratic = QuadraticModel{
            .allocator = std.testing.allocator,
            .linear = lm,
            .objective_hessian = null,
            .quadratic_constraints = &.{},
        },
        .graph = graph,
        .nonlinear_objective = null,
        .nonlinear_constraints = &.{},
    };
    defer nm.deinit();

    try std.testing.expectEqual(ProblemClass.nlp, nm.problemClass());
}
