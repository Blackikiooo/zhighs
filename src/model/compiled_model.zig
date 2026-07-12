//! Compiled model — solver‑internal snapshot of a fully validated model.
//!
//! `CompiledModel` is a tagged union that holds exactly one of the three
//! solver‑internal IR types.  Its purpose is to serve as the immutable,
//! fully‑validated input to presolve and solver dispatchers.
//!
//! ## Move‑only semantics
//!
//! `CompiledModel` owns its data.  Copying the union value produces two
//! handles to the same allocation; the second `deinit` would be a double
//! free.  In Zig this is prevented by the language (`union(enum)` cannot
//! be copied implicitly when it holds array fields), but the caller must
//! still avoid `std.mem.copy` or explicit field‑by‑field copies.
//!
//! Experimental API.  The union layout is stable; `clone` may be added
//! once the Model‑to‑Compiled‑Model compilation path is implemented.

const std = @import("std");
const problem_class_module = @import("problem_class.zig");
const linear_model_module = @import("linear_model.zig");
const quadratic_model_module = @import("quadratic_model.zig");
const nonlinear_model_module = @import("nonlinear_model.zig");

const ProblemClass = problem_class_module.ProblemClass;
const DomainClass = problem_class_module.DomainClass;
const ObjectiveClass = problem_class_module.ObjectiveClass;
const ConstraintClass = problem_class_module.ConstraintClass;
const LinearModel = linear_model_module.LinearModel;
const QuadraticModel = quadratic_model_module.QuadraticModel;
const NonlinearModel = nonlinear_model_module.NonlinearModel;

/// Experimental API: owning solver‑internal model snapshot.
pub const CompiledModel = union(enum) {
    linear: LinearModel,
    quadratic: QuadraticModel,
    nonlinear: NonlinearModel,

    const Self = @This();

    /// Returns the problem class derived from the active variant.
    pub fn problemClass(self: *const Self) ProblemClass {
        return switch (self.*) {
            inline .linear => |*m| m.problemClass(),
            inline .quadratic => |*m| m.problemClass(),
            inline .nonlinear => |*m| m.problemClass(),
        };
    }

    /// Returns the domain class (continuous or mixed‑integer).
    pub fn domainClass(self: *const Self) DomainClass {
        return switch (self.*) {
            inline .linear => |*m| m.domainClass(),
            inline .quadratic => |*m| m.domainClass(),
            inline .nonlinear => |*m| m.domainClass(),
        };
    }

    /// Returns a pointer to the underlying linear data.
    ///
    /// Every variant embeds a `LinearModel` (directly or transitively),
    /// so this accessor always succeeds.
    pub fn linearData(self: *const Self) *const LinearModel {
        return switch (self.*) {
            .linear => |*m| m,
            .quadratic => |*m| &m.linear,
            .nonlinear => |*m| &m.quadratic.linear,
        };
    }

    /// Releases all owned resources for the active variant.
    pub fn deinit(self: *Self) void {
        switch (self.*) {
            .linear => |*m| m.deinit(),
            .quadratic => |*m| m.deinit(),
            .nonlinear => |*m| m.deinit(),
        }
        self.* = undefined;
    }
};

// ── Tests ──────────────────────────────────────────────────────────────────

test "CompiledModel.linear classification" {
    const lm = try LinearModel.initEmpty(std.testing.allocator);
    var cm = CompiledModel{ .linear = lm };
    defer cm.deinit();

    try std.testing.expectEqual(ProblemClass.lp, cm.problemClass());
    try std.testing.expectEqual(DomainClass.continuous, cm.domainClass());
}

test "CompiledModel.linearData accessor" {
    const lm = try LinearModel.initEmpty(std.testing.allocator);
    var cm = CompiledModel{ .linear = lm };
    defer cm.deinit();

    const data = cm.linearData();
    try std.testing.expectEqual(@as(usize, 0), data.num_rows);
    try std.testing.expectEqual(@as(usize, 0), data.num_cols);
}

test "CompiledModel.quadratic classification" {
    const lm = try LinearModel.initEmpty(std.testing.allocator);
    const qm = QuadraticModel{
        .allocator = std.testing.allocator,
        .linear = lm,
        .objective_hessian = null,
        .quadratic_constraints = &.{},
    };
    var cm = CompiledModel{ .quadratic = qm };
    defer cm.deinit();

    try std.testing.expectEqual(ProblemClass.lp, cm.problemClass());
}

test "CompiledModel.nonlinear classification" {
    var bld = @import("expression_graph.zig").ExpressionGraphBuilder.init(std.testing.allocator);
    defer bld.deinit();
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
        .nonlinear_objective = null,
        .nonlinear_constraints = &.{},
    };
    var cm = CompiledModel{ .nonlinear = nm };
    defer cm.deinit();

    try std.testing.expectEqual(ProblemClass.nlp, cm.problemClass());
}

test "CompiledModel.linearData on quadratic variant" {
    const lm = try LinearModel.initEmpty(std.testing.allocator);
    const qm = QuadraticModel{
        .allocator = std.testing.allocator,
        .linear = lm,
        .objective_hessian = null,
        .quadratic_constraints = &.{},
    };
    var cm = CompiledModel{ .quadratic = qm };
    defer cm.deinit();

    const data = cm.linearData();
    try std.testing.expectEqual(@as(usize, 0), data.num_rows);
}
