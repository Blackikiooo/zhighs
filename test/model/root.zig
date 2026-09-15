//! Model core test entry point, independent of API, solver, and presolve modules.
//!
//! Run with:
//!   zig build test-model
//!   zig build test-model -Dhighs-int-width=w64

const std = @import("std");
const model = @import("model");

test "semi-variable zero branch survives the model solve pipeline" {
    for ([_]model.VarType{ .semicont, .semiint }) |kind| {
        for ([_]i64{ -1, 0 }) |presolve_setting| {
            var env = try model.Env.initSimple(std.testing.allocator);
            defer env.deinit();
            try env.setIntParam("Presolve", presolve_setting);
            var instance = try model.Model.init(std.testing.allocator, &env, "semi-zero");
            defer instance.deinit();
            // x is either zero or in [2, 4]. The row forces the zero branch.
            try instance.addVar(0, &.{}, &.{}, 1.0, 2.0, 4.0, kind, "x");
            try instance.addConstr(1, &.{0}, &.{1.0}, .less_equal, 1.0, "zero-only");
            try instance.optimize();
            try std.testing.expectEqual(model.Status.optimal, instance.status);
            try std.testing.expectApproxEqAbs(@as(f64, 0.0), instance.obj_val, 1e-8);
            try std.testing.expectApproxEqAbs(@as(f64, 0.0), instance.solution[0], 1e-8);
            try std.testing.expectApproxEqAbs(@as(f64, 1.0), instance.slack[0], 1e-8);
        }
    }
}

test "binary bounds survive empty-column presolve and postsolve" {
    for ([_]i64{ -1, 0 }) |presolve_setting| {
        var env = try model.Env.initSimple(std.testing.allocator);
        defer env.deinit();
        try env.setIntParam("Presolve", presolve_setting);
        var instance = try model.Model.init(std.testing.allocator, &env, "binary-domain");
        defer instance.deinit();
        try instance.addVar(0, &.{}, &.{}, -1.0, 0.0, 5.0, .binary, "x");
        try instance.optimize();
        try std.testing.expectEqual(model.Status.optimal, instance.status);
        try std.testing.expectApproxEqAbs(@as(f64, -1.0), instance.obj_val, 1e-8);
        try std.testing.expectApproxEqAbs(@as(f64, 1.0), instance.solution[0], 1e-8);
    }
}

test "MIP limit without incumbent publishes an unknown bound in either sense" {
    for ([_]i64{ 1, -1 }) |sense| {
        var env = try model.Env.initSimple(std.testing.allocator);
        defer env.deinit();
        try env.setIntParam("NodeLimit", 0);
        var instance = try model.Model.init(std.testing.allocator, &env, "no-incumbent");
        defer instance.deinit();
        try instance.addVar(0, &.{}, &.{}, 1.0, 0.0, 1.0, .binary, "x");
        try instance.addVar(0, &.{}, &.{}, 2.0, 0.0, 1.0, .binary, "y");
        try instance.addConstr(2, &.{ 0, 1 }, &.{ 2.0, 2.0 }, .less_equal, 3.0, "capacity");
        try instance.setIntAttr(.model_sense, sense);
        try instance.optimize();
        try std.testing.expectEqual(model.Status.node_limit, instance.status);
        try std.testing.expectEqual(@as(i64, 0), instance.sol_count);
        const bound = try instance.getDblAttr(.obj_bound);
        try std.testing.expect(if (sense == 1) std.math.isNegativeInf(bound) else std.math.isPositiveInf(bound));
    }
}

test {
    // Force compilation of all exported declarations (and their inline tests).
    std.testing.refAllDecls(model);
}

test "pending added row coefficients enter CSC and overlap added columns deterministically" {
    var env = try model.Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var instance = try model.Model.init(std.testing.allocator, &env, "structural_plan_rows");
    defer instance.deinit();

    // Both sides describe the same new coordinate; the later row payload is
    // the final set value and must not create a duplicate CSC entry.
    try instance.addConstr(1, &.{0}, &.{2.0}, .less_equal, 5.0, null);
    try instance.addVar(1, &.{0}, &.{1.0}, 0.0, 0.0, 10.0, .continuous, null);
    try instance.updateModel();
    try std.testing.expectEqual(@as(f64, 2.0), try instance.getCoeff(0, 0));
    try std.testing.expectEqual(@as(usize, 1), instance.matrix.csc().nnz());
    try instance.matrix.csc().validate();

    try instance.addConstr(1, &.{0}, &.{-3.0}, .greater_equal, -1.0, null);
    try instance.updateModel();
    try std.testing.expectEqual(@as(f64, -3.0), try instance.getCoeff(1, 0));
    try instance.matrix.csc().validate();
}

test "new object scalar edits fold into structural initialization" {
    var env = try model.Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var instance = try model.Model.init(std.testing.allocator, &env, "structural_plan_scalars");
    defer instance.deinit();

    try instance.addConstr(0, &.{}, &.{}, .less_equal, 1.0, "planned-row");
    try instance.addVar(1, &.{0}, &.{3.0}, 1.0, 0.0, 10.0, .continuous, "planned-column");
    try instance.chgBounds(0, -2.0, 6.0);
    try instance.chgObj(0, 9.0);
    try instance.chgVarType(0, .integer);
    try instance.chgRHS(0, 12.0);
    try instance.chgSense(0, .equal);
    try instance.updateModel();

    try std.testing.expectEqual(@as(f64, -2.0), instance.var_lb[0]);
    try std.testing.expectEqual(@as(f64, 6.0), instance.var_ub[0]);
    try std.testing.expectEqual(@as(f64, 9.0), instance.var_obj[0]);
    try std.testing.expectEqual(model.VarType.integer, instance.var_type[0]);
    try std.testing.expectEqual(@as(f64, 12.0), instance.constr_rhs[0]);
    try std.testing.expectEqual(model.Sense.equal, instance.constr_sense[0]);
    try std.testing.expectEqual(@as(usize, 0), try instance.getVarByName("planned-column"));
    try std.testing.expectEqual(@as(usize, 0), try instance.getConstrByName("planned-row"));
    try std.testing.expectEqual(@as(f64, 3.0), try instance.getCoeff(0, 0));
    try instance.matrix.csc().validate();
}

test "MPS file with integer variables solves to the known MIP optimum" {
    // test/model/data/tiny_mip.mps:
    //   maximize 3 X + 2 Y  s.t. 2 X + 3 Y <= 4, X, Y binary.
    // The LP relaxation is fractional (X=1, Y=2/3), so branch-and-bound must
    // run; the integer optimum is X=1, Y=0 with objective 3.
    var env = try model.Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var instance = try model.Model.init(std.testing.allocator, &env, "tiny_mip");
    defer instance.deinit();

    try instance.read("test/model/data/tiny_mip.mps");
    try std.testing.expectEqual(@as(i64, 1), try instance.getIntAttr(.is_mip));
    try instance.optimize();

    try std.testing.expectEqual(model.Status.optimal, instance.status);
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), try instance.getDblAttr(.obj_val), 1e-6);
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), try instance.getDblAttr(.obj_bound), 1e-6);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), try instance.getDblAttrElement(.x, 0), 1e-6);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), try instance.getDblAttrElement(.x, 1), 1e-6);
}

test "MIP presolve solves a reducible model and postsolve restores original coordinates" {
    // test/model/data/presolve_mip.mps:
    //   min X + Y + 5 Z  s.t. X + Y + Z >= 2, Z = 1, X, Y binary.
    // Presolve substitutes the fixed Z (=1) and removes the now-empty fixing
    // row, reducing the model to 2 columns and 1 row; the postsolved solution
    // must satisfy the original constraints with objective 6.
    var env = try model.Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var instance = try model.Model.init(std.testing.allocator, &env, "presolve_mip");
    defer instance.deinit();

    try instance.read("test/model/data/presolve_mip.mps");
    try std.testing.expectEqual(@as(i64, 1), try instance.getIntAttr(.is_mip));
    try instance.optimize();

    try std.testing.expectEqual(model.Status.optimal, instance.status);
    try std.testing.expectApproxEqAbs(@as(f64, 6.0), try instance.getDblAttr(.obj_val), 1e-6);
    try std.testing.expectApproxEqAbs(@as(f64, 6.0), try instance.getDblAttr(.obj_bound), 1e-6);
    // Original-coordinate solution: Z restored to its fixed value 1 and the
    // demand row X + Y + Z >= 2 is satisfied.
    const x = try instance.getDblAttrElement(.x, 0);
    const y = try instance.getDblAttrElement(.x, 1);
    const z = try instance.getDblAttrElement(.x, 2);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), z, 1e-6);
    try std.testing.expect(x + y + z >= 2.0 - 1e-6);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), @round(x), 1e-6);
    try std.testing.expectApproxEqAbs(x, @round(x), 1e-6);
    try std.testing.expectApproxEqAbs(y, @round(y), 1e-6);
}

test "MIP presolve=off matches the presolved optimum" {
    var env = try model.Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var instance = try model.Model.init(std.testing.allocator, &env, "presolve_mip_off");
    defer instance.deinit();

    try instance.read("test/model/data/presolve_mip.mps");
    try instance.setIntParam("Presolve", 0);
    try instance.optimize();

    try std.testing.expectEqual(model.Status.optimal, instance.status);
    try std.testing.expectApproxEqAbs(@as(f64, 6.0), try instance.getDblAttr(.obj_val), 1e-6);
    try std.testing.expectApproxEqAbs(@as(f64, 6.0), try instance.getDblAttr(.obj_bound), 1e-6);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), try instance.getDblAttrElement(.x, 2), 1e-6);
}

test "MIP presolve keeps a fractionally fixed integer column infeasible" {
    // An integer variable fixed at 0.5 is infeasible. MIP-mode presolve must
    // not substitute it away (that would silently drop the infeasibility).
    var env = try model.Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var instance = try model.Model.init(std.testing.allocator, &env, "fractional_fixed");
    defer instance.deinit();

    try instance.addVar(0, &.{}, &.{}, 1.0, 0.5, 0.5, .integer, "half");
    try instance.optimize();
    try std.testing.expectEqual(model.Status.infeasible, instance.status);
}

test "strong-branching probe interrupted mid Phase I cannot poison the warm basis" {
    // test/model/data/warm_restart_mip.lp: randomly generated 22 binary
    // variables / 14 covering constraints (seed 7). A strong-branching probe
    // hits its per-probe iteration limit inside primal Phase I and returns
    // with artificial columns still basic; the next warm reoptimization must
    // cold-start that basis instead of indexing original bounds through the
    // artificial column index (engine_dual.zig index-out-of-bounds panic).
    var env = try model.Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var instance = try model.Model.init(std.testing.allocator, &env, "warm_restart_mip");
    defer instance.deinit();

    try instance.read("test/model/data/warm_restart_mip.lp");
    try std.testing.expectEqual(@as(i64, 1), try instance.getIntAttr(.is_mip));
    try instance.optimize();

    try std.testing.expectEqual(model.Status.optimal, instance.status);
    try std.testing.expectApproxEqAbs(@as(f64, 63.0), try instance.getDblAttr(.obj_val), 1e-6);
    try std.testing.expectApproxEqAbs(@as(f64, 63.0), try instance.getDblAttr(.obj_bound), 1e-6);
}
