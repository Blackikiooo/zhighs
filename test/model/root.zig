//! Model core test entry point, independent of API, solver, and presolve modules.
//!
//! Run with:
//!   zig build test-model
//!   zig build test-model -Dhighs-int-width=w64

const std = @import("std");
const model = @import("model");

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
