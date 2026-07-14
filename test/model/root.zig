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
