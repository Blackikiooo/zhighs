//! Expression and column builder types.
//!
//! Re-exports the sub-modules:
//!
//! - [`Column`](column.zig) — sparse column for incremental model building
//! - [`LinExpr`](lin_expr.zig) — linear expression `c + Σ aᵢ xᵢ`
//! - [`QuadExpr`](quad_expr.zig) — quadratic expression
//! - [`TempConstr`](temp_constr.zig) — constraint builder from expressions

const std = @import("std");

pub const Column = @import("column.zig").Column;
pub const LinExpr = @import("lin_expr.zig").LinExpr;
pub const QuadExpr = @import("quad_expr.zig").QuadExpr;
pub const TempConstr = @import("temp_constr.zig").TempConstr;
pub const leExpr = @import("temp_constr.zig").leExpr;
pub const eqExpr = @import("temp_constr.zig").eqExpr;
pub const geExpr = @import("temp_constr.zig").geExpr;

// ── Tests ───────────────────────────────────────────────────────────────

test "LinExpr stores and evaluates terms" {
    const Env = @import("../env.zig").Env;
    const Model = @import("../model.zig").Model;
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, "c0");
    try model.addVar(1, &[_]usize{0}, &[_]f64{1.0}, 0.0, 0.0, 1.0, .continuous, "x0");
    try model.addVar(1, &[_]usize{0}, &[_]f64{1.0}, 0.0, 0.0, 1.0, .continuous, "x1");
    try model.updateModel();

    var expr = LinExpr.init(std.testing.allocator);
    defer expr.deinit();

    try expr.addTermByIndex(0, 2.0);
    try expr.addTermByIndex(1, 3.0);
    expr.addConstant(1.0);

    try std.testing.expectEqual(@as(usize, 2), expr.numTerms());
    try std.testing.expectEqual(@as(f64, 1.0), expr.getConstant());
}

test "Column stores terms" {
    var col = Column.init();
    defer col.deinit(std.testing.allocator);

    try col.addTermByIndex(std.testing.allocator, 0, 1.5);
    try col.addTermByIndex(std.testing.allocator, 1, -2.0);

    try std.testing.expectEqual(@as(usize, 2), col.numNz());
    try std.testing.expectEqual(@as(usize, 0), col.indices.items[0]);
    try std.testing.expectEqual(@as(f64, 1.5), col.values.items[0]);
    try std.testing.expectEqual(@as(usize, 1), col.indices.items[1]);
    try std.testing.expectEqual(@as(f64, -2.0), col.values.items[1]);
}

test "TempConstr from leExpr has correct sense" {
    const Env = @import("../env.zig").Env;
    const Model = @import("../model.zig").Model;
    const Sense = @import("../types.zig").Sense;
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, "c0");
    try model.addVar(0, &.{}, &.{}, 0.0, 0.0, 1.0, .continuous, null);
    try model.addVar(0, &.{}, &.{}, 0.0, 0.0, 1.0, .continuous, null);
    try model.updateModel();

    var expr = LinExpr.init(std.testing.allocator);
    defer expr.deinit();
    try expr.addTermByIndex(0, 1.0);
    try expr.addTermByIndex(1, 1.0);

    const tc = leExpr(&expr, 10.0);
    try std.testing.expectEqual(Sense.less_equal, tc.sense);
    try std.testing.expectEqual(@as(f64, 10.0), tc.rhs);
}
