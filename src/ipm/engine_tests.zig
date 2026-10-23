//! End-to-end IPM tests against analytic optima and failure paths.

const std = @import("std");
const foundation = @import("foundation");
const matrix = @import("matrix");
const lp = @import("lp");

const ipm = @import("root.zig");

const inf = std.math.inf(f64);
const RowId = foundation.RowId;

const check_tolerance = 1e-5;

fn solveExpectOptimal(problem: lp.ProblemView) !ipm.Result {
    var result = try ipm.solve(std.testing.allocator, problem, .{});
    errdefer result.deinit();
    var dense = try ipm.solve(std.testing.allocator, problem, .{ .dense_normal = true });
    defer dense.deinit();
    var scaled = try ipm.solve(std.testing.allocator, problem, .{ .equilibrate = true });
    defer scaled.deinit();
    try std.testing.expectEqual(dense.status, result.status);
    try std.testing.expectEqual(scaled.status, result.status);
    try std.testing.expectApproxEqAbs(dense.objective_value, result.objective_value, check_tolerance);
    try std.testing.expectApproxEqAbs(scaled.objective_value, result.objective_value, check_tolerance);
    return result;
}

fn oneColumnProblem(lower: f64, upper: f64, cost: f64, rhs: f64, storage: *[5]f64) lp.ProblemView {
    storage.* = .{ lower, upper, cost, rhs, rhs };
    return .{
        .num_rows = 1,
        .num_cols = 1,
        .col_lower = storage[0..1],
        .col_upper = storage[1..2],
        .col_cost = storage[2..3],
        .row_lower = storage[3..4],
        .row_upper = storage[4..5],
        .matrix = matrix.CscView.initAssumeValid(1, 1, &.{ 0, 1 }, &.{RowId.fromUsizeAssumeValid(0)}, &.{1}),
        .objective_sense = .minimize,
        .objective_offset = 4,
    };
}

test "ipm handles all-fixed models and contradictory constant rows" {
    var storage: [5]f64 = undefined;
    var result = try ipm.solve(std.testing.allocator, oneColumnProblem(2, 2, 3, 2, &storage), .{});
    defer result.deinit();
    try std.testing.expectEqual(.optimal, result.status);
    try std.testing.expectApproxEqAbs(10.0, result.objective_value, 1e-12);
    try std.testing.expectEqual(@as(f64, 2), result.primal[0]);
    try std.testing.expectEqual(@as(f64, 0), result.dual[0]);
    var infeasible = try ipm.solve(std.testing.allocator, oneColumnProblem(2, 2, 3, 3, &storage), .{});
    defer infeasible.deinit();
    try std.testing.expectEqual(.infeasible, infeasible.status);
}

test "ipm validates controls and obeys interruption time and work limits" {
    var storage: [5]f64 = undefined;
    const problem = oneColumnProblem(0, inf, 1, 2, &storage);
    const a = std.testing.allocator;
    try std.testing.expectError(error.InvalidOptions, ipm.solve(a, problem, .{ .step_fraction = 1 }));
    try std.testing.expectError(error.InvalidOptions, ipm.solve(a, problem, .{ .tolerances = .{ .regularization_growth = 1 } }));
    try std.testing.expectError(error.InvalidOptions, ipm.solve(a, problem, .{ .start_margin = std.math.nan(f64) }));
    try std.testing.expectError(error.ClockRequired, ipm.solve(a, problem, .{ .time_limit_ns = 0 }));
    var flag = std.atomic.Value(bool).init(true);
    var stopped = try ipm.solve(a, problem, .{ .interrupt_flag = &flag });
    defer stopped.deinit();
    try std.testing.expectEqual(.interrupted, stopped.status);
    var work = try ipm.solve(a, problem, .{ .work_limit = 0 });
    defer work.deinit();
    try std.testing.expectEqual(.work_limit, work.status);
    var timed = try ipm.engine.solveWithStats(a, problem, .{ .time_limit_ns = 0 }, std.testing.io, null);
    defer timed.deinit();
    try std.testing.expectEqual(.time_limit, timed.status);
}

test "profile Cholesky matches dense solves with fill and multiple right hand sides" {
    const a = std.testing.allocator;
    // A has non-contiguous row patterns; slack identity ensures full row rank.
    const problem = lp.ProblemView{
        .num_rows = 5,
        .num_cols = 3,
        .col_cost = &.{ 1, 2, 3 },
        .col_lower = &.{ 0, 0, 0 },
        .col_upper = &.{ inf, inf, inf },
        .row_lower = &.{ -inf, -inf, -inf, -inf, -inf },
        .row_upper = &.{ 1, 1, 1, 1, 1 },
        .matrix = matrix.CscView.initAssumeValid(5, 3, &.{ 0, 2, 4, 7 }, &.{
            RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(3),
            RowId.fromUsizeAssumeValid(1), RowId.fromUsizeAssumeValid(3),
            RowId.fromUsizeAssumeValid(1), RowId.fromUsizeAssumeValid(2),
            RowId.fromUsizeAssumeValid(4),
        }, &.{ 2, -1, 3, 1, -2, 4, 1 }),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var sf = try ipm.standard_form.build(a, problem);
    defer sf.deinit();
    const graph_mod = @import("graph.zig");
    var g = try graph_mod.build(a, sf.num_rows, sf.col_starts, sf.row_indices);
    defer g.deinit();
    var order = try graph_mod.Order.init(a, sf.num_rows, null);
    defer order.deinit();
    const Factor = @import("normal_factor.zig").NormalFactor;
    var profile = try Factor.init(a, &sf, false, &g, &order);
    defer profile.deinit();
    var dense = try Factor.init(a, &sf, true, null, null);
    defer dense.deinit();
    const profile_values = try a.alloc(f64, profile.storageLen());
    defer a.free(profile_values);
    var full: [25]f64 = @splat(0);
    @memset(profile_values, 0);
    for (0..sf.num_cols) |j| {
        for (sf.col_starts[j]..sf.col_starts[j + 1]) |p| {
            for (sf.col_starts[j]..sf.col_starts[j + 1]) |q| {
                const i = sf.row_indices[p];
                const k = sf.row_indices[q];
                full[i * 5 + k] += sf.values[p] * sf.values[q];
                if (k <= i) profile_values[profile.index(i, k)] += sf.values[p] * sf.values[q];
            }
        }
    }
    try profile.factor(profile_values, 0.0);
    try dense.factor(&full, 0.0);
    for (0..3) |repeat| {
        var rhs: [5]f64 = undefined;
        for (&rhs, 0..) |*v, i| v.* = @floatFromInt(i + repeat);
        const original = rhs;
        var reference = rhs;
        profile.solve(&rhs);
        dense.solve(&reference);
        for (rhs, reference) |x, y| try std.testing.expectApproxEqAbs(x, y, 1e-12);
        for (0..5) |i| {
            var value: f64 = 0;
            for (rhs, 0..) |x, j| value += full[i * 5 + j] * x;
            try std.testing.expectApproxEqAbs(original[i], value, 1e-11);
        }
    }
}

test "ipm certifies primal infeasibility from a stalled residual" {
    // min x1 + x2 s.t. x1 + x2 = 3 with x1, x2 ∈ [0,1]: the largest
    // attainable sum is 2, so the primal residual stalls while the dual
    // stays feasible.
    const problem = lp.ProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 1, 1 },
        .col_lower = &.{ 0, 0 },
        .col_upper = &.{ 1, 1 },
        .row_lower = &.{3},
        .row_upper = &.{3},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &.{ 0, 1, 2 }, &.{
            RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(0),
        }, &.{ 1, 1 }),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var result = try ipm.solve(std.testing.allocator, problem, .{});
    defer result.deinit();
    try std.testing.expectEqual(.infeasible, result.status);
}

test "ipm certifies unboundedness along an escaping ray" {
    // min -x1 s.t. x1 - x2 = 0, x >= 0: the ray x1 = x2 -> infinity keeps
    // the primal feasible while the objective plunges.
    const problem = lp.ProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ -1, 0 },
        .col_lower = &.{ 0, 0 },
        .col_upper = &.{ inf, inf },
        .row_lower = &.{0},
        .row_upper = &.{0},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &.{ 0, 1, 2 }, &.{
            RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(0),
        }, &.{ 1, -1 }),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var result = try ipm.solve(std.testing.allocator, problem, .{});
    defer result.deinit();
    try std.testing.expectEqual(.unbounded, result.status);
}

fn allocationProbe(a: std.mem.Allocator) !void {
    var storage: [5]f64 = undefined;
    var result = try ipm.solve(a, oneColumnProblem(0, 5, 2, 1, &storage), .{});
    defer result.deinit();
    try std.testing.expectEqual(.optimal, result.status);
}

test "ipm releases every partial allocation on failure" {
    try std.testing.checkAllAllocationFailures(std.testing.allocator, allocationProbe, .{});
}

test "ipm solves a bounded max LP to its vertex optimum" {
    // max x + y subject to x + 2y ≤ 4, 4x + 2y ≤ 12, x,y ≥ 0.
    // Optimum 10/3 at (8/3, 2/3).
    const col_cost = [_]f64{ 1.0, 1.0 };
    const col_lower = [_]f64{ 0.0, 0.0 };
    const col_upper = [_]f64{ inf, inf };
    const row_lower = [_]f64{ -inf, -inf };
    const row_upper = [_]f64{ 4.0, 12.0 };
    const col_starts = [_]usize{ 0, 2, 4 };
    const row_indices = [_]RowId{ RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(1), RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(1) };
    const values = [_]f64{ 1.0, 4.0, 2.0, 2.0 };
    const view = lp.ProblemView{
        .num_rows = 2,
        .num_cols = 2,
        .col_cost = &col_cost,
        .col_lower = &col_lower,
        .col_upper = &col_upper,
        .row_lower = &row_lower,
        .row_upper = &row_upper,
        .matrix = matrix.CscView.initAssumeValid(2, 2, &col_starts, &row_indices, &values),
        .objective_sense = .maximize,
        .objective_offset = 0.0,
    };
    var result = try solveExpectOptimal(view);
    defer result.deinit();
    try std.testing.expectEqual(.optimal, result.status);
    try std.testing.expectApproxEqAbs(10.0 / 3.0, result.objective_value, check_tolerance);
    try std.testing.expectApproxEqAbs(8.0 / 3.0, result.primal[0], check_tolerance);
    try std.testing.expectApproxEqAbs(2.0 / 3.0, result.primal[1], check_tolerance);
}

test "ipm solves an equality-constrained min LP" {
    // min 2x + y subject to x + y = 3, x,y ≥ 0. Optimum 3 at (0, 3).
    const col_cost = [_]f64{ 2.0, 1.0 };
    const col_lower = [_]f64{ 0.0, 0.0 };
    const col_upper = [_]f64{ inf, inf };
    const row_lower = [_]f64{3.0};
    const row_upper = [_]f64{3.0};
    const col_starts = [_]usize{ 0, 1, 2 };
    const row_indices = [_]RowId{ RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(0) };
    const values = [_]f64{ 1.0, 1.0 };
    const view = lp.ProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &col_cost,
        .col_lower = &col_lower,
        .col_upper = &col_upper,
        .row_lower = &row_lower,
        .row_upper = &row_upper,
        .matrix = matrix.CscView.initAssumeValid(1, 2, &col_starts, &row_indices, &values),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var result = try solveExpectOptimal(view);
    defer result.deinit();
    try std.testing.expectEqual(.optimal, result.status);
    try std.testing.expectApproxEqAbs(3.0, result.objective_value, check_tolerance);
    try std.testing.expectApproxEqAbs(0.0, result.primal[0], check_tolerance);
    try std.testing.expectApproxEqAbs(3.0, result.primal[1], check_tolerance);
}

test "ipm handles column upper bounds through bound rows" {
    // min x subject to x ≥ 1, 0 ≤ x ≤ 5. Optimum 1 at x = 1.
    const col_cost = [_]f64{1.0};
    const col_lower = [_]f64{0.0};
    const col_upper = [_]f64{5.0};
    const row_lower = [_]f64{1.0};
    const row_upper = [_]f64{inf};
    const col_starts = [_]usize{ 0, 1 };
    const row_indices = [_]RowId{RowId.fromUsizeAssumeValid(0)};
    const values = [_]f64{1.0};
    const view = lp.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &col_cost,
        .col_lower = &col_lower,
        .col_upper = &col_upper,
        .row_lower = &row_lower,
        .row_upper = &row_upper,
        .matrix = matrix.CscView.initAssumeValid(1, 1, &col_starts, &row_indices, &values),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var result = try solveExpectOptimal(view);
    defer result.deinit();
    try std.testing.expectEqual(.optimal, result.status);
    try std.testing.expectApproxEqAbs(1.0, result.objective_value, check_tolerance);
    try std.testing.expectApproxEqAbs(1.0, result.primal[0], check_tolerance);
}

test "ipm handles range rows" {
    // min x + 2y subject to 1 ≤ x + y ≤ 2, x,y ≥ 0. Optimum 1 at (1, 0).
    const col_cost = [_]f64{ 1.0, 2.0 };
    const col_lower = [_]f64{ 0.0, 0.0 };
    const col_upper = [_]f64{ inf, inf };
    const row_lower = [_]f64{1.0};
    const row_upper = [_]f64{2.0};
    const col_starts = [_]usize{ 0, 1, 2 };
    const row_indices = [_]RowId{ RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(0) };
    const values = [_]f64{ 1.0, 1.0 };
    const view = lp.ProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &col_cost,
        .col_lower = &col_lower,
        .col_upper = &col_upper,
        .row_lower = &row_lower,
        .row_upper = &row_upper,
        .matrix = matrix.CscView.initAssumeValid(1, 2, &col_starts, &row_indices, &values),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var result = try solveExpectOptimal(view);
    defer result.deinit();
    try std.testing.expectEqual(.optimal, result.status);
    try std.testing.expectApproxEqAbs(1.0, result.objective_value, check_tolerance);
    try std.testing.expectApproxEqAbs(1.0, result.primal[0], check_tolerance);
    try std.testing.expectApproxEqAbs(0.0, result.primal[1], check_tolerance);
}

test "ipm handles a free variable and objective offset" {
    // min x with x free, x ≥ 5, offset 10. Optimum 15 at x = 5.
    const col_cost = [_]f64{1.0};
    const col_lower = [_]f64{-inf};
    const col_upper = [_]f64{inf};
    const row_lower = [_]f64{5.0};
    const row_upper = [_]f64{inf};
    const col_starts = [_]usize{ 0, 1 };
    const row_indices = [_]RowId{RowId.fromUsizeAssumeValid(0)};
    const values = [_]f64{1.0};
    const view = lp.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &col_cost,
        .col_lower = &col_lower,
        .col_upper = &col_upper,
        .row_lower = &row_lower,
        .row_upper = &row_upper,
        .matrix = matrix.CscView.initAssumeValid(1, 1, &col_starts, &row_indices, &values),
        .objective_sense = .minimize,
        .objective_offset = 10.0,
    };
    var result = try solveExpectOptimal(view);
    defer result.deinit();
    try std.testing.expectEqual(.optimal, result.status);
    try std.testing.expectApproxEqAbs(15.0, result.objective_value, check_tolerance);
    try std.testing.expectApproxEqAbs(5.0, result.primal[0], check_tolerance);
}

test "ipm substitutes fixed columns" {
    // min 2x + y with x fixed at 3, y ≥ 0, row y ≥ 1. Optimum 7 at y = 1.
    const col_cost = [_]f64{ 2.0, 1.0 };
    const col_lower = [_]f64{ 3.0, 0.0 };
    const col_upper = [_]f64{ 3.0, inf };
    const row_lower = [_]f64{1.0};
    const row_upper = [_]f64{inf};
    const col_starts = [_]usize{ 0, 0, 1 };
    const row_indices = [_]RowId{RowId.fromUsizeAssumeValid(0)};
    const values = [_]f64{1.0};
    const view = lp.ProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &col_cost,
        .col_lower = &col_lower,
        .col_upper = &col_upper,
        .row_lower = &row_lower,
        .row_upper = &row_upper,
        .matrix = matrix.CscView.initAssumeValid(1, 2, &col_starts, &row_indices, &values),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var result = try solveExpectOptimal(view);
    defer result.deinit();
    try std.testing.expectEqual(.optimal, result.status);
    try std.testing.expectApproxEqAbs(7.0, result.objective_value, check_tolerance);
    try std.testing.expectApproxEqAbs(3.0, result.primal[0], check_tolerance);
    try std.testing.expectApproxEqAbs(1.0, result.primal[1], check_tolerance);
}

test "ipm duals and reduced costs satisfy stationarity" {
    // Same max LP as the vertex test; verify dual sign convention and that
    // primal feasibility plus stationarity hold in model coordinates.
    const col_cost = [_]f64{ 1.0, 1.0 };
    const col_lower = [_]f64{ 0.0, 0.0 };
    const col_upper = [_]f64{ inf, inf };
    const row_lower = [_]f64{ -inf, -inf };
    const row_upper = [_]f64{ 4.0, 12.0 };
    const col_starts = [_]usize{ 0, 2, 4 };
    const row_indices = [_]RowId{ RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(1), RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(1) };
    const values = [_]f64{ 1.0, 4.0, 2.0, 2.0 };
    const view = lp.ProblemView{
        .num_rows = 2,
        .num_cols = 2,
        .col_cost = &col_cost,
        .col_lower = &col_lower,
        .col_upper = &col_upper,
        .row_lower = &row_lower,
        .row_upper = &row_upper,
        .matrix = matrix.CscView.initAssumeValid(2, 2, &col_starts, &row_indices, &values),
        .objective_sense = .maximize,
        .objective_offset = 0.0,
    };
    var result = try solveExpectOptimal(view);
    defer result.deinit();
    try std.testing.expectEqual(.optimal, result.status);
    // Both rows bind at the optimum, so max-convention row duals are ≥ 0.
    try std.testing.expect(result.dual[0] > -check_tolerance);
    try std.testing.expect(result.dual[1] > -check_tolerance);
    // Stationarity: reduced cost c − Aᵀy vanishes for basic columns.
    try std.testing.expectApproxEqAbs(0.0, result.reduced_cost[0], check_tolerance);
    try std.testing.expectApproxEqAbs(0.0, result.reduced_cost[1], check_tolerance);
    // Complementary objective: cᵀx = yᵀ(row activity) = 10/3.
    const activity = [_]f64{
        result.primal[0] + 2.0 * result.primal[1],
        4.0 * result.primal[0] + 2.0 * result.primal[1],
    };
    const dual_obj = result.dual[0] * activity[0] + result.dual[1] * activity[1];
    try std.testing.expectApproxEqAbs(10.0 / 3.0, dual_obj, 1e-4);
}

test "ipm reports invalid_model for inconsistent bounds" {
    const col_cost = [_]f64{1.0};
    const col_lower = [_]f64{2.0};
    const col_upper = [_]f64{1.0};
    const view = lp.ProblemView{
        .num_rows = 0,
        .num_cols = 1,
        .col_cost = &col_cost,
        .col_lower = &col_lower,
        .col_upper = &col_upper,
        .row_lower = &.{},
        .row_upper = &.{},
        .matrix = matrix.CscView.initAssumeValid(0, 1, &[_]usize{ 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var result = try ipm.solve(std.testing.allocator, view, .{});
    defer result.deinit();
    try std.testing.expectEqual(.invalid_model, result.status);
}

test "ipm reports iteration_limit when the budget is zero" {
    const col_cost = [_]f64{1.0};
    const col_lower = [_]f64{0.0};
    const col_upper = [_]f64{inf};
    const row_lower = [_]f64{1.0};
    const row_upper = [_]f64{inf};
    const col_starts = [_]usize{ 0, 1 };
    const row_indices = [_]RowId{RowId.fromUsizeAssumeValid(0)};
    const values = [_]f64{1.0};
    const view = lp.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &col_cost,
        .col_lower = &col_lower,
        .col_upper = &col_upper,
        .row_lower = &row_lower,
        .row_upper = &row_upper,
        .matrix = matrix.CscView.initAssumeValid(1, 1, &col_starts, &row_indices, &values),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var result = try ipm.solve(std.testing.allocator, view, .{ .max_iterations = 0 });
    defer result.deinit();
    try std.testing.expectEqual(.iteration_limit, result.status);
}
