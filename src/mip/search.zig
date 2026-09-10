//! Plugin-driven serial branch-and-bound search over the LP solver boundary.
//!
//! The engine owns node management, LP reoptimization, incumbent bookkeeping,
//! and limit enforcement. Every algorithmic decision is delegated to
//! components registered through a `framework.PluginRegistry`:
//!
//! * branching rules pick the disjunction that splits a fractional node,
//! * primal heuristics propose incumbent candidates (validated here),
//! * a node selector picks the next open node to process.
//!
//! The default component set reproduces the historical behaviour exactly:
//! most-fractional branching, LP rounding, and depth-first search.

const std = @import("std");
const solver = @import("solver");
const plugin = @import("plugin");
const framework = @import("framework");
const builtin_plugins = @import("plugins_builtin");
const matrix = @import("matrix");
const foundation = @import("foundation");
const types = @import("types.zig");
const domain = @import("domain.zig");
const incumbent_candidate = @import("incumbent.zig");
const objective_bounds = @import("bounds.zig");
const NodeQueue = @import("node_queue.zig").NodeQueue;
const Relaxation = @import("relaxation.zig").Relaxation;

pub const Integrality = types.Integrality;
pub const ProblemView = types.ProblemView;
pub const SolveStatus = types.SolveStatus;
pub const SolveComponents = types.SolveComponents;
pub const SolveControl = types.SolveControl;
pub const SolveResult = types.SolveResult;

pub fn solve(allocator: std.mem.Allocator, problem: ProblemView, control: SolveControl) std.mem.Allocator.Error!SolveResult {
    problem.validate() catch return .{ .allocator = allocator, .status = .invalid_model };
    if (!types.validControl(control)) return .{ .allocator = allocator, .status = .invalid_model };

    var registry = framework.PluginRegistry.init(allocator);
    defer registry.deinit();

    var default_branching = builtin_plugins.branching.MostFractional{};
    var default_heuristic = builtin_plugins.heuristics.Rounding{};
    if (control.components.branching_rule) |rule| {
        try registry.addBranchingRule(rule, 0);
    } else {
        try registry.addBranchingRule(plugin.bindBranching(builtin_plugins.branching.MostFractional, &default_branching), 0);
    }
    if (control.components.node_selector) |selector| try registry.addNodeSelector(selector, 0);
    if (control.components.heuristics.len != 0) {
        for (control.components.heuristics) |heuristic| try registry.addHeuristic(heuristic, 0);
    } else if (control.rounding_heuristic) {
        try registry.addHeuristic(plugin.bindHeuristic(builtin_plugins.heuristics.Rounding, &default_heuristic), 0);
    }
    const node_selector = registry.primaryNodeSelector();

    var current_domain = (try domain.Domain.init(
        allocator,
        problem.integrality,
        problem.relaxation.col_lower,
        problem.relaxation.col_upper,
        control.integer_tolerance,
    )) orelse return .{ .allocator = allocator, .status = .infeasible };
    defer current_domain.deinit(allocator);
    var change_arena = domain.DomainChangeArena{};
    defer change_arena.deinit(allocator);

    var open_nodes = NodeQueue{};
    defer open_nodes.deinit(allocator);
    try open_nodes.append(allocator, .{});

    var incumbent: []f64 = &.{};
    errdefer allocator.free(incumbent);
    var incumbent_objective: f64 = 0.0;
    var have_incumbent = false;
    var nodes: usize = 0;
    var iterations: usize = 0;
    var work: u64 = 0;
    var solutions: usize = 0;
    var lp_relaxation = Relaxation.init(allocator);
    defer lp_relaxation.deinit();
    const clock_io = if (control.lp.time_limit_ns != null)
        control.lp.clock_io orelse std.Io.Threaded.global_single_threaded.io()
    else
        null;
    const start_ns = if (clock_io) |io| std.Io.Clock.awake.now(io).nanoseconds else 0;

    if (control.initial_solution) |start| {
        if (try incumbent_candidate.copyIfFeasible(allocator, problem, start, control.integer_tolerance)) |candidate| {
            incumbent = candidate.primal;
            incumbent_objective = candidate.objective;
            have_incumbent = true;
            solutions = 1;
            if (solutions >= control.solution_limit)
                return finish(allocator, .solution_limit, incumbent, incumbent_objective, objective_bounds.unknown(problem.relaxation.objective_sense), nodes, iterations, work, solutions, &lp_relaxation);
        }
    }

    while (open_nodes.len() != 0) {
        if (nodes >= control.node_limit)
            return finish(allocator, .node_limit, incumbent, incumbent_objective, objective_bounds.bestOpen(problem.relaxation.objective_sense, open_nodes.items(), incumbent_objective, have_incumbent), nodes, iterations, work, solutions, &lp_relaxation);
        if (iterations >= control.lp.max_iterations)
            return finish(allocator, .iteration_limit, incumbent, incumbent_objective, objective_bounds.bestOpen(problem.relaxation.objective_sense, open_nodes.items(), incumbent_objective, have_incumbent), nodes, iterations, work, solutions, &lp_relaxation);
        if (control.lp.work_limit) |limit| if (work >= limit)
            return finish(allocator, .work_limit, incumbent, incumbent_objective, objective_bounds.bestOpen(problem.relaxation.objective_sense, open_nodes.items(), incumbent_objective, have_incumbent), nodes, iterations, work, solutions, &lp_relaxation);

        var node = (try open_nodes.take(allocator, problem.relaxation.objective_sense, node_selector)) orelse break;
        defer node.deinit(allocator);
        if (!current_domain.load(&change_arena, node.change_head)) continue;
        nodes += 1;
        var relaxation = problem.relaxation;
        relaxation.col_lower = current_domain.lower;
        relaxation.col_upper = current_domain.upper;
        var node_control = control.lp;
        node_control.max_iterations = control.lp.max_iterations - iterations;
        node_control.work_limit = if (control.lp.work_limit) |limit| limit - work else null;
        node_control.initial_basis = if (nodes == 1) control.lp.initial_basis else null;
        node_control.time_limit_ns = remainingTime(control.lp.time_limit_ns, clock_io, start_ns);
        const lp_status = lp_relaxation.solve(relaxation, node_control, nodes, .{
            .structure = 0,
            .matrix_values = 0,
            .bounds = nodes,
            .objective = 0,
        });
        iterations = std.math.add(usize, iterations, lp_relaxation.lastIterations()) catch std.math.maxInt(usize);
        work = std.math.add(u64, work, lp_relaxation.lastWork()) catch std.math.maxInt(u64);

        switch (lp_status) {
            .infeasible => continue,
            .unbounded => return finish(allocator, if (domain.hasDiscreteColumns(problem.integrality)) .unbounded_or_infeasible else .unbounded, incumbent, incumbent_objective, 0.0, nodes, iterations, work, solutions, &lp_relaxation),
            .iteration_limit => return finish(allocator, .iteration_limit, incumbent, incumbent_objective, node.bound orelse 0.0, nodes, iterations, work, solutions, &lp_relaxation),
            .work_limit => return finish(allocator, .work_limit, incumbent, incumbent_objective, node.bound orelse 0.0, nodes, iterations, work, solutions, &lp_relaxation),
            .time_limit => return finish(allocator, .time_limit, incumbent, incumbent_objective, node.bound orelse 0.0, nodes, iterations, work, solutions, &lp_relaxation),
            .interrupted => return finish(allocator, .interrupted, incumbent, incumbent_objective, node.bound orelse 0.0, nodes, iterations, work, solutions, &lp_relaxation),
            .invalid_model => return finish(allocator, .invalid_model, incumbent, incumbent_objective, 0.0, nodes, iterations, work, solutions, &lp_relaxation),
            .numerical_failure, .not_implemented => return finish(allocator, .numerical_failure, incumbent, incumbent_objective, 0.0, nodes, iterations, work, solutions, &lp_relaxation),
            .optimal => {},
        }
        const view = lp_relaxation.resultView(relaxation, lp_status) orelse
            return finish(allocator, .numerical_failure, incumbent, incumbent_objective, 0.0, nodes, iterations, work, solutions, &lp_relaxation);
        if (have_incumbent and objective_bounds.cannotImprove(problem.relaxation.objective_sense, view.objective_value, incumbent_objective)) continue;

        const branch_context = plugin.BranchingContext{
            .problem = problem,
            .node_lower = current_domain.lower,
            .node_upper = current_domain.upper,
            .lp_primal = view.primal,
            .active_lower = current_domain.active_lower,
            .active_upper = current_domain.active_upper,
            .tolerance = control.integer_tolerance,
        };
        var violation: ?plugin.BranchingCandidate = null;
        for (registry.branchingRuleEntries()) |entry| {
            violation = entry.rule.select(branch_context);
            if (violation != null) break;
        }

        if (violation != null) {
            const heuristic_context = plugin.HeuristicContext{
                .allocator = allocator,
                .problem = problem,
                .node_lower = current_domain.lower,
                .node_upper = current_domain.upper,
                .lp_primal = view.primal,
                .active_lower = current_domain.active_lower,
                .active_upper = current_domain.active_upper,
                .tolerance = control.integer_tolerance,
            };
            for (registry.heuristicEntries()) |entry| {
                const proposed = entry.heuristic.trySolution(heuristic_context) catch return error.OutOfMemory;
                const raw = proposed orelse continue;
                defer allocator.free(raw);
                if (try incumbent_candidate.copyIfFeasible(allocator, problem, raw, control.integer_tolerance)) |candidate| {
                    if (!have_incumbent or objective_bounds.betterObjective(problem.relaxation.objective_sense, candidate.objective, incumbent_objective)) {
                        allocator.free(incumbent);
                        incumbent = candidate.primal;
                        incumbent_objective = candidate.objective;
                        have_incumbent = true;
                        solutions += 1;
                        if (solutions >= control.solution_limit) {
                            const pending_bound = objective_bounds.combine(problem.relaxation.objective_sense, objective_bounds.bestOpen(problem.relaxation.objective_sense, open_nodes.items(), view.objective_value, true), view.objective_value);
                            return finish(allocator, .solution_limit, incumbent, incumbent_objective, pending_bound, nodes, iterations, work, solutions, &lp_relaxation);
                        }
                    } else {
                        allocator.free(candidate.primal);
                    }
                }
            }
        }

        if (violation == null) {
            const candidate = try allocator.dupe(f64, view.primal);
            allocator.free(incumbent);
            incumbent = candidate;
            incumbent_objective = view.objective_value;
            have_incumbent = true;
            solutions += 1;
            if (solutions >= control.solution_limit)
                return finish(allocator, .solution_limit, incumbent, incumbent_objective, objective_bounds.bestOpen(problem.relaxation.objective_sense, open_nodes.items(), incumbent_objective, true), nodes, iterations, work, solutions, &lp_relaxation);
            const open_bound = objective_bounds.bestOpen(problem.relaxation.objective_sense, open_nodes.items(), incumbent_objective, true);
            if (objective_bounds.gapSatisfied(open_bound, incumbent_objective, control))
                return finish(allocator, .optimal, incumbent, incumbent_objective, open_bound, nodes, iterations, work, solutions, &lp_relaxation);
            continue;
        }

        const column = violation.?.column;
        const value = view.primal[column];
        switch (violation.?.kind) {
            .integer => {
                if (@floor(value) >= current_domain.lower[column] - control.integer_tolerance) {
                    var child = try node.childWith(&change_arena, allocator, view.objective_value, .{
                        .column = column,
                        .upper = @floor(value),
                    });
                    errdefer child.deinit(allocator);
                    try open_nodes.append(allocator, child);
                }
                if (@ceil(value) <= current_domain.upper[column] + control.integer_tolerance) {
                    var child = try node.childWith(&change_arena, allocator, view.objective_value, .{
                        .column = column,
                        .lower = @ceil(value),
                    });
                    errdefer child.deinit(allocator);
                    try open_nodes.append(allocator, child);
                }
            },
            .semi_activation => {
                if (current_domain.lower[column] <= 0 and current_domain.upper[column] >= 0) {
                    var zero = try node.childWith(&change_arena, allocator, view.objective_value, .{
                        .column = column,
                        .lower = 0,
                        .upper = 0,
                    });
                    errdefer zero.deinit(allocator);
                    try open_nodes.append(allocator, zero);
                }
                const active_lower = @max(current_domain.lower[column], current_domain.active_lower[column]);
                const active_upper = @min(current_domain.upper[column], current_domain.active_upper[column]);
                if (active_lower <= active_upper) {
                    var active = try node.childWith(&change_arena, allocator, view.objective_value, .{
                        .column = column,
                        .lower = active_lower,
                        .upper = active_upper,
                    });
                    errdefer active.deinit(allocator);
                    try open_nodes.append(allocator, active);
                }
            },
        }
    }
    if (!have_incumbent) {
        const statistics = lp_relaxation.statistics();
        return .{ .allocator = allocator, .status = .infeasible, .nodes = nodes, .lp_iterations = iterations, .lp_work = work, .lp_cold_solves = statistics.cold_solves, .lp_reoptimizations = statistics.reoptimizations };
    }
    return finish(allocator, .optimal, incumbent, incumbent_objective, incumbent_objective, nodes, iterations, work, solutions, &lp_relaxation);
}

fn remainingTime(limit: ?u64, clock_io: ?std.Io, start_ns: i96) ?u64 {
    const total = limit orelse return null;
    if (total == 0) return 0;
    const io = clock_io orelse return total;
    const now = std.Io.Clock.awake.now(io).nanoseconds;
    if (now <= start_ns) return total;
    const elapsed: u128 = @intCast(now - start_ns);
    if (elapsed >= total) return 0;
    return total - @as(u64, @intCast(elapsed));
}

fn finish(allocator: std.mem.Allocator, status: SolveStatus, primal: []f64, objective: f64, bound: f64, nodes: usize, iterations: usize, work: u64, solutions: usize, relaxation: *const Relaxation) SolveResult {
    const statistics = relaxation.statistics();
    return .{
        .allocator = allocator,
        .status = status,
        .primal = primal,
        .objective_value = objective,
        .objective_bound = bound,
        .nodes = nodes,
        .lp_iterations = iterations,
        .lp_work = work,
        .lp_cold_solves = statistics.cold_solves,
        .lp_reoptimizations = statistics.reoptimizations,
        .solution_count = solutions,
    };
}

fn fractionalBinaryRelaxation() solver.LpProblemView {
    return .{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 1.0, 1.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ 1.0, 1.0 },
        .row_lower = &.{1.5},
        .row_upper = &.{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &.{ 0, 1, 2 }, &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) }, &.{ 1.0, 1.0 }),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
}

test "branch-and-bound solves a fractional binary relaxation" {
    const relaxation = fractionalBinaryRelaxation();
    var result = try solve(std.testing.allocator, .{ .relaxation = relaxation, .integrality = &.{ .integer, .integer } }, .{});
    defer result.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, result.status);
    try std.testing.expectApproxEqAbs(@as(f64, 2), result.objective_value, 1e-9);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0 }, result.primal);
    try std.testing.expect(result.nodes >= 2);
    try std.testing.expectEqual(result.nodes, result.lp_cold_solves + result.lp_reoptimizations);
    try std.testing.expectEqual(@as(usize, 1), result.lp_cold_solves);
    try std.testing.expectEqual(result.nodes - 1, result.lp_reoptimizations);

    var limited = try solve(std.testing.allocator, .{ .relaxation = relaxation, .integrality = &.{ .integer, .integer } }, .{ .lp = .{ .max_iterations = 1 } });
    defer limited.deinit();
    try std.testing.expectEqual(SolveStatus.iteration_limit, limited.status);
    try std.testing.expect(limited.lp_iterations <= 1);

    var work_limited = try solve(std.testing.allocator, .{ .relaxation = relaxation, .integrality = &.{ .integer, .integer } }, .{ .lp = .{ .work_limit = 0 } });
    defer work_limited.deinit();
    try std.testing.expectEqual(SolveStatus.work_limit, work_limited.status);
    try std.testing.expectEqual(@as(u64, 0), work_limited.lp_work);
    try std.testing.expectEqual(@as(usize, 0), work_limited.nodes);

    var time_limited = try solve(std.testing.allocator, .{ .relaxation = relaxation, .integrality = &.{ .integer, .integer } }, .{ .lp = .{ .time_limit_ns = 0 } });
    defer time_limited.deinit();
    try std.testing.expectEqual(SolveStatus.time_limit, time_limited.status);

    var heuristic = try solve(std.testing.allocator, .{ .relaxation = relaxation, .integrality = &.{ .integer, .integer } }, .{ .solution_limit = 1 });
    defer heuristic.deinit();
    try std.testing.expectEqual(SolveStatus.solution_limit, heuristic.status);
    try std.testing.expectEqual(@as(usize, 1), heuristic.nodes);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0 }, heuristic.primal);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), heuristic.objective_value, 1e-9);
}

test "branch-and-bound reoptimizes stable node LP mappings" {
    const relaxation = solver.LpProblemView{
        .num_rows = 2,
        .num_cols = 2,
        .col_cost = &.{ 1.0, 0.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ 10.0, 10.0 },
        .row_lower = &.{ -std.math.inf(f64), -std.math.inf(f64) },
        .row_upper = &.{ 4.5, 9.0 },
        .matrix = matrix.CscView.initAssumeValid(2, 2, &.{ 0, 2, 4 }, &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1) }, &.{ 1.0, 1.0, 1.0, -1.0 }),
        .objective_sense = .maximize,
        .objective_offset = 0,
    };
    var result = try solve(std.testing.allocator, .{ .relaxation = relaxation, .integrality = &.{ .integer, .continuous } }, .{});
    defer result.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, result.status);
    try std.testing.expectApproxEqAbs(@as(f64, 4), result.objective_value, 1e-9);
    try std.testing.expect(result.lp_reoptimizations >= 1);
    try std.testing.expectEqual(result.nodes, result.lp_cold_solves + result.lp_reoptimizations);
}

test "branch-and-bound proves integer infeasibility" {
    const relaxation = solver.LpProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &.{0.0},
        .col_lower = &.{0.0},
        .col_upper = &.{1.0},
        .row_lower = &.{0.5},
        .row_upper = &.{0.5},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &.{ 0, 1 }, &.{foundation.RowId.fromUsizeAssumeValid(0)}, &.{1.0}),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var result = try solve(std.testing.allocator, .{ .relaxation = relaxation, .integrality = &.{.integer} }, .{});
    defer result.deinit();
    try std.testing.expectEqual(SolveStatus.infeasible, result.status);
}

test "plugin components can override node selection and branching" {
    const FirstFractional = struct {
        pub const name = "first-fractional-test";
        pub fn select(_: *const @This(), ctx: plugin.BranchingContext) ?plugin.BranchingCandidate {
            for (ctx.problem.integrality, ctx.lp_primal, 0..) |kind, value, column| {
                if (kind == .continuous or kind == .semi_continuous) continue;
                if (@abs(value - @round(value)) > ctx.tolerance) return .{ .column = column, .kind = .integer };
            }
            return null;
        }
    };
    const Fifo = struct {
        pub const name = "fifo-test";
        pub fn select(_: *const @This(), ctx: plugin.NodeSelectionContext) ?usize {
            return if (ctx.open.len == 0) null else 0;
        }
    };
    var branching = FirstFractional{};
    var selector = Fifo{};
    var result = try solve(std.testing.allocator, .{
        .relaxation = fractionalBinaryRelaxation(),
        .integrality = &.{ .integer, .integer },
    }, .{ .components = .{
        .branching_rule = plugin.bindBranching(FirstFractional, &branching),
        .node_selector = plugin.bindNodeSelector(Fifo, &selector),
    } });
    defer result.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, result.status);
    try std.testing.expectApproxEqAbs(@as(f64, 2), result.objective_value, 1e-9);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0 }, result.primal);
}

test "best-bound node selection preserves optimal objective" {
    var selector = builtin_plugins.node_selectors.BestBound{};
    var result = try solve(std.testing.allocator, .{
        .relaxation = fractionalBinaryRelaxation(),
        .integrality = &.{ .integer, .integer },
    }, .{ .components = .{
        .node_selector = plugin.bindNodeSelector(builtin_plugins.node_selectors.BestBound, &selector),
    } });
    defer result.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, result.status);
    try std.testing.expectApproxEqAbs(@as(f64, 2), result.objective_value, 1e-9);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0 }, result.primal);
}

test "disabling the rounding heuristic keeps the search correct" {
    var result = try solve(std.testing.allocator, .{
        .relaxation = fractionalBinaryRelaxation(),
        .integrality = &.{ .integer, .integer },
    }, .{ .rounding_heuristic = false });
    defer result.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, result.status);
    try std.testing.expectApproxEqAbs(@as(f64, 2), result.objective_value, 1e-9);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0 }, result.primal);
}
