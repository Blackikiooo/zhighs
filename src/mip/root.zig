//! Correctness-first serial branch-and-bound over the LP solver boundary.

const std = @import("std");
const solver = @import("solver");

pub const Integrality = enum { continuous, integer, binary, semi_continuous, semi_integer };
pub const SolveStatus = enum {
    optimal,
    infeasible,
    unbounded,
    node_limit,
    solution_limit,
    iteration_limit,
    work_limit,
    time_limit,
    interrupted,
    unbounded_or_infeasible,
    invalid_model,
    numerical_failure,
};

pub const ProblemView = struct {
    relaxation: solver.LpProblemView,
    integrality: []const Integrality,

    pub fn validate(self: ProblemView) !void {
        try self.relaxation.validate();
        if (self.integrality.len != self.relaxation.num_cols) return error.DimensionMismatch;
    }
};

pub const SolveControl = struct {
    node_limit: usize = std.math.maxInt(usize),
    solution_limit: usize = std.math.maxInt(usize),
    integer_tolerance: f64 = 1e-6,
    absolute_gap: f64 = 1e-10,
    relative_gap: f64 = 1e-4,
    /// Optional complete primal start in original column coordinates. NaN or
    /// non-finite entries make the candidate unavailable rather than invalid.
    initial_solution: ?[]const f64 = null,
    rounding_heuristic: bool = true,
    lp: solver.LpSolveControl = .{},
};

pub const SolveResult = struct {
    allocator: std.mem.Allocator,
    status: SolveStatus,
    primal: []f64 = &.{},
    objective_value: f64 = 0.0,
    objective_bound: f64 = 0.0,
    nodes: usize = 0,
    lp_iterations: usize = 0,
    lp_work: u64 = 0,
    lp_cold_solves: usize = 0,
    lp_reoptimizations: usize = 0,
    solution_count: usize = 0,

    pub fn deinit(self: *SolveResult) void {
        self.allocator.free(self.primal);
        self.* = undefined;
    }
};

const Node = struct {
    lower: []f64,
    upper: []f64,
    bound: ?f64 = null,

    fn deinit(self: *Node, allocator: std.mem.Allocator) void {
        allocator.free(self.lower);
        allocator.free(self.upper);
        self.* = undefined;
    }
};

pub fn solve(allocator: std.mem.Allocator, problem: ProblemView, control: SolveControl) std.mem.Allocator.Error!SolveResult {
    problem.validate() catch return .{ .allocator = allocator, .status = .invalid_model };
    if (!validControl(control)) return .{ .allocator = allocator, .status = .invalid_model };

    var stack: std.ArrayListUnmanaged(Node) = .empty;
    defer {
        for (stack.items) |*node| node.deinit(allocator);
        stack.deinit(allocator);
    }
    const root_lower = try allocator.dupe(f64, problem.relaxation.col_lower);
    errdefer allocator.free(root_lower);
    const root_upper = try allocator.dupe(f64, problem.relaxation.col_upper);
    errdefer allocator.free(root_upper);
    if (!normalizeRootBounds(problem.integrality, root_lower, root_upper, control.integer_tolerance))
        return .{ .allocator = allocator, .status = .infeasible };
    try stack.append(allocator, .{ .lower = root_lower, .upper = root_upper });

    var incumbent: []f64 = &.{};
    errdefer allocator.free(incumbent);
    var incumbent_objective: f64 = 0.0;
    var have_incumbent = false;
    var nodes: usize = 0;
    var iterations: usize = 0;
    var work: u64 = 0;
    var solutions: usize = 0;
    var lp_session = solver.LpSolveSession.init(allocator);
    defer lp_session.deinit();
    const clock_io = if (control.lp.time_limit_ns != null)
        control.lp.clock_io orelse std.Io.Threaded.global_single_threaded.io()
    else
        null;
    const start_ns = if (clock_io) |io| std.Io.Clock.awake.now(io).nanoseconds else 0;

    if (control.initial_solution) |start| {
        if (try copyFeasibleSolution(allocator, problem, start, control.integer_tolerance)) |candidate| {
            incumbent = candidate.primal;
            incumbent_objective = candidate.objective;
            have_incumbent = true;
            solutions = 1;
            if (solutions >= control.solution_limit)
                return finish(allocator, .solution_limit, incumbent, incumbent_objective, unknownBound(problem.relaxation.objective_sense), nodes, iterations, work, solutions, &lp_session);
        }
    }

    while (stack.items.len != 0) {
        if (nodes >= control.node_limit)
            return finish(allocator, .node_limit, incumbent, incumbent_objective, bestOpenBound(problem.relaxation.objective_sense, stack.items, incumbent_objective, have_incumbent), nodes, iterations, work, solutions, &lp_session);
        if (iterations >= control.lp.max_iterations)
            return finish(allocator, .iteration_limit, incumbent, incumbent_objective, bestOpenBound(problem.relaxation.objective_sense, stack.items, incumbent_objective, have_incumbent), nodes, iterations, work, solutions, &lp_session);
        if (control.lp.work_limit) |limit| if (work >= limit)
            return finish(allocator, .work_limit, incumbent, incumbent_objective, bestOpenBound(problem.relaxation.objective_sense, stack.items, incumbent_objective, have_incumbent), nodes, iterations, work, solutions, &lp_session);

        var node = stack.pop().?;
        defer node.deinit(allocator);
        nodes += 1;
        var relaxation = problem.relaxation;
        relaxation.col_lower = node.lower;
        relaxation.col_upper = node.upper;
        var node_control = control.lp;
        node_control.max_iterations = control.lp.max_iterations - iterations;
        node_control.work_limit = if (control.lp.work_limit) |limit| limit - work else null;
        node_control.initial_basis = if (nodes == 1) control.lp.initial_basis else null;
        node_control.time_limit_ns = remainingTime(control.lp.time_limit_ns, clock_io, start_ns);
        const lp_status = lp_session.solve(relaxation, node_control, .{
            .structure = 0,
            .matrix_values = 0,
            .bounds = nodes,
            .objective = 0,
        });
        iterations = std.math.add(usize, iterations, lp_session.engine.iterations) catch std.math.maxInt(usize);
        work = std.math.add(u64, work, lp_session.engine.work_used) catch std.math.maxInt(u64);

        switch (lp_status) {
            .infeasible => continue,
            .unbounded => return finish(allocator, if (hasIntegerColumns(problem.integrality)) .unbounded_or_infeasible else .unbounded, incumbent, incumbent_objective, 0.0, nodes, iterations, work, solutions, &lp_session),
            .iteration_limit => return finish(allocator, .iteration_limit, incumbent, incumbent_objective, node.bound orelse 0.0, nodes, iterations, work, solutions, &lp_session),
            .work_limit => return finish(allocator, .work_limit, incumbent, incumbent_objective, node.bound orelse 0.0, nodes, iterations, work, solutions, &lp_session),
            .time_limit => return finish(allocator, .time_limit, incumbent, incumbent_objective, node.bound orelse 0.0, nodes, iterations, work, solutions, &lp_session),
            .interrupted => return finish(allocator, .interrupted, incumbent, incumbent_objective, node.bound orelse 0.0, nodes, iterations, work, solutions, &lp_session),
            .invalid_model => return finish(allocator, .invalid_model, incumbent, incumbent_objective, 0.0, nodes, iterations, work, solutions, &lp_session),
            .numerical_failure, .not_implemented => return finish(allocator, .numerical_failure, incumbent, incumbent_objective, 0.0, nodes, iterations, work, solutions, &lp_session),
            .optimal => {},
        }
        const view = lp_session.resultView(relaxation, lp_status) orelse
            return finish(allocator, .numerical_failure, incumbent, incumbent_objective, 0.0, nodes, iterations, work, solutions, &lp_session);
        if (have_incumbent and cannotImprove(problem.relaxation.objective_sense, view.objective_value, incumbent_objective)) continue;

        const violation = chooseViolation(problem, node, view.primal, control.integer_tolerance);
        if (violation != null and control.rounding_heuristic) {
            if (try roundedFeasibleSolution(allocator, problem, view.primal, control.integer_tolerance)) |candidate| {
                if (!have_incumbent or betterObjective(problem.relaxation.objective_sense, candidate.objective, incumbent_objective)) {
                    allocator.free(incumbent);
                    incumbent = candidate.primal;
                    incumbent_objective = candidate.objective;
                    have_incumbent = true;
                    solutions += 1;
                    if (solutions >= control.solution_limit) {
                        const pending_bound = combineBound(problem.relaxation.objective_sense, bestOpenBound(problem.relaxation.objective_sense, stack.items, view.objective_value, true), view.objective_value);
                        return finish(allocator, .solution_limit, incumbent, incumbent_objective, pending_bound, nodes, iterations, work, solutions, &lp_session);
                    }
                } else {
                    allocator.free(candidate.primal);
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
                return finish(allocator, .solution_limit, incumbent, incumbent_objective, bestOpenBound(problem.relaxation.objective_sense, stack.items, incumbent_objective, true), nodes, iterations, work, solutions, &lp_session);
            const open_bound = bestOpenBound(problem.relaxation.objective_sense, stack.items, incumbent_objective, true);
            if (gapSatisfied(open_bound, incumbent_objective, control))
                return finish(allocator, .optimal, incumbent, incumbent_objective, open_bound, nodes, iterations, work, solutions, &lp_session);
            continue;
        }

        const column = violation.?.column;
        const value = view.primal[column];
        switch (violation.?.kind) {
            .integer => {
                if (@floor(value) >= node.lower[column] - control.integer_tolerance) {
                    var child = try cloneNode(allocator, node, view.objective_value);
                    errdefer child.deinit(allocator);
                    child.upper[column] = @min(child.upper[column], @floor(value));
                    if (child.lower[column] <= child.upper[column]) try stack.append(allocator, child) else child.deinit(allocator);
                }
                if (@ceil(value) <= node.upper[column] + control.integer_tolerance) {
                    var child = try cloneNode(allocator, node, view.objective_value);
                    errdefer child.deinit(allocator);
                    child.lower[column] = @max(child.lower[column], @ceil(value));
                    if (child.lower[column] <= child.upper[column]) try stack.append(allocator, child) else child.deinit(allocator);
                }
            },
            .semi_activation => {
                if (node.lower[column] <= 0 and node.upper[column] >= 0) {
                    var zero = try cloneNode(allocator, node, view.objective_value);
                    errdefer zero.deinit(allocator);
                    zero.lower[column] = 0;
                    zero.upper[column] = 0;
                    try stack.append(allocator, zero);
                }
                var active = try cloneNode(allocator, node, view.objective_value);
                errdefer active.deinit(allocator);
                active.lower[column] = @max(active.lower[column], activeLower(problem, column, control.integer_tolerance));
                active.upper[column] = @min(active.upper[column], activeUpper(problem, column, control.integer_tolerance));
                if (active.lower[column] <= active.upper[column]) try stack.append(allocator, active) else active.deinit(allocator);
            },
        }
    }
    if (!have_incumbent) return .{ .allocator = allocator, .status = .infeasible, .nodes = nodes, .lp_iterations = iterations, .lp_work = work, .lp_cold_solves = lp_session.cold_solves, .lp_reoptimizations = lp_session.reoptimizations };
    return finish(allocator, .optimal, incumbent, incumbent_objective, incumbent_objective, nodes, iterations, work, solutions, &lp_session);
}

fn validControl(control: SolveControl) bool {
    return std.math.isFinite(control.integer_tolerance) and control.integer_tolerance >= 0 and
        std.math.isFinite(control.absolute_gap) and control.absolute_gap >= 0 and
        std.math.isFinite(control.relative_gap) and control.relative_gap >= 0;
}

fn normalizeRootBounds(integrality: []const Integrality, lower: []f64, upper: []f64, tolerance: f64) bool {
    for (integrality, lower, upper) |kind, *lo, *up| {
        if (kind == .continuous) continue;
        switch (kind) {
            .continuous => unreachable,
            .binary => {
                lo.* = @max(lo.*, 0.0);
                up.* = @min(up.*, 1.0);
            },
            .integer => {},
            .semi_continuous, .semi_integer => {
                const active_lo = if (kind == .semi_integer and std.math.isFinite(lo.*)) @ceil(lo.* - tolerance) else lo.*;
                const active_up = if (kind == .semi_integer and std.math.isFinite(up.*)) @floor(up.* + tolerance) else up.*;
                lo.* = @min(0.0, active_lo);
                up.* = @max(0.0, active_up);
                continue;
            },
        }
        if (std.math.isFinite(lo.*)) lo.* = @ceil(lo.* - tolerance);
        if (std.math.isFinite(up.*)) up.* = @floor(up.* + tolerance);
        if (lo.* > up.*) return false;
    }
    return true;
}

fn hasIntegerColumns(integrality: []const Integrality) bool {
    for (integrality) |kind| if (kind != .continuous) return true;
    return false;
}

const Violation = struct {
    column: usize,
    kind: enum { integer, semi_activation },
};

fn chooseViolation(problem: ProblemView, node: Node, primal: []const f64, tolerance: f64) ?Violation {
    var best: ?Violation = null;
    var score: f64 = 0;
    for (problem.integrality, primal, 0..) |kind, value, column| {
        if (kind == .continuous) continue;
        if (kind == .semi_continuous or kind == .semi_integer) {
            const zero_allowed = node.lower[column] <= tolerance and node.upper[column] >= -tolerance;
            const active = value >= activeLower(problem, column, tolerance) - tolerance and
                value <= activeUpper(problem, column, tolerance) + tolerance;
            if (zero_allowed and @abs(value) <= tolerance) continue;
            if (!active) return .{ .column = column, .kind = .semi_activation };
            if (kind == .semi_continuous) continue;
        }
        const fractionality = @abs(value - @round(value));
        if (fractionality > tolerance and fractionality > score) {
            best = .{ .column = column, .kind = .integer };
            score = fractionality;
        }
    }
    return best;
}

fn activeLower(problem: ProblemView, column: usize, tolerance: f64) f64 {
    const lower = problem.relaxation.col_lower[column];
    return if (problem.integrality[column] == .semi_integer and std.math.isFinite(lower)) @ceil(lower - tolerance) else lower;
}

fn activeUpper(problem: ProblemView, column: usize, tolerance: f64) f64 {
    const upper = problem.relaxation.col_upper[column];
    return if (problem.integrality[column] == .semi_integer and std.math.isFinite(upper)) @floor(upper + tolerance) else upper;
}

const Candidate = struct { primal: []f64, objective: f64 };

fn copyFeasibleSolution(allocator: std.mem.Allocator, problem: ProblemView, values: []const f64, tolerance: f64) !?Candidate {
    if (values.len != problem.relaxation.num_cols) return null;
    const primal = try allocator.dupe(f64, values);
    errdefer allocator.free(primal);
    for (problem.integrality, primal, problem.relaxation.col_lower, problem.relaxation.col_upper) |kind, *value, lower, upper| {
        if (!std.math.isFinite(value.*)) {
            allocator.free(primal);
            return null;
        }
        switch (kind) {
            .continuous => {},
            .integer => {
                if (@abs(value.* - @round(value.*)) > tolerance) {
                    allocator.free(primal);
                    return null;
                }
                value.* = @round(value.*);
            },
            .binary => {
                if (@abs(value.* - @round(value.*)) > tolerance or value.* < -tolerance or value.* > 1 + tolerance) {
                    allocator.free(primal);
                    return null;
                }
                value.* = @round(value.*);
            },
            .semi_continuous, .semi_integer => {
                if (@abs(value.*) <= tolerance) {
                    value.* = 0;
                    continue;
                }
                if (kind == .semi_integer) {
                    if (@abs(value.* - @round(value.*)) > tolerance) {
                        allocator.free(primal);
                        return null;
                    }
                    value.* = @round(value.*);
                }
            },
        }
        if (value.* < lower - tolerance or value.* > upper + tolerance) {
            allocator.free(primal);
            return null;
        }
    }

    const activity = try allocator.alloc(f64, problem.relaxation.num_rows);
    defer allocator.free(activity);
    @memset(activity, 0);
    for (0..problem.relaxation.num_cols) |column| {
        const begin = problem.relaxation.matrix.col_starts[column];
        const end = problem.relaxation.matrix.col_starts[column + 1];
        for (problem.relaxation.matrix.row_indices[begin..end], problem.relaxation.matrix.values[begin..end]) |row, coefficient|
            activity[row.toUsize()] += coefficient * primal[column];
    }
    for (activity, problem.relaxation.row_lower, problem.relaxation.row_upper) |value, lower, upper| {
        if (value < lower - tolerance or value > upper + tolerance) {
            allocator.free(primal);
            return null;
        }
    }
    var objective = problem.relaxation.objective_offset;
    for (problem.relaxation.col_cost, primal) |cost, value| objective += cost * value;
    return .{ .primal = primal, .objective = objective };
}

fn roundedFeasibleSolution(allocator: std.mem.Allocator, problem: ProblemView, relaxation: []const f64, tolerance: f64) !?Candidate {
    const rounded = try allocator.dupe(f64, relaxation);
    defer allocator.free(rounded);
    for (problem.integrality, rounded, 0..) |kind, *value, column| switch (kind) {
        .continuous => {},
        .integer => value.* = std.math.clamp(@round(value.*), activeLower(problem, column, tolerance), activeUpper(problem, column, tolerance)),
        .binary => value.* = std.math.clamp(@round(value.*), @max(0.0, problem.relaxation.col_lower[column]), @min(1.0, problem.relaxation.col_upper[column])),
        .semi_continuous, .semi_integer => {
            const lower = activeLower(problem, column, tolerance);
            const upper = activeUpper(problem, column, tolerance);
            if (lower > upper) {
                value.* = 0;
                continue;
            }
            const active = if (kind == .semi_integer)
                std.math.clamp(@round(value.*), lower, upper)
            else
                std.math.clamp(value.*, lower, upper);
            value.* = if (@abs(value.*) <= @abs(value.* - active)) 0 else active;
        },
    };
    return copyFeasibleSolution(allocator, problem, rounded, tolerance);
}

fn betterObjective(sense: solver.LpObjectiveSense, candidate: f64, incumbent: f64) bool {
    return switch (sense) {
        .minimize => candidate < incumbent,
        .maximize => candidate > incumbent,
    };
}

fn unknownBound(sense: solver.LpObjectiveSense) f64 {
    return switch (sense) {
        .minimize => -std.math.inf(f64),
        .maximize => std.math.inf(f64),
    };
}

fn cloneNode(allocator: std.mem.Allocator, source: Node, bound: f64) !Node {
    const lower = try allocator.dupe(f64, source.lower);
    errdefer allocator.free(lower);
    return .{ .lower = lower, .upper = try allocator.dupe(f64, source.upper), .bound = bound };
}

fn cannotImprove(sense: solver.LpObjectiveSense, bound: f64, incumbent: f64) bool {
    return switch (sense) {
        .minimize => bound >= incumbent,
        .maximize => bound <= incumbent,
    };
}

fn gapSatisfied(bound: f64, incumbent: f64, control: SolveControl) bool {
    const gap = @abs(incumbent - bound);
    return gap <= control.absolute_gap or gap / @max(1.0, @abs(incumbent)) <= control.relative_gap;
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

fn bestOpenBound(sense: solver.LpObjectiveSense, nodes: []const Node, incumbent: f64, have_incumbent: bool) f64 {
    var result: ?f64 = null;
    for (nodes) |node| {
        if (node.bound) |bound| {
            result = if (result) |old| switch (sense) {
                .minimize => @min(old, bound),
                .maximize => @max(old, bound),
            } else bound;
        }
    }
    return result orelse if (have_incumbent) incumbent else 0.0;
}

fn combineBound(sense: solver.LpObjectiveSense, left: f64, right: f64) f64 {
    return switch (sense) {
        .minimize => @min(left, right),
        .maximize => @max(left, right),
    };
}

fn finish(allocator: std.mem.Allocator, status: SolveStatus, primal: []f64, objective: f64, bound: f64, nodes: usize, iterations: usize, work: u64, solutions: usize, session: *const solver.LpSolveSession) SolveResult {
    return .{
        .allocator = allocator,
        .status = status,
        .primal = primal,
        .objective_value = objective,
        .objective_bound = bound,
        .nodes = nodes,
        .lp_iterations = iterations,
        .lp_work = work,
        .lp_cold_solves = session.cold_solves,
        .lp_reoptimizations = session.reoptimizations,
        .solution_count = solutions,
    };
}

test "branch-and-bound solves a fractional binary relaxation" {
    const matrix = @import("matrix");
    const foundation = @import("foundation");
    const relaxation = solver.LpProblemView{
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
    var result = try solve(std.testing.allocator, .{ .relaxation = relaxation, .integrality = &.{ .integer, .integer } }, .{});
    defer result.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, result.status);
    try std.testing.expectApproxEqAbs(@as(f64, 2), result.objective_value, 1e-9);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0 }, result.primal);
    try std.testing.expect(result.nodes >= 2);
    try std.testing.expectEqual(result.nodes, result.lp_cold_solves + result.lp_reoptimizations);

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
    const matrix = @import("matrix");
    const foundation = @import("foundation");
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
    const matrix = @import("matrix");
    const foundation = @import("foundation");
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
