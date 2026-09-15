//! Bounded strong-branch LP probes for unreliable pseudocost candidates.
//!
//! Probes use a private LP session and reusable bound buffers, so exploratory
//! solves cannot disturb the main node-relaxation basis. Candidate visitation
//! uses generation stamps rather than clearing a boolean array at every node.

const std = @import("std");
const solver = @import("solver");
const plugin = @import("plugin");
const Pseudocost = @import("pseudocost.zig").Pseudocost;
const BranchOrigin = @import("node.zig").BranchOrigin;
const BranchDirection = @import("node.zig").BranchDirection;

pub const Control = struct {
    max_candidates: usize = 2,
    max_iterations_per_probe: usize = 20,
    max_work_per_node: u64 = 40,
};

pub const Statistics = struct {
    probes: usize = 0,
    iterations: usize = 0,
    work: u64 = 0,
};

pub const RefineResult = struct {
    statistics: Statistics = .{},
    limit_status: ?solver.LpSolveStatus = null,
};

const ProbeResult = struct {
    status: solver.LpSolveStatus,
    objective: ?f64 = null,
    iterations: usize = 0,
    work: u64 = 0,
};

pub const DirectionResult = struct {
    infeasible: bool = false,
    bound: ?f64 = null,
};

const CachedDirection = struct {
    generation: u32 = 0,
    result: DirectionResult = .{},
};

pub const StrongBranching = struct {
    allocator: std.mem.Allocator,
    session: solver.LpSolveSession,
    lower: []f64,
    upper: []f64,
    visited: []u32,
    down_results: []CachedDirection,
    up_results: []CachedDirection,
    generation: u32 = 0,
    revision: usize = 0,

    pub fn init(allocator: std.mem.Allocator, num_columns: usize) std.mem.Allocator.Error!StrongBranching {
        const lower = try allocator.alloc(f64, num_columns);
        errdefer allocator.free(lower);
        const upper = try allocator.alloc(f64, num_columns);
        errdefer allocator.free(upper);
        const visited = try allocator.alloc(u32, num_columns);
        errdefer allocator.free(visited);
        @memset(visited, 0);
        const down_results = try allocator.alloc(CachedDirection, num_columns);
        errdefer allocator.free(down_results);
        @memset(down_results, .{});
        const up_results = try allocator.alloc(CachedDirection, num_columns);
        @memset(up_results, .{});
        return .{
            .allocator = allocator,
            .session = solver.LpSolveSession.init(allocator),
            .lower = lower,
            .upper = upper,
            .visited = visited,
            .down_results = down_results,
            .up_results = up_results,
        };
    }

    pub fn deinit(self: *StrongBranching) void {
        self.session.deinit();
        self.allocator.free(self.lower);
        self.allocator.free(self.upper);
        self.allocator.free(self.visited);
        self.allocator.free(self.down_results);
        self.allocator.free(self.up_results);
        self.* = undefined;
    }

    pub fn refine(
        self: *StrongBranching,
        problem: solver.LpProblemView,
        integrality: []const plugin.Integrality,
        primal: []const f64,
        parent_objective: f64,
        pseudocost: *Pseudocost,
        lp_control: solver.LpSolveControl,
        control: Control,
        tolerance: f64,
    ) RefineResult {
        if (control.max_candidates == 0 or control.max_iterations_per_probe == 0 or control.max_work_per_node == 0)
            return .{};

        const scores = pseudocost.prepareScores(primal);
        var fractional_candidates: usize = 0;
        for (integrality, primal) |kind, value| {
            if ((kind == .integer or kind == .binary) and @abs(value - @round(value)) > tolerance)
                fractional_candidates += 1;
        }
        if (fractional_candidates < 2) return .{};

        self.nextGeneration();
        var result = RefineResult{};
        var candidates: usize = 0;
        var iterations_remaining = lp_control.max_iterations;
        var work_remaining = @min(lp_control.work_limit orelse std.math.maxInt(u64), control.max_work_per_node);
        const clock_io = if (lp_control.time_limit_ns != null)
            lp_control.clock_io orelse std.Io.Threaded.global_single_threaded.io()
        else
            null;
        const start_ns = if (clock_io) |io| std.Io.Clock.awake.now(io).nanoseconds else 0;

        while (candidates < control.max_candidates and iterations_remaining != 0 and work_remaining != 0) : (candidates += 1) {
            const column = self.bestUnreliableCandidate(integrality, primal, scores, tolerance) orelse break;
            self.visited[column] = self.generation;
            self.down_results[column] = .{ .generation = self.generation };
            self.up_results[column] = .{ .generation = self.generation };
            const value = primal[column];
            for ([_]BranchDirection{ .down, .up }) |direction| {
                const reliable = switch (direction) {
                    .down => scores[column].down_reliable,
                    .up => scores[column].up_reliable,
                };
                if (reliable or iterations_remaining == 0 or work_remaining == 0) continue;

                var probe_control = lp_control;
                probe_control.initial_basis = null;
                probe_control.max_iterations = @min(iterations_remaining, control.max_iterations_per_probe);
                probe_control.work_limit = work_remaining;
                probe_control.time_limit_ns = remainingTime(lp_control.time_limit_ns, clock_io, start_ns);
                if (probe_control.time_limit_ns) |remaining| if (remaining == 0) {
                    result.limit_status = .time_limit;
                    return result;
                };

                const probe_result = self.probe(problem, value, column, direction, probe_control);
                const cached = switch (direction) {
                    .down => &self.down_results[column],
                    .up => &self.up_results[column],
                };
                cached.result = .{
                    .infeasible = probe_result.status == .infeasible,
                    .bound = probe_result.objective,
                };
                result.statistics.probes += 1;
                result.statistics.iterations = std.math.add(usize, result.statistics.iterations, probe_result.iterations) catch std.math.maxInt(usize);
                result.statistics.work = std.math.add(u64, result.statistics.work, probe_result.work) catch std.math.maxInt(u64);
                iterations_remaining -|= probe_result.iterations;
                work_remaining -|= probe_result.work;

                if (probe_result.objective) |objective| {
                    const distance = switch (direction) {
                        .down => value - @floor(value),
                        .up => @ceil(value) - value,
                    };
                    pseudocost.observe(problem.objective_sense, .{
                        .column = column,
                        .direction = direction,
                        .distance = distance,
                        .parent_objective = parent_objective,
                    }, objective);
                }
                switch (probe_result.status) {
                    .time_limit, .interrupted => {
                        result.limit_status = probe_result.status;
                        return result;
                    },
                    else => {},
                }
            }
        }
        return result;
    }

    pub fn directionResult(self: *const StrongBranching, column: usize, direction: BranchDirection) ?DirectionResult {
        if (column >= self.visited.len) return null;
        const cached = switch (direction) {
            .down => self.down_results[column],
            .up => self.up_results[column],
        };
        return if (cached.generation == self.generation) cached.result else null;
    }

    fn probe(self: *StrongBranching, problem: solver.LpProblemView, value: f64, column: usize, direction: BranchDirection, control: solver.LpSolveControl) ProbeResult {
        @memcpy(self.lower, problem.col_lower);
        @memcpy(self.upper, problem.col_upper);
        switch (direction) {
            .down => self.upper[column] = @min(self.upper[column], @floor(value)),
            .up => self.lower[column] = @max(self.lower[column], @ceil(value)),
        }
        if (self.lower[column] > self.upper[column]) return .{ .status = .infeasible };

        var probe_problem = problem;
        probe_problem.col_lower = self.lower;
        probe_problem.col_upper = self.upper;
        self.revision +%= 1;
        const status = self.session.solveRelaxation(probe_problem, control, .{
            .structure = 0,
            .matrix_values = 0,
            .bounds = self.revision,
            .objective = 0,
        });
        const objective = if (status == .optimal)
            if (self.session.resultView(probe_problem, status)) |view| view.objective_value else null
        else
            null;
        return .{
            .status = status,
            .objective = objective,
            .iterations = self.session.engine.iterations,
            .work = self.session.engine.work_used,
        };
    }

    fn bestUnreliableCandidate(self: *const StrongBranching, integrality: []const plugin.Integrality, primal: []const f64, scores: []const plugin.PseudocostScore, tolerance: f64) ?usize {
        var best: ?usize = null;
        var best_score: f64 = -1.0;
        for (integrality, primal, scores, 0..) |kind, value, score, column| {
            if (self.visited[column] == self.generation or (kind != .integer and kind != .binary)) continue;
            if (@abs(value - @round(value)) <= tolerance or (score.down_reliable and score.up_reliable)) continue;
            const product = @max(score.down_gain, 1e-6) * @max(score.up_gain, 1e-6);
            if (product > best_score) {
                best = column;
                best_score = product;
            }
        }
        return best;
    }

    fn nextGeneration(self: *StrongBranching) void {
        self.generation +%= 1;
        if (self.generation == 0) {
            @memset(self.visited, 0);
            self.generation = 1;
        }
    }
};

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

test "strong branching refines both directions without touching caller bounds" {
    const foundation = @import("foundation");
    const matrix = @import("matrix");
    const problem = solver.LpProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 1.0, 0.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ 2.0, 2.0 },
        .row_lower = &.{0.5},
        .row_upper = &.{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &.{ 0, 1, 1 }, &.{foundation.RowId.fromUsizeAssumeValid(0)}, &.{1.0}),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var strong = try StrongBranching.init(std.testing.allocator, 2);
    defer strong.deinit();
    var pseudocost = try Pseudocost.init(std.testing.allocator, 2, 2);
    defer pseudocost.deinit(std.testing.allocator);
    const result = strong.refine(problem, &.{ .integer, .integer }, &.{ 0.5, 0.5 }, 0.5, &pseudocost, .{}, .{ .max_candidates = 1 }, 1e-6);
    try std.testing.expectEqual(@as(usize, 2), result.statistics.probes);
    try std.testing.expectEqual(@as(usize, 1), pseudocost.observations);
    try std.testing.expect(strong.directionResult(0, .down).?.infeasible);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), strong.directionResult(0, .up).?.bound.?, 1e-9);
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 0.0 }, problem.col_lower);
    try std.testing.expectEqualSlices(f64, &.{ 2.0, 2.0 }, problem.col_upper);
}
