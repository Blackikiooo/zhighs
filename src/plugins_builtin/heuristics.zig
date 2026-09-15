//! Built-in primal heuristics.
//!
//! Algorithm implementations live here; registration, ordering, timing, and
//! lifecycle remain the responsibility of `framework`.

const std = @import("std");
const plugin = @import("plugin");
const solver = @import("solver");

/// Simple LP-solution rounding. Snaps each integer-like column to its nearest
/// integral value inside the active column bounds and resolves semi-continuous
/// / semi-integer columns to either their zero branch or their active value,
/// whichever is closer to the LP point. The engine validates the candidate, so
/// an infeasible rounding is simply discarded.
pub const Rounding = struct {
    pub const name = "rounding";

    pub fn trySolution(_: *Rounding, ctx: plugin.HeuristicContext) std.mem.Allocator.Error!?[]f64 {
        const rounded = try ctx.allocator.dupe(f64, ctx.lp_primal);
        errdefer ctx.allocator.free(rounded);
        for (ctx.problem.integrality, rounded, 0..) |kind, *value, column| switch (kind) {
            .continuous => {},
            .integer => value.* = std.math.clamp(@round(value.*), ctx.active_lower[column], ctx.active_upper[column]),
            .binary => value.* = std.math.clamp(@round(value.*), @max(0.0, ctx.problem.relaxation.col_lower[column]), @min(1.0, ctx.problem.relaxation.col_upper[column])),
            .semi_continuous, .semi_integer => {
                const lower = ctx.active_lower[column];
                const upper = ctx.active_upper[column];
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
        return rounded;
    }
};

/// Fractionality of `value`: its distance to the nearest integer.
fn fractionality(value: f64) f64 {
    return @abs(value - @round(value));
}

/// Diving heuristic. Starting from the node LP point it repeatedly picks the
/// fractional integer-like column closest to integrality, rounds and fixes it,
/// and re-solves the LP relaxation until every integer-like column is
/// integral, a sub-LP fails, or the round cap is reached. Each fix is
/// permanent within one dive (no backtracking), so the loop always terminates
/// after at most `max_rounds` LP resolves. Semi-continuous / semi-integer
/// columns are never fixed during the dive; they are resolved to their zero
/// branch or active value at the end, exactly as in `Rounding`.
pub const Diving = struct {
    pub const name = "diving";

    /// Maximum number of bound-fixing LP resolves per invocation.
    max_rounds: usize = 64,
    /// Iteration budget for each sub-LP solve.
    max_lp_iterations: usize = 100_000,

    /// Reusable LP session kept across `trySolution` invocations so that the
    /// expensive engine state (basis, factorization, work arrays) survives
    /// from one dive to the next instead of being rebuilt per node. It is
    /// reset whenever the model matrix changes.
    session: ?solver.LpSolveSession = null,
    /// Matrix slices the retained session was last used with.
    matrix_identity: MatrixIdentity = .{},
    /// Whether the retained session already carries a solved state that
    /// supersedes any caller-provided basis.
    has_retained_state: bool = false,

    const MatrixIdentity = struct {
        col_starts: []const usize = &.{},
        row_indices: []const foundation.RowId = &.{},
        values: []const f64 = &.{},

        fn matches(self: MatrixIdentity, relaxation: solver.LpProblemView) bool {
            return self.col_starts.ptr == relaxation.matrix.col_starts.ptr and self.col_starts.len == relaxation.matrix.col_starts.len and
                self.row_indices.ptr == relaxation.matrix.row_indices.ptr and self.row_indices.len == relaxation.matrix.row_indices.len and
                self.values.ptr == relaxation.matrix.values.ptr and self.values.len == relaxation.matrix.values.len;
        }
    };

    /// Release the retained LP session. Must be called once after the final
    /// solve when the plugin instance outlives the search.
    pub fn deinit(self: *Diving) void {
        if (self.session) |*value| value.deinit();
        self.session = null;
        self.has_retained_state = false;
    }

    fn resetSession(self: *Diving, allocator: std.mem.Allocator) std.mem.Allocator.Error!void {
        self.deinit();
        self.session = solver.LpSolveSession.init(allocator);
        self.has_retained_state = false;
    }

    pub fn trySolution(self: *Diving, ctx: plugin.HeuristicContext) std.mem.Allocator.Error!?[]f64 {
        const relaxation = ctx.problem.relaxation;
        const num_cols = relaxation.num_cols;
        const node_lower = if (ctx.node_lower.len == num_cols) ctx.node_lower else relaxation.col_lower;
        const node_upper = if (ctx.node_upper.len == num_cols) ctx.node_upper else relaxation.col_upper;

        if (self.session == null or !self.matrix_identity.matches(relaxation))
            try self.resetSession(ctx.allocator);
        self.matrix_identity = .{
            .col_starts = relaxation.matrix.col_starts,
            .row_indices = relaxation.matrix.row_indices,
            .values = relaxation.matrix.values,
        };

        var lower = try ctx.allocator.dupe(f64, node_lower);
        defer ctx.allocator.free(lower);
        var upper = try ctx.allocator.dupe(f64, node_upper);
        defer ctx.allocator.free(upper);
        const current = try ctx.allocator.dupe(f64, ctx.lp_primal);
        defer ctx.allocator.free(current);

        var rounds: usize = 0;
        while (rounds < self.max_rounds) {
            var target: ?usize = null;
            var target_fractionality = std.math.inf(f64);
            for (ctx.problem.integrality, 0..) |kind, column| {
                switch (kind) {
                    .integer, .binary => {},
                    else => continue,
                }
                const frac = fractionality(current[column]);
                if (frac <= ctx.tolerance) continue;
                if (frac < target_fractionality) {
                    target = column;
                    target_fractionality = frac;
                }
            }
            const column = target orelse return try self.snap(ctx, current);
            rounds += 1;
            const fixed = std.math.clamp(@round(current[column]), lower[column], upper[column]);
            lower[column] = fixed;
            upper[column] = fixed;

            var sub_problem = relaxation;
            sub_problem.col_lower = lower;
            sub_problem.col_upper = upper;
            const status = self.session.?.solveRelaxation(sub_problem, .{
                .max_iterations = self.max_lp_iterations,
                .initial_basis = if (rounds == 1 and !self.has_retained_state) ctx.lp_basis else null,
            }, .{
                .structure = 0,
                .matrix_values = 0,
                .bounds = rounds,
                .objective = 0,
            });
            if (status == .numerical_failure) {
                // A numerically failed sub-LP leaves the engine in an
                // unreliable state; start the next dive from scratch.
                try self.resetSession(ctx.allocator);
                return null;
            }
            if (status != .optimal) return null;
            self.has_retained_state = true;
            const view = self.session.?.resultView(sub_problem, status) orelse return null;
            if (!view.primal_valid or view.primal.len != num_cols) return null;
            @memcpy(current, view.primal);
        }
        return try self.snap(ctx, current);
    }

    /// Snap every integer-like column of an (almost) integral point and return
    /// an owned candidate. Returns null when any integer/binary column is
    /// still fractional beyond tolerance — the dive simply failed.
    fn snap(self: *const Diving, ctx: plugin.HeuristicContext, current: []const f64) std.mem.Allocator.Error!?[]f64 {
        _ = self;
        const candidate = try ctx.allocator.dupe(f64, current);
        errdefer ctx.allocator.free(candidate);
        for (ctx.problem.integrality, candidate, 0..) |kind, *value, column| switch (kind) {
            .continuous => {},
            .integer => {
                if (fractionality(value.*) > ctx.tolerance) {
                    ctx.allocator.free(candidate);
                    return null;
                }
                value.* = std.math.clamp(@round(value.*), ctx.active_lower[column], ctx.active_upper[column]);
            },
            .binary => {
                if (fractionality(value.*) > ctx.tolerance) {
                    ctx.allocator.free(candidate);
                    return null;
                }
                value.* = std.math.clamp(@round(value.*), @max(0.0, ctx.problem.relaxation.col_lower[column]), @min(1.0, ctx.problem.relaxation.col_upper[column]));
            },
            .semi_continuous, .semi_integer => {
                const lower = ctx.active_lower[column];
                const upper = ctx.active_upper[column];
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
        return candidate;
    }
};

const matrix = @import("matrix");

fn contextFor(
    comptime integrality: []const plugin.Integrality,
    primal: []const f64,
    col_lower: []const f64,
    col_upper: []const f64,
    active_lower: []const f64,
    active_upper: []const f64,
) plugin.HeuristicContext {
    const num_cols = integrality.len;
    return .{
        .allocator = std.testing.allocator,
        .problem = .{
            .relaxation = .{
                .num_rows = 0,
                .num_cols = num_cols,
                .col_cost = &.{},
                .col_lower = col_lower,
                .col_upper = col_upper,
                .row_lower = &.{},
                .row_upper = &.{},
                .matrix = matrix.CscView.initAssumeValid(0, num_cols, &([_]usize{0} ** (num_cols + 1)), &.{}, &.{}),
                .objective_sense = .minimize,
                .objective_offset = 0.0,
            },
            .integrality = integrality,
        },
        .node_lower = &.{},
        .node_upper = &.{},
        .lp_primal = primal,
        .active_lower = active_lower,
        .active_upper = active_upper,
        .tolerance = 1e-6,
    };
}

test "rounding snaps integer columns to the nearest active integral value" {
    var heuristic = Rounding{};
    const ctx = contextFor(
        &.{ .integer, .integer },
        &.{ 0.4, 1.6 },
        &.{ 0.0, 0.0 },
        &.{ 5.0, 5.0 },
        &.{ 0.0, 0.0 },
        &.{ 5.0, 5.0 },
    );
    const candidate = try heuristic.trySolution(ctx);
    defer std.testing.allocator.free(candidate.?);
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 2.0 }, candidate.?);
}

test "rounding clamps to the active bounds when the rounded value escapes" {
    var heuristic = Rounding{};
    const ctx = contextFor(
        &.{.integer},
        &.{4.7},
        &.{0.0},
        &.{3.0},
        &.{0.0},
        &.{3.0},
    );
    const candidate = try heuristic.trySolution(ctx);
    defer std.testing.allocator.free(candidate.?);
    try std.testing.expectEqualSlices(f64, &.{3.0}, candidate.?);
}

test "rounding resolves semi-continuous columns to zero or the active value" {
    var heuristic = Rounding{};
    const ctx = contextFor(
        &.{ .semi_continuous, .semi_integer },
        &.{ 0.4, 4.6 },
        &.{ -1.0, -1.0 },
        &.{ 6.0, 6.0 },
        &.{ 2.0, 2.0 },
        &.{ 6.0, 6.0 },
    );
    const candidate = try heuristic.trySolution(ctx);
    defer std.testing.allocator.free(candidate.?);
    // 0.4 is closer to zero; 4.6 rounds to 5 inside the active range.
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 5.0 }, candidate.?);
}

test "rounding forces a semi column to zero when its active range is empty" {
    var heuristic = Rounding{};
    const ctx = contextFor(
        &.{.semi_continuous},
        &.{3.0},
        &.{
            -1.0,
        },
        &.{
            6.0,
        },
        &.{
            6.0,
        },
        &.{
            2.0,
        },
    );
    const candidate = try heuristic.trySolution(ctx);
    defer std.testing.allocator.free(candidate.?);
    try std.testing.expectEqualSlices(f64, &.{0.0}, candidate.?);
}

const foundation = @import("foundation");

fn divingContext(
    comptime integrality: []const plugin.Integrality,
    relaxation: solver.LpProblemView,
    primal: []const f64,
    lower: []const f64,
    upper: []const f64,
) plugin.HeuristicContext {
    return .{
        .allocator = std.testing.allocator,
        .problem = .{
            .relaxation = relaxation,
            .integrality = integrality,
        },
        .node_lower = lower,
        .node_upper = upper,
        .lp_primal = primal,
        .active_lower = lower,
        .active_upper = upper,
        .tolerance = 1e-6,
    };
}

// min x1 + x2 + x3 subject to x1 + x2 + x3 >= 2.5, all binary.
// LP relaxation: every xi = 5/6. Diving fixes each rounded column and
// re-solves until it reaches the incumbent (1, 1, 1) with objective 3.
fn divingFeasibleRelaxation() solver.LpProblemView {
    return .{
        .num_rows = 1,
        .num_cols = 3,
        .col_cost = &.{ 1.0, 1.0, 1.0 },
        .col_lower = &.{ 0.0, 0.0, 0.0 },
        .col_upper = &.{ 1.0, 1.0, 1.0 },
        .row_lower = &.{2.5},
        .row_upper = &.{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(
            1,
            3,
            &[_]usize{ 0, 1, 2, 3 },
            &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) },
            &.{ 1.0, 1.0, 1.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
}

test "diving re-solves the LP to reach an integer incumbent" {
    var heuristic = Diving{};
    defer heuristic.deinit();
    const relaxation = divingFeasibleRelaxation();
    const ctx = divingContext(
        &.{ .binary, .binary, .binary },
        relaxation,
        &.{ 5.0 / 6.0, 5.0 / 6.0, 5.0 / 6.0 },
        &.{ 0.0, 0.0, 0.0 },
        &.{ 1.0, 1.0, 1.0 },
    );
    const candidate = try heuristic.trySolution(ctx);
    defer std.testing.allocator.free(candidate.?);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0, 1.0 }, candidate.?);
}

test "diving borrows the parent basis without changing the node solution" {
    const relaxation = divingFeasibleRelaxation();
    var parent = solver.solveLp(std.testing.allocator, relaxation, .{});
    defer parent.deinit();
    try std.testing.expectEqual(solver.LpSolveStatus.optimal, parent.status);
    const view = parent.resultView(relaxation).?;
    const original_primal = try std.testing.allocator.dupe(f64, view.primal);
    defer std.testing.allocator.free(original_primal);
    var ctx = divingContext(&.{ .binary, .binary, .binary }, relaxation, view.primal, relaxation.col_lower, relaxation.col_upper);
    ctx.lp_basis = .{
        .structural_status = view.structural_status,
        .logical_status = view.logical_status,
        .basic_index = view.basic_index,
    };
    var heuristic = Diving{};
    defer heuristic.deinit();
    const candidate = (try heuristic.trySolution(ctx)).?;
    defer std.testing.allocator.free(candidate);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0, 1.0 }, candidate);
    try std.testing.expectEqualSlices(f64, original_primal, parent.resultView(relaxation).?.primal);
}

test "diving returns null when the dive hits an infeasible sub-LP" {
    var heuristic = Diving{};
    defer heuristic.deinit();
    // x + y == 1.5 with x,y binary has no integer feasible point; fixing the
    // rounded column drives the sub-LP infeasible.
    const relaxation = solver.LpProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 1.0, 1.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ 1.0, 1.0 },
        .row_lower = &.{1.5},
        .row_upper = &.{1.5},
        .matrix = matrix.CscView.initAssumeValid(
            1,
            2,
            &[_]usize{ 0, 1, 2 },
            &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) },
            &.{ 1.0, 1.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    const ctx = divingContext(
        &.{ .binary, .binary },
        relaxation,
        &.{ 0.75, 0.75 },
        &.{ 0.0, 0.0 },
        &.{ 1.0, 1.0 },
    );
    try std.testing.expectEqual(@as(?[]f64, null), try heuristic.trySolution(ctx));
}