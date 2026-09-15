//! Linear-row bound propagation for a reconstructed node domain.
//!
//! Tightenings are appended to the persistent domain-change arena so children
//! inherit them without copying full bound vectors. Semi-variable targets are
//! deliberately skipped until their disjunctive state is represented explicitly.

const std = @import("std");
const solver = @import("solver");
const types = @import("types.zig");
const domain = @import("domain.zig");

const Activity = struct {
    finite: f64 = 0.0,
    negative_infinite: usize = 0,
    positive_infinite: usize = 0,
    reliable: bool = true,

    fn add(self: *Activity, value: f64) void {
        if (value == -std.math.inf(f64)) {
            self.negative_infinite += 1;
        } else if (value == std.math.inf(f64)) {
            self.positive_infinite += 1;
        } else if (std.math.isFinite(value)) {
            const sum = self.finite + value;
            if (std.math.isFinite(sum)) {
                self.finite = sum;
            } else {
                self.reliable = false;
            }
        } else {
            self.reliable = false;
        }
    }

    fn excluding(self: Activity, value: f64) ?f64 {
        if (!self.reliable) return null;
        var negative = self.negative_infinite;
        var positive = self.positive_infinite;
        var finite = self.finite;
        if (value == -std.math.inf(f64)) {
            negative -= 1;
        } else if (value == std.math.inf(f64)) {
            positive -= 1;
        } else if (std.math.isFinite(value)) {
            finite -= value;
        } else {
            return null;
        }
        return if (negative == 0 and positive == 0 and std.math.isFinite(finite)) finite else null;
    }
};

pub const Result = struct {
    feasible: bool,
    change_head: ?domain.ChangeId,
    rounds: usize,
    tightenings: usize,
};

/// Solve-local scratch storage. Allocation cost is paid once, not per node.
pub const Workspace = struct {
    minimum: []Activity,
    maximum: []Activity,
    round_lower: []f64,
    round_upper: []f64,

    pub fn init(allocator: std.mem.Allocator, num_rows: usize, num_cols: usize) std.mem.Allocator.Error!Workspace {
        const minimum = try allocator.alloc(Activity, num_rows);
        errdefer allocator.free(minimum);
        const maximum = try allocator.alloc(Activity, num_rows);
        errdefer allocator.free(maximum);
        const round_lower = try allocator.alloc(f64, num_cols);
        errdefer allocator.free(round_lower);
        const round_upper = try allocator.alloc(f64, num_cols);
        return .{ .minimum = minimum, .maximum = maximum, .round_lower = round_lower, .round_upper = round_upper };
    }

    pub fn deinit(self: *Workspace, allocator: std.mem.Allocator) void {
        allocator.free(self.minimum);
        allocator.free(self.maximum);
        allocator.free(self.round_lower);
        allocator.free(self.round_upper);
        self.* = undefined;
    }
};

pub fn propagate(
    workspace: *Workspace,
    allocator: std.mem.Allocator,
    problem: solver.LpProblemView,
    integrality: []const types.Integrality,
    lower: []f64,
    upper: []f64,
    arena: *domain.DomainChangeArena,
    initial_head: ?domain.ChangeId,
    tolerance: f64,
    max_rounds: usize,
) std.mem.Allocator.Error!Result {
    if (max_rounds == 0 or problem.num_rows == 0 or problem.matrix.nnz() == 0)
        return .{ .feasible = true, .change_head = initial_head, .rounds = 0, .tightenings = 0 };

    std.debug.assert(workspace.minimum.len == problem.num_rows);
    std.debug.assert(workspace.maximum.len == problem.num_rows);
    std.debug.assert(workspace.round_lower.len == problem.num_cols);
    std.debug.assert(workspace.round_upper.len == problem.num_cols);
    const minimum = workspace.minimum;
    const maximum = workspace.maximum;
    const round_lower = workspace.round_lower;
    const round_upper = workspace.round_upper;

    var head = initial_head;
    var total_tightenings: usize = 0;
    var rounds: usize = 0;
    while (rounds < max_rounds) : (rounds += 1) {
        @memcpy(round_lower, lower);
        @memcpy(round_upper, upper);
        @memset(minimum, .{});
        @memset(maximum, .{});
        for (0..problem.num_cols) |column| {
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            for (begin..end) |position| {
                const row = problem.matrix.row_indices[position].toUsize();
                const coefficient = problem.matrix.values[position];
                minimum[row].add(coefficient * if (coefficient > 0.0) round_lower[column] else round_upper[column]);
                maximum[row].add(coefficient * if (coefficient > 0.0) round_upper[column] else round_lower[column]);
            }
        }

        for (0..problem.num_rows) |row| {
            const row_upper = problem.row_upper[row];
            if (minimum[row].reliable and std.math.isFinite(row_upper) and minimum[row].negative_infinite == 0 and
                (minimum[row].positive_infinite != 0 or minimum[row].finite > row_upper + tolerance))
                return .{ .feasible = false, .change_head = head, .rounds = rounds + 1, .tightenings = total_tightenings };
            const row_lower = problem.row_lower[row];
            if (maximum[row].reliable and std.math.isFinite(row_lower) and maximum[row].positive_infinite == 0 and
                (maximum[row].negative_infinite != 0 or maximum[row].finite < row_lower - tolerance))
                return .{ .feasible = false, .change_head = head, .rounds = rounds + 1, .tightenings = total_tightenings };
        }

        var changed = false;
        for (0..problem.num_cols) |column| {
            const kind = integrality[column];
            if (kind == .semi_continuous or kind == .semi_integer) continue;
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            for (begin..end) |position| {
                const row = problem.matrix.row_indices[position].toUsize();
                const coefficient = problem.matrix.values[position];
                const min_term = coefficient * if (coefficient > 0.0) round_lower[column] else round_upper[column];
                const max_term = coefficient * if (coefficient > 0.0) round_upper[column] else round_lower[column];

                if (std.math.isFinite(problem.row_upper[row])) {
                    if (minimum[row].excluding(min_term)) |other| {
                        const implied = (problem.row_upper[row] - other) / coefficient;
                        if (coefficient > 0.0) {
                            const candidate = roundedUpper(kind, implied, tolerance);
                            const outcome = try tightenUpper(allocator, arena, &head, column, lower, upper, candidate, row, tolerance);
                            if (!outcome.feasible) return .{ .feasible = false, .change_head = head, .rounds = rounds + 1, .tightenings = total_tightenings };
                            changed = changed or outcome.changed;
                            total_tightenings += @intFromBool(outcome.changed);
                        } else {
                            const candidate = roundedLower(kind, implied, tolerance);
                            const outcome = try tightenLower(allocator, arena, &head, column, lower, upper, candidate, row, tolerance);
                            if (!outcome.feasible) return .{ .feasible = false, .change_head = head, .rounds = rounds + 1, .tightenings = total_tightenings };
                            changed = changed or outcome.changed;
                            total_tightenings += @intFromBool(outcome.changed);
                        }
                    }
                }
                if (std.math.isFinite(problem.row_lower[row])) {
                    if (maximum[row].excluding(max_term)) |other| {
                        const implied = (problem.row_lower[row] - other) / coefficient;
                        if (coefficient > 0.0) {
                            const candidate = roundedLower(kind, implied, tolerance);
                            const outcome = try tightenLower(allocator, arena, &head, column, lower, upper, candidate, row, tolerance);
                            if (!outcome.feasible) return .{ .feasible = false, .change_head = head, .rounds = rounds + 1, .tightenings = total_tightenings };
                            changed = changed or outcome.changed;
                            total_tightenings += @intFromBool(outcome.changed);
                        } else {
                            const candidate = roundedUpper(kind, implied, tolerance);
                            const outcome = try tightenUpper(allocator, arena, &head, column, lower, upper, candidate, row, tolerance);
                            if (!outcome.feasible) return .{ .feasible = false, .change_head = head, .rounds = rounds + 1, .tightenings = total_tightenings };
                            changed = changed or outcome.changed;
                            total_tightenings += @intFromBool(outcome.changed);
                        }
                    }
                }
            }
        }
        if (!changed) return .{ .feasible = true, .change_head = head, .rounds = rounds + 1, .tightenings = total_tightenings };
    }
    return .{ .feasible = true, .change_head = head, .rounds = rounds, .tightenings = total_tightenings };
}

const TightenResult = struct { feasible: bool = true, changed: bool = false };

fn tightenLower(allocator: std.mem.Allocator, arena: *domain.DomainChangeArena, head: *?domain.ChangeId, column: usize, lower: []f64, upper: []const f64, raw_candidate: f64, row: usize, tolerance: f64) std.mem.Allocator.Error!TightenResult {
    if (!std.math.isFinite(raw_candidate) or raw_candidate <= lower[column] + tolerance) return .{};
    var candidate = raw_candidate;
    if (candidate > upper[column]) {
        if (candidate > upper[column] + tolerance) return .{ .feasible = false };
        candidate = upper[column];
    }
    lower[column] = candidate;
    head.* = try arena.appendWithReason(allocator, head.*, .{ .column = column, .lower = candidate }, .{ .row = row });
    return .{ .changed = true };
}

fn tightenUpper(allocator: std.mem.Allocator, arena: *domain.DomainChangeArena, head: *?domain.ChangeId, column: usize, lower: []const f64, upper: []f64, raw_candidate: f64, row: usize, tolerance: f64) std.mem.Allocator.Error!TightenResult {
    if (!std.math.isFinite(raw_candidate) or raw_candidate >= upper[column] - tolerance) return .{};
    var candidate = raw_candidate;
    if (candidate < lower[column]) {
        if (candidate < lower[column] - tolerance) return .{ .feasible = false };
        candidate = lower[column];
    }
    upper[column] = candidate;
    head.* = try arena.appendWithReason(allocator, head.*, .{ .column = column, .upper = candidate }, .{ .row = row });
    return .{ .changed = true };
}

fn roundedLower(kind: types.Integrality, value: f64, tolerance: f64) f64 {
    return switch (kind) {
        .binary, .integer => @ceil(value - tolerance),
        else => value,
    };
}

fn roundedUpper(kind: types.Integrality, value: f64, tolerance: f64) f64 {
    return switch (kind) {
        .binary, .integer => @floor(value + tolerance),
        else => value,
    };
}

test "row propagation rounds integer bounds and persists tightenings" {
    const foundation = @import("foundation");
    const matrix = @import("matrix");
    const problem = solver.LpProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 0.0, 0.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ 1.0, 1.0 },
        .row_lower = &.{1.5},
        .row_upper = &.{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &.{ 0, 1, 2 }, &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) }, &.{ 1.0, 1.0 }),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var lower = [_]f64{ 0.0, 0.0 };
    var upper = [_]f64{ 1.0, 1.0 };
    var arena = domain.DomainChangeArena{};
    defer arena.deinit(std.testing.allocator);
    var workspace = try Workspace.init(std.testing.allocator, problem.num_rows, problem.num_cols);
    defer workspace.deinit(std.testing.allocator);
    const result = try propagate(&workspace, std.testing.allocator, problem, &.{ .binary, .binary }, &lower, &upper, &arena, null, 1e-6, 8);
    try std.testing.expect(result.feasible);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0 }, &lower);
    try std.testing.expectEqual(@as(usize, 2), result.tightenings);
    try std.testing.expect(result.change_head != null);
    var replay = (try domain.Domain.init(std.testing.allocator, &.{ .binary, .binary }, problem.col_lower, problem.col_upper, 1e-6)).?;
    defer replay.deinit(std.testing.allocator);
    try std.testing.expect(replay.load(&arena, result.change_head));
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0 }, replay.lower);
}

test "row propagation detects infeasible domains without an LP solve" {
    const foundation = @import("foundation");
    const matrix = @import("matrix");
    const problem = solver.LpProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &.{0.0},
        .col_lower = &.{0.0},
        .col_upper = &.{1.0},
        .row_lower = &.{2.0},
        .row_upper = &.{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &.{ 0, 1 }, &.{foundation.RowId.fromUsizeAssumeValid(0)}, &.{1.0}),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var lower = [_]f64{0.0};
    var upper = [_]f64{1.0};
    var arena = domain.DomainChangeArena{};
    defer arena.deinit(std.testing.allocator);
    var workspace = try Workspace.init(std.testing.allocator, problem.num_rows, problem.num_cols);
    defer workspace.deinit(std.testing.allocator);
    const result = try propagate(&workspace, std.testing.allocator, problem, &.{.integer}, &lower, &upper, &arena, null, 1e-6, 8);
    try std.testing.expect(!result.feasible);
}

test "overflowed row activities are ignored conservatively" {
    const foundation = @import("foundation");
    const matrix = @import("matrix");
    const problem = solver.LpProblemView{
        .num_rows = 1,
        .num_cols = 3,
        .col_cost = &.{ 0.0, 0.0, 0.0 },
        .col_lower = &.{ 1.0, 1.0, 1.0 },
        .col_upper = &.{ 1.0, 1.0, 1.0 },
        .row_lower = &.{-std.math.inf(f64)},
        .row_upper = &.{1.1e308},
        .matrix = matrix.CscView.initAssumeValid(1, 3, &.{ 0, 1, 2, 3 }, &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) }, &.{ 1.0e308, 1.0e308, -1.0e308 }),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var lower = [_]f64{ 1.0, 1.0, 1.0 };
    var upper = [_]f64{ 1.0, 1.0, 1.0 };
    var arena = domain.DomainChangeArena{};
    defer arena.deinit(std.testing.allocator);
    var workspace = try Workspace.init(std.testing.allocator, problem.num_rows, problem.num_cols);
    defer workspace.deinit(std.testing.allocator);
    const result = try propagate(&workspace, std.testing.allocator, problem, &.{ .continuous, .continuous, .continuous }, &lower, &upper, &arena, null, 1e-6, 8);
    try std.testing.expect(result.feasible);
    try std.testing.expectEqual(@as(usize, 0), result.tightenings);
}
