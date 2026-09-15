//! Sparse binary-literal conflict graph extracted from globally valid rows.
//!
//! A graph edge means that its two literals cannot both be true. Both `x` and
//! `~x` are first-class vertices, so propagation can express implications in
//! either direction and follow them transitively without allocating per node.

const std = @import("std");
const solver = @import("solver");
const types = @import("types.zig");
const domain = @import("domain.zig");

pub const LiteralId = enum(usize) {
    _,

    pub fn init(column_index: usize, literal_value: bool) LiteralId {
        return @enumFromInt(2 * column_index + @intFromBool(literal_value));
    }
    pub fn column(self: LiteralId) usize {
        return @intFromEnum(self) / 2;
    }
    pub fn value(self: LiteralId) bool {
        return @intFromEnum(self) % 2 != 0;
    }
    pub fn complement(self: LiteralId) LiteralId {
        return @enumFromInt(@intFromEnum(self) ^ 1);
    }
};

const Pair = struct { first: LiteralId, second: LiteralId };
const Candidate = struct { literal: LiteralId, delta: f64 };

pub const PropagationResult = struct {
    feasible: bool,
    change_head: ?domain.ChangeId,
    tightenings: usize,
};

pub const CliqueTable = struct {
    allocator: std.mem.Allocator,
    starts: []usize,
    neighbors: []LiteralId,
    queue: []LiteralId,
    seen_generation: []u32,
    generation: u32 = 0,

    pub fn init(
        allocator: std.mem.Allocator,
        problem: solver.LpProblemView,
        integrality: []const types.Integrality,
        root_lower: []const f64,
        root_upper: []const f64,
        tolerance: f64,
        max_row_entries: usize,
        max_edges: usize,
    ) std.mem.Allocator.Error!CliqueTable {
        std.debug.assert(integrality.len == problem.num_cols);
        std.debug.assert(root_lower.len == problem.num_cols);
        std.debug.assert(root_upper.len == problem.num_cols);
        const literal_count = std.math.mul(usize, problem.num_cols, 2) catch return error.OutOfMemory;

        const row_starts = try allocator.alloc(usize, problem.num_rows + 1);
        defer allocator.free(row_starts);
        @memset(row_starts, 0);
        for (problem.matrix.row_indices) |row| row_starts[row.toUsize() + 1] += 1;
        for (1..row_starts.len) |index| row_starts[index] += row_starts[index - 1];
        const cursor = try allocator.dupe(usize, row_starts[0..problem.num_rows]);
        defer allocator.free(cursor);
        const row_columns = try allocator.alloc(usize, problem.matrix.nnz());
        defer allocator.free(row_columns);
        const row_values = try allocator.alloc(f64, problem.matrix.nnz());
        defer allocator.free(row_values);
        for (0..problem.num_cols) |column_index| {
            for (problem.matrix.col_starts[column_index]..problem.matrix.col_starts[column_index + 1]) |position| {
                const row = problem.matrix.row_indices[position].toUsize();
                const target = cursor[row];
                cursor[row] += 1;
                row_columns[target] = column_index;
                row_values[target] = problem.matrix.values[position];
            }
        }

        var pairs = std.ArrayListUnmanaged(Pair).empty;
        defer pairs.deinit(allocator);
        if (max_row_entries != 0 and max_edges != 0) row_loop: for (0..problem.num_rows) |row| {
            const begin = row_starts[row];
            const end = row_starts[row + 1];
            if (end - begin < 2 or end - begin > max_row_entries) continue;
            const columns = row_columns[begin..end];
            const values = row_values[begin..end];
            if (std.math.isFinite(problem.row_upper[row]))
                if (!try extractSide(allocator, &pairs, columns, values, integrality, root_lower, root_upper, problem.row_upper[row], 1.0, tolerance, max_edges)) break :row_loop;
            if (std.math.isFinite(problem.row_lower[row]))
                if (!try extractSide(allocator, &pairs, columns, values, integrality, root_lower, root_upper, -problem.row_lower[row], -1.0, tolerance, max_edges)) break :row_loop;
        };

        std.mem.sort(Pair, pairs.items, {}, pairBefore);
        var unique_count: usize = 0;
        for (pairs.items) |pair| {
            if (unique_count != 0 and samePair(pairs.items[unique_count - 1], pair)) continue;
            pairs.items[unique_count] = pair;
            unique_count += 1;
        }
        pairs.items.len = unique_count;

        const starts = try allocator.alloc(usize, literal_count + 1);
        errdefer allocator.free(starts);
        @memset(starts, 0);
        for (pairs.items) |pair| {
            starts[@intFromEnum(pair.first) + 1] += 1;
            starts[@intFromEnum(pair.second) + 1] += 1;
        }
        for (1..starts.len) |index| starts[index] += starts[index - 1];
        const neighbors = try allocator.alloc(LiteralId, pairs.items.len * 2);
        errdefer allocator.free(neighbors);
        const offsets = try allocator.dupe(usize, starts[0..literal_count]);
        defer allocator.free(offsets);
        for (pairs.items) |pair| {
            const first = @intFromEnum(pair.first);
            const second = @intFromEnum(pair.second);
            neighbors[offsets[first]] = pair.second;
            offsets[first] += 1;
            neighbors[offsets[second]] = pair.first;
            offsets[second] += 1;
        }
        const queue = try allocator.alloc(LiteralId, problem.num_cols);
        errdefer allocator.free(queue);
        const seen_generation = try allocator.alloc(u32, literal_count);
        errdefer allocator.free(seen_generation);
        @memset(seen_generation, 0);
        return .{ .allocator = allocator, .starts = starts, .neighbors = neighbors, .queue = queue, .seen_generation = seen_generation };
    }

    pub fn deinit(self: *CliqueTable) void {
        self.allocator.free(self.starts);
        self.allocator.free(self.neighbors);
        self.allocator.free(self.queue);
        self.allocator.free(self.seen_generation);
        self.* = undefined;
    }

    pub fn edgeCount(self: *const CliqueTable) usize {
        return self.neighbors.len / 2;
    }

    pub fn propagate(
        self: *CliqueTable,
        allocator: std.mem.Allocator,
        lower: []f64,
        upper: []f64,
        arena: *domain.DomainChangeArena,
        initial_head: ?domain.ChangeId,
        tolerance: f64,
    ) std.mem.Allocator.Error!PropagationResult {
        std.debug.assert(self.starts.len == lower.len * 2 + 1);
        std.debug.assert(upper.len == lower.len);
        self.nextGeneration();
        var queue_end: usize = 0;
        for (lower, upper, 0..) |lo, up, column_index| {
            const zero_true = up <= tolerance;
            const one_true = lo >= 1.0 - tolerance;
            if (zero_true and one_true) return .{ .feasible = false, .change_head = initial_head, .tightenings = 0 };
            if (zero_true or one_true) self.enqueue(LiteralId.init(column_index, one_true), &queue_end);
        }

        var head = initial_head;
        var tightenings: usize = 0;
        var queue_head: usize = 0;
        while (queue_head < queue_end) : (queue_head += 1) {
            const literal = self.queue[queue_head];
            const index = @intFromEnum(literal);
            for (self.neighbors[self.starts[index]..self.starts[index + 1]]) |conflict| {
                if (literalTrue(conflict, lower, upper, tolerance))
                    return .{ .feasible = false, .change_head = head, .tightenings = tightenings };
                if (literalFalse(conflict, lower, upper, tolerance)) continue;
                const column_index = conflict.column();
                if (conflict.value()) {
                    upper[column_index] = 0.0;
                    head = try arena.appendWithReason(allocator, head, .{ .column = column_index, .upper = 0.0 }, .clique);
                } else {
                    lower[column_index] = 1.0;
                    head = try arena.appendWithReason(allocator, head, .{ .column = column_index, .lower = 1.0 }, .clique);
                }
                tightenings += 1;
                self.enqueue(conflict.complement(), &queue_end);
            }
        }
        return .{ .feasible = true, .change_head = head, .tightenings = tightenings };
    }

    fn nextGeneration(self: *CliqueTable) void {
        self.generation +%= 1;
        if (self.generation == 0) {
            @memset(self.seen_generation, 0);
            self.generation = 1;
        }
    }

    fn enqueue(self: *CliqueTable, literal: LiteralId, queue_end: *usize) void {
        const index = @intFromEnum(literal);
        if (self.seen_generation[index] == self.generation) return;
        std.debug.assert(queue_end.* < self.queue.len);
        self.seen_generation[index] = self.generation;
        self.queue[queue_end.*] = literal;
        queue_end.* += 1;
    }
};

fn extractSide(
    allocator: std.mem.Allocator,
    pairs: *std.ArrayListUnmanaged(Pair),
    columns: []const usize,
    values: []const f64,
    integrality: []const types.Integrality,
    root_lower: []const f64,
    root_upper: []const f64,
    rhs: f64,
    multiplier: f64,
    tolerance: f64,
    max_edges: usize,
) std.mem.Allocator.Error!bool {
    var minimum: f64 = 0.0;
    for (columns, values) |column_index, raw_coefficient| {
        const coefficient = multiplier * raw_coefficient;
        if (coefficient == 0.0) continue;
        const bound = if (coefficient > 0.0) root_lower[column_index] else root_upper[column_index];
        const term = coefficient * bound;
        const next = minimum + term;
        if (!std.math.isFinite(term) or !std.math.isFinite(next)) return true;
        minimum = next;
    }
    for (columns, values, 0..) |left_column, left_raw, left_position| {
        const left = candidate(left_column, multiplier * left_raw, integrality, root_lower, root_upper, tolerance) orelse continue;
        for (left_position + 1..columns.len) |right_position| {
            const right_column = columns[right_position];
            if (left_column == right_column) continue;
            const right = candidate(right_column, multiplier * values[right_position], integrality, root_lower, root_upper, tolerance) orelse continue;
            if (minimum + left.delta + right.delta <= rhs + tolerance) continue;
            if (pairs.items.len >= max_edges) return false;
            try pairs.append(allocator, normalizedPair(left.literal, right.literal));
        }
    }
    return true;
}

fn candidate(column_index: usize, coefficient: f64, integrality: []const types.Integrality, lower: []const f64, upper: []const f64, tolerance: f64) ?Candidate {
    if (integrality[column_index] != .binary) return null;
    if (coefficient > 0.0 and lower[column_index] < 1.0 - tolerance and upper[column_index] >= 1.0 - tolerance)
        return .{ .literal = LiteralId.init(column_index, true), .delta = coefficient * (1.0 - lower[column_index]) };
    if (coefficient < 0.0 and lower[column_index] <= tolerance and upper[column_index] > tolerance)
        return .{ .literal = LiteralId.init(column_index, false), .delta = -coefficient * upper[column_index] };
    return null;
}

fn literalTrue(literal: LiteralId, lower: []const f64, upper: []const f64, tolerance: f64) bool {
    return if (literal.value()) lower[literal.column()] >= 1.0 - tolerance else upper[literal.column()] <= tolerance;
}
fn literalFalse(literal: LiteralId, lower: []const f64, upper: []const f64, tolerance: f64) bool {
    return if (literal.value()) upper[literal.column()] <= tolerance else lower[literal.column()] >= 1.0 - tolerance;
}
fn normalizedPair(first: LiteralId, second: LiteralId) Pair {
    return if (@intFromEnum(first) < @intFromEnum(second)) .{ .first = first, .second = second } else .{ .first = second, .second = first };
}
fn samePair(lhs: Pair, rhs: Pair) bool {
    return lhs.first == rhs.first and lhs.second == rhs.second;
}
fn pairBefore(_: void, lhs: Pair, rhs: Pair) bool {
    return @intFromEnum(lhs.first) < @intFromEnum(rhs.first) or (lhs.first == rhs.first and @intFromEnum(lhs.second) < @intFromEnum(rhs.second));
}

test "literal graph deduplicates positive conflicts and propagates" {
    const foundation = @import("foundation");
    const matrix = @import("matrix");
    const problem = solver.LpProblemView{
        .num_rows = 2,
        .num_cols = 3,
        .col_cost = &.{ 0.0, 0.0, 0.0 },
        .col_lower = &.{ 0.0, 0.0, 0.0 },
        .col_upper = &.{ 1.0, 1.0, 1.0 },
        .row_lower = &.{ -std.math.inf(f64), -std.math.inf(f64) },
        .row_upper = &.{ 1.0, 1.0 },
        .matrix = matrix.CscView.initAssumeValid(2, 3, &.{ 0, 2, 4, 4 }, &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1) }, &.{ 1.0, 1.0, 1.0, 1.0 }),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var table = try CliqueTable.init(std.testing.allocator, problem, &.{ .binary, .binary, .binary }, problem.col_lower, problem.col_upper, 1e-6, 16, 100);
    defer table.deinit();
    try std.testing.expectEqual(@as(usize, 1), table.edgeCount());
    var lower = [_]f64{ 1.0, 0.0, 0.0 };
    var upper = [_]f64{ 1.0, 1.0, 1.0 };
    var arena = domain.DomainChangeArena{};
    defer arena.deinit(std.testing.allocator);
    const result = try table.propagate(std.testing.allocator, &lower, &upper, &arena, null, 1e-6);
    try std.testing.expect(result.feasible);
    try std.testing.expectEqual(@as(usize, 1), result.tightenings);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 0.0, 1.0 }, &upper);
}

test "complemented implications propagate transitively" {
    const foundation = @import("foundation");
    const matrix = @import("matrix");
    const problem = solver.LpProblemView{
        .num_rows = 2,
        .num_cols = 3,
        .col_cost = &.{ 0.0, 0.0, 0.0 },
        .col_lower = &.{ 0.0, 0.0, 0.0 },
        .col_upper = &.{ 1.0, 1.0, 1.0 },
        .row_lower = &.{ -std.math.inf(f64), -std.math.inf(f64) },
        .row_upper = &.{ 0.0, 1.0 },
        .matrix = matrix.CscView.initAssumeValid(2, 3, &.{ 0, 1, 3, 4 }, &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(1) }, &.{ 1.0, -1.0, 1.0, 1.0 }),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var table = try CliqueTable.init(std.testing.allocator, problem, &.{ .binary, .binary, .binary }, problem.col_lower, problem.col_upper, 1e-6, 16, 100);
    defer table.deinit();
    try std.testing.expectEqual(@as(usize, 2), table.edgeCount());
    var lower = [_]f64{ 1.0, 0.0, 0.0 };
    var upper = [_]f64{ 1.0, 1.0, 1.0 };
    var arena = domain.DomainChangeArena{};
    defer arena.deinit(std.testing.allocator);
    const result = try table.propagate(std.testing.allocator, &lower, &upper, &arena, null, 1e-6);
    try std.testing.expect(result.feasible);
    try std.testing.expectEqual(@as(usize, 2), result.tightenings);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0, 0.0 }, &lower);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0, 0.0 }, &upper);
}

test "lower row side creates a complemented conflict" {
    const foundation = @import("foundation");
    const matrix = @import("matrix");
    const problem = solver.LpProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 0.0, 0.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ 1.0, 1.0 },
        .row_lower = &.{0.0},
        .row_upper = &.{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &.{ 0, 1, 2 }, &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) }, &.{ 1.0, -1.0 }),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var table = try CliqueTable.init(std.testing.allocator, problem, &.{ .binary, .binary }, problem.col_lower, problem.col_upper, 1e-6, 16, 100);
    defer table.deinit();
    var lower = [_]f64{ 0.0, 1.0 };
    var upper = [_]f64{ 1.0, 1.0 };
    var arena = domain.DomainChangeArena{};
    defer arena.deinit(std.testing.allocator);
    const result = try table.propagate(std.testing.allocator, &lower, &upper, &arena, null, 1e-6);
    try std.testing.expect(result.feasible);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0 }, &lower);
}
