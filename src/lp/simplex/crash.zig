//! Allocation-reusing sparse crash-basis planner.
//!
//! The planner operates only on borrowed CSC data. It peels singleton and
//! near-singleton structural columns, producing a row/column matching that is
//! installed transactionally by the simplex engine and validated by a fresh
//! numerical factorization before publication.

const std = @import("std");
const matrix = @import("matrix");

/// Output of a crash plan: parallel arrays of (row, column) pairs that the
/// engine should install as basic structural variables.
pub const CrashPlanView = struct {
    rows: []const u32,
    columns: []const u32,
};

/// Selection policy for ranking crash candidates.
/// - `ltssf`: Lexicographic tie-break: degree, bound, cost, row degree, pivot magnitude.
/// - `bixby`: Merit-based score that rewards strong pivots in sparse rows.
pub const CrashScoring = enum { ltssf, bixby };

pub const CrashWorkspace = struct {
    allocator: std.mem.Allocator,
    row_starts: []usize = &.{}, // CSR row pointers (length = num_rows + 1)
    row_columns: []u32 = &.{}, // Column index per retained nonzero
    row_cursor: []usize = &.{}, // Scatter cursor used while building CSR
    row_degree: []u32 = &.{}, // Live nonzero count per row
    column_degree: []u32 = &.{}, // Live nonzero count per column
    row_active: []bool = &.{}, // False once a row has been matched
    column_active: []bool = &.{}, // False once a column has been matched
    selected_rows: []u32 = &.{}, // Output: matched row for each selection
    selected_columns: []u32 = &.{}, // Output: matched column for each selection

    pub fn init(allocator: std.mem.Allocator) CrashWorkspace {
        return .{ .allocator = allocator };
    }

    /// Free all allocated buffers.
    pub fn deinit(self: *CrashWorkspace) void {
        self.allocator.free(self.row_starts);
        self.allocator.free(self.row_columns);
        self.allocator.free(self.row_cursor);
        self.allocator.free(self.row_degree);
        self.allocator.free(self.column_degree);
        self.allocator.free(self.row_active);
        self.allocator.free(self.column_active);
        self.allocator.free(self.selected_rows);
        self.allocator.free(self.selected_columns);
        self.* = .{ .allocator = self.allocator };
    }

    /// Build a deterministic sparse structural matching with the requested
    /// scoring policy. No allocation occurs after workspace capacity growth.
    pub fn plan(
        self: *CrashWorkspace,
        csc: matrix.CscView,
        cost: []const f64,
        lower: []const f64,
        upper: []const f64,
        row_scale: []const f64,
        column_scale: []const f64,
        near_singleton_limit: u32, // Max acceptable column degree for crash eligibility
        pivot_tolerance: f64, // Minimum scaled pivot magnitude
        scoring: CrashScoring,
    ) !CrashPlanView {
        if (cost.len != csc.num_cols or lower.len != csc.num_cols or upper.len != csc.num_cols or
            row_scale.len != csc.num_rows or column_scale.len != csc.num_cols)
            return error.DimensionMismatch;
        try self.ensureCapacity(csc.num_rows, csc.num_cols, csc.values.len);
        const row_starts = self.row_starts[0 .. csc.num_rows + 1];
        const row_cursor = self.row_cursor[0..csc.num_rows];
        const row_degree = self.row_degree[0..csc.num_rows];
        const column_degree = self.column_degree[0..csc.num_cols];
        const row_active = self.row_active[0..csc.num_rows];
        const column_active = self.column_active[0..csc.num_cols];
        @memset(row_starts, 0);
        @memset(column_degree, 0);
        @memset(row_active, true);
        @memset(column_active, true);

        // First pass: count retained nonzeros per row and per column.
        for (0..csc.num_cols) |column| {
            const begin = csc.col_starts[column];
            const end = csc.col_starts[column + 1];
            for (csc.row_indices[begin..end], csc.values[begin..end]) |row_id, coefficient| {
                if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                const row = row_id.toUsize();
                row_starts[row + 1] += 1;
                column_degree[column] += 1;
            }
        }
        // Prefix-sum per-row counts into CSR row pointers and copy them into the cursor.
        for (0..csc.num_rows) |row| row_starts[row + 1] += row_starts[row];
        @memcpy(row_cursor, row_starts[0..csc.num_rows]);
        // Second pass: scatter column indices into row-major order.
        for (0..csc.num_cols) |column| {
            const begin = csc.col_starts[column];
            const end = csc.col_starts[column + 1];
            for (csc.row_indices[begin..end], csc.values[begin..end]) |row_id, coefficient| {
                if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                const row = row_id.toUsize();
                self.row_columns[row_cursor[row]] = @intCast(column);
                row_cursor[row] += 1;
            }
        }
        for (0..csc.num_rows) |row| row_degree[row] = @intCast(row_starts[row + 1] - row_starts[row]);

        // Greedy matching loop: repeatedly pick the best eligible (column, row)
        // pair, mark both sides inactive, and decrement degrees of neighbors.
        var selected_count: usize = 0;
        while (selected_count < @min(csc.num_rows, csc.num_cols)) {
            var best: ?Candidate = null;
            for (0..csc.num_cols) |column| {
                if (!column_active[column]) continue;
                const degree = column_degree[column];
                if (degree == 0 or degree > near_singleton_limit) continue;
                const begin = csc.col_starts[column];
                const end = csc.col_starts[column + 1];
                // Find the strongest eligible pivot in this column (ties broken
                // by smaller row degree to keep fill-in low).
                var pivot_row: ?usize = null;
                var pivot_magnitude: f64 = 0.0;
                for (csc.row_indices[begin..end], csc.values[begin..end]) |row_id, coefficient| {
                    const row = row_id.toUsize();
                    if (!row_active[row] or !matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                    const magnitude = @abs(coefficient * row_scale[row] * column_scale[column]);
                    if (magnitude > pivot_magnitude or (magnitude == pivot_magnitude and
                        (pivot_row == null or row_degree[row] < row_degree[pivot_row.?])))
                    {
                        pivot_row = row;
                        pivot_magnitude = magnitude;
                    }
                }
                const row = pivot_row orelse continue;
                if (!std.math.isFinite(pivot_magnitude) or pivot_magnitude <= pivot_tolerance) continue;
                const candidate = Candidate{
                    .column = column,
                    .row = row,
                    .degree = degree,
                    .bound_penalty = boundPenalty(lower[column], upper[column]),
                    .cost = @abs(cost[column] * column_scale[column]),
                    .row_degree = row_degree[row],
                    .pivot_magnitude = pivot_magnitude,
                };
                if (best == null or candidate.betterThan(best.?, scoring)) best = candidate;
            }
            const chosen = best orelse break;
            self.selected_rows[selected_count] = @intCast(chosen.row);
            self.selected_columns[selected_count] = @intCast(chosen.column);
            selected_count += 1;
            row_active[chosen.row] = false;
            column_active[chosen.column] = false;

            // Decrement degree of every column touched by the matched row.
            for (self.row_columns[row_starts[chosen.row]..row_starts[chosen.row + 1]]) |column_u32| {
                const column: usize = @intCast(column_u32);
                if (column_active[column] and column_degree[column] != 0) column_degree[column] -= 1;
            }
            // Decrement degree of every row touched by the matched column.
            const begin = csc.col_starts[chosen.column];
            const end = csc.col_starts[chosen.column + 1];
            for (csc.row_indices[begin..end], csc.values[begin..end]) |row_id, coefficient| {
                const row = row_id.toUsize();
                if (row_active[row] and matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient) and row_degree[row] != 0)
                    row_degree[row] -= 1;
            }
        }
        return .{
            .rows = self.selected_rows[0..selected_count],
            .columns = self.selected_columns[0..selected_count],
        };
    }

    /// Total bytes currently held by all dynamic buffers.
    pub fn requestedBytes(self: *const CrashWorkspace) usize {
        return std.mem.sliceAsBytes(self.row_starts).len +
            std.mem.sliceAsBytes(self.row_columns).len +
            std.mem.sliceAsBytes(self.row_cursor).len +
            std.mem.sliceAsBytes(self.row_degree).len +
            std.mem.sliceAsBytes(self.column_degree).len +
            std.mem.sliceAsBytes(self.row_active).len +
            std.mem.sliceAsBytes(self.column_active).len +
            std.mem.sliceAsBytes(self.selected_rows).len +
            std.mem.sliceAsBytes(self.selected_columns).len;
    }

    /// Grow buffers as needed; only ever enlarges, never shrinks.
    fn ensureCapacity(self: *CrashWorkspace, rows: usize, columns: usize, nonzeros: usize) !void {
        if (self.row_starts.len < rows + 1) self.row_starts = try self.allocator.realloc(self.row_starts, rows + 1);
        if (self.row_columns.len < nonzeros) self.row_columns = try self.allocator.realloc(self.row_columns, nonzeros);
        if (self.row_cursor.len < rows) self.row_cursor = try self.allocator.realloc(self.row_cursor, rows);
        if (self.row_degree.len < rows) self.row_degree = try self.allocator.realloc(self.row_degree, rows);
        if (self.column_degree.len < columns) self.column_degree = try self.allocator.realloc(self.column_degree, columns);
        if (self.row_active.len < rows) self.row_active = try self.allocator.realloc(self.row_active, rows);
        if (self.column_active.len < columns) self.column_active = try self.allocator.realloc(self.column_active, columns);
        const selected = @min(rows, columns);
        if (self.selected_rows.len < selected) self.selected_rows = try self.allocator.realloc(self.selected_rows, selected);
        if (self.selected_columns.len < selected) self.selected_columns = try self.allocator.realloc(self.selected_columns, selected);
    }

    /// A single (column, row) crash candidate with the data needed for ranking.
    const Candidate = struct {
        column: usize,
        row: usize,
        degree: u32, // Column degree at selection time
        bound_penalty: u8, // 0 (free), 1 (boxed), 2 (fixed)
        cost: f64, // Scaled |objective coefficient|
        row_degree: u32,
        pivot_magnitude: f64, // Scaled |coefficient|

        /// Return true if `self` should be selected before `other` under `scoring`.
        fn betterThan(self: Candidate, other: Candidate, scoring: CrashScoring) bool {
            if (scoring == .ltssf) {
                // Lexicographic order on (degree, bound_penalty, cost, row_degree, -pivot).
                if (self.degree != other.degree) return self.degree < other.degree;
                if (self.bound_penalty != other.bound_penalty) return self.bound_penalty < other.bound_penalty;
                if (self.cost != other.cost) return self.cost < other.cost;
                if (self.row_degree != other.row_degree) return self.row_degree < other.row_degree;
                if (self.pivot_magnitude != other.pivot_magnitude) return self.pivot_magnitude > other.pivot_magnitude;
            } else {
                // Bixby-style merit rewards a strong pivot in a sparse row and
                // column, with bound/cost compatibility acting as a penalty.
                const self_fill = @as(f64, @floatFromInt(self.degree)) * @as(f64, @floatFromInt(self.row_degree));
                const other_fill = @as(f64, @floatFromInt(other.degree)) * @as(f64, @floatFromInt(other.row_degree));
                const self_merit = self.pivot_magnitude /
                    (@sqrt(self_fill) * (1.0 + self.cost) * (1.0 + @as(f64, @floatFromInt(self.bound_penalty))));
                const other_merit = other.pivot_magnitude /
                    (@sqrt(other_fill) * (1.0 + other.cost) * (1.0 + @as(f64, @floatFromInt(other.bound_penalty))));
                if (self_merit != other_merit) return self_merit > other_merit;
                if (self.degree != other.degree) return self.degree < other.degree;
                if (self.row_degree != other.row_degree) return self.row_degree < other.row_degree;
            }
            // Final deterministic tiebreak by position.
            if (self.column != other.column) return self.column < other.column;
            return self.row < other.row;
        }
    };

    /// Classify a column's bound tightness for the ltssf penalty:
    /// 0 = free (one-sided or unbounded), 1 = boxed, 2 = fixed.
    fn boundPenalty(lower: f64, upper: f64) u8 {
        if (lower == upper) return 2;
        if (std.math.isFinite(lower) and std.math.isFinite(upper)) return 1;
        return 0;
    }
};

test "crash planner peels a triangular structural basis" {
    const foundation = @import("foundation");
    const view = matrix.CscView.initAssumeValid(
        3,
        3,
        &[_]usize{ 0, 1, 3, 5 },
        &[_]foundation.RowId{ foundation.RowId.fromUsize(0), foundation.RowId.fromUsize(0), foundation.RowId.fromUsize(1), foundation.RowId.fromUsize(1), foundation.RowId.fromUsize(2) },
        &[_]f64{ 2, 1, 3, 1, 4 },
    );
    var workspace = CrashWorkspace.init(std.testing.allocator);
    defer workspace.deinit();
    const plan_view = try workspace.plan(view, &[_]f64{ 0, 0, 0 }, &[_]f64{ 0, 0, 0 }, &[_]f64{ std.math.inf(f64), std.math.inf(f64), std.math.inf(f64) }, &[_]f64{ 1, 1, 1 }, &[_]f64{ 1, 1, 1 }, 2, 1e-12, .ltssf);
    try std.testing.expectEqual(@as(usize, 3), plan_view.rows.len);
    try std.testing.expectEqualSlices(u32, &.{ 0, 1, 2 }, plan_view.rows);
    try std.testing.expectEqualSlices(u32, &.{ 0, 1, 2 }, plan_view.columns);
}

test "bixby scoring can prefer pivot quality over a singleton" {
    const foundation = @import("foundation");
    const view = matrix.CscView.initAssumeValid(
        2,
        2,
        &[_]usize{ 0, 1, 3 },
        &[_]foundation.RowId{ foundation.RowId.fromUsize(0), foundation.RowId.fromUsize(0), foundation.RowId.fromUsize(1) },
        &[_]f64{ 0.1, 0.2, 10.0 },
    );
    var workspace = CrashWorkspace.init(std.testing.allocator);
    defer workspace.deinit();
    const costs = [_]f64{ 0, 0 };
    const lower = [_]f64{ 0, 0 };
    const upper = [_]f64{ std.math.inf(f64), std.math.inf(f64) };
    const scales = [_]f64{ 1, 1 };

    const ltssf = try workspace.plan(view, &costs, &lower, &upper, &scales, &scales, 2, 1e-12, .ltssf);
    try std.testing.expectEqual(@as(u32, 0), ltssf.columns[0]);
    const bixby = try workspace.plan(view, &costs, &lower, &upper, &scales, &scales, 2, 1e-12, .bixby);
    try std.testing.expectEqual(@as(u32, 1), bixby.columns[0]);
    try std.testing.expectEqual(@as(u32, 1), bixby.rows[0]);
}
