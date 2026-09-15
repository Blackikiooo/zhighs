//! Root separation algorithms over the immutable model rows.

const std = @import("std");
const solver = @import("solver");
const plugin = @import("plugin");
const matrix = @import("matrix");
const foundation = @import("foundation");
const CutPool = @import("cut_pool.zig").CutPool;

const Candidate = struct {
    column: foundation.ColId,
    coefficient: f64,
    primal: f64,
};

pub const Statistics = struct {
    rows_considered: usize = 0,
    cuts_added: usize = 0,
    duplicates: usize = 0,
};

/// Minimal-cover separator for nonnegative pure-binary knapsack rows.
pub const CoverSeparator = struct {
    allocator: std.mem.Allocator,
    candidates: []Candidate,
    cut_columns: []foundation.ColId,
    cut_values: []f64,

    pub fn init(allocator: std.mem.Allocator, num_columns: usize) std.mem.Allocator.Error!CoverSeparator {
        const candidates = try allocator.alloc(Candidate, num_columns);
        errdefer allocator.free(candidates);
        const cut_columns = try allocator.alloc(foundation.ColId, num_columns);
        errdefer allocator.free(cut_columns);
        const cut_values = try allocator.alloc(f64, num_columns);
        return .{ .allocator = allocator, .candidates = candidates, .cut_columns = cut_columns, .cut_values = cut_values };
    }

    pub fn deinit(self: *CoverSeparator) void {
        self.allocator.free(self.candidates);
        self.allocator.free(self.cut_columns);
        self.allocator.free(self.cut_values);
        self.* = undefined;
    }

    pub fn separate(
        self: *CoverSeparator,
        problem: solver.LpProblemView,
        rows: matrix.CsrView,
        integrality: []const plugin.Integrality,
        primal: []const f64,
        pool: *CutPool,
        tolerance: f64,
        max_cuts: usize,
    ) !Statistics {
        var statistics = Statistics{};
        if (max_cuts == 0) return statistics;
        for (0..problem.num_rows) |row_index| {
            if (statistics.cuts_added >= max_cuts) break;
            const rhs = problem.row_upper[row_index];
            if (!std.math.isFinite(rhs) or rhs < 0.0) continue;
            const row = rows.rowAssumeValid(row_index);
            if (row.nnz() == 0 or row.nnz() > self.candidates.len) continue;
            statistics.rows_considered += 1;

            var eligible = true;
            for (row.indices, row.values, 0..) |column, coefficient, position| {
                const index = column.toUsize();
                if (coefficient <= 0.0 or integrality[index] != .binary or
                    problem.col_lower[index] < -tolerance or problem.col_upper[index] > 1.0 + tolerance)
                {
                    eligible = false;
                    break;
                }
                const value = primal[index];
                self.candidates[position] = .{ .column = column, .coefficient = coefficient, .primal = value };
            }
            if (!eligible) continue;

            const candidates = self.candidates[0..row.nnz()];
            std.mem.sort(Candidate, candidates, {}, candidateBefore);
            var cover_size: usize = 0;
            var coefficient_sum: f64 = 0.0;
            while (cover_size < candidates.len and coefficient_sum <= rhs + tolerance) : (cover_size += 1)
                coefficient_sum += candidates[cover_size].coefficient;
            if (coefficient_sum <= rhs + tolerance or cover_size == 0) continue;

            var cover_activity: f64 = 0.0;
            for (candidates[0..cover_size], 0..) |candidate, position| {
                self.cut_columns[position] = candidate.column;
                self.cut_values[position] = 1.0;
                cover_activity += candidate.primal;
            }
            const cover_rhs = @as(f64, @floatFromInt(cover_size - 1));
            const violation = cover_activity - cover_rhs;
            if (violation <= tolerance) continue;
            std.mem.sort(foundation.ColId, self.cut_columns[0..cover_size], {}, columnBefore);
            const cut = matrix.SparseVectorView(foundation.ColId).initAssumeValid(
                problem.num_cols,
                self.cut_columns[0..cover_size],
                self.cut_values[0..cover_size],
            );
            switch (try pool.add(self.allocator, cut, cover_rhs, violation / @sqrt(@as(f64, @floatFromInt(cover_size))))) {
                .added => statistics.cuts_added += 1,
                .existing => statistics.duplicates += 1,
            }
        }
        return statistics;
    }
};

fn candidateBefore(_: void, lhs: Candidate, rhs: Candidate) bool {
    if (lhs.primal != rhs.primal) return lhs.primal > rhs.primal;
    if (lhs.coefficient != rhs.coefficient) return lhs.coefficient > rhs.coefficient;
    return lhs.column.toUsize() < rhs.column.toUsize();
}

fn columnBefore(_: void, lhs: foundation.ColId, rhs: foundation.ColId) bool {
    return lhs.toUsize() < rhs.toUsize();
}

test "cover separator adds and deduplicates a violated binary cover" {
    const base = solver.LpProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 3.0, 2.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ 1.0, 1.0 },
        .row_lower = &.{-std.math.inf(f64)},
        .row_upper = &.{4.0},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &.{ 0, 1, 2 }, &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) }, &.{ 2.0, 3.0 }),
        .objective_sense = .maximize,
        .objective_offset = 0.0,
    };
    const starts = try std.testing.allocator.dupe(usize, base.matrix.col_starts);
    const row_ids = try std.testing.allocator.dupe(foundation.RowId, base.matrix.row_indices);
    const values = try std.testing.allocator.dupe(f64, base.matrix.values);
    var owned = matrix.CscMatrix.initOwnedSlicesAssumeValid(1, 2, starts, row_ids, values);
    defer owned.deinit(std.testing.allocator);
    var csr = try matrix.CsrCache.buildAssumeValid(std.testing.allocator, owned, 0);
    defer csr.deinit(std.testing.allocator);
    var pool = try CutPool.init(std.testing.allocator, 2);
    defer pool.deinit(std.testing.allocator);
    var separator = try CoverSeparator.init(std.testing.allocator, 2);
    defer separator.deinit();

    const first = try separator.separate(base, csr.viewAssumeCurrent(), &.{ .binary, .binary }, &.{ 1.0, 2.0 / 3.0 }, &pool, 1e-6, 4);
    try std.testing.expectEqual(@as(usize, 1), first.cuts_added);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0 }, pool.row(@enumFromInt(0)).values);
    try std.testing.expectEqual(@as(f64, 1.0), pool.rhs(@enumFromInt(0)));
    const second = try separator.separate(base, csr.viewAssumeCurrent(), &.{ .binary, .binary }, &.{ 1.0, 2.0 / 3.0 }, &pool, 1e-6, 4);
    try std.testing.expectEqual(@as(usize, 0), second.cuts_added);
    try std.testing.expectEqual(@as(usize, 1), second.duplicates);
}
