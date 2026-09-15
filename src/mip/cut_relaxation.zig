//! Materialized LP view combining the immutable model with active pool cuts.

const std = @import("std");
const solver = @import("solver");
const matrix = @import("matrix");
const CutPool = @import("cut_pool.zig").CutPool;

pub const CutRelaxation = struct {
    allocator: std.mem.Allocator,
    base: solver.LpProblemView,
    base_matrix: matrix.CscMatrix,
    base_rows: matrix.CsrCache,
    augmented_matrix: ?matrix.CscMatrix = null,
    row_lower: []f64 = &.{},
    row_upper: []f64 = &.{},
    pool_revision: ?usize = null,
    structure_revision: usize = 0,

    pub fn init(allocator: std.mem.Allocator, base: solver.LpProblemView) !CutRelaxation {
        const starts = try allocator.dupe(usize, base.matrix.col_starts);
        errdefer allocator.free(starts);
        const rows = try allocator.dupe(@TypeOf(base.matrix.row_indices[0]), base.matrix.row_indices);
        errdefer allocator.free(rows);
        const values = try allocator.dupe(f64, base.matrix.values);
        var base_matrix = matrix.CscMatrix.initOwnedSlicesAssumeValid(base.num_rows, base.num_cols, starts, rows, values);
        errdefer base_matrix.deinit(allocator);
        var base_rows = try matrix.CsrCache.buildAssumeValid(allocator, base_matrix, 0);
        errdefer base_rows.deinit(allocator);
        return .{ .allocator = allocator, .base = base, .base_matrix = base_matrix, .base_rows = base_rows, .pool_revision = 0 };
    }

    pub fn deinit(self: *CutRelaxation) void {
        if (self.augmented_matrix) |*value| value.deinit(self.allocator);
        self.allocator.free(self.row_lower);
        self.allocator.free(self.row_upper);
        self.base_rows.deinit(self.allocator);
        self.base_matrix.deinit(self.allocator);
        self.* = undefined;
    }

    pub fn originalRows(self: *const CutRelaxation) matrix.CsrView {
        return self.base_rows.viewAssumeCurrent();
    }

    /// Rebuild only when active pool membership changed.
    pub fn synchronize(self: *CutRelaxation, pool: *const CutPool) !bool {
        if (self.pool_revision == pool.revision) return false;
        var next_matrix = try pool.appendActiveToCsc(self.allocator, self.base_matrix);
        errdefer next_matrix.deinit(self.allocator);
        const row_count = self.base.num_rows + pool.activeCount();
        const next_lower = try self.allocator.alloc(f64, row_count);
        errdefer self.allocator.free(next_lower);
        const next_upper = try self.allocator.alloc(f64, row_count);
        @memcpy(next_lower[0..self.base.num_rows], self.base.row_lower);
        @memcpy(next_upper[0..self.base.num_rows], self.base.row_upper);
        @memset(next_lower[self.base.num_rows..], -std.math.inf(f64));
        pool.fillActiveUpper(next_upper[self.base.num_rows..]);

        if (self.augmented_matrix) |*value| value.deinit(self.allocator);
        self.allocator.free(self.row_lower);
        self.allocator.free(self.row_upper);
        self.augmented_matrix = next_matrix;
        self.row_lower = next_lower;
        self.row_upper = next_upper;
        self.pool_revision = pool.revision;
        self.structure_revision +%= 1;
        return true;
    }

    pub fn view(self: *const CutRelaxation, column_lower: []const f64, column_upper: []const f64) solver.LpProblemView {
        var result = self.base;
        result.col_lower = column_lower;
        result.col_upper = column_upper;
        if (self.augmented_matrix) |value| {
            result.num_rows = value.num_rows;
            result.row_lower = self.row_lower;
            result.row_upper = self.row_upper;
            result.matrix = value.view();
        }
        return result;
    }
};

test "cut relaxation appends active pool rows with stable base CSR" {
    const foundation = @import("foundation");
    const base = solver.LpProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 1.0, 1.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ 1.0, 1.0 },
        .row_lower = &.{-std.math.inf(f64)},
        .row_upper = &.{1.5},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &.{ 0, 1, 2 }, &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) }, &.{ 1.0, 1.0 }),
        .objective_sense = .maximize,
        .objective_offset = 0.0,
    };
    var cuts = try CutPool.init(std.testing.allocator, 2);
    defer cuts.deinit(std.testing.allocator);
    var relaxation = try CutRelaxation.init(std.testing.allocator, base);
    defer relaxation.deinit();
    const columns = [_]foundation.ColId{ try foundation.ColId.init(0), try foundation.ColId.init(1) };
    _ = try cuts.add(std.testing.allocator, matrix.SparseVectorView(foundation.ColId).initAssumeValid(2, &columns, &.{ 1.0, 1.0 }), 1.0, 0.5);
    try std.testing.expect(try relaxation.synchronize(&cuts));
    const view = relaxation.view(base.col_lower, base.col_upper);
    try std.testing.expectEqual(@as(usize, 2), view.num_rows);
    try std.testing.expectEqualSlices(f64, &.{ 1.5, 1.0 }, view.row_upper);
    try std.testing.expectEqual(@as(usize, 1), relaxation.originalRows().num_rows);
}
