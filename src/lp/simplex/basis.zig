//! Mutable simplex basis state.
//!
//! Hot numerical vectors are stored as separate contiguous arrays.  This
//! module owns basis membership and solution vectors, but not LU factors.

const std = @import("std");

/// Public-compatible values keep zero-copy status views possible at the model
/// adapter boundary. `free` and `fixed` remain solver-internal extensions.
pub const BasisStatus = enum(i8) { basic = 0, at_lower = -1, at_upper = 1, superbasic = 2, free = 3, fixed = 4 };
pub const BasisError = error{ RowOutOfRange, ColumnOutOfRange, ColumnAlreadyBasic };

pub const BasisState = struct {
    allocator: std.mem.Allocator,
    num_structural_cols: usize = 0,
    num_rows: usize = 0,
    row_status: []BasisStatus = &.{},
    col_status: []BasisStatus = &.{},
    basic_index: []u32 = &.{},
    basic_pos: []u32 = &.{},
    primal: []f64 = &.{},
    dual: []f64 = &.{},
    reduced_cost: []f64 = &.{},
    pivot_direction: []f64 = &.{},
    basic_value: []f64 = &.{},
    basic_lower: []f64 = &.{},
    basic_upper: []f64 = &.{},
    basic_margin: []f64 = &.{},
    ratio_direction: []f64 = &.{},
    row_scale: []f64 = &.{},
    column_scale: []f64 = &.{},
    row_rhs: []f64 = &.{},
    col_lower: []f64 = &.{},
    col_upper: []f64 = &.{},
    artificial_sign: []f64 = &.{},
    rhs_work: []f64 = &.{},
    /// Residual/correction workspace used by iterative refinement. Keeping it
    /// beside the other row vectors makes every FTRAN refinement allocation-free.
    residual_work: []f64 = &.{},
    /// Devex/steepest-edge reference weights. These are solver-owned because
    /// pricing mutates them, while pricing scans remain contiguous SoA loops.
    col_edge_weight: []f64 = &.{},
    row_edge_weight: []f64 = &.{},
    /// Frozen nonbasic reference set for a primal Devex framework. Bytes keep
    /// the hot weighted-norm loop branch-light and avoid packed-bit updates.
    devex_reference: []u8 = &.{},
    dual_row: []f64 = &.{},
    tableau: []f64 = &.{},
    dual_ratio: []f64 = &.{},
    dual_direction: []f64 = &.{},
    flip_columns: []u32 = &.{},
    dual_candidate_rows: []u32 = &.{},
    dual_candidate_score: []f64 = &.{},
    published_primal: []f64 = &.{},
    published_dual: []f64 = &.{},
    published_reduced_cost: []f64 = &.{},
    unbounded_ray: []f64 = &.{},

    pub fn init(allocator: std.mem.Allocator, rows: usize, cols: usize) !BasisState {
        const total_cols = cols + 2 * rows;
        var self = BasisState{ .allocator = allocator, .num_structural_cols = cols, .num_rows = rows };
        errdefer self.deinit();
        self.row_status = try allocator.alloc(BasisStatus, rows);
        self.col_status = try allocator.alloc(BasisStatus, total_cols);
        self.basic_index = try allocator.alloc(u32, rows);
        self.basic_pos = try allocator.alloc(u32, total_cols);
        self.primal = try allocator.alloc(f64, total_cols);
        self.dual = try allocator.alloc(f64, rows);
        self.reduced_cost = try allocator.alloc(f64, total_cols);
        self.pivot_direction = try allocator.alloc(f64, rows);
        self.basic_value = try allocator.alloc(f64, rows);
        self.basic_lower = try allocator.alloc(f64, rows);
        self.basic_upper = try allocator.alloc(f64, rows);
        self.basic_margin = try allocator.alloc(f64, rows);
        self.ratio_direction = try allocator.alloc(f64, rows);
        self.row_scale = try allocator.alloc(f64, rows);
        self.column_scale = try allocator.alloc(f64, total_cols);
        self.row_rhs = try allocator.alloc(f64, rows);
        self.col_lower = try allocator.alloc(f64, total_cols);
        self.col_upper = try allocator.alloc(f64, total_cols);
        self.artificial_sign = try allocator.alloc(f64, rows);
        self.rhs_work = try allocator.alloc(f64, rows);
        self.residual_work = try allocator.alloc(f64, rows);
        self.col_edge_weight = try allocator.alloc(f64, total_cols);
        self.row_edge_weight = try allocator.alloc(f64, rows);
        self.devex_reference = try allocator.alloc(u8, total_cols);
        self.dual_row = try allocator.alloc(f64, rows);
        self.tableau = try allocator.alloc(f64, total_cols);
        self.dual_ratio = try allocator.alloc(f64, total_cols);
        self.dual_direction = try allocator.alloc(f64, total_cols);
        self.flip_columns = try allocator.alloc(u32, total_cols);
        self.dual_candidate_rows = try allocator.alloc(u32, rows);
        self.dual_candidate_score = try allocator.alloc(f64, rows);
        self.published_primal = try allocator.alloc(f64, cols);
        self.published_dual = try allocator.alloc(f64, rows);
        self.published_reduced_cost = try allocator.alloc(f64, cols);
        self.unbounded_ray = try allocator.alloc(f64, cols);
        @memset(self.row_status, .basic);
        @memset(self.col_status, .at_lower);
        @memset(self.basic_index, 0);
        @memset(self.basic_pos, 0);
        @memset(self.primal, 0.0);
        @memset(self.dual, 0.0);
        @memset(self.reduced_cost, 0.0);
        @memset(self.pivot_direction, 0.0);
        @memset(self.basic_value, 0.0);
        @memset(self.basic_lower, 0.0);
        @memset(self.basic_upper, std.math.inf(f64));
        @memset(self.basic_margin, 0.0);
        @memset(self.ratio_direction, 0.0);
        @memset(self.row_scale, 1.0);
        @memset(self.column_scale, 1.0);
        @memset(self.row_rhs, 0.0);
        @memset(self.col_lower, 0.0);
        @memset(self.col_upper, std.math.inf(f64));
        @memset(self.artificial_sign, 0.0);
        @memset(self.rhs_work, 0.0);
        @memset(self.residual_work, 0.0);
        @memset(self.col_edge_weight, 1.0);
        @memset(self.row_edge_weight, 1.0);
        @memset(self.devex_reference, 0);
        @memset(self.dual_row, 0.0);
        @memset(self.tableau, 0.0);
        @memset(self.dual_ratio, std.math.inf(f64));
        @memset(self.dual_direction, 0.0);
        @memset(self.flip_columns, 0);
        @memset(self.dual_candidate_rows, 0);
        @memset(self.dual_candidate_score, 0.0);
        @memset(self.published_primal, 0.0);
        @memset(self.published_dual, 0.0);
        @memset(self.published_reduced_cost, 0.0);
        @memset(self.unbounded_ray, 0.0);
        return self;
    }

    pub fn initializeSlackBasis(self: *BasisState) void {
        @memset(self.basic_pos, std.math.maxInt(u32));
        for (self.basic_index, 0..) |*col, row| {
            const logical_col = self.num_structural_cols + row;
            col.* = @intCast(logical_col);
            self.basic_pos[logical_col] = @intCast(row);
        }
        @memset(self.row_status, .basic);
        @memset(self.col_status, .at_lower);
        for (self.basic_index) |col| self.col_status[col] = .basic;
        const artificial_begin = self.num_structural_cols + self.num_rows;
        @memset(self.col_status[artificial_begin..], .fixed);
        @memset(self.col_upper[artificial_begin..], 0.0);
    }

    /// Apply the combinatorial part of a simplex pivot. Numerical factor
    /// updates are deliberately handled by `Factorization`, keeping this hot
    /// state transition allocation-free.
    pub fn applyPivot(self: *BasisState, leaving_row: usize, entering_col: usize, leaving_status: BasisStatus) BasisError!void {
        if (leaving_row >= self.basic_index.len) return error.RowOutOfRange;
        if (entering_col >= self.col_status.len) return error.ColumnOutOfRange;
        if (self.col_status[entering_col] == .basic) return error.ColumnAlreadyBasic;
        const leaving_col = self.basic_index[leaving_row];
        if (leaving_col >= self.col_status.len) return error.ColumnOutOfRange;
        self.basic_index[leaving_row] = @intCast(entering_col);
        self.basic_pos[entering_col] = @intCast(leaving_row);
        self.col_status[entering_col] = .basic;
        self.col_status[leaving_col] = leaving_status;
        self.basic_pos[leaving_col] = std.math.maxInt(u32);
    }

    pub fn deinit(self: *BasisState) void {
        self.allocator.free(self.row_status);
        self.allocator.free(self.col_status);
        self.allocator.free(self.basic_index);
        self.allocator.free(self.basic_pos);
        self.allocator.free(self.primal);
        self.allocator.free(self.dual);
        self.allocator.free(self.reduced_cost);
        self.allocator.free(self.pivot_direction);
        self.allocator.free(self.basic_value);
        self.allocator.free(self.basic_lower);
        self.allocator.free(self.basic_upper);
        self.allocator.free(self.basic_margin);
        self.allocator.free(self.ratio_direction);
        self.allocator.free(self.row_scale);
        self.allocator.free(self.column_scale);
        self.allocator.free(self.row_rhs);
        self.allocator.free(self.col_lower);
        self.allocator.free(self.col_upper);
        self.allocator.free(self.artificial_sign);
        self.allocator.free(self.rhs_work);
        self.allocator.free(self.residual_work);
        self.allocator.free(self.col_edge_weight);
        self.allocator.free(self.row_edge_weight);
        self.allocator.free(self.devex_reference);
        self.allocator.free(self.dual_row);
        self.allocator.free(self.tableau);
        self.allocator.free(self.dual_ratio);
        self.allocator.free(self.dual_direction);
        self.allocator.free(self.flip_columns);
        self.allocator.free(self.dual_candidate_rows);
        self.allocator.free(self.dual_candidate_score);
        self.allocator.free(self.published_primal);
        self.allocator.free(self.published_dual);
        self.allocator.free(self.published_reduced_cost);
        self.allocator.free(self.unbounded_ray);
    }

    /// Bytes explicitly requested for retained basis SoA storage. Allocator
    /// metadata and the borrowed model matrix are intentionally excluded.
    pub fn requestedBytes(self: *const BasisState) usize {
        var total: usize = 0;
        total += std.mem.sliceAsBytes(self.row_status).len;
        total += std.mem.sliceAsBytes(self.col_status).len;
        total += std.mem.sliceAsBytes(self.basic_index).len;
        total += std.mem.sliceAsBytes(self.basic_pos).len;
        total += std.mem.sliceAsBytes(self.primal).len;
        total += std.mem.sliceAsBytes(self.dual).len;
        total += std.mem.sliceAsBytes(self.reduced_cost).len;
        total += std.mem.sliceAsBytes(self.pivot_direction).len;
        total += std.mem.sliceAsBytes(self.basic_value).len;
        total += std.mem.sliceAsBytes(self.basic_lower).len;
        total += std.mem.sliceAsBytes(self.basic_upper).len;
        total += std.mem.sliceAsBytes(self.basic_margin).len;
        total += std.mem.sliceAsBytes(self.ratio_direction).len;
        total += std.mem.sliceAsBytes(self.row_scale).len;
        total += std.mem.sliceAsBytes(self.column_scale).len;
        total += std.mem.sliceAsBytes(self.row_rhs).len;
        total += std.mem.sliceAsBytes(self.col_lower).len;
        total += std.mem.sliceAsBytes(self.col_upper).len;
        total += std.mem.sliceAsBytes(self.artificial_sign).len;
        total += std.mem.sliceAsBytes(self.rhs_work).len;
        total += std.mem.sliceAsBytes(self.residual_work).len;
        total += std.mem.sliceAsBytes(self.col_edge_weight).len;
        total += std.mem.sliceAsBytes(self.row_edge_weight).len;
        total += std.mem.sliceAsBytes(self.devex_reference).len;
        total += std.mem.sliceAsBytes(self.dual_row).len;
        total += std.mem.sliceAsBytes(self.tableau).len;
        total += std.mem.sliceAsBytes(self.dual_ratio).len;
        total += std.mem.sliceAsBytes(self.dual_direction).len;
        total += std.mem.sliceAsBytes(self.flip_columns).len;
        total += std.mem.sliceAsBytes(self.dual_candidate_rows).len;
        total += std.mem.sliceAsBytes(self.dual_candidate_score).len;
        total += std.mem.sliceAsBytes(self.published_primal).len;
        total += std.mem.sliceAsBytes(self.published_dual).len;
        total += std.mem.sliceAsBytes(self.published_reduced_cost).len;
        total += std.mem.sliceAsBytes(self.unbounded_ray).len;
        return total;
    }
};

test {
    std.testing.refAllDecls(@This());
}

test "BasisState pivot updates membership maps" {
    var basis = try BasisState.init(std.testing.allocator, 1, 2);
    defer basis.deinit();
    basis.initializeSlackBasis();
    try std.testing.expectEqual(@as(u32, 2), basis.basic_index[0]);
    try basis.applyPivot(0, 1, .at_lower);
    try std.testing.expectEqual(@as(u32, 1), basis.basic_index[0]);
    try std.testing.expectEqual(BasisStatus.basic, basis.col_status[1]);
    try std.testing.expectEqual(BasisStatus.at_lower, basis.col_status[2]);
    try std.testing.expectEqual(std.math.maxInt(u32), basis.basic_pos[2]);
}
