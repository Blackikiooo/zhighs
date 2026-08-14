//! Serial dual-simplex vectors whose support lives for one major iteration.
//!
//! The ownership mirrors `HEkkDual`: BTRAN produces `row_ep`, PRICE produces
//! `row_ap`, the entering FTRAN produces `col_aq`, and BFRT produces
//! `col_bfrt`.  CHUZC consumes one combined pack in producer order.  Keeping
//! these vectors out of `BasisState` prevents basis lifetime and iteration
//! scratch lifetime from being conflated.

const std = @import("std");
const indexed_vector = @import("indexed_vector.zig");

pub const DualIterationError = error{ DenseIndexState, PackCapacityExceeded, Overflow };

/// Allocation-free dense-addressable view of the pivotal row. HiGHS keeps
/// structural `row_ap` and logical `row_ep` as separate HVector instances;
/// this view preserves that ownership while giving scalar consumers one
/// global column namespace.
pub fn PivotalRowView(comptime Real: type) type {
    return struct {
        structural: []const Real,
        logical: []const Real,

        const Self = @This();

        pub inline fn init(structural: []const Real, logical: []const Real) Self {
            return .{ .structural = structural, .logical = logical };
        }

        pub inline fn fromDense(values: []const Real) Self {
            return .{ .structural = values, .logical = &.{} };
        }

        pub inline fn len(self: Self) usize {
            return self.structural.len + self.logical.len;
        }

        pub inline fn at(self: Self, column: usize) Real {
            return if (column < self.structural.len)
                self.structural[column]
            else
                self.logical[column - self.structural.len];
        }
    };
}

pub fn DualIterationWorkspace(comptime Real: type) type {
    const Vector = indexed_vector.IndexedVector(Real);
    return struct {
        allocator: std.mem.Allocator,
        num_rows: usize,
        num_cols: usize,

        /// Pivotal BTRAN row `B^-T e_p`, indexed by basis row.
        row_ep: Vector,
        /// Structural pivotal tableau row produced by PRICE.
        row_ap: Vector,
        /// Entering FTRAN column `B^-1 a_q`.
        col_aq: Vector,
        /// FTRAN of the accumulated bound-flip RHS.
        col_bfrt: Vector,
        /// DSE tau vector `B^-1 B^-T e_p` before the basis update.
        col_dse: Vector,

        /// CHUZC row: structural row_ap entries followed by logical row_ep.
        pack_indices: []u32,
        pack_values: []Real,
        pack_count: usize = 0,
        /// Mutable CHUZC work order and published BFRT flip prefix. This must
        /// not alias the producer-owned pivotal-row pack scanned by CHUZC2.
        candidate_indices: []u32,
        /// CHUZC ratio or breakpoint group keyed by global column.
        work_ratio: []Real,
        /// Tableau alpha times the leaving-row move, keyed by global column.
        /// Its sign is the selected nonbasic direction and its magnitude is
        /// the positive effective alpha consumed by later CHUZC stages.
        work_alpha: []Real,

        const Self = @This();
        pub const Scalar = Real;

        /// Allocate the four vectors and the combined CHUZC pack once.
        pub fn init(allocator: std.mem.Allocator, num_rows: usize, num_cols: usize) !Self {
            const total_columns = std.math.add(usize, num_cols, num_rows) catch return error.Overflow;
            if (total_columns > std.math.maxInt(u32)) return error.Overflow;
            var row_ep = try Vector.init(allocator, num_rows);
            errdefer row_ep.deinit();
            var row_ap = try Vector.init(allocator, num_cols);
            errdefer row_ap.deinit();
            var col_aq = try Vector.init(allocator, num_rows);
            errdefer col_aq.deinit();
            var col_bfrt = try Vector.init(allocator, num_rows);
            errdefer col_bfrt.deinit();
            var col_dse = try Vector.init(allocator, num_rows);
            errdefer col_dse.deinit();
            const pack_indices = try allocator.alloc(u32, total_columns);
            errdefer allocator.free(pack_indices);
            const pack_values = try allocator.alloc(Real, total_columns);
            errdefer allocator.free(pack_values);
            const candidate_indices = try allocator.alloc(u32, total_columns);
            errdefer allocator.free(candidate_indices);
            const work_ratio = try allocator.alloc(Real, total_columns);
            errdefer allocator.free(work_ratio);
            const work_alpha = try allocator.alloc(Real, total_columns);
            errdefer allocator.free(work_alpha);
            @memset(pack_indices, 0);
            @memset(pack_values, 0);
            @memset(candidate_indices, 0);
            @memset(work_ratio, std.math.inf(Real));
            @memset(work_alpha, 0);
            return .{
                .allocator = allocator,
                .num_rows = num_rows,
                .num_cols = num_cols,
                .row_ep = row_ep,
                .row_ap = row_ap,
                .col_aq = col_aq,
                .col_bfrt = col_bfrt,
                .col_dse = col_dse,
                .pack_indices = pack_indices,
                .pack_values = pack_values,
                .candidate_indices = candidate_indices,
                .work_ratio = work_ratio,
                .work_alpha = work_alpha,
            };
        }

        /// Release all iteration-owned storage.
        pub fn deinit(self: *Self) void {
            self.row_ep.deinit();
            self.row_ap.deinit();
            self.col_aq.deinit();
            self.col_bfrt.deinit();
            self.col_dse.deinit();
            self.allocator.free(self.pack_indices);
            self.allocator.free(self.pack_values);
            self.allocator.free(self.candidate_indices);
            self.allocator.free(self.work_ratio);
            self.allocator.free(self.work_alpha);
            self.* = undefined;
        }

        /// Begin a new iteration using each producer's retained support for
        /// sparse clear. Factorization and basis state are deliberately absent.
        pub fn beginIteration(self: *Self) void {
            self.row_ep.clear();
            self.row_ap.clear();
            self.col_aq.clear();
            self.col_bfrt.clear();
            self.col_dse.clear();
            self.pack_count = 0;
        }

        /// Pack PRICE structural entries first, then BTRAN logical entries,
        /// preserving the exact order published by both producers.
        pub fn buildPivotalRowPack(self: *Self) DualIterationError!void {
            const structural = self.row_ap.activeIndices() orelse return error.DenseIndexState;
            const logical = self.row_ep.activeIndices() orelse return error.DenseIndexState;
            if (structural.len + logical.len > self.pack_indices.len) return error.PackCapacityExceeded;
            var output: usize = 0;
            for (structural) |column| {
                self.pack_indices[output] = column;
                self.pack_values[output] = self.row_ap.values[column];
                output += 1;
            }
            for (logical) |row| {
                self.pack_indices[output] = @intCast(self.num_cols + row);
                self.pack_values[output] = self.row_ep.values[row];
                output += 1;
            }
            self.pack_count = output;
        }

        pub inline fn packedIndices(self: *const Self) []const u32 {
            return self.pack_indices[0..self.pack_count];
        }

        pub inline fn packedValues(self: *const Self) []const Real {
            return self.pack_values[0..self.pack_count];
        }

        pub inline fn pivotalRow(self: *const Self) PivotalRowView(Real) {
            return .init(self.row_ap.values, self.row_ep.values);
        }

        /// Bytes represented by the complete retained iteration workspace.
        pub fn retainedBytes(self: *const Self) usize {
            return self.row_ep.retainedBytes() + self.row_ap.retainedBytes() +
                self.col_aq.retainedBytes() + self.col_bfrt.retainedBytes() + self.col_dse.retainedBytes() +
                std.mem.sliceAsBytes(self.pack_indices).len +
                std.mem.sliceAsBytes(self.pack_values).len +
                std.mem.sliceAsBytes(self.candidate_indices).len +
                std.mem.sliceAsBytes(self.work_ratio).len +
                std.mem.sliceAsBytes(self.work_alpha).len;
        }
    };
}

test "dual iteration workspace preserves HiGHS CHUZC producer order" {
    const Workspace = DualIterationWorkspace(f64);
    var workspace = try Workspace.init(std.testing.allocator, 3, 4);
    defer workspace.deinit();

    workspace.row_ap.values[3] = 30;
    workspace.row_ap.values[1] = 10;
    try workspace.row_ap.appendIndex(3);
    try workspace.row_ap.appendIndex(1);
    workspace.row_ep.values[2] = 2;
    workspace.row_ep.values[0] = -1;
    try workspace.row_ep.appendIndex(2);
    try workspace.row_ep.appendIndex(0);

    try workspace.buildPivotalRowPack();
    try std.testing.expectEqualSlices(u32, &.{ 3, 1, 6, 4 }, workspace.packedIndices());
    try std.testing.expectEqualSlices(f64, &.{ 30, 10, 2, -1 }, workspace.packedValues());
    workspace.candidate_indices[0] = 6;
    workspace.candidate_indices[1] = 1;
    try std.testing.expectEqualSlices(u32, &.{ 3, 1, 6, 4 }, workspace.packedIndices());
}

test "dual iteration workspace clears only iteration-owned vector state" {
    const Workspace = DualIterationWorkspace(f128);
    var workspace = try Workspace.init(std.testing.allocator, 5, 7);
    defer workspace.deinit();
    workspace.col_aq.values[4] = 8;
    try workspace.col_aq.appendIndex(4);
    workspace.col_bfrt.values[1] = -3;
    try workspace.col_bfrt.appendIndex(1);
    workspace.col_dse.values[3] = 6;
    try workspace.col_dse.appendIndex(3);
    workspace.pack_count = 2;
    workspace.beginIteration();
    try std.testing.expectEqual(@as(f128, 0), workspace.col_aq.values[4]);
    try std.testing.expectEqual(@as(f128, 0), workspace.col_bfrt.values[1]);
    try std.testing.expectEqual(@as(f128, 0), workspace.col_dse.values[3]);
    try std.testing.expectEqual(@as(usize, 0), workspace.pack_count);
}

test "pivotal row view joins structural and logical values without copying" {
    const View = PivotalRowView(f64);
    const view = View.init(&.{ 2, -3 }, &.{ 5, 7, 11 });
    try std.testing.expectEqual(@as(usize, 5), view.len());
    try std.testing.expectEqual(@as(f64, 2), view.at(0));
    try std.testing.expectEqual(@as(f64, -3), view.at(1));
    try std.testing.expectEqual(@as(f64, 5), view.at(2));
    try std.testing.expectEqual(@as(f64, 11), view.at(4));
}
