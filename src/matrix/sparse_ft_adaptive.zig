//! Dimension-adaptive Forrest--Tomlin storage.
//!
//! Small factors retain the compact intrusive representation; larger factors
//! use contiguous U columns with reciprocal row slabs. Exactly one backend is
//! live and updated, so this dispatch does not double-write numerical state.

const std = @import("std");
const linked = @import("sparse_ft.zig");
const contiguous = @import("sparse_ft_packed.zig");

pub const FtError = linked.FtError;
pub const MutableUpperView = linked.MutableUpperView;

const packed_dimension_minimum = 128;

pub const SparseForrestTomlin = struct {
    allocator: std.mem.Allocator,
    backend: Backend,

    const Backend = union(enum) {
        linked: linked.SparseForrestTomlin,
        contiguous: contiguous.SparseForrestTomlin,
    };

    pub fn init(allocator: std.mem.Allocator) SparseForrestTomlin {
        return .{
            .allocator = allocator,
            .backend = .{ .linked = linked.SparseForrestTomlin.init(allocator) },
        };
    }

    pub fn deinit(self: *SparseForrestTomlin) void {
        switch (self.backend) {
            inline else => |*ft| ft.deinit(),
        }
        self.* = init(self.allocator);
    }

    pub fn reset(
        self: *SparseForrestTomlin,
        pivot_ids: []const u32,
        pivot_values: []const f64,
        column_starts: []const usize,
        column_rows: []const u32,
        column_values: []const f64,
    ) FtError!void {
        const use_contiguous = pivot_ids.len >= packed_dimension_minimum;
        switch (self.backend) {
            .linked => |*ft| {
                if (!use_contiguous) return ft.reset(pivot_ids, pivot_values, column_starts, column_rows, column_values);
                ft.deinit();
                self.backend = .{ .contiguous = contiguous.SparseForrestTomlin.init(self.allocator) };
                return self.backend.contiguous.reset(pivot_ids, pivot_values, column_starts, column_rows, column_values);
            },
            .contiguous => |*ft| {
                if (use_contiguous) return ft.reset(pivot_ids, pivot_values, column_starts, column_rows, column_values);
                ft.deinit();
                self.backend = .{ .linked = linked.SparseForrestTomlin.init(self.allocator) };
                return self.backend.linked.reset(pivot_ids, pivot_values, column_starts, column_rows, column_values);
            },
        }
    }

    pub fn mutableUpperView(self: *SparseForrestTomlin) MutableUpperView {
        return switch (self.backend) {
            inline else => |*ft| ft.mutableUpperView(),
        };
    }

    pub fn hasUpdates(self: *const SparseForrestTomlin) bool {
        return switch (self.backend) {
            inline else => |*ft| ft.hasUpdates(),
        };
    }

    pub fn updateCount(self: *const SparseForrestTomlin) usize {
        return switch (self.backend) {
            inline else => |ft| ft.update_count,
        };
    }

    pub fn retainedBytes(self: *const SparseForrestTomlin) usize {
        return switch (self.backend) {
            inline else => |*ft| ft.retainedBytes(),
        };
    }

    pub fn prepareFtran(self: *SparseForrestTomlin, values: []f64, capture: bool) FtError!void {
        return switch (self.backend) {
            inline else => |*ft| ft.prepareFtran(values, capture),
        };
    }

    pub fn solveUpper(self: *SparseForrestTomlin, values: []f64) FtError!void {
        return switch (self.backend) {
            inline else => |*ft| ft.solveUpper(values),
        };
    }

    pub fn lastFtranOutputIndices(self: *const SparseForrestTomlin) []const u32 {
        return switch (self.backend) {
            inline else => |*ft| ft.lastFtranOutputIndices(),
        };
    }

    pub fn solveUpperTranspose(self: *SparseForrestTomlin, values: []f64) FtError!void {
        return switch (self.backend) {
            inline else => |*ft| ft.solveUpperTranspose(values),
        };
    }

    pub fn captureEp(self: *SparseForrestTomlin, leaving_id: u32) FtError!void {
        return switch (self.backend) {
            inline else => |*ft| ft.captureEp(leaving_id),
        };
    }

    pub fn update(self: *SparseForrestTomlin, leaving_id: u32, alpha: f64, zero_tolerance: f64) FtError!void {
        return switch (self.backend) {
            inline else => |*ft| ft.update(leaving_id, alpha, zero_tolerance),
        };
    }
};

test "adaptive FT selects exactly one backend at the dimension boundary" {
    var ft = SparseForrestTomlin.init(std.testing.allocator);
    defer ft.deinit();

    var small_pivots: [packed_dimension_minimum - 1]u32 = undefined;
    var small_diagonal: [small_pivots.len]f64 = undefined;
    var small_starts: [small_pivots.len + 1]usize = @splat(0);
    for (&small_pivots, &small_diagonal, 0..) |*pivot, *diagonal, index| {
        pivot.* = @intCast(index);
        diagonal.* = 1;
    }
    try ft.reset(&small_pivots, &small_diagonal, &small_starts, &.{}, &.{});
    try std.testing.expect(ft.backend == .linked);

    var large_pivots: [packed_dimension_minimum]u32 = undefined;
    var large_diagonal: [large_pivots.len]f64 = undefined;
    var large_starts: [large_pivots.len + 1]usize = @splat(0);
    for (&large_pivots, &large_diagonal, 0..) |*pivot, *diagonal, index| {
        pivot.* = @intCast(index);
        diagonal.* = 1;
    }
    try ft.reset(&large_pivots, &large_diagonal, &large_starts, &.{}, &.{});
    try std.testing.expect(ft.backend == .contiguous);
    try ft.reset(&small_pivots, &small_diagonal, &small_starts, &.{}, &.{});
    try std.testing.expect(ft.backend == .linked);
}
