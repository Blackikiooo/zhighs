//! Reusable caller-owned output and scratch storage for CSC transformations.

const std = @import("std");
const foundation = @import("foundation");
const csc = @import("csc.zig");

/// Stable API: fixed-capacity buffers shared by slice, permutation, and edit.
///
/// Every `...Into` result is a borrowed `CscView`. It remains valid until the
/// next operation writing these buffers or `deinit`. The buffers retain their
/// high-water allocation across operations; callers choose when to release it.
pub const CscTransformBuffers = struct {
    col_starts: []align(64) usize,
    row_indices: []align(64) foundation.RowId,
    values: []align(64) f64,
    index_scratch: []usize,

    pub fn initCapacity(
        allocator: std.mem.Allocator,
        num_cols_capacity: usize,
        nnz_capacity: usize,
        index_scratch_capacity: usize,
    ) (std.mem.Allocator.Error || csc.MatrixError)!CscTransformBuffers {
        if (num_cols_capacity == std.math.maxInt(usize)) return error.DimensionTooLarge;
        const starts = try allocator.alignedAlloc(usize, .@"64", num_cols_capacity + 1);
        errdefer allocator.free(starts);
        const rows = try allocator.alignedAlloc(foundation.RowId, .@"64", nnz_capacity);
        errdefer allocator.free(rows);
        const vals = try allocator.alignedAlloc(f64, .@"64", nnz_capacity);
        errdefer allocator.free(vals);
        const scratch = try allocator.alloc(usize, index_scratch_capacity);
        return .{
            .col_starts = starts,
            .row_indices = rows,
            .values = vals,
            .index_scratch = scratch,
        };
    }

    pub fn deinit(self: *CscTransformBuffers, allocator: std.mem.Allocator) void {
        allocator.free(self.col_starts);
        allocator.free(self.row_indices);
        allocator.free(self.values);
        allocator.free(self.index_scratch);
        self.* = undefined;
    }

    pub fn requireCapacity(self: CscTransformBuffers, num_cols: usize, nnz: usize, scratch: usize) csc.MatrixError!void {
        if (num_cols == std.math.maxInt(usize)) return error.DimensionTooLarge;
        if (self.col_starts.len < num_cols + 1 or
            self.row_indices.len < nnz or
            self.values.len < nnz or
            self.index_scratch.len < scratch)
            return error.BufferTooSmall;
    }

    pub fn viewAssumeValid(self: CscTransformBuffers, num_rows: usize, num_cols: usize, nnz: usize) csc.CscView {
        std.debug.assert(self.col_starts.len >= num_cols + 1);
        std.debug.assert(self.row_indices.len >= nnz);
        std.debug.assert(self.values.len >= nnz);
        return csc.CscView.initAssumeValid(
            num_rows,
            num_cols,
            self.col_starts[0 .. num_cols + 1],
            self.row_indices[0..nnz],
            self.values[0..nnz],
        );
    }
};

test "transform buffers expose aligned reusable capacity" {
    var buffers = try CscTransformBuffers.initCapacity(std.testing.allocator, 4, 8, 6);
    defer buffers.deinit(std.testing.allocator);
    try buffers.requireCapacity(4, 8, 6);
    try std.testing.expectEqual(@as(usize, 0), @intFromPtr(buffers.col_starts.ptr) % 64);
    try std.testing.expectEqual(@as(usize, 0), @intFromPtr(buffers.row_indices.ptr) % 64);
    try std.testing.expectEqual(@as(usize, 0), @intFromPtr(buffers.values.ptr) % 64);
    try std.testing.expectError(error.BufferTooSmall, buffers.requireCapacity(5, 8, 6));
    try std.testing.expectError(error.BufferTooSmall, buffers.requireCapacity(4, 9, 6));
    try std.testing.expectError(error.BufferTooSmall, buffers.requireCapacity(4, 8, 7));
}
