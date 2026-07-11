//! Append-only dynamic sparse rows with cheap checkpoint and rollback.
//!
//! MIP cuts and presolve-generated rows should not repeatedly rebuild the base
//! CSC matrix. This store keeps canonical rows in CSR-like contiguous streams;
//! a batch merge creates a new CSC only at an explicit synchronization point.

const std = @import("std");
const foundation = @import("foundation");
const sparse_vector = @import("sparse_vector.zig");
const csc = @import("csc.zig");
const edit = @import("edit.zig");

pub const DynamicRowMatrix = struct {
    num_cols: usize,
    row_starts: std.ArrayList(usize) = .empty,
    col_indices: std.ArrayList(foundation.ColId) = .empty,
    values: std.ArrayList(f64) = .empty,

    const Self = @This();

    pub const Checkpoint = struct {
        num_rows: usize,
        nnz: usize,
    };

    pub fn init(allocator: std.mem.Allocator, num_cols: usize) (std.mem.Allocator.Error || csc.MatrixError)!Self {
        try csc.validateDimensions(0, num_cols);
        var self: Self = .{ .num_cols = num_cols };
        errdefer self.deinit(allocator);
        try self.row_starts.append(allocator, 0);
        return self;
    }

    pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
        self.row_starts.deinit(allocator);
        self.col_indices.deinit(allocator);
        self.values.deinit(allocator);
        self.* = undefined;
    }

    pub inline fn numRows(self: Self) usize {
        return self.row_starts.items.len - 1;
    }

    pub inline fn nnz(self: Self) usize {
        return self.values.items.len;
    }

    pub inline fn checkpoint(self: Self) Checkpoint {
        return .{ .num_rows = self.numRows(), .nnz = self.nnz() };
    }

    /// Appends one canonical row and returns its local RowId.
    pub fn appendRow(self: *Self, allocator: std.mem.Allocator, row_view: sparse_vector.SparseVectorView(foundation.ColId)) (std.mem.Allocator.Error || csc.MatrixError)!foundation.RowId {
        if (row_view.dimension != self.num_cols) return error.DimensionMismatch;
        try row_view.validate();
        const row_id = foundation.RowId.fromUsize(self.numRows()) catch return error.DimensionTooLarge;

        // Reserve all streams before mutation. Once reservation succeeds, the
        // following appends cannot fail and parallel lengths stay coherent.
        try self.col_indices.ensureUnusedCapacity(allocator, row_view.nnz());
        try self.values.ensureUnusedCapacity(allocator, row_view.nnz());
        try self.row_starts.ensureUnusedCapacity(allocator, 1);
        for (row_view.indices, row_view.values) |col_id, value| {
            self.col_indices.appendAssumeCapacity(col_id);
            self.values.appendAssumeCapacity(value);
        }
        self.row_starts.appendAssumeCapacity(self.nnz());
        return row_id;
    }

    pub fn row(self: Self, row_id: foundation.RowId) csc.MatrixError!sparse_vector.SparseVectorView(foundation.ColId) {
        const index = row_id.toUsize();
        if (index >= self.numRows()) return error.IndexOutOfBounds;
        return self.rowAssumeValid(index);
    }

    pub inline fn rowAssumeValid(self: Self, index: usize) sparse_vector.SparseVectorView(foundation.ColId) {
        const begin = self.row_starts.items[index];
        const end = self.row_starts.items[index + 1];
        return .{ .dimension = self.num_cols, .indices = self.col_indices.items[begin..end], .values = self.values.items[begin..end] };
    }

    /// Removes every row appended after checkpoint in O(1).
    pub fn rollback(self: *Self, target: Checkpoint) csc.MatrixError!void {
        if (target.num_rows > self.numRows() or target.nnz > self.nnz()) return error.InvalidCheckpoint;
        if (self.row_starts.items[target.num_rows] != target.nnz) return error.InvalidCheckpoint;
        self.row_starts.items.len = target.num_rows + 1;
        self.col_indices.items.len = target.nnz;
        self.values.items.len = target.nnz;
    }

    pub fn clearRetainingCapacity(self: *Self) void {
        self.row_starts.items.len = 1;
        self.row_starts.items[0] = 0;
        self.col_indices.clearRetainingCapacity();
        self.values.clearRetainingCapacity();
    }

    /// Batch-merges dynamic rows beneath a base CSC matrix.
    pub fn appendToCsc(self: Self, allocator: std.mem.Allocator, base: csc.CscMatrix) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
        if (base.num_cols != self.num_cols) return error.DimensionMismatch;
        const views = try allocator.alloc(sparse_vector.SparseVectorView(foundation.ColId), self.numRows());
        defer allocator.free(views);
        for (views, 0..) |*view, row_index| view.* = self.rowAssumeValid(row_index);
        return edit.appendRows(allocator, base, views);
    }
};

test "dynamic rows append access checkpoint and rollback" {
    var rows = try DynamicRowMatrix.init(std.testing.allocator, 4);
    defer rows.deinit(std.testing.allocator);
    var first_cols = [_]foundation.ColId{ try foundation.ColId.init(0), try foundation.ColId.init(3) };
    var first_values = [_]f64{ 2.0, -1.0 };
    const first_id = try rows.appendRow(std.testing.allocator, .{ .dimension = 4, .indices = &first_cols, .values = &first_values });
    try std.testing.expectEqual(@as(usize, 0), first_id.toUsize());
    const saved = rows.checkpoint();
    var second_cols = [_]foundation.ColId{try foundation.ColId.init(2)};
    var second_values = [_]f64{5.0};
    _ = try rows.appendRow(std.testing.allocator, .{ .dimension = 4, .indices = &second_cols, .values = &second_values });
    try std.testing.expectEqual(@as(usize, 2), rows.numRows());
    try rows.rollback(saved);
    try std.testing.expectEqual(@as(usize, 1), rows.numRows());
    try std.testing.expectEqual(@as(usize, 2), rows.nnz());
    const first = try rows.row(first_id);
    try std.testing.expectEqualSlices(f64, &.{ 2.0, -1.0 }, first.values);
}

test "dynamic rows batch merge into base CSC" {
    var base = try csc.CscMatrix.initZero(std.testing.allocator, 1, 2);
    defer base.deinit(std.testing.allocator);
    var rows = try DynamicRowMatrix.init(std.testing.allocator, 2);
    defer rows.deinit(std.testing.allocator);
    var cols = [_]foundation.ColId{ try foundation.ColId.init(0), try foundation.ColId.init(1) };
    var values = [_]f64{ 3.0, 4.0 };
    _ = try rows.appendRow(std.testing.allocator, .{ .dimension = 2, .indices = &cols, .values = &values });
    var merged = try rows.appendToCsc(std.testing.allocator, base);
    defer merged.deinit(std.testing.allocator);
    try merged.validate();
    try std.testing.expectEqual(@as(usize, 2), merged.num_rows);
    try std.testing.expectEqualSlices(usize, &.{ 0, 1, 2 }, merged.col_starts);
    try std.testing.expectEqual(@as(usize, 1), merged.row_indices[0].toUsize());
}

test "dynamic row validation and invalid checkpoints" {
    var rows = try DynamicRowMatrix.init(std.testing.allocator, 1);
    defer rows.deinit(std.testing.allocator);
    var no_cols = [_]foundation.ColId{};
    var no_values = [_]f64{};
    try std.testing.expectError(error.DimensionMismatch, rows.appendRow(std.testing.allocator, .{ .dimension = 2, .indices = &no_cols, .values = &no_values }));
    try std.testing.expectError(error.InvalidCheckpoint, rows.rollback(.{ .num_rows = 1, .nnz = 0 }));
}
