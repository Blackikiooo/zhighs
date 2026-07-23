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

const Entry = struct {
    /// Structural column of an appended row entry.
    col: foundation.ColId,
    /// Finite nonzero coefficient parallel to `col`.
    value: f64,
};

const EntryList = std.MultiArrayList(Entry);

pub const DynamicRowMatrix = struct {
    /// Fixed structural-column dimension shared by every appended row.
    num_cols: usize,
    /// CSR-style offsets; length is the number of dynamic rows plus one.
    row_starts: std.ArrayList(usize) = .empty,
    // One allocation with independent contiguous field streams. Column IDs and
    // values can no longer grow to inconsistent lengths, while numerical code
    // still receives cache-friendly SoA slices.
    entries: EntryList = .empty,

    const Self = @This();

    pub const Checkpoint = struct {
        /// Row count retained by a rollback to this checkpoint.
        num_rows: usize,
        /// Entry count retained by a rollback to this checkpoint.
        nnz: usize,
    };

    /// Construct an empty row store with the initial zero row offset.
    pub fn init(allocator: std.mem.Allocator, num_cols: usize) (std.mem.Allocator.Error || csc.MatrixError)!Self {
        try csc.validateDimensions(0, num_cols);
        var self: Self = .{ .num_cols = num_cols };
        errdefer self.deinit(allocator);
        try self.row_starts.append(allocator, 0);
        return self;
    }

    /// Release row offsets and the SoA entry allocation.
    pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
        self.row_starts.deinit(allocator);
        self.entries.deinit(allocator);
        self.* = undefined;
    }

    /// Number of currently appended rows.
    pub inline fn numRows(self: Self) usize {
        return self.row_starts.items.len - 1;
    }

    /// Number of stored dynamic-row entries.
    pub inline fn nnz(self: Self) usize {
        return self.entries.len;
    }

    /// Capture an O(1) rollback marker for the current logical lengths.
    pub inline fn checkpoint(self: Self) Checkpoint {
        return .{ .num_rows = self.numRows(), .nnz = self.nnz() };
    }

    /// Reserves both backing allocations for a batch of row appends. Callers
    /// that know the batch shape can hoist all growth out of the append loop.
    pub fn reserve(self: *Self, allocator: std.mem.Allocator, additional_rows: usize, additional_nnz: usize) (std.mem.Allocator.Error || csc.MatrixError)!void {
        const new_num_rows = std.math.add(usize, self.numRows(), additional_rows) catch return error.DimensionTooLarge;
        try csc.validateDimensions(new_num_rows, self.num_cols);
        try self.entries.ensureUnusedCapacity(allocator, additional_nnz);
        try self.row_starts.ensureUnusedCapacity(allocator, additional_rows);
    }

    /// Appends one canonical row and returns its local RowId.
    pub fn appendRow(self: *Self, allocator: std.mem.Allocator, row_view: sparse_vector.SparseVectorView(foundation.ColId)) (std.mem.Allocator.Error || csc.MatrixError)!foundation.RowId {
        if (row_view.dimension != self.num_cols) return error.DimensionMismatch;
        try row_view.validate();

        try self.reserve(allocator, 1, row_view.nnz());
        return self.appendRowPreReserved(row_view);
    }

    /// Trusted append for validated rows after a sufficient `reserve` call.
    /// The row view must be canonical, dimension-compatible, and its storage
    /// must not alias this matrix's entry allocation.
    pub fn appendRowPreReserved(self: *Self, row_view: sparse_vector.SparseVectorView(foundation.ColId)) foundation.RowId {
        std.debug.assert(row_view.dimension == self.num_cols);
        std.debug.assert(self.entries.capacity - self.entries.len >= row_view.nnz());
        std.debug.assert(self.row_starts.capacity - self.row_starts.items.len >= 1);
        const row_id = foundation.RowId.fromUsize(self.numRows()) catch unreachable;
        for (row_view.indices, row_view.values) |col_id, value| {
            self.entries.appendAssumeCapacity(.{ .col = col_id, .value = value });
        }
        self.row_starts.appendAssumeCapacity(self.nnz());
        return row_id;
    }

    /// Return a checked borrowed canonical view of one dynamic row.
    pub fn row(self: Self, row_id: foundation.RowId) csc.MatrixError!sparse_vector.SparseVectorView(foundation.ColId) {
        const index = row_id.toUsize();
        if (index >= self.numRows()) return error.IndexOutOfBounds;
        return self.rowAssumeValid(index);
    }

    /// Trusted row access for an in-range local row index.
    pub inline fn rowAssumeValid(self: Self, index: usize) sparse_vector.SparseVectorView(foundation.ColId) {
        const begin = self.row_starts.items[index];
        const end = self.row_starts.items[index + 1];
        const fields = self.entries.slice();
        return sparse_vector.SparseVectorView(foundation.ColId).initAssumeValid(self.num_cols, fields.items(.col)[begin..end], fields.items(.value)[begin..end]);
    }

    /// Removes every row appended after checkpoint in O(1).
    pub fn rollback(self: *Self, target: Checkpoint) csc.MatrixError!void {
        if (target.num_rows > self.numRows() or target.nnz > self.nnz()) return error.InvalidCheckpoint;
        if (self.row_starts.items[target.num_rows] != target.nnz) return error.InvalidCheckpoint;
        self.row_starts.items.len = target.num_rows + 1;
        self.entries.shrinkRetainingCapacity(target.nnz);
    }

    /// Remove every dynamic row while retaining both backing allocations.
    pub fn clearRetainingCapacity(self: *Self) void {
        self.row_starts.items.len = 1;
        self.row_starts.items[0] = 0;
        self.entries.clearRetainingCapacity();
    }

    /// Batch-merges dynamic rows beneath a base CSC matrix.
    pub fn appendToCsc(self: Self, allocator: std.mem.Allocator, base: csc.CscMatrix) (std.mem.Allocator.Error || csc.MatrixError)!csc.CscMatrix {
        if (base.num_cols != self.num_cols) return error.DimensionMismatch;
        try base.validate();
        const fields = self.entries.slice();
        return edit.appendRowsFromCsrAssumeValid(
            allocator,
            base,
            self.row_starts.items,
            fields.items(.col),
            fields.items(.value),
        );
    }
};

test "dynamic rows append access checkpoint and rollback" {
    var rows = try DynamicRowMatrix.init(std.testing.allocator, 4);
    defer rows.deinit(std.testing.allocator);
    var first_cols = [_]foundation.ColId{ try foundation.ColId.init(0), try foundation.ColId.init(3) };
    var first_values = [_]f64{ 2.0, -1.0 };
    const first_id = try rows.appendRow(std.testing.allocator, sparse_vector.SparseVectorView(foundation.ColId).initAssumeValid(4, &first_cols, &first_values));
    try std.testing.expectEqual(@as(usize, 0), first_id.toUsize());
    const saved = rows.checkpoint();
    var second_cols = [_]foundation.ColId{try foundation.ColId.init(2)};
    var second_values = [_]f64{5.0};
    _ = try rows.appendRow(std.testing.allocator, sparse_vector.SparseVectorView(foundation.ColId).initAssumeValid(4, &second_cols, &second_values));
    try std.testing.expectEqual(@as(usize, 2), rows.numRows());
    try rows.rollback(saved);
    try std.testing.expectEqual(@as(usize, 1), rows.numRows());
    try std.testing.expectEqual(@as(usize, 2), rows.nnz());
    const first = try rows.row(first_id);
    try std.testing.expectEqualSlices(f64, &.{ 2.0, -1.0 }, first.values);
}

test "dynamic rows batch reserve keeps append pointers stable" {
    var rows = try DynamicRowMatrix.init(std.testing.allocator, 3);
    defer rows.deinit(std.testing.allocator);
    try rows.reserve(std.testing.allocator, 2, 3);

    const reserved_fields = rows.entries.slice();
    const cols_ptr = reserved_fields.items(.col).ptr;
    const values_ptr = reserved_fields.items(.value).ptr;
    const starts_ptr = rows.row_starts.items.ptr;

    const first_cols = [_]foundation.ColId{ try foundation.ColId.init(0), try foundation.ColId.init(2) };
    const first_values = [_]f64{ 1.0, 2.0 };
    const second_cols = [_]foundation.ColId{try foundation.ColId.init(1)};
    const second_values = [_]f64{3.0};
    _ = rows.appendRowPreReserved(sparse_vector.SparseVectorView(foundation.ColId).initAssumeValid(3, &first_cols, &first_values));
    _ = rows.appendRowPreReserved(sparse_vector.SparseVectorView(foundation.ColId).initAssumeValid(3, &second_cols, &second_values));

    const appended_fields = rows.entries.slice();
    try std.testing.expectEqual(cols_ptr, appended_fields.items(.col).ptr);
    try std.testing.expectEqual(values_ptr, appended_fields.items(.value).ptr);
    try std.testing.expectEqual(starts_ptr, rows.row_starts.items.ptr);
    try std.testing.expectEqual(@as(usize, 2), rows.numRows());
    try std.testing.expectEqual(@as(usize, 3), rows.nnz());
}

test "dynamic rows batch merge into base CSC" {
    var base = try csc.CscMatrix.initZero(std.testing.allocator, 1, 2);
    defer base.deinit(std.testing.allocator);
    var rows = try DynamicRowMatrix.init(std.testing.allocator, 2);
    defer rows.deinit(std.testing.allocator);
    var cols = [_]foundation.ColId{ try foundation.ColId.init(0), try foundation.ColId.init(1) };
    var values = [_]f64{ 3.0, 4.0 };
    _ = try rows.appendRow(std.testing.allocator, sparse_vector.SparseVectorView(foundation.ColId).initAssumeValid(2, &cols, &values));
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
    try std.testing.expectError(error.DimensionMismatch, rows.appendRow(std.testing.allocator, sparse_vector.SparseVectorView(foundation.ColId).initAssumeValid(2, &no_cols, &no_values)));
    try std.testing.expectError(error.InvalidCheckpoint, rows.rollback(.{ .num_rows = 1, .nnz = 0 }));
}
