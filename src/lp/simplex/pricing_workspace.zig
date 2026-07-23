//! Reusable CSR companion view for row-oriented simplex pricing.
//!
//! The model remains borrowed in CSC form. This workspace is rebuilt once at
//! the solve boundary and then reused without allocation by every row-pricing
//! operation. Entries rejected by the common model dropping policy are not
//! copied, keeping the view consistent with factorization and column pricing.

const std = @import("std");
const matrix = @import("matrix");
const foundation = @import("foundation");

/// Row-major (CSR) transpose of the retained CSC entries. Built once per solve
/// and reused by every row-pricing operation without further allocation.
pub const RowView = struct {
    /// Allocator owning retained CSR buffers.
    allocator: std.mem.Allocator,
    /// CSR row offsets; active prefix length is `num_rows + 1`.
    row_starts: []usize = &.{},
    /// Structural column index parallel to `row_values`.
    row_columns: []u32 = &.{},
    /// Retained model coefficient parallel to `row_columns`.
    row_values: []f64 = &.{},
    /// Row insertion cursor used only while rebuilding CSR.
    row_cursor: []usize = &.{},
    /// Active row count of the last matrix passed to `build`.
    num_rows: usize = 0,
    /// Active structural column count of the last built matrix.
    num_cols: usize = 0,
    /// Number of retained entries in the active CSR prefix.
    nonzeros: usize = 0,

    /// Construct an empty reusable row workspace.
    pub fn init(allocator: std.mem.Allocator) RowView {
        return .{ .allocator = allocator };
    }

    /// Free all allocated buffers and reset to the empty state.
    pub fn deinit(self: *RowView) void {
        self.allocator.free(self.row_starts);
        self.allocator.free(self.row_columns);
        self.allocator.free(self.row_values);
        self.allocator.free(self.row_cursor);
        self.* = .{ .allocator = self.allocator };
    }

    /// Transpose the retained CSC entries into stable row-major order. Capacity
    /// growth is confined to this solve-boundary call.
    pub fn build(self: *RowView, csc: matrix.CscView) !void {
        try self.ensureCapacity(csc.num_rows, csc.values.len);
        const starts = self.row_starts[0 .. csc.num_rows + 1];
        const cursor = self.row_cursor[0..csc.num_rows];
        @memset(starts, 0);

        // First pass: count retained nonzeros per row to build the row pointers.
        for (csc.row_indices, csc.values) |row_id, coefficient| {
            if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
            starts[row_id.toUsize() + 1] += 1;
        }
        // Prefix-sum the per-row counts into CSR row pointers.
        for (0..csc.num_rows) |row| starts[row + 1] += starts[row];
        @memcpy(cursor, starts[0..csc.num_rows]);

        // Second pass: scatter (column, value) pairs into their row slots.
        for (0..csc.num_cols) |column| {
            const begin = csc.col_starts[column];
            const end = csc.col_starts[column + 1];
            for (csc.row_indices[begin..end], csc.values[begin..end]) |row_id, coefficient| {
                if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                const row = row_id.toUsize();
                const target = cursor[row];
                self.row_columns[target] = @intCast(column);
                self.row_values[target] = coefficient;
                cursor[row] += 1;
            }
        }
        self.num_rows = csc.num_rows;
        self.num_cols = csc.num_cols;
        self.nonzeros = starts[csc.num_rows];
    }

    /// Column indices for the nonzeros in `row`.
    pub fn rowColumns(self: *const RowView, row: usize) []const u32 {
        return self.row_columns[self.row_starts[row]..self.row_starts[row + 1]];
    }

    /// Values for the nonzeros in `row`.
    pub fn rowValues(self: *const RowView, row: usize) []const f64 {
        return self.row_values[self.row_starts[row]..self.row_starts[row + 1]];
    }

    /// Number of nonzeros in `row`.
    pub fn rowDegree(self: *const RowView, row: usize) usize {
        return self.row_starts[row + 1] - self.row_starts[row];
    }

    /// Total bytes currently held by all dynamic buffers (used for memory budgeting).
    pub fn requestedBytes(self: *const RowView) usize {
        return std.mem.sliceAsBytes(self.row_starts).len +
            std.mem.sliceAsBytes(self.row_columns).len +
            std.mem.sliceAsBytes(self.row_values).len +
            std.mem.sliceAsBytes(self.row_cursor).len;
    }

    /// Grow buffers as needed; only ever enlarges, never shrinks.
    fn ensureCapacity(self: *RowView, rows: usize, nonzeros: usize) !void {
        if (self.row_starts.len < rows + 1) self.row_starts = try self.allocator.realloc(self.row_starts, rows + 1);
        if (self.row_cursor.len < rows) self.row_cursor = try self.allocator.realloc(self.row_cursor, rows);
        if (self.row_columns.len < nonzeros) self.row_columns = try self.allocator.realloc(self.row_columns, nonzeros);
        if (self.row_values.len < nonzeros) self.row_values = try self.allocator.realloc(self.row_values, nonzeros);
    }
};

test "row view preserves retained entries in deterministic row order" {
    const row_indices = [_]foundation.RowId{ foundation.RowId.fromUsize(1), foundation.RowId.fromUsize(0), foundation.RowId.fromUsize(1) };
    const csc = matrix.CscView.initAssumeValid(2, 2, &.{ 0, 2, 3 }, &row_indices, &.{ 2.0, 3.0, 4.0 });
    var view = RowView.init(std.testing.allocator);
    defer view.deinit();
    try view.build(csc);
    try std.testing.expectEqualSlices(u32, &.{0}, view.rowColumns(0));
    try std.testing.expectEqualSlices(f64, &.{3.0}, view.rowValues(0));
    try std.testing.expectEqualSlices(u32, &.{ 0, 1 }, view.rowColumns(1));
    try std.testing.expectEqualSlices(f64, &.{ 2.0, 4.0 }, view.rowValues(1));
}
