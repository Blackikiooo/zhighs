//! Shared writer utilities: deterministic names and buffered formatting.

const std = @import("std");
const types = @import("types.zig");

pub const Names = struct {
    allocator: std.mem.Allocator,
    columns: [][]u8,
    rows: [][]u8,

    pub fn init(allocator: std.mem.Allocator, view: types.ModelView) types.IoError!Names {
        if ((view.col_names.len != 0 and view.col_names.len != view.numCols()) or
            (view.row_names.len != 0 and view.row_names.len != view.numRows())) return error.InvalidDimensions;
        const columns = allocator.alloc([]u8, view.numCols()) catch return error.OutOfMemory;
        errdefer allocator.free(columns);
        var columns_init: usize = 0;
        errdefer for (columns[0..columns_init]) |name| allocator.free(name);
        var seen_columns = std.StringHashMap(void).init(allocator);
        defer seen_columns.deinit();
        for (0..view.numCols()) |index| {
            const retained_name = if (view.col_names.len == 0) null else view.col_names[index];
            columns[index] = if (retained_name) |name|
                allocator.dupe(u8, name) catch return error.OutOfMemory
            else
                std.fmt.allocPrint(allocator, "x{d}", .{index}) catch return error.OutOfMemory;
            columns_init += 1;
            if (seen_columns.contains(columns[index])) return error.DuplicateName;
            seen_columns.put(columns[index], {}) catch return error.OutOfMemory;
        }

        const rows = allocator.alloc([]u8, view.numRows()) catch return error.OutOfMemory;
        errdefer allocator.free(rows);
        var rows_init: usize = 0;
        errdefer for (rows[0..rows_init]) |name| allocator.free(name);
        var seen_rows = std.StringHashMap(void).init(allocator);
        defer seen_rows.deinit();
        for (0..view.numRows()) |index| {
            const retained_name = if (view.row_names.len == 0) null else view.row_names[index];
            rows[index] = if (retained_name) |name|
                allocator.dupe(u8, name) catch return error.OutOfMemory
            else
                std.fmt.allocPrint(allocator, "c{d}", .{index}) catch return error.OutOfMemory;
            rows_init += 1;
            if (seen_rows.contains(rows[index])) return error.DuplicateName;
            seen_rows.put(rows[index], {}) catch return error.OutOfMemory;
        }
        return .{ .allocator = allocator, .columns = columns, .rows = rows };
    }

    pub fn deinit(self: *Names) void {
        for (self.columns) |name| self.allocator.free(name);
        for (self.rows) |name| self.allocator.free(name);
        self.allocator.free(self.columns);
        self.allocator.free(self.rows);
        self.* = undefined;
    }
};

pub fn print(writer: *std.Io.Writer, allocator: std.mem.Allocator, comptime format: []const u8, args: anytype) types.IoError!void {
    _ = allocator; // Kept in the signature for writer-policy evolution.
    writer.print(format, args) catch return error.WriteFailed;
}

pub fn write(writer: *std.Io.Writer, text: []const u8) types.IoError!void {
    writer.writeAll(text) catch return error.WriteFailed;
}
