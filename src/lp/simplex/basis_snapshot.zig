//! Stable borrowed and owning simplex-basis representations.
//!
//! `BasisView` is a zero-copy lifetime-bound projection. `BasisSnapshot` owns
//! one contiguous status block plus a dense basis-head array, making snapshots
//! cheap to validate, move, cache, and restore without per-entity objects.

const std = @import("std");
const basis = @import("basis.zig");

pub const BasisStatus = basis.BasisStatus;
pub const BasisViewError = error{
    DimensionMismatch,
    InvalidColumn,
    DuplicateBasicColumn,
    InconsistentStatus,
    OutOfMemory,
};

/// Borrowed basis description. All slices must outlive the consuming call.
pub const BasisView = struct {
    structural_status: []const BasisStatus,
    logical_status: []const BasisStatus,
    /// Global internal column index for every basis row. Imported snapshots
    /// may reference structural and logical columns, never artificial columns.
    basic_index: []const u32,

    pub fn validate(self: BasisView, num_cols: usize, num_rows: usize) BasisViewError!void {
        if (self.structural_status.len != num_cols or self.logical_status.len != num_rows or self.basic_index.len != num_rows)
            return error.DimensionMismatch;
        const total_cols = num_cols + num_rows;

        var status_basic_count: usize = 0;
        for (self.structural_status) |status| if (status == .basic) {
            status_basic_count += 1;
        };
        for (self.logical_status) |status| if (status == .basic) {
            status_basic_count += 1;
        };
        if (status_basic_count != num_rows) return error.InconsistentStatus;

        for (self.basic_index, 0..) |column, row| {
            if (column >= total_cols) return error.InvalidColumn;
            for (self.basic_index[0..row]) |previous| {
                if (previous == column) return error.DuplicateBasicColumn;
            }
            const status = if (column < num_cols)
                self.structural_status[column]
            else
                self.logical_status[column - num_cols];
            if (status != .basic) return error.InconsistentStatus;
        }
    }
};

/// Owning, move-only basis snapshot.
pub const BasisSnapshot = struct {
    allocator: std.mem.Allocator,
    status_storage: []BasisStatus,
    structural_status: []BasisStatus,
    logical_status: []BasisStatus,
    basic_index: []u32,

    pub fn initFromView(allocator: std.mem.Allocator, input_view: BasisView) BasisViewError!BasisSnapshot {
        try input_view.validate(input_view.structural_status.len, input_view.logical_status.len);
        const status_storage = allocator.alloc(BasisStatus, input_view.structural_status.len + input_view.logical_status.len) catch return error.OutOfMemory;
        errdefer allocator.free(status_storage);
        const basic_index = allocator.dupe(u32, input_view.basic_index) catch return error.OutOfMemory;
        errdefer allocator.free(basic_index);

        const structural_status = status_storage[0..input_view.structural_status.len];
        const logical_status = status_storage[input_view.structural_status.len..];
        @memcpy(structural_status, input_view.structural_status);
        @memcpy(logical_status, input_view.logical_status);
        return .{
            .allocator = allocator,
            .status_storage = status_storage,
            .structural_status = structural_status,
            .logical_status = logical_status,
            .basic_index = basic_index,
        };
    }

    pub fn deinit(self: *BasisSnapshot) void {
        self.allocator.free(self.status_storage);
        self.allocator.free(self.basic_index);
        self.* = undefined;
    }

    pub fn view(self: *const BasisSnapshot) BasisView {
        return .{
            .structural_status = self.structural_status,
            .logical_status = self.logical_status,
            .basic_index = self.basic_index,
        };
    }
};

test "BasisView validates and BasisSnapshot owns one status block" {
    const structural = [_]BasisStatus{ .basic, .at_lower };
    const logical = [_]BasisStatus{.basic};
    const head = [_]u32{0};
    const view = BasisView{ .structural_status = &structural, .logical_status = &logical, .basic_index = &head };
    try std.testing.expectError(error.InconsistentStatus, view.validate(2, 1));

    const valid_logical = [_]BasisStatus{.at_lower};
    const valid = BasisView{ .structural_status = &structural, .logical_status = &valid_logical, .basic_index = &head };
    try valid.validate(2, 1);
    var snapshot = try BasisSnapshot.initFromView(std.testing.allocator, valid);
    defer snapshot.deinit();
    try std.testing.expectEqual(BasisStatus.basic, snapshot.view().structural_status[0]);
}

test "BasisView rejects duplicate basis-head columns" {
    const structural = [_]BasisStatus{ .basic, .basic };
    const logical = [_]BasisStatus{ .at_lower, .at_lower };
    const duplicate_head = [_]u32{ 0, 0 };
    const view = BasisView{
        .structural_status = &structural,
        .logical_status = &logical,
        .basic_index = &duplicate_head,
    };
    try std.testing.expectError(error.DuplicateBasicColumn, view.validate(2, 2));
}
