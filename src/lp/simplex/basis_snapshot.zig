//! Stable borrowed and owning simplex-basis representations.
//!
//! `BasisView` is a zero-copy lifetime-bound projection. `BasisSnapshot` owns
//! one contiguous status block plus a dense basis-head array, making snapshots
//! cheap to validate, move, cache, and restore without per-entity objects.

const std = @import("std");
const basis = @import("basis.zig");

pub const BasisStatus = basis.BasisStatus;

/// Errors that can arise when validating or cloning a basis.
pub const BasisViewError = error{
    DimensionMismatch, // Slice lengths do not match the model shape
    InvalidColumn, // A basis-head entry is out of range
    DuplicateBasicColumn, // The same column appears twice in the basis head
    InconsistentStatus, // Status array and basis head disagree
    OutOfMemory,
};

/// Borrowed basis description. All slices must outlive the consuming call.
pub const BasisView = struct {
    structural_status: []const BasisStatus, // Per structural column (length = num_cols)
    logical_status: []const BasisStatus, // Per logical column / row (length = num_rows)
    /// Global internal column index for every basis row. Imported snapshots
    /// may reference structural and logical columns, never artificial columns.
    basic_index: []const u32,

    /// Cross-check slice lengths, basic count, head range, and status consistency.
    pub fn validate(self: BasisView, num_cols: usize, num_rows: usize) BasisViewError!void {
        if (self.structural_status.len != num_cols or self.logical_status.len != num_rows or self.basic_index.len != num_rows)
            return error.DimensionMismatch;
        const total_cols = num_cols + num_rows;

        // The number of `basic` statuses must equal the number of basis rows.
        var status_basic_count: usize = 0;
        for (self.structural_status) |status| if (status == .basic) {
            status_basic_count += 1;
        };
        for (self.logical_status) |status| if (status == .basic) {
            status_basic_count += 1;
        };
        if (status_basic_count != num_rows) return error.InconsistentStatus;

        // Each head entry must be in range, unique, and marked basic in the status arrays.
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

/// Owning, move-only basis snapshot. The status arrays are stored in one
/// contiguous allocation to keep cache traffic tight during MIP re-optimization.
pub const BasisSnapshot = struct {
    allocator: std.mem.Allocator,
    status_storage: []BasisStatus, // Backing buffer for structural_status + logical_status
    structural_status: []BasisStatus, // Alias into status_storage[0..num_cols]
    logical_status: []BasisStatus, // Alias into status_storage[num_cols..]
    basic_index: []u32, // Dense list of basic column indices (length = num_rows)

    /// Deep-copy a validated `BasisView` into owning storage.
    pub fn initFromView(allocator: std.mem.Allocator, input_view: BasisView) BasisViewError!BasisSnapshot {
        try input_view.validate(input_view.structural_status.len, input_view.logical_status.len);
        const status_storage = try allocator.alloc(BasisStatus, input_view.structural_status.len + input_view.logical_status.len);
        errdefer allocator.free(status_storage);
        const basic_index = try allocator.dupe(u32, input_view.basic_index);
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

    /// Release the backing buffers.
    pub fn deinit(self: *BasisSnapshot) void {
        self.allocator.free(self.status_storage);
        self.allocator.free(self.basic_index);
        self.* = undefined;
    }

    /// Project this owning snapshot as a borrowed `BasisView`.
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
