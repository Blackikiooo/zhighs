//! Mutable simplex basis state.
//!
//! Hot numerical vectors are stored as separate contiguous arrays.  This
//! module owns basis membership and solution vectors, but not LU factors.

const std = @import("std");

pub const BasisStatus = enum(u8) { basic, at_lower, at_upper, superbasic, free };

pub const BasisState = struct {
    allocator: std.mem.Allocator,
    row_status: []BasisStatus = &.{},
    col_status: []BasisStatus = &.{},
    basic_index: []u32 = &.{},
    basic_pos: []u32 = &.{},
    primal: []f64 = &.{},
    dual: []f64 = &.{},
    reduced_cost: []f64 = &.{},

    pub fn init(allocator: std.mem.Allocator, rows: usize, cols: usize) !BasisState {
        var self = BasisState{
            .allocator = allocator,
            .row_status = try allocator.alloc(BasisStatus, rows),
            .col_status = try allocator.alloc(BasisStatus, cols),
            .basic_index = try allocator.alloc(u32, rows),
            .basic_pos = try allocator.alloc(u32, cols),
            .primal = try allocator.alloc(f64, cols),
            .dual = try allocator.alloc(f64, rows),
            .reduced_cost = try allocator.alloc(f64, cols),
        };
        errdefer self.deinit();
        @memset(self.row_status, .basic);
        @memset(self.col_status, .at_lower);
        @memset(self.basic_index, 0);
        @memset(self.basic_pos, 0);
        @memset(self.primal, 0.0);
        @memset(self.dual, 0.0);
        @memset(self.reduced_cost, 0.0);
        return self;
    }

    pub fn initializeSlackBasis(self: *BasisState) void {
        for (self.basic_index, 0..) |*col, row| col.* = @intCast(row);
        for (self.basic_pos, 0..) |*row, col| row.* = if (col < self.basic_index.len) @intCast(col) else std.math.maxInt(u32);
        @memset(self.row_status, .basic);
        @memset(self.col_status, .at_lower);
    }

    pub fn deinit(self: *BasisState) void {
        self.allocator.free(self.row_status);
        self.allocator.free(self.col_status);
        self.allocator.free(self.basic_index);
        self.allocator.free(self.basic_pos);
        self.allocator.free(self.primal);
        self.allocator.free(self.dual);
        self.allocator.free(self.reduced_cost);
    }
};

test {
    std.testing.refAllDecls(@This());
}
