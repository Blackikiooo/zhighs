//! Per-quadratic-constraint data.
//!
//! Each quadratic constraint has a quadratic part (packed Q-triplets)
//! and an optional linear part (packed L-triplets), stored in
//! the parallel slices of [`QConstrArray`](array.zig).

const std = @import("std");
const types = @import("../types.zig");

const Sense = types.Sense;

/// Fixed fields for one quadratic constraint.
pub const QConstrData = struct {
    sense: Sense = .less_equal,
    rhs: f64 = 0.0,
    name: ?[]const u8 = null,

    pub fn deinit(data: QConstrData, allocator: std.mem.Allocator) void {
        if (data.name) |n| allocator.free(n);
    }

    pub fn clone(data: QConstrData, allocator: std.mem.Allocator) !QConstrData {
        return .{
            .sense = data.sense,
            .rhs = data.rhs,
            .name = if (data.name) |n| try allocator.dupe(u8, n) else null,
        };
    }
};
