//! Per-constraint data.

const std = @import("std");
const types = @import("../types.zig");

const Sense = types.Sense;

/// All per-constraint attributes stored in the model.
pub const ConstrData = struct {
    sense: Sense = .less_equal,
    rhs: f64 = 0.0,
    /// Heap-allocated name; freed by `ConstrArray` on removal / deinit.
    name: ?[]const u8 = null,

    pub fn deinit(data: ConstrData, allocator: std.mem.Allocator) void {
        if (data.name) |n| allocator.free(n);
    }

    pub fn clone(data: ConstrData, allocator: std.mem.Allocator) !ConstrData {
        return .{
            .sense = data.sense,
            .rhs = data.rhs,
            .name = if (data.name) |n| try allocator.dupe(u8, n) else null,
        };
    }
};
