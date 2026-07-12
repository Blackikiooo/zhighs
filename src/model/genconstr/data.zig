//! Per-general-constraint data.
//!
//! Each general constraint stores its type, the result variable index,
//! and a name.  The operand variable indices are held in the packed
//! `GenConstrArray` (nvars / indices).

const std = @import("std");
const types = @import("../types.zig");

const GenConstrType = types.GenConstrType;

/// Fixed fields for one general constraint.
pub const GenConstrData = struct {
    gc_type: GenConstrType = .max,
    resvar: usize = 0,
    name: ?[]const u8 = null,

    pub fn deinit(data: GenConstrData, allocator: std.mem.Allocator) void {
        if (data.name) |n| allocator.free(n);
    }

    pub fn clone(data: GenConstrData, allocator: std.mem.Allocator) !GenConstrData {
        return .{
            .gc_type = data.gc_type,
            .resvar = data.resvar,
            .name = if (data.name) |n| try allocator.dupe(u8, n) else null,
        };
    }
};
