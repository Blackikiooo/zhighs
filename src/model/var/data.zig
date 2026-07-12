//! Per-variable data.
//!
//! `VarData` holds the full state for one decision variable.
//! The owning `VarArray` stores these in SoA layout via
//! `std.MultiArrayList`.

const std = @import("std");
const types = @import("../types.zig");

const VarType = types.VarType;
const INFINITY = types.INFINITY;

/// All per-variable attributes stored in the model.
pub const VarData = struct {
    lb: f64 = 0.0,
    ub: f64 = INFINITY,
    obj: f64 = 0.0,
    v_type: VarType = .continuous,
    /// Heap-allocated name; freed by `VarArray` on removal / deinit.
    name: ?[]const u8 = null,

    /// Free the heap-allocated `name` if present.
    pub fn deinit(data: VarData, allocator: std.mem.Allocator) void {
        if (data.name) |n| allocator.free(n);
    }

    /// Deep-clone, duplicating the name string if present.
    pub fn clone(data: VarData, allocator: std.mem.Allocator) !VarData {
        return .{
            .lb = data.lb,
            .ub = data.ub,
            .obj = data.obj,
            .v_type = data.v_type,
            .name = if (data.name) |n| try allocator.dupe(u8, n) else null,
        };
    }
};
