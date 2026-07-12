//! Per-SOS constraint data.
//!
//! Each SOS entry stores its type and a name, while the variable
//! membership list is held in the packed `SosArray` (begin/indices/weights).

const std = @import("std");
const types = @import("../types.zig");

const SosType = types.SosType;

/// Fixed fields for one SOS constraint.
///
/// The member variable indices and weights are stored as packed
/// parallel slices in [`SosArray`](array.zig).
pub const SosData = struct {
    sos_type: SosType = .sos1,
    /// Heap-allocated name; freed by `SosArray` on removal / deinit.
    name: ?[]const u8 = null,

    pub fn deinit(data: SosData, allocator: std.mem.Allocator) void {
        if (data.name) |n| allocator.free(n);
    }

    pub fn clone(data: SosData, allocator: std.mem.Allocator) !SosData {
        return .{
            .sos_type = data.sos_type,
            .name = if (data.name) |n| try allocator.dupe(u8, n) else null,
        };
    }
};
