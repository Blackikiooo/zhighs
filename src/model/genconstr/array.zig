//! Managed container for general constraint data.
//!
//! SoA for fixed fields (type, resvar, name) + packed operand indices.

const std = @import("std");
const GenConstrData = @import("data.zig").GenConstrData;

/// Container for general constraints with packed operand index storage.
pub const GenConstrArray = struct {
    /// SoA for fixed per-GenConstr fields.
    inner: std.MultiArrayList(GenConstrData) = .{},

    /// Packed operand index lists (CSR-like: `begin[i] .. begin[i+1]` is
    /// the range into `indices` for GenConstr `i`).
    begin: std.ArrayListUnmanaged(usize) = .{},
    nvars: std.ArrayListUnmanaged(usize) = .{},
    indices: std.ArrayListUnmanaged(usize) = .{},

    pub inline fn len(self: GenConstrArray) usize {
        return self.inner.len;
    }

    pub inline fn get(self: GenConstrArray, index: usize) GenConstrData {
        return self.inner.get(index);
    }

    pub fn deinit(self: *GenConstrArray, allocator: std.mem.Allocator) void {
        for (self.inner.items(.name)) |n| if (n) |s| allocator.free(s);
        self.inner.deinit(allocator);
        self.begin.deinit(allocator);
        self.nvars.deinit(allocator);
        self.indices.deinit(allocator);
    }
};
