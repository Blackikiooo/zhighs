//! Managed SoA array of constraint data.

const std = @import("std");

const ConstrData = @import("data.zig").ConstrData;

/// SoA container for constraint data with automatic name lifecycle.
pub const ConstrArray = struct {
    inner: std.MultiArrayList(ConstrData) = .{},

    pub inline fn len(self: ConstrArray) usize { return self.inner.len; }
    pub inline fn capacity(self: ConstrArray) usize { return self.inner.capacity; }

    pub inline fn get(self: ConstrArray, index: usize) ConstrData {
        return self.inner.get(index);
    }

    pub inline fn set(self: *ConstrArray, index: usize, data: ConstrData) void {
        self.inner.set(index, data);
    }

    pub inline fn ensureUnusedCapacity(self: *ConstrArray, allocator: std.mem.Allocator, n: usize) !void {
        try self.inner.ensureUnusedCapacity(allocator, n);
    }

    pub inline fn ensureTotalCapacity(self: *ConstrArray, allocator: std.mem.Allocator, n: usize) !void {
        try self.inner.ensureTotalCapacity(allocator, n);
    }

    pub inline fn append(self: *ConstrArray, allocator: std.mem.Allocator, data: ConstrData) !void {
        try self.inner.append(allocator, data);
    }

    pub inline fn appendAssumeCapacity(self: *ConstrArray, data: ConstrData) void {
        self.inner.appendAssumeCapacity(data);
    }

    pub fn swapRemove(self: *ConstrArray, allocator: std.mem.Allocator, index: usize) void {
        const removed = self.inner.get(index);
        if (removed.name) |n| allocator.free(n);
        _ = self.inner.swapRemove(index);
    }

    pub fn clone(self: ConstrArray, allocator: std.mem.Allocator) !ConstrArray {
        var new: ConstrArray = .{};
        try new.inner.ensureTotalCapacity(allocator, self.inner.len);
        for (0..self.inner.len) |i| {
            const data = try self.inner.get(i).clone(allocator);
            new.inner.appendAssumeCapacity(data);
        }
        return new;
    }

    pub fn deinit(self: *ConstrArray, allocator: std.mem.Allocator) void {
        for (self.inner.items(.name)) |n| if (n) |s| allocator.free(s);
        self.inner.deinit(allocator);
    }
};
