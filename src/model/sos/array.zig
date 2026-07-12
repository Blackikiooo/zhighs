//! Managed container for SOS constraint data.
//!
//! Uses SoA layout for the per-SOS fields (type, name) and holds
//! packed member lists (begin, indices, weights) externally.

const std = @import("std");
const types = @import("../types.zig");
const SosData = @import("data.zig").SosData;
const SosType = types.SosType;

/// Container for SOS constraints with packed member lists.
pub const SosArray = struct {
    /// SoA for fixed per-SOS fields.
    inner: std.MultiArrayList(SosData) = .{},

    /// Packed member lists (CSR-like: `begin[i] .. begin[i+1]` is the range
    /// into `indices` / `weights` for SOS `i`).
    begin: std.ArrayListUnmanaged(usize) = .{},
    indices: std.ArrayListUnmanaged(usize) = .{},
    weights: std.ArrayListUnmanaged(f64) = .{},

    pub inline fn len(self: SosArray) usize {
        return self.inner.len;
    }

    pub inline fn get(self: SosArray, index: usize) SosData {
        return self.inner.get(index);
    }

    /// Returns the range of member indices for SOS `i`.
    pub inline fn memberRange(self: SosArray, index: usize) struct { start: usize, end: usize } {
        return .{ .start = self.begin.items[index], .end = self.begin.items[index + 1] };
    }

    /// Reserve capacity for `n` additional SOS constraints.
    /// Ensures both the inner SoA and the begin-offset array have room.
    pub fn ensureUnusedCapacity(self: *SosArray, allocator: std.mem.Allocator, n: usize) !void {
        try self.inner.ensureUnusedCapacity(allocator, n);
        try self.begin.ensureUnusedCapacity(allocator, n + 1);
    }

    /// Append one SOS constraint with its member list.
    pub fn append(
        self: *SosArray,
        allocator: std.mem.Allocator,
        data: SosData,
        member_indices: []const usize,
        member_weights: ?[]const f64,
    ) !void {
        const start = if (self.begin.items.len > 0) self.begin.getLast() else 0;
        try self.inner.append(allocator, data);
        try self.begin.append(allocator, start);
        try self.indices.appendSlice(allocator, member_indices);
        if (member_weights) |w| {
            try self.weights.appendSlice(allocator, w);
        } else {
            try self.weights.appendNTimes(allocator, 1.0, member_indices.len);
        }
    }

    /// Remove the SOS at `index` (swap-remove).
    ///
    /// **Note:** The CSR-packed begin/indices/weights invariants are not yet
    /// adjusted; this is currently a placeholder.
    pub fn swapRemove(self: *SosArray, allocator: std.mem.Allocator, index: usize) void {
        const removed = self.inner.get(index);
        if (removed.name) |n| allocator.free(n);
        _ = self.inner.swapRemove(index);
    }

    /// Deep-clone every SOS (names and packed lists are duplicated).
    pub fn clone(self: SosArray, allocator: std.mem.Allocator) !SosArray {
        var new: SosArray = .{};
        try new.inner.ensureTotalCapacity(allocator, self.inner.len);
        try new.begin.ensureTotalCapacity(allocator, self.begin.items.len);
        try new.indices.ensureTotalCapacity(allocator, self.indices.items.len);
        try new.weights.ensureTotalCapacity(allocator, self.weights.items.len);
        for (0..self.inner.len) |i| {
            const data = try self.inner.get(i).clone(allocator);
            new.inner.appendAssumeCapacity(data);
        }
        new.begin.appendSliceAssumeCapacity(self.begin.items);
        new.indices.appendSliceAssumeCapacity(self.indices.items);
        new.weights.appendSliceAssumeCapacity(self.weights.items);
        return new;
    }

    pub fn deinit(self: *SosArray, allocator: std.mem.Allocator) void {
        for (self.inner.items(.name)) |n| if (n) |s| allocator.free(s);
        self.inner.deinit(allocator);
        self.begin.deinit(allocator);
        self.indices.deinit(allocator);
        self.weights.deinit(allocator);
    }
};
