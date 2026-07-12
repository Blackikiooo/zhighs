//! Managed SoA array of variable data.
//!
//! Wraps [`std.MultiArrayList(VarData)`] and automatically frees
//! heap-allocated `.name` fields on removal or deinitialisation.

const std = @import("std");

const VarData = @import("data.zig").VarData;

/// SoA container for variable data with automatic name lifecycle.
pub const VarArray = struct {
    inner: std.MultiArrayList(VarData) = .{},

    // ── Query ─────────────────────────────────────────────────────────────

    pub inline fn len(self: VarArray) usize {
        return self.inner.len;
    }

    pub inline fn capacity(self: VarArray) usize {
        return self.inner.capacity;
    }

    /// Return the full `VarData` for one variable (heap copy of all fields).
    pub inline fn get(self: VarArray, index: usize) VarData {
        return self.inner.get(index);
    }

    /// Overwrite all fields of a variable.
    ///
    /// **Caller is responsible** for freeing the old `.name` before calling
    /// `set` if the name is being replaced.
    pub inline fn set(self: *VarArray, index: usize, data: VarData) void {
        self.inner.set(index, data);
    }

    // ── Growth ────────────────────────────────────────────────────────────

    /// Reserve capacity for at least `n` additional variables.
    pub inline fn ensureUnusedCapacity(self: *VarArray, allocator: std.mem.Allocator, n: usize) !void {
        try self.inner.ensureUnusedCapacity(allocator, n);
    }

    /// Reserve capacity for at least `n` variables total.
    pub inline fn ensureTotalCapacity(self: *VarArray, allocator: std.mem.Allocator, n: usize) !void {
        try self.inner.ensureTotalCapacity(allocator, n);
    }

    /// Append one variable (takes ownership of `.name`).
    pub inline fn append(self: *VarArray, allocator: std.mem.Allocator, data: VarData) !void {
        try self.inner.append(allocator, data);
    }

    /// Append one variable, assuming capacity is already reserved.
    pub inline fn appendAssumeCapacity(self: *VarArray, data: VarData) void {
        self.inner.appendAssumeCapacity(data);
    }

    // ── Removal ───────────────────────────────────────────────────────────

    /// Remove the variable at `index` by swapping in the last element.
    /// The removed variable's `.name` is freed automatically.
    pub fn swapRemove(self: *VarArray, allocator: std.mem.Allocator, index: usize) void {
        const removed = self.inner.get(index);
        if (removed.name) |n| allocator.free(n);
        _ = self.inner.swapRemove(index);
    }

    // ── Bulk ──────────────────────────────────────────────────────────────

    /// Deep-clone every variable (names are duplicated).
    pub fn clone(self: VarArray, allocator: std.mem.Allocator) !VarArray {
        var new: VarArray = .{};
        try new.inner.ensureTotalCapacity(allocator, self.inner.len);
        for (0..self.inner.len) |i| {
            const data = try self.inner.get(i).clone(allocator);
            new.inner.appendAssumeCapacity(data);
        }
        return new;
    }

    // ── Lifecycle ─────────────────────────────────────────────────────────

    /// Free all `.name` strings and release the backing memory.
    pub fn deinit(self: *VarArray, allocator: std.mem.Allocator) void {
        for (self.inner.items(.name)) |n| if (n) |s| allocator.free(s);
        self.inner.deinit(allocator);
    }
};
