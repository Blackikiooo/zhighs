//! Non-owning writer for an exact-size contiguous string pool.
//!
//! The arena never allocates and never moves bytes. Its backing storage is
//! owned by `ModelData`, so every returned slice remains stable for the model's
//! lifetime. Parsers can borrow names from their input and publish them here in
//! one final pass, avoiding one allocation per model object.

const std = @import("std");

/// Failure returned when a caller-provided arena cannot hold the next string.
pub const Error = error{CapacityExceeded};

/// Monotonic writer over caller-owned contiguous byte storage.
pub const StringArena = struct {
    /// Complete writable backing region; this type never frees it.
    storage: []u8,
    /// Byte offset immediately after the retained prefix and appended strings.
    cursor: usize = 0,

    /// Start an empty arena over `storage`.
    pub fn init(storage: []u8) StringArena {
        return .{ .storage = storage };
    }

    /// Continue writing an existing pool after `used_bytes` retained bytes.
    pub fn fromUsed(storage: []u8, used_bytes: usize) Error!StringArena {
        if (used_bytes > storage.len) return error.CapacityExceeded;
        return .{ .storage = storage, .cursor = used_bytes };
    }

    /// Copy one string into the pool and return its stable borrowed slice.
    pub fn append(self: *StringArena, source: []const u8) Error![]u8 {
        const end = std.math.add(usize, self.cursor, source.len) catch return error.CapacityExceeded;
        if (end > self.storage.len) return error.CapacityExceeded;
        const destination = self.storage[self.cursor..end];
        @memcpy(destination, source);
        self.cursor = end;
        return destination;
    }

    /// Number of bytes occupied by the retained prefix and appended strings.
    pub inline fn used(self: StringArena) usize {
        return self.cursor;
    }

    /// Number of bytes available before the next append would fail.
    pub inline fn remaining(self: StringArena) usize {
        return self.storage.len - self.cursor;
    }
};

test "string arena returns stable contiguous slices and fills exactly" {
    var storage: [9]u8 = undefined;
    var arena = StringArena.init(&storage);
    const alpha = try arena.append("alpha");
    const beta = try arena.append("beta");

    try std.testing.expectEqualStrings("alpha", alpha);
    try std.testing.expectEqualStrings("beta", beta);
    try std.testing.expectEqual(@as(usize, 9), arena.used());
    try std.testing.expectEqual(@as(usize, 0), arena.remaining());
    try std.testing.expectEqual(@intFromPtr(alpha.ptr) + alpha.len, @intFromPtr(beta.ptr));
    try std.testing.expectError(error.CapacityExceeded, arena.append("x"));
}

test "string arena resumes after retained prefix and handles empty strings" {
    var storage: [8]u8 = undefined;
    @memcpy(storage[0..5], "model");
    var arena = try StringArena.fromUsed(&storage, 5);
    const empty = try arena.append("");
    const name = try arena.append("row");

    try std.testing.expectEqual(@intFromPtr(storage[5..].ptr), @intFromPtr(empty.ptr));
    try std.testing.expectEqualStrings("row", name);
    try std.testing.expectEqual(@as(usize, 8), arena.used());
    try std.testing.expectError(error.CapacityExceeded, StringArena.fromUsed(&storage, 9));
}
