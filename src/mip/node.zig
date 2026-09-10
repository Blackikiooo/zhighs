//! Owned branch-and-bound node state.

const std = @import("std");
const BoundChange = @import("domain.zig").BoundChange;
const ChangeId = @import("domain.zig").ChangeId;
const DomainChangeArena = @import("domain.zig").DomainChangeArena;

pub const Node = struct {
    change_head: ?ChangeId = null,
    bound: ?f64 = null,
    depth: u32 = 0,

    pub fn childWith(self: Node, arena: *DomainChangeArena, allocator: std.mem.Allocator, bound: f64, change: BoundChange) !Node {
        return .{
            .change_head = try arena.append(allocator, self.change_head, change),
            .bound = bound,
            .depth = self.depth + 1,
        };
    }

    pub fn deinit(self: *Node, _: std.mem.Allocator) void {
        self.* = undefined;
    }
};

test "child nodes own compact accumulated domain changes" {
    var arena = DomainChangeArena{};
    defer arena.deinit(std.testing.allocator);
    var root = Node{};
    defer root.deinit(std.testing.allocator);
    var left = try root.childWith(&arena, std.testing.allocator, 1.5, .{ .column = 3, .upper = 7.0 });
    defer left.deinit(std.testing.allocator);
    var grandchild = try left.childWith(&arena, std.testing.allocator, 2.0, .{ .column = 8, .lower = 4.0 });
    defer grandchild.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 2), arena.entries.items.len);
    try std.testing.expect(root.change_head == null);
    try std.testing.expect(left.change_head != null);
    try std.testing.expect(grandchild.change_head != null);
}
