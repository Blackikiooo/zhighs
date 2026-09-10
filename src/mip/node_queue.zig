//! Ownership-aware open-node container with a zero-overhead default DFS path.
//!
//! The built-in depth-first policy removes the newest node directly. A plugin
//! selector pays for materializing the intentionally restricted `OpenNode`
//! view only when such a selector is actually installed.

const std = @import("std");
const solver = @import("solver");
const plugin = @import("plugin");
const Node = @import("node.zig").Node;

pub const NodeQueue = struct {
    nodes: std.ArrayListUnmanaged(Node) = .empty,
    plugin_views: std.ArrayListUnmanaged(plugin.OpenNode) = .empty,

    pub fn deinit(self: *NodeQueue, allocator: std.mem.Allocator) void {
        for (self.nodes.items) |*node| node.deinit(allocator);
        self.nodes.deinit(allocator);
        self.plugin_views.deinit(allocator);
        self.* = undefined;
    }

    pub fn append(self: *NodeQueue, allocator: std.mem.Allocator, node: Node) std.mem.Allocator.Error!void {
        try self.nodes.append(allocator, node);
    }

    pub fn len(self: *const NodeQueue) usize {
        return self.nodes.items.len;
    }

    pub fn items(self: *const NodeQueue) []const Node {
        return self.nodes.items;
    }

    /// Take the next node. `null` selects the allocation-free default DFS
    /// policy; a plugin sees metadata only and cannot mutate queue ownership.
    pub fn take(
        self: *NodeQueue,
        allocator: std.mem.Allocator,
        sense: solver.LpObjectiveSense,
        selector: ?plugin.NodeSelector,
    ) std.mem.Allocator.Error!?Node {
        if (self.nodes.items.len == 0) return null;
        const policy = selector orelse return self.nodes.orderedRemove(self.nodes.items.len - 1);

        self.plugin_views.clearRetainingCapacity();
        try self.plugin_views.ensureTotalCapacity(allocator, self.nodes.items.len);
        for (self.nodes.items) |node| self.plugin_views.appendAssumeCapacity(.{
            .bound = node.bound,
            .depth = node.depth,
        });
        const selected = policy.select(.{
            .objective_sense = sense,
            .open = self.plugin_views.items,
        }) orelse return self.nodes.orderedRemove(self.nodes.items.len - 1);
        if (selected >= self.nodes.items.len) return self.nodes.orderedRemove(self.nodes.items.len - 1);
        return self.nodes.orderedRemove(selected);
    }
};

test "default DFS takes the newest node without allocating plugin views" {
    var queue = NodeQueue{};
    defer queue.deinit(std.testing.allocator);
    try queue.append(std.testing.allocator, .{
        .depth = 1,
    });
    try queue.append(std.testing.allocator, .{
        .depth = 2,
    });
    var newest = (try queue.take(std.testing.allocator, .minimize, null)).?;
    defer newest.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(u32, 2), newest.depth);
    try std.testing.expectEqual(@as(usize, 0), queue.plugin_views.capacity);
}

test "invalid plugin selections fall back to correctness-preserving DFS" {
    const InvalidSelector = struct {
        pub const name = "invalid-test-selector";
        pub fn select(_: *@This(), _: plugin.NodeSelectionContext) ?usize {
            return std.math.maxInt(usize);
        }
    };
    var selector_impl = InvalidSelector{};
    var queue = NodeQueue{};
    defer queue.deinit(std.testing.allocator);
    try queue.append(std.testing.allocator, .{
        .depth = 7,
    });
    var node = (try queue.take(
        std.testing.allocator,
        .minimize,
        plugin.bindNodeSelector(InvalidSelector, &selector_impl),
    )).?;
    defer node.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(u32, 7), node.depth);
}
