//! Typed plugin registry with SCIP-style priority ordering.
//!
//! Plugins are registered with a signed priority; within each contract kind
//! the registry keeps entries sorted by descending priority (ties keep
//! registration order). Engines consult branching rules and heuristics in that
//! order and use the highest-priority node selector as primary.

const std = @import("std");
const plugin = @import("plugin");

pub const BranchingRuleEntry = struct {
    rule: plugin.BranchingRule,
    priority: i32,
};

pub const HeuristicEntry = struct {
    heuristic: plugin.PrimalHeuristic,
    priority: i32,
};

pub const NodeSelectorEntry = struct {
    selector: plugin.NodeSelector,
    priority: i32,
};

/// Registry of solver components. Owns only the entry lists; the plugin
/// instances themselves remain owned by their providers and must outlive the
/// registry.
pub const PluginRegistry = struct {
    allocator: std.mem.Allocator,
    branching_rules: std.ArrayListUnmanaged(BranchingRuleEntry) = .empty,
    heuristics: std.ArrayListUnmanaged(HeuristicEntry) = .empty,
    node_selectors: std.ArrayListUnmanaged(NodeSelectorEntry) = .empty,

    pub fn init(allocator: std.mem.Allocator) PluginRegistry {
        return .{ .allocator = allocator };
    }

    pub fn deinit(self: *PluginRegistry) void {
        self.branching_rules.deinit(self.allocator);
        self.heuristics.deinit(self.allocator);
        self.node_selectors.deinit(self.allocator);
        self.* = undefined;
    }

    /// Register a branching rule; higher-priority rules are consulted earlier.
    pub fn addBranchingRule(self: *PluginRegistry, rule: plugin.BranchingRule, priority: i32) std.mem.Allocator.Error!void {
        const slot = insertionIndex(BranchingRuleEntry, self.branching_rules.items, priority);
        try self.branching_rules.insert(self.allocator, slot, .{ .rule = rule, .priority = priority });
    }

    /// Register a primal heuristic; higher-priority heuristics run earlier.
    pub fn addHeuristic(self: *PluginRegistry, heuristic: plugin.PrimalHeuristic, priority: i32) std.mem.Allocator.Error!void {
        const slot = insertionIndex(HeuristicEntry, self.heuristics.items, priority);
        try self.heuristics.insert(self.allocator, slot, .{ .heuristic = heuristic, .priority = priority });
    }

    /// Register a node selector; the highest-priority entry becomes primary.
    pub fn addNodeSelector(self: *PluginRegistry, selector: plugin.NodeSelector, priority: i32) std.mem.Allocator.Error!void {
        const slot = insertionIndex(NodeSelectorEntry, self.node_selectors.items, priority);
        try self.node_selectors.insert(self.allocator, slot, .{ .selector = selector, .priority = priority });
    }

    pub fn branchingRuleEntries(self: *const PluginRegistry) []const BranchingRuleEntry {
        return self.branching_rules.items;
    }

    pub fn heuristicEntries(self: *const PluginRegistry) []const HeuristicEntry {
        return self.heuristics.items;
    }

    pub fn nodeSelectorEntries(self: *const PluginRegistry) []const NodeSelectorEntry {
        return self.node_selectors.items;
    }

    /// The node selector the engine should dispatch through, if any.
    pub fn primaryNodeSelector(self: *const PluginRegistry) ?plugin.NodeSelector {
        if (self.node_selectors.items.len == 0) return null;
        return self.node_selectors.items[0].selector;
    }

    fn insertionIndex(comptime E: type, entries: []const E, priority: i32) usize {
        var index: usize = 0;
        while (index < entries.len and entries[index].priority >= priority) index += 1;
        return index;
    }
};

const StubBranching = struct {
    id: usize = 0,
    pub const name = "stub-branching";
    pub fn select(_: *StubBranching, _: plugin.BranchingContext) ?plugin.BranchingCandidate {
        return null;
    }
};

const StubHeuristic = struct {
    id: usize = 0,
    pub const name = "stub-heuristic";
    pub fn trySolution(_: *StubHeuristic, _: plugin.HeuristicContext) std.mem.Allocator.Error!?[]f64 {
        return null;
    }
};

const StubSelector = struct {
    id: usize = 0,
    pub const name = "stub-selector";
    pub fn select(_: *StubSelector, _: plugin.NodeSelectionContext) ?usize {
        return null;
    }
};

test "registry orders plugins by descending priority" {
    var low = StubBranching{};
    var high = StubBranching{};
    var mid = StubBranching{};
    var registry = PluginRegistry.init(std.testing.allocator);
    defer registry.deinit();
    try registry.addBranchingRule(plugin.bindBranching(StubBranching, &low), 1);
    try registry.addBranchingRule(plugin.bindBranching(StubBranching, &high), 10);
    try registry.addBranchingRule(plugin.bindBranching(StubBranching, &mid), 5);
    const entries = registry.branchingRuleEntries();
    try std.testing.expectEqual(@as(i32, 10), entries[0].priority);
    try std.testing.expectEqual(@as(i32, 5), entries[1].priority);
    try std.testing.expectEqual(@as(i32, 1), entries[2].priority);
}

test "registry keeps registration order on priority ties" {
    var first = StubHeuristic{ .id = 1 };
    var second = StubHeuristic{ .id = 2 };
    var registry = PluginRegistry.init(std.testing.allocator);
    defer registry.deinit();
    try registry.addHeuristic(plugin.bindHeuristic(StubHeuristic, &first), 3);
    try registry.addHeuristic(plugin.bindHeuristic(StubHeuristic, &second), 3);
    const entries = registry.heuristicEntries();
    try std.testing.expectEqual(@as(i32, 3), entries[0].priority);
    try std.testing.expectEqualStrings("stub-heuristic", entries[0].heuristic.name());
    try std.testing.expectEqualStrings("stub-heuristic", entries[1].heuristic.name());
}

test "primary node selector is the highest-priority entry" {
    var secondary = StubSelector{ .id = 1 };
    var primary = StubSelector{ .id = 2 };
    var registry = PluginRegistry.init(std.testing.allocator);
    defer registry.deinit();
    try std.testing.expectEqual(@as(?plugin.NodeSelector, null), registry.primaryNodeSelector());
    try registry.addNodeSelector(plugin.bindNodeSelector(StubSelector, &secondary), -5);
    try registry.addNodeSelector(plugin.bindNodeSelector(StubSelector, &primary), 5);
    try std.testing.expectEqualStrings("stub-selector", registry.primaryNodeSelector().?.name());
}
