//! SCIP-inspired lifecycle and component scheduler.
//!
//! Owns stages, registries, scheduling metadata, events, explicit component
//! results, and restricted services exposed to plugins.

const std = @import("std");

pub const Stage = @import("stage.zig").Stage;
pub const registry = @import("registry.zig");

pub const PluginRegistry = registry.PluginRegistry;
pub const BranchingRuleEntry = registry.BranchingRuleEntry;
pub const HeuristicEntry = registry.HeuristicEntry;
pub const NodeSelectorEntry = registry.NodeSelectorEntry;

test {
    std.testing.refAllDecls(@This());
}
