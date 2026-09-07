//! Built-in implementations of the contracts in `plugin`.
//!
//! Algorithm implementations live here; registration, ordering, timing, and
//! lifecycle remain the responsibility of `framework`.

const std = @import("std");

pub const branching = @import("branching.zig");
pub const heuristics = @import("heuristics.zig");
pub const node_selectors = @import("node_selectors.zig");

test {
    std.testing.refAllDecls(@This());
}
