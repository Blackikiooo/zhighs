//! SCIP-inspired lifecycle and component scheduler.
//!
//! Owns stages, registries, scheduling metadata, events, explicit component
//! results, and restricted services exposed to plugins.

const std = @import("std");

pub const Stage = @import("stage.zig").Stage;

test {
    std.testing.refAllDecls(@This());
}
