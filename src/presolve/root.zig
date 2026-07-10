//! Model presolve and reversible postsolve.

const std = @import("std");

pub const rules = @import("rules/root.zig");

test {
    std.testing.refAllDecls(@This());
}
