//! Variable entity module.
//!
//! Exposes the typed [`Var`] handle. Variable data remains owned by `Model`;
//! this module deliberately does not define a second storage container.

const std = @import("std");

pub const Var = @import("index.zig").Var;

test {
    std.testing.refAllDecls(@This());
}
