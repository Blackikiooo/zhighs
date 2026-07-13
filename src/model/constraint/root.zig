//! Linear constraint entity module.
//!
//! Exposes the typed [`Constr`] handle. Constraint data remains owned by
//! `Model`; this module deliberately does not define a second storage container.

const std = @import("std");

pub const Constr = @import("index.zig").Constr;

test {
    std.testing.refAllDecls(@This());
}
