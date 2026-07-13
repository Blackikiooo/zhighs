//! Quadratic constraint entity module.
//!
//! Exposes the typed [`QConstr`] handle. Quadratic and linear term storage
//! remains owned by `Model`; this module defines no parallel container.

const std = @import("std");

pub const QConstr = @import("index.zig").QConstr;

test {
    std.testing.refAllDecls(@This());
}
