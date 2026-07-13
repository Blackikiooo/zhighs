//! General constraint entity module.
//!
//! Exposes the typed [`GenConstr`] handle. General-constraint fields and their
//! packed encoding remain owned by `Model` and `model_genconstr.zig`.

const std = @import("std");

pub const GenConstr = @import("index.zig").GenConstr;

test {
    std.testing.refAllDecls(@This());
}
