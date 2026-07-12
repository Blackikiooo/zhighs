//! Variable entity module.
//!
//! Provides the [`VarData`] type, the managed SoA container [`VarArray`],
//! and the typed index wrapper [`Var`].

const std = @import("std");

pub const VarData = @import("data.zig").VarData;
pub const VarArray = @import("array.zig").VarArray;
pub const Var = @import("index.zig").Var;

test {
    std.testing.refAllDecls(@This());
}
