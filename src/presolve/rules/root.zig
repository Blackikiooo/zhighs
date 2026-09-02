//! Built-in presolve reductions.
//!
//! Every reduction implemented here must also define how solution, basis,
//! certificates, and model status are restored during postsolve.

const std = @import("std");

pub const empty_row = @import("empty_row.zig");
pub const column_substitution = @import("column_substitution.zig");
pub const singleton_row = @import("singleton_row.zig");
pub const redundant_row = @import("redundant_row.zig");

test {
    std.testing.refAllDecls(@This());
}
