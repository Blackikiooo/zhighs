//! Built-in presolve reductions.
//!
//! Every reduction implemented here must also define how solution, basis,
//! certificates, and model status are restored during postsolve.

const std = @import("std");

test {
    std.testing.refAllDecls(@This());
}
