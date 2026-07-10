//! Top-level solve orchestration.
//!
//! Owns the incumbent model and coordinates validation, presolve, solver
//! selection, postsolve, callbacks, limits, and final status publication.

const std = @import("std");

test {
    std.testing.refAllDecls(@This());
}
