//! Optimization model and solution data.
//!
//! Covers LP/MIP/QP models, Hessians, multiple objectives, variable types,
//! names, solutions, bases, rays, and immutable model construction.

const std = @import("std");

test {
    std.testing.refAllDecls(@This());
}
