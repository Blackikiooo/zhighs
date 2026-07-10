//! Primal and dual revised simplex.
//!
//! This module will contain basis state, Phase I/II, pricing, ratio tests,
//! crash strategies, perturbation, reinversion, and warm-start support.

const std = @import("std");

test {
    std.testing.refAllDecls(@This());
}
