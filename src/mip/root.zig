//! Mixed-integer branch-and-cut framework.
//!
//! Covers LP relaxations, node search, domains, cut and conflict pools,
//! implications, pseudo-costs, primal heuristics, separation, and incumbents.

const std = @import("std");

test {
    std.testing.refAllDecls(@This());
}
