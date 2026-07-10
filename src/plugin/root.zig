//! Runtime component contracts.
//!
//! Contracts cover constraint handlers, presolvers, separators, propagators,
//! branching rules, primal heuristics, conflict handlers, node selectors, cut
//! selectors, variable pricers, Benders components, IIS finders, event
//! handlers, readers, and optional relaxation handlers.

const std = @import("std");

pub const Kind = @import("kind.zig").Kind;

test {
    std.testing.refAllDecls(@This());
}
