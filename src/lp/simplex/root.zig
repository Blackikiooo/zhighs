//! Revised-simplex algorithm family for sparse LPs.
//!
//! The package is intentionally split into policy/state modules so primal and
//! dual simplex, MIP re-optimization, and future parallel variants can share
//! the same basis and factorization infrastructure.

const std = @import("std");

pub const basis = @import("basis.zig");
pub const factorization = @import("factorization.zig");
pub const pricing = @import("pricing.zig");
pub const ratio_test = @import("ratio_test.zig");
pub const numerical = @import("numerical.zig");
pub const problem = @import("problem.zig");
pub const engine = @import("engine.zig");

test {
    std.testing.refAllDecls(@This());
}
