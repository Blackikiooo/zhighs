//! Interior-point LP solver support.
//!
//! Planned scope includes KKT systems, scaling, iteration control, crossover,
//! and translation between IPM results and the common solution/basis model.

const std = @import("std");

test {
    std.testing.refAllDecls(@This());
}
