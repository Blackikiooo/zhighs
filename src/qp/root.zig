//! Convex quadratic-programming solver support.
//!
//! Planned scope includes Hessian validation and an active-set solver with
//! reduced gradients, pricing, ratio tests, and reduced-Hessian operations.

const std = @import("std");

test {
    std.testing.refAllDecls(@This());
}
