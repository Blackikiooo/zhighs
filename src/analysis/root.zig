//! Solve-time and post-solve analysis.
//!
//! Covers KKT checks, primal/dual rays, IIS, ranging, sensitivity data,
//! ill-conditioning reports, feasibility relaxation, and solution assessment.

const std = @import("std");

test {
    std.testing.refAllDecls(@This());
}
