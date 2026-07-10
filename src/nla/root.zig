//! Numerical linear algebra used by optimization algorithms.
//!
//! Covers dense and sparse LU, basis factorization, FTRAN/BTRAN, update
//! strategies, singularity recovery, and reusable linear-system contracts.

const std = @import("std");

test {
    std.testing.refAllDecls(@This());
}
