//! First-order primal-dual LP solver support.
//!
//! Planned scope includes PDHG iterations, scaling, adaptive steps, restart,
//! termination checks, and optional accelerator backends.

const std = @import("std");

test {
    std.testing.refAllDecls(@This());
}
