//! Optional parallel execution primitives.
//!
//! Planned scope includes task scheduling, synchronization, cache-aware data,
//! deterministic controls, and cancellation. Serial execution remains valid.

const std = @import("std");

test {
    std.testing.refAllDecls(@This());
}
