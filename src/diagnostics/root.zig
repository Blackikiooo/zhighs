//! Logging, timing, statistics, tracing, and debug consistency checks.
//!
//! Diagnostics observe solver state through narrow views and must not own
//! algorithm state or change numerical decisions.

const std = @import("std");

test {
    std.testing.refAllDecls(@This());
}
