//! Built-in implementations of the contracts in `plugin`.
//!
//! Algorithm implementations live here; registration, ordering, timing, and
//! lifecycle remain the responsibility of `framework`.

const std = @import("std");

test {
    std.testing.refAllDecls(@This());
}
