//! Stable Zig-facing API.
//!
//! This module will expose solver construction, model mutation, options,
//! status, callbacks, and result access without leaking internal state.

const std = @import("std");

test {
    std.testing.refAllDecls(@This());
}
