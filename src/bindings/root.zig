//! Stable foreign-function interfaces.
//!
//! The first target is a C ABI over `api`; language-specific bindings should
//! build on that ABI instead of reaching into solver internals.

const std = @import("std");

test {
    std.testing.refAllDecls(@This());
}
