//! SOS constraint entity module.
//!
//! Exposes the typed [`SOS`] handle. SOS fields and packed membership data
//! remain owned by `Model`; this module defines no parallel storage container.

const std = @import("std");

pub const SOS = @import("index.zig").SOS;

test {
    std.testing.refAllDecls(@This());
}
