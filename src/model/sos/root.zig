//! SOS constraint entity module.

const std = @import("std");

pub const SosData = @import("data.zig").SosData;
pub const SosArray = @import("array.zig").SosArray;
pub const SOS = @import("index.zig").SOS;

test {
    std.testing.refAllDecls(@This());
}
