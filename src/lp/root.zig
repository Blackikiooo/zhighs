//! Linear-programming engines and shared LP solve state.

const std = @import("std");

pub const simplex = @import("simplex/root.zig");

test {
    std.testing.refAllDecls(@This());
}
