const std = @import("std");
const Int = @import("int.zig");
const Double = @import("double.zig");

pub const HInt = Int.HInt;
pub const HUInt = Int.HUInt;
pub const HD = Double.HD;

test {
    std.testing.refAllDecls(@This());
}
