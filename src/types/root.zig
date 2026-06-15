const std = @import("std");
const Int = @import("int.zig");
const Double = @import("double.zig");

pub const HighsInt = Int.HighsInt;
pub const HighsUInt = Int.HighsUInt;
pub const HighsDouble = Double.HighsDouble;

test {
    std.testing.refAllDecls(@This());
}
