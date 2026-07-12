//! Constraint entity module.

const std = @import("std");

pub const ConstrData = @import("data.zig").ConstrData;
pub const ConstrArray = @import("array.zig").ConstrArray;
pub const Constr = @import("index.zig").Constr;

test {
    std.testing.refAllDecls(@This());
}
