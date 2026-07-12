//! General constraint entity module.

const std = @import("std");

pub const GenConstrData = @import("data.zig").GenConstrData;
pub const GenConstrArray = @import("array.zig").GenConstrArray;
pub const GenConstr = @import("index.zig").GenConstr;

test {
    std.testing.refAllDecls(@This());
}
