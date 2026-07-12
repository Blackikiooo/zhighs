//! Quadratic constraint entity module.

const std = @import("std");

pub const QConstrData = @import("data.zig").QConstrData;
pub const QConstrArray = @import("array.zig").QConstrArray;
pub const QConstr = @import("index.zig").QConstr;

test {
    std.testing.refAllDecls(@This());
}
