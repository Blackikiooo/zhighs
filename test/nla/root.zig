//! NLA-only test entry point, independent of model and API modules.

const std = @import("std");
const nla = @import("nla");

test {
    std.testing.refAllDecls(nla);
}
