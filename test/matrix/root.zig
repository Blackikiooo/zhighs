//! Matrix-only test entry point, independent of model and API modules.

const std = @import("std");
const matrix = @import("matrix");

test {
    std.testing.refAllDecls(matrix);
}
