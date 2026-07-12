//! Matrix-only test entry point.
//!
//! Imports the full matrix module and runs all its inline tests, without
//! pulling in model, API, or the top-level zhighs module.
//!
//! This lets matrix performance work continue independently while higher
//! layers are being edited.

const matrix = @import("root.zig");
const std = @import("std");

test {
    std.testing.refAllDecls(matrix);
}
