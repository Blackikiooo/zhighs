//! Model core test entry point, independent of API, solver, and presolve modules.
//!
//! Run with:
//!   zig build test-model
//!   zig build test-model -Dhighs-int-width=w64

const std = @import("std");
const model = @import("model");

test {
    // Force compilation of all exported declarations (and their inline tests).
    std.testing.refAllDecls(model);
}
