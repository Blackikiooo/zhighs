//! Model, basis, option, and solution input/output.
//!
//! Planned formats include MPS and CPLEX LP, plus native solution and basis
//! formats. Parsing produces model builders rather than mutating solver state.

const std = @import("std");

test {
    std.testing.refAllDecls(@This());
}
