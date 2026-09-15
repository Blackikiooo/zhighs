//! Reversible LP presolve and postsolve.
//!
//! Operates at the ProblemView level, independent of the simplex engine.
//! Each reduction records transformation data so postsolve recovers exact
//! original-coordinate solutions.

const std = @import("std");

pub const rules = @import("rules/root.zig");
pub const presolve_mod = @import("presolve.zig");
pub const PresolvedProblem = presolve_mod.PresolvedProblem;
pub const PresolvePipeline = presolve_mod.PresolvePipeline;
pub const PresolveResult = presolve_mod.PresolveResult;
pub const ColumnNature = presolve_mod.ColumnNature;
pub const presolve = presolve_mod.presolve;
pub const presolveMip = presolve_mod.presolveMip;
pub const presolvePipeline = PresolvePipeline.init;
pub const postsolve = presolve_mod.postsolve;

test {
    std.testing.refAllDecls(@This());
}
