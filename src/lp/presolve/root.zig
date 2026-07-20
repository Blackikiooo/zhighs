//! Reversible LP presolve/postsolve module.
//!
//! Operates at the ProblemView level, independent of the simplex engine.
//! Each reduction records transformation data so postsolve can recover
//! exact original-coordinate solutions.

pub const presolve_mod = @import("presolve.zig");
pub const PresolvedProblem = presolve_mod.PresolvedProblem;
pub const PresolveResult = presolve_mod.PresolveResult;
pub const presolve = presolve_mod.presolve;
pub const postsolve = presolve_mod.postsolve;
