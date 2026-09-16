//! Interior-point LP solver.
//!
//! Mehrotra predictor-corrector primal-dual IPM over a reusable profile
//! normal-equations Cholesky kernel, with a dense reference backend.
//! `standard_form` rewrites a borrowed
//! `lp.ProblemView` into min ĉᵀx̂, Âx̂ = b̂, x̂ ≥ 0; `engine` iterates on that
//! form and maps the converged iterate back to model coordinates.

const std = @import("std");

pub const options = @import("options.zig");
pub const standard_form = @import("standard_form.zig");
pub const cholesky = @import("cholesky.zig");
pub const profile_cholesky = @import("profile_cholesky.zig");
pub const engine = @import("engine.zig");

pub const Options = options.Options;
pub const Tolerances = options.Tolerances;
pub const StandardForm = standard_form.StandardForm;
pub const Result = engine.Result;
pub const solve = engine.solve;
pub const solveWithStats = engine.solveWithStats;
pub const Stats = engine.Stats;

test {
    _ = @import("engine_tests.zig");
    std.testing.refAllDecls(@This());
}
