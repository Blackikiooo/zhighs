//! MIP correctness regression suite.
//!
//! Solves a fixed set of small/medium instances with known optima under
//! several solver profiles and asserts status=optimal with the objective
//! matching the reference value (tolerance 1e-6). Any correctness regression
//! in presolve, branching, propagation, cuts, or the underlying LP session
//! fails this suite immediately.
//!
//! Run with:
//!   zig build test-mip-regression

const std = @import("std");
const io = @import("io");
const mip = @import("mip");
const plugin = @import("plugin");
const plugins_builtin = @import("plugins_builtin");

const Case = struct {
    path: []const u8,
    objective: f64,
};

/// Reference objectives were fixed by hand/brute force for the handcrafted
/// instances and cross-validated across solver profiles for the generated
/// ones.
const cases = [_]Case{
    // Binary knapsack-style MPS: maximize 3X + 2Y, 2X + 3Y <= 4.
    .{ .path = "test/model/data/tiny_mip.mps", .objective = 3.0 },
    // Presolve-reducible MPS with a fixed column.
    .{ .path = "test/model/data/presolve_mip.mps", .objective = 6.0 },
    // Random 22 binary / 14 covering rows (seed 7).
    .{ .path = "test/model/data/warm_restart_mip.lp", .objective = 63.0 },
    // 16 binary / 10 rows branching comparison instance.
    .{ .path = "bench/mip/branching_compare.lp", .objective = 44.0 },
    // Same model as warm_restart_mip, kept under bench/ for profiling.
    .{ .path = "bench/mip/random_22x14_seed7.lp", .objective = 63.0 },
    // 6-item binary knapsack; brute-force optimum 12+20+15 with weight 16.
    .{ .path = "test/mip/data/knapsack_6.lp", .objective = 47.0 },
    // 5-set binary set cover; brute-force optimum s1+s5.
    .{ .path = "test/mip/data/set_cover_5.lp", .objective = 7.0 },
    // Mixed binary/integer/continuous; brute-force optimum x=0,y=1,z=1.5.
    .{ .path = "test/mip/data/mixed_int_3.lp", .objective = 5.5 },
};

fn solveCase(allocator: std.mem.Allocator, case: Case, control: mip.SolveControl) !mip.SolveResult {
    var model_data = try io.readFile(std.testing.io, allocator, case.path, .{ .keep_names = false });
    defer model_data.deinit();
    const integrality = try allocator.alloc(mip.Integrality, model_data.col_type.len);
    defer allocator.free(integrality);
    for (model_data.col_type, integrality) |source, *target| target.* = switch (source) {
        .continuous => .continuous,
        .binary => .binary,
        .integer => .integer,
        .semi_continuous => .semi_continuous,
        .semi_integer => .semi_integer,
    };
    return mip.solve(allocator, .{
        .relaxation = .{
            .num_rows = model_data.row_lower.len,
            .num_cols = model_data.col_cost.len,
            .col_cost = model_data.col_cost,
            .col_lower = model_data.col_lower,
            .col_upper = model_data.col_upper,
            .row_lower = model_data.row_lower,
            .row_upper = model_data.row_upper,
            .matrix = model_data.matrix.view(),
            .objective_sense = if (model_data.objective_sense == .minimize) .minimize else .maximize,
            .objective_offset = model_data.objective_offset,
        },
        .integrality = integrality,
    }, control);
}

fn expectSuite(control: mip.SolveControl) !void {
    for (cases) |case| {
        var result = try solveCase(std.testing.allocator, case, control);
        defer result.deinit();
        try std.testing.expectEqual(mip.SolveStatus.optimal, result.status);
        try std.testing.expectApproxEqAbs(case.objective, result.objective_value, 1e-6);
        try std.testing.expectApproxEqAbs(case.objective, result.objective_bound, 1e-6);
    }
}

test "MIP regression suite under the default profile" {
    try expectSuite(.{});
}

test "MIP regression suite without cut separation" {
    try expectSuite(.{ .cut_separation = false });
}

test "MIP regression suite without row propagation" {
    try expectSuite(.{ .propagation = false });
}

test "MIP regression suite with most-fractional branching and no propagation" {
    // Regression for the warm dual no-entering defect: this profile on
    // random_22x14_seed7.lp used to abort with numerical_failure at node 41
    // because the exact-cost dual loop refused to certify a genuine
    // infeasibility once its single rebuild-refresh had been spent.
    var most_fractional = plugins_builtin.branching.MostFractional{};
    try expectSuite(.{
        .propagation = false,
        .reliability_branching = false,
        .components = .{
            .branching_rule = plugin.bindBranching(plugins_builtin.branching.MostFractional, &most_fractional),
        },
    });
}
