//! Basis import, export, and repair methods for `SimplexEngine`.
//!
//! ## Responsibility
//!
//! Owns the `BasisView`/`BasisSnapshot` round trip, feasibility classification
//! of imported bases, and transactional repair of rank-deficient bases.

const std = @import("std");
const basis_module = @import("basis.zig");
const basis_snapshot_module = @import("basis_snapshot.zig");
const factorization_module = @import("factorization.zig");
const pricing_module = @import("pricing.zig");
const ratio_module = @import("ratio_test.zig");
const numerical_module = @import("numerical.zig");
const dual_phase_one_module = @import("dual_phase_one.zig");
const crash_module = @import("crash.zig");
const degeneracy_module = @import("degeneracy.zig");
const pricing_workspace_module = @import("pricing_workspace.zig");
const problem_module = @import("problem.zig");
const solution_module = @import("solution.zig");
const foundation = @import("foundation");
const matrix = @import("matrix");
const SimplexEngine = @import("engine.zig").SimplexEngine;
const Algorithm = @import("engine.zig").Algorithm;
const SolveStatus = @import("engine.zig").SolveStatus;
const BasisImportError = @import("engine.zig").BasisImportError;

/// Restore a validated borrowed basis and rebuild factorization/basic
/// values without retaining any caller memory.
pub fn importBasis(self: *SimplexEngine, problem: problem_module.ProblemView, view: basis_snapshot_module.BasisView) BasisImportError!void {
    try view.validate(problem.num_cols, problem.num_rows);
    const basis = if (self.basis) |*value| value else return error.NumericalFailure;
    @memcpy(basis.col_status[0..problem.num_cols], view.structural_status);
    @memcpy(basis.col_status[problem.num_cols..][0..problem.num_rows], view.logical_status);
    @memset(basis.col_status[problem.num_cols + problem.num_rows ..], .fixed);

    for (basis.col_status[0 .. problem.num_cols + problem.num_rows], basis.col_lower[0 .. problem.num_cols + problem.num_rows], basis.col_upper[0 .. problem.num_cols + problem.num_rows]) |status, lower, upper| {
        const valid = switch (status) {
            .basic => true,
            .at_lower => std.math.isFinite(lower),
            .at_upper => std.math.isFinite(upper),
            .fixed => std.math.isFinite(lower) and lower == upper,
            .free => !std.math.isFinite(lower) and !std.math.isFinite(upper),
            .superbasic => 0.0 >= lower and 0.0 <= upper,
        };
        if (!valid) return error.InvalidNonbasicStatus;
    }
    @memset(basis.basic_pos, std.math.maxInt(u32));
    @memcpy(basis.basic_index, view.basic_index);
    for (basis.basic_index, 0..) |column, row| basis.basic_pos[column] = @intCast(row);

    for (basis.col_status[0 .. problem.num_cols + problem.num_rows], 0..) |status, column| {
        basis.primal[column] = switch (status) {
            .at_lower, .fixed => basis.col_lower[column],
            .at_upper => basis.col_upper[column],
            .free, .superbasic, .basic => 0.0,
        };
    }
    for (basis.basic_index, 0..) |column, row| {
        basis.basic_lower[row] = basis.col_lower[column];
        basis.basic_upper[row] = basis.col_upper[column];
    }
    self.factorizeCurrentBasis(problem) catch |err| switch (err) {
        error.Singular => try self.repairRankDeficientBasis(problem),
        error.OutOfMemory => return error.OutOfMemory,
        error.DimensionMismatch, error.NotImplemented, error.NumericalFailure => return error.NumericalFailure,
    };
    if (self.finishRefactorization() != .optimal) return error.NumericalFailure;
    if (self.recomputeBasicValues(problem) != .optimal) {
        // A valid dual warm start is allowed to be primal infeasible, so
        // rebuild values without enforcing bounds.
        if (self.recomputeBasicValuesUnchecked(problem) != .optimal) return error.NumericalFailure;
    }
}

pub fn exportBasisView(self: *const SimplexEngine, problem: problem_module.ProblemView) ?basis_snapshot_module.BasisView {
    const basis = if (self.basis) |*value| value else return null;
    const artificial_begin = problem.num_cols + problem.num_rows;
    for (basis.basic_index) |column| if (column >= artificial_begin) return null;
    return .{
        .structural_status = basis.col_status[0..problem.num_cols],
        .logical_status = basis.col_status[problem.num_cols..][0..problem.num_rows],
        .basic_index = basis.basic_index,
    };
}

pub fn exportBasisSnapshot(self: *const SimplexEngine, allocator: std.mem.Allocator, problem: problem_module.ProblemView) BasisImportError!basis_snapshot_module.BasisSnapshot {
    const basis_view = self.exportBasisView(problem) orelse return error.NumericalFailure;
    return basis_snapshot_module.BasisSnapshot.initFromView(allocator, basis_view);
}

pub fn classifyFeasibility(self: *const SimplexEngine, problem: problem_module.ProblemView) struct { primal: bool, dual: bool } {
    const basis = if (self.basis) |*value| value else return .{ .primal = false, .dual = false };
    var primal = true;
    for (basis.basic_value, basis.basic_lower, basis.basic_upper) |value, lower, upper| {
        if (value < lower - self.numerical.primal_tolerance or value > upper + self.numerical.primal_tolerance) {
            primal = false;
            break;
        }
    }
    var dual = true;
    for (basis.reduced_cost[0 .. problem.num_cols + problem.num_rows], basis.col_status[0 .. problem.num_cols + problem.num_rows]) |reduced, status| {
        const infeasible = switch (status) {
            .at_lower => reduced < -self.numerical.dual_tolerance,
            .at_upper => reduced > self.numerical.dual_tolerance,
            .free, .superbasic => @abs(reduced) > self.numerical.dual_tolerance,
            .basic, .fixed => false,
        };
        if (infeasible) {
            dual = false;
            break;
        }
    }
    return .{ .primal = primal, .dual = dual };
}

/// Repair a singular imported basis by cumulatively replacing structural
/// basics with currently nonbasic logical columns. Unique logical columns
/// form an identity basis, so this deterministic process must recover a
/// nonsingular basis unless factorization fails for a non-rank reason.
pub fn repairRankDeficientBasis(self: *SimplexEngine, problem: problem_module.ProblemView) BasisImportError!void {
    const basis = if (self.basis) |*value| value else return error.NumericalFailure;
    const logical_begin = problem.num_cols;
    const logical_end = logical_begin + problem.num_rows;
    const maximum_incremental_trials = 8;
    var incremental_trials: usize = 0;
    for (0..problem.num_rows) |leaving_row| {
        const leaving_column: usize = basis.basic_index[leaving_row];
        if (leaving_column >= problem.num_cols) continue;

        var entering_column = logical_begin;
        while (entering_column < logical_end and basis.col_status[entering_column] == .basic) : (entering_column += 1) {}
        if (entering_column == logical_end) return error.SingularBasis;

        const leaving_status = nonbasicStatusForBounds(basis.col_lower[leaving_column], basis.col_upper[leaving_column]);
        basis.primal[leaving_column] = switch (leaving_status) {
            .at_upper => basis.col_upper[leaving_column],
            .at_lower, .fixed => basis.col_lower[leaving_column],
            .free, .superbasic => 0.0,
            .basic => unreachable,
        };
        basis.primal[entering_column] = 0.0;
        basis.applyPivot(leaving_row, entering_column, leaving_status) catch return error.NumericalFailure;
        basis.basic_lower[leaving_row] = basis.col_lower[entering_column];
        basis.basic_upper[leaving_row] = basis.col_upper[entering_column];
        self.rank_repair_count += 1;

        // Small deficiencies usually recover after one replacement. Cap
        // repeated INVERT trials; severe deficiency then falls through to
        // one final logical-basis factorization instead of O(n) retries.
        if (incremental_trials < maximum_incremental_trials) {
            incremental_trials += 1;
            self.factorizeCurrentBasis(problem) catch |err| switch (err) {
                error.Singular => continue,
                error.OutOfMemory => return error.OutOfMemory,
                error.DimensionMismatch, error.NotImplemented, error.NumericalFailure => return error.NumericalFailure,
            };
            return;
        }
    }
    self.factorizeCurrentBasis(problem) catch |err| return switch (err) {
        error.Singular => error.SingularBasis,
        error.OutOfMemory => error.OutOfMemory,
        error.DimensionMismatch, error.NotImplemented, error.NumericalFailure => error.NumericalFailure,
    };
}

pub fn nonbasicStatusForBounds(lower: f64, upper: f64) basis_module.BasisStatus {
    if (std.math.isFinite(lower) and std.math.isFinite(upper) and lower == upper) return .fixed;
    if (std.math.isFinite(lower)) return .at_lower;
    if (std.math.isFinite(upper)) return .at_upper;
    return .free;
}

test "imported dual-feasible basis reoptimizes a changed RHS" {
    const rows = [_]foundation.RowId{foundation.RowId.fromUsizeAssumeValid(0)};
    const original = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{1},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{2},
        .row_upper = &[_]f64{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &[_]usize{ 0, 1 }, &rows, &[_]f64{1}),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var first_engine = SimplexEngine.init(std.testing.allocator);
    defer first_engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, first_engine.solveProblem(original, .{}));
    const exported = first_engine.exportBasisView(original).?;
    var snapshot = try basis_snapshot_module.BasisSnapshot.initFromView(std.testing.allocator, exported);
    defer snapshot.deinit();

    const modified = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = original.col_cost,
        .col_lower = original.col_lower,
        .col_upper = original.col_upper,
        .row_lower = &[_]f64{-1},
        .row_upper = original.row_upper,
        .matrix = original.matrix,
        .objective_sense = original.objective_sense,
        .objective_offset = original.objective_offset,
    };
    var second_engine = SimplexEngine.init(std.testing.allocator);
    defer second_engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, second_engine.solveProblem(modified, .{ .initial_basis = snapshot.view() }));
    try std.testing.expectEqual(Algorithm.dual_revised, second_engine.algorithm);
    try std.testing.expectApproxEqAbs(@as(f64, 0), second_engine.basis.?.primal[0], 1e-12);
    try std.testing.expect(second_engine.stats.dual_reduced_cost_updates > 0);

    // Exact refresh replaces a deliberately drifted incremental value and
    // records the normalized discrepancy before terminal publication.
    second_engine.basis.?.reduced_cost[0] += 1e-4;
    try std.testing.expectEqual(SolveStatus.optimal, second_engine.recomputeReducedCostsWithDrift(modified));
    try std.testing.expect(second_engine.maximum_reduced_cost_drift >= 1e-4);
    try std.testing.expectEqual(@as(usize, 1), second_engine.stats.dual_exact_reprices);

    // The explicit steepest-devex lifecycle starts with exact DSE and uses a
    // one-update budget here to force a deterministic dual Devex fallback.
    var fallback_engine = SimplexEngine.init(std.testing.allocator);
    defer fallback_engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, fallback_engine.solveProblem(modified, .{
        .initial_basis = snapshot.view(),
        .dual_edge_weight_strategy = .steepest_devex,
        .dual_dse_update_budget = 1,
    }));
    try std.testing.expectEqual(Algorithm.dual_revised, fallback_engine.algorithm);
    try std.testing.expect(fallback_engine.stats.dual_dse_updates > 0);
    try std.testing.expectEqual(@as(usize, 1), fallback_engine.stats.dual_dse_budget_fallbacks);
}

test "singular imported basis is repaired without discarding independent structural basics" {
    const rows = [_]foundation.RowId{
        foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1),
        foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1),
    };
    const problem = problem_module.ProblemView{
        .num_rows = 2,
        .num_cols = 2,
        .col_cost = &[_]f64{ 0, 0 },
        .col_lower = &[_]f64{ 0, 0 },
        .col_upper = &[_]f64{ std.math.inf(f64), std.math.inf(f64) },
        .row_lower = &[_]f64{ 1, 1 },
        .row_upper = &[_]f64{ 1, 1 },
        .matrix = matrix.CscView.initAssumeValid(
            2,
            2,
            &[_]usize{ 0, 2, 4 },
            &rows,
            &[_]f64{ 1, 1, 1, 1 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    const initial_basis = basis_snapshot_module.BasisView{
        .structural_status = &[_]basis_module.BasisStatus{ .basic, .basic },
        .logical_status = &[_]basis_module.BasisStatus{ .fixed, .fixed },
        .basic_index = &[_]u32{ 0, 1 },
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();

    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{ .initial_basis = initial_basis }));
    try std.testing.expectEqual(@as(usize, 1), engine.rank_repair_count);
    var structural_basics: usize = 0;
    for (engine.basis.?.basic_index) |column| structural_basics += @intFromBool(column < problem.num_cols);
    try std.testing.expectEqual(@as(usize, 1), structural_basics);
    try std.testing.expect(engine.factorization.pivotConditionEstimate() < 1e12);
}

test "rank repair cumulatively replaces multiple dependent structural basics" {
    const rows = [_]foundation.RowId{
        foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(2),
        foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(2),
        foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(2),
    };
    const problem = problem_module.ProblemView{
        .num_rows = 3,
        .num_cols = 3,
        .col_cost = &[_]f64{ 0, 0, 0 },
        .col_lower = &[_]f64{ 0, 0, 0 },
        .col_upper = &[_]f64{ std.math.inf(f64), std.math.inf(f64), std.math.inf(f64) },
        .row_lower = &[_]f64{ 1, 1, 1 },
        .row_upper = &[_]f64{ 1, 1, 1 },
        .matrix = matrix.CscView.initAssumeValid(
            3,
            3,
            &[_]usize{ 0, 3, 6, 9 },
            &rows,
            &[_]f64{ 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    const initial_basis = basis_snapshot_module.BasisView{
        .structural_status = &[_]basis_module.BasisStatus{ .basic, .basic, .basic },
        .logical_status = &[_]basis_module.BasisStatus{ .fixed, .fixed, .fixed },
        .basic_index = &[_]u32{ 0, 1, 2 },
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();

    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{ .initial_basis = initial_basis }));
    try std.testing.expectEqual(@as(usize, 2), engine.rank_repair_count);
    var structural_basics: usize = 0;
    for (engine.basis.?.basic_index) |column| structural_basics += @intFromBool(column < problem.num_cols);
    try std.testing.expectEqual(@as(usize, 1), structural_basics);
}

test "rank repair restores a singular sparse backend basis" {
    const n = 64;
    const allocator = std.testing.allocator;
    const starts = try allocator.alloc(usize, n + 1);
    defer allocator.free(starts);
    const rows = try allocator.alloc(foundation.RowId, n);
    defer allocator.free(rows);
    const values = try allocator.alloc(f64, n);
    defer allocator.free(values);
    const zeros = try allocator.alloc(f64, n);
    defer allocator.free(zeros);
    const infinities = try allocator.alloc(f64, n);
    defer allocator.free(infinities);
    const structural_status = try allocator.alloc(basis_module.BasisStatus, n);
    defer allocator.free(structural_status);
    const logical_status = try allocator.alloc(basis_module.BasisStatus, n);
    defer allocator.free(logical_status);
    const basic_index = try allocator.alloc(u32, n);
    defer allocator.free(basic_index);
    for (0..n) |column| {
        starts[column] = column;
        rows[column] = foundation.RowId.fromUsizeAssumeValid(if (column == 1) 0 else column);
        values[column] = 1.0;
        basic_index[column] = @intCast(column);
    }
    starts[n] = n;
    @memset(zeros, 0.0);
    @memset(infinities, std.math.inf(f64));
    @memset(structural_status, .basic);
    @memset(logical_status, .fixed);

    const problem = problem_module.ProblemView{
        .num_rows = n,
        .num_cols = n,
        .col_cost = zeros,
        .col_lower = zeros,
        .col_upper = infinities,
        .row_lower = zeros,
        .row_upper = zeros,
        .matrix = matrix.CscView.initAssumeValid(n, n, starts, rows, values),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    const initial_basis = basis_snapshot_module.BasisView{
        .structural_status = structural_status,
        .logical_status = logical_status,
        .basic_index = basic_index,
    };
    var engine = SimplexEngine.init(allocator);
    defer engine.deinit();

    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{ .initial_basis = initial_basis }));
    try std.testing.expectEqual(factorization_module.BackendKind.sparse_lu, engine.factorization.backend_kind);
    try std.testing.expectEqual(@as(usize, 2), engine.rank_repair_count);
}
