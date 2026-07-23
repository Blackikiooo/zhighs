//! Problem-storage setup, scaling, and initial-basis installation for
//! `SimplexEngine`.
//!
//! ## Responsibility
//!
//! Translates the borrowed `ProblemView` into engine-owned basis storage,
//! computes row/column/objective scaling, and installs the logical, artificial
//! Phase-I, and sparse crash bases used at solve (re)start.

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
const SolveStatus = @import("engine.zig").SolveStatus;
const engine_basis_module = @import("engine_basis.zig");

/// Materialize the logical crash solution from current nonbasic structural
/// values. The basis remains the identity; values may violate logical
/// bounds so the caller can choose primal, dual, or artificial Phase I.
pub fn initializeLogicalBasicValues(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    @memcpy(basis.basic_value, basis.row_rhs);
    for (0..problem.num_cols) |column| {
        const initial = basis.primal[column];
        if (initial == 0.0) continue;
        const begin = problem.matrix.col_starts[column];
        const end = problem.matrix.col_starts[column + 1];
        for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
            if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
            const row_index = row.toUsize();
            basis.basic_value[row_index] -= basis.row_scale[row_index] * coefficient * basis.column_scale[column] * initial;
        }
    }
    for (basis.basic_index, basis.basic_value, 0..) |column, value, row| {
        if (!std.math.isFinite(value)) return .numerical_failure;
        basis.primal[column] = value;
        basis.basic_lower[row] = basis.col_lower[column];
        basis.basic_upper[row] = basis.col_upper[column];
    }
    return .optimal;
}

/// Replace each infeasible logical basic with its nonnegative artificial
/// column. Initial Phase I and every transactional restart must use this
/// exact state transition; rebuilding only the logical identity would
/// leave Phase-I costs disconnected from the violated rows.
pub fn installArtificialPhaseOneBasis(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    self.phase1_needed = false;
    for (basis.basic_value, basis.basic_lower, basis.basic_upper, 0..) |*value, lower, upper, row| {
        const logical_col = problem.num_cols + row;
        const artificial_col = problem.num_cols + problem.num_rows + row;
        if (!std.math.isFinite(value.*)) return .numerical_failure;
        if (value.* < lower - self.numerical.primal_tolerance) {
            basis.primal[logical_col] = lower;
            basis.col_status[logical_col] = .at_lower;
            basis.basic_pos[logical_col] = std.math.maxInt(u32);
            basis.artificial_sign[row] = -1.0;
            value.* = lower - value.*;
            basis.basic_index[row] = @intCast(artificial_col);
            basis.basic_pos[artificial_col] = @intCast(row);
            basis.col_status[artificial_col] = .basic;
            basis.col_upper[artificial_col] = std.math.inf(f64);
            basis.basic_lower[row] = 0.0;
            basis.basic_upper[row] = std.math.inf(f64);
            self.phase1_needed = true;
        } else if (value.* > upper + self.numerical.primal_tolerance) {
            basis.primal[logical_col] = upper;
            basis.col_status[logical_col] = .at_upper;
            basis.basic_pos[logical_col] = std.math.maxInt(u32);
            basis.artificial_sign[row] = 1.0;
            value.* -= upper;
            basis.basic_index[row] = @intCast(artificial_col);
            basis.basic_pos[artificial_col] = @intCast(row);
            basis.col_status[artificial_col] = .basic;
            basis.col_upper[artificial_col] = std.math.inf(f64);
            basis.basic_lower[row] = 0.0;
            basis.basic_upper[row] = std.math.inf(f64);
            self.phase1_needed = true;
        } else {
            basis.primal[logical_col] = value.*;
        }
    }
    return .optimal;
}

/// Install a sparse structural crash basis transactionally. The caller
/// restores the logical identity on false, so no partially matched basis
/// can escape numerical validation.
pub fn installSparseCrashBasis(
    self: *SimplexEngine,
    problem: problem_module.ProblemView,
    maximum_columns: ?usize,
    scoring: crash_module.CrashScoring,
) bool {
    const basis = if (self.basis) |*value| value else return false;
    const average_column_degree = if (problem.num_cols == 0) 0 else (problem.matrix.values.len + problem.num_cols - 1) / problem.num_cols;
    const near_singleton_limit: u32 = @intCast(@min(@max(average_column_degree * 2, 4), 32));
    const plan_view = self.crash.plan(
        problem.matrix,
        problem.col_cost,
        basis.col_lower[0..problem.num_cols],
        basis.col_upper[0..problem.num_cols],
        basis.row_scale,
        basis.column_scale[0..problem.num_cols],
        near_singleton_limit,
        self.numerical.pivot_tolerance,
        scoring,
    ) catch return false;
    if (plan_view.rows.len == 0) return false;
    self.stats.crash_planned_columns = plan_view.rows.len;

    // A near-singleton symbolic matching can still be numerically rank
    // deficient. Validate the longest deterministic prefix first, then
    // halve it until a safe basis is found. No partial basis is published.
    // Retain at least one quarter of the logical columns as stable anchors.
    // A fully structural crash can be perfectly invertible yet leave the
    // current dual Phase-I update path without a safe reinversion route.
    const anchored_limit = if (problem.num_rows == 0) 0 else @max(@as(usize, 1), problem.num_rows * 3 / 4);
    var prefix_count = @min(plan_view.rows.len, maximum_columns orelse anchored_limit);
    while (prefix_count != 0) : (prefix_count /= 2) {
        if (self.initializeProblemStorage(problem) != .optimal) return false;
        var basis_nonzeros = problem.num_rows - prefix_count;
        for (plan_view.rows[0..prefix_count], plan_view.columns[0..prefix_count]) |row_u32, column_u32| {
            const row: usize = @intCast(row_u32);
            const column: usize = @intCast(column_u32);
            const logical = problem.num_cols + row;
            const leaving_status = engine_basis_module.nonbasicStatusForBounds(basis.col_lower[logical], basis.col_upper[logical]);
            basis.primal[logical] = switch (leaving_status) {
                .at_lower, .fixed => basis.col_lower[logical],
                .at_upper => basis.col_upper[logical],
                .free, .superbasic => 0.0,
                .basic => unreachable,
            };
            basis.applyPivot(row, column, leaving_status) catch return false;
            basis.basic_lower[row] = basis.col_lower[column];
            basis.basic_upper[row] = basis.col_upper[column];
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            for (problem.matrix.values[begin..end]) |coefficient| {
                if (matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) basis_nonzeros += 1;
            }
        }
        self.factorizeCurrentBasis(problem) catch continue;
        if (self.finishRefactorization() != .optimal) continue;
        const condition = self.factorization.pivotConditionEstimate();
        if (!std.math.isFinite(condition) or condition > 1e12) continue;
        if (self.recomputeBasicValuesUnchecked(problem) != .optimal) continue;
        if (!std.math.isFinite(self.numerical.last_relative_residual) or
            self.numerical.last_relative_residual > self.numerical.residual_tolerance * 100.0)
            continue;
        self.stats.crash_structural_columns = prefix_count;
        self.stats.crash_basis_nonzeros = basis_nonzeros;
        self.stats.crash_condition_estimate = condition;
        return true;
    }
    return false;
}

/// Refresh bounds and nonbasic values while preserving the current basis/scaling.
///
/// Returns false when model dimensions changed or a retained nonbasic status is
/// incompatible with the new bounds, in which case the caller must cold-start.
pub fn refreshProblemStorage(self: *SimplexEngine, problem: problem_module.ProblemView) bool {
    const basis = if (self.basis) |*value| value else return false;
    if (basis.num_rows != problem.num_rows or basis.num_structural_cols != problem.num_cols) return false;
    for (problem.col_lower, problem.col_upper, 0..) |lower, upper, column| {
        if (lower > upper) return false;
        const scale = basis.column_scale[column];
        basis.col_lower[column] = lower / scale;
        basis.col_upper[column] = upper / scale;
        switch (basis.col_status[column]) {
            .basic => {},
            .at_lower => {
                if (!std.math.isFinite(lower)) return false;
                basis.primal[column] = lower / scale;
            },
            .at_upper => {
                if (!std.math.isFinite(upper)) return false;
                basis.primal[column] = upper / scale;
            },
            .fixed => {
                if (lower != upper) return false;
                basis.primal[column] = lower / scale;
            },
            .free, .superbasic => basis.primal[column] = @min(@max(0.0, lower), upper) / scale,
        }
    }
    for (problem.row_lower, problem.row_upper, 0..) |lower, upper, row| {
        const logical = problem.num_cols + row;
        if (self.active_dual_initialization_strategy == .highs) {
            const magnitude = basis.row_scale[row];
            if (!(magnitude > 0.0) or !std.math.isFinite(magnitude)) return false;
            basis.row_rhs[row] = 0.0;
            basis.col_lower[logical] = if (std.math.isFinite(upper)) -upper * magnitude else -std.math.inf(f64);
            basis.col_upper[logical] = if (std.math.isFinite(lower)) -lower * magnitude else std.math.inf(f64);
        } else {
            const sign: f64 = if (std.math.isFinite(upper)) 1.0 else if (std.math.isFinite(lower)) -1.0 else 0.0;
            if (std.math.sign(basis.row_scale[row]) != sign) return false;
            const magnitude = @abs(basis.row_scale[row]);
            basis.row_rhs[row] = if (sign > 0.0) upper * magnitude else if (sign < 0.0) -lower * magnitude else 0.0;
            basis.col_lower[logical] = 0.0;
            basis.col_upper[logical] = if (sign > 0.0 and std.math.isFinite(lower)) (upper - lower) * magnitude else std.math.inf(f64);
        }
        switch (basis.col_status[logical]) {
            .basic => {},
            .at_lower, .fixed => basis.primal[logical] = basis.col_lower[logical],
            .at_upper => {
                if (!std.math.isFinite(basis.col_upper[logical])) return false;
                basis.primal[logical] = basis.col_upper[logical];
            },
            .free, .superbasic => basis.primal[logical] = if (self.active_dual_initialization_strategy == .highs)
                @min(@max(0.0, basis.col_lower[logical]), basis.col_upper[logical])
            else
                0.0,
        }
    }
    for (basis.basic_index, 0..) |column, row| {
        basis.basic_lower[row] = basis.col_lower[column];
        basis.basic_upper[row] = basis.col_upper[column];
    }
    self.objective_scale = if (self.active_dual_initialization_strategy == .highs) 1.0 else objectiveScale(problem.col_cost);
    return true;
}

/// Initialize bounds, scaling and the logical basis using the selected policy.
pub fn initializeProblemStorage(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
    return if (self.active_dual_initialization_strategy == .highs)
        initializeProblemStorageHighs(self, problem)
    else
        initializeProblemStorageBaseline(self, problem);
}

/// Baseline initialization with row normalization and conditional column scaling.
fn initializeProblemStorageBaseline(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    basis.initializeSlackBasis();
    @memset(basis.primal, 0.0);
    @memset(basis.artificial_sign, 0.0);
    for (problem.col_lower, problem.col_upper, 0..) |lower, upper, col| {
        if (lower > upper) return .infeasible;
        if (lower == upper) {
            basis.primal[col] = lower;
            basis.col_status[col] = .fixed;
        } else if (std.math.isFinite(lower)) {
            basis.primal[col] = lower;
            basis.col_status[col] = .at_lower;
        } else if (std.math.isFinite(upper)) {
            basis.primal[col] = upper;
            basis.col_status[col] = .at_upper;
        } else {
            basis.primal[col] = 0.0;
            basis.col_status[col] = .free;
        }
        basis.col_lower[col] = lower;
        basis.col_upper[col] = upper;
    }
    @memset(basis.residual_work, 0.0);
    for (problem.matrix.row_indices, problem.matrix.values) |row, value| {
        if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(value)) continue;
        const row_index = row.toUsize();
        basis.residual_work[row_index] = @max(basis.residual_work[row_index], @abs(value));
    }
    for (problem.row_lower, problem.row_upper, 0..) |lower, upper, row| {
        if (lower > upper) return .infeasible;
        const maximum = basis.residual_work[row];
        const magnitude = if (maximum == 0.0)
            1.0
        else
            @exp2(std.math.clamp(@round(-@log2(maximum)), -20.0, 20.0));
        if (std.math.isFinite(upper)) {
            basis.row_scale[row] = magnitude;
            basis.row_rhs[row] = upper * magnitude;
            basis.col_lower[problem.num_cols + row] = 0.0;
            basis.col_upper[problem.num_cols + row] = if (std.math.isFinite(lower)) (upper - lower) * magnitude else std.math.inf(f64);
        } else if (std.math.isFinite(lower)) {
            basis.row_scale[row] = -magnitude;
            basis.row_rhs[row] = -lower * magnitude;
            basis.col_lower[problem.num_cols + row] = 0.0;
            basis.col_upper[problem.num_cols + row] = std.math.inf(f64);
        } else {
            basis.row_scale[row] = 0.0;
            basis.row_rhs[row] = 0.0;
            basis.col_lower[problem.num_cols + row] = 0.0;
            basis.col_upper[problem.num_cols + row] = std.math.inf(f64);
        }
        basis.basic_lower[row] = basis.col_lower[problem.num_cols + row];
        basis.basic_upper[row] = basis.col_upper[problem.num_cols + row];
    }
    var scaled_matrix_min = std.math.inf(f64);
    var scaled_matrix_max: f64 = 0.0;
    for (problem.matrix.row_indices, problem.matrix.values) |row, coefficient| {
        if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
        const magnitude = @abs(basis.row_scale[row.toUsize()] * coefficient);
        if (magnitude != 0.0) scaled_matrix_min = @min(scaled_matrix_min, magnitude);
        scaled_matrix_max = @max(scaled_matrix_max, magnitude);
    }
    const use_column_scaling = scaled_matrix_min < std.math.inf(f64) and scaled_matrix_max / scaled_matrix_min > 1e6;
    for (0..problem.num_cols) |column| {
        const begin = problem.matrix.col_starts[column];
        const end = problem.matrix.col_starts[column + 1];
        var maximum: f64 = 0.0;
        for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
            if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
            maximum = @max(maximum, @abs(basis.row_scale[row.toUsize()] * coefficient));
        }
        const scale = if (use_column_scaling) powerOfTwoScale(maximum, 20) else 1.0;
        basis.column_scale[column] = scale;
        basis.col_lower[column] = problem.col_lower[column] / scale;
        basis.col_upper[column] = problem.col_upper[column] / scale;
        basis.primal[column] /= scale;
    }
    self.objective_scale = objectiveScale(problem.col_cost);
    return .optimal;
}

/// HiGHS-compatible six-pass power-of-two equilibration and logical bounds.
fn initializeProblemStorageHighs(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    basis.initializeSlackBasis();
    @memset(basis.primal, 0.0);
    @memset(basis.artificial_sign, 0.0);
    for (problem.col_lower, problem.col_upper, 0..) |lower, upper, col| {
        if (lower > upper) return .infeasible;
        if (lower == upper) {
            basis.primal[col] = lower;
            basis.col_status[col] = .fixed;
        } else if (std.math.isFinite(lower)) {
            basis.primal[col] = lower;
            basis.col_status[col] = .at_lower;
        } else if (std.math.isFinite(upper)) {
            basis.primal[col] = upper;
            basis.col_status[col] = .at_upper;
        } else {
            basis.primal[col] = 0.0;
            basis.col_status[col] = .free;
        }
        basis.col_lower[col] = lower;
        basis.col_upper[col] = upper;
    }
    // Pinned HiGHS forced-equilibration scaling. Alternate column and row
    // equilibration six times, include small nonzero costs in column scaling,
    // clamp factors to 2^±20, then round to powers of two.
    @memset(basis.row_scale, 1.0);
    @memset(basis.column_scale[0..problem.num_cols], 1.0);
    var minimum_nonzero_cost = std.math.inf(f64);
    for (problem.col_cost) |cost| if (cost != 0.0) {
        minimum_nonzero_cost = @min(minimum_nonzero_cost, @abs(cost));
    };
    const include_cost = minimum_nonzero_cost < 0.1;
    const minimum_scale = @exp2(-20.0);
    const maximum_scale = @exp2(20.0);
    for (0..6) |_| {
        @memset(basis.basic_lower, 1e200);
        @memset(basis.basic_upper, 1e-200);
        for (0..problem.num_cols) |column| {
            var minimum: f64 = 1e200;
            var maximum: f64 = 1e-200;
            const cost = @abs(problem.col_cost[column]);
            if (include_cost and cost != 0.0) {
                minimum = @min(minimum, cost);
                maximum = @max(maximum, cost);
            }
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
                if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                const value = @abs(coefficient) * basis.row_scale[row.toUsize()];
                minimum = @min(minimum, value);
                maximum = @max(maximum, value);
            }
            const equilibration = 1.0 / @sqrt(minimum * maximum);
            basis.column_scale[column] = std.math.clamp(equilibration, minimum_scale, maximum_scale);
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
                if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                const row_index = row.toUsize();
                const value = @abs(coefficient) * basis.column_scale[column];
                basis.basic_lower[row_index] = @min(basis.basic_lower[row_index], value);
                basis.basic_upper[row_index] = @max(basis.basic_upper[row_index], value);
            }
        }
        for (basis.row_scale, basis.basic_lower, basis.basic_upper) |*scale, minimum, maximum| {
            scale.* = std.math.clamp(1.0 / @sqrt(minimum * maximum), minimum_scale, maximum_scale);
        }
    }
    for (basis.column_scale[0..problem.num_cols]) |*scale|
        scale.* = @exp2(@floor(@log2(scale.*) + 0.5));
    for (basis.row_scale) |*scale|
        scale.* = @exp2(@floor(@log2(scale.*) + 0.5));

    for (problem.row_lower, problem.row_upper, 0..) |lower, upper, row| {
        if (lower > upper) return .infeasible;
        const magnitude = basis.row_scale[row];
        basis.row_scale[row] = magnitude;
        basis.row_rhs[row] = 0.0;
        basis.col_lower[problem.num_cols + row] = if (std.math.isFinite(upper)) -upper * magnitude else -std.math.inf(f64);
        basis.col_upper[problem.num_cols + row] = if (std.math.isFinite(lower)) -lower * magnitude else std.math.inf(f64);
        basis.basic_lower[row] = basis.col_lower[problem.num_cols + row];
        basis.basic_upper[row] = basis.col_upper[problem.num_cols + row];
    }
    for (0..problem.num_cols) |column| {
        const scale = basis.column_scale[column];
        basis.column_scale[column] = scale;
        basis.col_lower[column] = problem.col_lower[column] / scale;
        basis.col_upper[column] = problem.col_upper[column] / scale;
        basis.primal[column] /= scale;
    }
    self.objective_scale = 1.0;
    return .optimal;
}

/// Return a bounded power-of-two scale that normalizes `maximum` near one.
fn powerOfTwoScale(maximum: f64, exponent_limit: comptime_int) f64 {
    if (maximum == 0.0 or !std.math.isFinite(maximum)) return 1.0;
    return @exp2(std.math.clamp(@round(-@log2(maximum)), -@as(f64, exponent_limit), @as(f64, exponent_limit)));
}

/// Compute a power-of-two objective scale from the largest absolute cost.
fn objectiveScale(cost: []const f64) f64 {
    var maximum: f64 = 0.0;
    for (cost) |value| maximum = @max(maximum, @abs(value));
    return powerOfTwoScale(maximum, 15);
}

/// Materialize a structural, logical or artificial internal column in scaled row space.
pub fn fillInternalColumn(self: *SimplexEngine, problem: problem_module.ProblemView, column: usize, output: []f64) SolveStatus {
    const basis = if (self.basis) |*value| value else return .numerical_failure;
    @memset(output, 0.0);
    if (column < problem.num_cols) {
        const column_scale = basis.column_scale[column];
        const begin = problem.matrix.col_starts[column];
        const end = problem.matrix.col_starts[column + 1];
        for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
            if (!matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
            const row_index = row.toUsize();
            output[row_index] = coefficient * basis.row_scale[row_index] * column_scale;
        }
    } else if (column < problem.num_cols + problem.num_rows) {
        output[column - problem.num_cols] = 1.0;
    } else if (column < problem.num_cols + 2 * problem.num_rows) {
        const row = column - problem.num_cols - problem.num_rows;
        output[row] = basis.artificial_sign[row];
    } else return .numerical_failure;
    return .optimal;
}
