//! Revised-simplex orchestration boundary; policies remain replaceable.
const std = @import("std");
const basis_module = @import("basis.zig");
const basis_snapshot_module = @import("basis_snapshot.zig");
const factorization_module = @import("factorization.zig");
const pricing_module = @import("pricing.zig");
const ratio_module = @import("ratio_test.zig");
const numerical_module = @import("numerical.zig");
const problem_module = @import("problem.zig");
const solution_module = @import("solution.zig");
const foundation = @import("foundation");

pub const Algorithm = enum { primal_revised, dual_revised };
pub const SolvePhase = enum { phase_one, dual_feasibility_repair, phase_two };
pub const CallbackAction = enum { continue_solve, stop };
/// Borrowed scalar progress snapshot. It owns no memory and is valid only for
/// the callback invocation.
pub const ProgressEventView = struct {
    phase: SolvePhase,
    algorithm: Algorithm,
    iterations: usize,
    work_used: u64,
    objective_value: f64,
    primal_infeasibility: f64,
    dual_infeasibility: f64,
};
pub const IterationCallback = *const fn (event: ProgressEventView, user_data: ?*anyopaque) CallbackAction;
pub const IterationLogCallback = *const fn (event: ProgressEventView, user_data: ?*anyopaque) void;
pub const LogLevel = enum { off, iterations };
pub const SolveStatus = solution_module.SolveStatus;
pub const BasisImportError = basis_snapshot_module.BasisViewError || error{
    InvalidNonbasicStatus,
    SingularBasis,
    NumericalFailure,
};
pub const SolveControl = struct {
    max_iterations: usize = 1_000_000,
    /// `null` disables the wall-clock limit; zero requests immediate stop.
    time_limit_ns: ?u64 = null,
    /// Optional caller-owned atomic flag. The flag must outlive `solveProblem`.
    interrupt_flag: ?*const std.atomic.Value(bool) = null,
    /// Optional borrowed warm-start basis.
    initial_basis: ?basis_snapshot_module.BasisView = null,
    /// Optional clock provider for deterministic tests or custom runtimes.
    /// The value is borrowed only for the duration of a solve call.
    clock_io: ?std.Io = null,
    /// Deterministic simplex work budget. One unit is charged for each
    /// attempted primal, dual, or Phase-I iteration.
    work_limit: ?u64 = null,
    iteration_callback: ?IterationCallback = null,
    callback_user_data: ?*anyopaque = null,
    callback_interval_work: u64 = 1,
    log_level: LogLevel = .off,
    log_callback: ?IterationLogCallback = null,
    log_user_data: ?*anyopaque = null,
    log_interval_work: u64 = 100,
};

pub const SimplexEngine = struct {
    allocator: std.mem.Allocator,
    algorithm: Algorithm = .dual_revised,
    basis: ?basis_module.BasisState = null,
    factorization: factorization_module.Factorization,
    pricing: pricing_module.Pricing = .{},
    ratio_test: ratio_module.RatioTest = .{},
    numerical: numerical_module.NumericalState = .{},
    iterations: usize = 0,
    objective_value: f64 = 0.0,
    phase1_needed: bool = false,
    solve_start_ns: ?i96 = null,
    solve_clock_io: ?std.Io = null,
    work_used: u64 = 0,
    dual_edge_weights_valid: bool = false,
    dual_row_index: ?u32 = null,
    dual_candidate_count: usize = 0,
    dual_candidate_cutoff: f64 = 0.0,
    dual_hyper_sparse_active: bool = false,

    pub fn init(a: std.mem.Allocator) SimplexEngine {
        return .{ .allocator = a, .factorization = factorization_module.Factorization.init(a) };
    }
    pub fn deinit(self: *SimplexEngine) void {
        if (self.basis) |*b| b.deinit();
        self.factorization.deinit();
    }
    pub fn solve(_: *SimplexEngine, _: usize, _: usize, _: SolveControl) SolveStatus {
        return .not_implemented;
    }

    /// Solve entry point consuming a borrowed LP `ProblemView`.
    ///
    /// The view must outlive this call; the engine never takes ownership of
    /// model arrays. Basis storage is owned by the engine instance.
    pub fn solveProblem(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
        self.startSolveClock(control);
        self.work_used = 0;
        if (self.controlledStop(control)) |status| return status;
        problem.validate() catch |err| return switch (err) {
            error.InvalidBounds => .infeasible,
            error.DimensionMismatch, error.InvalidMatrix => .numerical_failure,
        };
        if (self.basis) |*old| old.deinit();
        self.basis = basis_module.BasisState.init(self.allocator, problem.num_rows, problem.num_cols) catch return .numerical_failure;
        if (self.initializeProblemStorage(problem) != .optimal) return .infeasible;
        self.factorization.factorizeIdentity(problem.num_rows) catch return .numerical_failure;
        self.dual_edge_weights_valid = true;
        self.observeFactorizationStability();
        self.iterations = 0;

        if (control.initial_basis) |initial_basis| {
            if (self.importBasis(problem, initial_basis)) |_| {
                if (self.recomputeReducedCosts(problem) != .optimal) return .numerical_failure;
                const feasibility = self.classifyFeasibility(problem);
                if (feasibility.primal) {
                    self.algorithm = .primal_revised;
                    return self.solvePrimal(problem, control);
                }
                if (feasibility.dual) {
                    self.algorithm = .dual_revised;
                    return self.solveDual(problem, control);
                }
                const repair_status = self.repairWarmBasisWithDual(problem, control);
                if (repair_status == .optimal) {
                    if (self.recomputeReducedCosts(problem) != .optimal) return .numerical_failure;
                    const repaired = self.classifyFeasibility(problem);
                    if (repaired.primal) {
                        self.algorithm = .primal_revised;
                        return self.solvePrimal(problem, control);
                    }
                    if (repaired.dual) {
                        self.algorithm = .dual_revised;
                        return self.solveDual(problem, control);
                    }
                } else if (repair_status == .work_limit or repair_status == .time_limit or
                    repair_status == .iteration_limit or repair_status == .interrupted)
                {
                    return repair_status;
                }
            } else |_| {}
            // Invalid, singular, or neither-feasible warm starts fall back to
            // the deterministic logical crash basis.
            if (self.initializeProblemStorage(problem) != .optimal) return .infeasible;
            self.factorization.factorizeIdentity(problem.num_rows) catch return .numerical_failure;
            self.dual_edge_weights_valid = true;
            self.observeFactorizationStability();
        }

        self.algorithm = .primal_revised;
        if (self.basis) |*basis| {
            @memcpy(basis.basic_value, basis.row_rhs);
            for (0..problem.num_cols) |col| {
                const initial = basis.primal[col];
                const begin = problem.matrix.col_starts[col];
                const end = problem.matrix.col_starts[col + 1];
                for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, value| {
                    const row_index = row.toUsize();
                    basis.basic_value[row_index] -= basis.row_scale[row_index] * value * initial;
                }
            }
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
        }
        if (self.phase1_needed and self.refactorizeBasis(problem) != .optimal) return .numerical_failure;
        self.numerical.markRefactorized();
        self.iterations = 0;
        if (self.phase1_needed) {
            const phase1_status = self.solvePhaseOne(problem, control);
            if (phase1_status != .optimal) return phase1_status;
        }
        if (self.controlledStop(control)) |status| return status;
        return self.solvePrimal(problem, control);
    }

    /// Reoptimize with the existing basis and factorization. The caller must
    /// guarantee unchanged matrix structure and values. No model view is
    /// retained after this call.
    pub fn reoptimizeProblem(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
        self.startSolveClock(control);
        self.work_used = 0;
        if (self.controlledStop(control)) |status| return status;
        problem.validate() catch |err| return switch (err) {
            error.InvalidBounds => .infeasible,
            error.DimensionMismatch, error.InvalidMatrix => .numerical_failure,
        };
        if (!self.refreshProblemStorage(problem)) {
            var cold_control = control;
            cold_control.initial_basis = null;
            return self.solveProblem(problem, cold_control);
        }
        self.iterations = 0;
        if (self.recomputeBasicValuesUnchecked(problem) != .optimal or self.recomputeReducedCosts(problem) != .optimal)
            return .numerical_failure;
        const feasibility = self.classifyFeasibility(problem);
        if (feasibility.primal) {
            self.algorithm = .primal_revised;
            return self.solvePrimal(problem, control);
        }
        if (feasibility.dual) {
            self.algorithm = .dual_revised;
            return self.solveDual(problem, control);
        }
        const repair_status = self.repairWarmBasisWithDual(problem, control);
        if (repair_status == .optimal) {
            if (self.recomputeReducedCosts(problem) != .optimal) return .numerical_failure;
            if (self.classifyFeasibility(problem).primal) {
                self.algorithm = .primal_revised;
                return self.solvePrimal(problem, control);
            }
        } else if (repair_status == .work_limit or repair_status == .time_limit or
            repair_status == .iteration_limit or repair_status == .interrupted)
        {
            return repair_status;
        }
        var cold_control = control;
        cold_control.initial_basis = null;
        return self.solveProblem(problem, cold_control);
    }

    fn refreshProblemStorage(self: *SimplexEngine, problem: problem_module.ProblemView) bool {
        const basis = if (self.basis) |*value| value else return false;
        if (basis.num_rows != problem.num_rows or basis.num_structural_cols != problem.num_cols) return false;
        for (problem.col_lower, problem.col_upper, 0..) |lower, upper, column| {
            if (lower > upper) return false;
            basis.col_lower[column] = lower;
            basis.col_upper[column] = upper;
            switch (basis.col_status[column]) {
                .basic => {},
                .at_lower => {
                    if (!std.math.isFinite(lower)) return false;
                    basis.primal[column] = lower;
                },
                .at_upper => {
                    if (!std.math.isFinite(upper)) return false;
                    basis.primal[column] = upper;
                },
                .fixed => {
                    if (lower != upper) return false;
                    basis.primal[column] = lower;
                },
                .free, .superbasic => basis.primal[column] = @min(@max(0.0, lower), upper),
            }
        }
        for (problem.row_lower, problem.row_upper, 0..) |lower, upper, row| {
            const new_scale: f64 = if (std.math.isFinite(upper)) 1.0 else if (std.math.isFinite(lower)) -1.0 else 0.0;
            if (new_scale != basis.row_scale[row]) return false;
            basis.row_rhs[row] = if (new_scale > 0.0) upper else if (new_scale < 0.0) -lower else 0.0;
            const logical = problem.num_cols + row;
            basis.col_lower[logical] = 0.0;
            basis.col_upper[logical] = if (new_scale > 0.0 and std.math.isFinite(lower)) upper - lower else std.math.inf(f64);
            switch (basis.col_status[logical]) {
                .basic => {},
                .at_lower, .fixed => basis.primal[logical] = basis.col_lower[logical],
                .at_upper => {
                    if (!std.math.isFinite(basis.col_upper[logical])) return false;
                    basis.primal[logical] = basis.col_upper[logical];
                },
                .free, .superbasic => basis.primal[logical] = 0.0,
            }
        }
        for (basis.basic_index, 0..) |column, row| {
            basis.basic_lower[row] = basis.col_lower[column];
            basis.basic_upper[row] = basis.col_upper[column];
        }
        return true;
    }

    fn initializeProblemStorage(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
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
        for (problem.row_lower, problem.row_upper, 0..) |lower, upper, row| {
            if (lower > upper) return .infeasible;
            if (std.math.isFinite(upper)) {
                basis.row_scale[row] = 1.0;
                basis.row_rhs[row] = upper;
                basis.col_lower[problem.num_cols + row] = 0.0;
                basis.col_upper[problem.num_cols + row] = if (std.math.isFinite(lower)) upper - lower else std.math.inf(f64);
            } else if (std.math.isFinite(lower)) {
                basis.row_scale[row] = -1.0;
                basis.row_rhs[row] = -lower;
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
        return .optimal;
    }

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
        if (self.refactorizeBasis(problem) != .optimal) return error.SingularBasis;
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

    /// Compute the revised-simplex entering direction `B^-1 A_j` into the
    /// engine-owned SoA workspace without allocating or copying model data.
    pub fn computeDirection(self: *SimplexEngine, problem: problem_module.ProblemView, entering_col: usize) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        if (self.fillInternalColumn(problem, entering_col, basis.pivot_direction) != .optimal) return .numerical_failure;
        self.factorization.solve(basis.pivot_direction) catch return .numerical_failure;
        return .optimal;
    }

    /// Choose the leaving basic row for the currently materialized pivot
    /// direction. Returns `unbounded` when no positive pivot coefficient
    /// limits the entering variable.
    pub fn chooseLeaving(self: *SimplexEngine) struct { status: SolveStatus, row: ?u32 = null, step: f64 = 0.0, bound: basis_module.BasisStatus = .at_lower } {
        const basis = if (self.basis) |*value| value else return .{ .status = .numerical_failure };
        for (basis.basic_margin, basis.ratio_direction, basis.basic_value, basis.basic_lower, basis.basic_upper, basis.pivot_direction) |*margin, *ratio_direction, value, lower, upper, direction| {
            if (direction > self.numerical.zero_tolerance) {
                margin.* = value - lower;
                ratio_direction.* = direction;
            } else if (direction < -self.numerical.zero_tolerance) {
                margin.* = upper - value;
                ratio_direction.* = -direction;
            } else {
                margin.* = 0.0;
                ratio_direction.* = 0.0;
            }
        }
        const choice = self.ratio_test.chooseLeaving(basis.ratio_direction, basis.basic_margin);
        if (choice.row == null) return .{ .status = .unbounded };
        const row = choice.row.?;
        const bound: basis_module.BasisStatus = if (basis.pivot_direction[@intCast(row)] > 0) .at_lower else .at_upper;
        return .{ .status = .optimal, .row = row, .step = choice.step, .bound = bound };
    }

    /// Apply one primal pivot and rebuild the dense basis factorization in
    /// existing storage. This is allocation-free; update factorizations can
    /// replace the reinversion later without changing the state transition.
    pub fn performPivot(self: *SimplexEngine, problem: problem_module.ProblemView, entering_col: usize, entering_direction: f64, leaving_row: usize, leaving_bound: basis_module.BasisStatus, step: f64) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        if (entering_col >= basis.col_status.len or leaving_row >= problem.num_rows) return .numerical_failure;
        const leaving_col = basis.basic_index[leaving_row];
        const pivot = basis.pivot_direction[leaving_row];
        if (self.pricing.rule == .steepest_edge and self.dual_row_index == leaving_row) {
            if (self.updateDualSteepestEdgeWeights(leaving_row, pivot) != .optimal)
                self.dual_edge_weights_valid = false;
        } else {
            self.dual_edge_weights_valid = false;
            self.updateDevexWeights(entering_col, leaving_row, pivot);
        }
        self.dual_row_index = null;
        const update_succeeded = blk: {
            self.factorization.update(.{
                .leaving_row = @intCast(leaving_row),
                .entering_col = @intCast(entering_col),
                .direction = basis.pivot_direction,
                .column_scale = entering_direction,
            }) catch break :blk false;
            break :blk true;
        };
        self.numerical.observePivot(pivot);
        for (basis.basic_value, basis.pivot_direction) |*value, direction| value.* -= step * direction;
        const entering_value = basis.primal[entering_col] + entering_direction * step;
        basis.basic_value[leaving_row] = entering_value;
        basis.basic_lower[leaving_row] = basis.col_lower[entering_col];
        basis.basic_upper[leaving_row] = basis.col_upper[entering_col];
        basis.primal[entering_col] = entering_value;
        basis.primal[leaving_col] = if (leaving_bound == .at_upper) basis.col_upper[leaving_col] else basis.col_lower[leaving_col];
        basis.applyPivot(leaving_row, entering_col, leaving_bound) catch return .numerical_failure;
        for (basis.basic_index, basis.basic_value) |global_col, value| basis.primal[global_col] = value;
        if (!update_succeeded or self.factorization.needsRefactor(self.numerical.max_update_count) or self.numerical.needsRefactor())
            return self.refactorizeBasis(problem);
        return .optimal;
    }

    /// Lightweight Devex reference update. It deliberately uses the already
    /// hot FTRAN direction and does not allocate. Exact steepest-edge weights
    /// can replace this policy without changing basis storage or pivot code.
    fn updateDevexWeights(self: *SimplexEngine, entering_col: usize, leaving_row: usize, pivot: f64) void {
        const basis = if (self.basis) |*value| value else return;
        if (entering_col >= basis.col_edge_weight.len or leaving_row >= basis.row_edge_weight.len) return;
        var norm_squared: f64 = 1.0;
        for (basis.pivot_direction) |value| norm_squared += value * value;
        if (!std.math.isFinite(norm_squared)) norm_squared = 1.0;
        const entering_weight = @max(1.0, norm_squared);
        basis.col_edge_weight[entering_col] = entering_weight;
        basis.row_edge_weight[leaving_row] = if (@abs(pivot) > self.numerical.pivot_tolerance)
            @max(1.0, entering_weight / (pivot * pivot))
        else
            entering_weight;

        if (self.pricing.devex_reset_period != 0 and self.pricing.iterations % self.pricing.devex_reset_period == 0) {
            @memset(basis.col_edge_weight, 1.0);
            @memset(basis.row_edge_weight, 1.0);
        }
    }

    fn solvePrimal(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
        if (self.recomputeReducedCosts(problem) != .optimal) return .numerical_failure;
        while (self.iterations < control.max_iterations) : (self.iterations += 1) {
            if (self.beginIteration(problem, control, .phase_two)) |status| return status;
            const basis = if (self.basis) |*value| value else return .numerical_failure;
            const original_cols = problem.num_cols + problem.num_rows;
            const entering = self.pricing.choosePrimalEnteringWeighted(
                basis.reduced_cost[0..original_cols],
                basis.col_status[0..original_cols],
                basis.col_edge_weight[0..original_cols],
                self.numerical.dual_tolerance,
            ) orelse {
                return self.finishOptimal(problem);
            };
            if (self.computeDirection(problem, entering.column) != .optimal) return .numerical_failure;
            if (entering.direction < 0) {
                for (basis.pivot_direction) |*value| value.* = -value.*;
            }
            const leaving = self.chooseLeaving();
            const own_step = if (entering.direction > 0)
                basis.col_upper[entering.column] - basis.primal[entering.column]
            else
                basis.primal[entering.column] - basis.col_lower[entering.column];
            if (std.math.isFinite(own_step) and (leaving.status == .unbounded or own_step < leaving.step)) {
                for (basis.basic_value, basis.pivot_direction) |*value, direction| value.* -= own_step * direction;
                basis.primal[entering.column] = if (entering.direction > 0) basis.col_upper[entering.column] else basis.col_lower[entering.column];
                basis.col_status[entering.column] = if (entering.direction > 0) .at_upper else .at_lower;
            } else {
                if (leaving.status != .optimal) return leaving.status;
                if (self.performPivot(problem, entering.column, entering.direction, leaving.row.?, leaving.bound, leaving.step) != .optimal) return .numerical_failure;
            }
            if (self.recomputeReducedCosts(problem) != .optimal) return .numerical_failure;
        }
        return .iteration_limit;
    }

    fn solveDual(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
        if (self.recomputeReducedCosts(problem) != .optimal) return .numerical_failure;
        while (self.iterations < control.max_iterations) : (self.iterations += 1) {
            if (self.beginIteration(problem, control, .phase_two)) |status| return status;
            const basis = if (self.basis) |*value| value else return .numerical_failure;
            if (self.pricing.rule == .steepest_edge and self.ensureExactDualEdgeWeights() != .optimal)
                return .numerical_failure;
            const leaving = self.chooseDualLeavingRow() orelse return self.finishOptimal(problem);

            if (self.computeDualTableauRow(problem, leaving.row) != .optimal) return .numerical_failure;
            const original_cols = problem.num_cols + problem.num_rows;
            const entering = self.ratio_test.chooseDualEntering(
                basis.tableau[0..original_cols],
                basis.reduced_cost[0..original_cols],
                basis.col_status[0..original_cols],
                basis.col_lower[0..original_cols],
                basis.col_upper[0..original_cols],
                basis.primal[0..original_cols],
                leaving.bound,
                leaving.violation,
                basis.dual_ratio[0..original_cols],
                basis.dual_direction[0..original_cols],
                basis.flip_columns[0..original_cols],
            );

            if (self.applyBoundFlips(problem, entering.flip_count) != .optimal) return .numerical_failure;
            const entering_col = entering.column orelse return .infeasible;
            const entering_index: usize = @intCast(entering_col);
            if (self.computeDirection(problem, entering_index) != .optimal) return .numerical_failure;
            if (entering.direction < 0.0) {
                for (basis.pivot_direction) |*value| value.* = -value.*;
            }
            const leaving_row: usize = @intCast(leaving.row);
            const target = if (leaving.bound == .at_lower) basis.basic_lower[leaving_row] else basis.basic_upper[leaving_row];
            const pivot = basis.pivot_direction[leaving_row];
            if (@abs(pivot) <= self.numerical.pivot_tolerance) return .numerical_failure;
            const step = (basis.basic_value[leaving_row] - target) / pivot;
            if (!std.math.isFinite(step) or step < -self.numerical.primal_tolerance) return .numerical_failure;
            if (self.performPivot(
                problem,
                entering_index,
                entering.direction,
                leaving_row,
                leaving.bound,
                @max(step, 0.0),
            ) != .optimal) return .numerical_failure;
            if (self.recomputeReducedCosts(problem) != .optimal) return .numerical_failure;
            if (!self.classifyFeasibility(problem).dual) return .numerical_failure;
        }
        return .iteration_limit;
    }

    fn chooseDualLeavingRow(self: *SimplexEngine) ?pricing_module.DualLeavingChoice {
        const basis = if (self.basis) |*value| value else return null;
        if (self.pricing.rule != .hyper_sparse or !self.dual_hyper_sparse_active)
            return self.pricing.chooseDualLeavingWeighted(
                basis.basic_value,
                basis.basic_lower,
                basis.basic_upper,
                basis.row_edge_weight,
                self.numerical.primal_tolerance,
            );

        var best = self.bestDualCandidate();
        if (best == null or self.dualCandidateScore(best.?.row) + self.numerical.primal_tolerance < self.dual_candidate_cutoff) {
            self.rebuildDualCandidateList();
            best = self.bestDualCandidate();
        }
        return best;
    }

    fn dualCandidate(self: *const SimplexEngine, row: usize) ?pricing_module.DualLeavingChoice {
        const basis = if (self.basis) |*value| value else return null;
        const value = basis.basic_value[row];
        if (value < basis.basic_lower[row] - self.numerical.primal_tolerance)
            return .{ .row = @intCast(row), .bound = .at_lower, .violation = basis.basic_lower[row] - value };
        if (value > basis.basic_upper[row] + self.numerical.primal_tolerance)
            return .{ .row = @intCast(row), .bound = .at_upper, .violation = value - basis.basic_upper[row] };
        return null;
    }

    fn dualCandidateScore(self: *const SimplexEngine, row_u32: u32) f64 {
        const row: usize = @intCast(row_u32);
        const basis = if (self.basis) |*value| value else return 0.0;
        const candidate = self.dualCandidate(row) orelse return 0.0;
        return candidate.violation / @sqrt(@max(basis.row_edge_weight[row], 1.0));
    }

    fn bestDualCandidate(self: *SimplexEngine) ?pricing_module.DualLeavingChoice {
        const basis = if (self.basis) |*value| value else return null;
        var best: ?pricing_module.DualLeavingChoice = null;
        var best_score = self.numerical.primal_tolerance;
        for (basis.dual_candidate_rows[0..self.dual_candidate_count], basis.dual_candidate_score[0..self.dual_candidate_count]) |row, *stored_score| {
            const score = self.dualCandidateScore(row);
            stored_score.* = score;
            if (score > best_score) {
                best_score = score;
                best = self.dualCandidate(@intCast(row));
            }
        }
        return best;
    }

    fn rebuildDualCandidateList(self: *SimplexEngine) void {
        const basis = if (self.basis) |*value| value else return;
        const capacity = @min(basis.num_rows, 32);
        self.dual_candidate_count = 0;
        self.dual_candidate_cutoff = 0.0;
        if (capacity == 0) return;
        for (0..basis.num_rows) |row| {
            const score = self.dualCandidateScore(@intCast(row));
            if (score <= self.numerical.primal_tolerance) continue;
            if (self.dual_candidate_count < capacity) {
                const slot = self.dual_candidate_count;
                basis.dual_candidate_rows[slot] = @intCast(row);
                basis.dual_candidate_score[slot] = score;
                self.dual_candidate_count += 1;
                continue;
            }
            var weakest: usize = 0;
            for (basis.dual_candidate_score[1..capacity], 1..) |candidate_score, slot| {
                if (candidate_score < basis.dual_candidate_score[weakest]) weakest = slot;
            }
            if (score > basis.dual_candidate_score[weakest]) {
                basis.dual_candidate_rows[weakest] = @intCast(row);
                basis.dual_candidate_score[weakest] = score;
            }
        }
        if (self.dual_candidate_count > 0) {
            self.dual_candidate_cutoff = basis.dual_candidate_score[0];
            for (basis.dual_candidate_score[1..self.dual_candidate_count]) |score|
                self.dual_candidate_cutoff = @min(self.dual_candidate_cutoff, score);
        }
    }

    /// Initialize exact dual steepest-edge weights after a crash, imported
    /// basis, reinversion, or detected recurrence drift.
    fn ensureExactDualEdgeWeights(self: *SimplexEngine) SolveStatus {
        if (self.dual_edge_weights_valid) return .optimal;
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        for (basis.row_edge_weight, 0..) |*weight, row| {
            @memset(basis.dual_row, 0.0);
            basis.dual_row[row] = 1.0;
            self.factorization.solveTranspose(basis.dual_row) catch return .numerical_failure;
            var norm_squared: f64 = 0.0;
            for (basis.dual_row) |entry| norm_squared += entry * entry;
            if (!std.math.isFinite(norm_squared)) return .numerical_failure;
            weight.* = @max(norm_squared, self.numerical.zero_tolerance);
        }
        self.dual_edge_weights_valid = true;
        self.dual_row_index = null;
        return .optimal;
    }

    /// Forrest--Goldfarb DSE recurrence. `dual_row` is the freshly computed
    /// BTRAN result B^-T e_p. One additional FTRAN forms
    /// tau = B^-1 B^-T e_p before the factor update changes B.
    fn updateDualSteepestEdgeWeights(self: *SimplexEngine, leaving_row: usize, pivot: f64) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        if (!self.dual_edge_weights_valid or leaving_row >= basis.row_edge_weight.len or
            @abs(pivot) <= self.numerical.pivot_tolerance)
            return .numerical_failure;
        var exact_pivot_weight: f64 = 0.0;
        for (basis.dual_row) |entry| exact_pivot_weight += entry * entry;
        if (!std.math.isFinite(exact_pivot_weight) or exact_pivot_weight <= self.numerical.zero_tolerance)
            return .numerical_failure;
        basis.row_edge_weight[leaving_row] = exact_pivot_weight;
        @memcpy(basis.residual_work, basis.dual_row);
        self.factorization.solve(basis.residual_work) catch return .numerical_failure;

        for (basis.row_edge_weight, basis.pivot_direction, basis.residual_work, 0..) |*weight, alpha, tau, row| {
            if (row == leaving_row) continue;
            const ratio = alpha / pivot;
            const updated = weight.* - 2.0 * ratio * tau + ratio * ratio * exact_pivot_weight;
            if (!std.math.isFinite(updated)) return .numerical_failure;
            weight.* = @max(updated, 1e-4);
        }
        basis.row_edge_weight[leaving_row] = @max(exact_pivot_weight / (pivot * pivot), 1e-4);
        return .optimal;
    }

    /// Repair a warm basis that is neither primal nor dual feasible by using
    /// the zero auxiliary objective. Every reduced cost is then exactly zero,
    /// so dual pivots can restore primal feasibility without discarding the
    /// imported basis or factorization. The original objective is restored by
    /// the caller before entering Phase II.
    fn repairWarmBasisWithDual(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        self.algorithm = .dual_revised;
        @memset(basis.reduced_cost, 0.0);
        while (self.iterations < control.max_iterations) : (self.iterations += 1) {
            if (self.beginIteration(problem, control, .dual_feasibility_repair)) |status| return status;
            const leaving = self.pricing.chooseDualLeavingWeighted(
                basis.basic_value,
                basis.basic_lower,
                basis.basic_upper,
                basis.row_edge_weight,
                self.numerical.primal_tolerance,
            ) orelse return .optimal;
            if (self.computeDualTableauRow(problem, leaving.row) != .optimal) return .numerical_failure;
            const original_cols = problem.num_cols + problem.num_rows;
            const entering = self.ratio_test.chooseDualEntering(
                basis.tableau[0..original_cols],
                basis.reduced_cost[0..original_cols],
                basis.col_status[0..original_cols],
                basis.col_lower[0..original_cols],
                basis.col_upper[0..original_cols],
                basis.primal[0..original_cols],
                leaving.bound,
                leaving.violation,
                basis.dual_ratio[0..original_cols],
                basis.dual_direction[0..original_cols],
                basis.flip_columns[0..original_cols],
            );
            if (self.applyBoundFlips(problem, entering.flip_count) != .optimal) return .numerical_failure;
            const entering_col = entering.column orelse return .infeasible;
            const entering_index: usize = @intCast(entering_col);
            if (self.computeDirection(problem, entering_index) != .optimal) return .numerical_failure;
            if (entering.direction < 0.0) {
                for (basis.pivot_direction) |*value| value.* = -value.*;
            }
            const leaving_row: usize = @intCast(leaving.row);
            const target = if (leaving.bound == .at_lower) basis.basic_lower[leaving_row] else basis.basic_upper[leaving_row];
            const pivot = basis.pivot_direction[leaving_row];
            if (@abs(pivot) <= self.numerical.pivot_tolerance) return .numerical_failure;
            const step = (basis.basic_value[leaving_row] - target) / pivot;
            if (!std.math.isFinite(step) or step < -self.numerical.primal_tolerance) return .numerical_failure;
            if (self.performPivot(
                problem,
                entering_index,
                entering.direction,
                leaving_row,
                leaving.bound,
                @max(step, 0.0),
            ) != .optimal) return .numerical_failure;
            @memset(basis.reduced_cost, 0.0);
        }
        return .iteration_limit;
    }

    fn computeDualTableauRow(self: *SimplexEngine, problem: problem_module.ProblemView, leaving_row: u32) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        const row: usize = @intCast(leaving_row);
        if (row >= problem.num_rows) return .numerical_failure;
        @memset(basis.dual_row, 0.0);
        basis.dual_row[row] = 1.0;
        self.factorization.solveTranspose(basis.dual_row) catch return .numerical_failure;
        if (self.pricing.rule == .steepest_edge and self.dual_edge_weights_valid) {
            var exact_weight: f64 = 0.0;
            for (basis.dual_row) |entry| exact_weight += entry * entry;
            if (!std.math.isFinite(exact_weight) or exact_weight <= self.numerical.zero_tolerance)
                return .numerical_failure;
            const updated_weight = basis.row_edge_weight[row];
            const relative_error = @abs(updated_weight - exact_weight) / @max(exact_weight, self.numerical.zero_tolerance);
            if (relative_error > self.numerical.dual_edge_weight_error_tolerance) {
                basis.row_edge_weight[row] = exact_weight;
                self.numerical.dual_edge_weight_corrections += 1;
                // Huangfu's safeguard rejects a seriously underestimated
                // pivotal weight; force a full exact reset after this pivot.
                if (updated_weight < 0.5 * exact_weight) self.dual_edge_weights_valid = false;
            }
        }
        self.dual_row_index = leaving_row;
        @memset(basis.tableau, 0.0);
        for (0..problem.num_cols) |column| {
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            var alpha: f64 = 0.0;
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |matrix_row, coefficient| {
                const row_index = matrix_row.toUsize();
                alpha += basis.dual_row[row_index] * basis.row_scale[row_index] * coefficient;
            }
            basis.tableau[column] = alpha;
        }
        for (0..problem.num_rows) |logical_row| {
            basis.tableau[problem.num_cols + logical_row] = basis.dual_row[logical_row];
        }
        var nonzero_count: usize = 0;
        for (basis.tableau[0 .. problem.num_cols + problem.num_rows]) |value| {
            if (@abs(value) > self.numerical.zero_tolerance) nonzero_count += 1;
        }
        const tableau_count = problem.num_cols + problem.num_rows;
        self.dual_hyper_sparse_active = tableau_count > 0 and nonzero_count * 10 < tableau_count;
        if (!self.dual_hyper_sparse_active) self.dual_candidate_count = 0;
        return .optimal;
    }

    fn applyBoundFlips(self: *SimplexEngine, problem: problem_module.ProblemView, flip_count: usize) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        for (basis.flip_columns[0..flip_count]) |column_u32| {
            const column: usize = @intCast(column_u32);
            if (self.computeDirection(problem, column) != .optimal) return .numerical_failure;
            const delta = switch (basis.col_status[column]) {
                .at_lower => basis.col_upper[column] - basis.col_lower[column],
                .at_upper => -(basis.col_upper[column] - basis.col_lower[column]),
                else => return .numerical_failure,
            };
            if (!std.math.isFinite(delta)) return .numerical_failure;
            for (basis.basic_value, basis.pivot_direction) |*value, direction| value.* -= direction * delta;
            basis.primal[column] += delta;
            basis.col_status[column] = if (delta > 0.0) .at_upper else .at_lower;
            for (basis.basic_index, basis.basic_value) |basic_column, value| basis.primal[basic_column] = value;
        }
        return .optimal;
    }

    fn solvePhaseOne(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
        if (self.recomputePhaseOneReducedCosts(problem) != .optimal) return .numerical_failure;
        while (self.iterations < control.max_iterations) : (self.iterations += 1) {
            if (self.beginIteration(problem, control, .phase_one)) |status| return status;
            const basis = if (self.basis) |*value| value else return .numerical_failure;
            const entering = self.pricing.choosePrimalEnteringWeighted(
                basis.reduced_cost,
                basis.col_status,
                basis.col_edge_weight,
                self.numerical.dual_tolerance,
            ) orelse break;
            if (self.computeDirection(problem, entering.column) != .optimal) return .numerical_failure;
            if (entering.direction < 0) {
                for (basis.pivot_direction) |*value| value.* = -value.*;
            }
            const leaving = self.chooseLeaving();
            const own_step = if (entering.direction > 0)
                basis.col_upper[entering.column] - basis.primal[entering.column]
            else
                basis.primal[entering.column] - basis.col_lower[entering.column];
            if (std.math.isFinite(own_step) and (leaving.status == .unbounded or own_step < leaving.step)) {
                for (basis.basic_value, basis.pivot_direction) |*value, direction| value.* -= own_step * direction;
                basis.primal[entering.column] += entering.direction * own_step;
                basis.col_status[entering.column] = if (entering.direction > 0) .at_upper else .at_lower;
            } else {
                if (leaving.status != .optimal) return leaving.status;
                if (self.performPivot(problem, entering.column, entering.direction, leaving.row.?, leaving.bound, leaving.step) != .optimal)
                    return .numerical_failure;
            }
            if (self.recomputePhaseOneReducedCosts(problem) != .optimal) return .numerical_failure;
        }
        if (self.iterations >= control.max_iterations) return .iteration_limit;
        if (self.recomputeBasicValues(problem) != .optimal) return .numerical_failure;
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        const artificial_begin = problem.num_cols + problem.num_rows;
        var infeasibility: f64 = 0.0;
        for (basis.primal[artificial_begin..]) |value| infeasibility += @max(value, 0.0);
        if (infeasibility > self.numerical.primal_tolerance) return .infeasible;

        if (self.cleanupArtificialBasis(problem) != .optimal) return .numerical_failure;
        @memset(basis.col_upper[artificial_begin..], 0.0);
        for (basis.col_status[artificial_begin..]) |*status| {
            if (status.* != .basic) status.* = .fixed;
        }
        for (basis.basic_index, 0..) |column, row| {
            if (column >= artificial_begin) basis.basic_upper[row] = 0.0;
        }
        self.phase1_needed = false;
        return .optimal;
    }

    fn controlledStop(self: *SimplexEngine, control: SolveControl) ?SolveStatus {
        if (control.interrupt_flag) |flag| {
            if (flag.load(.acquire)) return .interrupted;
        }
        if (control.work_limit) |limit| {
            if (self.work_used >= limit) return .work_limit;
        }
        if (control.time_limit_ns) |limit| {
            if (limit == 0) return .time_limit;
            if (self.solve_start_ns) |start| {
                const io = self.solve_clock_io orelse return null;
                const now = std.Io.Clock.awake.now(io).nanoseconds;
                if (now >= start and @as(u128, @intCast(now - start)) >= limit) return .time_limit;
            }
        }
        return null;
    }

    fn beginIteration(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl, phase: SolvePhase) ?SolveStatus {
        if (self.controlledStop(control)) |status| return status;
        const callback_due = control.iteration_callback != null and
            self.work_used % @max(control.callback_interval_work, 1) == 0;
        const log_due = control.log_level == .iterations and control.log_callback != null and
            self.work_used % @max(control.log_interval_work, 1) == 0;
        if (callback_due or log_due) {
            const event = self.progressEvent(problem, phase);
            if (callback_due) {
                if (control.iteration_callback.?(event, control.callback_user_data) == .stop) return .interrupted;
            }
            if (log_due) control.log_callback.?(event, control.log_user_data);
        }
        self.work_used = std.math.add(u64, self.work_used, 1) catch std.math.maxInt(u64);
        return null;
    }

    fn progressEvent(self: *const SimplexEngine, problem: problem_module.ProblemView, phase: SolvePhase) ProgressEventView {
        var primal_infeasibility: f64 = 0.0;
        var dual_infeasibility: f64 = 0.0;
        var current_objective = problem.objective_offset;
        if (self.basis) |*basis| {
            for (problem.col_cost, basis.primal[0..problem.num_cols]) |cost, value| current_objective += cost * value;
            for (basis.basic_value, basis.basic_lower, basis.basic_upper) |value, lower, upper| {
                primal_infeasibility = @max(primal_infeasibility, @max(lower - value, value - upper));
            }
            for (basis.reduced_cost, basis.col_status) |reduced, status| {
                const violation: f64 = switch (status) {
                    .at_lower => -reduced,
                    .at_upper => reduced,
                    .free, .superbasic => @abs(reduced),
                    .basic, .fixed => 0.0,
                };
                dual_infeasibility = @max(dual_infeasibility, violation);
            }
        }
        return .{
            .phase = phase,
            .algorithm = self.algorithm,
            .iterations = self.iterations,
            .work_used = self.work_used,
            .objective_value = current_objective,
            .primal_infeasibility = @max(primal_infeasibility, 0.0),
            .dual_infeasibility = @max(dual_infeasibility, 0.0),
        };
    }

    fn startSolveClock(self: *SimplexEngine, control: SolveControl) void {
        if (control.time_limit_ns == null) {
            self.solve_start_ns = null;
            self.solve_clock_io = null;
            return;
        }
        const io = control.clock_io orelse std.Io.Threaded.global_single_threaded.io();
        self.solve_clock_io = io;
        self.solve_start_ns = std.Io.Clock.awake.now(io).nanoseconds;
    }

    /// Pivot zero-valued artificial basics out whenever a stable original or
    /// logical column is available. An artificial that remains basic denotes
    /// a rank-redundant row and is fixed at zero until presolve can remove that
    /// row explicitly.
    fn cleanupArtificialBasis(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        const artificial_begin = problem.num_cols + problem.num_rows;

        for (0..problem.num_rows) |leaving_row| {
            if (basis.basic_index[leaving_row] < artificial_begin) continue;

            var entering_col: ?usize = null;
            var best_pivot = self.numerical.pivot_tolerance;
            for (0..artificial_begin) |column| {
                if (basis.col_status[column] == .basic) continue;
                if (self.computeDirection(problem, column) != .optimal) return .numerical_failure;
                const candidate = @abs(basis.pivot_direction[leaving_row]);
                if (candidate > best_pivot) {
                    best_pivot = candidate;
                    entering_col = column;
                }
            }

            const column = entering_col orelse continue;
            if (self.computeDirection(problem, column) != .optimal) return .numerical_failure;
            if (self.performPivot(problem, column, 1.0, leaving_row, .fixed, 0.0) != .optimal)
                return .numerical_failure;
        }
        return self.recomputeBasicValues(problem);
    }

    fn recomputePhaseOneReducedCosts(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        const artificial_begin = problem.num_cols + problem.num_rows;
        for (basis.basic_index, 0..) |global_col, row| basis.dual[row] = if (global_col >= artificial_begin) 1.0 else 0.0;
        self.factorization.solveTranspose(basis.dual) catch return .numerical_failure;
        for (basis.reduced_cost, 0..) |*reduced, column| {
            reduced.* = if (column >= artificial_begin) 1.0 else 0.0;
            if (self.fillInternalColumn(problem, column, basis.pivot_direction) != .optimal) return .numerical_failure;
            for (basis.pivot_direction, basis.dual) |value, dual| reduced.* -= value * dual;
        }
        return .optimal;
    }

    fn recomputeReducedCosts(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        const maximize = problem.objective_sense == .maximize;
        for (basis.basic_index, 0..) |global_col, row| {
            basis.dual[row] = if (global_col < problem.num_cols)
                (if (maximize) -problem.col_cost[global_col] else problem.col_cost[global_col])
            else
                0.0;
        }
        self.factorization.solveTranspose(basis.dual) catch return .numerical_failure;
        for (0..problem.num_cols) |col| {
            var reduced = if (maximize) -problem.col_cost[col] else problem.col_cost[col];
            const begin = problem.matrix.col_starts[col];
            const end = problem.matrix.col_starts[col + 1];
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, value| {
                const row_index = row.toUsize();
                reduced -= basis.row_scale[row_index] * value * basis.dual[row_index];
            }
            basis.reduced_cost[col] = reduced;
        }
        for (0..problem.num_rows) |row| basis.reduced_cost[problem.num_cols + row] = -basis.dual[row];
        @memset(basis.reduced_cost[problem.num_cols + problem.num_rows ..], 0.0);
        return .optimal;
    }

    fn refactorizeBasis(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        const n = problem.num_rows;
        const buffer = self.factorization.mutableBasisBuffer(n) catch return .numerical_failure;
        @memset(buffer, 0.0);
        for (basis.basic_index, 0..) |global_col, basis_col| {
            if (global_col < basis.num_structural_cols) {
                const col: usize = @intCast(global_col);
                const begin = problem.matrix.col_starts[col];
                const end = problem.matrix.col_starts[col + 1];
                for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, value| {
                    const row_index = row.toUsize();
                    buffer[row_index * n + basis_col] = basis.row_scale[row_index] * value;
                }
            } else {
                const internal = @as(usize, @intCast(global_col)) - basis.num_structural_cols;
                if (internal < n) {
                    buffer[internal * n + basis_col] = 1.0;
                } else {
                    const artificial_row = internal - n;
                    if (artificial_row >= n) return .numerical_failure;
                    buffer[artificial_row * n + basis_col] = basis.artificial_sign[artificial_row];
                }
            }
        }
        self.factorization.refactorizeInPlace() catch return .numerical_failure;
        self.numerical.markRefactorized();
        self.observeFactorizationStability();
        self.dual_edge_weights_valid = false;
        self.dual_row_index = null;
        if (self.pricing.rule == .steepest_edge and self.ensureExactDualEdgeWeights() != .optimal)
            return .numerical_failure;
        return .optimal;
    }

    fn observeFactorizationStability(self: *SimplexEngine) void {
        self.numerical.pivot_condition_estimate = self.factorization.dense_lu.pivotConditionEstimate();
        if (!std.math.isFinite(self.numerical.pivot_condition_estimate) or self.numerical.pivot_condition_estimate > 1e12)
            self.numerical.numerical_warning = true;
    }

    fn fillInternalColumn(self: *SimplexEngine, problem: problem_module.ProblemView, column: usize, output: []f64) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        @memset(output, 0.0);
        if (column < problem.num_cols) {
            problem.fillColumn(column, output) catch return .numerical_failure;
            for (output, basis.row_scale) |*value, scale| value.* *= scale;
        } else if (column < problem.num_cols + problem.num_rows) {
            output[column - problem.num_cols] = 1.0;
        } else if (column < problem.num_cols + 2 * problem.num_rows) {
            const row = column - problem.num_cols - problem.num_rows;
            output[row] = basis.artificial_sign[row];
        } else return .numerical_failure;
        return .optimal;
    }

    /// Recompute the primal basic solution from the immutable problem and
    /// current nonbasic values. Uses engine-owned workspace and performs no
    /// allocation; this limits drift from long Eta update chains.
    fn recomputeBasicValues(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
        return self.recomputeBasicValuesImpl(problem, true);
    }

    fn recomputeBasicValuesUnchecked(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
        return self.recomputeBasicValuesImpl(problem, false);
    }

    fn recomputeBasicValuesImpl(self: *SimplexEngine, problem: problem_module.ProblemView, enforce_bounds: bool) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        @memcpy(basis.residual_work, basis.row_rhs);
        for (basis.col_status, basis.primal, 0..) |status, value, column| {
            if (status == .basic or value == 0.0) continue;
            if (column < problem.num_cols) {
                const begin = problem.matrix.col_starts[column];
                const end = problem.matrix.col_starts[column + 1];
                for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
                    const row_index = row.toUsize();
                    basis.residual_work[row_index] -= basis.row_scale[row_index] * coefficient * value;
                }
            } else if (column < problem.num_cols + problem.num_rows) {
                basis.residual_work[column - problem.num_cols] -= value;
            } else {
                const row = column - problem.num_cols - problem.num_rows;
                basis.residual_work[row] -= basis.artificial_sign[row] * value;
            }
        }
        @memcpy(basis.rhs_work, basis.residual_work);
        self.factorization.solve(basis.rhs_work) catch return .numerical_failure;

        var refinement_step: usize = 0;
        while (true) {
            @memcpy(basis.pivot_direction, basis.residual_work);
            self.subtractBasisProduct(problem, basis.rhs_work, basis.pivot_direction) catch return .numerical_failure;
            var residual_max: f64 = 0.0;
            var rhs_max: f64 = 0.0;
            for (basis.pivot_direction, basis.residual_work) |residual, rhs_value| {
                residual_max = @max(residual_max, @abs(residual));
                rhs_max = @max(rhs_max, @abs(rhs_value));
            }
            self.numerical.observeResidual(residual_max, rhs_max);
            if (!std.math.isFinite(residual_max)) return .numerical_failure;
            if (self.numerical.last_relative_residual <= self.numerical.residual_tolerance or
                refinement_step >= self.numerical.max_refinement_steps)
                break;
            self.factorization.solve(basis.pivot_direction) catch return .numerical_failure;
            for (basis.rhs_work, basis.pivot_direction) |*value, correction| value.* += correction;
            refinement_step += 1;
            self.numerical.refinement_count += 1;
        }
        for (basis.basic_index, basis.basic_value, basis.rhs_work, 0..) |column, *value, recomputed, row| {
            value.* = recomputed;
            basis.primal[column] = recomputed;
            if (enforce_bounds and (recomputed < basis.basic_lower[row] - self.numerical.primal_tolerance or
                recomputed > basis.basic_upper[row] + self.numerical.primal_tolerance))
                return .numerical_failure;
        }
        return .optimal;
    }

    /// Compute `residual -= B * x` directly from the borrowed CSC model and
    /// basis membership. No dense basis copy or temporary column is required.
    fn subtractBasisProduct(self: *SimplexEngine, problem: problem_module.ProblemView, x: []const f64, residual: []f64) !void {
        const basis = if (self.basis) |*value| value else return error.NumericalFailure;
        if (x.len != problem.num_rows or residual.len != problem.num_rows) return error.NumericalFailure;
        for (basis.basic_index, x) |global_col, value| {
            if (value == 0.0) continue;
            const column: usize = @intCast(global_col);
            if (column < problem.num_cols) {
                const begin = problem.matrix.col_starts[column];
                const end = problem.matrix.col_starts[column + 1];
                for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
                    const row_index = row.toUsize();
                    if (row_index >= residual.len) return error.NumericalFailure;
                    residual[row_index] -= basis.row_scale[row_index] * coefficient * value;
                }
            } else if (column < problem.num_cols + problem.num_rows) {
                residual[column - problem.num_cols] -= value;
            } else if (column < problem.num_cols + 2 * problem.num_rows) {
                const row = column - problem.num_cols - problem.num_rows;
                residual[row] -= basis.artificial_sign[row] * value;
            } else return error.NumericalFailure;
        }
    }

    fn finishOptimal(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
        if (self.recomputeBasicValues(problem) != .optimal) return .numerical_failure;
        if (self.recomputeReducedCosts(problem) != .optimal) return .numerical_failure;
        if (self.validateOptimalSolution(problem) != .optimal) return .numerical_failure;
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        self.objective_value = problem.objective_offset;
        for (problem.col_cost, basis.primal[0..problem.num_cols]) |cost, value| self.objective_value += cost * value;
        if (!std.math.isFinite(self.objective_value)) return .numerical_failure;
        return .optimal;
    }

    /// Allocation-free KKT feasibility check used before publishing an
    /// optimal result. This catches accumulated update drift and incorrect
    /// bound/status transitions at the solver boundary.
    fn validateOptimalSolution(self: *SimplexEngine, problem: problem_module.ProblemView) SolveStatus {
        const basis = if (self.basis) |*value| value else return .numerical_failure;
        const primal_tolerance = self.numerical.primal_tolerance;
        const dual_tolerance = self.numerical.dual_tolerance;

        for (basis.primal[0..problem.num_cols], problem.col_lower, problem.col_upper) |value, lower, upper| {
            if (!std.math.isFinite(value) or value < lower - primal_tolerance or value > upper + primal_tolerance)
                return .numerical_failure;
        }

        @memset(basis.rhs_work, 0.0);
        for (0..problem.num_cols) |column| {
            const value = basis.primal[column];
            if (value == 0.0) continue;
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
                basis.rhs_work[row.toUsize()] += coefficient * value;
            }
        }
        for (basis.rhs_work, problem.row_lower, problem.row_upper) |activity, lower, upper| {
            if (!std.math.isFinite(activity) or activity < lower - primal_tolerance or activity > upper + primal_tolerance)
                return .numerical_failure;
        }

        for (basis.reduced_cost[0 .. problem.num_cols + problem.num_rows], basis.col_status[0 .. problem.num_cols + problem.num_rows]) |reduced, status| {
            if (!std.math.isFinite(reduced)) return .numerical_failure;
            const infeasible = switch (status) {
                .at_lower => reduced < -dual_tolerance,
                .at_upper => reduced > dual_tolerance,
                .free, .superbasic, .basic => @abs(reduced) > dual_tolerance,
                .fixed => false,
            };
            if (infeasible) return .numerical_failure;
        }

        const artificial_begin = problem.num_cols + problem.num_rows;
        for (basis.primal[artificial_begin..]) |value| {
            if (!std.math.isFinite(value) or @abs(value) > primal_tolerance) return .numerical_failure;
        }
        return .optimal;
    }

    /// Borrow the current engine-owned solution arrays without copying.
    pub fn solutionView(self: *const SimplexEngine, problem: problem_module.ProblemView, status: SolveStatus) ?solution_module.SolutionView {
        const basis = if (self.basis) |*value| value else return null;
        return .{
            .status = status,
            .primal = basis.primal[0..problem.num_cols],
            .dual = basis.dual,
            .reduced_cost = basis.reduced_cost[0..problem.num_cols],
            .objective_value = self.objective_value,
            .iterations = self.iterations,
        };
    }
};

test {
    std.testing.refAllDecls(@This());
}

test "engine solves a standard-form LP with primal revised simplex" {
    const rows = [_]foundation.RowId{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) };
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &[_]f64{ -1, -1 },
        .col_lower = &[_]f64{ 0, 0 },
        .col_upper = &[_]f64{ std.math.inf(f64), std.math.inf(f64) },
        .row_lower = &[_]f64{-std.math.inf(f64)},
        .row_upper = &[_]f64{4},
        .matrix = .{ .num_rows = 1, .num_cols = 2, .col_starts = &[_]usize{ 0, 1, 2 }, .row_indices = &rows, .values = &[_]f64{ 1, 2 } },
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    const status = engine.solveProblem(problem, .{});
    try std.testing.expectEqual(SolveStatus.optimal, status);
    try std.testing.expectApproxEqAbs(@as(f64, 4), engine.basis.?.primal[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0), engine.basis.?.primal[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, -4), engine.objective_value, 1e-12);
    const solution = engine.solutionView(problem, status).?;
    try std.testing.expectEqual(@as(usize, 2), solution.primal.len);
    try std.testing.expectApproxEqAbs(@as(f64, -4), solution.objective_value, 1e-12);
}

test "engine solves with a nonzero structural lower bound" {
    const rows = [_]foundation.RowId{foundation.RowId.fromUsizeAssumeValid(0)};
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{1},
        .col_lower = &[_]f64{2},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{-std.math.inf(f64)},
        .row_upper = &[_]f64{5},
        .matrix = .{ .num_rows = 1, .num_cols = 1, .col_starts = &[_]usize{ 0, 1 }, .row_indices = &rows, .values = &[_]f64{1} },
        .objective_sense = .maximize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    const status = engine.solveProblem(problem, .{});
    try std.testing.expectEqual(SolveStatus.optimal, status);
    try std.testing.expectApproxEqAbs(@as(f64, 5), engine.basis.?.primal[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 5), engine.objective_value, 1e-12);
}

test "engine flips a nonbasic variable to its finite upper bound" {
    const rows = [_]foundation.RowId{foundation.RowId.fromUsizeAssumeValid(0)};
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{1},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{3},
        .row_lower = &[_]f64{-std.math.inf(f64)},
        .row_upper = &[_]f64{10},
        .matrix = .{ .num_rows = 1, .num_cols = 1, .col_starts = &[_]usize{ 0, 1 }, .row_indices = &rows, .values = &[_]f64{1} },
        .objective_sense = .maximize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{}));
    try std.testing.expectApproxEqAbs(@as(f64, 3), engine.basis.?.primal[0], 1e-12);
    try std.testing.expectEqual(basis_module.BasisStatus.at_upper, engine.basis.?.col_status[0]);
}

test "engine pivots a variable downward from its upper bound" {
    const rows = [_]foundation.RowId{foundation.RowId.fromUsizeAssumeValid(0)};
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{1},
        .col_lower = &[_]f64{-std.math.inf(f64)},
        .col_upper = &[_]f64{5},
        .row_lower = &[_]f64{2},
        .row_upper = &[_]f64{std.math.inf(f64)},
        .matrix = .{ .num_rows = 1, .num_cols = 1, .col_starts = &[_]usize{ 0, 1 }, .row_indices = &rows, .values = &[_]f64{1} },
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{}));
    try std.testing.expectApproxEqAbs(@as(f64, 2), engine.basis.?.primal[0], 1e-12);
    try std.testing.expectEqual(basis_module.BasisStatus.basic, engine.basis.?.col_status[0]);
}

test "engine normalizes and solves a ranged row" {
    const rows = [_]foundation.RowId{foundation.RowId.fromUsizeAssumeValid(0)};
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{1},
        .col_lower = &[_]f64{1},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{1},
        .row_upper = &[_]f64{5},
        .matrix = .{ .num_rows = 1, .num_cols = 1, .col_starts = &[_]usize{ 0, 1 }, .row_indices = &rows, .values = &[_]f64{1} },
        .objective_sense = .maximize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{}));
    try std.testing.expectApproxEqAbs(@as(f64, 5), engine.basis.?.primal[0], 1e-12);
}

test "Phase I restores feasibility for a greater-equal row" {
    const rows = [_]foundation.RowId{foundation.RowId.fromUsizeAssumeValid(0)};
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{1},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{2},
        .row_upper = &[_]f64{std.math.inf(f64)},
        .matrix = .{ .num_rows = 1, .num_cols = 1, .col_starts = &[_]usize{ 0, 1 }, .row_indices = &rows, .values = &[_]f64{1} },
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{}));
    try std.testing.expectApproxEqAbs(@as(f64, 2), engine.basis.?.primal[0], 1e-12);
}

test "Phase I removes zero artificial basics after redundant equalities" {
    const rows = [_]foundation.RowId{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1) };
    const problem = problem_module.ProblemView{
        .num_rows = 2,
        .num_cols = 1,
        .col_cost = &[_]f64{0},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{ 1, 1 },
        .row_upper = &[_]f64{ 1, 1 },
        .matrix = .{ .num_rows = 2, .num_cols = 1, .col_starts = &[_]usize{ 0, 2 }, .row_indices = &rows, .values = &[_]f64{ 1, 1 } },
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{}));
    try std.testing.expectApproxEqAbs(@as(f64, 1), engine.basis.?.primal[0], 1e-12);
    const artificial_begin = problem.num_cols + problem.num_rows;
    for (engine.basis.?.basic_index) |column| try std.testing.expect(column < artificial_begin);
}

test "Phase I detects an infeasible LP" {
    const rows = [_]foundation.RowId{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1) };
    const problem = problem_module.ProblemView{
        .num_rows = 2,
        .num_cols = 1,
        .col_cost = &[_]f64{0},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{ 2, -std.math.inf(f64) },
        .row_upper = &[_]f64{ std.math.inf(f64), 1 },
        .matrix = .{ .num_rows = 2, .num_cols = 1, .col_starts = &[_]usize{ 0, 2 }, .row_indices = &rows, .values = &[_]f64{ 1, 1 } },
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.infeasible, engine.solveProblem(problem, .{}));
}

test "engine detects an unbounded improving structural column" {
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{-1},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{-std.math.inf(f64)},
        .row_upper = &[_]f64{1},
        .matrix = .{ .num_rows = 1, .num_cols = 1, .col_starts = &[_]usize{ 0, 0 }, .row_indices = &.{}, .values = &.{} },
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.unbounded, engine.solveProblem(problem, .{}));
}

test "Phase I respects a zero iteration limit" {
    const rows = [_]foundation.RowId{foundation.RowId.fromUsizeAssumeValid(0)};
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{0},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{2},
        .row_upper = &[_]f64{std.math.inf(f64)},
        .matrix = .{ .num_rows = 1, .num_cols = 1, .col_starts = &[_]usize{ 0, 1 }, .row_indices = &rows, .values = &[_]f64{1} },
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.iteration_limit, engine.solveProblem(problem, .{ .max_iterations = 0 }));
}

test "deterministic work limit spans Phase I and Phase II" {
    const rows = [_]foundation.RowId{foundation.RowId.fromUsizeAssumeValid(0)};
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{1},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{2},
        .row_upper = &[_]f64{std.math.inf(f64)},
        .matrix = .{ .num_rows = 1, .num_cols = 1, .col_starts = &[_]usize{ 0, 1 }, .row_indices = &rows, .values = &[_]f64{1} },
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.work_limit, engine.solveProblem(problem, .{ .work_limit = 1 }));
    try std.testing.expectEqual(@as(u64, 1), engine.work_used);
}

test "iteration callback can stop without allocating callback state" {
    const Context = struct {
        calls: usize = 0,
        last_work: u64 = 0,

        fn callback(event: ProgressEventView, context_ptr: ?*anyopaque) CallbackAction {
            const self: *@This() = @ptrCast(@alignCast(context_ptr.?));
            self.calls += 1;
            self.last_work = event.work_used;
            return .stop;
        }
    };
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{-1},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{1},
        .row_lower = &[_]f64{-std.math.inf(f64)},
        .row_upper = &[_]f64{1},
        .matrix = .{ .num_rows = 1, .num_cols = 1, .col_starts = &[_]usize{ 0, 0 }, .row_indices = &.{}, .values = &.{} },
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var context = Context{};
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.interrupted, engine.solveProblem(problem, .{
        .iteration_callback = Context.callback,
        .callback_user_data = &context,
    }));
    try std.testing.expectEqual(@as(usize, 1), context.calls);
    try std.testing.expectEqual(@as(u64, 0), context.last_work);
    try std.testing.expectEqual(@as(u64, 0), engine.work_used);
}

test "structured iteration logging obeys its deterministic interval" {
    const Context = struct {
        calls: usize = 0,
        fn log(_: ProgressEventView, context_ptr: ?*anyopaque) void {
            const self: *@This() = @ptrCast(@alignCast(context_ptr.?));
            self.calls += 1;
        }
    };
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{-1},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{1},
        .row_lower = &[_]f64{-std.math.inf(f64)},
        .row_upper = &[_]f64{1},
        .matrix = .{ .num_rows = 1, .num_cols = 1, .col_starts = &[_]usize{ 0, 0 }, .row_indices = &.{}, .values = &.{} },
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var context = Context{};
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{
        .log_level = .iterations,
        .log_callback = Context.log,
        .log_user_data = &context,
        .log_interval_work = 1,
    }));
    try std.testing.expect(context.calls >= 1);
}

test "engine rejects malformed borrowed problem dimensions" {
    const problem = problem_module.ProblemView{
        .num_rows = 0,
        .num_cols = 1,
        .col_cost = &.{},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{1},
        .row_lower = &.{},
        .row_upper = &.{},
        .matrix = .{ .num_rows = 0, .num_cols = 1, .col_starts = &[_]usize{ 0, 0 }, .row_indices = &.{}, .values = &.{} },
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.numerical_failure, engine.solveProblem(problem, .{}));
}

test "engine honors an immediate time limit" {
    const problem = problem_module.ProblemView{
        .num_rows = 0,
        .num_cols = 1,
        .col_cost = &[_]f64{0},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{1},
        .row_lower = &.{},
        .row_upper = &.{},
        .matrix = .{ .num_rows = 0, .num_cols = 1, .col_starts = &[_]usize{ 0, 0 }, .row_indices = &.{}, .values = &.{} },
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.time_limit, engine.solveProblem(problem, .{ .time_limit_ns = 0 }));
}

test "engine honors a caller-owned atomic interrupt flag" {
    const problem = problem_module.ProblemView{
        .num_rows = 0,
        .num_cols = 1,
        .col_cost = &[_]f64{0},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{1},
        .row_lower = &.{},
        .row_upper = &.{},
        .matrix = .{ .num_rows = 0, .num_cols = 1, .col_starts = &[_]usize{ 0, 0 }, .row_indices = &.{}, .values = &.{} },
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var interrupted = std.atomic.Value(bool).init(true);
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.interrupted, engine.solveProblem(problem, .{ .interrupt_flag = &interrupted }));
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
        .matrix = .{ .num_rows = 1, .num_cols = 1, .col_starts = &[_]usize{ 0, 1 }, .row_indices = &rows, .values = &[_]f64{1} },
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
}

test "neither-feasible imported basis is repaired without a crash restart" {
    const Context = struct {
        saw_repair: bool = false,
        fn callback(event: ProgressEventView, context_ptr: ?*anyopaque) CallbackAction {
            const self: *@This() = @ptrCast(@alignCast(context_ptr.?));
            self.saw_repair = self.saw_repair or event.phase == .dual_feasibility_repair;
            return .continue_solve;
        }
    };
    const rows = [_]foundation.RowId{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) };
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &[_]f64{ 0, -1 },
        .col_lower = &[_]f64{ 0, 0 },
        .col_upper = &[_]f64{ 1, std.math.inf(f64) },
        .row_lower = &[_]f64{2},
        .row_upper = &[_]f64{2},
        .matrix = .{ .num_rows = 1, .num_cols = 2, .col_starts = &[_]usize{ 0, 1, 2 }, .row_indices = &rows, .values = &[_]f64{ 1, 1 } },
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    const initial_basis = basis_snapshot_module.BasisView{
        .structural_status = &[_]basis_module.BasisStatus{ .basic, .at_lower },
        .logical_status = &[_]basis_module.BasisStatus{.fixed},
        .basic_index = &[_]u32{0},
    };
    var context = Context{};
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{
        .initial_basis = initial_basis,
        .iteration_callback = Context.callback,
        .callback_user_data = &context,
    }));
    try std.testing.expect(context.saw_repair);
    try std.testing.expectEqual(Algorithm.primal_revised, engine.algorithm);
    try std.testing.expectApproxEqAbs(@as(f64, 0), engine.basis.?.primal[0], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 2), engine.basis.?.primal[1], 1e-10);
}

test "exact dual steepest-edge weights use BTRAN row norms" {
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 2, 0);
    engine.basis.?.basic_value[0] = -1.0;
    engine.basis.?.basic_value[1] = -2.0;
    @memset(engine.basis.?.basic_lower, 0.0);
    @memset(engine.basis.?.basic_upper, std.math.inf(f64));
    try engine.factorization.factorize(2, &[_]f64{ 2, 0, 0, 0.5 });
    try std.testing.expectEqual(SolveStatus.optimal, engine.ensureExactDualEdgeWeights());
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), engine.basis.?.row_edge_weight[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), engine.basis.?.row_edge_weight[1], 1e-12);

    engine.pricing.rule = .steepest_edge;
    const choice = engine.pricing.chooseDualLeavingWeighted(
        engine.basis.?.basic_value,
        engine.basis.?.basic_lower,
        engine.basis.?.basic_upper,
        engine.basis.?.row_edge_weight,
        engine.numerical.primal_tolerance,
    ).?;
    try std.testing.expectEqual(@as(u32, 0), choice.row);
}

test "incremental dual steepest-edge recurrence matches the new inverse rows" {
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 2, 0);
    try engine.factorization.factorizeIdentity(2);
    engine.pricing.rule = .steepest_edge;
    engine.dual_edge_weights_valid = true;
    engine.dual_row_index = 0;
    engine.basis.?.dual_row[0] = 1.0;
    engine.basis.?.dual_row[1] = 0.0;
    engine.basis.?.pivot_direction[0] = 2.0;
    engine.basis.?.pivot_direction[1] = 1.0;
    try std.testing.expectEqual(SolveStatus.optimal, engine.updateDualSteepestEdgeWeights(0, 2.0));
    try std.testing.expectApproxEqAbs(@as(f64, 0.25), engine.basis.?.row_edge_weight[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1.25), engine.basis.?.row_edge_weight[1], 1e-12);
}

test "hyper-sparse dual candidate list retains the most attractive rows" {
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    engine.basis = try basis_module.BasisState.init(std.testing.allocator, 64, 0);
    @memset(engine.basis.?.basic_lower, 0.0);
    @memset(engine.basis.?.basic_upper, std.math.inf(f64));
    @memset(engine.basis.?.basic_value, 0.0);
    for (24..64) |row| engine.basis.?.basic_value[row] = -@as(f64, @floatFromInt(row - 23));
    engine.pricing.rule = .hyper_sparse;
    engine.dual_hyper_sparse_active = true;
    const choice = engine.chooseDualLeavingRow().?;
    try std.testing.expectEqual(@as(u32, 63), choice.row);
    try std.testing.expectEqual(@as(usize, 32), engine.dual_candidate_count);
    try std.testing.expect(engine.dual_candidate_cutoff > 0.0);
}
