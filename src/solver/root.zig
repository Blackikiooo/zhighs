//! Top-level solve orchestration.
//!
//! Owns the incumbent model and coordinates validation, presolve, solver
//! selection, postsolve, callbacks, limits, and final status publication.

const std = @import("std");
const foundation = @import("foundation");
const lp = @import("lp");
const matrix = @import("matrix");
const presolve = @import("presolve");

pub const LpSolveControl = lp.simplex.engine.SolveControl;
pub const LpSolveStatus = lp.simplex.solution.SolveStatus;
pub const LpProblemView = lp.ProblemView;
pub const LpObjectiveSense = lp.ObjectiveSense;
pub const LpBasisStatus = lp.simplex.basis.BasisStatus;
pub const LpBasisView = lp.BasisView;
pub const LpBasisSnapshot = lp.BasisSnapshot;
pub const LpProgressEventView = lp.simplex.engine.ProgressEventView;
pub const LpCallbackAction = lp.simplex.engine.CallbackAction;
pub const LpPricingRule = lp.simplex.pricing.PricingRule;

pub const LpRevisionStamp = struct {
    structure: u64,
    matrix_values: u64,
    bounds: u64,
    objective: u64,
};

/// Borrowed LP result projected from an owning `LpSolveResult`.
///
/// Every slice is invalidated by `LpSolveResult.deinit`.
pub const LpSolveResultView = struct {
    status: LpSolveStatus,
    primal_valid: bool,
    dual_valid: bool,
    basis_valid: bool,
    primal: []const f64,
    dual: []const f64,
    reduced_cost: []const f64,
    unbounded_ray: []const f64,
    infeasibility_ray: []const f64,
    structural_status: []const LpBasisStatus,
    logical_status: []const LpBasisStatus,
    basic_index: []const u32,
    row_scale: []const f64,
    objective_value: f64,
    iterations: usize,
};

/// Owning lifetime boundary around one simplex invocation.
pub const LpSolveResult = struct {
    engine: lp.simplex.engine.SimplexEngine,
    status: LpSolveStatus,
    presolved_problem: ?presolve.PresolvePipeline = null,
    projection: ?PostsolveProjection = null,

    pub fn deinit(self: *LpSolveResult) void {
        if (self.projection) |*value| value.deinit();
        if (self.presolved_problem) |*value| value.deinit();
        self.engine.deinit();
        self.* = undefined;
    }

    pub fn resultView(self: *const LpSolveResult, problem: LpProblemView) ?LpSolveResultView {
        if (self.status == .invalid_model) return null;
        if (self.projection) |*projection| return projection.view(self.engine.iterations);
        if (self.presolved_problem != null) return null;
        const solution = self.engine.solutionView(problem, self.status) orelse return null;
        const basis = if (self.engine.basis) |*value| value else return null;
        return .{
            .status = self.status,
            .primal_valid = self.status == .optimal,
            .dual_valid = self.status == .optimal,
            .basis_valid = basisIsValid(problem, basis),
            .primal = solution.primal,
            .dual = solution.dual,
            .reduced_cost = solution.reduced_cost,
            .unbounded_ray = solution.unbounded_ray,
            .infeasibility_ray = solution.infeasibility_ray,
            .structural_status = basis.col_status[0..problem.num_cols],
            .logical_status = basis.col_status[problem.num_cols..][0..problem.num_rows],
            .basic_index = basis.basic_index,
            .row_scale = basis.row_scale,
            .objective_value = solution.objective_value,
            .iterations = solution.iterations,
        };
    }
};

/// Execute the LP backend while borrowing all immutable problem arrays.
pub fn solveLp(allocator: std.mem.Allocator, problem: LpProblemView, control: LpSolveControl) LpSolveResult {
    var engine = lp.simplex.engine.SimplexEngine.init(allocator);
    problem.validate() catch {
        const status = engine.solveProblem(problem, control);
        return .{ .engine = engine, .status = status };
    };
    if (!hasPresolveReductions(allocator, problem)) {
        const status = engine.solveProblem(problem, control);
        return .{ .engine = engine, .status = status };
    }

    var pp = presolve.presolvePipeline(allocator, problem) catch
        return .{ .engine = engine, .status = .numerical_failure };
    if (!pp.wasApplied()) {
        pp.deinit();
        const status = engine.solveProblem(problem, control);
        return .{ .engine = engine, .status = status };
    }
    var reduced_basis: ?ReducedBasisChain = null;
    defer if (reduced_basis) |*value| value.deinit();
    var reduced_control = control;
    if (control.initial_basis) |initial_basis|
        reduced_basis = ReducedBasisChain.init(allocator, &pp, initial_basis) catch null;
    reduced_control.initial_basis = if (reduced_basis) |*value| value.view() else null;
    const reduced_problem = pp.reducedView();
    const status = engine.solveProblem(reduced_problem, reduced_control);
    const projection = PostsolveProjection.init(allocator, problem, &pp, &engine, status) catch null;
    return .{
        .engine = engine,
        .status = if (projection == null) .numerical_failure else status,
        .presolved_problem = pp,
        .projection = projection,
    };
}

/// Persistent LP solve lifetime. It owns all mutable numerical state but never
/// retains borrowed problem arrays between calls.
pub const LpSolveSession = struct {
    allocator: std.mem.Allocator,
    engine: lp.simplex.engine.SimplexEngine,
    presolved_problem: ?presolve.PresolvePipeline = null,
    projection: ?PostsolveProjection = null,
    last_revision: ?LpRevisionStamp = null,
    last_rows: usize = 0,
    last_cols: usize = 0,
    cold_solves: usize = 0,
    reoptimizations: usize = 0,

    pub fn init(allocator: std.mem.Allocator) LpSolveSession {
        return .{ .allocator = allocator, .engine = lp.simplex.engine.SimplexEngine.init(allocator) };
    }

    pub fn deinit(self: *LpSolveSession) void {
        self.clearPresolveState();
        self.engine.deinit();
        self.* = undefined;
    }

    pub fn solve(self: *LpSolveSession, problem: LpProblemView, control: LpSolveControl, revision: LpRevisionStamp) LpSolveStatus {
        if (self.projection) |*value| value.deinit();
        self.projection = null;

        problem.validate() catch {
            self.clearPresolveState();
            self.last_revision = null;
            self.last_rows = 0;
            self.last_cols = 0;
            return self.engine.solveProblem(problem, control);
        };

        if (!hasPresolveReductions(self.allocator, problem)) {
            if (self.presolved_problem) |*value| value.deinit();
            self.presolved_problem = null;
            return self.solveReduced(problem, control, revision, false);
        }

        var next_pp = presolve.presolvePipeline(self.allocator, problem) catch {
            self.last_revision = revision;
            return .numerical_failure;
        };
        if (!next_pp.wasApplied()) {
            next_pp.deinit();
            if (self.presolved_problem) |*value| value.deinit();
            self.presolved_problem = null;
            return self.solveReduced(problem, control, revision, false);
        }
        const mapping_changed = if (self.presolved_problem) |*previous|
            !previous.sameMapping(&next_pp)
        else
            true;
        if (self.presolved_problem) |*value| value.deinit();
        self.presolved_problem = next_pp;

        var reduced_control = control;
        var reduced_basis: ?ReducedBasisChain = null;
        defer if (reduced_basis) |*value| value.deinit();
        if (control.initial_basis) |initial_basis|
            reduced_basis = ReducedBasisChain.init(self.allocator, &self.presolved_problem.?, initial_basis) catch null;
        reduced_control.initial_basis = if (reduced_basis) |*value| value.view() else null;
        const reduced_problem = self.presolved_problem.?.reducedView();
        const status = self.solveReduced(reduced_problem, reduced_control, revision, mapping_changed);
        self.projection = PostsolveProjection.init(
            self.allocator,
            problem,
            &self.presolved_problem.?,
            &self.engine,
            status,
        ) catch {
            return .numerical_failure;
        };
        return status;
    }

    fn solveReduced(self: *LpSolveSession, problem: LpProblemView, control: LpSolveControl, revision: LpRevisionStamp, force_cold: bool) LpSolveStatus {
        const reusable = !force_cold and if (self.last_revision) |last|
            last.structure == revision.structure and last.matrix_values == revision.matrix_values and
                self.last_rows == problem.num_rows and self.last_cols == problem.num_cols
        else
            false;
        const status = if (reusable) blk: {
            self.reoptimizations += 1;
            break :blk self.engine.reoptimizeProblem(problem, control);
        } else blk: {
            self.cold_solves += 1;
            break :blk self.engine.solveProblem(problem, control);
        };
        self.last_revision = revision;
        self.last_rows = problem.num_rows;
        self.last_cols = problem.num_cols;
        return status;
    }

    pub fn resultView(self: *const LpSolveSession, problem: LpProblemView, status: LpSolveStatus) ?LpSolveResultView {
        if (status == .invalid_model) return null;
        if (self.projection) |*projection| return projection.view(self.engine.iterations);
        if (self.presolved_problem != null) return null;
        const solution = self.engine.solutionView(problem, status) orelse return null;
        const basis = if (self.engine.basis) |*value| value else return null;
        return .{
            .status = status,
            .primal_valid = status == .optimal,
            .dual_valid = status == .optimal,
            .basis_valid = basisIsValid(problem, basis),
            .primal = solution.primal,
            .dual = solution.dual,
            .reduced_cost = solution.reduced_cost,
            .unbounded_ray = solution.unbounded_ray,
            .infeasibility_ray = solution.infeasibility_ray,
            .structural_status = basis.col_status[0..problem.num_cols],
            .logical_status = basis.col_status[problem.num_cols..][0..problem.num_rows],
            .basic_index = basis.basic_index,
            .row_scale = basis.row_scale,
            .objective_value = solution.objective_value,
            .iterations = solution.iterations,
        };
    }

    fn clearPresolveState(self: *LpSolveSession) void {
        if (self.projection) |*value| value.deinit();
        if (self.presolved_problem) |*value| value.deinit();
        self.projection = null;
        self.presolved_problem = null;
    }
};

fn hasPresolveReductions(allocator: std.mem.Allocator, problem: LpProblemView) bool {
    const objective_sign: f64 = @floatFromInt(@intFromEnum(problem.objective_sense));
    for (0..problem.num_cols) |column| {
        const begin = problem.matrix.col_starts[column];
        const end = problem.matrix.col_starts[column + 1];
        switch (presolve.rules.column_substitution.classify(
            end - begin,
            problem.col_cost[column],
            problem.col_lower[column],
            problem.col_upper[column],
            objective_sign,
        )) {
            .keep => {},
            .substitute => return true,
        }
        if (presolve.rules.column_substitution.isFreeSingleton(
            end - begin,
            problem.col_cost[column],
            problem.col_lower[column],
            problem.col_upper[column],
        )) return true;
    }

    // CSC has no row counts. Count them once to select empty- and singleton-row
    // presolve. Allocation failure safely leaves either row type to simplex.
    const row_counts = allocator.alloc(u32, problem.num_rows) catch return false;
    defer allocator.free(row_counts);
    @memset(row_counts, 0);
    for (problem.matrix.row_indices) |row| row_counts[row.toUsize()] += 1;
    for (row_counts, problem.row_lower, problem.row_upper) |count, lower, upper| {
        if (count == 1) return true;
        if (presolve.rules.empty_row.classify(count, lower, upper) == .remove) return true;
    }
    if (hasExactDuplicateRows(allocator, problem, row_counts)) return true;

    const row_activity = allocator.alloc(presolve.rules.redundant_row.ActivityBounds, problem.num_rows) catch return false;
    defer allocator.free(row_activity);
    for (row_activity) |*activity| activity.* = .{};
    for (0..problem.num_cols) |column| {
        const begin = problem.matrix.col_starts[column];
        const end = problem.matrix.col_starts[column + 1];
        for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient|
            row_activity[row.toUsize()].add(coefficient, problem.col_lower[column], problem.col_upper[column]);
    }
    for (row_activity, problem.row_lower, problem.row_upper) |activity, lower, upper| {
        if (presolve.rules.redundant_row.classify(activity, lower, upper) == .remove) return true;
    }
    return false;
}

fn hasExactDuplicateRows(allocator: std.mem.Allocator, problem: LpProblemView, row_counts: []const u32) bool {
    const row_starts = allocator.alloc(usize, problem.num_rows + 1) catch return false;
    defer allocator.free(row_starts);
    row_starts[0] = 0;
    for (row_counts, 0..) |count, row|
        row_starts[row + 1] = row_starts[row] + count;
    const row_columns = allocator.alloc(u32, problem.matrix.values.len) catch return false;
    defer allocator.free(row_columns);
    const row_values = allocator.alloc(f64, problem.matrix.values.len) catch return false;
    defer allocator.free(row_values);
    const cursor = allocator.dupe(usize, row_starts[0..problem.num_rows]) catch return false;
    defer allocator.free(cursor);
    for (0..problem.num_cols) |column| {
        const begin = problem.matrix.col_starts[column];
        const end = problem.matrix.col_starts[column + 1];
        for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row_id, value| {
            const row = row_id.toUsize();
            row_columns[cursor[row]] = @intCast(column);
            row_values[cursor[row]] = value;
            cursor[row] += 1;
        }
    }

    var heads: std.AutoHashMapUnmanaged(u64, u32) = .empty;
    defer heads.deinit(allocator);
    const next = allocator.alloc(u32, problem.num_rows) catch return false;
    defer allocator.free(next);
    @memset(next, std.math.maxInt(u32));
    for (0..problem.num_rows) |row| {
        const begin = row_starts[row];
        const end = row_starts[row + 1];
        const signature = presolve.rules.duplicate_row.hash(
            problem.row_lower[row],
            problem.row_upper[row],
            row_columns[begin..end],
            row_values[begin..end],
        );
        var candidate = heads.get(signature) orelse std.math.maxInt(u32);
        while (candidate != std.math.maxInt(u32)) : (candidate = next[candidate]) {
            const candidate_begin = row_starts[candidate];
            const candidate_end = row_starts[candidate + 1];
            if (presolve.rules.duplicate_row.exactlyEqual(
                problem.row_lower[row],
                problem.row_upper[row],
                row_columns[begin..end],
                row_values[begin..end],
                problem.row_lower[candidate],
                problem.row_upper[candidate],
                row_columns[candidate_begin..candidate_end],
                row_values[candidate_begin..candidate_end],
            )) return true;
        }
        next[row] = heads.get(signature) orelse std.math.maxInt(u32);
        heads.put(allocator, signature, @intCast(row)) catch return false;
    }
    return false;
}

fn basisIsValid(problem: LpProblemView, basis: *const lp.simplex.basis.BasisState) bool {
    const view = LpBasisView{
        .structural_status = basis.col_status[0..problem.num_cols],
        .logical_status = basis.col_status[problem.num_cols..][0..problem.num_rows],
        .basic_index = basis.basic_index,
    };
    view.validate(problem.num_cols, problem.num_rows) catch return false;
    return true;
}

/// Temporary original-to-reduced basis mapping. It is borrowed by simplex
/// only for the duration of one solve call.
const ReducedBasisProjection = struct {
    allocator: std.mem.Allocator,
    status_storage: []LpBasisStatus,
    structural_status: []LpBasisStatus,
    logical_status: []LpBasisStatus,
    basic_index: []u32,

    fn init(
        allocator: std.mem.Allocator,
        pp: *const presolve.PresolvedProblem,
        original: LpBasisView,
    ) !ReducedBasisProjection {
        try original.validate(pp.original_num_cols, pp.original_num_rows);
        var result = ReducedBasisProjection{
            .allocator = allocator,
            .status_storage = &.{},
            .structural_status = &.{},
            .logical_status = &.{},
            .basic_index = &.{},
        };
        errdefer result.deinit();
        result.status_storage = try allocator.alloc(LpBasisStatus, pp.num_cols + pp.num_rows);
        result.structural_status = result.status_storage[0..pp.num_cols];
        result.logical_status = result.status_storage[pp.num_cols..];
        result.basic_index = try allocator.alloc(u32, pp.num_rows);

        for (pp.reduced_to_original_col, 0..) |original_col, reduced_col|
            result.structural_status[reduced_col] = original.structural_status[original_col];
        for (pp.reduced_to_original_row, 0..) |original_row, reduced_row|
            result.logical_status[reduced_row] = original.logical_status[original_row];
        for (pp.original_to_reduced_row, 0..) |reduced_row, original_row| {
            if (reduced_row != std.math.maxInt(u32)) continue;
            var free_singleton_column: ?u32 = null;
            for (pp.free_singleton_columns[0..pp.free_singleton_column_count]) |record| {
                if (record.original_row == original_row) {
                    free_singleton_column = record.original_column;
                    break;
                }
            }
            if (free_singleton_column) |column| {
                if (original.structural_status[column] != .basic or
                    original.logical_status[original_row] == .basic or
                    original.basic_index[original_row] != column)
                    return error.RemovedFreeSingletonBasisMismatch;
            } else if (original.logical_status[original_row] != .basic or
                original.basic_index[original_row] != pp.original_num_cols + original_row)
                return error.RemovedRowNotLogicalBasic;
        }
        for (pp.reduced_to_original_row, 0..) |original_row, row| {
            const original_col = original.basic_index[original_row];
            if (original_col < pp.original_num_cols) {
                const reduced_col = pp.original_to_reduced_col[original_col];
                if (reduced_col == std.math.maxInt(u32)) return error.SubstitutedBasicColumn;
                result.basic_index[row] = reduced_col;
            } else {
                const original_logical_row = original_col - pp.original_num_cols;
                const reduced_logical_row = pp.original_to_reduced_row[original_logical_row];
                if (reduced_logical_row == std.math.maxInt(u32)) return error.RemovedLogicalBasicColumn;
                result.basic_index[row] = pp.num_cols + reduced_logical_row;
            }
        }
        try result.view().validate(pp.num_cols, pp.num_rows);
        return result;
    }

    fn deinit(self: *ReducedBasisProjection) void {
        self.allocator.free(self.status_storage);
        self.allocator.free(self.basic_index);
        self.* = undefined;
    }

    fn view(self: *const ReducedBasisProjection) LpBasisView {
        return .{
            .structural_status = self.structural_status,
            .logical_status = self.logical_status,
            .basic_index = self.basic_index,
        };
    }
};

/// Owns every forward basis projection so each borrowed intermediate view
/// remains alive until simplex has consumed the final reduced basis.
const ReducedBasisChain = struct {
    allocator: std.mem.Allocator,
    projections: []ReducedBasisProjection,
    initialized: usize = 0,

    fn init(
        allocator: std.mem.Allocator,
        pipeline: *const presolve.PresolvePipeline,
        original: LpBasisView,
    ) !ReducedBasisChain {
        var result = ReducedBasisChain{
            .allocator = allocator,
            .projections = try allocator.alloc(ReducedBasisProjection, pipeline.passes.items.len),
        };
        errdefer result.deinit();
        var current = original;
        for (pipeline.passes.items) |*pass| {
            result.projections[result.initialized] = try ReducedBasisProjection.init(allocator, pass, current);
            result.initialized += 1;
            current = result.projections[result.initialized - 1].view();
        }
        return result;
    }

    fn deinit(self: *ReducedBasisChain) void {
        for (self.projections[0..self.initialized]) |*projection| projection.deinit();
        self.allocator.free(self.projections);
        self.* = undefined;
    }

    fn view(self: *const ReducedBasisChain) LpBasisView {
        std.debug.assert(self.initialized != 0);
        return self.projections[self.initialized - 1].view();
    }
};

/// One owning reverse basis projection across a relative presolve pass.
const ExpandedBasis = struct {
    allocator: std.mem.Allocator,
    structural_status: []LpBasisStatus = &.{},
    logical_status: []LpBasisStatus = &.{},
    basic_index: []u32 = &.{},
    row_scale: []f64 = &.{},
    basis_valid: bool = true,

    fn init(
        allocator: std.mem.Allocator,
        pp: *const presolve.PresolvedProblem,
        reduced_structural_status: []const LpBasisStatus,
        reduced_logical_status: []const LpBasisStatus,
        reduced_basic_index: []const u32,
        reduced_row_scale: []const f64,
        input_valid: bool,
    ) !ExpandedBasis {
        var result = ExpandedBasis{ .allocator = allocator, .basis_valid = input_valid };
        errdefer result.deinit();
        result.structural_status = try allocator.alloc(LpBasisStatus, pp.original_num_cols);
        result.logical_status = try allocator.alloc(LpBasisStatus, pp.original_num_rows);
        result.basic_index = try allocator.alloc(u32, pp.original_num_rows);
        result.row_scale = try allocator.alloc(f64, pp.original_num_rows);

        @memset(result.structural_status, .at_lower);
        for (pp.column_substitutions[0..pp.column_substitution_count]) |record|
            result.structural_status[record.original_column] = record.basis_status;
        for (pp.reduced_to_original_col, 0..) |original_col, reduced_col|
            result.structural_status[original_col] = reduced_structural_status[reduced_col];

        @memset(result.logical_status, .basic);
        @memset(result.row_scale, 1.0);
        for (pp.reduced_to_original_row, 0..) |original_row, reduced_row| {
            result.logical_status[original_row] = reduced_logical_status[reduced_row];
            result.row_scale[original_row] = reduced_row_scale[reduced_row];
        }
        for (0..pp.original_num_rows) |row| {
            if (pp.original_to_reduced_row[row] == std.math.maxInt(u32))
                result.basic_index[row] = pp.original_num_cols + @as(u32, @intCast(row));
        }
        for (pp.free_singleton_columns[0..pp.free_singleton_column_count]) |record| {
            result.structural_status[record.original_column] = .basic;
            result.logical_status[record.original_row] = switch (record.side) {
                .fixed => .fixed,
                .lower => .at_upper,
                .upper => .at_lower,
                .free => .free,
            };
            result.basic_index[record.original_row] = record.original_column;
        }
        for (reduced_basic_index, 0..) |reduced_col, reduced_row| {
            const original_row = pp.reduced_to_original_row[reduced_row];
            result.basic_index[original_row] = if (reduced_col < pp.num_cols)
                pp.reduced_to_original_col[reduced_col]
            else if (reduced_col < pp.num_cols + pp.num_rows)
                pp.original_num_cols + pp.reduced_to_original_row[reduced_col - pp.num_cols]
            else blk: {
                result.basis_valid = false;
                break :blk std.math.maxInt(u32);
            };
        }
        var singleton_index: usize = pp.singleton_row_count;
        while (singleton_index > 0) {
            singleton_index -= 1;
            const record = pp.singleton_rows[singleton_index];
            const structural_status = result.structural_status[record.original_column];
            if (structural_status == .basic) continue;
            const source = switch (structural_status) {
                .at_lower => record.lower_source,
                .at_upper => record.upper_source,
                .fixed => record.lower_source orelse record.upper_source,
                .basic, .free, .superbasic => null,
            } orelse continue;
            result.structural_status[record.original_column] = .basic;
            result.logical_status[record.original_row] = switch (source) {
                .lower => .at_upper,
                .upper => .at_lower,
            };
            result.basic_index[record.original_row] = record.original_column;
        }
        const mapped_basis = LpBasisView{
            .structural_status = result.structural_status,
            .logical_status = result.logical_status,
            .basic_index = result.basic_index,
        };
        mapped_basis.validate(pp.original_num_cols, pp.original_num_rows) catch {
            result.basis_valid = false;
        };
        return result;
    }

    fn deinit(self: *ExpandedBasis) void {
        self.allocator.free(self.structural_status);
        self.allocator.free(self.logical_status);
        self.allocator.free(self.basic_index);
        self.allocator.free(self.row_scale);
        self.* = undefined;
    }
};

/// Owning original-coordinate result produced after solving a reduced LP.
const PostsolveProjection = struct {
    allocator: std.mem.Allocator,
    solution: presolve.PresolveResult,
    structural_status: []LpBasisStatus,
    logical_status: []LpBasisStatus,
    basic_index: []u32,
    row_scale: []f64,
    basis_valid: bool = true,

    fn init(
        allocator: std.mem.Allocator,
        original_problem: LpProblemView,
        pipeline: *const presolve.PresolvePipeline,
        engine: *const lp.simplex.engine.SimplexEngine,
        status: LpSolveStatus,
    ) !PostsolveProjection {
        const reduced_problem = pipeline.reducedView();
        const reduced_solution = engine.solutionView(reduced_problem, status) orelse return error.SolutionUnavailable;
        var result = PostsolveProjection{
            .allocator = allocator,
            .solution = try pipeline.postsolveAll(reduced_solution, allocator),
            .structural_status = &.{},
            .logical_status = &.{},
            .basic_index = &.{},
            .row_scale = &.{},
        };
        errdefer result.deinit();
        const basis = if (engine.basis) |*value| value else return error.BasisUnavailable;
        var current: ?ExpandedBasis = null;
        errdefer if (current) |*value| value.deinit();
        var structural_status: []const LpBasisStatus = basis.col_status[0..reduced_problem.num_cols];
        var logical_status: []const LpBasisStatus = basis.col_status[reduced_problem.num_cols..][0..reduced_problem.num_rows];
        var basic_index: []const u32 = basis.basic_index[0..reduced_problem.num_rows];
        var row_scale: []const f64 = basis.row_scale[0..reduced_problem.num_rows];
        var basis_valid = basisIsValid(reduced_problem, basis);
        var pass_index = pipeline.passes.items.len;
        while (pass_index > 0) {
            pass_index -= 1;
            const next = try ExpandedBasis.init(
                allocator,
                &pipeline.passes.items[pass_index],
                structural_status,
                logical_status,
                basic_index,
                row_scale,
                basis_valid,
            );
            if (current) |*value| value.deinit();
            current = next;
            structural_status = current.?.structural_status;
            logical_status = current.?.logical_status;
            basic_index = current.?.basic_index;
            row_scale = current.?.row_scale;
            basis_valid = current.?.basis_valid;
        }
        var expanded = current.?;
        current = null;
        result.structural_status = expanded.structural_status;
        result.logical_status = expanded.logical_status;
        result.basic_index = expanded.basic_index;
        result.row_scale = expanded.row_scale;
        result.basis_valid = expanded.basis_valid;
        expanded.structural_status = &.{};
        expanded.logical_status = &.{};
        expanded.basic_index = &.{};
        expanded.row_scale = &.{};
        std.debug.assert(result.structural_status.len == original_problem.num_cols);
        std.debug.assert(result.logical_status.len == original_problem.num_rows);
        return result;
    }

    fn deinit(self: *PostsolveProjection) void {
        self.solution.deinit();
        self.allocator.free(self.structural_status);
        self.allocator.free(self.logical_status);
        self.allocator.free(self.basic_index);
        self.allocator.free(self.row_scale);
        self.* = undefined;
    }

    fn view(self: *const PostsolveProjection, iterations: usize) LpSolveResultView {
        return .{
            .status = self.solution.status,
            .primal_valid = self.solution.status == .optimal,
            .dual_valid = self.solution.status == .optimal,
            .basis_valid = self.basis_valid,
            .primal = self.solution.primal,
            .dual = self.solution.dual,
            .reduced_cost = self.solution.reduced_cost,
            .unbounded_ray = self.solution.unbounded_ray,
            .infeasibility_ray = self.solution.infeasibility_ray,
            .structural_status = self.structural_status,
            .logical_status = self.logical_status,
            .basic_index = self.basic_index,
            .row_scale = self.row_scale,
            .objective_value = self.solution.objective_value,
            .iterations = iterations,
        };
    }
};

fn fixedColumnProblem(fixed_first: bool) LpProblemView {
    return .{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 1.0, 2.0 },
        .col_lower = if (fixed_first) &.{ 2.0, 0.0 } else &.{ 0.0, 3.0 },
        .col_upper = if (fixed_first) &.{ 2.0, std.math.inf(f64) } else &.{ std.math.inf(f64), 3.0 },
        .row_lower = &.{5.0},
        .row_upper = &.{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(
            1,
            2,
            &[_]usize{ 0, 1, 2 },
            &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) },
            &.{ 1.0, 1.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
}

fn expectValidFarkas(problem: LpProblemView, ray: []const f64) !void {
    try std.testing.expectEqual(problem.num_rows, ray.len);
    var proof_lower: f64 = 0.0;
    for (ray, problem.row_lower, problem.row_upper) |multiplier, lower, upper| {
        if (multiplier > 0.0) {
            try std.testing.expect(std.math.isFinite(lower));
            proof_lower += multiplier * lower;
        } else if (multiplier < 0.0) {
            try std.testing.expect(std.math.isFinite(upper));
            proof_lower += multiplier * upper;
        }
    }
    var implied_upper: f64 = 0.0;
    for (0..problem.num_cols) |column| {
        const begin = problem.matrix.col_starts[column];
        const end = problem.matrix.col_starts[column + 1];
        var coefficient: f64 = 0.0;
        for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, value|
            coefficient += ray[row.toUsize()] * value;
        if (coefficient > 0.0) {
            try std.testing.expect(std.math.isFinite(problem.col_upper[column]));
            implied_upper += coefficient * problem.col_upper[column];
        } else if (coefficient < 0.0) {
            try std.testing.expect(std.math.isFinite(problem.col_lower[column]));
            implied_upper += coefficient * problem.col_lower[column];
        }
    }
    try std.testing.expect(proof_lower - implied_upper > 1e-9);
}

test "solveLp returns fixed-column solutions and basis in original coordinates" {
    const problem = fixedColumnProblem(false);
    var result = solveLp(std.testing.allocator, problem, .{});
    defer result.deinit();
    try std.testing.expectEqual(LpSolveStatus.optimal, result.status);
    const view = result.resultView(problem).?;
    try std.testing.expect(view.primal_valid);
    try std.testing.expect(view.dual_valid);
    try std.testing.expect(view.basis_valid);
    try std.testing.expectEqualSlices(f64, &.{ 2.0, 3.0 }, view.primal);
    try std.testing.expectApproxEqAbs(@as(f64, 8.0), view.objective_value, 1e-9);
    try std.testing.expectEqual(LpBasisStatus.fixed, view.structural_status[1]);
    try std.testing.expectEqual(@as(usize, 2), view.structural_status.len);
    try std.testing.expectEqual(@as(usize, 1), view.logical_status.len);
    try std.testing.expect(view.basic_index[0] < problem.num_cols + problem.num_rows);
}

test "solveLp maps an exported original basis back through presolve" {
    const problem = fixedColumnProblem(false);
    var first = solveLp(std.testing.allocator, problem, .{});
    defer first.deinit();
    const first_view = first.resultView(problem).?;
    var snapshot = try LpBasisSnapshot.initFromView(std.testing.allocator, .{
        .structural_status = first_view.structural_status,
        .logical_status = first_view.logical_status,
        .basic_index = first_view.basic_index,
    });
    defer snapshot.deinit();

    var second = solveLp(std.testing.allocator, problem, .{ .initial_basis = snapshot.view() });
    defer second.deinit();
    try std.testing.expectEqual(LpSolveStatus.optimal, second.status);
    const second_view = second.resultView(problem).?;
    try std.testing.expect(second_view.basis_valid);
    try std.testing.expectEqualSlices(f64, first_view.primal, second_view.primal);
    try std.testing.expectEqual(@as(usize, 0), second_view.iterations);
}

test "presolve restores and warm-starts a redundant row exposed by a fixed column" {
    const problem = LpProblemView{
        .num_rows = 2,
        .num_cols = 2,
        .col_cost = &.{ 1.0, 2.0 },
        .col_lower = &.{ 0.0, 2.0 },
        .col_upper = &.{ std.math.inf(f64), 2.0 },
        .row_lower = &.{ 2.0, 1.0 },
        .row_upper = &.{ 2.0, std.math.inf(f64) },
        .matrix = matrix.CscView.initAssumeValid(
            2,
            2,
            &[_]usize{ 0, 1, 2 },
            &.{ foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(0) },
            &.{ 1.0, 1.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var first = solveLp(std.testing.allocator, problem, .{});
    defer first.deinit();
    try std.testing.expectEqual(LpSolveStatus.optimal, first.status);
    const first_view = first.resultView(problem).?;
    try std.testing.expect(first_view.basis_valid);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 2.0 }, first_view.primal);
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 1.0 }, first_view.dual);
    try std.testing.expectEqual(LpBasisStatus.basic, first_view.logical_status[0]);
    try std.testing.expectEqual(@as(u32, 2), first_view.basic_index[0]);

    var snapshot = try LpBasisSnapshot.initFromView(std.testing.allocator, .{
        .structural_status = first_view.structural_status,
        .logical_status = first_view.logical_status,
        .basic_index = first_view.basic_index,
    });
    defer snapshot.deinit();
    var second = solveLp(std.testing.allocator, problem, .{ .initial_basis = snapshot.view() });
    defer second.deinit();
    try std.testing.expectEqual(LpSolveStatus.optimal, second.status);
    const second_view = second.resultView(problem).?;
    try std.testing.expect(second_view.basis_valid);
    try std.testing.expectEqual(@as(usize, 0), second_view.iterations);
    try std.testing.expectEqualSlices(f64, first_view.primal, second_view.primal);
}

test "solver routes an original redundant empty row through reversible presolve" {
    const problem = LpProblemView{
        .num_rows = 2,
        .num_cols = 1,
        .col_cost = &.{1.0},
        .col_lower = &.{0.0},
        .col_upper = &.{std.math.inf(f64)},
        .row_lower = &.{ -1.0, 2.0 },
        .row_upper = &.{ 1.0, std.math.inf(f64) },
        .matrix = matrix.CscView.initAssumeValid(
            2,
            1,
            &[_]usize{ 0, 1 },
            &.{foundation.RowId.fromUsizeAssumeValid(1)},
            &.{1.0},
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var result = solveLp(std.testing.allocator, problem, .{});
    defer result.deinit();
    try std.testing.expectEqual(LpSolveStatus.optimal, result.status);
    try std.testing.expect(result.presolved_problem.?.removed_row_count == 2);
    try std.testing.expect(result.presolved_problem.?.singleton_row_count == 1);
    const view = result.resultView(problem).?;
    try std.testing.expect(view.basis_valid);
    try std.testing.expectEqualSlices(f64, &.{2.0}, view.primal);
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 1.0 }, view.dual);
    try std.testing.expectEqual(LpBasisStatus.basic, view.logical_status[0]);
    try std.testing.expectEqual(@as(u32, 1), view.basic_index[0]);
    try std.testing.expectEqual(LpBasisStatus.basic, view.structural_status[0]);
    try std.testing.expectEqual(LpBasisStatus.at_upper, view.logical_status[1]);
    try std.testing.expectEqual(@as(u32, 0), view.basic_index[1]);
}

test "activity presolve removes a general redundant row and warm-starts postsolve" {
    const problem = LpProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 1.0, 2.0 },
        .col_lower = &.{ 0.0, -2.0 },
        .col_upper = &.{ 1.0, 3.0 },
        .row_lower = &.{-3.0},
        .row_upper = &.{6.0},
        .matrix = matrix.CscView.initAssumeValid(
            1,
            2,
            &[_]usize{ 0, 1, 2 },
            &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) },
            &.{ 2.0, 1.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var first = solveLp(std.testing.allocator, problem, .{});
    defer first.deinit();
    try std.testing.expectEqual(LpSolveStatus.optimal, first.status);
    try std.testing.expectEqual(@as(u32, 1), first.presolved_problem.?.redundant_row_count);
    try std.testing.expectEqual(@as(usize, 2), first.presolved_problem.?.passes.items.len);
    const first_view = first.resultView(problem).?;
    try std.testing.expect(first_view.basis_valid);
    try std.testing.expectEqualSlices(f64, &.{ 0.0, -2.0 }, first_view.primal);
    try std.testing.expectEqualSlices(f64, &.{0.0}, first_view.dual);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 2.0 }, first_view.reduced_cost);
    try std.testing.expectEqual(LpBasisStatus.basic, first_view.logical_status[0]);
    try std.testing.expectEqual(@as(u32, 2), first_view.basic_index[0]);

    var snapshot = try LpBasisSnapshot.initFromView(std.testing.allocator, .{
        .structural_status = first_view.structural_status,
        .logical_status = first_view.logical_status,
        .basic_index = first_view.basic_index,
    });
    defer snapshot.deinit();
    var second = solveLp(std.testing.allocator, problem, .{ .initial_basis = snapshot.view() });
    defer second.deinit();
    try std.testing.expectEqual(LpSolveStatus.optimal, second.status);
    const second_view = second.resultView(problem).?;
    try std.testing.expect(second_view.basis_valid);
    try std.testing.expectEqual(@as(usize, 0), second_view.iterations);
    try std.testing.expectEqualSlices(f64, first_view.primal, second_view.primal);

    var session = LpSolveSession.init(std.testing.allocator);
    defer session.deinit();
    const revision = LpRevisionStamp{ .structure = 1, .matrix_values = 1, .bounds = 1, .objective = 1 };
    try std.testing.expectEqual(LpSolveStatus.optimal, session.solve(problem, .{}, revision));
    const warm_status = session.solve(problem, .{}, revision);
    try std.testing.expectEqual(LpSolveStatus.optimal, warm_status);
    try std.testing.expectEqual(@as(usize, 1), session.cold_solves);
    try std.testing.expectEqual(@as(usize, 1), session.reoptimizations);
    const warm_view = session.resultView(problem, warm_status).?;
    try std.testing.expect(warm_view.basis_valid);
    try std.testing.expectEqual(@as(usize, 0), warm_view.iterations);
    try std.testing.expectEqualSlices(f64, first_view.primal, warm_view.primal);
}

test "activity-redundant rows preserve unbounded and infeasible certificates" {
    const unbounded_problem = LpProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ -1.0, 0.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ std.math.inf(f64), std.math.inf(f64) },
        .row_lower = &.{0.0},
        .row_upper = &.{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(
            1,
            2,
            &[_]usize{ 0, 1, 2 },
            &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) },
            &.{ 1.0, 1.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var unbounded = solveLp(std.testing.allocator, unbounded_problem, .{});
    defer unbounded.deinit();
    try std.testing.expectEqual(LpSolveStatus.unbounded, unbounded.status);
    try std.testing.expectEqual(@as(u32, 1), unbounded.presolved_problem.?.redundant_row_count);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 0.0 }, unbounded.resultView(unbounded_problem).?.unbounded_ray);

    const infeasible_problem = LpProblemView{
        .num_rows = 2,
        .num_cols = 2,
        .col_cost = &.{ 0.0, 0.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ 1.0, 1.0 },
        .row_lower = &.{ 0.0, 3.0 },
        .row_upper = &.{ 2.0, std.math.inf(f64) },
        .matrix = matrix.CscView.initAssumeValid(
            2,
            2,
            &[_]usize{ 0, 2, 4 },
            &.{
                foundation.RowId.fromUsizeAssumeValid(0),
                foundation.RowId.fromUsizeAssumeValid(1),
                foundation.RowId.fromUsizeAssumeValid(0),
                foundation.RowId.fromUsizeAssumeValid(1),
            },
            &.{ 1.0, 1.0, 1.0, 1.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var infeasible = solveLp(std.testing.allocator, infeasible_problem, .{ .phase_one_strategy = .dual });
    defer infeasible.deinit();
    try std.testing.expectEqual(LpSolveStatus.infeasible, infeasible.status);
    try std.testing.expectEqual(@as(u32, 1), infeasible.presolved_problem.?.redundant_row_count);
    const infeasible_view = infeasible.resultView(infeasible_problem).?;
    try expectValidFarkas(infeasible_problem, infeasible_view.infeasibility_ray);
    try std.testing.expectEqual(@as(f64, 0.0), infeasible_view.infeasibility_ray[0]);
    try std.testing.expect(infeasible_view.infeasibility_ray[1] > 0.0);
}

test "exact duplicate rows restore dual basis and warm sessions" {
    const row0 = foundation.RowId.fromUsizeAssumeValid(0);
    const row1 = foundation.RowId.fromUsizeAssumeValid(1);
    const problem = LpProblemView{
        .num_rows = 2,
        .num_cols = 2,
        .col_cost = &.{ 1.0, 2.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ 10.0, 10.0 },
        .row_lower = &.{ 1.0, 1.0 },
        .row_upper = &.{ std.math.inf(f64), std.math.inf(f64) },
        .matrix = matrix.CscView.initAssumeValid(
            2,
            2,
            &[_]usize{ 0, 2, 4 },
            &.{ row0, row1, row0, row1 },
            &.{ 1.0, 1.0, 1.0, 1.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var first = solveLp(std.testing.allocator, problem, .{});
    defer first.deinit();
    try std.testing.expectEqual(LpSolveStatus.optimal, first.status);
    try std.testing.expectEqual(@as(u32, 1), first.presolved_problem.?.duplicate_row_count);
    try std.testing.expectEqual(@as(usize, 1), first.presolved_problem.?.passes.items.len);
    const first_view = first.resultView(problem).?;
    try std.testing.expect(first_view.basis_valid);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 0.0 }, first_view.primal);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 0.0 }, first_view.dual);
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 1.0 }, first_view.reduced_cost);
    try std.testing.expectEqualSlices(LpBasisStatus, &.{ .at_upper, .basic }, first_view.logical_status);
    try std.testing.expectEqual(@as(u32, 3), first_view.basic_index[1]);

    var snapshot = try LpBasisSnapshot.initFromView(std.testing.allocator, .{
        .structural_status = first_view.structural_status,
        .logical_status = first_view.logical_status,
        .basic_index = first_view.basic_index,
    });
    defer snapshot.deinit();
    var reduced_chain = try ReducedBasisChain.init(std.testing.allocator, &first.presolved_problem.?, snapshot.view());
    defer reduced_chain.deinit();
    try reduced_chain.view().validate(2, 1);
    var warm = solveLp(std.testing.allocator, problem, .{ .initial_basis = snapshot.view() });
    defer warm.deinit();
    try std.testing.expectEqual(LpSolveStatus.optimal, warm.status);
    try std.testing.expectEqual(@as(usize, 0), warm.resultView(problem).?.iterations);

    var session = LpSolveSession.init(std.testing.allocator);
    defer session.deinit();
    const revision = LpRevisionStamp{ .structure = 1, .matrix_values = 1, .bounds = 1, .objective = 1 };
    try std.testing.expectEqual(LpSolveStatus.optimal, session.solve(problem, .{}, revision));
    try std.testing.expectEqual(LpSolveStatus.optimal, session.solve(problem, .{}, revision));
    try std.testing.expectEqual(@as(usize, 1), session.cold_solves);
    try std.testing.expectEqual(@as(usize, 1), session.reoptimizations);
    try std.testing.expect(session.resultView(problem, .optimal).?.basis_valid);
}

test "duplicate row postsolve leaves the removed Farkas multiplier zero" {
    const row0 = foundation.RowId.fromUsizeAssumeValid(0);
    const row1 = foundation.RowId.fromUsizeAssumeValid(1);
    const row2 = foundation.RowId.fromUsizeAssumeValid(2);
    const problem = LpProblemView{
        .num_rows = 3,
        .num_cols = 2,
        .col_cost = &.{ 0.0, 0.0 },
        .col_lower = &.{ -std.math.inf(f64), -std.math.inf(f64) },
        .col_upper = &.{ std.math.inf(f64), std.math.inf(f64) },
        .row_lower = &.{ 0.0, 0.0, -std.math.inf(f64) },
        .row_upper = &.{ std.math.inf(f64), std.math.inf(f64), -1.0 },
        .matrix = matrix.CscView.initAssumeValid(
            3,
            2,
            &[_]usize{ 0, 3, 6 },
            &.{ row0, row1, row2, row0, row1, row2 },
            &.{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var result = solveLp(std.testing.allocator, problem, .{ .phase_one_strategy = .dual });
    defer result.deinit();
    try std.testing.expectEqual(LpSolveStatus.infeasible, result.status);
    try std.testing.expectEqual(@as(u32, 1), result.presolved_problem.?.duplicate_row_count);
    const ray = result.resultView(problem).?.infeasibility_ray;
    try expectValidFarkas(problem, ray);
    try std.testing.expect(ray[0] > 0.0);
    try std.testing.expectEqual(@as(f64, 0.0), ray[1]);
    try std.testing.expect(ray[2] < 0.0);
}

test "fixed substitution can expose duplicate adjusted rows in the same pass" {
    const row0 = foundation.RowId.fromUsizeAssumeValid(0);
    const row1 = foundation.RowId.fromUsizeAssumeValid(1);
    const problem = LpProblemView{
        .num_rows = 2,
        .num_cols = 3,
        .col_cost = &.{ 1.0, 2.0, 0.0 },
        .col_lower = &.{ 0.0, 0.0, 1.0 },
        .col_upper = &.{ 10.0, 10.0, 1.0 },
        .row_lower = &.{ 2.0, 3.0 },
        .row_upper = &.{ 5.0, 6.0 },
        .matrix = matrix.CscView.initAssumeValid(
            2,
            3,
            &[_]usize{ 0, 2, 4, 6 },
            &.{ row0, row1, row0, row1, row0, row1 },
            &.{ 1.0, 1.0, 1.0, 1.0, 1.0, 2.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var result = solveLp(std.testing.allocator, problem, .{});
    defer result.deinit();
    try std.testing.expectEqual(LpSolveStatus.optimal, result.status);
    try std.testing.expectEqual(@as(u32, 1), result.presolved_problem.?.column_substitution_count);
    try std.testing.expectEqual(@as(u32, 1), result.presolved_problem.?.duplicate_row_count);
    const view = result.resultView(problem).?;
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 0.0, 1.0 }, view.primal);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 0.0 }, view.dual);
    try std.testing.expect(view.basis_valid);
}

test "presolve substitutes empty columns at finite objective-optimal bounds" {
    const problem = LpProblemView{
        .num_rows = 1,
        .num_cols = 3,
        .col_cost = &.{ 1.0, 2.0, -4.0 },
        .col_lower = &.{ 0.0, -3.0, 0.0 },
        .col_upper = &.{ std.math.inf(f64), 7.0, 5.0 },
        .row_lower = &.{1.0},
        .row_upper = &.{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(
            1,
            3,
            &[_]usize{ 0, 1, 1, 1 },
            &.{foundation.RowId.fromUsizeAssumeValid(0)},
            &.{1.0},
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var first = solveLp(std.testing.allocator, problem, .{});
    defer first.deinit();
    try std.testing.expectEqual(LpSolveStatus.optimal, first.status);
    try std.testing.expectEqual(@as(usize, 2), first.presolved_problem.?.passes.items.len);
    try std.testing.expectEqual(@as(u32, 3), first.presolved_problem.?.column_substitution_count);
    try std.testing.expectEqual(@as(u32, 1), first.presolved_problem.?.singleton_row_count);
    try std.testing.expectEqual(@as(usize, 0), first.presolved_problem.?.reducedView().num_cols);
    try std.testing.expectEqual(@as(usize, 0), first.presolved_problem.?.reducedView().num_rows);
    const first_view = first.resultView(problem).?;
    try std.testing.expect(first_view.basis_valid);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, -3.0, 5.0 }, first_view.primal);
    try std.testing.expectEqualSlices(f64, &.{1.0}, first_view.dual);
    try std.testing.expectApproxEqAbs(@as(f64, -25.0), first_view.objective_value, 1e-12);
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 2.0, -4.0 }, first_view.reduced_cost);
    try std.testing.expectEqual(LpBasisStatus.basic, first_view.structural_status[0]);
    try std.testing.expectEqual(LpBasisStatus.at_lower, first_view.structural_status[1]);
    try std.testing.expectEqual(LpBasisStatus.at_upper, first_view.structural_status[2]);
    try std.testing.expectEqual(LpBasisStatus.at_upper, first_view.logical_status[0]);
    try std.testing.expectEqual(@as(u32, 0), first_view.basic_index[0]);

    var snapshot = try LpBasisSnapshot.initFromView(std.testing.allocator, .{
        .structural_status = first_view.structural_status,
        .logical_status = first_view.logical_status,
        .basic_index = first_view.basic_index,
    });
    defer snapshot.deinit();
    var second = solveLp(std.testing.allocator, problem, .{ .initial_basis = snapshot.view() });
    defer second.deinit();
    try std.testing.expectEqual(LpSolveStatus.optimal, second.status);
    const second_view = second.resultView(problem).?;
    try std.testing.expect(second_view.basis_valid);
    try std.testing.expectEqual(@as(usize, 0), second_view.iterations);
    try std.testing.expectEqualSlices(f64, first_view.primal, second_view.primal);
}

test "zero-cost empty columns publish canonical values basis and warm sessions" {
    const problem = LpProblemView{
        .num_rows = 0,
        .num_cols = 4,
        .col_cost = &.{ 0.0, 0.0, -0.0, 0.0 },
        .col_lower = &.{ -std.math.inf(f64), -std.math.inf(f64), 2.0, -2.0 },
        .col_upper = &.{ std.math.inf(f64), 3.0, std.math.inf(f64), 4.0 },
        .row_lower = &.{},
        .row_upper = &.{},
        .matrix = matrix.CscView.initAssumeValid(0, 4, &[_]usize{ 0, 0, 0, 0, 0 }, &.{}, &.{}),
        .objective_sense = .maximize,
        .objective_offset = 7.0,
    };
    var first = solveLp(std.testing.allocator, problem, .{});
    defer first.deinit();
    try std.testing.expectEqual(LpSolveStatus.optimal, first.status);
    try std.testing.expectEqual(@as(usize, 1), first.presolved_problem.?.passes.items.len);
    try std.testing.expectEqual(@as(u32, 4), first.presolved_problem.?.column_substitution_count);
    const first_view = first.resultView(problem).?;
    try std.testing.expect(first_view.basis_valid);
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 3.0, 2.0, -2.0 }, first_view.primal);
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 0.0, 0.0, 0.0 }, first_view.reduced_cost);
    try std.testing.expectEqualSlices(LpBasisStatus, &.{ .free, .at_upper, .at_lower, .at_lower }, first_view.structural_status);
    try std.testing.expectApproxEqAbs(@as(f64, 7.0), first_view.objective_value, 1e-12);

    var snapshot = try LpBasisSnapshot.initFromView(std.testing.allocator, .{
        .structural_status = first_view.structural_status,
        .logical_status = first_view.logical_status,
        .basic_index = first_view.basic_index,
    });
    defer snapshot.deinit();
    var warm = solveLp(std.testing.allocator, problem, .{ .initial_basis = snapshot.view() });
    defer warm.deinit();
    try std.testing.expectEqual(LpSolveStatus.optimal, warm.status);
    const warm_view = warm.resultView(problem).?;
    try std.testing.expect(warm_view.basis_valid);
    try std.testing.expectEqual(@as(usize, 0), warm_view.iterations);
    try std.testing.expectEqualSlices(f64, first_view.primal, warm_view.primal);

    var lower = [_]f64{ -std.math.inf(f64), 2.0 };
    var upper = [_]f64{ std.math.inf(f64), 5.0 };
    const session_problem = LpProblemView{
        .num_rows = 0,
        .num_cols = 2,
        .col_cost = &.{ 0.0, 0.0 },
        .col_lower = &lower,
        .col_upper = &upper,
        .row_lower = &.{},
        .row_upper = &.{},
        .matrix = matrix.CscView.initAssumeValid(0, 2, &[_]usize{ 0, 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var session = LpSolveSession.init(std.testing.allocator);
    defer session.deinit();
    const first_revision = LpRevisionStamp{ .structure = 1, .matrix_values = 1, .bounds = 1, .objective = 1 };
    try std.testing.expectEqual(LpSolveStatus.optimal, session.solve(session_problem, .{}, first_revision));
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 2.0 }, session.resultView(session_problem, .optimal).?.primal);
    lower[1] = 3.0;
    upper[0] = 4.0;
    const second_revision = LpRevisionStamp{ .structure = 1, .matrix_values = 1, .bounds = 2, .objective = 1 };
    try std.testing.expectEqual(LpSolveStatus.optimal, session.solve(session_problem, .{}, second_revision));
    const session_view = session.resultView(session_problem, .optimal).?;
    try std.testing.expectEqualSlices(f64, &.{ 4.0, 3.0 }, session_view.primal);
    try std.testing.expectEqualSlices(LpBasisStatus, &.{ .at_upper, .at_lower }, session_view.structural_status);
    try std.testing.expectEqual(@as(usize, 1), session.cold_solves);
    try std.testing.expectEqual(@as(usize, 1), session.reoptimizations);
}

test "zero-cost free-column substitution preserves a Farkas certificate" {
    const problem = LpProblemView{
        .num_rows = 2,
        .num_cols = 2,
        .col_cost = &.{ 0.0, 0.0 },
        .col_lower = &.{ 0.0, -std.math.inf(f64) },
        .col_upper = &.{ std.math.inf(f64), std.math.inf(f64) },
        .row_lower = &.{ 2.0, -std.math.inf(f64) },
        .row_upper = &.{ std.math.inf(f64), 1.0 },
        .matrix = matrix.CscView.initAssumeValid(
            2,
            2,
            &[_]usize{ 0, 2, 2 },
            &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1) },
            &.{ 1.0, 1.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var result = solveLp(std.testing.allocator, problem, .{ .phase_one_strategy = .dual });
    defer result.deinit();
    try std.testing.expectEqual(LpSolveStatus.infeasible, result.status);
    const view = result.resultView(problem).?;
    try expectValidFarkas(problem, view.infeasibility_ray);
    try std.testing.expectEqual(LpBasisStatus.free, view.structural_status[1]);
    try std.testing.expectEqual(@as(f64, 0.0), view.primal[1]);
}

test "free zero-cost singleton column restores primal basis and warm sessions" {
    const problem = LpProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 0.0, 1.0 },
        .col_lower = &.{ -std.math.inf(f64), 0.0 },
        .col_upper = &.{ std.math.inf(f64), 10.0 },
        .row_lower = &.{4.0},
        .row_upper = &.{8.0},
        .matrix = matrix.CscView.initAssumeValid(
            1,
            2,
            &[_]usize{ 0, 1, 2 },
            &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) },
            &.{ 2.0, 3.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var first = solveLp(std.testing.allocator, problem, .{});
    defer first.deinit();
    try std.testing.expectEqual(LpSolveStatus.optimal, first.status);
    try std.testing.expectEqual(@as(u32, 1), first.presolved_problem.?.free_singleton_column_count);
    const first_view = first.resultView(problem).?;
    try std.testing.expect(first_view.basis_valid);
    try std.testing.expectEqualSlices(f64, &.{ 2.0, 0.0 }, first_view.primal);
    try std.testing.expectEqualSlices(f64, &.{0.0}, first_view.dual);
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 1.0 }, first_view.reduced_cost);
    try std.testing.expectEqualSlices(LpBasisStatus, &.{ .basic, .at_lower }, first_view.structural_status);
    try std.testing.expectEqual(LpBasisStatus.at_upper, first_view.logical_status[0]);
    try std.testing.expectEqual(@as(u32, 0), first_view.basic_index[0]);

    var snapshot = try LpBasisSnapshot.initFromView(std.testing.allocator, .{
        .structural_status = first_view.structural_status,
        .logical_status = first_view.logical_status,
        .basic_index = first_view.basic_index,
    });
    defer snapshot.deinit();
    var reduced_chain = try ReducedBasisChain.init(std.testing.allocator, &first.presolved_problem.?, snapshot.view());
    defer reduced_chain.deinit();
    try reduced_chain.view().validate(0, 0);
    var warm = solveLp(std.testing.allocator, problem, .{ .initial_basis = snapshot.view() });
    defer warm.deinit();
    try std.testing.expectEqual(LpSolveStatus.optimal, warm.status);
    try std.testing.expect(warm.resultView(problem).?.basis_valid);
    try std.testing.expectEqualSlices(f64, first_view.primal, warm.resultView(problem).?.primal);

    var session = LpSolveSession.init(std.testing.allocator);
    defer session.deinit();
    const revision = LpRevisionStamp{ .structure = 1, .matrix_values = 1, .bounds = 1, .objective = 1 };
    try std.testing.expectEqual(LpSolveStatus.optimal, session.solve(problem, .{}, revision));
    try std.testing.expectEqual(LpSolveStatus.optimal, session.solve(problem, .{}, revision));
    try std.testing.expectEqual(@as(usize, 1), session.cold_solves);
    try std.testing.expectEqual(@as(usize, 1), session.reoptimizations);
    try std.testing.expect(session.resultView(problem, .optimal).?.basis_valid);
}

test "free singleton postsolve cancels an unbounded row direction" {
    const problem = LpProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 0.0, -1.0 },
        .col_lower = &.{ -std.math.inf(f64), 0.0 },
        .col_upper = &.{ std.math.inf(f64), std.math.inf(f64) },
        .row_lower = &.{5.0},
        .row_upper = &.{5.0},
        .matrix = matrix.CscView.initAssumeValid(
            1,
            2,
            &[_]usize{ 0, 1, 2 },
            &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) },
            &.{ 2.0, 3.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var result = solveLp(std.testing.allocator, problem, .{});
    defer result.deinit();
    try std.testing.expectEqual(LpSolveStatus.unbounded, result.status);
    const ray = result.resultView(problem).?.unbounded_ray;
    try std.testing.expect(result.resultView(problem).?.basis_valid);
    try std.testing.expectEqual(LpBasisStatus.basic, result.resultView(problem).?.structural_status[0]);
    try std.testing.expectEqual(LpBasisStatus.fixed, result.resultView(problem).?.logical_status[0]);
    try std.testing.expectEqualSlices(f64, &.{ -1.5, 1.0 }, ray);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), 2.0 * ray[0] + 3.0 * ray[1], 1e-12);
    try std.testing.expect(problem.col_cost[0] * ray[0] + problem.col_cost[1] * ray[1] < 0.0);
}

test "free singleton postsolve leaves its Farkas row multiplier zero" {
    const problem = LpProblemView{
        .num_rows = 2,
        .num_cols = 2,
        .col_cost = &.{ 0.0, 0.0 },
        .col_lower = &.{ -std.math.inf(f64), 0.0 },
        .col_upper = &.{ std.math.inf(f64), 1.0 },
        .row_lower = &.{ 0.0, 2.0 },
        .row_upper = &.{ 0.0, std.math.inf(f64) },
        .matrix = matrix.CscView.initAssumeValid(
            2,
            2,
            &[_]usize{ 0, 1, 3 },
            &.{
                foundation.RowId.fromUsizeAssumeValid(0),
                foundation.RowId.fromUsizeAssumeValid(0),
                foundation.RowId.fromUsizeAssumeValid(1),
            },
            &.{ 1.0, 1.0, 1.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var result = solveLp(std.testing.allocator, problem, .{ .phase_one_strategy = .dual });
    defer result.deinit();
    try std.testing.expectEqual(LpSolveStatus.infeasible, result.status);
    const ray = result.resultView(problem).?.infeasibility_ray;
    try expectValidFarkas(problem, ray);
    try std.testing.expectEqual(@as(f64, 0.0), ray[0]);
    try std.testing.expect(ray[1] > 0.0);
}

test "multiple free singleton columns on one row converge across passes" {
    const row = foundation.RowId.fromUsizeAssumeValid(0);
    const problem = LpProblemView{
        .num_rows = 1,
        .num_cols = 3,
        .col_cost = &.{ 0.0, 0.0, 0.0 },
        .col_lower = &.{ -std.math.inf(f64), -std.math.inf(f64), 2.0 },
        .col_upper = &.{ std.math.inf(f64), std.math.inf(f64), 4.0 },
        .row_lower = &.{7.0},
        .row_upper = &.{9.0},
        .matrix = matrix.CscView.initAssumeValid(
            1,
            3,
            &[_]usize{ 0, 1, 2, 3 },
            &.{ row, row, row },
            &.{ 1.0, 2.0, 3.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var result = solveLp(std.testing.allocator, problem, .{});
    defer result.deinit();
    try std.testing.expectEqual(LpSolveStatus.optimal, result.status);
    try std.testing.expectEqual(@as(usize, 2), result.presolved_problem.?.passes.items.len);
    try std.testing.expectEqual(@as(u32, 1), result.presolved_problem.?.free_singleton_column_count);
    try std.testing.expectEqual(@as(u32, 2), result.presolved_problem.?.column_substitution_count);
    const view = result.resultView(problem).?;
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 0.0, 2.0 }, view.primal);
    try std.testing.expectEqualSlices(LpBasisStatus, &.{ .basic, .free, .at_lower }, view.structural_status);
    try std.testing.expectEqual(LpBasisStatus.at_upper, view.logical_status[0]);
    try std.testing.expect(view.basis_valid);
}

test "maximize empty-column substitution publishes normalized reduced cost" {
    const problem = LpProblemView{
        .num_rows = 0,
        .num_cols = 1,
        .col_cost = &.{3.0},
        .col_lower = &.{-2.0},
        .col_upper = &.{5.0},
        .row_lower = &.{},
        .row_upper = &.{},
        .matrix = matrix.CscView.initAssumeValid(0, 1, &[_]usize{ 0, 0 }, &.{}, &.{}),
        .objective_sense = .maximize,
        .objective_offset = 0.0,
    };
    var result = solveLp(std.testing.allocator, problem, .{});
    defer result.deinit();
    try std.testing.expectEqual(LpSolveStatus.optimal, result.status);
    const view = result.resultView(problem).?;
    try std.testing.expect(view.primal_valid);
    try std.testing.expect(view.dual_valid);
    try std.testing.expect(view.basis_valid);
    try std.testing.expectEqualSlices(f64, &.{5.0}, view.primal);
    try std.testing.expectEqualSlices(f64, &.{-3.0}, view.reduced_cost);
    try std.testing.expectEqual(LpBasisStatus.at_upper, view.structural_status[0]);
    try std.testing.expectApproxEqAbs(@as(f64, 15.0), view.objective_value, 1e-12);
}

test "singleton upper bound postsolves dual reduced cost and basis" {
    const problem = LpProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &.{-1.0},
        .col_lower = &.{0.0},
        .col_upper = &.{10.0},
        .row_lower = &.{-std.math.inf(f64)},
        .row_upper = &.{4.0},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &[_]usize{ 0, 1 }, &.{foundation.RowId.fromUsizeAssumeValid(0)}, &.{1.0}),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var result = solveLp(std.testing.allocator, problem, .{});
    defer result.deinit();
    try std.testing.expectEqual(LpSolveStatus.optimal, result.status);
    const view = result.resultView(problem).?;
    try std.testing.expect(view.primal_valid);
    try std.testing.expect(view.dual_valid);
    try std.testing.expect(view.basis_valid);
    try std.testing.expectEqualSlices(f64, &.{4.0}, view.primal);
    try std.testing.expectEqualSlices(f64, &.{-1.0}, view.dual);
    try std.testing.expectEqualSlices(f64, &.{0.0}, view.reduced_cost);
    try std.testing.expectEqual(LpBasisStatus.basic, view.structural_status[0]);
    try std.testing.expectEqual(LpBasisStatus.at_lower, view.logical_status[0]);
    try std.testing.expectEqual(@as(u32, 0), view.basic_index[0]);
}

test "singleton postsolve restores Farkas rays for every coefficient and bound orientation" {
    const cases = [_]struct {
        coefficients: [2]f64,
        row_lower: [2]f64,
        row_upper: [2]f64,
        restored_sign: f64,
    }{
        .{
            .coefficients = .{ 1.0, 1.0 },
            .row_lower = .{ 2.0, -std.math.inf(f64) },
            .row_upper = .{ std.math.inf(f64), 1.0 },
            .restored_sign = 1.0,
        },
        .{
            .coefficients = .{ 1.0, 1.0 },
            .row_lower = .{ -std.math.inf(f64), 2.0 },
            .row_upper = .{ 1.0, std.math.inf(f64) },
            .restored_sign = -1.0,
        },
        .{
            .coefficients = .{ -1.0, 1.0 },
            .row_lower = .{ -std.math.inf(f64), -std.math.inf(f64) },
            .row_upper = .{ -2.0, 1.0 },
            .restored_sign = -1.0,
        },
        .{
            .coefficients = .{ -1.0, 1.0 },
            .row_lower = .{ -1.0, 2.0 },
            .row_upper = .{ std.math.inf(f64), std.math.inf(f64) },
            .restored_sign = 1.0,
        },
    };
    for (cases) |case| {
        const problem = LpProblemView{
            .num_rows = 2,
            .num_cols = 1,
            .col_cost = &.{0.0},
            .col_lower = &.{0.0},
            .col_upper = &.{std.math.inf(f64)},
            .row_lower = &case.row_lower,
            .row_upper = &case.row_upper,
            .matrix = matrix.CscView.initAssumeValid(
                2,
                1,
                &[_]usize{ 0, 2 },
                &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1) },
                &case.coefficients,
            ),
            .objective_sense = .minimize,
            .objective_offset = 0.0,
        };
        var result = solveLp(std.testing.allocator, problem, .{ .phase_one_strategy = .dual });
        defer result.deinit();
        try std.testing.expectEqual(LpSolveStatus.infeasible, result.status);
        const view = result.resultView(problem).?;
        try std.testing.expectEqual(@as(u32, 1), result.presolved_problem.?.singleton_row_count);
        try expectValidFarkas(problem, view.infeasibility_ray);
        try std.testing.expect(view.infeasibility_ray[0] * case.restored_sign > 0.0);
        try std.testing.expectApproxEqAbs(
            @as(f64, 0.0),
            view.infeasibility_ray[0] * case.coefficients[0] + view.infeasibility_ray[1] * case.coefficients[1],
            1e-9,
        );
    }
}

test "singleton Farkas postsolve selects the last active bound source" {
    const problem = LpProblemView{
        .num_rows = 3,
        .num_cols = 1,
        .col_cost = &.{0.0},
        .col_lower = &.{0.0},
        .col_upper = &.{std.math.inf(f64)},
        .row_lower = &.{ 1.0, 2.0, -std.math.inf(f64) },
        .row_upper = &.{ std.math.inf(f64), std.math.inf(f64), 0.0 },
        .matrix = matrix.CscView.initAssumeValid(
            3,
            1,
            &[_]usize{ 0, 3 },
            &.{
                foundation.RowId.fromUsizeAssumeValid(0),
                foundation.RowId.fromUsizeAssumeValid(1),
                foundation.RowId.fromUsizeAssumeValid(2),
            },
            &.{ 1.0, 1.0, 1.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var result = solveLp(std.testing.allocator, problem, .{ .phase_one_strategy = .dual });
    defer result.deinit();
    try std.testing.expectEqual(LpSolveStatus.infeasible, result.status);
    try std.testing.expectEqual(@as(u32, 2), result.presolved_problem.?.singleton_row_count);
    const view = result.resultView(problem).?;
    try expectValidFarkas(problem, view.infeasibility_ray);
    try std.testing.expectEqual(@as(f64, 0.0), view.infeasibility_ray[0]);
    try std.testing.expect(view.infeasibility_ray[1] > 0.0);
    try std.testing.expect(view.infeasibility_ray[2] < 0.0);
}

test "multi-pass postsolve preserves a singleton Farkas certificate" {
    const problem = LpProblemView{
        .num_rows = 3,
        .num_cols = 2,
        .col_cost = &.{ 0.0, 1.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ std.math.inf(f64), 4.0 },
        .row_lower = &.{ 2.0, -std.math.inf(f64), -std.math.inf(f64) },
        .row_upper = &.{ std.math.inf(f64), 1.0, 3.0 },
        .matrix = matrix.CscView.initAssumeValid(
            3,
            2,
            &[_]usize{ 0, 2, 3 },
            &.{
                foundation.RowId.fromUsizeAssumeValid(0),
                foundation.RowId.fromUsizeAssumeValid(1),
                foundation.RowId.fromUsizeAssumeValid(2),
            },
            &.{ 1.0, 1.0, 1.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var result = solveLp(std.testing.allocator, problem, .{ .phase_one_strategy = .dual });
    defer result.deinit();
    try std.testing.expectEqual(LpSolveStatus.infeasible, result.status);
    try std.testing.expectEqual(@as(usize, 2), result.presolved_problem.?.passes.items.len);
    try std.testing.expectEqual(@as(u32, 2), result.presolved_problem.?.singleton_row_count);
    const view = result.resultView(problem).?;
    try expectValidFarkas(problem, view.infeasibility_ray);
    try std.testing.expect(view.infeasibility_ray[0] > 0.0);
    try std.testing.expect(view.infeasibility_ray[1] < 0.0);
    try std.testing.expectEqual(@as(f64, 0.0), view.infeasibility_ray[2]);
}

test "reduced basis mapping rejects a basic substituted column" {
    const problem = LpProblemView{
        .num_rows = 1,
        .num_cols = 3,
        .col_cost = &.{ 1.0, 1.0, 0.0 },
        .col_lower = &.{ 0.0, 0.0, 2.0 },
        .col_upper = &.{ std.math.inf(f64), std.math.inf(f64), 2.0 },
        .row_lower = &.{-std.math.inf(f64)},
        .row_upper = &.{10.0},
        .matrix = matrix.CscView.initAssumeValid(
            1,
            3,
            &[_]usize{ 0, 1, 2, 2 },
            &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) },
            &.{ 1.0, 1.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var pp = try presolve.presolve(std.testing.allocator, problem);
    defer pp.deinit();
    const original = LpBasisView{
        .structural_status = &.{ .at_lower, .at_lower, .basic },
        .logical_status = &.{.at_lower},
        .basic_index = &.{2},
    };
    try std.testing.expectError(
        error.SubstitutedBasicColumn,
        ReducedBasisProjection.init(std.testing.allocator, &pp, original),
    );
}

test "multi-pass basis chain projects a logical warm basis to the fixed point" {
    const problem = LpProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 1.0, 0.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ 1.0, 0.0 },
        .row_lower = &.{-std.math.inf(f64)},
        .row_upper = &.{2.0},
        .matrix = matrix.CscView.initAssumeValid(
            1,
            2,
            &[_]usize{ 0, 1, 2 },
            &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) },
            &.{ 1.0, 1.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var pipeline = try presolve.presolvePipeline(std.testing.allocator, problem);
    defer pipeline.deinit();
    try std.testing.expectEqual(@as(usize, 2), pipeline.passes.items.len);
    const original = LpBasisView{
        .structural_status = &.{ .at_lower, .fixed },
        .logical_status = &.{.basic},
        .basic_index = &.{2},
    };
    var chain = try ReducedBasisChain.init(std.testing.allocator, &pipeline, original);
    defer chain.deinit();
    try std.testing.expectEqual(@as(usize, 0), chain.view().structural_status.len);
    try std.testing.expectEqual(@as(usize, 0), chain.view().logical_status.len);
    try std.testing.expectEqual(@as(usize, 0), chain.view().basic_index.len);
}

test "postsolve preserves unbounded certificates and validity semantics" {
    const problem = LpProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ -1.0, 3.0 },
        .col_lower = &.{ 0.0, 2.0 },
        .col_upper = &.{ std.math.inf(f64), 2.0 },
        .row_lower = &.{-std.math.inf(f64)},
        .row_upper = &.{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &[_]usize{ 0, 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var result = solveLp(std.testing.allocator, problem, .{});
    defer result.deinit();
    try std.testing.expectEqual(LpSolveStatus.unbounded, result.status);
    const view = result.resultView(problem).?;
    try std.testing.expect(!view.primal_valid);
    try std.testing.expect(!view.dual_valid);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 0.0 }, view.unbounded_ray);
    try std.testing.expectEqual(@as(usize, 0), view.infeasibility_ray.len);
}

test "postsolve preserves infeasibility certificates independently of basis validity" {
    const problem = LpProblemView{
        .num_rows = 2,
        .num_cols = 3,
        .col_cost = &.{ 0.0, 0.0, 5.0 },
        .col_lower = &.{ 0.0, 0.0, 0.0 },
        .col_upper = &.{ std.math.inf(f64), std.math.inf(f64), 0.0 },
        .row_lower = &.{ 2.0, -std.math.inf(f64) },
        .row_upper = &.{ std.math.inf(f64), 1.0 },
        .matrix = matrix.CscView.initAssumeValid(
            2,
            3,
            &[_]usize{ 0, 2, 4, 4 },
            &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1) },
            &.{ 1.0, 1.0, 1.0, 1.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var result = solveLp(std.testing.allocator, problem, .{ .phase_one_strategy = .dual });
    defer result.deinit();
    try std.testing.expectEqual(LpSolveStatus.infeasible, result.status);
    const view = result.resultView(problem).?;
    try std.testing.expect(!view.primal_valid);
    try std.testing.expect(!view.dual_valid);
    try std.testing.expectEqual(@as(usize, 0), view.unbounded_ray.len);
    try std.testing.expectEqual(@as(usize, 2), view.infeasibility_ray.len);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), view.infeasibility_ray[0] + view.infeasibility_ray[1], 1e-9);
}

test "presolve handles an all-fixed feasible model" {
    const problem = LpProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &.{4.0},
        .col_lower = &.{2.0},
        .col_upper = &.{2.0},
        .row_lower = &.{2.0},
        .row_upper = &.{2.0},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &[_]usize{ 0, 1 }, &.{foundation.RowId.fromUsizeAssumeValid(0)}, &.{1.0}),
        .objective_sense = .maximize,
        .objective_offset = 3.0,
    };
    var result = solveLp(std.testing.allocator, problem, .{});
    defer result.deinit();
    try std.testing.expectEqual(LpSolveStatus.optimal, result.status);
    const view = result.resultView(problem).?;
    try std.testing.expectEqualSlices(f64, &.{2.0}, view.primal);
    try std.testing.expectApproxEqAbs(@as(f64, 11.0), view.objective_value, 1e-12);
    try std.testing.expect(view.basis_valid);
}

test "presolve handles an all-fixed infeasible model" {
    const problem = LpProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &.{0.0},
        .col_lower = &.{2.0},
        .col_upper = &.{2.0},
        .row_lower = &.{3.0},
        .row_upper = &.{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &[_]usize{ 0, 1 }, &.{foundation.RowId.fromUsizeAssumeValid(0)}, &.{1.0}),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var result = solveLp(std.testing.allocator, problem, .{ .phase_one_strategy = .dual });
    defer result.deinit();
    try std.testing.expectEqual(LpSolveStatus.infeasible, result.status);
    const view = result.resultView(problem).?;
    try std.testing.expect(!view.primal_valid);
    try expectValidFarkas(problem, view.infeasibility_ray);
    try std.testing.expectEqualSlices(f64, &.{1.0}, view.infeasibility_ray);
}

test "simplex publishes both orientations of an original empty-row Farkas ray" {
    const cases = [_]struct { lower: f64, upper: f64, multiplier: f64 }{
        .{ .lower = 2.0, .upper = std.math.inf(f64), .multiplier = 1.0 },
        .{ .lower = -std.math.inf(f64), .upper = -3.0, .multiplier = -1.0 },
    };
    for (cases) |case| {
        const problem = LpProblemView{
            .num_rows = 1,
            .num_cols = 0,
            .col_cost = &.{},
            .col_lower = &.{},
            .col_upper = &.{},
            .row_lower = &.{case.lower},
            .row_upper = &.{case.upper},
            .matrix = matrix.CscView.initAssumeValid(1, 0, &[_]usize{0}, &.{}, &.{}),
            .objective_sense = .minimize,
            .objective_offset = 0.0,
        };
        var result = solveLp(std.testing.allocator, problem, .{ .phase_one_strategy = .dual });
        defer result.deinit();
        try std.testing.expectEqual(LpSolveStatus.infeasible, result.status);
        try std.testing.expect(result.presolved_problem == null);
        const ray = result.resultView(problem).?.infeasibility_ray;
        try std.testing.expectEqualSlices(f64, &.{case.multiplier}, ray);
        try expectValidFarkas(problem, ray);
    }
}

test "empty-row Farkas terminal is refreshed during direct reoptimization" {
    var row_lower = [_]f64{-1.0};
    var row_upper = [_]f64{1.0};
    const problem = LpProblemView{
        .num_rows = 1,
        .num_cols = 0,
        .col_cost = &.{},
        .col_lower = &.{},
        .col_upper = &.{},
        .row_lower = &row_lower,
        .row_upper = &row_upper,
        .matrix = matrix.CscView.initAssumeValid(1, 0, &[_]usize{0}, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var engine = lp.simplex.engine.SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(LpSolveStatus.optimal, engine.solveProblem(problem, .{}));
    row_lower[0] = -std.math.inf(f64);
    row_upper[0] = -2.0;
    const status = engine.reoptimizeProblem(problem, .{});
    try std.testing.expectEqual(LpSolveStatus.infeasible, status);
    const ray = engine.solutionView(problem, status).?.infeasibility_ray;
    try std.testing.expectEqualSlices(f64, &.{-1.0}, ray);
    try expectValidFarkas(problem, ray);
}

test "empty-row Farkas ray survives a two-pass presolve chain" {
    const problem = LpProblemView{
        .num_rows = 2,
        .num_cols = 2,
        .col_cost = &.{ 0.0, 1.0 },
        .col_lower = &.{ 2.0, 0.0 },
        .col_upper = &.{ 2.0, 1.0 },
        .row_lower = &.{ 3.0, 0.0 },
        .row_upper = &.{ std.math.inf(f64), 1.0 },
        .matrix = matrix.CscView.initAssumeValid(
            2,
            2,
            &[_]usize{ 0, 1, 2 },
            &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1) },
            &.{ 1.0, 1.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var result = solveLp(std.testing.allocator, problem, .{ .phase_one_strategy = .dual });
    defer result.deinit();
    try std.testing.expectEqual(LpSolveStatus.infeasible, result.status);
    try std.testing.expectEqual(@as(usize, 2), result.presolved_problem.?.passes.items.len);
    const ray = result.resultView(problem).?.infeasibility_ray;
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 0.0 }, ray);
    try expectValidFarkas(problem, ray);
}

test "session invalidates reuse when the fixed-column mapping changes" {
    var session = LpSolveSession.init(std.testing.allocator);
    defer session.deinit();
    const revision = LpRevisionStamp{ .structure = 1, .matrix_values = 1, .bounds = 1, .objective = 1 };

    const first = fixedColumnProblem(false);
    try std.testing.expectEqual(LpSolveStatus.optimal, session.solve(first, .{}, revision));
    try std.testing.expectEqual(@as(usize, 1), session.cold_solves);
    try std.testing.expectEqualSlices(f64, &.{ 2.0, 3.0 }, session.resultView(first, .optimal).?.primal);

    const second = fixedColumnProblem(true);
    const changed_revision = LpRevisionStamp{ .structure = 1, .matrix_values = 1, .bounds = 2, .objective = 1 };
    try std.testing.expectEqual(LpSolveStatus.optimal, session.solve(second, .{}, changed_revision));
    try std.testing.expectEqual(@as(usize, 2), session.cold_solves);
    try std.testing.expectEqualSlices(f64, &.{ 2.0, 3.0 }, session.resultView(second, .optimal).?.primal);
}

test "session reoptimizes when an empty column switches optimal bound" {
    var costs = [_]f64{2.0};
    const problem = LpProblemView{
        .num_rows = 0,
        .num_cols = 1,
        .col_cost = &costs,
        .col_lower = &.{-3.0},
        .col_upper = &.{7.0},
        .row_lower = &.{},
        .row_upper = &.{},
        .matrix = matrix.CscView.initAssumeValid(0, 1, &[_]usize{ 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var session = LpSolveSession.init(std.testing.allocator);
    defer session.deinit();
    const first_revision = LpRevisionStamp{ .structure = 1, .matrix_values = 1, .bounds = 1, .objective = 1 };
    try std.testing.expectEqual(LpSolveStatus.optimal, session.solve(problem, .{}, first_revision));
    const first = session.resultView(problem, .optimal).?;
    try std.testing.expectEqualSlices(f64, &.{-3.0}, first.primal);
    try std.testing.expectApproxEqAbs(@as(f64, -6.0), first.objective_value, 1e-12);
    try std.testing.expectEqual(LpBasisStatus.at_lower, first.structural_status[0]);
    try std.testing.expectEqual(@as(usize, 1), session.cold_solves);

    costs[0] = -2.0;
    const second_revision = LpRevisionStamp{ .structure = 1, .matrix_values = 1, .bounds = 1, .objective = 2 };
    try std.testing.expectEqual(LpSolveStatus.optimal, session.solve(problem, .{}, second_revision));
    const second = session.resultView(problem, .optimal).?;
    try std.testing.expectEqualSlices(f64, &.{7.0}, second.primal);
    try std.testing.expectApproxEqAbs(@as(f64, -14.0), second.objective_value, 1e-12);
    try std.testing.expectEqual(LpBasisStatus.at_upper, second.structural_status[0]);
    try std.testing.expectEqual(@as(usize, 1), session.cold_solves);
    try std.testing.expectEqual(@as(usize, 1), session.reoptimizations);
}

test "invalid models are rejected before presolve and do not publish results" {
    const invalid_matrix_problem = LpProblemView{
        .num_rows = 0,
        .num_cols = 1,
        .col_cost = &.{0.0},
        .col_lower = &.{0.0},
        .col_upper = &.{1.0},
        .row_lower = &.{},
        .row_upper = &.{},
        .matrix = matrix.CscView.initAssumeValid(0, 1, &[_]usize{ 1, 1 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var malformed = solveLp(std.testing.allocator, invalid_matrix_problem, .{ .time_limit_ns = 0 });
    defer malformed.deinit();
    try std.testing.expectEqual(LpSolveStatus.invalid_model, malformed.status);
    try std.testing.expect(malformed.presolved_problem == null);
    try std.testing.expect(malformed.resultView(invalid_matrix_problem) == null);
    try std.testing.expectEqual(@as(usize, 0), malformed.engine.iterations);

    const invalid_bounds_problem = LpProblemView{
        .num_rows = 0,
        .num_cols = 1,
        .col_cost = &.{0.0},
        .col_lower = &.{2.0},
        .col_upper = &.{1.0},
        .row_lower = &.{},
        .row_upper = &.{},
        .matrix = matrix.CscView.initAssumeValid(0, 1, &[_]usize{ 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var reversed = solveLp(std.testing.allocator, invalid_bounds_problem, .{});
    defer reversed.deinit();
    try std.testing.expectEqual(LpSolveStatus.invalid_model, reversed.status);
    try std.testing.expect(reversed.resultView(invalid_bounds_problem) == null);
}

test "session clears stale state after an invalid model and resumes cold" {
    var lower = [_]f64{0.0};
    const problem = LpProblemView{
        .num_rows = 0,
        .num_cols = 1,
        .col_cost = &.{0.0},
        .col_lower = &lower,
        .col_upper = &.{1.0},
        .row_lower = &.{},
        .row_upper = &.{},
        .matrix = matrix.CscView.initAssumeValid(0, 1, &[_]usize{ 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var session = LpSolveSession.init(std.testing.allocator);
    defer session.deinit();
    const revision = LpRevisionStamp{ .structure = 1, .matrix_values = 1, .bounds = 1, .objective = 1 };

    try std.testing.expectEqual(LpSolveStatus.optimal, session.solve(problem, .{}, revision));
    try std.testing.expect(session.presolved_problem != null);
    try std.testing.expectEqual(@as(usize, 1), session.cold_solves);

    lower[0] = std.math.nan(f64);
    const invalid_status = session.solve(problem, .{}, revision);
    try std.testing.expectEqual(LpSolveStatus.invalid_model, invalid_status);
    try std.testing.expect(session.presolved_problem == null);
    try std.testing.expect(session.last_revision == null);
    try std.testing.expect(session.resultView(problem, invalid_status) == null);

    lower[0] = 0.0;
    try std.testing.expectEqual(LpSolveStatus.optimal, session.solve(problem, .{}, revision));
    try std.testing.expectEqual(@as(usize, 2), session.cold_solves);
    try std.testing.expectEqual(@as(usize, 0), session.reoptimizations);
    try std.testing.expectEqualSlices(f64, &.{0.0}, session.resultView(problem, .optimal).?.primal);
}

test {
    std.testing.refAllDecls(@This());
}
