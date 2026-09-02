//! Reversible LP presolve operating at the ProblemView level.
//!
//! Reduces a problem by removing redundant structure, then postsolves
//! the simplex engine output back to original coordinates. Each reduction
//! records the transformation so postsolve is exact (not approximate).
//!
//! Phase 1: Exact column substitution.
//! Fixed columns and structurally empty columns with a finite objective-optimal
//! bound are substituted. Their contribution is moved into the row bounds and
//! objective offset before the column is deleted.
//!
//! Phase 2: Redundant empty-row elimination.
//! Rows with no retained coefficients and bounds containing zero are removed.
//! Original row coordinates are retained so duals, certificates, and bases can
//! be reconstructed exactly during postsolve.
//!
//! Phase 3: Singleton-row bound transfer.
//! A row containing one retained coefficient is converted into bounds on that
//! column and removed. Postsolve records the originating row side so dual and
//! basis coordinates can be restored.

const std = @import("std");
const foundation = @import("foundation");
const matrix = @import("matrix");
const lp = @import("lp");
const rules = @import("rules/root.zig");
const ProblemView = lp.ProblemView;
const SolveStatus = lp.SolveStatus;
const SolutionView = lp.SolutionView;

/// Record for one exactly substituted column.
pub const ColumnSubstitutionRecord = struct {
    original_column: u32,
    value: f64,
    col_cost: f64,
    basis_status: lp.simplex.basis.BasisStatus,
    row_start: u32,
    row_end: u32,
};

pub const SingletonRowRecord = struct {
    original_row: u32,
    original_column: u32,
    coefficient: f64,
    previous_lower: f64,
    previous_upper: f64,
    tightened_lower: f64,
    tightened_upper: f64,
    lower_source: ?rules.singleton_row.RowSide,
    upper_source: ?rules.singleton_row.RowSide,
};

/// Owns reduced problem data and transformation records for postsolve.
pub const PresolvedProblem = struct {
    allocator: std.mem.Allocator,

    /// Reduced dimensions.
    num_rows: u32,
    num_cols: u32,

    /// Original dimensions required to restore rectangular models.
    original_num_rows: u32,
    original_num_cols: u32,

    /// Reduced problem data (owning slices).
    col_cost: []f64,
    col_lower: []f64,
    col_upper: []f64,
    row_lower: []f64,
    row_upper: []f64,
    objective_offset: f64,
    objective_sense: lp.ObjectiveSense,

    /// Reduced constraint matrix (owning CSC).
    matrix_storage: matrix.CscMatrix = undefined,
    _matrix_valid: bool = false,

    // ── Index mapping ──
    reduced_to_original_col: []u32,
    original_to_reduced_col: []u32,
    reduced_to_original_row: []u32,
    original_to_reduced_row: []u32,

    // ── Exact column-substitution records ──
    substitution_row_indices: []u32,
    substitution_row_values: []f64,
    column_substitutions: []ColumnSubstitutionRecord,
    column_substitution_count: u32,
    singleton_rows: []SingletonRowRecord,
    singleton_row_count: u32,
    removed_row_count: u32,

    pub fn wasApplied(self: *const PresolvedProblem) bool {
        return self.column_substitution_count > 0 or self.removed_row_count > 0;
    }

    pub fn reducedView(self: *const PresolvedProblem) ProblemView {
        const csc = self.matrix_storage.view();
        return .{
            .num_rows = self.num_rows,
            .num_cols = self.num_cols,
            .col_cost = self.col_cost,
            .col_lower = self.col_lower,
            .col_upper = self.col_upper,
            .row_lower = self.row_lower,
            .row_upper = self.row_upper,
            .matrix = csc,
            .objective_sense = self.objective_sense,
            .objective_offset = self.objective_offset,
        };
    }

    pub fn deinit(self: *PresolvedProblem) void {
        if (self._matrix_valid) self.matrix_storage.deinit(self.allocator);
        self.allocator.free(self.col_cost);
        self.allocator.free(self.col_lower);
        self.allocator.free(self.col_upper);
        self.allocator.free(self.row_lower);
        self.allocator.free(self.row_upper);
        self.allocator.free(self.reduced_to_original_col);
        self.allocator.free(self.original_to_reduced_col);
        self.allocator.free(self.reduced_to_original_row);
        self.allocator.free(self.original_to_reduced_row);
        self.allocator.free(self.substitution_row_indices);
        self.allocator.free(self.substitution_row_values);
        self.allocator.free(self.column_substitutions);
        self.allocator.free(self.singleton_rows);
    }
};

/// Result of presolve → solve → postsolve round-trip.
pub const PresolveResult = struct {
    allocator: std.mem.Allocator,
    status: SolveStatus,
    objective_value: f64,
    primal: []f64,
    dual: []f64,
    reduced_cost: []f64,
    unbounded_ray: []f64,
    infeasibility_ray: []f64,

    pub fn deinit(self: *PresolveResult) void {
        self.allocator.free(self.primal);
        self.allocator.free(self.dual);
        self.allocator.free(self.reduced_cost);
        self.allocator.free(self.unbounded_ray);
        self.allocator.free(self.infeasibility_ray);
        self.* = undefined;
    }

    pub fn solutionView(self: *const PresolveResult, iterations: usize) SolutionView {
        return .{
            .status = self.status,
            .primal = self.primal,
            .dual = self.dual,
            .reduced_cost = self.reduced_cost,
            .unbounded_ray = self.unbounded_ray,
            .infeasibility_ray = self.infeasibility_ray,
            .objective_value = self.objective_value,
            .iterations = iterations,
        };
    }
};

/// Owns a converged sequence of relative-coordinate presolve passes.
/// Capacity is reserved before the first pass so every intermediate
/// `reducedView` remains valid while the next pass is constructed.
pub const PresolvePipeline = struct {
    allocator: std.mem.Allocator,
    passes: std.ArrayListUnmanaged(PresolvedProblem) = .empty,
    column_substitution_count: u32 = 0,
    singleton_row_count: u32 = 0,
    removed_row_count: u32 = 0,

    pub fn init(allocator: std.mem.Allocator, problem: ProblemView) !PresolvePipeline {
        var result = PresolvePipeline{ .allocator = allocator };
        errdefer result.deinit();
        const maximum_passes = problem.num_rows + problem.num_cols + 1;
        try result.passes.ensureTotalCapacity(allocator, maximum_passes);

        var current = problem;
        while (result.passes.items.len < maximum_passes) {
            var pass = try presolve(allocator, current);
            if (!pass.wasApplied()) {
                pass.deinit();
                break;
            }
            std.debug.assert(pass.num_rows < current.num_rows or pass.num_cols < current.num_cols);
            result.column_substitution_count += pass.column_substitution_count;
            result.singleton_row_count += pass.singleton_row_count;
            result.removed_row_count += pass.removed_row_count;
            result.passes.appendAssumeCapacity(pass);
            current = result.passes.items[result.passes.items.len - 1].reducedView();
        }
        return result;
    }

    pub fn deinit(self: *PresolvePipeline) void {
        for (self.passes.items) |*pass| pass.deinit();
        self.passes.deinit(self.allocator);
        self.* = undefined;
    }

    pub fn wasApplied(self: *const PresolvePipeline) bool {
        return self.passes.items.len != 0;
    }

    pub fn reducedView(self: *const PresolvePipeline) ProblemView {
        std.debug.assert(self.passes.items.len != 0);
        return self.passes.items[self.passes.items.len - 1].reducedView();
    }

    pub fn sameMapping(self: *const PresolvePipeline, other: *const PresolvePipeline) bool {
        if (self.passes.items.len != other.passes.items.len) return false;
        for (self.passes.items, other.passes.items) |left, right| {
            if (!std.mem.eql(u32, left.reduced_to_original_col, right.reduced_to_original_col) or
                !std.mem.eql(u32, left.reduced_to_original_row, right.reduced_to_original_row)) return false;
        }
        return true;
    }

    pub fn postsolveAll(
        self: *const PresolvePipeline,
        reduced_solution: SolutionView,
        allocator: std.mem.Allocator,
    ) !PresolveResult {
        std.debug.assert(self.passes.items.len != 0);
        var current: ?PresolveResult = null;
        errdefer if (current) |*value| value.deinit();
        var view = reduced_solution;
        var index = self.passes.items.len;
        while (index > 0) {
            index -= 1;
            const next = try postsolve(&self.passes.items[index], view, allocator);
            if (current) |*value| value.deinit();
            current = next;
            view = current.?.solutionView(reduced_solution.iterations);
        }
        return current.?;
    }
};

pub const PresolveError = error{OutOfMemory};

const SENTINEL = std.math.maxInt(u32);

/// Apply exact column substitutions and redundant empty-row elimination.
pub fn presolve(allocator: std.mem.Allocator, problem: ProblemView) (PresolveError || std.mem.Allocator.Error)!PresolvedProblem {
    // ── Count substitutable columns and their nonzeros ──
    var substitution_count: u32 = 0;
    var substitution_nnz: usize = 0;
    const objective_sign: f64 = @floatFromInt(@intFromEnum(problem.objective_sense));
    for (0..problem.num_cols) |j| {
        const begin = problem.matrix.col_starts[j];
        const end = problem.matrix.col_starts[j + 1];
        switch (rules.column_substitution.classify(
            end - begin,
            problem.col_cost[j],
            problem.col_lower[j],
            problem.col_upper[j],
            objective_sign,
        )) {
            .keep => {},
            .substitute => {
                substitution_count += 1;
                substitution_nnz += end - begin;
            },
        }
    }
    const reduced_cols = problem.num_cols - substitution_count;

    var adjusted_row_lower = try allocator.dupe(f64, problem.row_lower);
    defer allocator.free(adjusted_row_lower);
    var adjusted_row_upper = try allocator.dupe(f64, problem.row_upper);
    defer allocator.free(adjusted_row_upper);
    const retained_row_nnz = try allocator.alloc(u32, problem.num_rows);
    defer allocator.free(retained_row_nnz);
    @memset(retained_row_nnz, 0);
    const retained_row_column = try allocator.alloc(u32, problem.num_rows);
    defer allocator.free(retained_row_column);
    @memset(retained_row_column, SENTINEL);
    const retained_row_value = try allocator.alloc(f64, problem.num_rows);
    defer allocator.free(retained_row_value);
    @memset(retained_row_value, 0.0);

    // ── Allocate result ──
    var pp = PresolvedProblem{
        .allocator = allocator,
        .num_rows = 0,
        .num_cols = @intCast(reduced_cols),
        .original_num_rows = @intCast(problem.num_rows),
        .original_num_cols = @intCast(problem.num_cols),
        .col_cost = &.{},
        .col_lower = &.{},
        .col_upper = &.{},
        .row_lower = &.{},
        .row_upper = &.{},
        .objective_offset = problem.objective_offset,
        .objective_sense = problem.objective_sense,
        .reduced_to_original_col = &.{},
        .original_to_reduced_col = &.{},
        .reduced_to_original_row = &.{},
        .original_to_reduced_row = &.{},
        .substitution_row_indices = &.{},
        .substitution_row_values = &.{},
        .column_substitutions = &.{},
        .column_substitution_count = substitution_count,
        .singleton_rows = &.{},
        .singleton_row_count = 0,
        .removed_row_count = 0,
    };
    errdefer pp.deinit();

    pp.col_cost = try allocator.alloc(f64, reduced_cols);
    pp.col_lower = try allocator.alloc(f64, reduced_cols);
    pp.col_upper = try allocator.alloc(f64, reduced_cols);
    pp.reduced_to_original_col = try allocator.alloc(u32, reduced_cols);
    pp.original_to_reduced_col = try allocator.alloc(u32, problem.num_cols);
    pp.substitution_row_indices = try allocator.alloc(u32, substitution_nnz);
    pp.substitution_row_values = try allocator.alloc(f64, substitution_nnz);
    pp.column_substitutions = try allocator.alloc(ColumnSubstitutionRecord, substitution_count);
    pp.singleton_rows = try allocator.alloc(SingletonRowRecord, problem.num_rows);

    @memset(pp.original_to_reduced_col, SENTINEL);

    // ── Substitute columns ──
    var rc: u32 = 0;
    var fi: u32 = 0;
    var fnz: u32 = 0;

    for (0..problem.num_cols) |j| {
        const begin = problem.matrix.col_starts[j];
        const end = problem.matrix.col_starts[j + 1];
        switch (rules.column_substitution.classify(
            end - begin,
            problem.col_cost[j],
            problem.col_lower[j],
            problem.col_upper[j],
            objective_sign,
        )) {
            .substitute => |substitution| {
                pp.column_substitutions[fi] = .{
                    .original_column = @intCast(j),
                    .value = substitution.value,
                    .col_cost = problem.col_cost[j],
                    .basis_status = switch (substitution.bound) {
                        .fixed => .fixed,
                        .lower => .at_lower,
                        .upper => .at_upper,
                    },
                    .row_start = fnz,
                    .row_end = fnz + @as(u32, @intCast(end - begin)),
                };

                for (begin..end) |k| {
                    const row = problem.matrix.row_indices[k].toUsize();
                    const coeff = problem.matrix.values[k];
                    adjusted_row_lower[row] -= coeff * substitution.value;
                    adjusted_row_upper[row] -= coeff * substitution.value;
                    pp.substitution_row_indices[fnz] = @intCast(row);
                    pp.substitution_row_values[fnz] = coeff;
                    fnz += 1;
                }
                pp.objective_offset += problem.col_cost[j] * substitution.value;
                fi += 1;
            },
            .keep => {
                pp.col_cost[rc] = problem.col_cost[j];
                pp.col_lower[rc] = problem.col_lower[j];
                pp.col_upper[rc] = problem.col_upper[j];
                pp.reduced_to_original_col[rc] = @intCast(j);
                pp.original_to_reduced_col[j] = rc;
                for (begin..end) |k| {
                    const row = problem.matrix.row_indices[k].toUsize();
                    if (retained_row_nnz[row] == 0) {
                        retained_row_column[row] = @intCast(j);
                        retained_row_value[row] = problem.matrix.values[k];
                    }
                    retained_row_nnz[row] += 1;
                }
                rc += 1;
            },
        }
    }
    std.debug.assert(fi == substitution_count);
    std.debug.assert(rc == reduced_cols);

    // Remove feasible empty rows and transfer feasible singleton rows into the
    // retained column bounds. Infeasible rows remain for simplex to certify.
    const row_removed = try allocator.alloc(bool, problem.num_rows);
    defer allocator.free(row_removed);
    @memset(row_removed, false);
    var singleton_count: u32 = 0;
    for (0..problem.num_rows) |row| {
        if (retained_row_nnz[row] == 0) {
            row_removed[row] = rules.empty_row.classify(0, adjusted_row_lower[row], adjusted_row_upper[row]) == .remove;
            continue;
        }
        if (retained_row_nnz[row] != 1) continue;
        const original_column = retained_row_column[row];
        const reduced_column = pp.original_to_reduced_col[original_column];
        std.debug.assert(reduced_column != SENTINEL);
        const previous_lower = pp.col_lower[reduced_column];
        const previous_upper = pp.col_upper[reduced_column];
        switch (rules.singleton_row.classify(
            retained_row_value[row],
            adjusted_row_lower[row],
            adjusted_row_upper[row],
            previous_lower,
            previous_upper,
        )) {
            .keep, .infeasible => {},
            .reduce => |reduction| {
                pp.singleton_rows[singleton_count] = .{
                    .original_row = @intCast(row),
                    .original_column = original_column,
                    .coefficient = retained_row_value[row],
                    .previous_lower = previous_lower,
                    .previous_upper = previous_upper,
                    .tightened_lower = reduction.lower,
                    .tightened_upper = reduction.upper,
                    .lower_source = reduction.lower_source,
                    .upper_source = reduction.upper_source,
                };
                singleton_count += 1;
                pp.col_lower[reduced_column] = reduction.lower;
                pp.col_upper[reduced_column] = reduction.upper;
                row_removed[row] = true;
            },
        }
    }
    pp.singleton_row_count = singleton_count;

    var reduced_rows: u32 = 0;
    for (row_removed) |removed| reduced_rows += @intFromBool(!removed);
    pp.num_rows = reduced_rows;
    pp.removed_row_count = @intCast(problem.num_rows - reduced_rows);
    pp.row_lower = try allocator.alloc(f64, reduced_rows);
    pp.row_upper = try allocator.alloc(f64, reduced_rows);
    pp.reduced_to_original_row = try allocator.alloc(u32, reduced_rows);
    pp.original_to_reduced_row = try allocator.alloc(u32, problem.num_rows);
    @memset(pp.original_to_reduced_row, SENTINEL);
    var rr: u32 = 0;
    for (0..problem.num_rows) |row| {
        if (row_removed[row]) continue;
        pp.row_lower[rr] = adjusted_row_lower[row];
        pp.row_upper[rr] = adjusted_row_upper[row];
        pp.reduced_to_original_row[rr] = @intCast(row);
        pp.original_to_reduced_row[row] = rr;
        rr += 1;
    }
    std.debug.assert(rr == reduced_rows);

    // ── Build compact CSC ──
    var retained_nnz: usize = 0;
    for (retained_row_nnz, 0..) |count, row| {
        if (pp.original_to_reduced_row[row] != SENTINEL) retained_nnz += count;
    }
    var csc = matrix.CscMatrix.initPackedUninitialized(allocator, reduced_rows, reduced_cols, retained_nnz) catch return error.OutOfMemory;
    errdefer csc.deinit(allocator);
    pp._matrix_valid = true;

    // Count preserved nonzeros
    var offset: usize = 0;
    var col: u32 = 0;
    for (0..problem.num_cols) |j| {
        if (pp.original_to_reduced_col[j] == SENTINEL) continue;
        csc.col_starts[col] = offset;
        const begin = problem.matrix.col_starts[j];
        const end = problem.matrix.col_starts[j + 1];
        for (begin..end) |k| {
            const reduced_row = pp.original_to_reduced_row[problem.matrix.row_indices[k].toUsize()];
            if (reduced_row == SENTINEL) continue;
            csc.row_indices[offset] = foundation.RowId.fromUsizeAssumeValid(reduced_row);
            csc.values[offset] = problem.matrix.values[k];
            offset += 1;
        }
        col += 1;
    }
    csc.col_starts[reduced_cols] = offset;
    pp.matrix_storage = csc;

    return pp;
}

// ═══════════════════════════════════════════════════════════════════
// Postsolve
// ═══════════════════════════════════════════════════════════════════

/// Recover original-coordinate solution after solving the reduced problem.
/// Caller owns the returned arrays (free with their own allocator).
pub fn postsolve(
    pp: *const PresolvedProblem,
    solution: SolutionView,
    allocator: std.mem.Allocator,
) !PresolveResult {
    const orig_rows = pp.original_num_rows;
    const orig_cols = pp.original_num_cols;

    var result = PresolveResult{
        .allocator = allocator,
        .status = solution.status,
        .objective_value = solution.objective_value,
        .primal = &.{},
        .dual = &.{},
        .reduced_cost = &.{},
        .unbounded_ray = &.{},
        .infeasibility_ray = &.{},
    };
    errdefer result.deinit();
    result.primal = try allocator.alloc(f64, orig_cols);
    result.dual = try allocator.alloc(f64, orig_rows);
    result.reduced_cost = try allocator.alloc(f64, orig_cols);
    if (solution.unbounded_ray.len > 0)
        result.unbounded_ray = try allocator.alloc(f64, orig_cols);
    if (solution.infeasibility_ray.len > 0)
        result.infeasibility_ray = try allocator.alloc(f64, orig_rows);

    // ── Primal: copy preserved, set fixed ──
    @memset(result.primal, 0.0);
    for (0..pp.num_cols) |rc| {
        result.primal[pp.reduced_to_original_col[rc]] = solution.primal[rc];
    }
    for (pp.column_substitutions[0..pp.column_substitution_count]) |rec| {
        result.primal[rec.original_column] = rec.value;
    }

    // ── Dual: restore retained rows; redundant empty rows have multiplier 0 ──
    @memset(result.dual, 0.0);
    for (0..pp.num_rows) |row| {
        result.dual[pp.reduced_to_original_row[row]] = if (row < solution.dual.len) solution.dual[row] else 0.0;
    }

    // ── Reduced cost: copy preserved, then undo singleton bound transfers ──
    @memset(result.reduced_cost, 0.0);
    for (0..pp.num_cols) |rc| {
        result.reduced_cost[pp.reduced_to_original_col[rc]] = solution.reduced_cost[rc];
    }
    var singleton_index: usize = pp.singleton_row_count;
    while (singleton_index > 0) {
        singleton_index -= 1;
        const rec = pp.singleton_rows[singleton_index];
        const column = rec.original_column;
        const value = result.primal[column];
        const reduced = result.reduced_cost[column];
        const tolerance = 1e-9 * @max(1.0, @abs(value));
        const lower_active = rec.lower_source != null and @abs(value - rec.tightened_lower) <= tolerance;
        const upper_active = rec.upper_source != null and @abs(value - rec.tightened_upper) <= tolerance;
        const source = if (reduced > 0.0 and lower_active)
            rec.lower_source
        else if (reduced < 0.0 and upper_active)
            rec.upper_source
        else
            null;
        if (source != null) {
            result.dual[rec.original_row] = reduced / rec.coefficient;
            result.reduced_cost[column] = 0.0;
        }
    }

    // Substituted-column reduced costs use the engine's normalized
    // minimization sign, matching every reduced cost published by simplex.
    const objective_sign: f64 = @floatFromInt(@intFromEnum(pp.objective_sense));
    for (pp.column_substitutions[0..pp.column_substitution_count]) |rec| {
        var rc = objective_sign * rec.col_cost;
        const begin = rec.row_start;
        const end = rec.row_end;
        for (begin..end) |k| {
            rc -= result.dual[pp.substitution_row_indices[k]] * pp.substitution_row_values[k];
        }
        result.reduced_cost[rec.original_column] = rc;
    }

    // ── Ray: zero for fixed columns ──
    if (result.unbounded_ray.len > 0) {
        @memset(result.unbounded_ray, 0.0);
        for (0..pp.num_cols) |rc| {
            result.unbounded_ray[pp.reduced_to_original_col[rc]] = solution.unbounded_ray[rc];
        }
    }

    // Restore retained row coordinates in a Farkas certificate. Removed
    // feasible empty rows contribute zero. A singleton row replaced one of
    // the structural column bounds by `a*x >= b` or `a*x <= b`. If the
    // reduced proof uses that transferred bound, cancel the current column
    // coefficient q with the deleted-row multiplier y = -q/a. Its row-bound
    // contribution is exactly the removed column-bound contribution, so the
    // Farkas gap is preserved.
    if (result.infeasibility_ray.len > 0) {
        @memset(result.infeasibility_ray, 0.0);
        const count = @min(pp.num_rows, solution.infeasibility_ray.len);
        for (0..count) |row|
            result.infeasibility_ray[pp.reduced_to_original_row[row]] = solution.infeasibility_ray[row];

        const ray_column_coefficient = try allocator.alloc(f64, orig_cols);
        defer allocator.free(ray_column_coefficient);
        @memset(ray_column_coefficient, 0.0);
        for (0..pp.num_cols) |reduced_column| {
            const begin = pp.matrix_storage.col_starts[reduced_column];
            const end = pp.matrix_storage.col_starts[reduced_column + 1];
            var coefficient: f64 = 0.0;
            for (pp.matrix_storage.row_indices[begin..end], pp.matrix_storage.values[begin..end]) |row, value| {
                coefficient += result.infeasibility_ray[pp.reduced_to_original_row[row.toUsize()]] * value;
            }
            ray_column_coefficient[pp.reduced_to_original_col[reduced_column]] = coefficient;
        }

        singleton_index = pp.singleton_row_count;
        while (singleton_index > 0) {
            singleton_index -= 1;
            const rec = pp.singleton_rows[singleton_index];
            const coefficient = ray_column_coefficient[rec.original_column];
            const uses_transferred_bound = if (coefficient > 0.0)
                rec.upper_source != null
            else if (coefficient < 0.0)
                rec.lower_source != null
            else
                false;
            if (!uses_transferred_bound) continue;
            const multiplier = -coefficient / rec.coefficient;
            result.infeasibility_ray[rec.original_row] = multiplier;
            ray_column_coefficient[rec.original_column] += multiplier * rec.coefficient;
        }
    }

    return result;
}

// ═══════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════

test "problem without reducible rows or columns remains structurally unchanged" {
    const a = std.testing.allocator;
    const problem = ProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 1.0, 2.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ 10.0, 10.0 },
        .row_lower = &.{-std.math.inf(f64)},
        .row_upper = &.{5.0},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &[_]usize{ 0, 1, 2 }, &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) }, &.{ 1.0, 1.0 }),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var pp = try presolve(a, problem);
    defer pp.deinit();
    try std.testing.expect(!pp.wasApplied());
    try std.testing.expectEqual(@as(u32, 2), pp.num_cols);
    try std.testing.expectEqual(@as(u32, 1), pp.num_rows);
}

test "single fixed column elimination" {
    const a = std.testing.allocator;
    // min x1 + 2*x2  s.t.  x1 + x2 <= 5,  0 <= x1, x2 = 3
    const problem = ProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 1.0, 2.0 },
        .col_lower = &.{ 0.0, 3.0 },
        .col_upper = &.{ std.math.inf(f64), 3.0 },
        .row_lower = &.{-std.math.inf(f64)},
        .row_upper = &.{5.0},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &[_]usize{ 0, 1, 2 }, &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) }, &.{ 1.0, 1.0 }),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var pp = try presolve(a, problem);
    defer pp.deinit();

    try std.testing.expect(pp.wasApplied());
    try std.testing.expectEqual(@as(u32, 1), pp.num_cols);
    try std.testing.expectEqual(@as(u32, 1), pp.column_substitution_count);
    // offset = 2*3 = 6; the remaining singleton transfers x1 <= 2.
    try std.testing.expectApproxEqAbs(@as(f64, 6.0), pp.objective_offset, 1e-12);
    try std.testing.expectEqual(@as(u32, 0), pp.num_rows);
    try std.testing.expectEqual(@as(u32, 1), pp.singleton_row_count);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), pp.col_upper[0], 1e-12);

    // Postsolve
    const sol = SolutionView{
        .status = .optimal,
        .objective_value = 6.0,
        .iterations = 1,
        .primal = &.{0.0},
        .dual = &.{},
        .reduced_cost = &.{1.0},
        .unbounded_ray = &.{},
        .infeasibility_ray = &.{},
    };
    var res = try postsolve(&pp, sol, a);
    defer res.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), res.primal[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), res.primal[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 6.0), res.objective_value, 1e-12);
}

test "rectangular fixed-column postsolve restores every coordinate" {
    const a = std.testing.allocator;
    // Three rows deliberately differ from the two columns. The retained
    // column has one nonzero while the fixed column has two, exercising both
    // original dimension tracking and retained-nnz allocation.
    const problem = ProblemView{
        .num_rows = 3,
        .num_cols = 2,
        .col_cost = &.{ 2.0, 7.0 },
        .col_lower = &.{ 0.0, 4.0 },
        .col_upper = &.{ 10.0, 4.0 },
        .row_lower = &.{ -std.math.inf(f64), 1.0, -3.0 },
        .row_upper = &.{ 8.0, 9.0, 5.0 },
        .matrix = matrix.CscView.initAssumeValid(
            3,
            2,
            &[_]usize{ 0, 1, 3 },
            &.{
                foundation.RowId.fromUsizeAssumeValid(1),
                foundation.RowId.fromUsizeAssumeValid(0),
                foundation.RowId.fromUsizeAssumeValid(2),
            },
            &.{ 3.0, 2.0, -1.0 },
        ),
        .objective_sense = .maximize,
        .objective_offset = 11.0,
    };
    var pp = try presolve(a, problem);
    defer pp.deinit();

    try std.testing.expectEqual(@as(u32, 3), pp.original_num_rows);
    try std.testing.expectEqual(@as(u32, 2), pp.original_num_cols);
    try std.testing.expectEqual(@as(u32, 1), pp.num_rows);
    try std.testing.expectEqual(@as(u32, 2), pp.removed_row_count);
    try std.testing.expectEqual(@as(u32, 1), pp.singleton_row_count);
    try std.testing.expectEqualSlices(u32, &.{2}, pp.reduced_to_original_row);
    try std.testing.expectEqual(@as(usize, 0), pp.matrix_storage.nnz());
    try std.testing.expectApproxEqAbs(@as(f64, 39.0), pp.objective_offset, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 9.0), pp.row_upper[0], 1e-12);

    const sol = SolutionView{
        .status = .infeasible,
        .objective_value = 39.0,
        .iterations = 2,
        .primal = &.{1.5},
        .dual = &.{3.0},
        .reduced_cost = &.{8.0},
        .unbounded_ray = &.{},
        .infeasibility_ray = &.{2.0},
    };
    var res = try postsolve(&pp, sol, a);
    defer res.deinit();

    try std.testing.expectEqual(@as(usize, 2), res.primal.len);
    try std.testing.expectEqual(@as(usize, 3), res.dual.len);
    try std.testing.expectApproxEqAbs(@as(f64, 1.5), res.primal[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), res.primal[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 8.0), res.reduced_cost[0], 1e-12);
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 0.0, 3.0 }, res.dual);
    // Maximize reduced costs use the normalized minimization sign:
    // -c_j - A_j^T y = -7 - (-1*3) = -4.
    try std.testing.expectApproxEqAbs(@as(f64, -4.0), res.reduced_cost[1], 1e-12);
    // The retained infeasible empty row is restored while the unrelated
    // singleton row contributes zero.
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 0.0, 2.0 }, res.infeasibility_ray);
}

test "postsolve expands an unbounded ray with zero fixed direction" {
    const a = std.testing.allocator;
    const problem = ProblemView{
        .num_rows = 0,
        .num_cols = 2,
        .col_cost = &.{ -1.0, 3.0 },
        .col_lower = &.{ 0.0, 2.0 },
        .col_upper = &.{ std.math.inf(f64), 2.0 },
        .row_lower = &.{},
        .row_upper = &.{},
        .matrix = matrix.CscView.initAssumeValid(0, 2, &[_]usize{ 0, 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var pp = try presolve(a, problem);
    defer pp.deinit();
    const sol = SolutionView{
        .status = .unbounded,
        .objective_value = 6.0,
        .iterations = 0,
        .primal = &.{0.0},
        .dual = &.{},
        .reduced_cost = &.{-1.0},
        .unbounded_ray = &.{1.0},
        .infeasibility_ray = &.{},
    };
    var res = try postsolve(&pp, sol, a);
    defer res.deinit();
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 0.0 }, res.unbounded_ray);
}

test "two fixed columns" {
    const a = std.testing.allocator;
    const problem = ProblemView{
        .num_rows = 1,
        .num_cols = 3,
        .col_cost = &.{ 1.0, 2.0, 3.0 },
        .col_lower = &.{ 1.0, 2.0, 0.0 },
        .col_upper = &.{ 1.0, 2.0, 10.0 },
        .row_lower = &.{-std.math.inf(f64)},
        .row_upper = &.{20.0},
        .matrix = matrix.CscView.initAssumeValid(1, 3, &[_]usize{ 0, 1, 2, 3 }, &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) }, &.{ 2.0, 3.0, 1.0 }),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var pp = try presolve(a, problem);
    defer pp.deinit();
    try std.testing.expectEqual(@as(u32, 2), pp.column_substitution_count);
    // offset = 1+4+0=5; the remaining singleton row is redundant with x3 <= 10.
    try std.testing.expectApproxEqAbs(@as(f64, 5.0), pp.objective_offset, 1e-12);
    try std.testing.expectEqual(@as(u32, 0), pp.num_rows);
    try std.testing.expectEqual(@as(u32, 1), pp.singleton_row_count);
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), pp.col_upper[0], 1e-12);
}

test "finite-inf bounds not fixed" {
    const a = std.testing.allocator;
    const problem = ProblemView{
        .num_rows = 0,
        .num_cols = 1,
        .col_cost = &.{1.0},
        .col_lower = &.{-std.math.inf(f64)},
        .col_upper = &.{std.math.inf(f64)},
        .row_lower = &.{},
        .row_upper = &.{},
        .matrix = matrix.CscView.initAssumeValid(0, 1, &[_]usize{ 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var pp = try presolve(a, problem);
    defer pp.deinit();
    try std.testing.expect(!pp.wasApplied());
}
