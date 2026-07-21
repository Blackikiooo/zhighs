//! Reversible LP presolve operating at the ProblemView level.
//!
//! Reduces a problem by removing redundant structure, then postsolves
//! the simplex engine output back to original coordinates. Each reduction
//! records the transformation so postsolve is exact (not approximate).
//!
//! Phase 1: Fixed column elimination.
//! A column where col_lower == col_upper is deterministic — its value is
//! known. The column's contribution is moved into the RHS and objective
//! offset, and the column is deleted from the reduced problem.

const std = @import("std");
const foundation = @import("foundation");
const matrix = @import("matrix");
const lp = @import("lp");
const ProblemView = lp.ProblemView;
const SolveStatus = lp.SolveStatus;
const SolutionView = lp.SolutionView;

/// Record for one eliminated fixed column.
pub const FixedColumnRecord = struct {
    original_column: u32,
    fixed_value: f64,
    row_start: u32,
    row_end: u32,
};

/// Owns reduced problem data and transformation records for postsolve.
pub const PresolvedProblem = struct {
    allocator: std.mem.Allocator,

    /// Reduced dimensions.
    num_rows: u32,
    num_cols: u32,

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

    // ── Fixed column records ──
    fixed_row_indices: []u32,
    fixed_row_values: []f64,
    fixed_columns: []FixedColumnRecord,
    fixed_column_count: u32,

    pub fn wasApplied(self: *const PresolvedProblem) bool {
        return self.fixed_column_count > 0;
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
        self.allocator.free(self.fixed_row_indices);
        self.allocator.free(self.fixed_row_values);
        self.allocator.free(self.fixed_columns);
    }
};

/// Result of presolve → solve → postsolve round-trip.
pub const PresolveResult = struct {
    status: SolveStatus,
    objective_value: f64,
    primal: []f64,
    dual: []f64,
    reduced_cost: []f64,
    unbounded_ray: []f64,
};

pub const PresolveError = error{ OutOfMemory };

const SENTINEL = std.math.maxInt(u32);

/// Apply fixed-column elimination. Returns a trivial wrapper when no fixed
/// columns exist (so callers don't need a separate code path).
pub fn presolve(allocator: std.mem.Allocator, problem: ProblemView) (PresolveError || std.mem.Allocator.Error)!PresolvedProblem {
    // ── Count fixed columns and their nonzeros ──
    var fixed: u32 = 0;
    var fixed_nnz: usize = 0;
    for (0..problem.num_cols) |j| {
        if (problem.col_lower[j] == problem.col_upper[j] and std.math.isFinite(problem.col_lower[j])) {
            fixed += 1;
            fixed_nnz += problem.matrix.col_starts[j + 1] - problem.matrix.col_starts[j];
        }
    }
    if (fixed == 0) return triviallyPresolved(allocator, problem);

    const reduced_cols = problem.num_cols - fixed;

    // ── Allocate result ──
    var pp = PresolvedProblem{
        .allocator = allocator,
        .num_rows = @intCast(problem.num_rows),
        .num_cols = @intCast(reduced_cols),
        .col_cost = try allocator.alloc(f64, reduced_cols),
        .col_lower = try allocator.alloc(f64, reduced_cols),
        .col_upper = try allocator.alloc(f64, reduced_cols),
        .row_lower = try allocator.alloc(f64, problem.num_rows),
        .row_upper = try allocator.alloc(f64, problem.num_rows),
        .objective_offset = problem.objective_offset,
        .objective_sense = problem.objective_sense,
        .reduced_to_original_col = try allocator.alloc(u32, reduced_cols),
        .original_to_reduced_col = try allocator.alloc(u32, problem.num_cols),
        .fixed_row_indices = try allocator.alloc(u32, fixed_nnz),
        .fixed_row_values = try allocator.alloc(f64, fixed_nnz),
        .fixed_columns = try allocator.alloc(FixedColumnRecord, fixed),
        .fixed_column_count = fixed,
    };
    errdefer pp.deinit();

    @memset(pp.original_to_reduced_col, SENTINEL);
    @memcpy(pp.row_lower, problem.row_lower);
    @memcpy(pp.row_upper, problem.row_upper);

    // ── Eliminate fixed columns ──
    var rc: u32 = 0;
    var fi: u32 = 0;
    var fnz: u32 = 0;

    for (0..problem.num_cols) |j| {
        if (problem.col_lower[j] == problem.col_upper[j] and std.math.isFinite(problem.col_lower[j])) {
            const value = problem.col_lower[j];
            const begin = problem.matrix.col_starts[j];
            const end = problem.matrix.col_starts[j + 1];

            pp.fixed_columns[fi] = .{
                .original_column = @as(u32, @intCast(j)),
                .fixed_value = value,
                .row_start = fnz,
                .row_end = fnz + @as(u32, @intCast(end - begin)),
            };

            for (begin..end) |k| {
                const row = problem.matrix.row_indices[k].toUsize();
                const coeff = problem.matrix.values[k];
                pp.row_lower[row] -= coeff * value;
                pp.row_upper[row] -= coeff * value;
                pp.fixed_row_indices[fnz] = @intCast(row);
                pp.fixed_row_values[fnz] = coeff;
                fnz += 1;
            }
            pp.objective_offset += problem.col_cost[j] * value;
            fi += 1;
        } else {
            pp.col_cost[rc] = problem.col_cost[j];
            pp.col_lower[rc] = problem.col_lower[j];
            pp.col_upper[rc] = problem.col_upper[j];
            pp.reduced_to_original_col[rc] = @intCast(j);
            pp.original_to_reduced_col[j] = rc;
            rc += 1;
        }
    }
    std.debug.assert(fi == fixed);
    std.debug.assert(rc == reduced_cols);

    // ── Build compact CSC ──
    var csc = matrix.CscMatrix.initPackedUninitialized(allocator, problem.num_rows, reduced_cols, fixed_nnz) catch return error.OutOfMemory;
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
            csc.row_indices[offset] = problem.matrix.row_indices[k];
            csc.values[offset] = problem.matrix.values[k];
            offset += 1;
        }
        col += 1;
    }
    csc.col_starts[reduced_cols] = offset;
    pp.matrix_storage = csc;

    return pp;
}

/// Build a trivial presolved problem that clones the original.
fn triviallyPresolved(allocator: std.mem.Allocator, problem: ProblemView) (PresolveError || std.mem.Allocator.Error)!PresolvedProblem {
    const nc = problem.num_cols;
    const nr = problem.num_rows;
    var pp = PresolvedProblem{
        .allocator = allocator,
        .num_rows = @intCast(nr),
        .num_cols = @intCast(nc),
        .col_cost = try allocator.dupe(f64, problem.col_cost),
        .col_lower = try allocator.dupe(f64, problem.col_lower),
        .col_upper = try allocator.dupe(f64, problem.col_upper),
        .row_lower = try allocator.dupe(f64, problem.row_lower),
        .row_upper = try allocator.dupe(f64, problem.row_upper),
        .objective_offset = problem.objective_offset,
        .objective_sense = problem.objective_sense,
        .reduced_to_original_col = try allocator.alloc(u32, nc),
        .original_to_reduced_col = try allocator.alloc(u32, nc),
        .fixed_row_indices = &.{},
        .fixed_row_values = &.{},
        .fixed_columns = &.{},
        .fixed_column_count = 0,
    };
    errdefer pp.deinit();
    for (0..nc) |j| {
        pp.reduced_to_original_col[j] = @intCast(j);
        pp.original_to_reduced_col[j] = @intCast(j);
    }
    pp.matrix_storage = matrix.CscMatrix.initPackedUninitialized(allocator, nr, nc, problem.matrix.values.len) catch return error.OutOfMemory;
    pp._matrix_valid = true;
    for (0..nc + 1) |j| pp.matrix_storage.col_starts[j] = problem.matrix.col_starts[j];
    @memcpy(pp.matrix_storage.row_indices[0..problem.matrix.values.len], problem.matrix.row_indices);
    @memcpy(pp.matrix_storage.values[0..problem.matrix.values.len], problem.matrix.values);
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
    const orig_rows = pp.original_to_reduced_col.len; // also num rows (identity mapping in Phase 1)
    const orig_cols = pp.original_to_reduced_col.len;

    var result = PresolveResult{
        .status = solution.status,
        .objective_value = solution.objective_value,
        .primal = try allocator.alloc(f64, orig_cols),
        .dual = try allocator.alloc(f64, orig_rows),
        .reduced_cost = try allocator.alloc(f64, orig_cols),
        .unbounded_ray = try allocator.alloc(f64, orig_cols),
    };
    errdefer {
        allocator.free(result.primal);
        allocator.free(result.dual);
        allocator.free(result.reduced_cost);
        allocator.free(result.unbounded_ray);
    }

    // ── Primal: copy preserved, set fixed ──
    @memset(result.primal, 0.0);
    for (0..pp.num_cols) |rc| {
        result.primal[pp.reduced_to_original_col[rc]] = solution.primal[rc];
    }
    for (pp.fixed_columns[0..pp.fixed_column_count]) |rec| {
        result.primal[rec.original_column] = rec.fixed_value;
    }

    // ── Dual: rows unchanged in Phase 1 ──
    @memset(result.dual, 0.0);
    for (0..orig_rows) |i| {
        result.dual[i] = if (i < solution.dual.len) solution.dual[i] else 0.0;
    }

    // ── Reduced cost: copy preserved, compute for fixed ──
    @memset(result.reduced_cost, 0.0);
    for (0..pp.num_cols) |rc| {
        result.reduced_cost[pp.reduced_to_original_col[rc]] = solution.reduced_cost[rc];
    }
    for (pp.fixed_columns[0..pp.fixed_column_count]) |rec| {
        var rc: f64 = 0.0;
        const begin = rec.row_start;
        const end = rec.row_end;
        for (begin..end) |k| {
            rc -= result.dual[pp.fixed_row_indices[k]] * pp.fixed_row_values[k];
        }
        result.reduced_cost[rec.original_column] = rc;
    }

    // ── Ray: zero for fixed columns ──
    @memset(result.unbounded_ray, 0.0);
    if (solution.status == .unbounded and solution.unbounded_ray.len > 0) {
        for (0..pp.num_cols) |rc| {
            result.unbounded_ray[pp.reduced_to_original_col[rc]] = solution.unbounded_ray[rc];
        }
    }

    return result;
}

// ═══════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════

test "no fixed columns → trivial" {
    const a = std.testing.allocator;
    const problem = ProblemView{
        .num_rows = 1, .num_cols = 1,
        .col_cost = &.{1.0}, .col_lower = &.{0.0}, .col_upper = &.{10.0},
        .row_lower = &.{-std.math.inf(f64)}, .row_upper = &.{5.0},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &[_]usize{ 0, 1 }, &.{foundation.RowId.fromUsizeAssumeValid(0)}, &.{1.0}),
        .objective_sense = .minimize, .objective_offset = 0.0,
    };
    var pp = try presolve(a, problem);
    defer pp.deinit();
    try std.testing.expect(!pp.wasApplied());
    try std.testing.expectEqual(@as(u32, 1), pp.num_cols);
}

test "single fixed column elimination" {
    const a = std.testing.allocator;
    // min x1 + 2*x2  s.t.  x1 + x2 <= 5,  0 <= x1, x2 = 3
    const problem = ProblemView{
        .num_rows = 1, .num_cols = 2,
        .col_cost = &.{ 1.0, 2.0 },
        .col_lower = &.{ 0.0, 3.0 },
        .col_upper = &.{ std.math.inf(f64), 3.0 },
        .row_lower = &.{-std.math.inf(f64)},
        .row_upper = &.{5.0},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &[_]usize{ 0, 1, 2 },
            &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) },
            &.{ 1.0, 1.0 }),
        .objective_sense = .minimize, .objective_offset = 0.0,
    };
    var pp = try presolve(a, problem);
    defer pp.deinit();

    try std.testing.expect(pp.wasApplied());
    try std.testing.expectEqual(@as(u32, 1), pp.num_cols);
    try std.testing.expectEqual(@as(u32, 1), pp.fixed_column_count);
    // offset = 2*3 = 6, row upper = 5 - 1*3 = 2
    try std.testing.expectApproxEqAbs(@as(f64, 6.0), pp.objective_offset, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), pp.row_upper[0], 1e-12);

    // Postsolve
    const sol = SolutionView{
        .status = .optimal, .objective_value = 6.0, .iterations = 1,
        .primal = &.{0.0}, .dual = &.{0.0}, .reduced_cost = &.{1.0}, .unbounded_ray = &.{},
        .infeasibility_ray = &.{},
    };
    const res = try postsolve(&pp, sol, a);
    defer {
        a.free(res.primal); a.free(res.dual); a.free(res.reduced_cost); a.free(res.unbounded_ray);
    }
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), res.primal[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), res.primal[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 6.0), res.objective_value, 1e-12);
}

test "two fixed columns" {
    const a = std.testing.allocator;
    const problem = ProblemView{
        .num_rows = 1, .num_cols = 3,
        .col_cost = &.{ 1.0, 2.0, 3.0 },
        .col_lower = &.{ 1.0, 2.0, 0.0 },
        .col_upper = &.{ 1.0, 2.0, 10.0 },
        .row_lower = &.{-std.math.inf(f64)},
        .row_upper = &.{20.0},
        .matrix = matrix.CscView.initAssumeValid(1, 3, &[_]usize{ 0, 1, 2, 3 },
            &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) },
            &.{ 2.0, 3.0, 1.0 }),
        .objective_sense = .minimize, .objective_offset = 0.0,
    };
    var pp = try presolve(a, problem);
    defer pp.deinit();
    try std.testing.expectEqual(@as(u32, 2), pp.fixed_column_count);
    // offset = 1+4+0=5, row upper = 20-2-6=12
    try std.testing.expectApproxEqAbs(@as(f64, 5.0), pp.objective_offset, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 12.0), pp.row_upper[0], 1e-12);
}

test "finite-inf bounds not fixed" {
    const a = std.testing.allocator;
    const problem = ProblemView{
        .num_rows = 0, .num_cols = 1,
        .col_cost = &.{1.0}, .col_lower = &.{-std.math.inf(f64)}, .col_upper = &.{std.math.inf(f64)},
        .row_lower = &.{}, .row_upper = &.{},
        .matrix = matrix.CscView.initAssumeValid(0, 1, &[_]usize{ 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize, .objective_offset = 0.0,
    };
    var pp = try presolve(a, problem);
    defer pp.deinit();
    try std.testing.expect(!pp.wasApplied());
}
