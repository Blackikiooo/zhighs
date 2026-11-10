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
//! Phase 2: Free zero-cost singleton-column elimination.
//! One such column can satisfy its only row at any activity, so the row and
//! column are removed together and reconstructed as a basic pair.
//!
//! Phase 3: Redundant empty-row elimination.
//! Rows with no retained coefficients and bounds containing zero are removed.
//! Original row coordinates are retained so duals, certificates, and bases can
//! be reconstructed exactly during postsolve.
//!
//! Phase 4: Singleton-row bound transfer.
//! A row containing one retained coefficient is converted into bounds on that
//! column and removed. Postsolve records the originating row side so dual and
//! basis coordinates can be restored.
//!
//! Phase 5: Activity-implied redundant-row elimination.
//! Outward-rounded row activity intervals prove that every point in the
//! retained column box satisfies a row before it is removed.
//!
//! Phase 6: Exact duplicate-row elimination.
//! Rows with identical adjusted bounds and coefficient sequences share one
//! retained representative; removed copies recover zero multipliers.
//!
//! MIP mode (`presolveMip`, `PresolvePipeline.initMip`) restricts column-level
//! reductions to ones that remain exact under integrality: free-singleton
//! elimination applies to continuous columns only, and a discrete column is
//! substituted only at an integer-feasible value (empty-column bounds are
//! rounded into the feasible domain). Row-level rules are purely
//! feasibility-based and apply unchanged.

const std = @import("std");
const foundation = @import("foundation");
const matrix = @import("matrix");
const lp = @import("lp");
const rules = @import("rules/root.zig");
const ProblemView = lp.ProblemView;
const SolveStatus = lp.SolveStatus;
const SolutionView = lp.SolutionView;

/// Integrality summary of one column as presolve needs it. MIP mode only
/// distinguishes binary and general discrete columns, whose
/// substituted values must be integer-feasible, from continuous ones.
pub const ColumnNature = enum {
    continuous,
    discrete,
    binary,
};

/// Rule switches used by benchmark/ablation callers. The regular `presolve`,
/// `presolveMip`, and pipeline entry points keep every rule enabled.
pub const PresolveOptions = struct {
    /// Row-activity propagation that commits only forced (collapsed) fixings.
    implied_fixing: bool = true,
    /// Remove rows that are exact non-zero scalar multiples of another row.
    parallel_rows: bool = true,
};

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

pub const FreeSingletonColumnSide = enum {
    fixed,
    lower,
    upper,
    free,
};

/// One free zero-cost singleton column paired with its deleted row. Row-entry
/// storage excludes the singleton itself and is used to reconstruct primal and
/// unbounded-ray coordinates.
pub const FreeSingletonColumnRecord = struct {
    original_column: u32,
    original_row: u32,
    coefficient: f64,
    side: FreeSingletonColumnSide,
    target_activity: f64,
    row_start: u32,
    row_end: u32,
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
    free_singleton_row_columns: []u32,
    free_singleton_row_values: []f64,
    free_singleton_columns: []FreeSingletonColumnRecord,
    free_singleton_column_count: u32,
    singleton_rows: []SingletonRowRecord,
    singleton_row_count: u32,
    redundant_row_count: u32,
    duplicate_row_count: u32,
    implied_fixed_column_count: u32,
    parallel_row_count: u32,
    removed_row_count: u32,

    pub fn wasApplied(self: *const PresolvedProblem) bool {
        return self.column_substitution_count > 0 or self.free_singleton_column_count > 0 or self.removed_row_count > 0;
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
        self.allocator.free(self.free_singleton_row_columns);
        self.allocator.free(self.free_singleton_row_values);
        self.allocator.free(self.free_singleton_columns);
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
    free_singleton_column_count: u32 = 0,
    singleton_row_count: u32 = 0,
    redundant_row_count: u32 = 0,
    duplicate_row_count: u32 = 0,
    implied_fixed_column_count: u32 = 0,
    parallel_row_count: u32 = 0,
    removed_row_count: u32 = 0,

    pub fn init(allocator: std.mem.Allocator, problem: *const ProblemView) !PresolvePipeline {
        return initWithOptions(allocator, problem, .{});
    }

    pub fn initWithOptions(allocator: std.mem.Allocator, problem: *const ProblemView, options: PresolveOptions) !PresolvePipeline {
        return initInternal(allocator, problem, null, options);
    }

    /// MIP-safe pipeline: `nature` marks which original columns are discrete.
    /// Column-nature information is remapped through every pass so each pass
    /// applies only integrality-safe reductions.
    pub fn initMip(allocator: std.mem.Allocator, problem: *const ProblemView, nature: []const ColumnNature) !PresolvePipeline {
        return initMipWithOptions(allocator, problem, nature, .{});
    }

    pub fn initMipWithOptions(
        allocator: std.mem.Allocator,
        problem: *const ProblemView,
        nature: []const ColumnNature,
        options: PresolveOptions,
    ) !PresolvePipeline {
        std.debug.assert(nature.len == problem.num_cols);
        return initInternal(allocator, problem, nature, options);
    }

    fn initInternal(
        allocator: std.mem.Allocator,
        problem: *const ProblemView,
        nature: ?[]const ColumnNature,
        options: PresolveOptions,
    ) !PresolvePipeline {
        var result = PresolvePipeline{ .allocator = allocator };
        errdefer result.deinit();
        const maximum_passes = problem.num_rows + problem.num_cols + 1;
        try result.passes.ensureTotalCapacity(allocator, maximum_passes);

        var current_nature: ?[]ColumnNature = if (nature) |source| try allocator.dupe(ColumnNature, source) else null;
        defer if (current_nature) |slice| allocator.free(slice);

        var current: ProblemView = problem.*;
        while (result.passes.items.len < maximum_passes) {
            var pass = try presolveImpl(allocator, &current, current_nature, options);
            if (!pass.wasApplied()) {
                pass.deinit();
                break;
            }
            std.debug.assert(pass.num_rows < current.num_rows or pass.num_cols < current.num_cols);
            if (current_nature) |slice| {
                const next_nature = try allocator.alloc(ColumnNature, pass.num_cols);
                for (next_nature, pass.reduced_to_original_col) |*target, source| target.* = slice[source];
                allocator.free(slice);
                current_nature = next_nature;
            }
            result.column_substitution_count += pass.column_substitution_count;
            result.free_singleton_column_count += pass.free_singleton_column_count;
            result.singleton_row_count += pass.singleton_row_count;
            result.redundant_row_count += pass.redundant_row_count;
            result.duplicate_row_count += pass.duplicate_row_count;
            result.implied_fixed_column_count += pass.implied_fixed_column_count;
            result.parallel_row_count += pass.parallel_row_count;
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

const INTEGER_TOLERANCE: f64 = 1e-6;

/// MIP-safe column-substitution classification. Continuous columns use the
/// LP rule unchanged. A discrete column is substituted only at an
/// integer-feasible value: a fractional fixed value would hide an infeasible
/// domain, so the column is kept; an empty column's objective-optimal bound is
/// rounded to the nearest feasible integer, which stays optimal.
fn classifySubstitution(
    nonzero_count: usize,
    cost: f64,
    lower: f64,
    upper: f64,
    objective_sign: f64,
    nature: ?ColumnNature,
) rules.column_substitution.Action {
    if (nature == .binary) {
        const binary_lower = @max(lower, 0.0);
        const binary_upper = @min(upper, 1.0);
        // Keep an empty binary domain visible to the MIP infeasibility check.
        if (binary_lower > binary_upper) return .keep;
        return classifySubstitution(nonzero_count, cost, binary_lower, binary_upper, objective_sign, .discrete);
    }
    const action = rules.column_substitution.classify(nonzero_count, cost, lower, upper, objective_sign);
    if (nature != .discrete) return action;
    switch (action) {
        .keep => return .keep,
        .substitute => |substitution| {
            switch (substitution.bound) {
                .fixed => {
                    const rounded = @round(substitution.value);
                    if (@abs(rounded - substitution.value) > INTEGER_TOLERANCE) return .keep;
                    return .{ .substitute = .{ .value = rounded, .bound = .fixed } };
                },
                .lower => {
                    const value = @ceil(lower - INTEGER_TOLERANCE);
                    if (value > upper) return .keep;
                    return .{ .substitute = .{ .value = value, .bound = .lower } };
                },
                .upper => {
                    const value = @floor(upper + INTEGER_TOLERANCE);
                    if (value < lower) return .keep;
                    return .{ .substitute = .{ .value = value, .bound = .upper } };
                },
                .free => return .{ .substitute = .{ .value = 0.0, .bound = .free } },
            }
        },
    }
}

/// Apply the forced fixing override first, then the original column rule.
fn classifyColumn(
    implied_fixings: ?[]const ?f64,
    column: usize,
    nonzero_count: usize,
    cost: f64,
    lower: f64,
    upper: f64,
    objective_sign: f64,
    nature: ?ColumnNature,
) rules.column_substitution.Action {
    if (implied_fixings) |values| {
        if (values[column]) |value| return .{ .substitute = .{ .value = value, .bound = .fixed } };
    }
    return classifySubstitution(nonzero_count, cost, lower, upper, objective_sign, nature);
}

const MAX_PROPAGATION_PASSES: u32 = 8;

/// Propagate row activities until no bound or fixing changes, and return the
/// forced value of every column whose implied interval collapses. Only fixings
/// are returned: temporary bound tightenings are internal to propagation and
/// never leak into the reduced model, which keeps dual postsolve exact.
fn computeImpliedFixings(
    allocator: std.mem.Allocator,
    problem: *const ProblemView,
    nature: ?[]const ColumnNature,
) ![]?f64 {
    const result = try allocator.alloc(?f64, problem.num_cols);
    errdefer allocator.free(result);
    @memset(result, null);
    if (problem.num_rows == 0 or problem.num_cols == 0) return result;

    const row_counts = try allocator.alloc(u32, problem.num_rows);
    defer allocator.free(row_counts);
    @memset(row_counts, 0);
    for (problem.matrix.row_indices) |row| row_counts[row.toUsize()] += 1;

    const row_starts = try allocator.alloc(usize, problem.num_rows + 1);
    defer allocator.free(row_starts);
    row_starts[0] = 0;
    for (row_counts, 0..) |count, row| row_starts[row + 1] = row_starts[row] + count;

    const row_columns = try allocator.alloc(u32, problem.matrix.values.len);
    defer allocator.free(row_columns);
    const row_values = try allocator.alloc(f64, problem.matrix.values.len);
    defer allocator.free(row_values);
    const row_cursor = try allocator.dupe(usize, row_starts[0..problem.num_rows]);
    defer allocator.free(row_cursor);
    for (0..problem.num_cols) |column| {
        const begin = problem.matrix.col_starts[column];
        const end = problem.matrix.col_starts[column + 1];
        for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row_id, value| {
            const row = row_id.toUsize();
            row_columns[row_cursor[row]] = @intCast(column);
            row_values[row_cursor[row]] = value;
            row_cursor[row] += 1;
        }
    }

    const lower = try allocator.dupe(f64, problem.col_lower);
    defer allocator.free(lower);
    const upper = try allocator.dupe(f64, problem.col_upper);
    defer allocator.free(upper);
    const fixed = try allocator.alloc(bool, problem.num_cols);
    defer allocator.free(fixed);
    @memset(fixed, false);
    const fixed_values = try allocator.alloc(f64, problem.num_cols);
    defer allocator.free(fixed_values);
    @memset(fixed_values, 0.0);
    const activity = try allocator.alloc(rules.redundant_row.ActivityBounds, problem.num_rows);
    defer allocator.free(activity);
    // Integer domains can be rounded before the first pass. This is the only
    // bound information that is safe to propagate without a postsolve record,
    // because it is never published: it only strengthens the fixing test.
    for (0..problem.num_cols) |column| {
        const kind = columnKind(nature, column);
        if (kind == .continuous) continue;
        if (kind == .binary) {
            lower[column] = @max(lower[column], 0.0);
            upper[column] = @min(upper[column], 1.0);
        }
        lower[column] = @ceil(lower[column]);
        upper[column] = @floor(upper[column]);
    }

    var pass: u32 = 0;
    var changed = true;
    while (changed and pass < MAX_PROPAGATION_PASSES) : (pass += 1) {
        changed = false;
        for (activity) |*bounds| bounds.* = .{};
        for (0..problem.num_cols) |column| {
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row_id, coefficient|
                activity[row_id.toUsize()].add(coefficient, lower[column], upper[column]);
        }

        for (0..problem.num_rows) |row| {
            const total = activity[row];
            const row_lower = problem.row_lower[row];
            const row_upper = problem.row_upper[row];
            if (std.math.isNan(total.lower) or std.math.isNan(total.upper) or
                std.math.isNan(row_lower) or std.math.isNan(row_upper)) continue;
            // An infeasible row under the current box is left to simplex to
            // certify; it must not drive fixings.
            if (total.upper < row_lower or total.lower > row_upper) continue;

            for (row_starts[row]..row_starts[row + 1]) |entry| {
                const column = row_columns[entry];
                if (fixed[column]) continue;
                const coefficient = row_values[entry];
                // total is the sum of per-term interval endpoints. Excluding
                // this term uses its minimum rounded up (respectively maximum
                // rounded down), which keeps both remaining bounds rigorous.
                const others = rules.implied_fixing.Bounds{
                    .lower = rules.redundant_row.outwardAdd(total.lower, -rules.implied_fixing.minimumUpper(coefficient, lower[column], upper[column]), -std.math.inf(f64)),
                    .upper = rules.redundant_row.outwardAdd(total.upper, -rules.implied_fixing.maximumLower(coefficient, lower[column], upper[column]), std.math.inf(f64)),
                };
                const implied = rules.implied_fixing.impliedBounds(coefficient, row_lower, row_upper, others);
                if (std.math.isNan(implied.lower) or std.math.isNan(implied.upper)) continue;

                var new_lower = @max(lower[column], implied.lower);
                var new_upper = @min(upper[column], implied.upper);
                const kind = columnKind(nature, column);
                if (kind == .binary) {
                    new_lower = @max(new_lower, 0.0);
                    new_upper = @min(new_upper, 1.0);
                }
                if (kind != .continuous) {
                    new_lower = @ceil(new_lower);
                    new_upper = @floor(new_upper);
                }
                // Outward rounding can leave a small gap even when the true
                // interval collapses; only commit a fixing when it is gone.
                if (new_lower > new_upper) continue;

                if (new_lower != lower[column] or new_upper != upper[column]) {
                    lower[column] = new_lower;
                    upper[column] = new_upper;
                    changed = true;
                }
                if (rules.implied_fixing.fixValue(.{ .lower = new_lower, .upper = new_upper }, kind)) |value| {
                    fixed[column] = true;
                    fixed_values[column] = value;
                    lower[column] = value;
                    upper[column] = value;
                    changed = true;
                }
            }
        }
    }

    for (0..problem.num_cols) |column| {
        if (fixed[column]) result[column] = fixed_values[column];
    }
    return result;
}

fn columnKind(nature: ?[]const ColumnNature, column: usize) rules.implied_fixing.Kind {
    const kinds = nature orelse return .continuous;
    return switch (kinds[column]) {
        .continuous => .continuous,
        .discrete => .integer,
        .binary => .binary,
    };
}

/// Apply exact column substitutions and redundant empty-row elimination.


pub fn presolve(allocator: std.mem.Allocator, problem: *const ProblemView) (PresolveError || std.mem.Allocator.Error)!PresolvedProblem {
    return presolveImpl(allocator, problem, null, .{});
}

/// LP variant with explicit rule switches for ablation benchmarks.
pub fn presolveWithOptions(
    allocator: std.mem.Allocator,
    problem: *const ProblemView,
    options: PresolveOptions,
) (PresolveError || std.mem.Allocator.Error)!PresolvedProblem {
    return presolveImpl(allocator, problem, null, options);
}

/// MIP-safe variant: `nature[column]` marks discrete columns. Discrete
/// columns are substituted only at integer-feasible values, and
/// free-singleton elimination applies to continuous columns only.
pub fn presolveMip(allocator: std.mem.Allocator, problem: *const ProblemView, nature: []const ColumnNature) (PresolveError || std.mem.Allocator.Error)!PresolvedProblem {
    std.debug.assert(nature.len == problem.num_cols);
    return presolveImpl(allocator, problem, nature, .{});
}

/// MIP-safe variant with explicit rule switches for ablation benchmarks.
pub fn presolveMipWithOptions(
    allocator: std.mem.Allocator,
    problem: *const ProblemView,
    nature: []const ColumnNature,
    options: PresolveOptions,
) (PresolveError || std.mem.Allocator.Error)!PresolvedProblem {
    std.debug.assert(nature.len == problem.num_cols);
    return presolveImpl(allocator, problem, nature, options);
}

fn presolveImpl(
    allocator: std.mem.Allocator,
    problem: *const ProblemView,
    nature: ?[]const ColumnNature,
    options: PresolveOptions,
) (PresolveError || std.mem.Allocator.Error)!PresolvedProblem {
    // Select at most one free zero-cost singleton column for each row. If a
    // row has several, later passes see the remaining columns as empty.
    const raw_row_nnz = try allocator.alloc(u32, problem.num_rows);
    defer allocator.free(raw_row_nnz);
    @memset(raw_row_nnz, 0);
    for (problem.matrix.row_indices) |row| raw_row_nnz[row.toUsize()] += 1;

    // Row-activity propagation is intentionally committing-only: it returns
    // values for columns whose implied interval collapses, never looser output
    // bounds. The exact fixed-column substitution below therefore remains the
    // sole postsolve record and dual reconstruction is unchanged.
    const implied_fixings: ?[]?f64 = if (options.implied_fixing)
        try computeImpliedFixings(allocator, problem, nature)
    else
        null;
    defer if (implied_fixings) |values| allocator.free(values);

    const free_singleton_by_row = try allocator.alloc(u32, problem.num_rows);
    defer allocator.free(free_singleton_by_row);
    @memset(free_singleton_by_row, SENTINEL);
    const free_singleton_row_by_col = try allocator.alloc(u32, problem.num_cols);
    defer allocator.free(free_singleton_row_by_col);
    @memset(free_singleton_row_by_col, SENTINEL);
    var free_singleton_count: u32 = 0;
    var free_singleton_row_nnz: usize = 0;
    for (0..problem.num_cols) |column| {
        if (nature) |kinds| {
            // A free singleton reconstructs its value from the row activity,
            // which may be fractional; only continuous columns are MIP-safe.
            if (kinds[column] != .continuous) continue;
        }
        const begin = problem.matrix.col_starts[column];
        const end = problem.matrix.col_starts[column + 1];
        if (!rules.column_substitution.isFreeSingleton(
            end - begin,
            problem.col_cost[column],
            problem.col_lower[column],
            problem.col_upper[column],
        )) continue;
        const row = problem.matrix.row_indices[begin].toUsize();
        if (free_singleton_by_row[row] != SENTINEL) continue;
        free_singleton_by_row[row] = @intCast(column);
        free_singleton_row_by_col[column] = @intCast(row);
        free_singleton_count += 1;
        free_singleton_row_nnz += raw_row_nnz[row] - 1;
    }

    // ── Count substitutable columns and their nonzeros ──
    var substitution_count: u32 = 0;
    var substitution_nnz: usize = 0;
    const objective_sign: f64 = @floatFromInt(@intFromEnum(problem.objective_sense));
    for (0..problem.num_cols) |j| {
        const begin = problem.matrix.col_starts[j];
        const end = problem.matrix.col_starts[j + 1];
        if (free_singleton_row_by_col[j] != SENTINEL) continue;
        switch (classifyColumn(
            implied_fixings,
            j,
            end - begin,
            problem.col_cost[j],
            problem.col_lower[j],
            problem.col_upper[j],
            objective_sign,
            if (nature) |kinds| kinds[j] else null,
        )) {
            .keep => {},
            .substitute => {
                substitution_count += 1;
                substitution_nnz += end - begin;
            },
        }
    }
    const reduced_cols = problem.num_cols - substitution_count - free_singleton_count;

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
    var implied_fixed_count: u32 = 0;
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
        .free_singleton_row_columns = &.{},
        .free_singleton_row_values = &.{},
        .free_singleton_columns = &.{},
        .free_singleton_column_count = free_singleton_count,
        .singleton_rows = &.{},
        .singleton_row_count = 0,
        .redundant_row_count = 0,
        .duplicate_row_count = 0,
        .implied_fixed_column_count = 0,
        .parallel_row_count = 0,
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
    pp.free_singleton_row_columns = try allocator.alloc(u32, free_singleton_row_nnz);
    pp.free_singleton_row_values = try allocator.alloc(f64, free_singleton_row_nnz);
    pp.free_singleton_columns = try allocator.alloc(FreeSingletonColumnRecord, free_singleton_count);
    pp.singleton_rows = try allocator.alloc(SingletonRowRecord, problem.num_rows);

    @memset(pp.original_to_reduced_col, SENTINEL);

    const free_singleton_record_by_row = try allocator.alloc(u32, problem.num_rows);
    defer allocator.free(free_singleton_record_by_row);
    @memset(free_singleton_record_by_row, SENTINEL);
    var free_record_index: u32 = 0;
    var free_entry_offset: u32 = 0;
    for (free_singleton_by_row, 0..) |column, row| {
        if (column == SENTINEL) continue;
        const begin = problem.matrix.col_starts[column];
        const side: FreeSingletonColumnSide = if (std.math.isFinite(problem.row_lower[row]) and
            problem.row_lower[row] == problem.row_upper[row])
            .fixed
        else if (std.math.isFinite(problem.row_lower[row]))
            .lower
        else if (std.math.isFinite(problem.row_upper[row]))
            .upper
        else
            .free;
        const entry_count = raw_row_nnz[row] - 1;
        pp.free_singleton_columns[free_record_index] = .{
            .original_column = column,
            .original_row = @intCast(row),
            .coefficient = problem.matrix.values[begin],
            .side = side,
            .target_activity = switch (side) {
                .fixed => problem.row_lower[row],
                .lower => problem.row_lower[row],
                .upper => problem.row_upper[row],
                .free => 0.0,
            },
            .row_start = free_entry_offset,
            .row_end = free_entry_offset + entry_count,
        };
        free_singleton_record_by_row[row] = free_record_index;
        free_entry_offset += entry_count;
        free_record_index += 1;
    }
    std.debug.assert(free_record_index == free_singleton_count);
    std.debug.assert(free_entry_offset == free_singleton_row_nnz);

    const free_entry_cursor = try allocator.alloc(u32, free_singleton_count);
    defer allocator.free(free_entry_cursor);
    for (pp.free_singleton_columns, 0..) |record, index|
        free_entry_cursor[index] = record.row_start;
    for (0..problem.num_cols) |column| {
        const begin = problem.matrix.col_starts[column];
        const end = problem.matrix.col_starts[column + 1];
        for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row_id, value| {
            const record_index = free_singleton_record_by_row[row_id.toUsize()];
            if (record_index == SENTINEL) continue;
            const record = pp.free_singleton_columns[record_index];
            if (record.original_column == column) continue;
            const cursor = free_entry_cursor[record_index];
            pp.free_singleton_row_columns[cursor] = @intCast(column);
            pp.free_singleton_row_values[cursor] = value;
            free_entry_cursor[record_index] += 1;
        }
    }
    for (pp.free_singleton_columns, free_entry_cursor) |record, cursor|
        std.debug.assert(cursor == record.row_end);

    // ── Substitute columns ──
    var rc: u32 = 0;
    var fi: u32 = 0;
    var fnz: u32 = 0;

    for (0..problem.num_cols) |j| {
        const begin = problem.matrix.col_starts[j];
        const end = problem.matrix.col_starts[j + 1];
        if (free_singleton_row_by_col[j] != SENTINEL) continue;
        const implied_fixed = if (implied_fixings) |values| values[j] != null else false;
        switch (classifyColumn(
            implied_fixings,
            j,
            end - begin,
            problem.col_cost[j],
            problem.col_lower[j],
            problem.col_upper[j],
            objective_sign,
            if (nature) |kinds| kinds[j] else null,
        )) {
            .substitute => |substitution| {
                if (implied_fixed) implied_fixed_count += 1;
                pp.column_substitutions[fi] = .{
                    .original_column = @intCast(j),
                    .value = substitution.value,
                    .col_cost = problem.col_cost[j],
                    .basis_status = switch (substitution.bound) {
                        .fixed => .fixed,
                        .lower => .at_lower,
                        .upper => .at_upper,
                        .free => .free,
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
    pp.implied_fixed_column_count = implied_fixed_count;

    // Remove feasible empty rows and transfer feasible singleton rows into the
    // retained column bounds. Infeasible rows remain for simplex to certify.
    const row_removed = try allocator.alloc(bool, problem.num_rows);
    defer allocator.free(row_removed);
    @memset(row_removed, false);
    for (pp.free_singleton_columns) |record| row_removed[record.original_row] = true;
    var singleton_count: u32 = 0;
    for (0..problem.num_rows) |row| {
        if (row_removed[row]) continue;
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

    // Prove general retained rows redundant from an outward-rounded interval
    // of every possible row activity under the current column bounds. Bounds
    // transferred by singleton rows above participate immediately.
    const row_activity = try allocator.alloc(rules.redundant_row.ActivityBounds, problem.num_rows);
    defer allocator.free(row_activity);
    for (row_activity) |*activity| activity.* = .{};
    for (0..problem.num_cols) |original_column| {
        const reduced_column = pp.original_to_reduced_col[original_column];
        if (reduced_column == SENTINEL) continue;
        const begin = problem.matrix.col_starts[original_column];
        const end = problem.matrix.col_starts[original_column + 1];
        for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row_id, coefficient| {
            const row = row_id.toUsize();
            if (row_removed[row]) continue;
            row_activity[row].add(coefficient, pp.col_lower[reduced_column], pp.col_upper[reduced_column]);
        }
    }
    var redundant_count: u32 = 0;
    for (0..problem.num_rows) |row| {
        if (row_removed[row]) continue;
        if (rules.redundant_row.classify(
            row_activity[row],
            adjusted_row_lower[row],
            adjusted_row_upper[row],
        ) == .remove) {
            row_removed[row] = true;
            redundant_count += 1;
        }
    }
    pp.redundant_row_count = redundant_count;

    // Materialize retained rows once in column order. Hashes only select
    // candidates; exact comparison makes collisions unable to remove a row.
    const retained_row_starts = try allocator.alloc(usize, problem.num_rows + 1);
    defer allocator.free(retained_row_starts);
    retained_row_starts[0] = 0;
    for (retained_row_nnz, 0..) |count, row|
        retained_row_starts[row + 1] = retained_row_starts[row] + count;
    const retained_entry_count = retained_row_starts[problem.num_rows];
    const retained_row_columns = try allocator.alloc(u32, retained_entry_count);
    defer allocator.free(retained_row_columns);
    const retained_row_values = try allocator.alloc(f64, retained_entry_count);
    defer allocator.free(retained_row_values);
    const retained_row_cursor = try allocator.dupe(usize, retained_row_starts[0..problem.num_rows]);
    defer allocator.free(retained_row_cursor);
    for (0..problem.num_cols) |original_column| {
        if (pp.original_to_reduced_col[original_column] == SENTINEL) continue;
        const begin = problem.matrix.col_starts[original_column];
        const end = problem.matrix.col_starts[original_column + 1];
        for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row_id, value| {
            const row = row_id.toUsize();
            const cursor = retained_row_cursor[row];
            retained_row_columns[cursor] = @intCast(original_column);
            retained_row_values[cursor] = value;
            retained_row_cursor[row] += 1;
        }
    }

    var duplicate_heads: std.AutoHashMapUnmanaged(u64, u32) = .empty;
    defer duplicate_heads.deinit(allocator);
    const duplicate_next = try allocator.alloc(u32, problem.num_rows);
    defer allocator.free(duplicate_next);
    @memset(duplicate_next, SENTINEL);
    var duplicate_count: u32 = 0;
    for (0..problem.num_rows) |row| {
        if (row_removed[row]) continue;
        const begin = retained_row_starts[row];
        const end = retained_row_starts[row + 1];
        const columns = retained_row_columns[begin..end];
        const values = retained_row_values[begin..end];
        const signature = rules.duplicate_row.hash(
            adjusted_row_lower[row],
            adjusted_row_upper[row],
            columns,
            values,
        );
        var candidate = duplicate_heads.get(signature) orelse SENTINEL;
        while (candidate != SENTINEL) : (candidate = duplicate_next[candidate]) {
            const candidate_begin = retained_row_starts[candidate];
            const candidate_end = retained_row_starts[candidate + 1];
            if (!rules.duplicate_row.exactlyEqual(
                adjusted_row_lower[row],
                adjusted_row_upper[row],
                columns,
                values,
                adjusted_row_lower[candidate],
                adjusted_row_upper[candidate],
                retained_row_columns[candidate_begin..candidate_end],
                retained_row_values[candidate_begin..candidate_end],
            )) continue;
            row_removed[row] = true;
            duplicate_count += 1;
            break;
        }
        if (row_removed[row]) continue;
        duplicate_next[row] = duplicate_heads.get(signature) orelse SENTINEL;
        try duplicate_heads.put(allocator, signature, @intCast(row));
    }
    pp.duplicate_row_count = duplicate_count;

    // Remove rows that are exact non-zero scalar multiples of a retained row.
    // A scaled copy defines the same half-spaces, so its dual multiplier can be
    // zero during postsolve. Approximate proportionality is deliberately not
    // accepted: every coefficient and both bounds must be bit-exact after
    // scaling, which keeps the reduction certificate exact.
    var parallel_count: u32 = 0;
    if (options.parallel_rows) {
        var parallel_heads: std.AutoHashMapUnmanaged(u64, u32) = .empty;
        defer parallel_heads.deinit(allocator);
        const parallel_next = try allocator.alloc(u32, problem.num_rows);
        defer allocator.free(parallel_next);
        @memset(parallel_next, SENTINEL);
        for (0..problem.num_rows) |row| {
            if (row_removed[row]) continue;
            const begin = retained_row_starts[row];
            const end = retained_row_starts[row + 1];
            const columns = retained_row_columns[begin..end];
            const signature = rules.duplicate_row.patternHash(columns);
            var candidate = parallel_heads.get(signature) orelse SENTINEL;
            while (candidate != SENTINEL) : (candidate = parallel_next[candidate]) {
                const candidate_begin = retained_row_starts[candidate];
                const candidate_end = retained_row_starts[candidate + 1];
                const candidate_columns = retained_row_columns[candidate_begin..candidate_end];
                if (!std.mem.eql(u32, columns, candidate_columns)) continue;
                const scale = rules.duplicate_row.parallelScale(
                    candidate_columns,
                    retained_row_values[candidate_begin..candidate_end],
                    columns,
                    retained_row_values[begin..end],
                ) orelse continue;
                if (!rules.duplicate_row.boundsMatchScale(
                    scale,
                    adjusted_row_lower[candidate],
                    adjusted_row_upper[candidate],
                    adjusted_row_lower[row],
                    adjusted_row_upper[row],
                )) continue;
                row_removed[row] = true;
                parallel_count += 1;
                break;
            }
            if (row_removed[row]) continue;
            parallel_next[row] = parallel_heads.get(signature) orelse SENTINEL;
            try parallel_heads.put(allocator, signature, @intCast(row));
        }
    }
    pp.parallel_row_count = parallel_count;

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
    for (pp.free_singleton_columns[0..pp.free_singleton_column_count]) |rec| {
        var other_activity: f64 = 0.0;
        for (rec.row_start..rec.row_end) |entry| {
            other_activity += pp.free_singleton_row_values[entry] *
                result.primal[pp.free_singleton_row_columns[entry]];
        }
        result.primal[rec.original_column] = (rec.target_activity - other_activity) / rec.coefficient;
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
        for (pp.free_singleton_columns[0..pp.free_singleton_column_count]) |rec| {
            var other_direction: f64 = 0.0;
            for (rec.row_start..rec.row_end) |entry| {
                other_direction += pp.free_singleton_row_values[entry] *
                    result.unbounded_ray[pp.free_singleton_row_columns[entry]];
            }
            result.unbounded_ray[rec.original_column] = -other_direction / rec.coefficient;
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
    var pp = try presolve(a, &problem);
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
    var pp = try presolve(a, &problem);
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

test "singleton bounds immediately expose a general redundant row" {
    const problem = ProblemView{
        .num_rows = 2,
        .num_cols = 2,
        .col_cost = &.{ 1.0, 2.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ 10.0, 1.0 },
        .row_lower = &.{ -std.math.inf(f64), 0.0 },
        .row_upper = &.{ 1.0, 2.0 },
        .matrix = matrix.CscView.initAssumeValid(
            2,
            2,
            &[_]usize{ 0, 2, 3 },
            &.{
                foundation.RowId.fromUsizeAssumeValid(0),
                foundation.RowId.fromUsizeAssumeValid(1),
                foundation.RowId.fromUsizeAssumeValid(1),
            },
            &.{ 1.0, 1.0, 1.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var pp = try presolve(std.testing.allocator, &problem);
    defer pp.deinit();
    try std.testing.expectEqual(@as(u32, 1), pp.singleton_row_count);
    try std.testing.expectEqual(@as(u32, 1), pp.redundant_row_count);
    try std.testing.expectEqual(@as(u32, 0), pp.num_rows);
    try std.testing.expectEqual(@as(f64, 1.0), pp.col_upper[0]);
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
    var pp = try presolve(a, &problem);
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
    var pp = try presolve(a, &problem);
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
    var pp = try presolve(a, &problem);
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
    var pp = try presolve(a, &problem);
    defer pp.deinit();
    try std.testing.expect(!pp.wasApplied());
}

test "mip mode keeps a discrete column fixed at a fractional value" {
    const a = std.testing.allocator;
    // An integer variable fixed at 2.5 makes the model infeasible; presolve
    // must not substitute it away and hide the infeasibility.
    const problem = ProblemView{
        .num_rows = 0,
        .num_cols = 1,
        .col_cost = &.{1.0},
        .col_lower = &.{2.5},
        .col_upper = &.{2.5},
        .row_lower = &.{},
        .row_upper = &.{},
        .matrix = matrix.CscView.initAssumeValid(0, 1, &[_]usize{ 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var pp = try presolveMip(a, &problem, &.{.discrete});
    defer pp.deinit();
    try std.testing.expect(!pp.wasApplied());
    try std.testing.expectEqual(@as(u32, 1), pp.num_cols);

    var lp_pp = try presolve(a, &problem);
    defer lp_pp.deinit();
    try std.testing.expectEqual(@as(u32, 1), lp_pp.column_substitution_count);
}

test "mip mode substitutes discrete fixed and rounded empty columns" {
    const a = std.testing.allocator;
    // Column 0: integer fixed at 3 → substituted exactly.
    // Column 1: integer, empty, positive cost, bounds [0.3, 5] → rounded up to 1.
    // Column 2: continuous, empty, positive cost, bounds [0.3, 5] → kept at 0.3.
    const problem = ProblemView{
        .num_rows = 0,
        .num_cols = 3,
        .col_cost = &.{ 1.0, 2.0, 2.0 },
        .col_lower = &.{ 3.0, 0.3, 0.3 },
        .col_upper = &.{ 3.0, 5.0, 5.0 },
        .row_lower = &.{},
        .row_upper = &.{},
        .matrix = matrix.CscView.initAssumeValid(0, 3, &[_]usize{ 0, 0, 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var pp = try presolveMip(a, &problem, &.{ .discrete, .discrete, .continuous });
    defer pp.deinit();
    try std.testing.expectEqual(@as(u32, 3), pp.column_substitution_count);
    try std.testing.expectEqual(@as(u32, 0), pp.num_cols);
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), pp.column_substitutions[0].value, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), pp.column_substitutions[1].value, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), pp.column_substitutions[2].value, 1e-12);
    // offset = 1*3 + 2*1 + 2*0.3
    try std.testing.expectApproxEqAbs(@as(f64, 5.6), pp.objective_offset, 1e-12);
}

test "mip mode skips free singleton elimination for discrete columns only" {
    const a = std.testing.allocator;
    // Row: x + y = 4 where x is a free zero-cost singleton and y is boxed.
    const row_id = foundation.RowId.fromUsizeAssumeValid(0);
    const problem = ProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 0.0, 1.0 },
        .col_lower = &.{ -std.math.inf(f64), 0.0 },
        .col_upper = &.{ std.math.inf(f64), 10.0 },
        .row_lower = &.{4.0},
        .row_upper = &.{4.0},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &[_]usize{ 0, 1, 2 }, &.{ row_id, row_id }, &.{ 1.0, 1.0 }),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    // Discrete singleton: reconstruction could be fractional, so keep both.
    var discrete_pp = try presolveMip(a, &problem, &.{ .discrete, .discrete });
    defer discrete_pp.deinit();
    try std.testing.expectEqual(@as(u32, 0), discrete_pp.free_singleton_column_count);

    // Continuous singleton: eliminated and reconstructed exactly.
    var continuous_pp = try presolveMip(a, &problem, &.{ .continuous, .discrete });
    defer continuous_pp.deinit();
    try std.testing.expectEqual(@as(u32, 1), continuous_pp.free_singleton_column_count);
    try std.testing.expectEqual(@as(u32, 0), continuous_pp.num_rows);
    const sol = SolutionView{
        .status = .optimal,
        .objective_value = 0.0,
        .iterations = 0,
        .primal = &.{0.0},
        .dual = &.{},
        .reduced_cost = &.{1.0},
        .unbounded_ray = &.{},
        .infeasibility_ray = &.{},
    };
    var res = try postsolve(&continuous_pp, sol, a);
    defer res.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), res.primal[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), res.primal[1], 1e-12);
}

test "mip pipeline threads column nature through passes" {
    const a = std.testing.allocator;
    // x fixed at 2 (discrete) exposes y's row as a singleton: y >= 1.
    // Pass 2 then sees y as a discrete column with tightened bounds.
    const problem = ProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 1.0, 1.0 },
        .col_lower = &.{ 2.0, 0.0 },
        .col_upper = &.{ 2.0, 5.0 },
        .row_lower = &.{3.0},
        .row_upper = &.{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &[_]usize{ 0, 1, 2 }, &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) }, &.{ 1.0, 1.0 }),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var pipeline = try PresolvePipeline.initMip(a, &problem, &.{ .discrete, .discrete });
    defer pipeline.deinit();
    try std.testing.expect(pipeline.wasApplied());
    // Pass 1 substitutes x (fixed at 2) and transfers the singleton row into
    // y >= 1; pass 2 substitutes the now-empty y at its rounded lower bound 1.
    try std.testing.expectEqual(@as(u32, 2), pipeline.column_substitution_count);
    const reduced = pipeline.reducedView();
    try std.testing.expectEqual(@as(u32, 0), reduced.num_rows);
    try std.testing.expectEqual(@as(u32, 0), reduced.num_cols);
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), reduced.objective_offset, 1e-12);
}



test "row activity forcing fixes both continuous columns" {
    const a = std.testing.allocator;
    // min x + y  s.t. x + y = 20, 0 <= x,y <= 10.
    const problem = ProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 1.0, 1.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ 10.0, 10.0 },
        .row_lower = &.{20.0},
        .row_upper = &.{20.0},
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
    var pp = try presolve(a, &problem);
    defer pp.deinit();
    try std.testing.expectEqual(@as(u32, 2), pp.implied_fixed_column_count);
    try std.testing.expectEqual(@as(u32, 0), pp.num_cols);
    try std.testing.expectEqual(@as(u32, 0), pp.num_rows);
    try std.testing.expectApproxEqAbs(@as(f64, 20.0), pp.objective_offset, 1e-12);

    const sol = SolutionView{
        .status = .optimal,
        .objective_value = 20.0,
        .iterations = 0,
        .primal = &.{},
        .dual = &.{},
        .reduced_cost = &.{},
        .unbounded_ray = &.{},
        .infeasibility_ray = &.{},
    };
    var res = try postsolve(&pp, sol, a);
    defer res.deinit();
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), res.primal[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), res.primal[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 20.0), res.objective_value, 1e-12);
}

test "mip row activity fixes a forced integer pair" {
    const a = std.testing.allocator;
    // x + y = 7 with integer x in [6.5,10], y in [0,10]: integer rounding
    // gives x = 7 and then y = 0.
    const problem = ProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 1.0, 1.0 },
        .col_lower = &.{ 6.5, 0.0 },
        .col_upper = &.{ 10.0, 10.0 },
        .row_lower = &.{7.0},
        .row_upper = &.{7.0},
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
    var pp = try presolveMip(a, &problem, &.{ .discrete, .discrete });
    defer pp.deinit();
    try std.testing.expectEqual(@as(u32, 2), pp.implied_fixed_column_count);
    try std.testing.expectEqual(@as(u32, 0), pp.num_cols);
}

test "exact scaled rows are removed and weaker options keep them" {
    const a = std.testing.allocator;
    // r1: x + 2y <= 5; r2 = 2*r1; r3 = -1*r1. Every row defines the same
    // half-space, so two copies can be removed with zero postsolve duals.
    const problem = ProblemView{
        .num_rows = 3,
        .num_cols = 2,
        .col_cost = &.{ 1.0, 1.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ 10.0, 10.0 },
        .row_lower = &.{ -std.math.inf(f64), -std.math.inf(f64), -5.0 },
        .row_upper = &.{ 5.0, 10.0, std.math.inf(f64) },
        .matrix = matrix.CscView.initAssumeValid(
            3,
            2,
            &[_]usize{ 0, 3, 6 },
            &.{
                foundation.RowId.fromUsizeAssumeValid(0),
                foundation.RowId.fromUsizeAssumeValid(1),
                foundation.RowId.fromUsizeAssumeValid(2),
                foundation.RowId.fromUsizeAssumeValid(0),
                foundation.RowId.fromUsizeAssumeValid(1),
                foundation.RowId.fromUsizeAssumeValid(2),
            },
            &.{ 1.0, 2.0, -1.0, 2.0, 4.0, -2.0 },
        ),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var pp = try presolve(a, &problem);
    defer pp.deinit();
    try std.testing.expectEqual(@as(u32, 1), pp.num_rows);
    try std.testing.expectEqual(@as(u32, 2), pp.parallel_row_count);

    var disabled = try presolveWithOptions(a, &problem, .{ .parallel_rows = false });
    defer disabled.deinit();
    try std.testing.expectEqual(@as(u32, 3), disabled.num_rows);
    try std.testing.expectEqual(@as(u32, 0), disabled.parallel_row_count);
}
