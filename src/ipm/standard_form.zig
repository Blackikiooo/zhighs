//! Transformation of a borrowed `ProblemView` into IPM standard form
//! min ĉᵀx̂ subject to Â·x̂ = b̂, x̂ ≥ 0.
//!
//! The conversion shifts finite lower bounds into the right-hand side, splits
//! free columns into a positive/negative difference, negates upper-only
//! columns, and turns column upper bounds and row bounds into non-negative
//! slack columns with unit entries. Fixed columns are folded into the
//! right-hand side and the objective shift. The engine minimizes with the
//! objective sense baked into ĉ, so the sense only reappears when the result
//! is mapped back to model coordinates.

const std = @import("std");
const Allocator = std.mem.Allocator;
const lp = @import("lp");

const inf = std.math.inf(f64);

pub const BuildError = error{InvalidModel} || Allocator.Error;

/// How one original column appears in standard form.
pub const ColKind = enum {
    /// l == u: substituted out, contributes only to rhs and objective shift.
    fixed,
    /// finite l: x = l + x̂ with x̂ ≥ 0; a finite u adds a bound row.
    shifted,
    /// l = −∞, finite u: x = u − x̂ with x̂ ≥ 0 and negated coefficients.
    flipped,
    /// l = −∞, u = +∞: x = x̂⁺ − x̂⁻ with both parts ≥ 0.
    free_split,
};

/// Per-original-column mapping record used to translate results back.
pub const ColRecord = struct {
    kind: ColKind,
    /// Constant part of x (lower bound, upper bound for flipped, 0 for free).
    base: f64,
    /// Standard-form column of the primary part; invalid for fixed columns.
    primary: usize,
    /// Standard-form column of the negative part for free splits.
    negative: ?usize,
};

/// Standard-form equality row origin: `null` marks an auxiliary bound row
/// introduced for a column upper bound or a range-row slack.
pub const RowOrigin = ?usize;

/// Owning standard-form model plus the records required to map back.
pub const StandardForm = struct {
    allocator: Allocator,
    /// Number of equality rows in Â.
    num_rows: usize,
    /// Number of non-negative standard-form columns.
    num_cols: usize,
    /// CSC storage for Â (row indices are plain usize, canonical order).
    col_starts: []usize,
    row_indices: []usize,
    values: []f64,
    /// Standard-form objective (minimization, sense already applied).
    cost: []f64,
    /// Standard-form right-hand side b̂.
    rhs: []f64,
    /// One record per original structural column.
    col_records: []ColRecord,
    /// Original row index for each standard-form equality row.
    row_origin: []RowOrigin,
    /// Standard-form row index for each original row; `null` for free rows.
    row_to_std: []?usize,
    /// Objective contribution of substituted bound constants (min convention).
    objective_shift: f64,

    pub fn deinit(self: *StandardForm) void {
        const a = self.allocator;
        a.free(self.col_starts);
        a.free(self.row_indices);
        a.free(self.values);
        a.free(self.cost);
        a.free(self.rhs);
        a.free(self.col_records);
        a.free(self.row_origin);
        a.free(self.row_to_std);
    }

    /// Objective value in the original sense and coordinates, including the
    /// model offset, from the standard-form objective `costᵀx̂`.
    pub fn originalObjective(self: *const StandardForm, problem: lp.ProblemView, std_objective: f64) f64 {
        const sense = senseSign(problem.objective_sense);
        return sense * (std_objective + self.objective_shift) + problem.objective_offset;
    }

    /// Min-convention reduced cost of original column j from min-convention
    /// row multipliers `row_dual`: c̃_j − a_jᵀy.
    pub fn columnReducedCostMin(self: *const StandardForm, problem: lp.ProblemView, column: usize, row_dual: []const f64) f64 {
        _ = self;
        const sense = senseSign(problem.objective_sense);
        var rc = sense * problem.col_cost[column];
        const m = problem.matrix;
        for (m.col_starts[column]..m.col_starts[column + 1]) |k| {
            rc -= m.values[k] * row_dual[m.row_indices[k].toUsize()];
        }
        return rc;
    }
};

fn senseSign(sense: lp.ObjectiveSense) f64 {
    return @floatFromInt(@intFromEnum(sense));
}

const ColClass = struct {
    kind: ColKind,
    base: f64,
    has_bound_row: bool,
};

const RowClass = enum { free, equality, lower, upper, range };

fn classifyColumn(lower: f64, upper: f64) ColClass {
    if (lower == upper) return .{ .kind = .fixed, .base = lower, .has_bound_row = false };
    if (lower == -inf) {
        if (upper == inf) return .{ .kind = .free_split, .base = 0.0, .has_bound_row = false };
        return .{ .kind = .flipped, .base = upper, .has_bound_row = false };
    }
    return .{ .kind = .shifted, .base = lower, .has_bound_row = upper != inf };
}

fn classifyRow(lower: f64, upper: f64) RowClass {
    if (lower == -inf and upper == inf) return .free;
    if (lower == upper) return .equality;
    if (upper == inf) return .lower;
    if (lower == -inf) return .upper;
    return .range;
}

/// Build the owning standard-form model from a validated borrowed view.
pub fn build(allocator: Allocator, problem: lp.ProblemView) BuildError!StandardForm {
    problem.validate() catch return error.InvalidModel;

    const m_in = problem.num_rows;
    const n_in = problem.num_cols;
    const sense = senseSign(problem.objective_sense);
    const csc = problem.matrix;

    const col_class = try allocator.alloc(ColClass, n_in);
    defer allocator.free(col_class);
    const row_class = try allocator.alloc(RowClass, m_in);
    defer allocator.free(row_class);
    const row_to_std = try allocator.alloc(?usize, m_in);
    errdefer allocator.free(row_to_std);

    // Structural-row constants from substituted bound offsets: for x = base + d·x̂
    // the term a·base moves to the right-hand side of every row touching a.
    const row_const = try allocator.alloc(f64, m_in);
    defer allocator.free(row_const);
    @memset(row_const, 0.0);
    var objective_shift: f64 = 0.0;
    for (0..n_in) |j| {
        col_class[j] = classifyColumn(problem.col_lower[j], problem.col_upper[j]);
        const base = col_class[j].base;
        if (base != 0.0) {
            for (csc.col_starts[j]..csc.col_starts[j + 1]) |k| {
                row_const[csc.row_indices[k].toUsize()] += csc.values[k] * base;
            }
            objective_shift += sense * problem.col_cost[j] * base;
        }
    }

    // Assign standard-form row indices: original non-free rows first, then one
    // bound row per upper-bounded shifted column and per range-row slack.
    var m_std: usize = 0;
    for (0..m_in) |i| {
        row_class[i] = classifyRow(problem.row_lower[i], problem.row_upper[i]);
        if (row_class[i] == .free) {
            row_to_std[i] = null;
        } else {
            row_to_std[i] = m_std;
            m_std += 1;
        }
    }
    const bound_row_base = m_std;
    for (0..n_in) |j| {
        if (col_class[j].has_bound_row) m_std += 1;
    }
    const range_row_base = m_std;
    for (0..m_in) |i| {
        if (row_class[i] == .range) m_std += 1;
    }

    // Column counts per original column / slack.
    var n_std: usize = 0;
    for (0..n_in) |j| {
        n_std += switch (col_class[j].kind) {
            .fixed => 0,
            .shifted, .flipped => 1,
            .free_split => 2,
        };
    }
    for (0..m_in) |i| {
        n_std += switch (row_class[i]) {
            .free, .equality => 0,
            .lower, .upper => 1,
            .range => 2,
        };
    }
    // One slack column per column-bound row (x̂ + t = u − l).
    for (0..n_in) |j| {
        if (col_class[j].has_bound_row) n_std += 1;
    }

    const col_starts = try allocator.alloc(usize, n_std + 1);
    errdefer allocator.free(col_starts);
    const col_records = try allocator.alloc(ColRecord, n_in);
    errdefer allocator.free(col_records);
    const row_origin = try allocator.alloc(RowOrigin, m_std);
    errdefer allocator.free(row_origin);
    const cost = try allocator.alloc(f64, n_std);
    errdefer allocator.free(cost);
    const rhs = try allocator.alloc(f64, m_std);
    errdefer allocator.free(rhs);

    // Row origins and right-hand sides.
    for (0..m_in) |i| {
        if (row_to_std[i]) |r| {
            row_origin[r] = i;
            const bound = switch (row_class[i]) {
                .equality, .lower => problem.row_lower[i],
                .upper, .range => problem.row_upper[i],
                .free => unreachable,
            };
            rhs[r] = bound - row_const[i];
        }
    }
    var bound_cursor = bound_row_base;
    for (0..n_in) |j| {
        if (!col_class[j].has_bound_row) continue;
        row_origin[bound_cursor] = null;
        rhs[bound_cursor] = problem.col_upper[j] - col_class[j].base;
        bound_cursor += 1;
    }
    var range_cursor = range_row_base;
    for (0..m_in) |i| {
        if (row_class[i] != .range) continue;
        row_origin[range_cursor] = null;
        rhs[range_cursor] = problem.row_upper[i] - problem.row_lower[i];
        range_cursor += 1;
    }

    // Pass 1: column lengths (CSC entries per standard-form column).
    const lengths = col_starts[0..n_std];
    var cursor: usize = 0;
    bound_cursor = bound_row_base;
    range_cursor = range_row_base;
    for (0..n_in) |j| {
        const class = col_class[j];
        const entries = blk: {
            var count: usize = 0;
            for (csc.col_starts[j]..csc.col_starts[j + 1]) |k| {
                if (row_to_std[csc.row_indices[k].toUsize()] != null) count += 1;
            }
            break :blk count;
        };
        switch (class.kind) {
            .fixed => col_records[j] = .{ .kind = .fixed, .base = class.base, .primary = 0, .negative = null },
            .shifted, .flipped => {
                col_records[j] = .{ .kind = class.kind, .base = class.base, .primary = cursor, .negative = null };
                lengths[cursor] = entries + @intFromBool(class.has_bound_row);
                const dir: f64 = if (class.kind == .flipped) -1.0 else 1.0;
                cost[cursor] = sense * problem.col_cost[j] * dir;
                cursor += 1;
                if (class.has_bound_row) bound_cursor += 1;
            },
            .free_split => {
                col_records[j] = .{ .kind = .free_split, .base = 0.0, .primary = cursor, .negative = cursor + 1 };
                lengths[cursor] = entries;
                lengths[cursor + 1] = entries;
                cost[cursor] = sense * problem.col_cost[j];
                cost[cursor + 1] = -sense * problem.col_cost[j];
                cursor += 2;
            },
        }
    }
    for (0..m_in) |i| {
        switch (row_class[i]) {
            .free, .equality => {},
            .lower, .upper => {
                lengths[cursor] = 1;
                cost[cursor] = 0.0;
                cursor += 1;
            },
            .range => {
                lengths[cursor] = 2;
                cost[cursor] = 0.0;
                lengths[cursor + 1] = 1;
                cost[cursor + 1] = 0.0;
                cursor += 2;
                range_cursor += 1;
            },
        }
    }
    // Column-bound slacks trail the row slacks, one per bounded column.
    for (0..n_in) |j| {
        if (!col_class[j].has_bound_row) continue;
        lengths[cursor] = 1;
        cost[cursor] = 0.0;
        cursor += 1;
        bound_cursor += 1;
    }
    std.debug.assert(cursor == n_std);

    var start: usize = 0;
    for (0..n_std) |j| {
        // `lengths` aliases col_starts[0..n_std]; read before overwriting.
        const len = lengths[j];
        col_starts[j] = start;
        start += len;
    }
    col_starts[n_std] = start;
    const nnz = start;
    const row_indices = try allocator.alloc(usize, nnz);
    errdefer allocator.free(row_indices);
    const values = try allocator.alloc(f64, nnz);
    errdefer allocator.free(values);

    // Pass 2: fill CSC payload. `next` walks the write cursor per column.
    const next = try allocator.alloc(usize, n_std);
    defer allocator.free(next);
    @memcpy(next, col_starts[0..n_std]);

    bound_cursor = bound_row_base;
    range_cursor = range_row_base;
    for (0..n_in) |j| {
        const class = col_class[j];
        if (class.kind == .fixed) continue;
        const sign: f64 = if (class.kind == .flipped) -1.0 else 1.0;
        const primary = col_records[j].primary;
        for (csc.col_starts[j]..csc.col_starts[j + 1]) |k| {
            const std_row = row_to_std[csc.row_indices[k].toUsize()] orelse continue;
            row_indices[next[primary]] = std_row;
            values[next[primary]] = sign * csc.values[k];
            next[primary] += 1;
            if (class.kind == .free_split) {
                const neg = col_records[j].negative.?;
                row_indices[next[neg]] = std_row;
                values[next[neg]] = -csc.values[k];
                next[neg] += 1;
            }
        }
        if (class.has_bound_row) {
            row_indices[next[primary]] = bound_cursor;
            values[next[primary]] = 1.0;
            next[primary] += 1;
            bound_cursor += 1;
        }
    }
    // Slack columns live after all structural parts; recompute their start.
    var slack_cursor: usize = 0;
    for (0..n_in) |j| {
        slack_cursor += switch (col_class[j].kind) {
            .fixed => 0,
            .shifted, .flipped => 1,
            .free_split => 2,
        };
    }
    for (0..m_in) |i| {
        const r = row_to_std[i] orelse continue;
        switch (row_class[i]) {
            .free, .equality => {},
            .lower => {
                row_indices[next[slack_cursor]] = r;
                values[next[slack_cursor]] = -1.0;
                next[slack_cursor] += 1;
                slack_cursor += 1;
            },
            .upper => {
                row_indices[next[slack_cursor]] = r;
                values[next[slack_cursor]] = 1.0;
                next[slack_cursor] += 1;
                slack_cursor += 1;
            },
            .range => {
                // Ax + s = u with 0 ≤ s ≤ u − l, so s also enters its bound row.
                row_indices[next[slack_cursor]] = r;
                values[next[slack_cursor]] = 1.0;
                next[slack_cursor] += 1;
                row_indices[next[slack_cursor]] = range_cursor;
                values[next[slack_cursor]] = 1.0;
                next[slack_cursor] += 1;
                slack_cursor += 1;
                row_indices[next[slack_cursor]] = range_cursor;
                values[next[slack_cursor]] = 1.0;
                next[slack_cursor] += 1;
                slack_cursor += 1;
                range_cursor += 1;
            },
        }
    }
    // Column-bound slacks trail the row slacks, matching the pass-1 layout.
    bound_cursor = bound_row_base;
    for (0..n_in) |j| {
        if (!col_class[j].has_bound_row) continue;
        row_indices[next[slack_cursor]] = bound_cursor;
        values[next[slack_cursor]] = 1.0;
        next[slack_cursor] += 1;
        slack_cursor += 1;
        bound_cursor += 1;
    }
    std.debug.assert(slack_cursor == n_std);
    std.debug.assert(std.mem.eql(usize, next, col_starts[1..]));

    return .{
        .allocator = allocator,
        .num_rows = m_std,
        .num_cols = n_std,
        .col_starts = col_starts,
        .row_indices = row_indices,
        .values = values,
        .cost = cost,
        .rhs = rhs,
        .col_records = col_records,
        .row_origin = row_origin,
        .row_to_std = row_to_std,
        .objective_shift = objective_shift,
    };
}
