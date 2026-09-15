//! Borrowed LP input view for simplex.
//!
//! `ProblemView` never owns or copies model data. It is deliberately
//! independent of the public model layer so the numerical kernel can sit below
//! solver dispatch and model adapters in the dependency graph.

const std = @import("std");
const foundation = @import("foundation");
const matrix = @import("matrix");

/// Optimization direction. Encoded as +1/-1 so flipping the sense is just a
/// sign change on the objective vector.
pub const ObjectiveSense = enum(i2) { minimize = 1, maximize = -1 };

/// Non-owning view over a sparse LP in CSC form. All slices are borrowed from
/// the caller; the view performs no allocation.
pub const ProblemView = struct {
    /// Number of constraint rows and therefore the dimension of every basis.
    num_rows: usize,
    /// Number of structural columns, excluding logical and artificial columns.
    num_cols: usize,
    /// Objective coefficient for each structural column.
    col_cost: []const f64,
    /// Inclusive lower bound for each structural column.
    col_lower: []const f64,
    /// Inclusive upper bound for each structural column.
    col_upper: []const f64,
    /// Inclusive lower activity bound for each constraint row.
    row_lower: []const f64,
    /// Inclusive upper activity bound for each constraint row.
    row_upper: []const f64,
    /// Borrowed structural constraint matrix in compressed-column form.
    matrix: matrix.CscView,
    /// Direction of optimization before internal minimization normalization.
    objective_sense: ObjectiveSense,
    /// Constant added to the reported objective value.
    objective_offset: f64,

    /// Structural validation failures detected before simplex initialization.
    pub const ViewError = error{ DimensionMismatch, InvalidBounds, InvalidObjective, InvalidMatrix };

    /// Validate borrowed dimensions, numerical values, and canonical CSC
    /// storage without taking ownership or allocating. Intended for the
    /// compile/dispatch boundary, not for repeated simplex iterations.
    pub fn validate(self: ProblemView) ViewError!void {
        try self.validateExceptMatrix();
        self.matrix.validate() catch return error.InvalidMatrix;
    }

    /// Validate everything except the O(nnz) CSC scan. Callers that reuse the
    /// exact matrix slices validated previously (bound-only reoptimizations)
    /// can skip that scan through slice-identity caching.
    pub fn validateExceptMatrix(self: ProblemView) ViewError!void {
        if (self.col_cost.len != self.num_cols or self.col_lower.len != self.num_cols or self.col_upper.len != self.num_cols or
            self.row_lower.len != self.num_rows or self.row_upper.len != self.num_rows)
            return error.DimensionMismatch;
        if (self.matrix.num_rows != self.num_rows or self.matrix.num_cols != self.num_cols)
            return error.InvalidMatrix;
        if (!std.math.isFinite(self.objective_offset)) return error.InvalidObjective;
        for (self.col_cost) |cost| {
            if (!std.math.isFinite(cost)) return error.InvalidObjective;
        }
        for (self.col_lower, self.col_upper) |lower, upper| {
            if (std.math.isNan(lower) or std.math.isNan(upper) or
                lower == std.math.inf(f64) or upper == -std.math.inf(f64) or lower > upper)
                return error.InvalidBounds;
        }
        for (self.row_lower, self.row_upper) |lower, upper| {
            if (std.math.isNan(lower) or std.math.isNan(upper) or
                lower == std.math.inf(f64) or upper == -std.math.inf(f64) or lower > upper)
                return error.InvalidBounds;
        }
    }

    /// Materialize one borrowed CSC column into caller-owned dense storage.
    /// This hot-path helper performs no allocation.
    pub fn fillColumn(self: ProblemView, column: usize, output: []f64) ViewError!void {
        if (column >= self.num_cols or output.len != self.num_rows) return error.DimensionMismatch;
        @memset(output, 0.0);
        const begin = self.matrix.col_starts[column];
        const end = self.matrix.col_starts[column + 1];
        for (self.matrix.row_indices[begin..end], self.matrix.values[begin..end]) |row, value| {
            const index = row.toUsize();
            if (index >= output.len) return error.InvalidMatrix;
            output[index] = value;
        }
    }
};

test "ProblemView validates borrowed dimensions" {
    const view = ProblemView{
        .num_rows = 0,
        .num_cols = 0,
        .col_cost = &.{},
        .col_lower = &.{},
        .col_upper = &.{},
        .row_lower = &.{},
        .row_upper = &.{},
        .matrix = matrix.CscView.initAssumeValid(0, 0, &[_]usize{0}, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    try view.validate();
}

test "ProblemView rejects nonfinite objectives and invalid infinite bounds" {
    var view = ProblemView{
        .num_rows = 0,
        .num_cols = 1,
        .col_cost = &.{0.0},
        .col_lower = &.{0.0},
        .col_upper = &.{1.0},
        .row_lower = &.{},
        .row_upper = &.{},
        .matrix = matrix.CscView.initAssumeValid(0, 1, &[_]usize{ 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };

    view.col_cost = &.{std.math.nan(f64)};
    try std.testing.expectError(error.InvalidObjective, view.validate());
    view.col_cost = &.{0.0};
    view.objective_offset = std.math.inf(f64);
    try std.testing.expectError(error.InvalidObjective, view.validate());
    view.objective_offset = 0.0;
    view.col_lower = &.{std.math.nan(f64)};
    try std.testing.expectError(error.InvalidBounds, view.validate());
    view.col_lower = &.{std.math.inf(f64)};
    try std.testing.expectError(error.InvalidBounds, view.validate());
    view.col_lower = &.{-std.math.inf(f64)};
    view.col_upper = &.{-std.math.inf(f64)};
    try std.testing.expectError(error.InvalidBounds, view.validate());
}

test "ProblemView validates the full canonical CSC contract" {
    var view = ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &.{0.0},
        .col_lower = &.{0.0},
        .col_upper = &.{1.0},
        .row_lower = &.{-std.math.inf(f64)},
        .row_upper = &.{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &[_]usize{ 1, 1 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    try std.testing.expectError(error.InvalidMatrix, view.validate());

    view.matrix = matrix.CscView.initAssumeValid(
        1,
        1,
        &[_]usize{ 0, 1 },
        &.{foundation.RowId.fromUsizeAssumeValid(0)},
        &.{0.0},
    );
    try std.testing.expectError(error.InvalidMatrix, view.validate());
}