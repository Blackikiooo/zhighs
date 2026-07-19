//! Borrowed LP input view for simplex.
//!
//! `ProblemView` never owns or copies model data. It is deliberately
//! independent of the public model layer so the numerical kernel can sit below
//! solver dispatch and model adapters in the dependency graph.

const matrix = @import("matrix");

/// Optimization direction. Encoded as +1/-1 so flipping the sense is just a
/// sign change on the objective vector.
pub const ObjectiveSense = enum(i2) { minimize = 1, maximize = -1 };

/// Non-owning view over a sparse LP in CSC form. All slices are borrowed from
/// the caller; the view performs no allocation.
pub const ProblemView = struct {
    num_rows: usize,
    num_cols: usize,
    col_cost: []const f64, // Objective coefficients (length = num_cols)
    col_lower: []const f64, // Column lower bounds (length = num_cols)
    col_upper: []const f64, // Column upper bounds (length = num_cols)
    row_lower: []const f64, // Row lower bounds (length = num_rows)
    row_upper: []const f64, // Row upper bounds (length = num_rows)
    matrix: matrix.CscView, // Constraint matrix in compressed sparse column form
    objective_sense: ObjectiveSense,
    objective_offset: f64, // Constant added to the objective value

    pub const ViewError = error{ DimensionMismatch, InvalidBounds, InvalidMatrix };

    /// Validate borrowed dimensions and finite bound ordering without taking
    /// ownership or allocating. Intended for the compile/dispatch boundary,
    /// not for repeated simplex iterations.
    pub fn validate(self: ProblemView) ViewError!void {
        if (self.col_cost.len != self.num_cols or self.col_lower.len != self.num_cols or self.col_upper.len != self.num_cols or
            self.row_lower.len != self.num_rows or self.row_upper.len != self.num_rows)
            return error.DimensionMismatch;
        if (self.matrix.num_rows != self.num_rows or self.matrix.num_cols != self.num_cols or
            self.matrix.col_starts.len != self.num_cols + 1 or self.matrix.row_indices.len != self.matrix.values.len)
            return error.InvalidMatrix;
        for (self.col_lower, self.col_upper) |lower, upper| {
            if (lower > upper) return error.InvalidBounds;
        }
        for (self.row_lower, self.row_upper) |lower, upper| {
            if (lower > upper) return error.InvalidBounds;
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
