//! Borrowed LP input view for simplex.
//!
//! `ProblemView` never owns or copies model data. It is deliberately
//! independent of the public model layer so the numerical kernel can sit below
//! solver dispatch and model adapters in the dependency graph.

const matrix = @import("matrix");

pub const ObjectiveSense = enum(i2) { minimize = 1, maximize = -1 };

pub const ProblemView = struct {
    num_rows: usize,
    num_cols: usize,
    col_cost: []const f64,
    col_lower: []const f64,
    col_upper: []const f64,
    row_lower: []const f64,
    row_upper: []const f64,
    matrix: matrix.CscView,
    objective_sense: ObjectiveSense,
    objective_offset: f64,

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
        .matrix = .{ .num_rows = 0, .num_cols = 0, .col_starts = &[_]usize{0}, .row_indices = &.{}, .values = &.{} },
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    try view.validate();
}
