//! Borrowed LP input view for simplex.
//!
//! `ProblemView` never owns or copies model data. Its lifetime must be nested
//! inside the immutable `CompiledModel` snapshot that supplied it.

const model = @import("model");
const matrix = @import("matrix");

pub const ProblemView = struct {
    num_rows: usize,
    num_cols: usize,
    col_cost: []const f64,
    col_lower: []const f64,
    col_upper: []const f64,
    row_lower: []const f64,
    row_upper: []const f64,
    matrix: matrix.CscView,
    objective_sense: model.ObjectiveSense,
    objective_offset: f64,

    pub fn fromLinearModel(model_data: *const model.linear_model.LinearModel) ProblemView {
        const csc = model_data.matrix.csc();
        return .{
            .num_rows = model_data.num_rows,
            .num_cols = model_data.num_cols,
            .col_cost = model_data.col_cost,
            .col_lower = model_data.col_lower,
            .col_upper = model_data.col_upper,
            .row_lower = model_data.row_lower,
            .row_upper = model_data.row_upper,
            .matrix = .{
                .num_rows = csc.num_rows,
                .num_cols = csc.num_cols,
                .col_starts = csc.col_starts,
                .row_indices = csc.row_indices,
                .values = csc.values,
            },
            .objective_sense = model_data.objective_sense,
            .objective_offset = model_data.objective_offset,
        };
    }
};
