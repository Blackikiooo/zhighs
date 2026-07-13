//! Revision-cached, borrowed compilation view for continuous LP models.
//!
//! ## Responsibility
//!
//! Borrows the committed variable SoA arrays and CSC matrix without copying,
//! while owning only the row-bound arrays derived from public `Sense + RHS`
//! storage. The cache is rebuilt exactly when the source model revision
//! changes. Solver workspaces must not retain a view after the source model is
//! mutated or destroyed.

const std = @import("std");
const types = @import("types.zig");
const matrix = @import("matrix");

const Sense = types.Sense;
const ObjectiveSense = types.ObjectiveSense;

pub const CompileViewError = error{ OutOfMemory, DimensionMismatch };

/// Non-owning source projection of committed public model storage.
pub const ModelSourceView = struct {
    revision: u64,
    objective_sense: ObjectiveSense,
    objective_offset: f64,
    col_cost: []const f64,
    col_lower: []const f64,
    col_upper: []const f64,
    row_sense: []const Sense,
    row_rhs: []const f64,
    matrix: matrix.CscView,
};

/// Non-owning canonical LP projection consumed by solver adapters.
pub const CompiledLinearModelView = struct {
    source_revision: u64,
    objective_sense: ObjectiveSense,
    objective_offset: f64,
    num_rows: usize,
    num_cols: usize,
    col_cost: []const f64,
    col_lower: []const f64,
    col_upper: []const f64,
    row_lower: []const f64,
    row_upper: []const f64,
    matrix: matrix.CscView,
};

/// Owning storage for data that cannot be borrowed in canonical form.
pub const CompiledModelViewCache = struct {
    row_lower: []f64 = &.{},
    row_upper: []f64 = &.{},
    source_revision: ?u64 = null,

    pub fn deinit(self: *CompiledModelViewCache, allocator: std.mem.Allocator) void {
        allocator.free(self.row_lower);
        allocator.free(self.row_upper);
        self.* = .{};
    }

    /// Compile a canonical LP view, reusing derived row bounds when possible.
    pub fn compileLinearView(self: *CompiledModelViewCache, allocator: std.mem.Allocator, source: ModelSourceView) CompileViewError!CompiledLinearModelView {
        const num_cols = source.col_cost.len;
        const num_rows = source.row_rhs.len;
        if (source.col_lower.len != num_cols or source.col_upper.len != num_cols or
            source.row_sense.len != num_rows or source.matrix.num_cols != num_cols or source.matrix.num_rows != num_rows)
            return error.DimensionMismatch;

        if (self.source_revision == null or self.source_revision.? != source.revision or
            self.row_lower.len != num_rows or self.row_upper.len != num_rows)
        {
            const row_lower = allocator.alloc(f64, num_rows) catch return error.OutOfMemory;
            errdefer allocator.free(row_lower);
            const row_upper = allocator.alloc(f64, num_rows) catch return error.OutOfMemory;
            errdefer allocator.free(row_upper);

            for (source.row_sense, source.row_rhs, row_lower, row_upper) |sense, rhs, *lower, *upper| {
                switch (sense) {
                    .less_equal => {
                        lower.* = -std.math.inf(f64);
                        upper.* = rhs;
                    },
                    .equal => {
                        lower.* = rhs;
                        upper.* = rhs;
                    },
                    .greater_equal => {
                        lower.* = rhs;
                        upper.* = std.math.inf(f64);
                    },
                }
            }

            allocator.free(self.row_lower);
            allocator.free(self.row_upper);
            self.row_lower = row_lower;
            self.row_upper = row_upper;
            self.source_revision = source.revision;
        }

        return .{
            .source_revision = source.revision,
            .objective_sense = source.objective_sense,
            .objective_offset = source.objective_offset,
            .num_rows = num_rows,
            .num_cols = num_cols,
            .col_cost = source.col_cost,
            .col_lower = source.col_lower,
            .col_upper = source.col_upper,
            .row_lower = self.row_lower,
            .row_upper = self.row_upper,
            .matrix = source.matrix,
        };
    }
};

test {
    std.testing.refAllDecls(@This());
}
