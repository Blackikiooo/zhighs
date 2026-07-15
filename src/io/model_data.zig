//! Owning model interchange representation returned by parsers.
//!
//! Numeric arrays and the canonical CSC matrix can be moved into a solver
//! model without reparsing. Names are individually owned for simple, explicit
//! lifetime management; a pooled representation can replace this internally
//! without changing the public view.

const std = @import("std");
const matrix = @import("matrix");
const types = @import("types.zig");

pub const ModelData = struct {
    allocator: std.mem.Allocator,
    name: []u8,
    objective_sense: types.ObjectiveSense = .minimize,
    objective_offset: f64 = 0.0,
    col_cost: []f64,
    col_lower: []f64,
    col_upper: []f64,
    col_type: []types.VariableType,
    col_names: []?[]u8,
    row_lower: []f64,
    row_upper: []f64,
    row_names: []?[]u8,
    matrix: matrix.CscMatrix,

    pub fn deinit(self: *ModelData) void {
        const allocator = self.allocator;
        allocator.free(self.name);
        allocator.free(self.col_cost);
        allocator.free(self.col_lower);
        allocator.free(self.col_upper);
        allocator.free(self.col_type);
        for (self.col_names) |name| if (name) |value| allocator.free(value);
        allocator.free(self.col_names);
        allocator.free(self.row_lower);
        allocator.free(self.row_upper);
        for (self.row_names) |name| if (name) |value| allocator.free(value);
        allocator.free(self.row_names);
        self.matrix.deinit(allocator);
        self.* = undefined;
    }

    pub fn view(self: *const ModelData) types.ModelView {
        return .{
            .name = self.name,
            .objective_sense = self.objective_sense,
            .objective_offset = self.objective_offset,
            .col_cost = self.col_cost,
            .col_lower = self.col_lower,
            .col_upper = self.col_upper,
            .col_type = self.col_type,
            .col_names = @ptrCast(self.col_names),
            .row_lower = self.row_lower,
            .row_upper = self.row_upper,
            .row_names = @ptrCast(self.row_names),
            .matrix = self.matrix.view(),
        };
    }
};
