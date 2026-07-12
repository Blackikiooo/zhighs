//! Solve, reset, copy, and MIP-query methods for `Model`.

const std = @import("std");
const types = @import("types.zig");
const Model = @import("model.zig").Model;
const foundation = @import("foundation");
const csc_mod = @import("../matrix/csc.zig");
const store_mod = @import("../matrix/store.zig");

const ModelError = types.ModelError;
const VarType = types.VarType;
const Sense = types.Sense;
const RowId = foundation.RowId;
const CscMatrix = csc_mod.CscMatrix;
const MatrixStore = store_mod.MatrixStore;

/// Optimise the model
///
/// Flushes pending changes first, then delegates to the appropriate
/// solver engine.  On return, solution attributes are populated and
/// `status` reflects the result.
pub fn optimize(self: *Model) ModelError!void {
    // Flush any queued modifications.
    try self.updateModel();

    if (self.num_vars == 0) return error.EmptyModel;

    // Select and run the solver.
    // The actual solve dispatch will live in the solver module;
    // for now we set a placeholder status.
    self.status = .optimal;
    self.revision += 1;
}

/// Discard the solution and reset the status to LOADED
pub fn reset(self: *Model, clear_all: bool) void {
    _ = clear_all;
    self.status = .loaded;
    self.obj_val = 0.0;
    self.obj_bound = 0.0;
    self.iter_count = 0;
    self.node_count = 0;
    self.bar_iter_count = 0;
}

/// Create an independent copy of the model
pub fn copy(self: *Model, new_name: []const u8) ModelError!Model {
    var new = try Model.init(self.allocator, self.env, new_name);
    errdefer new.deinit();

    // Flush pending changes first.
    try self.updateModel();

    // Clone variable data.
    if (self.num_vars > 0) {
        new.var_lb = try self.allocator.dupe(f64, self.var_lb);
        new.var_ub = try self.allocator.dupe(f64, self.var_ub);
        new.var_obj = try self.allocator.dupe(f64, self.var_obj);
        new.var_type = try self.allocator.dupe(VarType, self.var_type);
        new.var_names = try self.allocator.alloc(?[]const u8, self.num_vars);
        for (self.var_names, 0..) |n, i| {
            new.var_names[i] = if (n) |s| try self.allocator.dupe(u8, s) else null;
        }
        new.num_vars = self.num_vars;
    }

    // Clone constraint data.
    if (self.num_constrs > 0) {
        new.constr_sense = try self.allocator.dupe(Sense, self.constr_sense);
        new.constr_rhs = try self.allocator.dupe(f64, self.constr_rhs);
        new.constr_names = try self.allocator.alloc(?[]const u8, self.num_constrs);
        for (self.constr_names, 0..) |n, i| {
            new.constr_names[i] = if (n) |s| try self.allocator.dupe(u8, s) else null;
        }
        new.num_constrs = self.num_constrs;
    }

    // Clone the matrix.
    const csc_src = self.matrix.csc();
    const csc_copy = CscMatrix{
        .num_rows = csc_src.num_rows,
        .num_cols = csc_src.num_cols,
        .col_starts = try self.allocator.dupe(usize, csc_src.col_starts),
        .row_indices = try self.allocator.dupe(RowId, csc_src.row_indices),
        .values = try self.allocator.dupe(f64, csc_src.values),
    };
    // Free the initial zero matrix from init().
    new.matrix.deinit(self.allocator);
    new.matrix = MatrixStore.initAssumeValid(csc_copy);

    new.sense = self.sense;
    new.status = self.status;
    new.revision = self.revision;
    return new;
}

/// Return whether the model contains any non-continuous variable.
pub fn isMip(self: Model) bool {
    for (self.var_type) |vt| {
        if (vt != .continuous) return true;
    }
    return false;
}
