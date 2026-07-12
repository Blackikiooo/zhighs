//! Lazy update system for `Model`.
//!
//! All queued modifications are flushed to the committed state by
//! `updateModel` / `applyPending`.

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

/// Apply all queued modifications.
///
/// Until this is called, queries reflect the last committed state.
pub fn updateModel(self: *Model) ModelError!void {
    if (!self.has_pending) return;
    try self.applyPending();
}

/// Apply all queued modifications to the committed state.
pub fn applyPending(self: *Model) ModelError!void {
    const alloc = self.allocator;

    // Count how many vars and constrs will be added.
    var new_vars: usize = 0;
    var new_constrs: usize = 0;
    for (self.pending.items) |chg| {
        switch (chg) {
            .add_var => new_vars += 1,
            .add_constr => new_constrs += 1,
            else => {},
        }
    }

    // Allocate space for the expanded arrays.
    const total_vars = self.num_vars + new_vars;
    const total_constrs = self.num_constrs + new_constrs;

    // Commit variable additions.
    if (new_vars > 0) {
        var lb = try alloc.alloc(f64, total_vars);
        var ub = try alloc.alloc(f64, total_vars);
        var obj = try alloc.alloc(f64, total_vars);
        var vtype = try alloc.alloc(VarType, total_vars);
        var vnames = try alloc.alloc(?[]const u8, total_vars);

        @memcpy(lb[0..self.num_vars], self.var_lb);
        @memcpy(ub[0..self.num_vars], self.var_ub);
        @memcpy(obj[0..self.num_vars], self.var_obj);
        @memcpy(vtype[0..self.num_vars], self.var_type);
        @memcpy(vnames[0..self.num_vars], self.var_names);

        var idx = self.num_vars;
        for (self.pending.items) |chg| {
            switch (chg) {
                .add_var => |v| {
                    lb[idx] = v.lb;
                    ub[idx] = v.ub;
                    obj[idx] = v.obj;
                    vtype[idx] = v.vtype;
                    vnames[idx] = v.name;
                    idx += 1;
                },
                else => {},
            }
        }

        alloc.free(self.var_lb);
        alloc.free(self.var_ub);
        alloc.free(self.var_obj);
        alloc.free(self.var_type);
        alloc.free(self.var_names);
        self.var_lb = lb;
        self.var_ub = ub;
        self.var_obj = obj;
        self.var_type = vtype;
        self.var_names = vnames;
    }

    // Commit constraint additions.
    if (new_constrs > 0) {
        var sense = try alloc.alloc(Sense, total_constrs);
        var rhs = try alloc.alloc(f64, total_constrs);
        var cnames = try alloc.alloc(?[]const u8, total_constrs);

        @memcpy(sense[0..self.num_constrs], self.constr_sense);
        @memcpy(rhs[0..self.num_constrs], self.constr_rhs);
        @memcpy(cnames[0..self.num_constrs], self.constr_names);

        var idx = self.num_constrs;
        for (self.pending.items) |chg| {
            switch (chg) {
                .add_constr => |c| {
                    sense[idx] = c.sense;
                    rhs[idx] = c.rhs;
                    cnames[idx] = c.name;
                    idx += 1;
                },
                else => {},
            }
        }

        alloc.free(self.constr_sense);
        alloc.free(self.constr_rhs);
        alloc.free(self.constr_names);
        self.constr_sense = sense;
        self.constr_rhs = rhs;
        self.constr_names = cnames;
    }

    // Build the constraint matrix by collecting all column data
    // from pending add_var entries.
    if (new_vars > 0 or new_constrs > 0) {
        // Count column non-zero entries from pending add_var changes.
        var col_nnz: usize = 0;
        var col_nz_counts: std.ArrayListUnmanaged(usize) = .empty;
        // Initialise with existing columns.
        const existing_csc = self.matrix.csc();
        for (0..self.num_vars) |col| {
            const cnt = existing_csc.col_starts[col + 1] - existing_csc.col_starts[col];
            col_nz_counts.append(alloc, cnt) catch return error.OutOfMemory;
            col_nnz += cnt;
        }
        // Add pending columns.
        for (self.pending.items) |chg| {
            switch (chg) {
                .add_var => |v| {
                    col_nz_counts.append(alloc, v.num_nz) catch return error.OutOfMemory;
                    col_nnz += v.num_nz;
                },
                else => {},
            }
        }

        // Build the new CSC matrix.
        const num_cols = total_vars;
        const num_rows = total_constrs;
        var new_col_starts = try alloc.alloc(usize, num_cols + 1);

        var new_row_indices = try alloc.alloc(RowId, col_nnz);
        var new_values = try alloc.alloc(f64, col_nnz);

        new_col_starts[0] = 0;
        var pos: usize = 0;

        // Copy existing columns.
        var col: usize = 0;
        while (col < self.num_vars) : (col += 1) {
            const begin = existing_csc.col_starts[col];
            const end = existing_csc.col_starts[col + 1];
            const len = end - begin;
            @memcpy(new_row_indices[pos .. pos + len], existing_csc.row_indices[begin..end]);
            @memcpy(new_values[pos .. pos + len], existing_csc.values[begin..end]);
            pos += len;
            new_col_starts[col + 1] = pos;
        }

        // Add pending columns.
        for (self.pending.items) |chg| {
            switch (chg) {
                .add_var => |v| {
                    for (v.vind[0..v.num_nz], 0..) |ci, j| {
                        new_row_indices[pos + j] = RowId.fromUsizeAssumeValid(ci);
                    }
                    @memcpy(new_values[pos .. pos + v.num_nz], v.vval);
                    pos += v.num_nz;
                    new_col_starts[col + 1] = pos;
                    col += 1;
                },
                else => {},
            }
        }

        // Validate row indices are within the total constraint count.
        for (new_row_indices) |ri| {
            if (ri.toUsize() >= num_rows) {
                alloc.free(new_col_starts);
                alloc.free(new_row_indices);
                alloc.free(new_values);
                col_nz_counts.deinit(alloc);
                return error.InvalidArgument;
            }
        }

        // Create the new CSC matrix and replace the store.
        var new_csc = CscMatrix{
            .num_rows = num_rows,
            .num_cols = num_cols,
            .col_starts = new_col_starts,
            .row_indices = new_row_indices,
            .values = new_values,
        };
        // Validate (debug mode catches structural bugs).
        new_csc.validate() catch {
            alloc.free(new_col_starts);
            alloc.free(new_row_indices);
            alloc.free(new_values);
            col_nz_counts.deinit(alloc);
            return error.InvalidArgument;
        };
        self.matrix.replaceMatrixAssumeValid(alloc, new_csc) catch {
            alloc.free(new_col_starts);
            alloc.free(new_row_indices);
            alloc.free(new_values);
            col_nz_counts.deinit(alloc);
            return error.RevisionOverflow;
        };
        col_nz_counts.deinit(alloc);
    }

    // Apply coefficient changes, deletions, etc.
    for (self.pending.items) |chg| {
        switch (chg) {
            .chg_coeff => |c| {
                if (c.var_idx < total_vars and c.constr_idx < total_constrs) {
                    const cur = self.matrix.csc();
                    const cs = cur.col_starts[c.var_idx];
                    const ce = cur.col_starts[c.var_idx + 1];
                    if (ce > cs) {
                        for (cur.row_indices[cs..ce], 0..) |rid, off| {
                            if (rid.toUsize() == c.constr_idx) {
                                var new_vals = try alloc.dupe(f64, cur.values);
                                new_vals[cs + off] = c.new_val;
                                const new_csc = CscMatrix{
                                    .num_rows = cur.num_rows,
                                    .num_cols = cur.num_cols,
                                    .col_starts = try alloc.dupe(usize, cur.col_starts),
                                    .row_indices = try alloc.dupe(RowId, cur.row_indices),
                                    .values = new_vals,
                                };
                                self.matrix.replaceMatrixAssumeValid(alloc, new_csc) catch {};
                                break;
                            }
                        }
                    }
                }
            },
            .chg_bounds => |b| {
                if (b.var_idx < self.var_lb.len) {
                    self.var_lb[b.var_idx] = b.lb;
                    self.var_ub[b.var_idx] = b.ub;
                }
            },
            .chg_obj => |o| {
                if (o.var_idx < self.var_obj.len) {
                    self.var_obj[o.var_idx] = o.obj;
                }
            },
            .chg_rhs => |r| {
                if (r.constr_idx < self.constr_rhs.len) {
                    self.constr_rhs[r.constr_idx] = r.rhs;
                }
            },
            .chg_sense => |s| {
                if (s.constr_idx < self.constr_sense.len) {
                    self.constr_sense[s.constr_idx] = s.sense;
                }
            },
            .chg_type => |t| {
                if (t.var_idx < self.var_type.len) {
                    self.var_type[t.var_idx] = t.vtype;
                }
            },
            else => {},
        }
    }

    // Free the pending change records (owned slices were moved or freed).
    for (self.pending.items) |*chg| {
        switch (chg.*) {
            .add_var => |*v| {
                alloc.free(v.vind);
                alloc.free(v.vval);
            },
            .add_constr => |*c| {
                alloc.free(c.cind);
                alloc.free(c.cval);
            },
            .del_vars => |d| alloc.free(d.indices),
            .del_constrs => |d| alloc.free(d.indices),
            else => {},
        }
    }
    self.pending.clearRetainingCapacity();
    self.has_pending = false;

    self.num_vars = total_vars;
    self.num_constrs = total_constrs;
    self.revision += 1;
}
