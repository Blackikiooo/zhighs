//! Lazy update system for `Model`.
//!
//! All queued modifications are flushed to the committed state by
//! `updateModel` / `applyPending`.
//!
//! ## Responsibility
//!
//! Interprets `PendingChange` records and atomically updates committed arrays,
//! names, dimensions, and the constraint matrix.  Queue ownership and cleanup
//! are defined in `model_pending.zig`; public edit requests are created by
//! `model_linear.zig`.

const std = @import("std");
const types = @import("types.zig");
const Model = @import("model.zig").Model;
const genconstr = @import("model_genconstr.zig");
const foundation = @import("foundation");
const matrix = @import("matrix");

const ModelError = types.ModelError;
const VarType = types.VarType;
const Sense = types.Sense;
const GenConstrType = types.GenConstrType;
const RowId = foundation.RowId;
const CscMatrix = matrix.CscMatrix;
const MatrixStore = matrix.MatrixStore;

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

    var structure_changed = false;
    var matrix_values_changed = false;
    var bounds_changed = false;
    var objective_changed = false;
    for (self.pending.items) |change| switch (change) {
        .add_var, .add_constr, .del_vars, .del_constrs, .chg_type => structure_changed = true,
        .chg_coeff => matrix_values_changed = true,
        .chg_bounds, .chg_rhs, .chg_sense => bounds_changed = true,
        .chg_obj => objective_changed = true,
    };

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
    var total_vars = self.num_vars + new_vars;
    var total_constrs = self.num_constrs + new_constrs;

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

    // Deletion compacts the linear core and remaps every packed store that
    // references variables (quadratic objective/constraints, PWL and SOS).
    // General constraints still require a dedicated remapping pass; refusing
    // that mixed case is safer than leaving stale indices.
    var delete_var_count: usize = 0;
    var delete_constr_count: usize = 0;
    for (self.pending.items) |chg| switch (chg) {
        .del_vars => |d| delete_var_count += d.indices.len,
        .del_constrs => |d| delete_constr_count += d.indices.len,
        else => {},
    };
    if (delete_var_count > 0 or delete_constr_count > 0) {
        var deleted_vars = try alloc.alloc(bool, total_vars);
        defer alloc.free(deleted_vars);
        var deleted_constrs = try alloc.alloc(bool, total_constrs);
        defer alloc.free(deleted_constrs);
        @memset(deleted_vars, false);
        @memset(deleted_constrs, false);

        for (self.pending.items) |chg| switch (chg) {
            .del_vars => |d| for (d.indices) |idx| {
                if (idx >= total_vars) return error.InvalidArgument;
                deleted_vars[idx] = true;
            },
            .del_constrs => |d| for (d.indices) |idx| {
                if (idx >= total_constrs) return error.InvalidArgument;
                deleted_constrs[idx] = true;
            },
            else => {},
        };

        while (self.var_handles.liveLen() < total_vars) {
            const id = self.var_handles.allocate(alloc) catch return error.OutOfMemory;
            self.var_handles.bindDenseWithAllocator(alloc, id, @intCast(self.var_handles.liveLen())) catch return error.OutOfMemory;
        }
        const var_remap = self.var_handles.compact(alloc, deleted_vars) catch |err| return switch (err) {
            error.DenseIndexOutOfRange, error.InvalidHandle => error.InvalidArgument,
            error.HandleExhausted => error.OutOfMemory,
        };
        defer alloc.free(var_remap);
        while (self.genconstr_handles.liveLen() < self.genconstr_count) {
            const id = self.genconstr_handles.allocate(alloc) catch return error.OutOfMemory;
            self.genconstr_handles.bindDenseWithAllocator(alloc, id, @intCast(self.genconstr_handles.liveLen())) catch return error.OutOfMemory;
        }

        // General constraints use a packed representation with a small
        // type-specific set of variable slots.  Constraints whose result or
        // any operand is deleted are removed; surviving constraints are
        // copied with their variable slots remapped in place.
        if (self.genconstr_count > 0 and delete_var_count > 0) {
            var new_types: std.ArrayListUnmanaged(GenConstrType) = .empty;
            var new_resvars: std.ArrayListUnmanaged(usize) = .empty;
            var new_nvars: std.ArrayListUnmanaged(usize) = .empty;
            var new_indices: std.ArrayListUnmanaged(usize) = .empty;
            var new_begin: std.ArrayListUnmanaged(usize) = .empty;
            var new_names: std.ArrayListUnmanaged(?[]const u8) = .empty;
            new_begin.append(alloc, 0) catch return error.OutOfMemory;
            var gen_deleted = try alloc.alloc(bool, self.genconstr_count);
            defer alloc.free(gen_deleted);
            @memset(gen_deleted, false);

            for (0..self.genconstr_count) |i| {
                const t = self.genconstr_types[i];
                const nvars = self.genconstr_nvars[i];
                const off = genconstr.genConstrOffset(self.*, i);
                const len = genconstr.genConstrDataLen(self.*, i);
                const data = self.genconstr_indices[off .. off + len];
                var drop = deleted_vars[self.genconstr_resvar[i]];
                var kept_nvars = nvars;
                switch (t) {
                    .max, .min, .abs, .and_, .or_, .norm, .nl, .exp, .log, .sin, .cos, .tan, .logistic => {
                        kept_nvars = 0;
                        for (data) |v| {
                            if (!deleted_vars[v]) kept_nvars += 1;
                        }
                        if (kept_nvars == 0) drop = true;
                    },
                    .indicator => {
                        kept_nvars = 0;
                        for (0..nvars) |j| {
                            if (!deleted_vars[data[4 + 2 * j]]) kept_nvars += 1;
                        }
                    },
                    .pwl, .expa, .loga, .pow, .poly => {
                        if (deleted_vars[data[0]]) drop = true;
                    },
                }
                if (drop) {
                    gen_deleted[i] = true;
                    if (self.genconstr_names[i]) |name| alloc.free(name);
                    continue;
                }

                new_types.append(alloc, t) catch return error.OutOfMemory;
                new_resvars.append(alloc, var_remap[self.genconstr_resvar[i]]) catch return error.OutOfMemory;
                new_nvars.append(alloc, kept_nvars) catch return error.OutOfMemory;
                new_names.append(alloc, self.genconstr_names[i]) catch return error.OutOfMemory;
                switch (t) {
                    .max, .min, .abs, .and_, .or_, .norm, .nl, .exp, .log, .sin, .cos, .tan, .logistic => {
                        for (data) |value| {
                            if (!deleted_vars[value]) {
                                new_indices.append(alloc, var_remap[value]) catch return error.OutOfMemory;
                            }
                        }
                    },
                    .indicator => {
                        new_indices.appendSlice(alloc, data[0..4]) catch return error.OutOfMemory;
                        new_indices.items[new_indices.items.len - 3] = @as(usize, @bitCast(@as(i64, @intCast(kept_nvars))));
                        for (0..nvars) |j| {
                            if (!deleted_vars[data[4 + 2 * j]]) {
                                new_indices.append(alloc, var_remap[data[4 + 2 * j]]) catch return error.OutOfMemory;
                                new_indices.append(alloc, data[5 + 2 * j]) catch return error.OutOfMemory;
                            }
                        }
                    },
                    .pwl, .expa, .loga, .pow, .poly => {
                        for (data, 0..) |value, j| {
                            const is_var = j == 0;
                            new_indices.append(alloc, if (is_var) var_remap[value] else value) catch return error.OutOfMemory;
                        }
                    },
                }
                new_begin.append(alloc, new_indices.items.len) catch return error.OutOfMemory;
            }
            alloc.free(self.genconstr_types);
            alloc.free(self.genconstr_resvar);
            alloc.free(self.genconstr_nvars);
            alloc.free(self.genconstr_indices);
            alloc.free(self.genconstr_begin);
            alloc.free(self.genconstr_names);
            self.genconstr_types = try new_types.toOwnedSlice(alloc);
            self.genconstr_resvar = try new_resvars.toOwnedSlice(alloc);
            self.genconstr_nvars = try new_nvars.toOwnedSlice(alloc);
            self.genconstr_indices = try new_indices.toOwnedSlice(alloc);
            self.genconstr_begin = try new_begin.toOwnedSlice(alloc);
            self.genconstr_names = try new_names.toOwnedSlice(alloc);
            self.genconstr_count = self.genconstr_types.len;
            const gen_remap = self.genconstr_handles.compact(alloc, gen_deleted) catch return error.OutOfMemory;
            alloc.free(gen_remap);
        }

        // Count unique deletions, not command entries: callers may submit the
        // same index more than once in one batch.
        delete_var_count = 0;
        for (deleted_vars) |marked| {
            if (marked) delete_var_count += 1;
        }
        delete_constr_count = 0;
        for (deleted_constrs) |marked| {
            if (marked) delete_constr_count += 1;
        }

        // Ensure handle tables contain every pre-compaction dense entry.
        while (self.constr_handles.liveLen() < total_constrs) {
            const id = self.constr_handles.allocate(alloc) catch return error.OutOfMemory;
            self.constr_handles.bindDenseWithAllocator(alloc, id, @intCast(self.constr_handles.liveLen())) catch return error.OutOfMemory;
        }

        const constr_remap = self.constr_handles.compact(alloc, deleted_constrs) catch |err| return switch (err) {
            error.DenseIndexOutOfRange, error.InvalidHandle => error.InvalidArgument,
            error.HandleExhausted => error.OutOfMemory,
        };
        defer alloc.free(constr_remap);

        // Quadratic objective terms use dense variable indices. Drop terms
        // touching deleted variables and remap surviving endpoints in one
        // linear pass over the packed triples.
        if (self.q_nz > 0 and delete_var_count > 0) {
            var kept_q: usize = 0;
            for (0..self.q_nz) |i| {
                if (self.q_row[i] < 0 or self.q_col[i] < 0) return error.InvalidArgument;
                const row = @as(usize, @intCast(self.q_row[i]));
                const col = @as(usize, @intCast(self.q_col[i]));
                if (row < deleted_vars.len and col < deleted_vars.len and !deleted_vars[row] and !deleted_vars[col]) {
                    self.q_row[kept_q] = @intCast(var_remap[row]);
                    self.q_col[kept_q] = @intCast(var_remap[col]);
                    self.q_val[kept_q] = self.q_val[i];
                    kept_q += 1;
                }
            }
            self.q_row = try alloc.realloc(self.q_row, kept_q);
            self.q_col = try alloc.realloc(self.q_col, kept_q);
            self.q_val = try alloc.realloc(self.q_val, kept_q);
            self.q_nz = kept_q;
        }

        if (self.qconstr_count > 0 and delete_var_count > 0) {
            if (self.qconstr_qbegin.len != self.qconstr_count + 1 or self.qconstr_lbegin.len != self.qconstr_count + 1)
                return error.InvalidArgument;
            var new_qbegin: std.ArrayListUnmanaged(usize) = .empty;
            var new_qrow: std.ArrayListUnmanaged(i32) = .empty;
            var new_qcol: std.ArrayListUnmanaged(i32) = .empty;
            var new_qval: std.ArrayListUnmanaged(f64) = .empty;
            var new_lbegin: std.ArrayListUnmanaged(usize) = .empty;
            var new_lind: std.ArrayListUnmanaged(usize) = .empty;
            var new_lval: std.ArrayListUnmanaged(f64) = .empty;
            new_qbegin.append(alloc, 0) catch return error.OutOfMemory;
            new_lbegin.append(alloc, 0) catch return error.OutOfMemory;
            for (0..self.qconstr_count) |constraint_index| {
                const qb = self.qconstr_qbegin[constraint_index];
                const qe = self.qconstr_qbegin[constraint_index + 1];
                for (self.qconstr_qrow[qb..qe], self.qconstr_qcol[qb..qe], self.qconstr_qval[qb..qe]) |row, col, value| {
                    if (row < 0 or col < 0) return error.InvalidArgument;
                    const old_row: usize = @intCast(row);
                    const old_col: usize = @intCast(col);
                    if (old_row < deleted_vars.len and old_col < deleted_vars.len and !deleted_vars[old_row] and !deleted_vars[old_col]) {
                        new_qrow.append(alloc, @intCast(var_remap[old_row])) catch return error.OutOfMemory;
                        new_qcol.append(alloc, @intCast(var_remap[old_col])) catch return error.OutOfMemory;
                        new_qval.append(alloc, value) catch return error.OutOfMemory;
                    }
                }
                new_qbegin.append(alloc, new_qrow.items.len) catch return error.OutOfMemory;

                const lb = self.qconstr_lbegin[constraint_index];
                const le = self.qconstr_lbegin[constraint_index + 1];
                for (self.qconstr_lind[lb..le], self.qconstr_lval[lb..le]) |old_var, value| {
                    if (old_var >= deleted_vars.len) return error.InvalidArgument;
                    if (!deleted_vars[old_var]) {
                        new_lind.append(alloc, var_remap[old_var]) catch return error.OutOfMemory;
                        new_lval.append(alloc, value) catch return error.OutOfMemory;
                    }
                }
                new_lbegin.append(alloc, new_lind.items.len) catch return error.OutOfMemory;
            }
            alloc.free(self.qconstr_qrow);
            alloc.free(self.qconstr_qcol);
            alloc.free(self.qconstr_qval);
            alloc.free(self.qconstr_qbegin);
            alloc.free(self.qconstr_lind);
            alloc.free(self.qconstr_lval);
            alloc.free(self.qconstr_lbegin);
            self.qconstr_qrow = try new_qrow.toOwnedSlice(alloc);
            self.qconstr_qcol = try new_qcol.toOwnedSlice(alloc);
            self.qconstr_qval = try new_qval.toOwnedSlice(alloc);
            self.qconstr_qbegin = try new_qbegin.toOwnedSlice(alloc);
            self.qconstr_lind = try new_lind.toOwnedSlice(alloc);
            self.qconstr_lval = try new_lval.toOwnedSlice(alloc);
            self.qconstr_lbegin = try new_lbegin.toOwnedSlice(alloc);
        }

        if (self.pwlobj_count > 0 and delete_var_count > 0) {
            var new_pwl_var: std.ArrayListUnmanaged(usize) = .empty;
            var new_pwl_npts: std.ArrayListUnmanaged(usize) = .empty;
            var new_pwl_x: std.ArrayListUnmanaged(f64) = .empty;
            var new_pwl_y: std.ArrayListUnmanaged(f64) = .empty;
            var offset: usize = 0;
            for (0..self.pwlobj_count) |i| {
                const old_var = self.pwlobj_var[i];
                const npts = self.pwlobj_npts[i];
                if (old_var >= deleted_vars.len) return error.InvalidArgument;
                if (!deleted_vars[old_var]) {
                    new_pwl_var.append(alloc, var_remap[old_var]) catch return error.OutOfMemory;
                    new_pwl_npts.append(alloc, npts) catch return error.OutOfMemory;
                    new_pwl_x.appendSlice(alloc, self.pwlobj_xdata[offset .. offset + npts]) catch return error.OutOfMemory;
                    new_pwl_y.appendSlice(alloc, self.pwlobj_ydata[offset .. offset + npts]) catch return error.OutOfMemory;
                }
                offset += npts;
            }
            self.allocator.free(self.pwlobj_var);
            self.allocator.free(self.pwlobj_npts);
            self.allocator.free(self.pwlobj_xdata);
            self.allocator.free(self.pwlobj_ydata);
            self.pwlobj_var = try new_pwl_var.toOwnedSlice(alloc);
            self.pwlobj_npts = try new_pwl_npts.toOwnedSlice(alloc);
            self.pwlobj_xdata = try new_pwl_x.toOwnedSlice(alloc);
            self.pwlobj_ydata = try new_pwl_y.toOwnedSlice(alloc);
            self.pwlobj_count = self.pwlobj_var.len;
        }

        if (self.sos_count > 0 and delete_var_count > 0) {
            var new_begin: std.ArrayListUnmanaged(usize) = .empty;
            var new_indices: std.ArrayListUnmanaged(usize) = .empty;
            var new_weights: std.ArrayListUnmanaged(f64) = .empty;
            new_begin.append(alloc, 0) catch return error.OutOfMemory;
            for (0..self.sos_count) |sos_index| {
                const begin = self.sos_begin[sos_index];
                const end = self.sos_begin[sos_index + 1];
                for (self.sos_indices[begin..end], self.sos_weights[begin..end]) |old_var, weight| {
                    if (old_var >= deleted_vars.len) return error.InvalidArgument;
                    if (!deleted_vars[old_var]) {
                        new_indices.append(alloc, var_remap[old_var]) catch return error.OutOfMemory;
                        new_weights.append(alloc, weight) catch return error.OutOfMemory;
                    }
                }
                new_begin.append(alloc, new_indices.items.len) catch return error.OutOfMemory;
            }
            self.allocator.free(self.sos_begin);
            self.allocator.free(self.sos_indices);
            self.allocator.free(self.sos_weights);
            self.sos_begin = try new_begin.toOwnedSlice(alloc);
            self.sos_indices = try new_indices.toOwnedSlice(alloc);
            self.sos_weights = try new_weights.toOwnedSlice(alloc);
        }

        const kept_vars = total_vars - delete_var_count;
        const kept_constrs = total_constrs - delete_constr_count;

        var new_lb = try alloc.alloc(f64, kept_vars);
        var new_ub = try alloc.alloc(f64, kept_vars);
        var new_obj = try alloc.alloc(f64, kept_vars);
        var new_type = try alloc.alloc(VarType, kept_vars);
        var new_names = try alloc.alloc(?[]const u8, kept_vars);
        var vp: usize = 0;
        for (0..total_vars) |old| if (!deleted_vars[old]) {
            new_lb[vp] = self.var_lb[old];
            new_ub[vp] = self.var_ub[old];
            new_obj[vp] = self.var_obj[old];
            new_type[vp] = self.var_type[old];
            new_names[vp] = self.var_names[old];
            vp += 1;
        } else if (self.var_names[old]) |name| alloc.free(name);
        alloc.free(self.var_lb);
        alloc.free(self.var_ub);
        alloc.free(self.var_obj);
        alloc.free(self.var_type);
        alloc.free(self.var_names);
        self.var_lb = new_lb;
        self.var_ub = new_ub;
        self.var_obj = new_obj;
        self.var_type = new_type;
        self.var_names = new_names;

        var new_sense = try alloc.alloc(Sense, kept_constrs);
        var new_rhs = try alloc.alloc(f64, kept_constrs);
        var new_cnames = try alloc.alloc(?[]const u8, kept_constrs);
        var cp: usize = 0;
        for (0..total_constrs) |old| if (!deleted_constrs[old]) {
            new_sense[cp] = self.constr_sense[old];
            new_rhs[cp] = self.constr_rhs[old];
            new_cnames[cp] = self.constr_names[old];
            cp += 1;
        } else if (self.constr_names[old]) |name| alloc.free(name);
        alloc.free(self.constr_sense);
        alloc.free(self.constr_rhs);
        alloc.free(self.constr_names);
        self.constr_sense = new_sense;
        self.constr_rhs = new_rhs;
        self.constr_names = new_cnames;

        const old_matrix = self.matrix.csc();
        var col_starts = try alloc.alloc(usize, kept_vars + 1);
        var row_list: std.ArrayListUnmanaged(RowId) = .empty;
        var val_list: std.ArrayListUnmanaged(f64) = .empty;
        col_starts[0] = 0;
        var new_col: usize = 0;
        for (0..total_vars) |old_col| if (!deleted_vars[old_col]) {
            const begin = old_matrix.col_starts[old_col];
            const end = old_matrix.col_starts[old_col + 1];
            for (old_matrix.row_indices[begin..end], old_matrix.values[begin..end]) |row, value| {
                const old_row = row.toUsize();
                if (!deleted_constrs[old_row]) {
                    row_list.append(alloc, RowId.fromUsizeAssumeValid(constr_remap[old_row])) catch return error.OutOfMemory;
                    val_list.append(alloc, value) catch return error.OutOfMemory;
                }
            }
            new_col += 1;
            col_starts[new_col] = row_list.items.len;
        };
        const new_csc = CscMatrix{
            .num_rows = kept_constrs,
            .num_cols = kept_vars,
            .col_starts = col_starts,
            .row_indices = try alloc.dupe(RowId, row_list.items),
            .values = try alloc.dupe(f64, val_list.items),
        };
        row_list.deinit(alloc);
        val_list.deinit(alloc);
        self.matrix.replaceMatrixAssumeValid(alloc, new_csc) catch return error.RevisionOverflow;
        total_vars = kept_vars;
        total_constrs = kept_constrs;
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

    // Materialize stable IDs for every committed dense entity.  This keeps
    // handle resolution O(1) after updateModel and gives later compaction a
    // complete dense-to-slot map to update.
    while (self.var_handles.liveLen() < self.num_vars) {
        const id = self.var_handles.allocate(alloc) catch return error.OutOfMemory;
        self.var_handles.bindDenseWithAllocator(alloc, id, @intCast(self.var_handles.liveLen())) catch return error.OutOfMemory;
    }
    while (self.constr_handles.liveLen() < self.num_constrs) {
        const id = self.constr_handles.allocate(alloc) catch return error.OutOfMemory;
        self.constr_handles.bindDenseWithAllocator(alloc, id, @intCast(self.constr_handles.liveLen())) catch return error.OutOfMemory;
    }
    while (self.qconstr_handles.liveLen() < self.qconstr_count) {
        const id = self.qconstr_handles.allocate(alloc) catch return error.OutOfMemory;
        self.qconstr_handles.bindDenseWithAllocator(alloc, id, @intCast(self.qconstr_handles.liveLen())) catch return error.OutOfMemory;
    }
    while (self.sos_handles.liveLen() < self.sos_count) {
        const id = self.sos_handles.allocate(alloc) catch return error.OutOfMemory;
        self.sos_handles.bindDenseWithAllocator(alloc, id, @intCast(self.sos_handles.liveLen())) catch return error.OutOfMemory;
    }
    while (self.genconstr_handles.liveLen() < self.genconstr_count) {
        const id = self.genconstr_handles.allocate(alloc) catch return error.OutOfMemory;
        self.genconstr_handles.bindDenseWithAllocator(alloc, id, @intCast(self.genconstr_handles.liveLen())) catch return error.OutOfMemory;
    }
    self.revision += 1;
    if (structure_changed) try self.markRevision(.structure);
    if (matrix_values_changed) try self.markRevision(.matrix_values);
    if (bounds_changed) try self.markRevision(.bounds);
    if (objective_changed) try self.markRevision(.objective);
}
