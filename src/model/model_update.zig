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
const edit_plan_module = @import("model_edit_plan.zig");
const foundation = @import("foundation");
const matrix = @import("matrix");

const ModelError = types.ModelError;
const VarType = types.VarType;
const Sense = types.Sense;
const GenConstrType = types.GenConstrType;
const RowId = foundation.RowId;
const CscMatrix = matrix.CscMatrix;
const MatrixStore = matrix.MatrixStore;
const ModelEditPlan = edit_plan_module.ModelEditPlan;

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

    if (edit_plan_module.isDirectScalarSegment(self.pending.items))
        return applyDirectScalarSegment(self);

    var edit_plan = ModelEditPlan.build(alloc, self.pending.items) catch return error.OutOfMemory;
    defer edit_plan.deinit(alloc);
    const structure_changed = edit_plan.has_structure;
    const matrix_values_changed = edit_plan.coefficients.len != 0;
    const bounds_changed = edit_plan.bounds.len != 0 or edit_plan.rhs.len != 0 or edit_plan.senses.len != 0;
    const objective_changed = edit_plan.objective.len != 0;

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
        // Count both appended CSC columns and appended CSR rows. The latter
        // must be scattered into existing/new columns during the same packed
        // rebuild; previously add_constr coefficient payloads were ignored.
        var col_nnz: usize = 0;
        const existing_csc = self.matrix.csc();
        for (0..self.num_vars) |col| {
            const cnt = existing_csc.col_starts[col + 1] - existing_csc.col_starts[col];
            col_nnz += cnt;
        }
        for (self.pending.items) |chg| {
            switch (chg) {
                .add_var => |v| {
                    col_nnz += v.num_nz;
                },
                else => {},
            }
        }
        var counted_added_row: usize = 0;
        for (self.pending.items) |chg| switch (chg) {
            .add_constr => |c| {
                const target_row = self.num_constrs + counted_added_row;
                for (c.cind[0..c.num_nz]) |target_col| {
                    if (addedColumnRowOffset(self.pending.items, self.num_vars, target_col, target_row) == null)
                        col_nnz = std.math.add(usize, col_nnz, 1) catch return error.InvalidArgument;
                }
                counted_added_row += 1;
            },
            else => {},
        };

        // Build the new CSC matrix.
        const num_cols = total_vars;
        const num_rows = total_constrs;
        var new_col_starts = try alloc.alloc(usize, num_cols + 1);

        var new_row_indices = try alloc.alloc(RowId, col_nnz);
        var new_values = try alloc.alloc(f64, col_nnz);

        @memset(new_col_starts, 0);
        for (0..self.num_vars) |col|
            new_col_starts[col + 1] = existing_csc.col_starts[col + 1] - existing_csc.col_starts[col];
        var col = self.num_vars;
        for (self.pending.items) |chg| {
            switch (chg) {
                .add_var => |v| {
                    new_col_starts[col + 1] = v.num_nz;
                    col += 1;
                },
                else => {},
            }
        }
        var offset_counted_row: usize = 0;
        for (self.pending.items) |chg| switch (chg) {
            .add_constr => |c| {
                const target_row = self.num_constrs + offset_counted_row;
                for (c.cind[0..c.num_nz]) |target_col| {
                    if (target_col >= num_cols) return error.InvalidArgument;
                    if (addedColumnRowOffset(self.pending.items, self.num_vars, target_col, target_row) == null)
                        new_col_starts[target_col + 1] += 1;
                }
                offset_counted_row += 1;
            },
            else => {},
        };
        for (0..num_cols) |column| new_col_starts[column + 1] += new_col_starts[column];
        const column_begins = try alloc.dupe(usize, new_col_starts[0..num_cols]);
        defer alloc.free(column_begins);

        // Use starts as per-column cursors. Existing entries precede appended
        // rows; added constraints are visited in increasing new row order.
        for (0..self.num_vars) |column| {
            const begin = existing_csc.col_starts[column];
            const end = existing_csc.col_starts[column + 1];
            const count = end - begin;
            const destination = new_col_starts[column];
            @memcpy(new_row_indices[destination..][0..count], existing_csc.row_indices[begin..end]);
            @memcpy(new_values[destination..][0..count], existing_csc.values[begin..end]);
            new_col_starts[column] += count;
        }
        col = self.num_vars;
        for (self.pending.items) |chg| switch (chg) {
            .add_var => |v| {
                const destination = new_col_starts[col];
                for (v.vind[0..v.num_nz], 0..) |row, index|
                    new_row_indices[destination + index] = RowId.fromUsizeAssumeValid(row);
                @memcpy(new_values[destination..][0..v.num_nz], v.vval);
                new_col_starts[col] += v.num_nz;
                col += 1;
            },
            else => {},
        };
        var added_row: usize = 0;
        for (self.pending.items) |chg| switch (chg) {
            .add_constr => |c| {
                const row_id = RowId.fromUsizeAssumeValid(self.num_constrs + added_row);
                for (c.cind[0..c.num_nz], c.cval[0..c.num_nz]) |target_col, value| {
                    if (addedColumnRowOffset(self.pending.items, self.num_vars, target_col, row_id.toUsize())) |offset| {
                        new_values[column_begins[target_col] + offset] = value;
                        continue;
                    }
                    const destination = new_col_starts[target_col];
                    new_row_indices[destination] = row_id;
                    new_values[destination] = value;
                    new_col_starts[target_col] += 1;
                }
                added_row += 1;
            },
            else => {},
        };
        col = num_cols;
        while (col > 0) {
            new_col_starts[col] = new_col_starts[col - 1];
            col -= 1;
        }
        new_col_starts[0] = 0;

        // Validate row indices are within the total constraint count.
        for (new_row_indices) |ri| {
            if (ri.toUsize() >= num_rows) {
                alloc.free(new_col_starts);
                alloc.free(new_row_indices);
                alloc.free(new_values);
                return error.InvalidArgument;
            }
        }

        // Create the new CSC matrix and replace the store.
        var new_csc = CscMatrix.initOwnedSlicesAssumeValid(num_rows, num_cols, new_col_starts, new_row_indices, new_values);
        // Validate (debug mode catches structural bugs).
        new_csc.validate() catch {
            alloc.free(new_col_starts);
            alloc.free(new_row_indices);
            alloc.free(new_values);
            return error.InvalidArgument;
        };
        self.matrix.replaceMatrixAssumeValid(alloc, new_csc) catch {
            alloc.free(new_col_starts);
            alloc.free(new_row_indices);
            alloc.free(new_values);
            return error.RevisionOverflow;
        };
    }

    // Apply every coefficient delta in one canonical merge. This avoids one
    // complete CSC copy and one cache invalidation per queued change.
    if (matrix_values_changed)
        try applyCoefficientChanges(self, &edit_plan, total_vars, total_constrs);

    // Apply normalized scalar streams once per target. Deletions are handled
    // below after scalar attributes have been committed.
    {
        const fields = edit_plan.bounds.slice();
        for (fields.items(.index), fields.items(.lower), fields.items(.upper)) |index, lower, upper| {
            if (index < self.var_lb.len) {
                self.var_lb[index] = lower;
                self.var_ub[index] = upper;
            }
        }
    }
    {
        const fields = edit_plan.objective.slice();
        for (fields.items(.index), fields.items(.value)) |index, value| {
            if (index < self.var_obj.len) self.var_obj[index] = value;
        }
    }
    {
        const fields = edit_plan.rhs.slice();
        for (fields.items(.index), fields.items(.value)) |index, value| {
            if (index < self.constr_rhs.len) self.constr_rhs[index] = value;
        }
    }
    {
        const fields = edit_plan.senses.slice();
        for (fields.items(.index), fields.items(.value)) |index, value| {
            if (index < self.constr_sense.len) self.constr_sense[index] = value;
        }
    }
    {
        const fields = edit_plan.types.slice();
        for (fields.items(.index), fields.items(.value)) |index, value| {
            if (index < self.var_type.len) self.var_type[index] = value;
        }
    }

    // Deletion compacts the linear core and remaps every packed store that
    // references variables (quadratic objective/constraints, PWL and SOS).
    // General constraints still require a dedicated remapping pass; refusing
    // that mixed case is safer than leaving stale indices.
    const delete_var_count = edit_plan.deleted_vars.items.len;
    const delete_constr_count = edit_plan.deleted_constraints.items.len;
    if (delete_var_count > 0 or delete_constr_count > 0) {
        var deleted_vars = try alloc.alloc(bool, total_vars);
        defer alloc.free(deleted_vars);
        var deleted_constrs = try alloc.alloc(bool, total_constrs);
        defer alloc.free(deleted_constrs);
        @memset(deleted_vars, false);
        @memset(deleted_constrs, false);

        for (edit_plan.deleted_vars.items) |idx| {
            if (idx >= total_vars) return error.InvalidArgument;
            deleted_vars[idx] = true;
        }
        for (edit_plan.deleted_constraints.items) |idx| {
            if (idx >= total_constrs) return error.InvalidArgument;
            deleted_constrs[idx] = true;
        }

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

        // Pass 1 determines the exact surviving nonzero count. The previous
        // ArrayList path repeatedly grew two temporary streams and then copied
        // both into separate owning arrays, temporarily retaining two copies
        // of every surviving entry.
        var kept_nnz: usize = 0;
        for (0..total_vars) |old_col| if (!deleted_vars[old_col]) {
            const begin = old_matrix.col_starts[old_col];
            const end = old_matrix.col_starts[old_col + 1];
            for (old_matrix.row_indices[begin..end]) |row| {
                if (!deleted_constrs[row.toUsize()]) kept_nnz += 1;
            }
        };

        var new_csc = CscMatrix.initPackedUninitialized(alloc, kept_constrs, kept_vars, kept_nnz) catch |err| {
            if (err == error.OutOfMemory) return error.OutOfMemory;
            return error.InvalidArgument;
        };
        var owns_new_csc = true;
        errdefer if (owns_new_csc) new_csc.deinit(alloc);

        // Pass 2 streams surviving columns directly into their final packed
        // row/value spans. Filtering rows preserves the existing canonical
        // order, while constr_remap is monotonic over surviving rows.
        var output_pos: usize = 0;
        var new_col: usize = 0;
        for (0..total_vars) |old_col| if (!deleted_vars[old_col]) {
            new_csc.col_starts[new_col] = output_pos;
            const begin = old_matrix.col_starts[old_col];
            const end = old_matrix.col_starts[old_col + 1];
            for (old_matrix.row_indices[begin..end], old_matrix.values[begin..end]) |row, value| {
                const old_row = row.toUsize();
                if (!deleted_constrs[old_row]) {
                    new_csc.row_indices[output_pos] = RowId.fromUsizeAssumeValid(constr_remap[old_row]);
                    new_csc.values[output_pos] = value;
                    output_pos += 1;
                }
            }
            new_col += 1;
        };
        new_csc.col_starts[kept_vars] = output_pos;
        std.debug.assert(output_pos == kept_nnz);
        std.debug.assert(new_col == kept_vars);
        self.matrix.replaceMatrixAssumeValid(alloc, new_csc) catch return error.RevisionOverflow;
        owns_new_csc = false;
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

fn addedColumnRowOffset(pending: []const @import("model_pending.zig").PendingChange, existing_cols: usize, target_col: usize, target_row: usize) ?usize {
    if (target_col < existing_cols) return null;
    const added_col = target_col - existing_cols;
    var ordinal: usize = 0;
    for (pending) |change| switch (change) {
        .add_var => |column| {
            if (ordinal == added_col) {
                for (column.vind[0..column.num_nz], 0..) |row, offset| {
                    if (row == target_row) return offset;
                }
                return null;
            }
            ordinal += 1;
        },
        else => {},
    };
    return null;
}

fn applyDirectScalarSegment(self: *Model) ModelError!void {
    var structure_changed = false;
    var bounds_changed = false;
    var objective_changed = false;
    for (self.pending.items) |change| switch (change) {
        .chg_bounds => |edit| {
            bounds_changed = true;
            if (edit.var_idx < self.var_lb.len) {
                self.var_lb[edit.var_idx] = edit.lb;
                self.var_ub[edit.var_idx] = edit.ub;
            }
        },
        .chg_obj => |edit| {
            objective_changed = true;
            if (edit.var_idx < self.var_obj.len) self.var_obj[edit.var_idx] = edit.obj;
        },
        .chg_rhs => |edit| {
            bounds_changed = true;
            if (edit.constr_idx < self.constr_rhs.len) self.constr_rhs[edit.constr_idx] = edit.rhs;
        },
        .chg_sense => |edit| {
            bounds_changed = true;
            if (edit.constr_idx < self.constr_sense.len) self.constr_sense[edit.constr_idx] = edit.sense;
        },
        .chg_type => |edit| {
            structure_changed = true;
            if (edit.var_idx < self.var_type.len) self.var_type[edit.var_idx] = edit.vtype;
        },
        else => unreachable,
    };
    self.pending.clearRetainingCapacity();
    self.has_pending = false;
    self.revision += 1;
    if (structure_changed) try self.markRevision(.structure);
    if (bounds_changed) try self.markRevision(.bounds);
    if (objective_changed) try self.markRevision(.objective);
}

/// Sorts and coalesces queued coefficient sets, then merges them with the
/// canonical CSC matrix in two linear passes. Repeated coordinates use
/// last-write-wins semantics. A zero delta removes an existing nonzero; a
/// nonzero delta inserts a previously absent coordinate.
fn applyCoefficientChanges(self: *Model, plan: *ModelEditPlan, num_cols: usize, num_rows: usize) ModelError!void {
    const allocator = self.allocator;
    var fields = plan.coefficients.slice();
    const input_rows = fields.items(.row);
    const input_cols = fields.items(.col);
    var write: usize = 0;
    for (0..plan.coefficients.len) |read| {
        if (input_cols[read] >= num_cols or input_rows[read] >= num_rows) continue;
        fields.set(write, fields.get(read));
        write += 1;
    }
    plan.coefficients.shrinkRetainingCapacity(write);
    if (write == 0) return;
    fields = plan.coefficients.slice();
    const rows = fields.items(.row);
    const cols = fields.items(.col);
    const values = fields.items(.value);
    for (values) |value| if (!std.math.isFinite(value)) return error.InvalidArgument;
    const current = self.matrix.csc();

    // The overwhelmingly common case changes values at coordinates that
    // already exist. Validate every position first, then mutate only the value
    // stream in place: no O(nnz) copy and no structural allocation.
    if (try applyExistingCoefficientValuesFast(self, current.*, rows, cols, values)) return;

    var changed = false;
    var delta_pos: usize = 0;
    var output_nnz: usize = 0;

    // Pass 1 determines exact output sizes and whether the batch has any
    // observable effect, avoiding a replacement for no-op sets.
    for (0..num_cols) |col| {
        var old_pos = current.col_starts[col];
        const old_end = current.col_starts[col + 1];
        var count: usize = 0;
        while (old_pos < old_end or (delta_pos < values.len and cols[delta_pos] == col)) {
            const old_row = if (old_pos < old_end) current.row_indices[old_pos].toUsize() else std.math.maxInt(usize);
            const delta_row = if (delta_pos < values.len and cols[delta_pos] == col) rows[delta_pos] else std.math.maxInt(usize);
            if (old_row < delta_row) {
                count += 1;
                old_pos += 1;
            } else if (delta_row < old_row) {
                if (values[delta_pos] != 0.0) {
                    count += 1;
                    changed = true;
                }
                delta_pos += 1;
            } else {
                if (values[delta_pos] != 0.0) {
                    count += 1;
                    if (values[delta_pos] != current.values[old_pos]) changed = true;
                } else {
                    changed = true;
                }
                old_pos += 1;
                delta_pos += 1;
            }
        }
        output_nnz = std.math.add(usize, output_nnz, count) catch return error.InvalidArgument;
    }
    if (!changed) return;

    var replacement = CscMatrix.initPackedUninitialized(allocator, num_rows, num_cols, output_nnz) catch |err| {
        if (err == error.OutOfMemory) return error.OutOfMemory;
        return error.InvalidArgument;
    };
    errdefer replacement.deinit(allocator);

    // Pass 2 writes canonical rows directly into their final column spans.
    delta_pos = 0;
    var output_pos: usize = 0;
    for (0..num_cols) |col| {
        replacement.col_starts[col] = output_pos;
        var old_pos = current.col_starts[col];
        const old_end = current.col_starts[col + 1];
        while (old_pos < old_end or (delta_pos < values.len and cols[delta_pos] == col)) {
            const old_row = if (old_pos < old_end) current.row_indices[old_pos].toUsize() else std.math.maxInt(usize);
            const delta_row = if (delta_pos < values.len and cols[delta_pos] == col) rows[delta_pos] else std.math.maxInt(usize);
            if (old_row < delta_row) {
                replacement.row_indices[output_pos] = current.row_indices[old_pos];
                replacement.values[output_pos] = current.values[old_pos];
                output_pos += 1;
                old_pos += 1;
            } else if (delta_row < old_row) {
                if (values[delta_pos] != 0.0) {
                    replacement.row_indices[output_pos] = RowId.fromUsizeAssumeValid(rows[delta_pos]);
                    replacement.values[output_pos] = values[delta_pos];
                    output_pos += 1;
                }
                delta_pos += 1;
            } else {
                if (values[delta_pos] != 0.0) {
                    replacement.row_indices[output_pos] = current.row_indices[old_pos];
                    replacement.values[output_pos] = values[delta_pos];
                    output_pos += 1;
                }
                old_pos += 1;
                delta_pos += 1;
            }
        }
    }
    replacement.col_starts[num_cols] = output_pos;
    std.debug.assert(output_pos == output_nnz);
    std.debug.assert(delta_pos == values.len);
    self.matrix.replaceMatrixAssumeValid(allocator, replacement) catch return error.RevisionOverflow;
}

fn applyExistingCoefficientValuesFast(self: *Model, current: CscMatrix, rows: []const usize, cols: []const usize, new_values: []const f64) ModelError!bool {
    const allocator = self.allocator;
    const positions = try allocator.alloc(usize, new_values.len);
    defer allocator.free(positions);
    const values = try allocator.alloc(f64, new_values.len);
    defer allocator.free(values);

    var changed = false;
    for (rows, cols, new_values, 0..) |row, col, value, index| {
        if (value == 0.0) return false;
        const position = findRowPosition(current, col, row) orelse return false;
        positions[index] = position;
        values[index] = value;
        if (current.values[position] != value) changed = true;
    }
    if (!changed) return true;
    self.matrix.updateValuesAtPositionsAssumeValid(allocator, positions, values) catch return error.RevisionOverflow;
    return true;
}

fn findRowPosition(matrix_value: CscMatrix, col: usize, target_row: usize) ?usize {
    var low = matrix_value.col_starts[col];
    var high = matrix_value.col_starts[col + 1];
    while (low < high) {
        const middle = low + (high - low) / 2;
        const row = matrix_value.row_indices[middle].toUsize();
        if (row < target_row) {
            low = middle + 1;
        } else {
            high = middle;
        }
    }
    if (low < matrix_value.col_starts[col + 1] and matrix_value.row_indices[low].toUsize() == target_row)
        return low;
    return null;
}
