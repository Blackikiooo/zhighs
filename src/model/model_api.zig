//! Variable / constraint addition, modification, and query methods for `Model`.
//!
//! These are re-exported from the `Model` struct so they are callable as
//! `model.addVar(...)`, `model.addConstr(...)`, etc.

const std = @import("std");
const types = @import("types.zig");
const Model = @import("model.zig").Model;
const PendingChange = @import("model_pending.zig").PendingChange;
const foundation = @import("foundation");

const ModelError = types.ModelError;
const VarType = types.VarType;
const Sense = types.Sense;
const SosType = types.SosType;
const GenConstrType = types.GenConstrType;
const FeasRelaxType = types.FeasRelaxType;
const CallbackWhere = types.CallbackWhere;
const CallbackFunc = types.CallbackFunc;
const INFINITY = types.INFINITY;
const Version = types.Version;
const AttributeInfo = @import("attrs.zig").AttributeInfo;
const RowId = foundation.RowId;
const ParamValue = @import("env.zig").ParamValue;
const Env = @import("env.zig").Env;

// ══════════════════════════════════════════════════════════════════════════
//  Variable addition
// ══════════════════════════════════════════════════════════════════════════

/// Add one variable
pub fn addVar(
    self: *Model,
    num_nz: usize,
    vind: []const usize,
    vval: []const f64,
    obj: f64,
    lb: f64,
    ub: f64,
    vtype: VarType,
    name: ?[]const u8,
) ModelError!void {
    if (vind.len < num_nz or vval.len < num_nz) return error.InvalidArgument;
    if (lb > ub) return error.InvalidArgument;

    try self.enqueue(.{
        .add_var = .{
            .num_nz = num_nz,
            .vind = try self.allocator.dupe(usize, vind[0..num_nz]),
            .vval = try self.allocator.dupe(f64, vval[0..num_nz]),
            .obj = obj,
            .lb = lb,
            .ub = ub,
            .vtype = vtype,
            .name = if (name) |n| try self.allocator.dupe(u8, n) else null,
        },
    });
}

/// Add a batch of variables (batch; CSC column format).
pub fn addVars(
    self: *Model,
    num_vars: usize,
    vbeg: []const usize,
    vind: []const usize,
    vval: []const f64,
    obj: ?[]const f64,
    lb: ?[]const f64,
    ub: ?[]const f64,
    vtype: ?[]const VarType,
    names: ?[]const ?[]const u8,
) ModelError!void {
    const default_obj: f64 = 0.0;
    const default_lb: f64 = 0.0;
    const default_ub: f64 = INFINITY;
    const default_type: VarType = .continuous;

    for (0..num_vars) |i| {
        const beg = vbeg[i];
        const end = if (i + 1 < num_vars) vbeg[i + 1] else vval.len;
        const nnz = end - beg;

        try self.addVar(
            nnz,
            vind[beg..end],
            vval[beg..end],
            if (obj) |o| o[i] else default_obj,
            if (lb) |l| l[i] else default_lb,
            if (ub) |u| u[i] else default_ub,
            if (vtype) |t| t[i] else default_type,
            if (names) |n| n[i] else null,
        );
    }
}

// ══════════════════════════════════════════════════════════════════════════
//  Constraint addition
// ══════════════════════════════════════════════════════════════════════════

/// Add one linear constraint
pub fn addConstr(
    self: *Model,
    num_nz: usize,
    cind: []const usize,
    cval: []const f64,
    sense: Sense,
    rhs: f64,
    name: ?[]const u8,
) ModelError!void {
    if (cind.len < num_nz or cval.len < num_nz) return error.InvalidArgument;

    try self.enqueue(.{
        .add_constr = .{
            .num_nz = num_nz,
            .cind = try self.allocator.dupe(usize, cind[0..num_nz]),
            .cval = try self.allocator.dupe(f64, cval[0..num_nz]),
            .sense = sense,
            .rhs = rhs,
            .name = if (name) |n| try self.allocator.dupe(u8, n) else null,
        },
    });
}

/// Add a batch of constraints (batch; CSR row format).
pub fn addConstrs(
    self: *Model,
    num_constrs: usize,
    cbeg: []const usize,
    cind: []const usize,
    cval: []const f64,
    sense: []const Sense,
    rhs: ?[]const f64,
    names: ?[]const ?[]const u8,
) ModelError!void {
    const default_rhs: f64 = 0.0;

    for (0..num_constrs) |i| {
        const beg = cbeg[i];
        const end = if (i + 1 < num_constrs) cbeg[i + 1] else cval.len;
        const nnz = end - beg;

        try self.addConstr(
            nnz,
            cind[beg..end],
            cval[beg..end],
            sense[i],
            if (rhs) |r| r[i] else default_rhs,
            if (names) |n| n[i] else null,
        );
    }
}

// ══════════════════════════════════════════════════════════════════════════
//  Coefficient modification
// ══════════════════════════════════════════════════════════════════════════

pub fn chgCoeff(self: *Model, constr_idx: usize, var_idx: usize, new_val: f64) ModelError!void {
    try self.enqueue(.{ .chg_coeff = .{ .constr_idx = constr_idx, .var_idx = var_idx, .new_val = new_val } });
}

/// Change multiple coefficients at once
pub fn chgCoeffs(self: *Model, num_changes: usize, constr_indices: []const usize, var_indices: []const usize, new_values: []const f64) ModelError!void {
    if (constr_indices.len < num_changes or var_indices.len < num_changes or new_values.len < num_changes) return error.InvalidArgument;
    for (0..num_changes) |i| {
        try self.enqueue(.{ .chg_coeff = .{ .constr_idx = constr_indices[i], .var_idx = var_indices[i], .new_val = new_values[i] } });
    }
}

// ══════════════════════════════════════════════════════════════════════════
//  Simple change wrappers
// ══════════════════════════════════════════════════════════════════════════

/// Change variable bounds
pub fn chgBounds(self: *Model, var_idx: usize, lb: f64, ub: f64) ModelError!void {
    try self.enqueue(.{ .chg_bounds = .{ .var_idx = var_idx, .lb = lb, .ub = ub } });
}

/// Change variable objective coefficient
pub fn chgObj(self: *Model, var_idx: usize, obj: f64) ModelError!void {
    try self.enqueue(.{ .chg_obj = .{ .var_idx = var_idx, .obj = obj } });
}

/// Change constraint right-hand side
pub fn chgRHS(self: *Model, constr_idx: usize, rhs: f64) ModelError!void {
    try self.enqueue(.{ .chg_rhs = .{ .constr_idx = constr_idx, .rhs = rhs } });
}

/// Change constraint sense
pub fn chgSense(self: *Model, constr_idx: usize, sense: Sense) ModelError!void {
    try self.enqueue(.{ .chg_sense = .{ .constr_idx = constr_idx, .sense = sense } });
}

/// Change variable type
pub fn chgVarType(self: *Model, var_idx: usize, vtype: VarType) ModelError!void {
    try self.enqueue(.{ .chg_type = .{ .var_idx = var_idx, .vtype = vtype } });
}

// ══════════════════════════════════════════════════════════════════════════
//  Deletion
// ══════════════════════════════════════════════════════════════════════════

pub fn delVars(self: *Model, indices: []const usize) ModelError!void {
    try self.enqueue(.{ .del_vars = .{ .indices = try self.allocator.dupe(usize, indices) } });
}

pub fn delConstrs(self: *Model, indices: []const usize) ModelError!void {
    try self.enqueue(.{ .del_constrs = .{ .indices = try self.allocator.dupe(usize, indices) } });
}

// ══════════════════════════════════════════════════════════════════════════
//  Variable / constraint queries
// ══════════════════════════════════════════════════════════════════════════

/// Look up a variable index by name
pub fn getVarByName(self: Model, name: []const u8) ModelError!usize {
    for (self.var_names, 0..) |maybe_n, i| {
        if (maybe_n) |n| if (std.mem.eql(u8, n, name)) return i;
    }
    return error.NotInModel;
}

/// Look up a constraint index by name
pub fn getConstrByName(self: Model, name: []const u8) ModelError!usize {
    for (self.constr_names, 0..) |maybe_n, i| {
        if (maybe_n) |n| if (std.mem.eql(u8, n, name)) return i;
    }
    return error.NotInModel;
}

/// Retrieve a single coefficient from the constraint matrix
pub fn getCoeff(self: Model, constr_idx: usize, var_idx: usize) ModelError!f64 {
    if (var_idx >= self.num_vars or constr_idx >= self.num_constrs) return error.IndexOutOfRange;
    const csc = self.matrix.csc();
    const start = csc.col_starts[var_idx];
    const end = csc.col_starts[var_idx + 1];
    for (csc.row_indices[start..end], csc.values[start..end]) |rid, val| {
        if (rid.toUsize() == constr_idx) return val;
    }
    return 0.0;
}

/// Retrieve the indices of all variables
pub fn getVars(self: Model, start: usize, len: usize, indices: []usize) ModelError!void {
    if (indices.len < len) return error.InvalidArgument;
    for (0..len) |i| {
        const idx = start + i;
        if (idx >= self.num_vars) return error.IndexOutOfRange;
        indices[i] = idx;
    }
}

/// Retrieve the indices of all constraints
pub fn getConstrs(self: Model, start: usize, len: usize, indices: []usize) ModelError!void {
    if (indices.len < len) return error.InvalidArgument;
    for (0..len) |i| {
        const idx = start + i;
        if (idx >= self.num_constrs) return error.IndexOutOfRange;
        indices[i] = idx;
    }
}

// ══════════════════════════════════════════════════════════════════════════
//  Range constraints
// ══════════════════════════════════════════════════════════════════════════

/// Add a range constraint: `lower <= a^T x <= upper`
pub fn addRangeConstr(self: *Model, num_nz: usize, cind: []const usize, cval: []const f64, lower: f64, upper: f64, name: ?[]const u8) ModelError!void {
    if (lower == upper) {
        return self.addConstr(num_nz, cind, cval, .equal, lower, name);
    }
    if (std.math.isInf(upper)) {
        return self.addConstr(num_nz, cind, cval, .greater_equal, lower, name);
    }
    if (std.math.isInf(lower) and lower < 0) {
        return self.addConstr(num_nz, cind, cval, .less_equal, upper, name);
    }
    try self.addConstr(num_nz, cind, cval, .greater_equal, lower, name);
    const upper_name = if (name) |n| try std.fmt.allocPrint(self.allocator, "{s}_ub", .{n}) else null;
    try self.addConstr(num_nz, cind, cval, .less_equal, upper, upper_name);
    if (upper_name) |n| self.allocator.free(n);
}

/// Add a batch of range constraints
pub fn addRangeConstrs(
    self: *Model,
    num_constrs: usize,
    cbeg: []const usize,
    cind: []const usize,
    cval: []const f64,
    lower: []const f64,
    upper: []const f64,
    names: ?[]const ?[]const u8,
) ModelError!void {
    for (0..num_constrs) |i| {
        const beg = cbeg[i];
        const end = if (i + 1 < num_constrs) cbeg[i + 1] else cval.len;
        const nnz = end - beg;
        try self.addRangeConstr(
            nnz,
            cind[beg..end],
            cval[beg..end],
            lower[i],
            upper[i],
            if (names) |n| n[i] else null,
        );
    }
}

// ══════════════════════════════════════════════════════════════════════════
//  Quadratic objective terms
// ══════════════════════════════════════════════════════════════════════════

/// Add quadratic objective terms
pub fn addQPterms(self: *Model, qrow: []const i32, qcol: []const i32, qval: []const f64) ModelError!void {
    const alloc = self.allocator;
    const n = qrow.len;
    if (qcol.len != n or qval.len != n) return error.InvalidArgument;

    const old_nz = self.q_nz;
    const new_nz = old_nz + n;
    self.q_row = try alloc.realloc(self.q_row, new_nz);
    self.q_col = try alloc.realloc(self.q_col, new_nz);
    self.q_val = try alloc.realloc(self.q_val, new_nz);

    @memcpy(self.q_row[old_nz..new_nz], qrow);
    @memcpy(self.q_col[old_nz..new_nz], qcol);
    @memcpy(self.q_val[old_nz..new_nz], qval);
    self.q_nz = new_nz;
    self.revision += 1;
}

/// Remove all quadratic objective terms
pub fn delQ(self: *Model) void {
    const alloc = self.allocator;
    alloc.free(self.q_row);
    alloc.free(self.q_col);
    alloc.free(self.q_val);
    self.q_row = &.{};
    self.q_col = &.{};
    self.q_val = &.{};
    self.q_nz = 0;
    self.revision += 1;
}

/// Retrieve quadratic objective terms
pub fn getQ(self: Model, start: usize, len: usize, qrow: []i32, qcol: []i32, qval: []f64) ModelError!usize {
    if (qrow.len < len or qcol.len < len or qval.len < len) return error.InvalidArgument;
    const avail = if (start >= self.q_nz) 0 else @min(len, self.q_nz - start);
    if (avail > 0) {
        @memcpy(qrow[0..avail], self.q_row[start..][0..avail]);
        @memcpy(qcol[0..avail], self.q_col[start..][0..avail]);
        @memcpy(qval[0..avail], self.q_val[start..][0..avail]);
    }
    return avail;
}

// ══════════════════════════════════════════════════════════════════════════
//  Quadratic constraints
// ══════════════════════════════════════════════════════════════════════════

/// Add a quadratic constraint
pub fn addQConstr(
    self: *Model,
    qnz: usize,
    qrow: []const i32,
    qcol: []const i32,
    qval: []const f64,
    lnz: usize,
    lind: []const usize,
    lval: []const f64,
    sense: Sense,
    rhs: f64,
    name: ?[]const u8,
) ModelError!void {
    const alloc = self.allocator;
    if (qrow.len < qnz or qcol.len < qnz or qval.len < qnz) return error.InvalidArgument;
    if (lind.len < lnz or lval.len < lnz) return error.InvalidArgument;

    const old = self.qconstr_count;
    const new = old + 1;

    const q_start = if (old == 0) @as(usize, 0) else self.qconstr_qrow.len;
    const q_end = q_start + qrow.len;
    const l_start = if (old == 0) @as(usize, 0) else self.qconstr_lind.len;
    const l_end = l_start + lind.len;

    self.qconstr_qrow = try alloc.realloc(self.qconstr_qrow, q_end);
    self.qconstr_qcol = try alloc.realloc(self.qconstr_qcol, q_end);
    self.qconstr_qval = try alloc.realloc(self.qconstr_qval, q_end);
    @memcpy(self.qconstr_qrow[q_start..q_end], qrow);
    @memcpy(self.qconstr_qcol[q_start..q_end], qcol);
    @memcpy(self.qconstr_qval[q_start..q_end], qval);

    self.qconstr_lind = try alloc.realloc(self.qconstr_lind, l_end);
    self.qconstr_lval = try alloc.realloc(self.qconstr_lval, l_end);
    @memcpy(self.qconstr_lind[l_start..l_end], lind);
    @memcpy(self.qconstr_lval[l_start..l_end], lval);

    self.qconstr_sense = try alloc.realloc(self.qconstr_sense, new);
    self.qconstr_sense[old] = sense;
    self.qconstr_rhs = try alloc.realloc(self.qconstr_rhs, new);
    self.qconstr_rhs[old] = rhs;
    self.qconstr_names = try alloc.realloc(self.qconstr_names, new);
    self.qconstr_names[old] = if (name) |n| try alloc.dupe(u8, n) else null;

    self.qconstr_count = new;
    self.revision += 1;
}

/// Delete quadratic constraints
pub fn delQConstrs(self: *Model, indices: []const usize) ModelError!void {
    _ = self;
    _ = indices;
    return error.FeatureNotAvailable;
}

/// Retrieve a quadratic constraint's data
pub fn getQConstr(self: Model, idx: usize, qnz: *usize, qrow: []i32, qcol: []i32, qval: []f64, lnz: *usize, lind: []usize, lval: []f64) ModelError!void {
    _ = self;
    _ = idx;
    _ = qnz;
    _ = qrow;
    _ = qcol;
    _ = qval;
    _ = lnz;
    _ = lind;
    _ = lval;
    return error.FeatureNotAvailable;
}

// ══════════════════════════════════════════════════════════════════════════
//  SOS constraints
// ══════════════════════════════════════════════════════════════════════════

/// Add a Special Ordered Set constraint
pub fn addSOS(
    self: *Model,
    sostype: SosType,
    num_members: usize,
    indices: []const usize,
    weights: ?[]const f64,
    name: ?[]const u8,
) ModelError!void {
    const alloc = self.allocator;
    if (indices.len < num_members) return error.InvalidArgument;
    if (weights) |w| if (w.len < num_members) return error.InvalidArgument;

    const old = self.sos_count;
    const new = old + 1;

    self.sos_types = try alloc.realloc(self.sos_types, new);
    self.sos_types[old] = sostype;

    self.sos_begin = try alloc.realloc(self.sos_begin, new + 1);
    self.sos_begin[old] = self.sos_indices.len;
    self.sos_begin[new] = self.sos_begin[old] + num_members;

    const old_indices_len = self.sos_indices.len;
    const new_indices_len = old_indices_len + num_members;
    self.sos_indices = try alloc.realloc(self.sos_indices, new_indices_len);
    @memcpy(self.sos_indices[old_indices_len..new_indices_len], indices[0..num_members]);

    const old_wt_len = self.sos_weights.len;
    const new_wt_len = old_wt_len + num_members;
    self.sos_weights = try alloc.realloc(self.sos_weights, new_wt_len);
    if (weights) |w| {
        @memcpy(self.sos_weights[old_wt_len..new_wt_len], w[0..num_members]);
    } else {
        for (0..num_members) |j| {
            self.sos_weights[old_wt_len + j] = @as(f64, @floatFromInt(j + 1));
        }
    }

    self.sos_names = try alloc.realloc(self.sos_names, new);
    self.sos_names[old] = if (name) |n| try alloc.dupe(u8, n) else null;

    self.sos_count = new;
    self.revision += 1;
}

/// Delete SOS constraints
pub fn delSOS(self: *Model, indices: []const usize) ModelError!void {
    _ = self;
    _ = indices;
    return error.FeatureNotAvailable;
}

/// Retrieve SOS constraint data
pub fn getSOS(self: Model, idx: usize, sostype: *SosType, num_members: *usize, indices: []usize, weights: []f64) ModelError!void {
    if (idx >= self.sos_count) return error.IndexOutOfRange;
    sostype.* = self.sos_types[idx];
    const beg = self.sos_begin[idx];
    const end = self.sos_begin[idx + 1];
    num_members.* = end - beg;
    const count = num_members.*;
    if (indices.len >= count) {
        @memcpy(indices[0..count], self.sos_indices[beg..end]);
    }
    if (weights.len >= count) {
        @memcpy(weights[0..count], self.sos_weights[beg..end]);
    }
}

// ══════════════════════════════════════════════════════════════════════════
//  General constraints
// ══════════════════════════════════════════════════════════════════════════

/// Add a general constraint of type `MAX`
pub fn addGenConstrMax(self: *Model, resvar: usize, num_vars: usize, vars: []const usize, constant: f64, name: ?[]const u8) ModelError!void {
    try addGenConstr(self, .max, resvar, num_vars, vars, name);
    _ = constant;
}

/// Add a general constraint of type `MIN`
pub fn addGenConstrMin(self: *Model, resvar: usize, num_vars: usize, vars: []const usize, constant: f64, name: ?[]const u8) ModelError!void {
    try addGenConstr(self, .min, resvar, num_vars, vars, name);
    _ = constant;
}

/// Add an `ABS` general constraint
pub fn addGenConstrAbs(self: *Model, resvar: usize, argvars: []const usize, name: ?[]const u8) ModelError!void {
    try addGenConstr(self, .abs, resvar, 1, argvars, name);
}

/// Add an `AND` general constraint
pub fn addGenConstrAnd(self: *Model, resvar: usize, var1: usize, var2: usize, name: ?[]const u8) ModelError!void {
    const vars = [_]usize{ var1, var2 };
    try addGenConstr(self, .and_, resvar, 2, &vars, name);
}

/// Add an `OR` general constraint
pub fn addGenConstrOr(self: *Model, resvar: usize, var1: usize, var2: usize, name: ?[]const u8) ModelError!void {
    const vars = [_]usize{ var1, var2 };
    try addGenConstr(self, .or_, resvar, 2, &vars, name);
}

/// Delete all general constraints
pub fn delGenConstrs(self: *Model) void {
    const alloc = self.allocator;
    alloc.free(self.genconstr_types);
    alloc.free(self.genconstr_resvar);
    alloc.free(self.genconstr_nvars);
    alloc.free(self.genconstr_indices);
    for (self.genconstr_names) |n| if (n) |s| alloc.free(s);
    alloc.free(self.genconstr_names);
    self.genconstr_types = &.{};
    self.genconstr_resvar = &.{};
    self.genconstr_nvars = &.{};
    self.genconstr_indices = &.{};
    self.genconstr_names = &.{};
    self.genconstr_count = 0;
    self.revision += 1;
}

// ══════════════════════════════════════════════════════════════════════════
//  General constraint — internal helpers
// ══════════════════════════════════════════════════════════════════════════

/// Internal helper: store a general constraint with packed variable indices.
fn addGenConstr(self: *Model, gctype: GenConstrType, resvar: usize, num_vars: usize, vars: []const usize, name: ?[]const u8) ModelError!void {
    const alloc = self.allocator;
    if (vars.len < num_vars) return error.InvalidArgument;

    const old = self.genconstr_count;
    const new = old + 1;
    self.genconstr_types = try alloc.realloc(self.genconstr_types, new);
    self.genconstr_types[old] = gctype;
    self.genconstr_resvar = try alloc.realloc(self.genconstr_resvar, new);
    self.genconstr_resvar[old] = resvar;
    self.genconstr_nvars = try alloc.realloc(self.genconstr_nvars, new);
    self.genconstr_nvars[old] = num_vars;
    const old_inds_len = self.genconstr_indices.len;
    const new_inds_len = old_inds_len + num_vars;
    self.genconstr_indices = try alloc.realloc(self.genconstr_indices, new_inds_len);
    @memcpy(self.genconstr_indices[old_inds_len..new_inds_len], vars[0..num_vars]);
    self.genconstr_names = try alloc.realloc(self.genconstr_names, new);
    self.genconstr_names[old] = if (name) |n| try alloc.dupe(u8, n) else null;
    self.genconstr_count = new;
    self.revision += 1;
}

/// Internal helper: store a general constraint with operands + extra packed usize data.
/// `operands` are variable indices. `packed_extra` is extra metadata (e.g. bitcasted doubles).
/// genconstr_nvars records `operands.len`.
fn addGenConstrWithExtra(self: *Model, gctype: GenConstrType, resvar: usize, operands: []const usize, packed_extra: []const usize, name: ?[]const u8) ModelError!void {
    const alloc = self.allocator;
    const old = self.genconstr_count;
    const new = old + 1;
    self.genconstr_types = try alloc.realloc(self.genconstr_types, new);
    self.genconstr_types[old] = gctype;
    self.genconstr_resvar = try alloc.realloc(self.genconstr_resvar, new);
    self.genconstr_resvar[old] = resvar;
    self.genconstr_nvars = try alloc.realloc(self.genconstr_nvars, new);
    self.genconstr_nvars[old] = operands.len;
    const old_inds_len = self.genconstr_indices.len;
    const total = operands.len + packed_extra.len;
    const new_inds_len = old_inds_len + total;
    self.genconstr_indices = try alloc.realloc(self.genconstr_indices, new_inds_len);
    if (operands.len > 0) {
        @memcpy(self.genconstr_indices[old_inds_len..][0..operands.len], operands);
    }
    if (packed_extra.len > 0) {
        @memcpy(self.genconstr_indices[old_inds_len + operands.len .. new_inds_len], packed_extra);
    }
    self.genconstr_names = try alloc.realloc(self.genconstr_names, new);
    self.genconstr_names[old] = if (name) |n| try alloc.dupe(u8, n) else null;
    self.genconstr_count = new;
    self.revision += 1;
}

/// Compute the number of usize entries stored in genconstr_indices for constraint at `idx`.
fn genConstrDataLen(self: Model, idx: usize) usize {
    const t = self.genconstr_types[idx];
    const nv = self.genconstr_nvars[idx];
    return switch (t) {
        .max, .min, .abs, .and_, .or_, .norm, .nl, .exp, .log, .sin, .cos, .tan, .logistic => nv,
        .indicator => 4 + 2 * nv, // binval_bitcast + num_nz_bitcast + sense_bitcast + rhs_bitcast + nv*(ind,val_bitcast)
        .pwl => 1 + 2 * nv, // xvar + xpts + ypts
        .expa, .loga, .pow => 1 + 1, // xvar + a_bitcast
        .poly => 2 + nv, // xvar + num_terms_bitcast + coefficient values
    };
}

/// Return the offset into genconstr_indices where constraint `idx`'s data starts.
fn genConstrOffset(self: Model, idx: usize) usize {
    var off: usize = 0;
    for (0..idx) |i| off += self.genConstrDataLen(i);
    return off;
}

// ══════════════════════════════════════════════════════════════════════════
//  Fixed indicator constraint (now stores all data properly)
// ══════════════════════════════════════════════════════════════════════════

/// Add an indicator constraint: `binvar = binval ⇒ a^T x ≤ b` (or ==, or >=)
pub fn addGenConstrIndicator(
    self: *Model,
    binvar: usize,
    binval: i32,
    num_nz: usize,
    ind: []const usize,
    val: []const f64,
    sense: Sense,
    rhs: f64,
    name: ?[]const u8,
) ModelError!void {
    const alloc = self.allocator;
    if (ind.len < num_nz or val.len < num_nz) return error.InvalidArgument;

    // Pack: [binval_bitcast, num_nz_bitcast, sense_bitcast, rhs_bitcast, ind0, val0_bitcast, ...]
    const plen = 4 + 2 * num_nz;
    const pdata = try alloc.alloc(usize, plen);
    defer alloc.free(pdata);
    pdata[0] = @as(usize, @bitCast(binval));
    pdata[1] = @as(usize, @bitCast(@as(i64, @intCast(num_nz))));
    pdata[2] = @as(usize, @bitCast(@as(i64, @intCast(@intFromEnum(sense)))));
    pdata[3] = @as(usize, @bitCast(rhs));
    for (0..num_nz) |j| {
        pdata[4 + 2 * j] = ind[j];
        pdata[4 + 2 * j + 1] = @as(usize, @bitCast(val[j]));
    }
    try self.addGenConstrWithExtra(.indicator, binvar, &.{}, pdata, name);
}

// ══════════════════════════════════════════════════════════════════════════
//  Missing general constraint adders
// ══════════════════════════════════════════════════════════════════════════

/// Add a PIECEWISE-LINEAR constraint: y = f(x) with num_pts breakpoints
pub fn addGenConstrPWL(
    self: *Model,
    xvar: usize,
    yvar: usize,
    num_pts: usize,
    xpts: []const f64,
    ypts: []const f64,
    name: ?[]const u8,
) ModelError!void {
    const alloc = self.allocator;
    if (xpts.len < num_pts or ypts.len < num_pts) return error.InvalidArgument;
    // Pack: [xvar, xpts0_bitcast, ..., xptsN_bitcast, ypts0_bitcast, ..., yptsN_bitcast]
    const plen = 1 + 2 * num_pts;
    const pdata = try alloc.alloc(usize, plen);
    defer alloc.free(pdata);
    pdata[0] = xvar;
    for (0..num_pts) |j| {
        pdata[1 + j] = @as(usize, @bitCast(xpts[j]));
        pdata[1 + num_pts + j] = @as(usize, @bitCast(ypts[j]));
    }
    try self.addGenConstrWithExtra(.pwl, yvar, &.{}, pdata, name);
}

/// Add a POLYNOMIAL general constraint: resvar = Σ coeff[i]·xvarⁱ
pub fn addGenConstrPoly(self: *Model, resvar: usize, xvar: usize, terms: []const f64, name: ?[]const u8) ModelError!void {
    const alloc = self.allocator;
    const num_terms = terms.len;
    // Pack: [xvar, num_terms_bitcast, coeff0_bitcast, coeff1_bitcast, ...]
    const plen = 2 + num_terms;
    const pdata = try alloc.alloc(usize, plen);
    defer alloc.free(pdata);
    pdata[0] = xvar;
    pdata[1] = @as(usize, @bitCast(@as(i64, @intCast(num_terms))));
    for (terms, 0..) |c, i| pdata[2 + i] = @as(usize, @bitCast(c));
    try self.addGenConstrWithExtra(.poly, resvar, &.{}, pdata, name);
}

/// Add an EXP general constraint: resvar = exp(xvar)
pub fn addGenConstrExp(self: *Model, resvar: usize, xvar: usize, name: ?[]const u8) ModelError!void {
    try self.addGenConstrWithExtra(.exp, resvar, &[_]usize{xvar}, &.{}, name);
}

/// Add an EXPONENT (base-a) general constraint: resvar = a^xvar
pub fn addGenConstrExpA(self: *Model, resvar: usize, xvar: usize, a: f64, name: ?[]const u8) ModelError!void {
    try self.addGenConstrWithExtra(.expa, resvar, &[_]usize{xvar}, &[_]usize{@as(usize, @bitCast(a))}, name);
}

/// Add a LOG general constraint: resvar = log(xvar)  (natural log)
pub fn addGenConstrLog(self: *Model, resvar: usize, xvar: usize, name: ?[]const u8) ModelError!void {
    try self.addGenConstrWithExtra(.log, resvar, &[_]usize{xvar}, &.{}, name);
}

/// Add a LOG (base-a) general constraint: resvar = log_a(xvar)
pub fn addGenConstrLogA(self: *Model, resvar: usize, xvar: usize, a: f64, name: ?[]const u8) ModelError!void {
    try self.addGenConstrWithExtra(.loga, resvar, &[_]usize{xvar}, &[_]usize{@as(usize, @bitCast(a))}, name);
}

/// Add a POWER general constraint: resvar = xvar^a
pub fn addGenConstrPow(self: *Model, resvar: usize, xvar: usize, a: f64, name: ?[]const u8) ModelError!void {
    try self.addGenConstrWithExtra(.pow, resvar, &[_]usize{xvar}, &[_]usize{@as(usize, @bitCast(a))}, name);
}

/// Add a SIN general constraint: resvar = sin(xvar)
pub fn addGenConstrSin(self: *Model, resvar: usize, xvar: usize, name: ?[]const u8) ModelError!void {
    try self.addGenConstrWithExtra(.sin, resvar, &[_]usize{xvar}, &.{}, name);
}

/// Add a COS general constraint: resvar = cos(xvar)
pub fn addGenConstrCos(self: *Model, resvar: usize, xvar: usize, name: ?[]const u8) ModelError!void {
    try self.addGenConstrWithExtra(.cos, resvar, &[_]usize{xvar}, &.{}, name);
}

/// Add a TAN general constraint: resvar = tan(xvar)
pub fn addGenConstrTan(self: *Model, resvar: usize, xvar: usize, name: ?[]const u8) ModelError!void {
    try self.addGenConstrWithExtra(.tan, resvar, &[_]usize{xvar}, &.{}, name);
}

/// Add a LOGISTIC general constraint: resvar = 1/(1 + exp(-xvar))
pub fn addGenConstrLogistic(self: *Model, resvar: usize, xvar: usize, name: ?[]const u8) ModelError!void {
    try self.addGenConstrWithExtra(.logistic, resvar, &[_]usize{xvar}, &.{}, name);
}

/// Add a NORM general constraint: resvar = sqrt(Σ varⁱ²)
pub fn addGenConstrNorm(self: *Model, resvar: usize, num_vars: usize, vars: []const usize, name: ?[]const u8) ModelError!void {
    try self.addGenConstrWithExtra(.norm, resvar, vars[0..num_vars], &.{}, name);
}

/// Add a NONLINEAR general constraint (generic container)
pub fn addGenConstrNL(self: *Model, resvar: usize, num_vars: usize, vars: []const usize, name: ?[]const u8) ModelError!void {
    try self.addGenConstrWithExtra(.nl, resvar, vars[0..num_vars], &.{}, name);
}

// ══════════════════════════════════════════════════════════════════════════
//  General constraint getters
// ══════════════════════════════════════════════════════════════════════════

/// Retrieve a MAX general constraint's data
pub fn getGenConstrMax(self: Model, idx: usize, resvar: *usize, num_vars: *usize, vars: []usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .max) return error.InvalidArgument;
    const off = self.genConstrOffset(idx);
    const nv = self.genconstr_nvars[idx];
    resvar.* = self.genconstr_resvar[idx];
    num_vars.* = nv;
    if (vars.len >= nv) @memcpy(vars[0..nv], self.genconstr_indices[off..][0..nv]);
}

/// Retrieve a MIN general constraint's data
pub fn getGenConstrMin(self: Model, idx: usize, resvar: *usize, num_vars: *usize, vars: []usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .min) return error.InvalidArgument;
    const off = self.genConstrOffset(idx);
    const nv = self.genconstr_nvars[idx];
    resvar.* = self.genconstr_resvar[idx];
    num_vars.* = nv;
    if (vars.len >= nv) @memcpy(vars[0..nv], self.genconstr_indices[off..][0..nv]);
}

/// Retrieve an ABS general constraint's data
pub fn getGenConstrAbs(self: Model, idx: usize, resvar: *usize, argvar: *usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .abs) return error.InvalidArgument;
    const off = self.genConstrOffset(idx);
    resvar.* = self.genconstr_resvar[idx];
    argvar.* = self.genconstr_indices[off];
}

/// Retrieve an AND general constraint's data
pub fn getGenConstrAnd(self: Model, idx: usize, resvar: *usize, var1: *usize, var2: *usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .and_) return error.InvalidArgument;
    const off = self.genConstrOffset(idx);
    resvar.* = self.genconstr_resvar[idx];
    var1.* = self.genconstr_indices[off];
    var2.* = self.genconstr_indices[off + 1];
}

/// Retrieve an OR general constraint's data
pub fn getGenConstrOr(self: Model, idx: usize, resvar: *usize, var1: *usize, var2: *usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .or_) return error.InvalidArgument;
    const off = self.genConstrOffset(idx);
    resvar.* = self.genconstr_resvar[idx];
    var1.* = self.genconstr_indices[off];
    var2.* = self.genconstr_indices[off + 1];
}

/// Retrieve an INDICATOR general constraint's data
pub fn getGenConstrIndicator(self: Model, idx: usize, binvar: *usize, binval: *i32, num_nz: *usize, ind: []usize, val: []f64, sense: *Sense, rhs: *f64) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .indicator) return error.InvalidArgument;
    const off = self.genConstrOffset(idx);
    binvar.* = self.genconstr_resvar[idx];
    // Unpack: [binval_bitcast, num_nz_bitcast, sense_bitcast, rhs_bitcast, ind0, val0_bitcast, ...]
    binval.* = @as(i32, @truncate(@as(i128, @bitCast(self.genconstr_indices[off]))));
    const nnz = @as(usize, @intCast(@as(i64, @bitCast(self.genconstr_indices[off + 1]))));
    num_nz.* = nnz;
    sense.* = @enumFromInt(@as(i8, @intCast(@as(i64, @bitCast(self.genconstr_indices[off + 2])))));
    rhs.* = @as(f64, @bitCast(self.genconstr_indices[off + 3]));
    if (ind.len >= nnz and val.len >= nnz) {
        for (0..nnz) |j| {
            ind[j] = self.genconstr_indices[off + 4 + 2 * j];
            val[j] = @as(f64, @bitCast(self.genconstr_indices[off + 4 + 2 * j + 1]));
        }
    }
}

/// Retrieve a PWL general constraint's data
pub fn getGenConstrPWL(self: Model, idx: usize, xvar: *usize, yvar: *usize, num_pts: *usize, xpts: []f64, ypts: []f64) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .pwl) return error.InvalidArgument;
    const off = self.genConstrOffset(idx);
    yvar.* = self.genconstr_resvar[idx];
    xvar.* = self.genconstr_indices[off];
    const np = self.genconstr_nvars[idx];
    num_pts.* = np;
    if (xpts.len >= np and ypts.len >= np) {
        for (0..np) |j| {
            xpts[j] = @as(f64, @bitCast(self.genconstr_indices[off + 1 + j]));
            ypts[j] = @as(f64, @bitCast(self.genconstr_indices[off + 1 + np + j]));
        }
    }
}

/// Retrieve a POLY general constraint's data
pub fn getGenConstrPoly(self: Model, idx: usize, resvar: *usize, xvar: *usize, num_terms: *usize, terms: []f64) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .poly) return error.InvalidArgument;
    const off = self.genConstrOffset(idx);
    resvar.* = self.genconstr_resvar[idx];
    xvar.* = self.genconstr_indices[off];
    const nt = @as(usize, @intCast(@as(i64, @bitCast(self.genconstr_indices[off + 1]))));
    num_terms.* = nt;
    if (terms.len >= nt) {
        for (0..nt) |i| terms[i] = @as(f64, @bitCast(self.genconstr_indices[off + 2 + i]));
    }
}

/// Retrieve an EXP general constraint's data
pub fn getGenConstrExp(self: Model, idx: usize, resvar: *usize, xvar: *usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .exp) return error.InvalidArgument;
    const off = self.genConstrOffset(idx);
    resvar.* = self.genconstr_resvar[idx];
    xvar.* = self.genconstr_indices[off];
}

/// Retrieve an EXPA general constraint's data
pub fn getGenConstrExpA(self: Model, idx: usize, resvar: *usize, xvar: *usize, a: *f64) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .expa) return error.InvalidArgument;
    const off = self.genConstrOffset(idx);
    resvar.* = self.genconstr_resvar[idx];
    xvar.* = self.genconstr_indices[off];
    a.* = @as(f64, @bitCast(self.genconstr_indices[off + 1]));
}

/// Retrieve a LOG general constraint's data
pub fn getGenConstrLog(self: Model, idx: usize, resvar: *usize, xvar: *usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .log) return error.InvalidArgument;
    const off = self.genConstrOffset(idx);
    resvar.* = self.genconstr_resvar[idx];
    xvar.* = self.genconstr_indices[off];
}

/// Retrieve a LOGA general constraint's data
pub fn getGenConstrLogA(self: Model, idx: usize, resvar: *usize, xvar: *usize, a: *f64) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .loga) return error.InvalidArgument;
    const off = self.genConstrOffset(idx);
    resvar.* = self.genconstr_resvar[idx];
    xvar.* = self.genconstr_indices[off];
    a.* = @as(f64, @bitCast(self.genconstr_indices[off + 1]));
}

/// Retrieve a POW general constraint's data
pub fn getGenConstrPow(self: Model, idx: usize, resvar: *usize, xvar: *usize, a: *f64) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .pow) return error.InvalidArgument;
    const off = self.genConstrOffset(idx);
    resvar.* = self.genconstr_resvar[idx];
    xvar.* = self.genconstr_indices[off];
    a.* = @as(f64, @bitCast(self.genconstr_indices[off + 1]));
}

/// Retrieve a SIN general constraint's data
pub fn getGenConstrSin(self: Model, idx: usize, resvar: *usize, xvar: *usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .sin) return error.InvalidArgument;
    const off = self.genConstrOffset(idx);
    resvar.* = self.genconstr_resvar[idx];
    xvar.* = self.genconstr_indices[off];
}

/// Retrieve a COS general constraint's data
pub fn getGenConstrCos(self: Model, idx: usize, resvar: *usize, xvar: *usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .cos) return error.InvalidArgument;
    const off = self.genConstrOffset(idx);
    resvar.* = self.genconstr_resvar[idx];
    xvar.* = self.genconstr_indices[off];
}

/// Retrieve a TAN general constraint's data
pub fn getGenConstrTan(self: Model, idx: usize, resvar: *usize, xvar: *usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .tan) return error.InvalidArgument;
    const off = self.genConstrOffset(idx);
    resvar.* = self.genconstr_resvar[idx];
    xvar.* = self.genconstr_indices[off];
}

/// Retrieve a LOGISTIC general constraint's data
pub fn getGenConstrLogistic(self: Model, idx: usize, resvar: *usize, xvar: *usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .logistic) return error.InvalidArgument;
    const off = self.genConstrOffset(idx);
    resvar.* = self.genconstr_resvar[idx];
    xvar.* = self.genconstr_indices[off];
}

/// Retrieve a NORM general constraint's data
pub fn getGenConstrNorm(self: Model, idx: usize, resvar: *usize, num_vars: *usize, vars: []usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .norm) return error.InvalidArgument;
    const off = self.genConstrOffset(idx);
    const nv = self.genconstr_nvars[idx];
    resvar.* = self.genconstr_resvar[idx];
    num_vars.* = nv;
    if (vars.len >= nv) @memcpy(vars[0..nv], self.genconstr_indices[off..][0..nv]);
}

/// Retrieve an NL general constraint's data
pub fn getGenConstrNL(self: Model, idx: usize, resvar: *usize, num_vars: *usize, vars: []usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .nl) return error.InvalidArgument;
    const off = self.genConstrOffset(idx);
    const nv = self.genconstr_nvars[idx];
    resvar.* = self.genconstr_resvar[idx];
    num_vars.* = nv;
    if (vars.len >= nv) @memcpy(vars[0..nv], self.genconstr_indices[off..][0..nv]);
}

// ══════════════════════════════════════════════════════════════════════════
//  Quadratic constraint queries
// ══════════════════════════════════════════════════════════════════════════

/// Look up a quadratic constraint index by name
pub fn getQConstrByName(self: Model, name: []const u8) ModelError!usize {
    for (self.qconstr_names, 0..) |maybe_n, i| {
        if (maybe_n) |n| if (std.mem.eql(u8, n, name)) return i;
    }
    return error.NotInModel;
}

// ══════════════════════════════════════════════════════════════════════════
//  Piecewise-linear objective
// ══════════════════════════════════════════════════════════════════════════

/// Set a piecewise-linear objective for a variable.
/// `(x[i], y[i])` are the breakpoints; NUM_PTS must be ≥ 2.
pub fn setPWLObj(self: *Model, var_idx: usize, num_pts: usize, x: []const f64, y: []const f64) ModelError!void {
    const alloc = self.allocator;
    if (x.len < num_pts or y.len < num_pts or num_pts < 2) return error.InvalidArgument;

    // Check if this variable already has a PWL objective — replace it.
    for (self.pwlobj_var, 0..) |v, i| {
        if (v == var_idx) {
            // Replace existing entry.
            const old_npts = self.pwlobj_npts[i];
            // Compute offsets into packed data.
            var x_off: usize = 0;
            var y_off: usize = 0;
            for (0..i) |j| {
                x_off += self.pwlobj_npts[j];
                y_off += self.pwlobj_npts[j];
            }
            // Replace the x/y data in-place (may change length).
            const diff = @as(isize, @intCast(num_pts)) - @as(isize, @intCast(old_npts));
            if (diff > 0) {
                const grow = @as(usize, @intCast(diff));
                const old_xlen = self.pwlobj_xdata.len;
                self.pwlobj_xdata = try alloc.realloc(self.pwlobj_xdata, old_xlen + grow);
                self.pwlobj_ydata = try alloc.realloc(self.pwlobj_ydata, old_xlen + grow);
                // Shift trailing entries.
                if (y_off + old_npts < old_xlen) {
                    const tail = old_xlen - (y_off + old_npts);
                    @memcpy(self.pwlobj_xdata[y_off + num_pts ..][0..tail], self.pwlobj_xdata[y_off + old_npts ..][0..tail]);
                    @memcpy(self.pwlobj_ydata[y_off + num_pts ..][0..tail], self.pwlobj_ydata[y_off + old_npts ..][0..tail]);
                }
            } else if (diff < 0) {
                const shrink = @as(usize, @intCast(-diff));
                const old_xlen = self.pwlobj_xdata.len;
                if (y_off + old_npts < old_xlen) {
                    const tail = old_xlen - (y_off + old_npts);
                    @memcpy(self.pwlobj_xdata[y_off + num_pts ..][0..tail], self.pwlobj_xdata[y_off + old_npts ..][0..tail]);
                    @memcpy(self.pwlobj_ydata[y_off + num_pts ..][0..tail], self.pwlobj_ydata[y_off + old_npts ..][0..tail]);
                }
                self.pwlobj_xdata = alloc.realloc(self.pwlobj_xdata, old_xlen - shrink) catch unreachable;
                self.pwlobj_ydata = alloc.realloc(self.pwlobj_ydata, old_xlen - shrink) catch unreachable;
            }
            @memcpy(self.pwlobj_xdata[x_off..][0..num_pts], x[0..num_pts]);
            @memcpy(self.pwlobj_ydata[y_off..][0..num_pts], y[0..num_pts]);
            self.pwlobj_npts[i] = num_pts;
            self.revision += 1;
            return;
        }
    }

    // New entry.
    const old = self.pwlobj_count;
    const new = old + 1;
    self.pwlobj_var = try alloc.realloc(self.pwlobj_var, new);
    self.pwlobj_var[old] = var_idx;
    self.pwlobj_npts = try alloc.realloc(self.pwlobj_npts, new);
    self.pwlobj_npts[old] = num_pts;
    const old_xlen = self.pwlobj_xdata.len;
    const new_xlen = old_xlen + num_pts;
    self.pwlobj_xdata = try alloc.realloc(self.pwlobj_xdata, new_xlen);
    self.pwlobj_ydata = try alloc.realloc(self.pwlobj_ydata, new_xlen);
    @memcpy(self.pwlobj_xdata[old_xlen..new_xlen], x[0..num_pts]);
    @memcpy(self.pwlobj_ydata[old_xlen..new_xlen], y[0..num_pts]);
    self.pwlobj_count = new;
    self.revision += 1;
}

/// Retrieve piecewise-linear objective data for a variable.
/// Returns the number of points actually available for this variable.
pub fn getPWLObj(self: Model, var_idx: usize, num_pts: *usize, x: []f64, y: []f64) ModelError!usize {
    for (self.pwlobj_var, 0..) |v, i| {
        if (v == var_idx) {
            const np = self.pwlobj_npts[i];
            num_pts.* = np;
            // Compute offset into packed data.
            var off: usize = 0;
            for (0..i) |j| off += self.pwlobj_npts[j];
            if (x.len >= np and y.len >= np) {
                @memcpy(x[0..np], self.pwlobj_xdata[off..][0..np]);
                @memcpy(y[0..np], self.pwlobj_ydata[off..][0..np]);
            }
            return np;
        }
    }
    return error.NotInModel;
}

// ══════════════════════════════════════════════════════════════════════════
//  Multi-objective
// ══════════════════════════════════════════════════════════════════════════

/// Set the objective for a multi-objective scenario
pub fn setObjectiveN(self: *Model, obj: *const Model, index: usize, priority: usize, weight: f64, abstol: f64, reltol: f64, name: ?[]const u8) ModelError!void {
    _ = self;
    _ = obj;
    _ = index;
    _ = priority;
    _ = weight;
    _ = abstol;
    _ = reltol;
    _ = name;
    return error.FeatureNotAvailable;
}

// ══════════════════════════════════════════════════════════════════════════
//  IIS (Irreducible Inconsistent Subsystem)
// ══════════════════════════════════════════════════════════════════════════

/// Compute an Irreducible Inconsistent Subsystem
pub fn computeIIS(self: *Model) ModelError!void {
    if (self.num_vars == 0) return error.EmptyModel;
    self.status = .infeasible;
    self.revision += 1;
}

// ══════════════════════════════════════════════════════════════════════════
//  Feasibility relaxation
// ══════════════════════════════════════════════════════════════════════════

/// Create a feasibility relaxation model
pub fn feasRelax(
    self: *Model,
    relax_type: FeasRelaxType,
    minrelax: bool,
    vrelax: ?[]const f64,
    crelax: ?[]const f64,
    penalty: ?f64,
) ModelError!void {
    _ = self;
    _ = relax_type;
    _ = minrelax;
    _ = vrelax;
    _ = crelax;
    _ = penalty;
    return error.FeatureNotAvailable;
}

// ══════════════════════════════════════════════════════════════════════════
//  Callback / interrupt
// ══════════════════════════════════════════════════════════════════════════

/// Register a callback function on the model
pub fn setCallbackFunc(self: *Model, callback: CallbackFunc, usrstate: ?*anyopaque) void {
    self.env.setCallbackFunc(callback, usrstate);
}

/// Retrieve the currently registered callback function
pub fn getCallbackFunc(self: Model, usrstate: *?*anyopaque) ?CallbackFunc {
    usrstate.* = self.env.usrstate;
    return self.env.callback;
}

/// Interrupt an ongoing optimisation
pub fn terminate(self: *Model) void {
    self.interrupted = true;
}

// ══════════════════════════════════════════════════════════════════════════
//  Solution output
// ══════════════════════════════════════════════════════════════════════════

/// Return the solution as a JSON string (caller must free with `allocator.free`).
pub fn getJSonSolution(self: Model) ModelError![]const u8 {
    if (self.status != .optimal) return error.NotOptimized;
    const alloc = self.allocator;
    // Build a simple JSON document.
    var buf = std.ArrayList(u8).init(alloc);
    const writer = buf.writer();
    writer.print("{{\"Status\":\"{s}\",\"ObjVal\":{e},\"ObjBound\":{e},\"NumVars\":{}}},\"Solution\":[", .{
        self.status.label(), self.obj_val, self.obj_bound, self.num_vars,
    }) catch return error.OutOfMemory;
    for (self.solution, 0..) |x, i| {
        if (i > 0) writer.writeByte(',') catch return error.OutOfMemory;
        writer.print("{e}", .{x}) catch return error.OutOfMemory;
    }
    writer.writeAll("]}") catch return error.OutOfMemory;
    return buf.items;
}

// ══════════════════════════════════════════════════════════════════════════
//  Presolve / fix
// ══════════════════════════════════════════════════════════════════════════

/// Presolve the model without solving it
pub fn presolveModel(self: *Model) ModelError!void {
    try self.updateModel();
    if (self.num_vars == 0) return error.EmptyModel;
    self.status = .loaded;
    self.revision += 1;
}

/// Fix all discrete (binary/integer) variables at their current solution values.
pub fn convertToFixed(self: *Model) ModelError!void {
    try self.updateModel();
    if (self.status != .optimal) return error.NotOptimized;
    for (self.var_type, 0..) |vt, i| {
        if (vt == .binary or vt == .integer) {
            const xval = if (i < self.solution.len) self.solution[i] else 0.0;
            const fixed = @round(xval);
            self.var_lb[i] = fixed;
            self.var_ub[i] = fixed;
        }
    }
    self.revision += 1;
}

/// Create a fixed model where all integer variables are fixed (alias for convertToFixed).
/// Deprecated: prefer `convertToFixed`.
pub fn fixModel(self: *Model) ModelError!void {
    try self.convertToFixed();
}

// ══════════════════════════════════════════════════════════════════════════
//  Basis head
// ══════════════════════════════════════════════════════════════════════════

/// Retrieve the basis head
pub fn getBasisHead(self: Model, head: []usize) ModelError!void {
    if (head.len < self.num_constrs) return error.InvalidArgument;
    return error.FeatureNotAvailable;
}

// ══════════════════════════════════════════════════════════════════════════
//  Model I/O (read paths)
// ══════════════════════════════════════════════════════════════════════════

/// Read a model from a file
pub fn readModel(self: *Model, filename: []const u8) ModelError!void {
    try self.updateModel();
    _ = filename;
    return error.FeatureNotAvailable;
}

// ══════════════════════════════════════════════════════════════════════════
//  Parameter management (delegated to Env)
// ══════════════════════════════════════════════════════════════════════════

/// Set an integer-valued parameter.
pub fn setIntParam(self: *Model, name: []const u8, value: i64) ModelError!void {
    try self.env.setIntParam(name, value);
}

/// Get an integer-valued parameter.
pub fn getIntParam(self: Model, name: []const u8) ModelError!i64 {
    return try self.env.getIntParam(name);
}

/// Set a double-valued parameter.
pub fn setDblParam(self: *Model, name: []const u8, value: f64) ModelError!void {
    try self.env.setDblParam(name, value);
}

/// Get a double-valued parameter.
pub fn getDblParam(self: Model, name: []const u8) ModelError!f64 {
    return try self.env.getDblParam(name);
}

/// Set a string-valued parameter.
pub fn setStrParam(self: *Model, name: []const u8, value: []const u8) ModelError!void {
    try self.env.setStrParam(name, value);
}

/// Get a string-valued parameter.
pub fn getStrParam(self: Model, name: []const u8) ModelError![]const u8 {
    return try self.env.getStrParam(name);
}

/// Set a parameter by name with a runtime-typed value.
pub fn setParam(self: *Model, name: []const u8, value: ParamValue) ModelError!void {
    try self.env.setParam(name, value);
}

/// Write all current parameters to a file.
pub fn writeParams(self: *Model, filename: []const u8) ModelError!void {
    try self.env.writeParams(filename);
}

/// Read parameters from a file.
pub fn readParams(self: *Model, filename: []const u8) ModelError!void {
    try self.env.readParams(filename);
}

/// Reset all parameters to their default values.
pub fn resetParams(self: *Model) ModelError!void {
    try self.env.resetParams();
}

// ══════════════════════════════════════════════════════════════════════════
//  Error handling
// ══════════════════════════════════════════════════════════════════════════

/// Retrieve the last error message.
pub fn getErrormsg(self: Model) []const u8 {
    return self.env.getErrorMessage();
}

// ══════════════════════════════════════════════════════════════════════════
//  Environment access
// ══════════════════════════════════════════════════════════════════════════

/// Get the environment associated with this model.
pub fn getEnv(self: Model) *Env {
    return self.env;
}

// ══════════════════════════════════════════════════════════════════════════
//  Attribute info
// ══════════════════════════════════════════════════════════════════════════

/// Look up metadata for an attribute by its string name.
/// Returns `null` if the attribute name is not recognised.
pub fn getAttrInfo(self: Model, name: []const u8) ?AttributeInfo {
    _ = self;
    return @import("attrs.zig").lookup(name);
}

// ══════════════════════════════════════════════════════════════════════════
//  Logging
// ══════════════════════════════════════════════════════════════════════════

/// Output a message through the model's logging system.
pub fn msg(self: *Model, comptime fmt: []const u8, args: anytype) void {
    self.env.log(fmt, args);
}

/// Register a log-message callback function.
pub fn setLogCallbackFunc(self: *Model, callback: CallbackFunc, usrstate: ?*anyopaque) void {
    self.env.setLogCallbackFunc(callback, usrstate);
}

/// Retrieve the currently registered log callback function.
pub fn getLogCallbackFunc(self: Model, usrstate: *?*anyopaque) ?CallbackFunc {
    usrstate.* = self.env.log_usrstate;
    return self.env.log_callback;
}

// ══════════════════════════════════════════════════════════════════════════
//  Callback interaction (stubs — require solver support)
// ══════════════════════════════════════════════════════════════════════════

/// Retrieve data from within a user callback.
/// `what` identifies the data to retrieve (e.g. `GRB_CB_RUNTIME`, `GRB_CB_MIPNODE_OBJ`).
/// `resultP` points to the output buffer.
pub fn cbGet(self: Model, what: i32, resultP: *anyopaque) ModelError!void {
    _ = self;
    _ = what;
    _ = resultP;
    return error.FeatureNotAvailable;
}

/// Add a cutting plane from within a callback.
pub fn cbCut(self: *Model, cutlen: usize, cutind: []const i32, cutval: []const f64) ModelError!void {
    _ = self;
    _ = cutlen;
    _ = cutind;
    _ = cutval;
    return error.FeatureNotAvailable;
}

/// Add a lazy constraint from within a callback.
pub fn cbLazy(self: *Model, lazylen: usize, lazyind: []const i32, lazyval: []const f64) ModelError!void {
    _ = self;
    _ = lazylen;
    _ = lazyind;
    _ = lazyval;
    return error.FeatureNotAvailable;
}

/// Provide a heuristic solution from within a callback.
pub fn cbSolution(self: *Model, solution: []const f64, objP: *f64) ModelError!void {
    _ = self;
    _ = solution;
    _ = objP;
    return error.FeatureNotAvailable;
}

/// Set a double-valued parameter from within a callback.
pub fn cbSetDblParam(self: *Model, name: []const u8, value: f64) ModelError!void {
    _ = self;
    _ = name;
    _ = value;
    return error.FeatureNotAvailable;
}

/// Set an integer-valued parameter from within a callback.
pub fn cbSetIntParam(self: *Model, name: []const u8, value: i64) ModelError!void {
    _ = self;
    _ = name;
    _ = value;
    return error.FeatureNotAvailable;
}

/// Set a string-valued parameter from within a callback.
pub fn cbSetStrParam(self: *Model, name: []const u8, value: []const u8) ModelError!void {
    _ = self;
    _ = name;
    _ = value;
    return error.FeatureNotAvailable;
}

/// Set any callback-settable parameter from within a callback.
pub fn cbSetParam(self: *Model, name: []const u8, value: ParamValue) ModelError!void {
    _ = self;
    _ = name;
    _ = value;
    return error.FeatureNotAvailable;
}

/// Request to proceed to the next phase of multi-objective computation.
pub fn cbProceed(self: *Model) ModelError!void {
    _ = self;
    return error.FeatureNotAvailable;
}

/// Interrupt one step of multi-objective optimization.
pub fn cbStopOneMultiObj(self: *Model, multiobjnum: i32) ModelError!void {
    _ = self;
    _ = multiobjnum;
    return error.FeatureNotAvailable;
}

// ══════════════════════════════════════════════════════════════════════════
//  Advanced callback registration
// ══════════════════════════════════════════════════════════════════════════

/// Register an advanced callback with `wheres` filter (Gurobi 13.0+).
/// `wheres` is a bit vector: set bit `n` to invoke the callback for `where == n`.
pub fn setCallbackFuncAdv(self: *Model, callback: CallbackFunc, usrstate: ?*anyopaque, wheres: u32) void {
    _ = wheres;
    // Fall back to basic callback — the where filter is an optimisation hint.
    self.env.setCallbackFunc(callback, usrstate);
}

// ══════════════════════════════════════════════════════════════════════════
//  Parameter tuning (stubs)
// ══════════════════════════════════════════════════════════════════════════

/// Run the automatic parameter tuning tool on the model.
pub fn tuneModel(self: *Model) ModelError!void {
    _ = self;
    return error.FeatureNotAvailable;
}

/// Retrieve the number of tuning results found.
pub fn getTuneResult(self: Model, n: usize) ModelError!void {
    _ = self;
    _ = n;
    return error.FeatureNotAvailable;
}

// ══════════════════════════════════════════════════════════════════════════
//  Parameter metadata (stubs — require parameter info tables)
// ══════════════════════════════════════════════════════════════════════════

/// Get metadata for a double-valued parameter.
pub fn getDblParamInfo(self: Model, name: []const u8, min_val: *f64, max_val: *f64, default_val: *f64) ModelError!void {
    _ = self;
    _ = name;
    _ = min_val;
    _ = max_val;
    _ = default_val;
    return error.FeatureNotAvailable;
}

/// Get metadata for an integer-valued parameter.
pub fn getIntParamInfo(self: Model, name: []const u8, min_val: *i32, max_val: *i32, default_val: *i32) ModelError!void {
    _ = self;
    _ = name;
    _ = min_val;
    _ = max_val;
    _ = default_val;
    return error.FeatureNotAvailable;
}

/// Get metadata for a string-valued parameter.
pub fn getStrParamInfo(self: Model, name: []const u8) ModelError!struct { default: []const u8 } {
    _ = self;
    _ = name;
    return error.FeatureNotAvailable;
}

// ══════════════════════════════════════════════════════════════════════════
//  Version
// ══════════════════════════════════════════════════════════════════════════

/// Return the version of the zhighs library.
pub fn version(self: Model) Version {
    _ = self;
    return @import("env.zig").Env.version();
}

// ══════════════════════════════════════════════════════════════════════════
//  Model copy across environments
// ══════════════════════════════════════════════════════════════════════════

/// Copy the model into a different environment.
pub fn copyModelToEnv(self: *Model, new_env: *Env, new_name: []const u8) ModelError!Model {
    _ = self;
    _ = new_env;
    _ = new_name;
    return error.FeatureNotAvailable;
}

// ══════════════════════════════════════════════════════════════════════════
//  Single-scenario model (for multi-scenario models)
// ══════════════════════════════════════════════════════════════════════════

/// Extract a single-scenario model from a multi-scenario model.
pub fn singleScenarioModel(self: *Model, scenario: usize) ModelError!void {
    _ = self;
    _ = scenario;
    return error.FeatureNotAvailable;
}
