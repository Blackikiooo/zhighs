//! Linear variable and constraint construction, mutation, and query methods.
//!
//! ## Responsibility
//!
//! Owns variable and linear/range-constraint creation, queued coefficient and
//! bound/type/RHS changes, deletion requests, and basic linear-model queries.
//! Applying queued changes is the responsibility of `model_update.zig`.

const std = @import("std");
const types = @import("types.zig");
const Model = @import("model.zig").Model;
const RowId = @import("foundation").RowId;

const ModelError = types.ModelError;
const VarType = types.VarType;
const Sense = types.Sense;
const INFINITY = types.INFINITY;
const Column = @import("expr/column.zig").Column;

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

/// Add a variable from a sparse column. Stable constraint handles are
/// resolved once at the API boundary; the pending/model layers keep dense
/// indices for cache-friendly matrix storage.
pub fn addVarColumn(
    self: *Model,
    column: Column,
    obj: f64,
    lb: f64,
    ub: f64,
    vtype: VarType,
    name: ?[]const u8,
) ModelError!void {
    const indices = try column.resolveIndices(self.allocator, self.*);
    defer self.allocator.free(indices);
    return self.addVar(indices.len, indices, column.values.items, obj, lb, ub, vtype, name);
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
