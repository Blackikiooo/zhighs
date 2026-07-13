//! General constraint storage and API methods for `Model`.
//!
//! ## Responsibility
//!
//! Owns all general-constraint add/get/delete APIs and the internal packed-data
//! layout used by indicator, PWL, polynomial, transcendental, norm, and related
//! constraints.  Changes to packed encoding, length accounting, and offsets
//! must remain centralized here.

const std = @import("std");
const types = @import("types.zig");
const Model = @import("model.zig").Model;

const ModelError = types.ModelError;
const Sense = types.Sense;
const GenConstrType = types.GenConstrType;

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
    alloc.free(self.genconstr_begin);
    for (self.genconstr_names) |n| if (n) |s| alloc.free(s);
    alloc.free(self.genconstr_names);
    self.genconstr_types = &.{};
    self.genconstr_resvar = &.{};
    self.genconstr_nvars = &.{};
    self.genconstr_indices = &.{};
    self.genconstr_begin = &.{};
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
    if (old == 0) {
        self.genconstr_begin = try alloc.realloc(self.genconstr_begin, 1);
        self.genconstr_begin[0] = 0;
    }
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
    self.genconstr_begin = try alloc.realloc(self.genconstr_begin, new + 1);
    self.genconstr_begin[new] = new_inds_len;
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
    if (old == 0) {
        self.genconstr_begin = try alloc.realloc(self.genconstr_begin, 1);
        self.genconstr_begin[0] = 0;
    }
    self.genconstr_types = try alloc.realloc(self.genconstr_types, new);
    self.genconstr_types[old] = gctype;
    self.genconstr_resvar = try alloc.realloc(self.genconstr_resvar, new);
    self.genconstr_resvar[old] = resvar;
    self.genconstr_nvars = try alloc.realloc(self.genconstr_nvars, new);
    // For packed constraint kinds this field records the logical repeated-item
    // count used by genConstrDataLen (nonzeros, points, or coefficients).
    self.genconstr_nvars[old] = switch (gctype) {
        .indicator => (packed_extra.len - 4) / 2,
        .pwl => (packed_extra.len - 1) / 2,
        .poly => packed_extra.len - 2,
        else => operands.len,
    };
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
    self.genconstr_begin = try alloc.realloc(self.genconstr_begin, new + 1);
    self.genconstr_begin[new] = new_inds_len;
    self.genconstr_count = new;
    self.revision += 1;
}

/// Compute the number of usize entries stored in genconstr_indices for constraint at `idx`.
pub fn genConstrDataLen(self: Model, idx: usize) usize {
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
pub fn genConstrOffset(self: Model, idx: usize) usize {
    return self.genconstr_begin[idx];
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
    try addGenConstrWithExtra(self, .indicator, binvar, &.{}, pdata, name);
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
    try addGenConstrWithExtra(self, .pwl, yvar, &.{}, pdata, name);
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
    try addGenConstrWithExtra(self, .poly, resvar, &.{}, pdata, name);
}

/// Add an EXP general constraint: resvar = exp(xvar)
pub fn addGenConstrExp(self: *Model, resvar: usize, xvar: usize, name: ?[]const u8) ModelError!void {
    try addGenConstrWithExtra(self, .exp, resvar, &[_]usize{xvar}, &.{}, name);
}

/// Add an EXPONENT (base-a) general constraint: resvar = a^xvar
pub fn addGenConstrExpA(self: *Model, resvar: usize, xvar: usize, a: f64, name: ?[]const u8) ModelError!void {
    try addGenConstrWithExtra(self, .expa, resvar, &[_]usize{xvar}, &[_]usize{@as(usize, @bitCast(a))}, name);
}

/// Add a LOG general constraint: resvar = log(xvar)  (natural log)
pub fn addGenConstrLog(self: *Model, resvar: usize, xvar: usize, name: ?[]const u8) ModelError!void {
    try addGenConstrWithExtra(self, .log, resvar, &[_]usize{xvar}, &.{}, name);
}

/// Add a LOG (base-a) general constraint: resvar = log_a(xvar)
pub fn addGenConstrLogA(self: *Model, resvar: usize, xvar: usize, a: f64, name: ?[]const u8) ModelError!void {
    try addGenConstrWithExtra(self, .loga, resvar, &[_]usize{xvar}, &[_]usize{@as(usize, @bitCast(a))}, name);
}

/// Add a POWER general constraint: resvar = xvar^a
pub fn addGenConstrPow(self: *Model, resvar: usize, xvar: usize, a: f64, name: ?[]const u8) ModelError!void {
    try addGenConstrWithExtra(self, .pow, resvar, &[_]usize{xvar}, &[_]usize{@as(usize, @bitCast(a))}, name);
}

/// Add a SIN general constraint: resvar = sin(xvar)
pub fn addGenConstrSin(self: *Model, resvar: usize, xvar: usize, name: ?[]const u8) ModelError!void {
    try addGenConstrWithExtra(self, .sin, resvar, &[_]usize{xvar}, &.{}, name);
}

/// Add a COS general constraint: resvar = cos(xvar)
pub fn addGenConstrCos(self: *Model, resvar: usize, xvar: usize, name: ?[]const u8) ModelError!void {
    try addGenConstrWithExtra(self, .cos, resvar, &[_]usize{xvar}, &.{}, name);
}

/// Add a TAN general constraint: resvar = tan(xvar)
pub fn addGenConstrTan(self: *Model, resvar: usize, xvar: usize, name: ?[]const u8) ModelError!void {
    try addGenConstrWithExtra(self, .tan, resvar, &[_]usize{xvar}, &.{}, name);
}

/// Add a LOGISTIC general constraint: resvar = 1/(1 + exp(-xvar))
pub fn addGenConstrLogistic(self: *Model, resvar: usize, xvar: usize, name: ?[]const u8) ModelError!void {
    try addGenConstrWithExtra(self, .logistic, resvar, &[_]usize{xvar}, &.{}, name);
}

/// Add a NORM general constraint: resvar = sqrt(Σ varⁱ²)
pub fn addGenConstrNorm(self: *Model, resvar: usize, num_vars: usize, vars: []const usize, name: ?[]const u8) ModelError!void {
    try addGenConstrWithExtra(self, .norm, resvar, vars[0..num_vars], &.{}, name);
}

/// Add a NONLINEAR general constraint (generic container)
pub fn addGenConstrNL(self: *Model, resvar: usize, num_vars: usize, vars: []const usize, name: ?[]const u8) ModelError!void {
    try addGenConstrWithExtra(self, .nl, resvar, vars[0..num_vars], &.{}, name);
}

// ══════════════════════════════════════════════════════════════════════════
//  General constraint getters
// ══════════════════════════════════════════════════════════════════════════

/// Retrieve a MAX general constraint's data
pub fn getGenConstrMax(self: Model, idx: usize, resvar: *usize, num_vars: *usize, vars: []usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .max) return error.InvalidArgument;
    const off = genConstrOffset(self, idx);
    const nv = self.genconstr_nvars[idx];
    resvar.* = self.genconstr_resvar[idx];
    num_vars.* = nv;
    if (vars.len >= nv) @memcpy(vars[0..nv], self.genconstr_indices[off..][0..nv]);
}

/// Retrieve a MIN general constraint's data
pub fn getGenConstrMin(self: Model, idx: usize, resvar: *usize, num_vars: *usize, vars: []usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .min) return error.InvalidArgument;
    const off = genConstrOffset(self, idx);
    const nv = self.genconstr_nvars[idx];
    resvar.* = self.genconstr_resvar[idx];
    num_vars.* = nv;
    if (vars.len >= nv) @memcpy(vars[0..nv], self.genconstr_indices[off..][0..nv]);
}

/// Retrieve an ABS general constraint's data
pub fn getGenConstrAbs(self: Model, idx: usize, resvar: *usize, argvar: *usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .abs) return error.InvalidArgument;
    const off = genConstrOffset(self, idx);
    resvar.* = self.genconstr_resvar[idx];
    argvar.* = self.genconstr_indices[off];
}

/// Retrieve an AND general constraint's data
pub fn getGenConstrAnd(self: Model, idx: usize, resvar: *usize, var1: *usize, var2: *usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .and_) return error.InvalidArgument;
    const off = genConstrOffset(self, idx);
    resvar.* = self.genconstr_resvar[idx];
    var1.* = self.genconstr_indices[off];
    var2.* = self.genconstr_indices[off + 1];
}

/// Retrieve an OR general constraint's data
pub fn getGenConstrOr(self: Model, idx: usize, resvar: *usize, var1: *usize, var2: *usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .or_) return error.InvalidArgument;
    const off = genConstrOffset(self, idx);
    resvar.* = self.genconstr_resvar[idx];
    var1.* = self.genconstr_indices[off];
    var2.* = self.genconstr_indices[off + 1];
}

/// Retrieve an INDICATOR general constraint's data
pub fn getGenConstrIndicator(self: Model, idx: usize, binvar: *usize, binval: *i32, num_nz: *usize, ind: []usize, val: []f64, sense: *Sense, rhs: *f64) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .indicator) return error.InvalidArgument;
    const off = genConstrOffset(self, idx);
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
    const off = genConstrOffset(self, idx);
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
    const off = genConstrOffset(self, idx);
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
    const off = genConstrOffset(self, idx);
    resvar.* = self.genconstr_resvar[idx];
    xvar.* = self.genconstr_indices[off];
}

/// Retrieve an EXPA general constraint's data
pub fn getGenConstrExpA(self: Model, idx: usize, resvar: *usize, xvar: *usize, a: *f64) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .expa) return error.InvalidArgument;
    const off = genConstrOffset(self, idx);
    resvar.* = self.genconstr_resvar[idx];
    xvar.* = self.genconstr_indices[off];
    a.* = @as(f64, @bitCast(self.genconstr_indices[off + 1]));
}

/// Retrieve a LOG general constraint's data
pub fn getGenConstrLog(self: Model, idx: usize, resvar: *usize, xvar: *usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .log) return error.InvalidArgument;
    const off = genConstrOffset(self, idx);
    resvar.* = self.genconstr_resvar[idx];
    xvar.* = self.genconstr_indices[off];
}

/// Retrieve a LOGA general constraint's data
pub fn getGenConstrLogA(self: Model, idx: usize, resvar: *usize, xvar: *usize, a: *f64) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .loga) return error.InvalidArgument;
    const off = genConstrOffset(self, idx);
    resvar.* = self.genconstr_resvar[idx];
    xvar.* = self.genconstr_indices[off];
    a.* = @as(f64, @bitCast(self.genconstr_indices[off + 1]));
}

/// Retrieve a POW general constraint's data
pub fn getGenConstrPow(self: Model, idx: usize, resvar: *usize, xvar: *usize, a: *f64) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .pow) return error.InvalidArgument;
    const off = genConstrOffset(self, idx);
    resvar.* = self.genconstr_resvar[idx];
    xvar.* = self.genconstr_indices[off];
    a.* = @as(f64, @bitCast(self.genconstr_indices[off + 1]));
}

/// Retrieve a SIN general constraint's data
pub fn getGenConstrSin(self: Model, idx: usize, resvar: *usize, xvar: *usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .sin) return error.InvalidArgument;
    const off = genConstrOffset(self, idx);
    resvar.* = self.genconstr_resvar[idx];
    xvar.* = self.genconstr_indices[off];
}

/// Retrieve a COS general constraint's data
pub fn getGenConstrCos(self: Model, idx: usize, resvar: *usize, xvar: *usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .cos) return error.InvalidArgument;
    const off = genConstrOffset(self, idx);
    resvar.* = self.genconstr_resvar[idx];
    xvar.* = self.genconstr_indices[off];
}

/// Retrieve a TAN general constraint's data
pub fn getGenConstrTan(self: Model, idx: usize, resvar: *usize, xvar: *usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .tan) return error.InvalidArgument;
    const off = genConstrOffset(self, idx);
    resvar.* = self.genconstr_resvar[idx];
    xvar.* = self.genconstr_indices[off];
}

/// Retrieve a LOGISTIC general constraint's data
pub fn getGenConstrLogistic(self: Model, idx: usize, resvar: *usize, xvar: *usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .logistic) return error.InvalidArgument;
    const off = genConstrOffset(self, idx);
    resvar.* = self.genconstr_resvar[idx];
    xvar.* = self.genconstr_indices[off];
}

/// Retrieve a NORM general constraint's data
pub fn getGenConstrNorm(self: Model, idx: usize, resvar: *usize, num_vars: *usize, vars: []usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .norm) return error.InvalidArgument;
    const off = genConstrOffset(self, idx);
    const nv = self.genconstr_nvars[idx];
    resvar.* = self.genconstr_resvar[idx];
    num_vars.* = nv;
    if (vars.len >= nv) @memcpy(vars[0..nv], self.genconstr_indices[off..][0..nv]);
}

/// Retrieve an NL general constraint's data
pub fn getGenConstrNL(self: Model, idx: usize, resvar: *usize, num_vars: *usize, vars: []usize) ModelError!void {
    if (idx >= self.genconstr_count or self.genconstr_types[idx] != .nl) return error.InvalidArgument;
    const off = genConstrOffset(self, idx);
    const nv = self.genconstr_nvars[idx];
    resvar.* = self.genconstr_resvar[idx];
    num_vars.* = nv;
    if (vars.len >= nv) @memcpy(vars[0..nv], self.genconstr_indices[off..][0..nv]);
}
