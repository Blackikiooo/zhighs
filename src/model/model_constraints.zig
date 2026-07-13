//! Quadratic objective, quadratic constraint, and SOS methods for `Model`.
//!
//! ## Responsibility
//!
//! Owns storage-facing add/delete/query operations for quadratic objective
//! terms, quadratic constraints, and special ordered sets.  Linear constraints
//! belong to `model_linear.zig`; nonlinear/general constraints belong to
//! `model_genconstr.zig`.

const std = @import("std");
const types = @import("types.zig");
const Model = @import("model.zig").Model;

const ModelError = types.ModelError;
const Sense = types.Sense;
const SosType = types.SosType;

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

pub fn addQConstr(self: *Model, qnz: usize, qrow: []const i32, qcol: []const i32, qval: []const f64, lnz: usize, lind: []const usize, lval: []const f64, sense: Sense, rhs: f64, name: ?[]const u8) ModelError!void {
    const alloc = self.allocator;
    if (qrow.len < qnz or qcol.len < qnz or qval.len < qnz) return error.InvalidArgument;
    if (lind.len < lnz or lval.len < lnz) return error.InvalidArgument;

    const old = self.qconstr_count;
    const new = old + 1;
    const q_start = self.qconstr_qrow.len;
    const q_end = q_start + qnz;
    const l_start = self.qconstr_lind.len;
    const l_end = l_start + lnz;

    self.qconstr_qrow = try alloc.realloc(self.qconstr_qrow, q_end);
    self.qconstr_qcol = try alloc.realloc(self.qconstr_qcol, q_end);
    self.qconstr_qval = try alloc.realloc(self.qconstr_qval, q_end);
    @memcpy(self.qconstr_qrow[q_start..q_end], qrow[0..qnz]);
    @memcpy(self.qconstr_qcol[q_start..q_end], qcol[0..qnz]);
    @memcpy(self.qconstr_qval[q_start..q_end], qval[0..qnz]);

    self.qconstr_lind = try alloc.realloc(self.qconstr_lind, l_end);
    self.qconstr_lval = try alloc.realloc(self.qconstr_lval, l_end);
    @memcpy(self.qconstr_lind[l_start..l_end], lind[0..lnz]);
    @memcpy(self.qconstr_lval[l_start..l_end], lval[0..lnz]);

    self.qconstr_sense = try alloc.realloc(self.qconstr_sense, new);
    self.qconstr_sense[old] = sense;
    self.qconstr_rhs = try alloc.realloc(self.qconstr_rhs, new);
    self.qconstr_rhs[old] = rhs;
    self.qconstr_names = try alloc.realloc(self.qconstr_names, new);
    self.qconstr_names[old] = if (name) |n| try alloc.dupe(u8, n) else null;
    self.qconstr_count = new;
    self.revision += 1;
}

pub fn delQConstrs(self: *Model, indices: []const usize) ModelError!void {
    _ = self;
    _ = indices;
    return error.FeatureNotAvailable;
}

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

pub fn addSOS(self: *Model, sostype: SosType, num_members: usize, indices: []const usize, weights: ?[]const f64, name: ?[]const u8) ModelError!void {
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

    const old_weights_len = self.sos_weights.len;
    const new_weights_len = old_weights_len + num_members;
    self.sos_weights = try alloc.realloc(self.sos_weights, new_weights_len);
    if (weights) |w| {
        @memcpy(self.sos_weights[old_weights_len..new_weights_len], w[0..num_members]);
    } else {
        for (0..num_members) |i| self.sos_weights[old_weights_len + i] = @floatFromInt(i + 1);
    }

    self.sos_names = try alloc.realloc(self.sos_names, new);
    self.sos_names[old] = if (name) |n| try alloc.dupe(u8, n) else null;
    self.sos_count = new;
    self.revision += 1;
}

pub fn delSOS(self: *Model, indices: []const usize) ModelError!void {
    _ = self;
    _ = indices;
    return error.FeatureNotAvailable;
}

pub fn getSOS(self: Model, idx: usize, sostype: *SosType, num_members: *usize, indices: []usize, weights: []f64) ModelError!void {
    if (idx >= self.sos_count) return error.IndexOutOfRange;
    sostype.* = self.sos_types[idx];
    const begin = self.sos_begin[idx];
    const end = self.sos_begin[idx + 1];
    num_members.* = end - begin;
    const count = num_members.*;
    if (indices.len >= count) @memcpy(indices[0..count], self.sos_indices[begin..end]);
    if (weights.len >= count) @memcpy(weights[0..count], self.sos_weights[begin..end]);
}

pub fn getQConstrByName(self: Model, name: []const u8) ModelError!usize {
    for (self.qconstr_names, 0..) |maybe_name, i| {
        if (maybe_name) |candidate| if (std.mem.eql(u8, candidate, name)) return i;
    }
    return error.NotInModel;
}
