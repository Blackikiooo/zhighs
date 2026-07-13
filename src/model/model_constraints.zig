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
const QuadExpr = @import("expr/quad_expr.zig").QuadExpr;

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

/// Add quadratic objective terms from an expression, resolving stable VarIds
/// once before appending the dense lower-triangle triples.
pub fn addQPtermsExpr(self: *Model, expr: QuadExpr) ModelError!void {
    const resolved = try expr.resolveQTerms(self.allocator, self.*);
    defer self.allocator.free(resolved);
    const qrow = try self.allocator.alloc(i32, resolved.len);
    defer self.allocator.free(qrow);
    const qcol = try self.allocator.alloc(i32, resolved.len);
    defer self.allocator.free(qcol);
    const qval = try self.allocator.alloc(f64, resolved.len);
    defer self.allocator.free(qval);
    for (resolved, 0..) |term, i| {
        qrow[i] = std.math.cast(i32, term.row) orelse return error.InvalidArgument;
        qcol[i] = std.math.cast(i32, term.col) orelse return error.InvalidArgument;
        qval[i] = term.coeff;
    }
    return self.addQPterms(qrow, qcol, qval);
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

    if (old == 0) {
        self.qconstr_qbegin = try alloc.realloc(self.qconstr_qbegin, 1);
        self.qconstr_qbegin[0] = 0;
        self.qconstr_lbegin = try alloc.realloc(self.qconstr_lbegin, 1);
        self.qconstr_lbegin[0] = 0;
    }

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

    self.qconstr_qbegin = try alloc.realloc(self.qconstr_qbegin, new + 1);
    self.qconstr_qbegin[new] = q_end;
    self.qconstr_lbegin = try alloc.realloc(self.qconstr_lbegin, new + 1);
    self.qconstr_lbegin[new] = l_end;

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
    const alloc = self.allocator;
    if (indices.len == 0) return;
    var deleted = try alloc.alloc(bool, self.qconstr_count);
    defer alloc.free(deleted);
    @memset(deleted, false);
    for (indices) |idx| {
        if (idx >= self.qconstr_count) return error.IndexOutOfRange;
        deleted[idx] = true;
    }
    while (self.qconstr_handles.liveLen() < self.qconstr_count) {
        const id = self.qconstr_handles.allocate(alloc) catch return error.OutOfMemory;
        self.qconstr_handles.bindDenseWithAllocator(alloc, id, @intCast(self.qconstr_handles.liveLen())) catch return error.OutOfMemory;
    }

    var qbegin: std.ArrayListUnmanaged(usize) = .empty;
    var qrow: std.ArrayListUnmanaged(i32) = .empty;
    var qcol: std.ArrayListUnmanaged(i32) = .empty;
    var qval: std.ArrayListUnmanaged(f64) = .empty;
    var lbegin: std.ArrayListUnmanaged(usize) = .empty;
    var lind: std.ArrayListUnmanaged(usize) = .empty;
    var lval: std.ArrayListUnmanaged(f64) = .empty;
    var sense: std.ArrayListUnmanaged(Sense) = .empty;
    var rhs: std.ArrayListUnmanaged(f64) = .empty;
    var names: std.ArrayListUnmanaged(?[]const u8) = .empty;
    qbegin.append(alloc, 0) catch return error.OutOfMemory;
    lbegin.append(alloc, 0) catch return error.OutOfMemory;
    for (0..self.qconstr_count) |i| {
        if (deleted[i]) {
            if (self.qconstr_names[i]) |name| alloc.free(name);
            continue;
        }
        const qb = self.qconstr_qbegin[i];
        const qe = self.qconstr_qbegin[i + 1];
        qrow.appendSlice(alloc, self.qconstr_qrow[qb..qe]) catch return error.OutOfMemory;
        qcol.appendSlice(alloc, self.qconstr_qcol[qb..qe]) catch return error.OutOfMemory;
        qval.appendSlice(alloc, self.qconstr_qval[qb..qe]) catch return error.OutOfMemory;
        qbegin.append(alloc, qrow.items.len) catch return error.OutOfMemory;
        const lb = self.qconstr_lbegin[i];
        const le = self.qconstr_lbegin[i + 1];
        lind.appendSlice(alloc, self.qconstr_lind[lb..le]) catch return error.OutOfMemory;
        lval.appendSlice(alloc, self.qconstr_lval[lb..le]) catch return error.OutOfMemory;
        lbegin.append(alloc, lind.items.len) catch return error.OutOfMemory;
        sense.append(alloc, self.qconstr_sense[i]) catch return error.OutOfMemory;
        rhs.append(alloc, self.qconstr_rhs[i]) catch return error.OutOfMemory;
        names.append(alloc, self.qconstr_names[i]) catch return error.OutOfMemory;
    }
    alloc.free(self.qconstr_qrow);
    alloc.free(self.qconstr_qcol);
    alloc.free(self.qconstr_qval);
    alloc.free(self.qconstr_qbegin);
    alloc.free(self.qconstr_lind);
    alloc.free(self.qconstr_lval);
    alloc.free(self.qconstr_lbegin);
    alloc.free(self.qconstr_sense);
    alloc.free(self.qconstr_rhs);
    alloc.free(self.qconstr_names);
    self.qconstr_qrow = try qrow.toOwnedSlice(alloc);
    self.qconstr_qcol = try qcol.toOwnedSlice(alloc);
    self.qconstr_qval = try qval.toOwnedSlice(alloc);
    self.qconstr_qbegin = try qbegin.toOwnedSlice(alloc);
    self.qconstr_lind = try lind.toOwnedSlice(alloc);
    self.qconstr_lval = try lval.toOwnedSlice(alloc);
    self.qconstr_lbegin = try lbegin.toOwnedSlice(alloc);
    self.qconstr_sense = try sense.toOwnedSlice(alloc);
    self.qconstr_rhs = try rhs.toOwnedSlice(alloc);
    self.qconstr_names = try names.toOwnedSlice(alloc);
    self.qconstr_count = self.qconstr_sense.len;
    const remap = self.qconstr_handles.compact(alloc, deleted) catch return error.OutOfMemory;
    alloc.free(remap);
    self.revision += 1;
}

pub fn getQConstr(self: Model, idx: usize, qnz: *usize, qrow: []i32, qcol: []i32, qval: []f64, lnz: *usize, lind: []usize, lval: []f64) ModelError!void {
    if (idx >= self.qconstr_count) return error.IndexOutOfRange;
    const qb = self.qconstr_qbegin[idx];
    const qe = self.qconstr_qbegin[idx + 1];
    const lb = self.qconstr_lbegin[idx];
    const le = self.qconstr_lbegin[idx + 1];
    qnz.* = qe - qb;
    lnz.* = le - lb;
    if (qrow.len >= qnz.*) @memcpy(qrow[0..qnz.*], self.qconstr_qrow[qb..qe]);
    if (qcol.len >= qnz.*) @memcpy(qcol[0..qnz.*], self.qconstr_qcol[qb..qe]);
    if (qval.len >= qnz.*) @memcpy(qval[0..qnz.*], self.qconstr_qval[qb..qe]);
    if (lind.len >= lnz.*) @memcpy(lind[0..lnz.*], self.qconstr_lind[lb..le]);
    if (lval.len >= lnz.*) @memcpy(lval[0..lnz.*], self.qconstr_lval[lb..le]);
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
    const alloc = self.allocator;
    if (indices.len == 0) return;
    var deleted = try alloc.alloc(bool, self.sos_count);
    defer alloc.free(deleted);
    @memset(deleted, false);
    for (indices) |idx| {
        if (idx >= self.sos_count) return error.IndexOutOfRange;
        deleted[idx] = true;
    }
    while (self.sos_handles.liveLen() < self.sos_count) {
        const id = self.sos_handles.allocate(alloc) catch return error.OutOfMemory;
        self.sos_handles.bindDenseWithAllocator(alloc, id, @intCast(self.sos_handles.liveLen())) catch return error.OutOfMemory;
    }

    var sos_types: std.ArrayListUnmanaged(SosType) = .empty;
    var begin: std.ArrayListUnmanaged(usize) = .empty;
    var members: std.ArrayListUnmanaged(usize) = .empty;
    var weights: std.ArrayListUnmanaged(f64) = .empty;
    var names: std.ArrayListUnmanaged(?[]const u8) = .empty;
    begin.append(alloc, 0) catch return error.OutOfMemory;
    for (0..self.sos_count) |i| {
        if (deleted[i]) {
            if (self.sos_names[i]) |name| alloc.free(name);
            continue;
        }
        const start = self.sos_begin[i];
        const end = self.sos_begin[i + 1];
        sos_types.append(alloc, self.sos_types[i]) catch return error.OutOfMemory;
        members.appendSlice(alloc, self.sos_indices[start..end]) catch return error.OutOfMemory;
        weights.appendSlice(alloc, self.sos_weights[start..end]) catch return error.OutOfMemory;
        begin.append(alloc, members.items.len) catch return error.OutOfMemory;
        names.append(alloc, self.sos_names[i]) catch return error.OutOfMemory;
    }
    alloc.free(self.sos_types);
    alloc.free(self.sos_begin);
    alloc.free(self.sos_indices);
    alloc.free(self.sos_weights);
    alloc.free(self.sos_names);
    self.sos_types = try sos_types.toOwnedSlice(alloc);
    self.sos_begin = try begin.toOwnedSlice(alloc);
    self.sos_indices = try members.toOwnedSlice(alloc);
    self.sos_weights = try weights.toOwnedSlice(alloc);
    self.sos_names = try names.toOwnedSlice(alloc);
    self.sos_count = self.sos_types.len;
    const remap = self.sos_handles.compact(alloc, deleted) catch return error.OutOfMemory;
    alloc.free(remap);
    self.revision += 1;
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
