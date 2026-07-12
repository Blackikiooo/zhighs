//! Attribute get/set methods for `Model`.
//!
//! Problem and solution data are read/written through uniform
//! `get*Attr` / `set*Attr` methods keyed by an `Attr` enum value.

const std = @import("std");
const types = @import("types.zig");
const Model = @import("model.zig").Model;
const Attr = @import("attrs.zig").Attr;

const ModelError = types.ModelError;
const VarType = types.VarType;
const Sense = types.Sense;
const BasisStatus = types.BasisStatus;
const ObjectiveSense = types.ObjectiveSense;
const Status = types.Status;

// ══════════════════════════════════════════════════════════════════════════
//  Scalar attributes
// ══════════════════════════════════════════════════════════════════════════

// ── Integer scalar attributes ───────────────────────────────────────────

pub fn getIntAttr(self: Model, attr: Attr) ModelError!i64 {
    return switch (attr) {
        .num_vars => @intCast(self.num_vars),
        .num_constrs => @intCast(self.num_constrs),
        .num_nz => @intCast(self.numNz()),
        .status => @intFromEnum(self.status),
        .model_sense => self.sense.toModelSenseValue(),
        .is_mip => if (self.isMip()) 1 else 0,
        .iter_count => self.iter_count,
        .node_count => self.node_count,
        .bar_iter_count => self.bar_iter_count,
        .num_qnz => @intCast(self.q_nz),
        .num_sos => @intCast(self.sos_count),
        .num_gen_constrs => @intCast(self.genconstr_count),
        .num_q_constrs => @intCast(self.qconstr_count),
        .num_bin_vars => @intCast(countVarType(self, .binary)),
        .num_int_vars => @intCast(countVarType(self, .integer)),
        .sol_count => if (self.status == .optimal) 1 else 0,
        else => error.InvalidAttribute,
    };
}

pub fn setIntAttr(self: *Model, attr: Attr, value: i64) ModelError!void {
    switch (attr) {
        .model_sense => {
            self.sense = try ObjectiveSense.fromModelSenseValue(@intCast(value));
        },
        else => return error.InvalidAttribute,
    }
}

// ── Double scalar attributes ────────────────────────────────────────────

pub fn getDblAttr(self: Model, attr: Attr) ModelError!f64 {
    return switch (attr) {
        .obj_val => self.obj_val,
        .obj_bound => self.obj_bound,
        .obj_con => self.obj_con,
        else => error.InvalidAttribute,
    };
}

pub fn setDblAttr(self: *Model, attr: Attr, value: f64) ModelError!void {
    switch (attr) {
        .obj_con => {
            self.obj_con = value;
        },
        else => return error.InvalidAttribute,
    }
}

// ── String scalar attributes ────────────────────────────────────────────

pub fn getStrAttr(self: Model, attr: Attr) ModelError![]const u8 {
    return switch (attr) {
        .model_name => self.name,
        .status_label => self.status.label(),
        else => error.InvalidAttribute,
    };
}

pub fn setStrAttr(self: *Model, attr: Attr, value: []const u8) ModelError!void {
    switch (attr) {
        .model_name => {
            self.allocator.free(self.name);
            self.name = self.allocator.dupe(u8, value) catch return error.OutOfMemory;
        },
        else => return error.InvalidAttribute,
    }
}

// ══════════════════════════════════════════════════════════════════════════
//  Element attributes
// ══════════════════════════════════════════════════════════════════════════

// ── Integer element attributes ──────────────────────────────────────────

pub fn getIntAttrElement(self: Model, attr: Attr, element: usize) ModelError!i64 {
    switch (attr) {
        .v_basis => {
            if (element >= self.vbasis.len) return error.IndexOutOfRange;
            return @intFromEnum(self.vbasis[element]);
        },
        .c_basis => {
            if (element >= self.cbasis.len) return error.IndexOutOfRange;
            return @intFromEnum(self.cbasis[element]);
        },
        else => return error.InvalidAttribute,
    }
}

pub fn setIntAttrElement(self: *Model, attr: Attr, element: usize, value: i64) ModelError!void {
    switch (attr) {
        .v_basis => {
            if (element >= self.vbasis.len) return error.IndexOutOfRange;
            self.vbasis[element] = @enumFromInt(@as(i8, @intCast(value)));
        },
        .c_basis => {
            if (element >= self.cbasis.len) return error.IndexOutOfRange;
            self.cbasis[element] = @enumFromInt(@as(i8, @intCast(value)));
        },
        else => return error.InvalidAttribute,
    }
}

// ── Double element attributes ───────────────────────────────────────────

pub fn getDblAttrElement(self: Model, attr: Attr, element: usize) ModelError!f64 {
    switch (attr) {
        .lb => {
            if (element >= self.var_lb.len) return error.IndexOutOfRange;
            return self.var_lb[element];
        },
        .ub => {
            if (element >= self.var_ub.len) return error.IndexOutOfRange;
            return self.var_ub[element];
        },
        .obj => {
            if (element >= self.var_obj.len) return error.IndexOutOfRange;
            return self.var_obj[element];
        },
        .x => {
            if (element >= self.solution.len) return error.IndexOutOfRange;
            return self.solution[element];
        },
        .rc => {
            if (element >= self.reduced_cost.len) return error.IndexOutOfRange;
            return self.reduced_cost[element];
        },
        .rhs => {
            if (element >= self.constr_rhs.len) return error.IndexOutOfRange;
            return self.constr_rhs[element];
        },
        .pi => {
            if (element >= self.pi.len) return error.IndexOutOfRange;
            return self.pi[element];
        },
        .slack => {
            if (element >= self.slack.len) return error.IndexOutOfRange;
            return self.slack[element];
        },
        .start => {
            if (element >= self.mip_start.len) return error.IndexOutOfRange;
            return self.mip_start[element];
        },
        .p_start => {
            if (element >= self.p_start.len) return error.IndexOutOfRange;
            return self.p_start[element];
        },
        .d_start => {
            if (element >= self.d_start.len) return error.IndexOutOfRange;
            return self.d_start[element];
        },
        else => return error.InvalidAttribute,
    }
}

pub fn setDblAttrElement(self: *Model, attr: Attr, element: usize, value: f64) ModelError!void {
    switch (attr) {
        .lb => {
            if (element >= self.var_lb.len) return error.IndexOutOfRange;
            self.var_lb[element] = value;
        },
        .ub => {
            if (element >= self.var_ub.len) return error.IndexOutOfRange;
            self.var_ub[element] = value;
        },
        .obj => {
            if (element >= self.var_obj.len) return error.IndexOutOfRange;
            self.var_obj[element] = value;
        },
        .rhs => {
            if (element >= self.constr_rhs.len) return error.IndexOutOfRange;
            self.constr_rhs[element] = value;
        },
        .start => {
            if (element >= self.mip_start.len) return error.IndexOutOfRange;
            self.mip_start[element] = value;
        },
        .p_start => {
            if (element >= self.p_start.len) return error.IndexOutOfRange;
            self.p_start[element] = value;
        },
        .d_start => {
            if (element >= self.d_start.len) return error.IndexOutOfRange;
            self.d_start[element] = value;
        },
        else => return error.InvalidAttribute,
    }
}

// ── Char element attributes ─────────────────────────────────────────────

pub fn getCharAttrElement(self: Model, attr: Attr, element: usize) ModelError!u8 {
    switch (attr) {
        .v_type => {
            if (element >= self.var_type.len) return error.IndexOutOfRange;
            return @intFromEnum(self.var_type[element]);
        },
        .sense => {
            if (element >= self.constr_sense.len) return error.IndexOutOfRange;
            return @intFromEnum(self.constr_sense[element]);
        },
        else => return error.InvalidAttribute,
    }
}

pub fn setCharAttrElement(self: *Model, attr: Attr, element: usize, value: u8) ModelError!void {
    switch (attr) {
        .v_type => {
            if (element >= self.var_type.len) return error.IndexOutOfRange;
            self.var_type[element] = try VarType.fromCode(value);
        },
        .sense => {
            if (element >= self.constr_sense.len) return error.IndexOutOfRange;
            self.constr_sense[element] = try Sense.fromCode(value);
        },
        else => return error.InvalidAttribute,
    }
}

// ── String element attributes ───────────────────────────────────────────

pub fn getStrAttrElement(self: Model, attr: Attr, element: usize) ModelError![]const u8 {
    switch (attr) {
        .var_name => {
            if (element >= self.var_names.len) return error.IndexOutOfRange;
            return self.var_names[element] orelse "";
        },
        .constr_name => {
            if (element >= self.constr_names.len) return error.IndexOutOfRange;
            return self.constr_names[element] orelse "";
        },
        else => return error.InvalidAttribute,
    }
}

pub fn setStrAttrElement(self: *Model, attr: Attr, element: usize, value: []const u8) ModelError!void {
    const alloc = self.allocator;
    switch (attr) {
        .var_name => {
            if (element >= self.var_names.len) return error.IndexOutOfRange;
            if (self.var_names[element]) |old| alloc.free(old);
            self.var_names[element] = try alloc.dupe(u8, value);
        },
        .constr_name => {
            if (element >= self.constr_names.len) return error.IndexOutOfRange;
            if (self.constr_names[element]) |old| alloc.free(old);
            self.constr_names[element] = try alloc.dupe(u8, value);
        },
        else => return error.InvalidAttribute,
    }
}

// ══════════════════════════════════════════════════════════════════════════
//  Array attribute setters (contiguous ranges)
// ══════════════════════════════════════════════════════════════════════════

/// Set a range of double attribute values starting at `start`
pub fn setDblAttrArray(self: *Model, attr: Attr, start: usize, values: []const f64) ModelError!void {
    for (values, 0..) |v, i| {
        try self.setDblAttrElement(attr, start + i, v);
    }
}

/// Set a range of integer attribute values starting at `start`
pub fn setIntAttrArray(self: *Model, attr: Attr, start: usize, values: []const i64) ModelError!void {
    for (values, 0..) |v, i| {
        try self.setIntAttrElement(attr, start + i, v);
    }
}

/// Set a range of string attribute values starting at `start`
pub fn setStrAttrArray(self: *Model, attr: Attr, start: usize, values: []const []const u8) ModelError!void {
    for (values, 0..) |v, i| {
        try self.setStrAttrElement(attr, start + i, v);
    }
}

/// Set a range of char attribute values starting at `start`
pub fn setCharAttrArray(self: *Model, attr: Attr, start: usize, values: []const u8) ModelError!void {
    for (values, 0..) |v, i| {
        try self.setCharAttrElement(attr, start + i, v);
    }
}

/// Retrieve a range of string attribute values
pub fn getStrAttrArray(self: Model, attr: Attr, start: usize, len: usize, values: [][]const u8) ModelError!void {
    if (values.len < len) return error.InvalidArgument;
    for (0..len) |i| {
        values[i] = try self.getStrAttrElement(attr, start + i);
    }
}

/// Retrieve a range of double attribute values
pub fn getDblAttrArray(self: Model, attr: Attr, start: usize, len: usize, values: []f64) ModelError!void {
    if (values.len < len) return error.InvalidArgument;
    for (0..len) |i| {
        values[i] = try self.getDblAttrElement(attr, start + i);
    }
}

/// Retrieve a range of integer attribute values
pub fn getIntAttrArray(self: Model, attr: Attr, start: usize, len: usize, values: []i64) ModelError!void {
    if (values.len < len) return error.InvalidArgument;
    for (0..len) |i| {
        values[i] = try self.getIntAttrElement(attr, start + i);
    }
}

/// Retrieve a range of char attribute values
pub fn getCharAttrArray(self: Model, attr: Attr, start: usize, len: usize, values: []u8) ModelError!void {
    if (values.len < len) return error.InvalidArgument;
    for (0..len) |i| {
        values[i] = try self.getCharAttrElement(attr, start + i);
    }
}

// ══════════════════════════════════════════════════════════════════════════
//  List attribute operations (by explicit index list)
// ══════════════════════════════════════════════════════════════════════════

pub fn setDblAttrList(self: *Model, attr: Attr, indices: []const usize, values: []const f64) ModelError!void {
    if (indices.len != values.len) return error.InvalidArgument;
    for (indices, values) |idx, v| {
        try self.setDblAttrElement(attr, idx, v);
    }
}

pub fn setIntAttrList(self: *Model, attr: Attr, indices: []const usize, values: []const i64) ModelError!void {
    if (indices.len != values.len) return error.InvalidArgument;
    for (indices, values) |idx, v| {
        try self.setIntAttrElement(attr, idx, v);
    }
}

pub fn setStrAttrList(self: *Model, attr: Attr, indices: []const usize, values: []const []const u8) ModelError!void {
    if (indices.len != values.len) return error.InvalidArgument;
    for (indices, values) |idx, v| {
        try self.setStrAttrElement(attr, idx, v);
    }
}

pub fn setCharAttrList(self: *Model, attr: Attr, indices: []const usize, values: []const u8) ModelError!void {
    if (indices.len != values.len) return error.InvalidArgument;
    for (indices, values) |idx, v| {
        try self.setCharAttrElement(attr, idx, v);
    }
}

pub fn getDblAttrList(self: Model, attr: Attr, indices: []const usize, values: []f64) ModelError!void {
    if (indices.len != values.len) return error.InvalidArgument;
    for (indices, values) |idx, *v| {
        v.* = try self.getDblAttrElement(attr, idx);
    }
}

pub fn getIntAttrList(self: Model, attr: Attr, indices: []const usize, values: []i64) ModelError!void {
    if (indices.len != values.len) return error.InvalidArgument;
    for (indices, values) |idx, *v| {
        v.* = try self.getIntAttrElement(attr, idx);
    }
}

pub fn getStrAttrList(self: Model, attr: Attr, indices: []const usize, values: [][]const u8) ModelError!void {
    if (indices.len != values.len) return error.InvalidArgument;
    for (indices, values) |idx, *v| {
        v.* = try self.getStrAttrElement(attr, idx);
    }
}

pub fn getCharAttrList(self: Model, attr: Attr, indices: []const usize, values: []u8) ModelError!void {
    if (indices.len != values.len) return error.InvalidArgument;
    for (indices, values) |idx, *v| {
        v.* = try self.getCharAttrElement(attr, idx);
    }
}

// ── Helpers ──────────────────────────────────────────────────────────────

/// Count how many variables have the given type.
fn countVarType(self: Model, vt: VarType) usize {
    var count: usize = 0;
    for (self.var_type) |t| {
        if (t == vt) count += 1;
    }
    return count;
}
