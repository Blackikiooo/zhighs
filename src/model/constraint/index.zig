//! Typed index wrapper for a linear constraint.

const std = @import("std");
const Model = @import("../model.zig").Model;
const types = @import("../types.zig");
const Attr = @import("../attrs.zig").Attr;

const Sense = types.Sense;
const BasisStatus = types.BasisStatus;
const ModelError = types.ModelError;

/// A linear constraint, identified by index within a model.
pub const Constr = struct {
    model: *Model,
    index: usize,

    const Self = @This();

    // ── Generic attribute access ──────────────────────────────────────────

    pub fn get(self: Self, attr: Attr) ModelError!f64 {
        return self.model.getDblAttrElement(attr, self.index);
    }

    pub fn set(self: *Self, attr: Attr, value: f64) ModelError!void {
        try self.model.setDblAttrElement(attr, self.index, value);
    }

    pub fn getInt(self: Self, attr: Attr) ModelError!i64 {
        return self.model.getIntAttrElement(attr, self.index);
    }

    pub fn setInt(self: *Self, attr: Attr, value: i64) ModelError!void {
        try self.model.setIntAttrElement(attr, self.index, value);
    }

    pub fn getChar(self: Self, attr: Attr) ModelError!u8 {
        return self.model.getCharAttrElement(attr, self.index);
    }

    pub fn setChar(self: *Self, attr: Attr, value: u8) ModelError!void {
        try self.model.setCharAttrElement(attr, self.index, value);
    }

    pub fn getStr(self: Self, attr: Attr) ModelError![]const u8 {
        return self.model.getStrAttrElement(attr, self.index);
    }

    pub fn setStr(self: *Self, attr: Attr, value: []const u8) ModelError!void {
        try self.model.setStrAttrElement(attr, self.index, value);
    }

    // ── Convenience accessors ─────────────────────────────────────────────

    pub fn getSense(self: Self) ModelError!Sense {
        const code = try self.model.getCharAttrElement(.sense, self.index);
        return Sense.fromCode(code);
    }
    pub fn setSense(self: *Self, sense: Sense) ModelError!void {
        try self.model.setCharAttrElement(.sense, self.index, @intFromEnum(sense));
    }

    pub fn getRHS(self: Self) ModelError!f64 {
        return self.model.getDblAttrElement(.rhs, self.index);
    }
    pub fn setRHS(self: *Self, rhs: f64) ModelError!void {
        try self.model.setDblAttrElement(.rhs, self.index, rhs);
    }

    pub fn getPi(self: Self) ModelError!f64 {
        return self.model.getDblAttrElement(.pi, self.index);
    }

    pub fn getSlack(self: Self) ModelError!f64 {
        return self.model.getDblAttrElement(.slack, self.index);
    }

    pub fn getConstrName(self: Self) ModelError![]const u8 {
        return self.model.getStrAttrElement(.constr_name, self.index);
    }
    pub fn setConstrName(self: *Self, name: []const u8) ModelError!void {
        try self.model.setStrAttrElement(.constr_name, self.index, name);
    }

    pub fn getBasis(self: Self) ModelError!BasisStatus {
        const v = try self.model.getIntAttrElement(.c_basis, self.index);
        return @enumFromInt(@as(i8, @intCast(v)));
    }
    pub fn setBasis(self: *Self, status: BasisStatus) ModelError!void {
        try self.model.setIntAttrElement(.c_basis, self.index, @intFromEnum(status));
    }

    // ── Comparison ────────────────────────────────────────────────────────

    pub fn sameAs(self: Self, other: Self) bool {
        return self.model == other.model and self.index == other.index;
    }
};
