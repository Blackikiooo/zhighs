//! Typed index wrapper for a linear constraint.

const std = @import("std");
const Model = @import("../model.zig").Model;
const types = @import("../types.zig");
const Attr = @import("../attrs.zig").Attr;

const Sense = types.Sense;
const BasisStatus = types.BasisStatus;
const ModelError = types.ModelError;
const ConstrId = @import("../entity_handle.zig").ConstrId;

/// A linear constraint, identified by index within a model.
pub const Constr = struct {
    model: *Model,
    index: usize,
    id: ?ConstrId = null,

    const Self = @This();

    pub fn at(model: *Model, index: usize) ModelError!Self {
        return .{ .model = model, .index = index, .id = try model.constrIdAt(index) };
    }

    fn denseIndex(self: Self) ModelError!usize {
        if (self.id) |id| return self.model.resolveConstrId(id);
        return self.index;
    }

    // ── Generic attribute access ──────────────────────────────────────────

    pub fn get(self: Self, attr: Attr) ModelError!f64 {
        return self.model.getDblAttrElement(attr, (try self.denseIndex()));
    }

    pub fn set(self: *Self, attr: Attr, value: f64) ModelError!void {
        try self.model.setDblAttrElement(attr, (try self.denseIndex()), value);
    }

    pub fn getInt(self: Self, attr: Attr) ModelError!i64 {
        return self.model.getIntAttrElement(attr, (try self.denseIndex()));
    }

    pub fn setInt(self: *Self, attr: Attr, value: i64) ModelError!void {
        try self.model.setIntAttrElement(attr, (try self.denseIndex()), value);
    }

    pub fn getChar(self: Self, attr: Attr) ModelError!u8 {
        return self.model.getCharAttrElement(attr, (try self.denseIndex()));
    }

    pub fn setChar(self: *Self, attr: Attr, value: u8) ModelError!void {
        try self.model.setCharAttrElement(attr, (try self.denseIndex()), value);
    }

    pub fn getStr(self: Self, attr: Attr) ModelError![]const u8 {
        return self.model.getStrAttrElement(attr, (try self.denseIndex()));
    }

    pub fn setStr(self: *Self, attr: Attr, value: []const u8) ModelError!void {
        try self.model.setStrAttrElement(attr, (try self.denseIndex()), value);
    }

    // ── Convenience accessors ─────────────────────────────────────────────

    pub fn getSense(self: Self) ModelError!Sense {
        const code = try self.model.getCharAttrElement(.sense, (try self.denseIndex()));
        return Sense.fromCode(code);
    }
    pub fn setSense(self: *Self, sense: Sense) ModelError!void {
        try self.model.setCharAttrElement(.sense, (try self.denseIndex()), @intFromEnum(sense));
    }

    pub fn getRHS(self: Self) ModelError!f64 {
        return self.model.getDblAttrElement(.rhs, (try self.denseIndex()));
    }
    pub fn setRHS(self: *Self, rhs: f64) ModelError!void {
        try self.model.setDblAttrElement(.rhs, (try self.denseIndex()), rhs);
    }

    pub fn getPi(self: Self) ModelError!f64 {
        return self.model.getDblAttrElement(.pi, (try self.denseIndex()));
    }

    pub fn getSlack(self: Self) ModelError!f64 {
        return self.model.getDblAttrElement(.slack, (try self.denseIndex()));
    }

    pub fn getConstrName(self: Self) ModelError![]const u8 {
        return self.model.getStrAttrElement(.constr_name, (try self.denseIndex()));
    }
    pub fn setConstrName(self: *Self, name: []const u8) ModelError!void {
        try self.model.setStrAttrElement(.constr_name, (try self.denseIndex()), name);
    }

    pub fn getBasis(self: Self) ModelError!BasisStatus {
        const v = try self.model.getIntAttrElement(.c_basis, (try self.denseIndex()));
        return @enumFromInt(@as(i8, @intCast(v)));
    }
    pub fn setBasis(self: *Self, status: BasisStatus) ModelError!void {
        try self.model.setIntAttrElement(.c_basis, (try self.denseIndex()), @intFromEnum(status));
    }

    // ── Comparison ────────────────────────────────────────────────────────

    pub fn sameAs(self: Self, other: Self) bool {
        return self.model == other.model and self.index == other.index;
    }
};
