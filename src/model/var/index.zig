//! Typed index wrapper for a decision variable.
//!
//! A `Var` pairs a `*Model` with a `usize` index, providing strongly-typed
//! attribute accessors without exposing the string-keyed attribute system.

const std = @import("std");
const Model = @import("../model.zig").Model;
const types = @import("../types.zig");
const Attr = @import("../attrs.zig").Attr;

const VarType = types.VarType;
const BasisStatus = types.BasisStatus;
const ModelError = types.ModelError;

/// A decision variable, identified by index within a model.
pub const Var = struct {
    model: *Model,
    index: usize,

    const Self = @This();

    // ── Generic attribute access ──────────────────────────────────────────

    /// Get a double-typed variable attribute.
    pub fn get(self: Self, attr: Attr) ModelError!f64 {
        return self.model.getDblAttrElement(attr, self.index);
    }

    /// Set a double-typed variable attribute.
    pub fn set(self: *Self, attr: Attr, value: f64) ModelError!void {
        try self.model.setDblAttrElement(attr, self.index, value);
    }

    /// Get an integer-typed variable attribute.
    pub fn getInt(self: Self, attr: Attr) ModelError!i64 {
        return self.model.getIntAttrElement(attr, self.index);
    }

    /// Set an integer-typed variable attribute.
    pub fn setInt(self: *Self, attr: Attr, value: i64) ModelError!void {
        try self.model.setIntAttrElement(attr, self.index, value);
    }

    /// Get a char-typed variable attribute.
    pub fn getChar(self: Self, attr: Attr) ModelError!u8 {
        return self.model.getCharAttrElement(attr, self.index);
    }

    /// Set a char-typed variable attribute.
    pub fn setChar(self: *Self, attr: Attr, value: u8) ModelError!void {
        try self.model.setCharAttrElement(attr, self.index, value);
    }

    /// Get a string-typed variable attribute.
    pub fn getStr(self: Self, attr: Attr) ModelError![]const u8 {
        return self.model.getStrAttrElement(attr, self.index);
    }

    /// Set a string-typed variable attribute.
    pub fn setStr(self: *Self, attr: Attr, value: []const u8) ModelError!void {
        try self.model.setStrAttrElement(attr, self.index, value);
    }

    // ── Convenience accessors (most common attributes) ────────────────────

    pub fn getLB(self: Self) ModelError!f64 { return self.model.getDblAttrElement(.lb, self.index); }
    pub fn setLB(self: *Self, lb: f64) ModelError!void { try self.model.setDblAttrElement(.lb, self.index, lb); }

    pub fn getUB(self: Self) ModelError!f64 { return self.model.getDblAttrElement(.ub, self.index); }
    pub fn setUB(self: *Self, ub: f64) ModelError!void { try self.model.setDblAttrElement(.ub, self.index, ub); }

    pub fn getObj(self: Self) ModelError!f64 { return self.model.getDblAttrElement(.obj, self.index); }
    pub fn setObj(self: *Self, obj: f64) ModelError!void { try self.model.setDblAttrElement(.obj, self.index, obj); }

    pub fn getType(self: Self) ModelError!VarType {
        const code = try self.model.getCharAttrElement(.v_type, self.index);
        return VarType.fromCode(code);
    }
    pub fn setType(self: *Self, vtype: VarType) ModelError!void {
        try self.model.setCharAttrElement(.v_type, self.index, @intFromEnum(vtype));
    }

    /// Primal solution value (attribute `X`), available after `optimize`.
    pub fn getX(self: Self) ModelError!f64 { return self.model.getDblAttrElement(.x, self.index); }

    /// Reduced cost (attribute `RC`), available after `optimize`.
    pub fn getRC(self: Self) ModelError!f64 { return self.model.getDblAttrElement(.rc, self.index); }

    pub fn getVarName(self: Self) ModelError![]const u8 {
        return self.model.getStrAttrElement(.var_name, self.index);
    }
    pub fn setVarName(self: *Self, name: []const u8) ModelError!void {
        try self.model.setStrAttrElement(.var_name, self.index, name);
    }

    pub fn getBasis(self: Self) ModelError!BasisStatus {
        const v = try self.model.getIntAttrElement(.v_basis, self.index);
        return @enumFromInt(@as(i8, @intCast(v)));
    }
    pub fn setBasis(self: *Self, status: BasisStatus) ModelError!void {
        try self.model.setIntAttrElement(.v_basis, self.index, @intFromEnum(status));
    }

    pub fn getStart(self: Self) ModelError!f64 { return self.model.getDblAttrElement(.start, self.index); }
    pub fn setStart(self: *Self, value: f64) ModelError!void { try self.model.setDblAttrElement(.start, self.index, value); }

    pub fn getPStart(self: Self) ModelError!f64 { return self.model.getDblAttrElement(.p_start, self.index); }
    pub fn setPStart(self: *Self, value: f64) ModelError!void { try self.model.setDblAttrElement(.p_start, self.index, value); }

    pub fn getDStart(self: Self) ModelError!f64 { return self.model.getDblAttrElement(.d_start, self.index); }
    pub fn setDStart(self: *Self, value: f64) ModelError!void { try self.model.setDblAttrElement(.d_start, self.index, value); }

    // ── Comparison ────────────────────────────────────────────────────────

    /// Returns `true` when two `Var` values refer to the same model variable.
    pub fn sameAs(self: Self, other: Self) bool {
        return self.model == other.model and self.index == other.index;
    }
};
