//! Typed index wrapper for a general (non-linear) constraint.

const std = @import("std");
const Model = @import("../model.zig").Model;
const types = @import("../types.zig");
const Var = @import("../var/index.zig").Var;

const GenConstrType = types.GenConstrType;
const ModelError = types.ModelError;
const GenConstrId = @import("../entity_handle.zig").GenConstrId;

/// A general (non-linear) constraint, identified by index within a model.
pub const GenConstr = struct {
    model: *Model,
    index: usize,
    id: ?GenConstrId = null,

    const Self = @This();

    pub fn at(model: *Model, index: usize) ModelError!Self {
        return .{ .model = model, .index = index, .id = try model.genconstrIdAt(index) };
    }

    fn denseIndex(self: Self) ModelError!usize {
        if (self.id) |id| return self.model.resolveGenConstrId(id);
        return self.index;
    }

    pub fn getType(self: Self) ModelError!GenConstrType {
        const index = try self.denseIndex();
        if (index >= self.model.genconstr_count) return error.IndexOutOfRange;
        return self.model.genconstr_types[index];
    }

    pub fn getResultVar(self: Self) ModelError!Var {
        const index = try self.denseIndex();
        if (index >= self.model.genconstr_count) return error.IndexOutOfRange;
        return Var{ .model = self.model, .index = self.model.genconstr_resvar[index] };
    }

    pub fn getGenConstrName(self: Self) ModelError![]const u8 {
        const index = try self.denseIndex();
        if (index >= self.model.genconstr_count) return error.IndexOutOfRange;
        return self.model.genconstr_names[index] orelse "";
    }

    pub fn sameAs(self: Self, other: Self) bool {
        return self.model == other.model and self.index == other.index;
    }
};
