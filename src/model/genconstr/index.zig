//! Typed index wrapper for a general (non-linear) constraint.

const std = @import("std");
const Model = @import("../model.zig").Model;
const types = @import("../types.zig");
const Var = @import("../var/index.zig").Var;

const GenConstrType = types.GenConstrType;
const ModelError = types.ModelError;

/// A general (non-linear) constraint, identified by index within a model.
pub const GenConstr = struct {
    model: *Model,
    index: usize,

    const Self = @This();

    pub fn getType(self: Self) ModelError!GenConstrType {
        if (self.index >= self.model.genconstr_count) return error.IndexOutOfRange;
        return self.model.genconstr_types[self.index];
    }

    pub fn getResultVar(self: Self) ModelError!Var {
        if (self.index >= self.model.genconstr_count) return error.IndexOutOfRange;
        return Var{ .model = self.model, .index = self.model.genconstr_resvar[self.index] };
    }

    pub fn getGenConstrName(self: Self) ModelError![]const u8 {
        if (self.index >= self.model.genconstr_count) return error.IndexOutOfRange;
        return self.model.genconstr_names[self.index] orelse "";
    }

    pub fn sameAs(self: Self, other: Self) bool {
        return self.model == other.model and self.index == other.index;
    }
};
