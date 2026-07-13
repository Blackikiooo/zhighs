//! Typed index wrapper for a quadratic constraint.

const std = @import("std");
const Model = @import("../model.zig").Model;
const types = @import("../types.zig");

const Sense = types.Sense;
const ModelError = types.ModelError;
const QConstrId = @import("../entity_handle.zig").QConstrId;

/// A quadratic constraint, identified by index within a model.
pub const QConstr = struct {
    model: *Model,
    index: usize,
    id: ?QConstrId = null,

    const Self = @This();

    pub fn at(model: *Model, index: usize) ModelError!Self {
        return .{ .model = model, .index = index, .id = try model.qconstrIdAt(index) };
    }

    fn denseIndex(self: Self) ModelError!usize {
        if (self.id) |id| return self.model.resolveQConstrId(id);
        return self.index;
    }

    pub fn getSense(self: Self) ModelError!Sense {
        const index = try self.denseIndex();
        if (index >= self.model.qconstr_count) return error.IndexOutOfRange;
        return self.model.qconstr_sense[index];
    }

    pub fn getRHS(self: Self) ModelError!f64 {
        const index = try self.denseIndex();
        if (index >= self.model.qconstr_count) return error.IndexOutOfRange;
        return self.model.qconstr_rhs[index];
    }

    pub fn getQConstrName(self: Self) ModelError![]const u8 {
        const index = try self.denseIndex();
        if (index >= self.model.qconstr_count) return error.IndexOutOfRange;
        return self.model.qconstr_names[index] orelse "";
    }

    pub fn sameAs(self: Self, other: Self) bool {
        return self.model == other.model and self.index == other.index;
    }
};
