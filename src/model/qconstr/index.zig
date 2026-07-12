//! Typed index wrapper for a quadratic constraint.

const std = @import("std");
const Model = @import("../model.zig").Model;
const types = @import("../types.zig");

const Sense = types.Sense;
const ModelError = types.ModelError;

/// A quadratic constraint, identified by index within a model.
pub const QConstr = struct {
    model: *Model,
    index: usize,

    const Self = @This();

    pub fn getSense(self: Self) ModelError!Sense {
        if (self.index >= self.model.qconstr_count) return error.IndexOutOfRange;
        return self.model.qconstr_sense[self.index];
    }

    pub fn getRHS(self: Self) ModelError!f64 {
        if (self.index >= self.model.qconstr_count) return error.IndexOutOfRange;
        return self.model.qconstr_rhs[self.index];
    }

    pub fn getQConstrName(self: Self) ModelError![]const u8 {
        if (self.index >= self.model.qconstr_count) return error.IndexOutOfRange;
        return self.model.qconstr_names[self.index] orelse "";
    }

    pub fn sameAs(self: Self, other: Self) bool {
        return self.model == other.model and self.index == other.index;
    }
};
