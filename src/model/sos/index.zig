//! Typed index wrapper for a Special Ordered Set constraint.

const std = @import("std");
const Model = @import("../model.zig").Model;
const types = @import("../types.zig");

const SosType = types.SosType;
const ModelError = types.ModelError;

/// A Special Ordered Set constraint, identified by index within a model.
pub const SOS = struct {
    model: *Model,
    index: usize,

    const Self = @This();

    pub fn getType(self: Self) ModelError!SosType {
        if (self.index >= self.model.sos_count) return error.IndexOutOfRange;
        return self.model.sos_types[self.index];
    }

    pub fn getSOSName(self: Self) ModelError![]const u8 {
        if (self.index >= self.model.sos_count) return error.IndexOutOfRange;
        return self.model.sos_names[self.index] orelse "";
    }

    pub fn numMembers(self: Self) ModelError!usize {
        if (self.index >= self.model.sos_count) return error.IndexOutOfRange;
        const beg = self.model.sos_begin[self.index];
        const end = self.model.sos_begin[self.index + 1];
        return end - beg;
    }

    pub fn sameAs(self: Self, other: Self) bool {
        return self.model == other.model and self.index == other.index;
    }
};
