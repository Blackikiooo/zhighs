//! Typed index wrapper for a Special Ordered Set constraint.

const std = @import("std");
const Model = @import("../model.zig").Model;
const types = @import("../types.zig");

const SosType = types.SosType;
const ModelError = types.ModelError;
const SosId = @import("../entity_handle.zig").SosId;

/// A Special Ordered Set constraint, identified by index within a model.
pub const SOS = struct {
    model: *Model,
    index: usize,
    id: ?SosId = null,

    const Self = @This();

    pub fn at(model: *Model, index: usize) ModelError!Self {
        return .{ .model = model, .index = index, .id = try model.sosIdAt(index) };
    }

    fn denseIndex(self: Self) ModelError!usize {
        if (self.id) |id| return self.model.resolveSosId(id);
        return self.index;
    }

    pub fn getType(self: Self) ModelError!SosType {
        const index = try self.denseIndex();
        if (index >= self.model.sos_count) return error.IndexOutOfRange;
        return self.model.sos_types[index];
    }

    pub fn getSOSName(self: Self) ModelError![]const u8 {
        const index = try self.denseIndex();
        if (index >= self.model.sos_count) return error.IndexOutOfRange;
        return self.model.sos_names[index] orelse "";
    }

    pub fn numMembers(self: Self) ModelError!usize {
        const index = try self.denseIndex();
        if (index >= self.model.sos_count) return error.IndexOutOfRange;
        const beg = self.model.sos_begin[index];
        const end = self.model.sos_begin[index + 1];
        return end - beg;
    }

    pub fn sameAs(self: Self, other: Self) bool {
        return self.model == other.model and self.index == other.index;
    }
};
