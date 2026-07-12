//! Pending-change type and queue management for `Model`.
//!
//! Modifications are queued and applied in batch by `updateModel`.  This module
//! owns the change record type and the low-level queue primitives.

const std = @import("std");
const types = @import("types.zig");
const Model = @import("model.zig").Model;

const ModelError = types.ModelError;
const VarType = types.VarType;
const Sense = types.Sense;

/// Upper bound on pending change count before an auto-flush is forced.
pub const PENDING_FLUSH_THRESHOLD: usize = std.math.maxInt(u16) + 1; // 65536

/// One atomic change queued for the next `updateModel` call.
pub const PendingChange = union(enum) {
    add_var: struct {
        num_nz: usize,
        vind: []const usize,
        vval: []const f64,
        obj: f64,
        lb: f64,
        ub: f64,
        vtype: VarType,
        name: ?[]const u8,
    },
    add_constr: struct {
        num_nz: usize,
        cind: []const usize,
        cval: []const f64,
        sense: Sense,
        rhs: f64,
        name: ?[]const u8,
    },
    del_vars: struct { indices: []const usize },
    del_constrs: struct { indices: []const usize },
    chg_coeff: struct { constr_idx: usize, var_idx: usize, new_val: f64 },
    chg_bounds: struct { var_idx: usize, lb: f64, ub: f64 },
    chg_obj: struct { var_idx: usize, obj: f64 },
    chg_rhs: struct { constr_idx: usize, rhs: f64 },
    chg_sense: struct { constr_idx: usize, sense: Sense },
    chg_type: struct { var_idx: usize, vtype: VarType },
};

/// Queue a change for the next `updateModel` / `applyPending` batch.
pub fn enqueue(self: *Model, change: PendingChange) ModelError!void {
    self.pending.append(self.allocator, change) catch return error.OutOfMemory;
    self.has_pending = true;

    // Auto-flush if the queue grows too large.
    if (self.pending.items.len >= PENDING_FLUSH_THRESHOLD) {
        try self.applyPending();
    }
}

/// Discard all queued modifications without applying them.
pub fn discardPending(self: *Model) void {
    for (self.pending.items) |*chg| {
        switch (chg.*) {
            .add_var => |v| {
                self.allocator.free(v.vind);
                self.allocator.free(v.vval);
                if (v.name) |n| self.allocator.free(n);
            },
            .add_constr => |c| {
                self.allocator.free(c.cind);
                self.allocator.free(c.cval);
                if (c.name) |n| self.allocator.free(n);
            },
            .del_vars => |d| self.allocator.free(d.indices),
            .del_constrs => |d| self.allocator.free(d.indices),
            else => {},
        }
    }
    self.pending.clearAndFree(self.allocator);
    self.has_pending = false;
}

/// Count how many `add_var` entries are queued.
pub fn countPendingAddVar(self: Model) usize {
    var count: usize = 0;
    for (self.pending.items) |chg| {
        switch (chg) {
            .add_var => count += 1,
            else => {},
        }
    }
    return count;
}

/// Count how many `add_constr` entries are queued.
pub fn countPendingAddConstr(self: Model) usize {
    var count: usize = 0;
    for (self.pending.items) |chg| {
        switch (chg) {
            .add_constr => count += 1,
            else => {},
        }
    }
    return count;
}
