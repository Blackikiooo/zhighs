//! Managed container for quadratic constraint data.
//!
//! SoA for fixed fields (sense, rhs, name) + packed triplet storage
//! for the quadratic and linear terms.

const std = @import("std");
const QConstrData = @import("data.zig").QConstrData;

/// Container for quadratic constraints with packed term storage.
pub const QConstrArray = struct {
    /// SoA for fixed per-QConstr fields.
    inner: std.MultiArrayList(QConstrData) = .{},

    /// Packed quadratic triplets (CSR-like offsets into qrow/qcol/qval).
    q_begin: std.ArrayListUnmanaged(usize) = .{},
    qrow: std.ArrayListUnmanaged(i32) = .{},
    qcol: std.ArrayListUnmanaged(i32) = .{},
    qval: std.ArrayListUnmanaged(f64) = .{},

    /// Packed linear terms.
    l_begin: std.ArrayListUnmanaged(usize) = .{},
    lind: std.ArrayListUnmanaged(usize) = .{},
    lval: std.ArrayListUnmanaged(f64) = .{},

    pub inline fn len(self: QConstrArray) usize {
        return self.inner.len;
    }

    pub inline fn get(self: QConstrArray, index: usize) QConstrData {
        return self.inner.get(index);
    }

    pub fn deinit(self: *QConstrArray, allocator: std.mem.Allocator) void {
        for (self.inner.items(.name)) |n| if (n) |s| allocator.free(s);
        self.inner.deinit(allocator);
        self.q_begin.deinit(allocator);
        self.qrow.deinit(allocator);
        self.qcol.deinit(allocator);
        self.qval.deinit(allocator);
        self.l_begin.deinit(allocator);
        self.lind.deinit(allocator);
        self.lval.deinit(allocator);
    }
};
