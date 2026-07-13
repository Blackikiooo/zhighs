//! Column builder — a sparse column of non-zero coefficients.

const std = @import("std");
const Constr = @import("../constraint/index.zig").Constr;
const ModelError = @import("../types.zig").ModelError;
const ConstrId = @import("../entity_handle.zig").ConstrId;
const Model = @import("../model.zig").Model;

/// A column describing the non-zero entries of a variable in the constraint
/// matrix.  Used in conjunction with `Model.addVar`:
///
/// ```zig
/// var col = Column.init();
/// defer col.deinit(allocator);
/// try col.addTerm(constr_0, 2.5);
/// try col.addTerm(constr_1, -1.0);
/// try model.addVar(col, 0.0, 1.0, 1.0, .continuous, "x");
/// ```
pub const Column = struct {
    /// Constraint indices.
    indices: std.ArrayListUnmanaged(usize) = .empty,
    stable_ids: std.ArrayListUnmanaged(?ConstrId) = .empty,
    /// Coefficient values.
    values: std.ArrayListUnmanaged(f64) = .empty,

    const Self = @This();

    /// Create an empty column.
    pub fn init() Self {
        return .{};
    }

    /// Release memory held by the column.
    pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
        self.indices.deinit(allocator);
        self.stable_ids.deinit(allocator);
        self.values.deinit(allocator);
    }

    /// Add one non-zero `(constraint, coefficient)` pair to the column.
    pub fn addTerm(self: *Self, allocator: std.mem.Allocator, constr: Constr, coeff: f64) ModelError!void {
        self.indices.append(allocator, constr.index) catch return error.OutOfMemory;
        self.stable_ids.append(allocator, constr.id) catch return error.OutOfMemory;
        self.values.append(allocator, coeff) catch return error.OutOfMemory;
    }

    /// Add a term by raw constraint index.
    pub fn addTermByIndex(self: *Self, allocator: std.mem.Allocator, constr_idx: usize, coeff: f64) ModelError!void {
        self.indices.append(allocator, constr_idx) catch return error.OutOfMemory;
        self.stable_ids.append(allocator, null) catch return error.OutOfMemory;
        self.values.append(allocator, coeff) catch return error.OutOfMemory;
    }

    /// Return the number of non-zeros in this column.
    pub fn numNz(self: Self) usize {
        return self.indices.items.len;
    }

    /// Resolve stable constraint handles to the current dense row indices.
    /// The returned array is owned by the caller and is suitable for the
    /// matrix-building API.
    pub fn resolveIndices(self: Self, allocator: std.mem.Allocator, model: Model) ModelError![]usize {
        if (self.stable_ids.items.len != self.indices.items.len) return error.InvalidArgument;
        const resolved = allocator.alloc(usize, self.indices.items.len) catch return error.OutOfMemory;
        errdefer allocator.free(resolved);
        for (self.stable_ids.items, 0..) |maybe_id, i| {
            resolved[i] = if (maybe_id) |id| try model.resolveConstrId(id) else self.indices.items[i];
        }
        return resolved;
    }
};
