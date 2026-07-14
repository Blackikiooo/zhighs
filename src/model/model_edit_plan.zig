//! Normalized data-oriented execution plan for one pending model-edit segment.

const std = @import("std");
const pending_module = @import("model_pending.zig");
const types = @import("types.zig");

const PendingChange = pending_module.PendingChange;

/// Small scalar segments stay on the allocation-free direct path. Larger
/// segments amortize DOD plan construction through duplicate coalescing.
pub const small_direct_edit_threshold: usize = 8;

pub fn isDirectScalarSegment(pending: []const PendingChange) bool {
    if (pending.len == 0 or pending.len > small_direct_edit_threshold) return false;
    for (pending) |change| switch (change) {
        .chg_bounds, .chg_obj, .chg_rhs, .chg_sense, .chg_type => {},
        else => return false,
    };
    return true;
}

pub const PlanKind = enum {
    empty,
    scalar_only,
    coefficients_only,
    mixed_nonstructural,
    structural,
};

pub const CoefficientEdit = struct { row: usize, col: usize, value: f64, sequence: usize };
pub const BoundsEdit = struct { index: usize, lower: f64, upper: f64, sequence: usize };
pub const ObjectiveEdit = struct { index: usize, value: f64, sequence: usize };
pub const RhsEdit = struct { index: usize, value: f64, sequence: usize };
pub const SenseEdit = struct { index: usize, value: types.Sense, sequence: usize };
pub const TypeEdit = struct { index: usize, value: types.VarType, sequence: usize };

/// First-stage edit plan. Structural payloads remain owned by `PendingChange`
/// for now; all coefficient and scalar streams are normalized into contiguous
/// SoA fields and duplicate targets use last-write-wins.
pub const ModelEditPlan = struct {
    kind: PlanKind = .empty,
    coefficients: std.MultiArrayList(CoefficientEdit) = .empty,
    bounds: std.MultiArrayList(BoundsEdit) = .empty,
    objective: std.MultiArrayList(ObjectiveEdit) = .empty,
    rhs: std.MultiArrayList(RhsEdit) = .empty,
    senses: std.MultiArrayList(SenseEdit) = .empty,
    types: std.MultiArrayList(TypeEdit) = .empty,
    deleted_vars: std.ArrayListUnmanaged(usize) = .empty,
    deleted_constraints: std.ArrayListUnmanaged(usize) = .empty,
    has_structure: bool = false,

    pub fn build(allocator: std.mem.Allocator, pending: []const PendingChange) std.mem.Allocator.Error!ModelEditPlan {
        var result = ModelEditPlan{};
        errdefer result.deinit(allocator);

        var coefficient_count: usize = 0;
        var bounds_count: usize = 0;
        var objective_count: usize = 0;
        var rhs_count: usize = 0;
        var sense_count: usize = 0;
        var type_count: usize = 0;
        var deleted_var_count: usize = 0;
        var deleted_constraint_count: usize = 0;
        for (pending) |change| switch (change) {
            .chg_coeff => coefficient_count += 1,
            .chg_bounds => bounds_count += 1,
            .chg_obj => objective_count += 1,
            .chg_rhs => rhs_count += 1,
            .chg_sense => sense_count += 1,
            .chg_type => {
                type_count += 1;
                result.has_structure = true;
            },
            .del_vars => |edit| {
                deleted_var_count += edit.indices.len;
                result.has_structure = true;
            },
            .del_constrs => |edit| {
                deleted_constraint_count += edit.indices.len;
                result.has_structure = true;
            },
            .add_var, .add_constr => result.has_structure = true,
        };
        try result.coefficients.ensureUnusedCapacity(allocator, coefficient_count);
        try result.bounds.ensureUnusedCapacity(allocator, bounds_count);
        try result.objective.ensureUnusedCapacity(allocator, objective_count);
        try result.rhs.ensureUnusedCapacity(allocator, rhs_count);
        try result.senses.ensureUnusedCapacity(allocator, sense_count);
        try result.types.ensureUnusedCapacity(allocator, type_count);
        try result.deleted_vars.ensureUnusedCapacity(allocator, deleted_var_count);
        try result.deleted_constraints.ensureUnusedCapacity(allocator, deleted_constraint_count);

        for (pending, 0..) |change, sequence| switch (change) {
            .chg_coeff => |edit| result.coefficients.appendAssumeCapacity(.{
                .row = edit.constr_idx,
                .col = edit.var_idx,
                .value = edit.new_val,
                .sequence = sequence,
            }),
            .chg_bounds => |edit| result.bounds.appendAssumeCapacity(.{
                .index = edit.var_idx,
                .lower = edit.lb,
                .upper = edit.ub,
                .sequence = sequence,
            }),
            .chg_obj => |edit| result.objective.appendAssumeCapacity(.{ .index = edit.var_idx, .value = edit.obj, .sequence = sequence }),
            .chg_rhs => |edit| result.rhs.appendAssumeCapacity(.{ .index = edit.constr_idx, .value = edit.rhs, .sequence = sequence }),
            .chg_sense => |edit| result.senses.appendAssumeCapacity(.{ .index = edit.constr_idx, .value = edit.sense, .sequence = sequence }),
            .chg_type => |edit| result.types.appendAssumeCapacity(.{ .index = edit.var_idx, .value = edit.vtype, .sequence = sequence }),
            .del_vars => |edit| for (edit.indices) |index| result.deleted_vars.appendAssumeCapacity(index),
            .del_constrs => |edit| for (edit.indices) |index| result.deleted_constraints.appendAssumeCapacity(index),
            else => {},
        };

        normalizeCoefficients(&result.coefficients);
        normalizeByIndex(BoundsEdit, &result.bounds);
        normalizeByIndex(ObjectiveEdit, &result.objective);
        normalizeByIndex(RhsEdit, &result.rhs);
        normalizeByIndex(SenseEdit, &result.senses);
        normalizeByIndex(TypeEdit, &result.types);
        normalizeIds(&result.deleted_vars);
        normalizeIds(&result.deleted_constraints);
        result.kind = result.classify();
        return result;
    }

    pub fn deinit(self: *ModelEditPlan, allocator: std.mem.Allocator) void {
        self.coefficients.deinit(allocator);
        self.bounds.deinit(allocator);
        self.objective.deinit(allocator);
        self.rhs.deinit(allocator);
        self.senses.deinit(allocator);
        self.types.deinit(allocator);
        self.deleted_vars.deinit(allocator);
        self.deleted_constraints.deinit(allocator);
        self.* = undefined;
    }

    fn classify(self: ModelEditPlan) PlanKind {
        if (self.has_structure) return .structural;
        const has_coefficients = self.coefficients.len != 0;
        const has_scalars = self.bounds.len != 0 or self.objective.len != 0 or self.rhs.len != 0 or self.senses.len != 0 or self.types.len != 0;
        if (has_coefficients and has_scalars) return .mixed_nonstructural;
        if (has_coefficients) return .coefficients_only;
        if (has_scalars) return .scalar_only;
        return .empty;
    }
};

fn normalizeIds(ids: *std.ArrayListUnmanaged(usize)) void {
    if (ids.items.len < 2) return;
    std.sort.pdq(usize, ids.items, {}, lessThanUsize);
    var write: usize = 1;
    for (ids.items[1..]) |index| {
        if (index == ids.items[write - 1]) continue;
        ids.items[write] = index;
        write += 1;
    }
    ids.shrinkRetainingCapacity(write);
}

fn lessThanUsize(_: void, lhs: usize, rhs: usize) bool {
    return lhs < rhs;
}

const CoefficientSortContext = struct {
    rows: []usize,
    cols: []usize,
    sequences: []usize,

    pub fn lessThan(self: @This(), lhs: usize, rhs: usize) bool {
        if (self.cols[lhs] != self.cols[rhs]) return self.cols[lhs] < self.cols[rhs];
        if (self.rows[lhs] != self.rows[rhs]) return self.rows[lhs] < self.rows[rhs];
        return self.sequences[lhs] < self.sequences[rhs];
    }
};

fn normalizeCoefficients(edits: *std.MultiArrayList(CoefficientEdit)) void {
    if (edits.len < 2) return;
    var fields = edits.slice();
    edits.sort(CoefficientSortContext{
        .rows = fields.items(.row),
        .cols = fields.items(.col),
        .sequences = fields.items(.sequence),
    });
    fields = edits.slice();
    const rows = fields.items(.row);
    const cols = fields.items(.col);
    var read: usize = 0;
    var write: usize = 0;
    while (read < edits.len) {
        var last = read;
        read += 1;
        while (read < edits.len and rows[read] == rows[last] and cols[read] == cols[last]) : (read += 1) last = read;
        fields.set(write, fields.get(last));
        write += 1;
    }
    edits.shrinkRetainingCapacity(write);
}

fn normalizeByIndex(comptime Entry: type, edits: *std.MultiArrayList(Entry)) void {
    if (edits.len < 2) return;
    var fields = edits.slice();
    edits.sort(IndexSortContext{
        .indices = fields.items(.index),
        .sequences = fields.items(.sequence),
    });
    fields = edits.slice();
    const indices = fields.items(.index);
    var read: usize = 0;
    var write: usize = 0;
    while (read < edits.len) {
        var last = read;
        read += 1;
        while (read < edits.len and indices[read] == indices[last]) : (read += 1) last = read;
        fields.set(write, fields.get(last));
        write += 1;
    }
    edits.shrinkRetainingCapacity(write);
}

const IndexSortContext = struct {
    indices: []usize,
    sequences: []usize,

    pub fn lessThan(self: @This(), lhs: usize, rhs: usize) bool {
        if (self.indices[lhs] != self.indices[rhs]) return self.indices[lhs] < self.indices[rhs];
        return self.sequences[lhs] < self.sequences[rhs];
    }
};

test "edit plan normalizes coefficient and scalar streams with last write wins" {
    const pending = [_]PendingChange{
        .{ .chg_coeff = .{ .constr_idx = 2, .var_idx = 1, .new_val = 3.0 } },
        .{ .chg_bounds = .{ .var_idx = 4, .lb = 0.0, .ub = 1.0 } },
        .{ .chg_coeff = .{ .constr_idx = 2, .var_idx = 1, .new_val = 7.0 } },
        .{ .chg_obj = .{ .var_idx = 3, .obj = 2.0 } },
        .{ .chg_bounds = .{ .var_idx = 4, .lb = -1.0, .ub = 5.0 } },
    };
    var plan = try ModelEditPlan.build(std.testing.allocator, &pending);
    defer plan.deinit(std.testing.allocator);
    try std.testing.expectEqual(PlanKind.mixed_nonstructural, plan.kind);
    try std.testing.expectEqual(@as(usize, 1), plan.coefficients.len);
    try std.testing.expectEqual(@as(f64, 7.0), plan.coefficients.items(.value)[0]);
    try std.testing.expectEqual(@as(usize, 1), plan.bounds.len);
    try std.testing.expectEqual(@as(f64, -1.0), plan.bounds.items(.lower)[0]);
    try std.testing.expectEqual(@as(f64, 5.0), plan.bounds.items(.upper)[0]);
}

test "small direct scalar segment threshold rejects matrix and structural edits" {
    const scalar = [_]PendingChange{
        .{ .chg_obj = .{ .var_idx = 0, .obj = 1.0 } },
        .{ .chg_rhs = .{ .constr_idx = 0, .rhs = 2.0 } },
    };
    try std.testing.expect(isDirectScalarSegment(&scalar));
    const coefficient = [_]PendingChange{.{ .chg_coeff = .{ .constr_idx = 0, .var_idx = 0, .new_val = 1.0 } }};
    try std.testing.expect(!isDirectScalarSegment(&coefficient));
    var large = [_]PendingChange{.{ .chg_obj = .{ .var_idx = 0, .obj = 1.0 } }} ** (small_direct_edit_threshold + 1);
    try std.testing.expect(!isDirectScalarSegment(&large));
}

test "edit plan sorts and deduplicates deleted dense IDs" {
    var first = [_]usize{ 4, 1, 4 };
    var second = [_]usize{ 3, 1 };
    const pending = [_]PendingChange{
        .{ .del_vars = .{ .indices = &first } },
        .{ .del_vars = .{ .indices = &second } },
    };
    var plan = try ModelEditPlan.build(std.testing.allocator, &pending);
    defer plan.deinit(std.testing.allocator);
    try std.testing.expectEqualSlices(usize, &.{ 1, 3, 4 }, plan.deleted_vars.items);
    try std.testing.expectEqual(PlanKind.structural, plan.kind);
}
