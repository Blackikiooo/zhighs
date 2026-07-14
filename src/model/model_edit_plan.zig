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
pub const AddedColumn = struct {
    objective: f64,
    lower: f64,
    upper: f64,
    var_type: types.VarType,
    name: ?[]const u8,
};
pub const AddedRow = struct { sense: types.Sense, rhs: f64, name: ?[]const u8 };

/// Data-oriented execution plan for one pending edit segment. Added columns
/// and rows are flattened into offset/index/value streams; coefficient and
/// scalar streams are normalized with last-write-wins semantics.
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
    added_columns: std.MultiArrayList(AddedColumn) = .empty,
    added_column_starts: std.ArrayListUnmanaged(usize) = .empty,
    added_column_rows: std.ArrayListUnmanaged(usize) = .empty,
    added_column_values: std.ArrayListUnmanaged(f64) = .empty,
    added_rows: std.MultiArrayList(AddedRow) = .empty,
    added_row_starts: std.ArrayListUnmanaged(usize) = .empty,
    added_row_columns: std.ArrayListUnmanaged(usize) = .empty,
    added_row_values: std.ArrayListUnmanaged(f64) = .empty,
    has_structure: bool = false,
    changes_matrix_values: bool = false,
    changes_bounds: bool = false,
    changes_objective: bool = false,

    pub fn build(allocator: std.mem.Allocator, pending: []const PendingChange, existing_vars: usize, existing_rows: usize) std.mem.Allocator.Error!ModelEditPlan {
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
        var added_column_count: usize = 0;
        var added_column_nnz: usize = 0;
        var added_row_count: usize = 0;
        var added_row_nnz: usize = 0;
        for (pending) |change| switch (change) {
            .chg_coeff => {
                coefficient_count += 1;
                result.changes_matrix_values = true;
            },
            .chg_bounds => {
                bounds_count += 1;
                result.changes_bounds = true;
            },
            .chg_obj => {
                objective_count += 1;
                result.changes_objective = true;
            },
            .chg_rhs => {
                rhs_count += 1;
                result.changes_bounds = true;
            },
            .chg_sense => {
                sense_count += 1;
                result.changes_bounds = true;
            },
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
            .add_var => |edit| {
                added_column_count += 1;
                added_column_nnz += edit.num_nz;
                result.has_structure = true;
            },
            .add_constr => |edit| {
                added_row_count += 1;
                added_row_nnz += edit.num_nz;
                result.has_structure = true;
            },
        };
        try result.coefficients.ensureUnusedCapacity(allocator, coefficient_count);
        try result.bounds.ensureUnusedCapacity(allocator, bounds_count);
        try result.objective.ensureUnusedCapacity(allocator, objective_count);
        try result.rhs.ensureUnusedCapacity(allocator, rhs_count);
        try result.senses.ensureUnusedCapacity(allocator, sense_count);
        try result.types.ensureUnusedCapacity(allocator, type_count);
        try result.deleted_vars.ensureUnusedCapacity(allocator, deleted_var_count);
        try result.deleted_constraints.ensureUnusedCapacity(allocator, deleted_constraint_count);
        try result.added_columns.ensureUnusedCapacity(allocator, added_column_count);
        try result.added_column_starts.ensureUnusedCapacity(allocator, added_column_count + @intFromBool(added_column_count != 0));
        try result.added_column_rows.ensureUnusedCapacity(allocator, added_column_nnz);
        try result.added_column_values.ensureUnusedCapacity(allocator, added_column_nnz);
        try result.added_rows.ensureUnusedCapacity(allocator, added_row_count);
        try result.added_row_starts.ensureUnusedCapacity(allocator, added_row_count + @intFromBool(added_row_count != 0));
        try result.added_row_columns.ensureUnusedCapacity(allocator, added_row_nnz);
        try result.added_row_values.ensureUnusedCapacity(allocator, added_row_nnz);
        if (added_column_count != 0) result.added_column_starts.appendAssumeCapacity(0);
        if (added_row_count != 0) result.added_row_starts.appendAssumeCapacity(0);

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
            .add_var => |edit| {
                result.added_columns.appendAssumeCapacity(.{
                    .objective = edit.obj,
                    .lower = edit.lb,
                    .upper = edit.ub,
                    .var_type = edit.vtype,
                    .name = edit.name,
                });
                for (edit.vind[0..edit.num_nz], edit.vval[0..edit.num_nz]) |row, value| {
                    result.added_column_rows.appendAssumeCapacity(row);
                    result.added_column_values.appendAssumeCapacity(value);
                }
                result.added_column_starts.appendAssumeCapacity(result.added_column_rows.items.len);
            },
            .add_constr => |edit| {
                result.added_rows.appendAssumeCapacity(.{ .sense = edit.sense, .rhs = edit.rhs, .name = edit.name });
                for (edit.cind[0..edit.num_nz], edit.cval[0..edit.num_nz]) |column, value| {
                    result.added_row_columns.appendAssumeCapacity(column);
                    result.added_row_values.appendAssumeCapacity(value);
                }
                result.added_row_starts.appendAssumeCapacity(result.added_row_columns.items.len);
            },
        };

        normalizeCoefficients(&result.coefficients);
        normalizeByIndex(BoundsEdit, &result.bounds);
        normalizeByIndex(ObjectiveEdit, &result.objective);
        normalizeByIndex(RhsEdit, &result.rhs);
        normalizeByIndex(SenseEdit, &result.senses);
        normalizeByIndex(TypeEdit, &result.types);
        normalizeIds(&result.deleted_vars);
        normalizeIds(&result.deleted_constraints);
        result.foldNewObjectScalars(existing_vars, existing_rows);
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
        self.added_columns.deinit(allocator);
        self.added_column_starts.deinit(allocator);
        self.added_column_rows.deinit(allocator);
        self.added_column_values.deinit(allocator);
        self.added_rows.deinit(allocator);
        self.added_row_starts.deinit(allocator);
        self.added_row_columns.deinit(allocator);
        self.added_row_values.deinit(allocator);
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

    fn foldNewObjectScalars(self: *ModelEditPlan, existing_vars: usize, existing_rows: usize) void {
        if (self.added_columns.len != 0) {
            var column_fields = self.added_columns.slice();
            self.foldBounds(existing_vars, &column_fields);
            self.foldObjective(existing_vars, &column_fields);
            self.foldTypes(existing_vars, &column_fields);
        }
        if (self.added_rows.len != 0) {
            var row_fields = self.added_rows.slice();
            self.foldRhs(existing_rows, &row_fields);
            self.foldSenses(existing_rows, &row_fields);
        }
    }

    fn foldBounds(self: *ModelEditPlan, existing: usize, added: *std.MultiArrayList(AddedColumn).Slice) void {
        var fields = self.bounds.slice();
        var write: usize = 0;
        for (0..self.bounds.len) |read| {
            const edit = fields.get(read);
            if (edit.index >= existing and edit.index - existing < self.added_columns.len) {
                const target = edit.index - existing;
                added.items(.lower)[target] = edit.lower;
                added.items(.upper)[target] = edit.upper;
            } else {
                fields.set(write, edit);
                write += 1;
            }
        }
        self.bounds.shrinkRetainingCapacity(write);
    }

    fn foldObjective(self: *ModelEditPlan, existing: usize, added: *std.MultiArrayList(AddedColumn).Slice) void {
        var fields = self.objective.slice();
        var write: usize = 0;
        for (0..self.objective.len) |read| {
            const edit = fields.get(read);
            if (edit.index >= existing and edit.index - existing < self.added_columns.len)
                added.items(.objective)[edit.index - existing] = edit.value
            else {
                fields.set(write, edit);
                write += 1;
            }
        }
        self.objective.shrinkRetainingCapacity(write);
    }

    fn foldTypes(self: *ModelEditPlan, existing: usize, added: *std.MultiArrayList(AddedColumn).Slice) void {
        var fields = self.types.slice();
        var write: usize = 0;
        for (0..self.types.len) |read| {
            const edit = fields.get(read);
            if (edit.index >= existing and edit.index - existing < self.added_columns.len)
                added.items(.var_type)[edit.index - existing] = edit.value
            else {
                fields.set(write, edit);
                write += 1;
            }
        }
        self.types.shrinkRetainingCapacity(write);
    }

    fn foldRhs(self: *ModelEditPlan, existing: usize, added: *std.MultiArrayList(AddedRow).Slice) void {
        var fields = self.rhs.slice();
        var write: usize = 0;
        for (0..self.rhs.len) |read| {
            const edit = fields.get(read);
            if (edit.index >= existing and edit.index - existing < self.added_rows.len)
                added.items(.rhs)[edit.index - existing] = edit.value
            else {
                fields.set(write, edit);
                write += 1;
            }
        }
        self.rhs.shrinkRetainingCapacity(write);
    }

    fn foldSenses(self: *ModelEditPlan, existing: usize, added: *std.MultiArrayList(AddedRow).Slice) void {
        var fields = self.senses.slice();
        var write: usize = 0;
        for (0..self.senses.len) |read| {
            const edit = fields.get(read);
            if (edit.index >= existing and edit.index - existing < self.added_rows.len)
                added.items(.sense)[edit.index - existing] = edit.value
            else {
                fields.set(write, edit);
                write += 1;
            }
        }
        self.senses.shrinkRetainingCapacity(write);
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
    var plan = try ModelEditPlan.build(std.testing.allocator, &pending, 5, 3);
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
    var plan = try ModelEditPlan.build(std.testing.allocator, &pending, 5, 0);
    defer plan.deinit(std.testing.allocator);
    try std.testing.expectEqualSlices(usize, &.{ 1, 3, 4 }, plan.deleted_vars.items);
    try std.testing.expectEqual(PlanKind.structural, plan.kind);
}

test "edit plan flattens structural payloads and folds new object scalars" {
    const column_rows = [_]usize{ 0, 3 };
    const column_values = [_]f64{ 1.5, -2.0 };
    const row_columns = [_]usize{ 1, 2 };
    const row_values = [_]f64{ 4.0, 8.0 };
    const pending = [_]PendingChange{
        .{ .add_var = .{
            .num_nz = column_rows.len,
            .vind = &column_rows,
            .vval = &column_values,
            .obj = 1.0,
            .lb = 0.0,
            .ub = 10.0,
            .vtype = .continuous,
            .name = "new-column",
        } },
        .{ .add_constr = .{
            .num_nz = row_columns.len,
            .cind = &row_columns,
            .cval = &row_values,
            .sense = .less_equal,
            .rhs = 5.0,
            .name = "new-row",
        } },
        .{ .chg_bounds = .{ .var_idx = 2, .lb = -3.0, .ub = 7.0 } },
        .{ .chg_obj = .{ .var_idx = 2, .obj = 9.0 } },
        .{ .chg_type = .{ .var_idx = 2, .vtype = .integer } },
        .{ .chg_rhs = .{ .constr_idx = 3, .rhs = 12.0 } },
        .{ .chg_sense = .{ .constr_idx = 3, .sense = .equal } },
    };

    var plan = try ModelEditPlan.build(std.testing.allocator, &pending, 2, 3);
    defer plan.deinit(std.testing.allocator);

    try std.testing.expectEqual(PlanKind.structural, plan.kind);
    try std.testing.expectEqualSlices(usize, &.{ 0, 2 }, plan.added_column_starts.items);
    try std.testing.expectEqualSlices(usize, &column_rows, plan.added_column_rows.items);
    try std.testing.expectEqualSlices(f64, &column_values, plan.added_column_values.items);
    try std.testing.expectEqualSlices(usize, &.{ 0, 2 }, plan.added_row_starts.items);
    try std.testing.expectEqualSlices(usize, &row_columns, plan.added_row_columns.items);
    try std.testing.expectEqualSlices(f64, &row_values, plan.added_row_values.items);

    const columns = plan.added_columns.slice();
    try std.testing.expectEqual(@as(f64, -3.0), columns.items(.lower)[0]);
    try std.testing.expectEqual(@as(f64, 7.0), columns.items(.upper)[0]);
    try std.testing.expectEqual(@as(f64, 9.0), columns.items(.objective)[0]);
    try std.testing.expectEqual(types.VarType.integer, columns.items(.var_type)[0]);
    try std.testing.expectEqualStrings("new-column", columns.items(.name)[0].?);
    const rows = plan.added_rows.slice();
    try std.testing.expectEqual(@as(f64, 12.0), rows.items(.rhs)[0]);
    try std.testing.expectEqual(types.Sense.equal, rows.items(.sense)[0]);
    try std.testing.expectEqualStrings("new-row", rows.items(.name)[0].?);

    try std.testing.expectEqual(@as(usize, 0), plan.bounds.len);
    try std.testing.expectEqual(@as(usize, 0), plan.objective.len);
    try std.testing.expectEqual(@as(usize, 0), plan.types.len);
    try std.testing.expectEqual(@as(usize, 0), plan.rhs.len);
    try std.testing.expectEqual(@as(usize, 0), plan.senses.len);
    try std.testing.expect(plan.changes_bounds);
    try std.testing.expect(plan.changes_objective);
}
