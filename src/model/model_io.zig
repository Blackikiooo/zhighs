//! Thin adapter between the public `Model` API and the independent `io` module.
//!
//! File grammars, buffering, suffix dispatch, and canonical sparse assembly
//! belong to `io`. This file only maps types, preserves lazy-update semantics,
//! and atomically publishes a successfully imported model.

const std = @import("std");
const io = @import("io");
const types = @import("types.zig");
const Model = @import("model.zig").Model;

const ModelError = types.ModelError;

pub fn writeModel(self: *Model, filename: []const u8) ModelError!void {
    try self.updateModel();
    const col_types = self.allocator.alloc(io.VariableType, self.num_vars) catch return error.OutOfMemory;
    defer self.allocator.free(col_types);
    for (self.var_type, col_types) |source, *target| target.* = switch (source) {
        .continuous => .continuous,
        .binary => .binary,
        .integer => .integer,
        .semicont => .semi_continuous,
        .semiint => .semi_integer,
    };
    const row_lower = self.allocator.alloc(f64, self.num_constrs) catch return error.OutOfMemory;
    defer self.allocator.free(row_lower);
    const row_upper = self.allocator.alloc(f64, self.num_constrs) catch return error.OutOfMemory;
    defer self.allocator.free(row_upper);
    for (self.constr_sense, self.constr_rhs, row_lower, row_upper) |sense, rhs, *lower, *upper| switch (sense) {
        .less_equal => {
            lower.* = -std.math.inf(f64);
            upper.* = rhs;
        },
        .equal => {
            lower.* = rhs;
            upper.* = rhs;
        },
        .greater_equal => {
            lower.* = rhs;
            upper.* = std.math.inf(f64);
        },
    };
    const view = io.ModelView{
        .name = self.name,
        .objective_sense = if (self.sense == .minimize) .minimize else .maximize,
        .objective_offset = self.obj_con,
        .col_cost = self.var_obj,
        .col_lower = self.var_lb,
        .col_upper = self.var_ub,
        .col_type = col_types,
        .col_names = self.var_names,
        .row_lower = row_lower,
        .row_upper = row_upper,
        .row_names = self.constr_names,
        .matrix = self.matrix.csc().view(),
    };
    const io_context = std.Io.Threaded.global_single_threaded.io();
    io.writeFile(io_context, self.allocator, filename, view, .{}) catch |err| return mapIoError(err);
}

pub const write = writeModel;

pub fn readModel(self: *Model, filename: []const u8) ModelError!void {
    try self.updateModel();
    const io_context = std.Io.Threaded.global_single_threaded.io();
    var imported = io.readFile(io_context, self.allocator, filename, .{}) catch |err| return mapIoError(err);
    defer imported.deinit();
    var replacement = Model.init(self.allocator, self.env, imported.name) catch |err| return err;
    errdefer replacement.deinit();
    replacement.sense = if (imported.objective_sense == .minimize) .minimize else .maximize;
    replacement.obj_con = imported.objective_offset;

    const row_first = self.allocator.alloc(usize, imported.row_lower.len) catch return error.OutOfMemory;
    defer self.allocator.free(row_first);
    const row_second = self.allocator.alloc(?usize, imported.row_lower.len) catch return error.OutOfMemory;
    defer self.allocator.free(row_second);
    @memset(row_second, null);
    for (0..imported.row_lower.len) |row| {
        row_first[row] = replacement.numConstrs();
        const lower = imported.row_lower[row];
        const upper = imported.row_upper[row];
        const name = if (imported.row_names[row]) |value| value else null;
        if (lower == upper) {
            try replacement.addConstr(0, &.{}, &.{}, .equal, lower, name);
        } else if (std.math.isInf(lower)) {
            try replacement.addConstr(0, &.{}, &.{}, .less_equal, upper, name);
        } else if (std.math.isInf(upper)) {
            try replacement.addConstr(0, &.{}, &.{}, .greater_equal, lower, name);
        } else {
            try replacement.addConstr(0, &.{}, &.{}, .greater_equal, lower, name);
            row_second[row] = replacement.numConstrs();
            const upper_name = if (name) |value| std.fmt.allocPrint(self.allocator, "{s}_range_upper", .{value}) catch return error.OutOfMemory else null;
            defer if (upper_name) |value| self.allocator.free(value);
            try replacement.addConstr(0, &.{}, &.{}, .less_equal, upper, upper_name);
        }
    }
    const scratch_capacity = std.math.mul(usize, maxColumnLength(imported.matrix), 2) catch return error.InvalidArgument;
    const scratch_rows = self.allocator.alloc(usize, scratch_capacity) catch return error.OutOfMemory;
    defer self.allocator.free(scratch_rows);
    const scratch_values = self.allocator.alloc(f64, scratch_capacity) catch return error.OutOfMemory;
    defer self.allocator.free(scratch_values);
    for (0..imported.col_cost.len) |column| {
        const begin = imported.matrix.col_starts[column];
        const end = imported.matrix.col_starts[column + 1];
        var output: usize = 0;
        for (imported.matrix.row_indices[begin..end], imported.matrix.values[begin..end]) |row_id, value| {
            const row = row_id.toUsize();
            scratch_rows[output] = row_first[row];
            scratch_values[output] = value;
            output += 1;
            if (row_second[row]) |second| {
                scratch_rows[output] = second;
                scratch_values[output] = value;
                output += 1;
            }
        }
        try replacement.addVar(
            output,
            scratch_rows[0..output],
            scratch_values[0..output],
            imported.col_cost[column],
            imported.col_lower[column],
            imported.col_upper[column],
            switch (imported.col_type[column]) {
                .continuous => .continuous,
                .binary => .binary,
                .integer => .integer,
                .semi_continuous => .semicont,
                .semi_integer => .semiint,
            },
            if (imported.col_names[column]) |name| name else null,
        );
    }
    try replacement.updateModel();
    self.deinit();
    self.* = replacement;
}

pub const read = readModel;

fn maxColumnLength(matrix: @import("matrix").CscMatrix) usize {
    var maximum: usize = 0;
    for (0..matrix.num_cols) |column| maximum = @max(maximum, matrix.col_starts[column + 1] - matrix.col_starts[column]);
    return maximum;
}

fn mapIoError(err: io.IoError) ModelError {
    return switch (err) {
        error.OutOfMemory => error.OutOfMemory,
        error.FileNotFound, error.PermissionDenied, error.ReadFailed, error.WriteFailed, error.FileTooLarge => error.IoError,
        error.UnsupportedFormat, error.UnsupportedCompression, error.UnsupportedFeature => error.FeatureNotAvailable,
        error.DuplicateName => error.DuplicateName,
        else => error.InvalidArgument,
    };
}

test "model I/O adapter stays a thin type and ownership boundary" {
    try std.testing.expect(@sizeOf(io.ModelView) > 0);
}

test "model adapter atomically imports MPS and exports LP" {
    const Env = @import("env.zig").Env;
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "old");
    defer model.deinit();
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const source =
        \\NAME ADAPTER
        \\ROWS
        \\ N OBJ
        \\ L CAP
        \\COLUMNS
        \\ X OBJ -2 CAP 3
        \\RHS
        \\ RHS1 CAP 7
        \\BOUNDS
        \\ UP BND X 4
        \\ENDATA
    ;
    try tmp.dir.writeFile(std.testing.io, "adapter.mps", source);
    const input_path = try std.fmt.allocPrint(std.testing.allocator, ".zig-cache/tmp/{s}/adapter.mps", .{tmp.sub_path});
    defer std.testing.allocator.free(input_path);
    const output_path = try std.fmt.allocPrint(std.testing.allocator, ".zig-cache/tmp/{s}/adapter.lp", .{tmp.sub_path});
    defer std.testing.allocator.free(output_path);
    try model.read(input_path);
    try std.testing.expectEqualStrings("ADAPTER", model.name);
    try std.testing.expectEqual(@as(usize, 1), model.num_vars);
    try std.testing.expectEqual(@as(usize, 1), model.num_constrs);
    try std.testing.expectEqual(@as(f64, 3.0), try model.getCoeff(0, 0));
    try model.write(output_path);
}
