//! Shared semantic builder used by text frontends.
//!
//! LP benefits from its duplicate-merging triplet path. MPS may switch to a
//! direct CSC builder while keeping the same final `ModelData` contract.

const std = @import("std");
const foundation = @import("foundation");
const matrix = @import("matrix");
const types = @import("types.zig");
const ModelData = @import("model_data.zig").ModelData;

pub const Column = struct {
    name: []const u8,
    cost: f64 = 0.0,
    lower: f64 = 0.0,
    upper: f64 = std.math.inf(f64),
    kind: types.VariableType = .continuous,
};

pub const Row = struct {
    name: ?[]const u8,
    lower: f64 = -std.math.inf(f64),
    upper: f64 = std.math.inf(f64),
};

pub const Builder = struct {
    allocator: std.mem.Allocator,
    name: []const u8 = "",
    objective_sense: types.ObjectiveSense = .minimize,
    objective_offset: f64 = 0.0,
    columns: std.ArrayListUnmanaged(Column) = .empty,
    rows: std.ArrayListUnmanaged(Row) = .empty,
    terms: std.ArrayListUnmanaged(Term) = .empty,

    const Term = struct { row: usize, col: usize, value: f64 };

    pub fn init(allocator: std.mem.Allocator) Builder {
        return .{ .allocator = allocator };
    }

    pub fn deinit(self: *Builder) void {
        self.columns.deinit(self.allocator);
        self.rows.deinit(self.allocator);
        self.terms.deinit(self.allocator);
        self.* = undefined;
    }

    pub fn addColumn(self: *Builder, column: Column) types.IoError!usize {
        if (column.lower > column.upper) return error.InvalidBounds;
        const index = self.columns.items.len;
        self.columns.append(self.allocator, column) catch return error.OutOfMemory;
        return index;
    }

    pub fn addRow(self: *Builder, row: Row) types.IoError!usize {
        const index = self.rows.items.len;
        self.rows.append(self.allocator, row) catch return error.OutOfMemory;
        return index;
    }

    pub fn addTerm(self: *Builder, row: usize, col: usize, value: f64) types.IoError!void {
        if (row >= self.rows.items.len or col >= self.columns.items.len) return error.InvalidDimensions;
        if (!std.math.isFinite(value)) return error.NonFiniteValue;
        if (value == 0.0) return;
        self.terms.append(self.allocator, .{ .row = row, .col = col, .value = value }) catch return error.OutOfMemory;
    }

    pub fn finish(self: *Builder, options: types.ReadOptions) types.IoError!ModelData {
        const allocator = self.allocator;
        var matrix_builder = matrix.MatrixBuilder.init(self.rows.items.len, self.columns.items.len) catch return error.InvalidDimensions;
        defer matrix_builder.deinit(allocator);
        matrix_builder.reserve(allocator, self.terms.items.len) catch return error.OutOfMemory;
        for (self.terms.items) |term| {
            const row = foundation.RowId.fromUsize(term.row) catch return error.InvalidDimensions;
            const col = foundation.ColId.fromUsize(term.col) catch return error.InvalidDimensions;
            matrix_builder.appendPreReserved(row, col, term.value);
        }
        const csc = matrix_builder.freeze(allocator, options.zero_tolerance) catch |err| return switch (err) {
            error.OutOfMemory => error.OutOfMemory,
            error.NonFiniteValue => error.NonFiniteValue,
            else => error.InvalidDimensions,
        };
        return self.finishWithMatrix(options, csc);
    }

    /// Fast path for MPS COLUMNS records. When coordinates arrive ordered by
    /// `(column,row)`, canonical CSC is created in one exact-size allocation
    /// without copying triplets into a second builder or running a global sort.
    /// Unordered input transparently falls back to the general deterministic
    /// path, preserving format compatibility.
    pub fn finishColumnOrdered(self: *Builder, options: types.ReadOptions) types.IoError!ModelData {
        var previous_col: usize = 0;
        var previous_row: usize = 0;
        var first = true;
        for (self.terms.items) |term| {
            if (!first and (term.col < previous_col or (term.col == previous_col and term.row < previous_row))) return self.finish(options);
            previous_col = term.col;
            previous_row = term.row;
            first = false;
        }
        var nonzeros: usize = 0;
        var index: usize = 0;
        while (index < self.terms.items.len) {
            const coordinate = self.terms.items[index];
            var sum: f64 = 0.0;
            while (index < self.terms.items.len and self.terms.items[index].col == coordinate.col and self.terms.items[index].row == coordinate.row) : (index += 1) sum += self.terms.items[index].value;
            if (!std.math.isFinite(sum)) return error.NonFiniteValue;
            if (@abs(sum) > options.zero_tolerance) nonzeros += 1;
        }
        var csc = matrix.CscMatrix.initPackedUninitialized(self.allocator, self.rows.items.len, self.columns.items.len, nonzeros) catch |err| return switch (err) {
            error.OutOfMemory => error.OutOfMemory,
            else => error.InvalidDimensions,
        };
        @memset(csc.col_starts, 0);
        index = 0;
        var output_index: usize = 0;
        while (index < self.terms.items.len) {
            const coordinate = self.terms.items[index];
            var sum: f64 = 0.0;
            while (index < self.terms.items.len and self.terms.items[index].col == coordinate.col and self.terms.items[index].row == coordinate.row) : (index += 1) sum += self.terms.items[index].value;
            if (@abs(sum) <= options.zero_tolerance) continue;
            csc.row_indices[output_index] = foundation.RowId.fromUsize(coordinate.row) catch {
                csc.deinit(self.allocator);
                return error.InvalidDimensions;
            };
            csc.values[output_index] = sum;
            csc.col_starts[coordinate.col + 1] += 1;
            output_index += 1;
        }
        for (0..self.columns.items.len) |column| csc.col_starts[column + 1] += csc.col_starts[column];
        return self.finishWithMatrix(options, csc);
    }

    fn finishWithMatrix(self: *Builder, options: types.ReadOptions, csc: matrix.CscMatrix) types.IoError!ModelData {
        const allocator = self.allocator;
        var owned_csc = csc;
        errdefer owned_csc.deinit(allocator);

        const name = allocator.dupe(u8, self.name) catch return error.OutOfMemory;
        errdefer allocator.free(name);
        const col_cost = allocator.alloc(f64, self.columns.items.len) catch return error.OutOfMemory;
        errdefer allocator.free(col_cost);
        const col_lower = allocator.alloc(f64, self.columns.items.len) catch return error.OutOfMemory;
        errdefer allocator.free(col_lower);
        const col_upper = allocator.alloc(f64, self.columns.items.len) catch return error.OutOfMemory;
        errdefer allocator.free(col_upper);
        const col_type = allocator.alloc(types.VariableType, self.columns.items.len) catch return error.OutOfMemory;
        errdefer allocator.free(col_type);
        const col_names = allocator.alloc(?[]u8, self.columns.items.len) catch return error.OutOfMemory;
        errdefer allocator.free(col_names);
        @memset(col_names, null);
        var col_names_init: usize = 0;
        errdefer for (col_names[0..col_names_init]) |maybe| if (maybe) |value| allocator.free(value);
        for (self.columns.items, 0..) |column, i| {
            col_cost[i] = column.cost;
            col_lower[i] = column.lower;
            col_upper[i] = column.upper;
            col_type[i] = column.kind;
            if (options.keep_names) col_names[i] = allocator.dupe(u8, column.name) catch return error.OutOfMemory;
            col_names_init += 1;
        }

        const row_lower = allocator.alloc(f64, self.rows.items.len) catch return error.OutOfMemory;
        errdefer allocator.free(row_lower);
        const row_upper = allocator.alloc(f64, self.rows.items.len) catch return error.OutOfMemory;
        errdefer allocator.free(row_upper);
        const row_names = allocator.alloc(?[]u8, self.rows.items.len) catch return error.OutOfMemory;
        errdefer allocator.free(row_names);
        @memset(row_names, null);
        var row_names_init: usize = 0;
        errdefer for (row_names[0..row_names_init]) |maybe| if (maybe) |value| allocator.free(value);
        for (self.rows.items, 0..) |row, i| {
            row_lower[i] = row.lower;
            row_upper[i] = row.upper;
            if (options.keep_names and row.name != null) row_names[i] = allocator.dupe(u8, row.name.?) catch return error.OutOfMemory;
            row_names_init += 1;
        }

        return .{
            .allocator = allocator,
            .name = name,
            .objective_sense = self.objective_sense,
            .objective_offset = self.objective_offset,
            .col_cost = col_cost,
            .col_lower = col_lower,
            .col_upper = col_upper,
            .col_type = col_type,
            .col_names = col_names,
            .row_lower = row_lower,
            .row_upper = row_upper,
            .row_names = row_names,
            .matrix = owned_csc,
        };
    }
};
