//! Shared semantic builder used by text frontends.
//!
//! General terms are canonicalized in place and frozen directly into final
//! CSC storage. LP additionally uses column chains: row-ordered terms enter a
//! shared node pool, duplicates merge at append time, and columns freeze to CSC
//! without a global coordinate sort or per-column heap allocations.

const std = @import("std");
const foundation = @import("foundation");
const matrix = @import("matrix");
const types = @import("types.zig");
const ModelData = @import("model_data.zig").ModelData;
const StringArena = @import("string_arena.zig").StringArena;

pub const Column = struct {
    /// Borrowed unique symbol used while resolving parser references.
    name: []const u8,
    /// Linear objective coefficient.
    cost: f64 = 0.0,
    /// Inclusive variable lower bound.
    lower: f64 = 0.0,
    /// Inclusive variable upper bound.
    upper: f64 = std.math.inf(f64),
    /// Variable domain/integrality classification.
    kind: types.VariableType = .continuous,
};

pub const Row = struct {
    /// Optional borrowed row symbol.
    name: ?[]const u8,
    /// Inclusive row-activity lower bound.
    lower: f64 = -std.math.inf(f64),
    /// Inclusive row-activity upper bound.
    upper: f64 = std.math.inf(f64),
};

pub const Builder = struct {
    /// Owner of all temporary semantic arrays.
    allocator: std.mem.Allocator,
    /// Borrowed model name copied only during final publication.
    name: []const u8 = "",
    /// Objective direction parsed from the source.
    objective_sense: types.ObjectiveSense = .minimize,
    /// Constant objective term.
    objective_offset: f64 = 0.0,
    /// Temporary structural-column metadata.
    columns: std.ArrayListUnmanaged(Column) = .empty,
    /// Temporary row metadata.
    rows: std.ArrayListUnmanaged(Row) = .empty,
    /// General unordered coordinate stream used by MPS and generic callers.
    terms: std.ArrayListUnmanaged(Term) = .empty,
    /// LP-specific linked-node pool grouped by structural column.
    column_terms: std.ArrayListUnmanaged(ColumnTerm) = .empty,
    /// Head/tail descriptor for each LP column chain.
    column_chains: std.ArrayListUnmanaged(ColumnChain) = .empty,
    /// Whether terms must use the LP column-chain fast path.
    column_term_mode: bool = false,
    /// Configured semantic row limit.
    max_rows: usize = std.math.maxInt(usize),
    /// Configured semantic column limit.
    max_columns: usize = std.math.maxInt(usize),
    /// Configured temporary coordinate/node limit.
    max_matrix_terms: usize = std.math.maxInt(usize),
    /// Configured sum limit for retained unique name bytes.
    max_name_bytes: usize = std.math.maxInt(usize),
    /// Name bytes accepted so far.
    name_bytes: usize = 0,
    /// Cooperative cancellation poller for long finalization loops.
    control: types.ParseControl = types.ParseControl.init(.{}),

    /// Strong IDs make invalid dimensions fail at append time and keep each
    /// term at 16 bytes with the default 32-bit index configuration.
    const Term = struct {
        /// Temporary row coordinate.
        row: foundation.RowId,
        /// Temporary structural-column coordinate.
        col: foundation.ColId,
        /// Finite nonzero coefficient.
        value: f64,
    };
    const no_term = std.math.maxInt(foundation.HUInt);
    const ColumnTerm = struct {
        /// Row coordinate within this column's ordered chain.
        row: foundation.RowId,
        /// Next node index, or `no_term` at chain end.
        next: foundation.HUInt = no_term,
        /// Merged finite coefficient for this coordinate.
        value: f64,
    };
    const ColumnChain = struct {
        /// First node index, or `no_term` for an empty column.
        head: foundation.HUInt = no_term,
        /// Last node index, used for O(1) append and duplicate merging.
        tail: foundation.HUInt = no_term,
    };

    /// Construct an empty semantic builder with unrestricted default limits.
    pub fn init(allocator: std.mem.Allocator) Builder {
        return .{ .allocator = allocator };
    }

    /// Copy resource limits and cancellation policy from public read options.
    pub fn configureLimits(self: *Builder, options: types.ReadOptions) void {
        self.max_rows = options.max_rows;
        self.max_columns = options.max_columns;
        self.max_matrix_terms = options.max_matrix_terms;
        self.max_name_bytes = options.max_name_bytes;
        self.control = types.ParseControl.init(options);
    }

    /// Release every temporary semantic and term array.
    pub fn deinit(self: *Builder) void {
        self.columns.deinit(self.allocator);
        self.rows.deinit(self.allocator);
        self.terms.deinit(self.allocator);
        self.column_terms.deinit(self.allocator);
        self.column_chains.deinit(self.allocator);
        self.* = undefined;
    }

    /// Select the LP fast path before adding columns or generic terms.
    pub fn enableColumnTermStorage(self: *Builder) void {
        std.debug.assert(self.columns.items.len == 0 and self.terms.items.len == 0);
        self.column_term_mode = true;
    }

    /// Validate and append column metadata, returning its stable temporary index.
    pub fn addColumn(self: *Builder, column: Column) types.IoError!usize {
        if (column.lower > column.upper) return error.InvalidBounds;
        if (self.columns.items.len >= self.max_columns) return error.ResourceLimitExceeded;
        const next_name_bytes = std.math.add(usize, self.name_bytes, column.name.len) catch return error.ResourceLimitExceeded;
        if (next_name_bytes > self.max_name_bytes) return error.ResourceLimitExceeded;
        const index = self.columns.items.len;
        self.columns.ensureUnusedCapacity(self.allocator, 1) catch return error.OutOfMemory;
        if (self.column_term_mode) self.column_chains.ensureUnusedCapacity(self.allocator, 1) catch return error.OutOfMemory;
        self.columns.appendAssumeCapacity(column);
        if (self.column_term_mode) self.column_chains.appendAssumeCapacity(.{});
        self.name_bytes = next_name_bytes;
        return index;
    }

    /// Append row metadata after enforcing row and name-byte limits.
    pub fn addRow(self: *Builder, row: Row) types.IoError!usize {
        if (self.rows.items.len >= self.max_rows) return error.ResourceLimitExceeded;
        const name_len = if (row.name) |name| name.len else 0;
        const next_name_bytes = std.math.add(usize, self.name_bytes, name_len) catch return error.ResourceLimitExceeded;
        if (next_name_bytes > self.max_name_bytes) return error.ResourceLimitExceeded;
        const index = self.rows.items.len;
        self.rows.append(self.allocator, row) catch return error.OutOfMemory;
        self.name_bytes = next_name_bytes;
        return index;
    }

    /// Append one checked generic matrix coordinate.
    pub fn addTerm(self: *Builder, row: usize, col: usize, value: f64) types.IoError!void {
        if (self.column_term_mode) return error.InvalidDimensions;
        if (row >= self.rows.items.len or col >= self.columns.items.len) return error.InvalidDimensions;
        if (!std.math.isFinite(value)) return error.NonFiniteValue;
        if (value == 0.0) return;
        if (self.terms.items.len >= self.max_matrix_terms) return error.ResourceLimitExceeded;
        const row_id = foundation.RowId.fromUsize(row) catch return error.InvalidDimensions;
        const col_id = foundation.ColId.fromUsize(col) catch return error.InvalidDimensions;
        self.terms.append(self.allocator, .{ .row = row_id, .col = col_id, .value = value }) catch return error.OutOfMemory;
    }

    /// Append an LP term directly to its column chain. Rows for each column
    /// must be nondecreasing, which follows naturally from row-wise LP input.
    /// A repeated `(column,row)` is merged into the tail node immediately.
    pub fn addColumnTerm(self: *Builder, row: usize, col: usize, value: f64) types.IoError!void {
        if (!self.column_term_mode or row >= self.rows.items.len or col >= self.columns.items.len) return error.InvalidDimensions;
        if (!std.math.isFinite(value)) return error.NonFiniteValue;
        if (value == 0.0) return;
        const row_id = foundation.RowId.fromUsize(row) catch return error.InvalidDimensions;
        const chain = &self.column_chains.items[col];
        if (chain.tail != no_term) {
            const tail = &self.column_terms.items[@intCast(chain.tail)];
            const tail_row = tail.row.toUsize();
            if (tail_row > row) return error.InvalidDimensions;
            if (tail_row == row) {
                const sum = tail.value + value;
                if (!std.math.isFinite(sum)) return error.NonFiniteValue;
                tail.value = sum;
                return;
            }
        }

        if (self.column_terms.items.len >= self.max_matrix_terms) return error.ResourceLimitExceeded;
        if (self.column_terms.items.len >= no_term) return error.InvalidDimensions;
        const node_index: foundation.HUInt = @intCast(self.column_terms.items.len);
        self.column_terms.append(self.allocator, .{ .row = row_id, .value = value }) catch return error.OutOfMemory;
        if (chain.tail == no_term)
            chain.head = node_index
        else
            self.column_terms.items[@intCast(chain.tail)].next = node_index;
        chain.tail = node_index;
    }

    /// Canonicalize the active term representation and publish owned model data.
    pub fn finish(self: *Builder, options: types.ReadOptions) types.IoError!ModelData {
        if (self.column_term_mode) return self.finishColumnTerms(options);
        return self.finishTerms(options, false);
    }

    /// Freeze LP column chains directly into canonical CSC. The first pass
    /// counts surviving values for an exact allocation; the second writes each
    /// already row-ordered chain without sorting or cursor scratch storage.
    pub fn finishColumnTerms(self: *Builder, options: types.ReadOptions) types.IoError!ModelData {
        if (!self.column_term_mode or self.terms.items.len != 0) return error.InvalidDimensions;
        if (!std.math.isFinite(options.zero_tolerance) or options.zero_tolerance < 0.0) return error.InvalidTolerance;

        var nonzeros: usize = 0;
        for (self.column_chains.items) |chain| {
            var node_index = chain.head;
            while (node_index != no_term) {
                try self.control.tick();
                const node = self.column_terms.items[@intCast(node_index)];
                if (!std.math.isFinite(node.value)) return error.NonFiniteValue;
                if (@abs(node.value) > options.zero_tolerance) nonzeros = std.math.add(usize, nonzeros, 1) catch return error.InvalidDimensions;
                node_index = node.next;
            }
        }

        var csc = matrix.CscMatrix.initPackedUninitialized(self.allocator, self.rows.items.len, self.columns.items.len, nonzeros) catch |err| return switch (err) {
            error.OutOfMemory => error.OutOfMemory,
            else => error.InvalidDimensions,
        };
        var owns_csc = true;
        errdefer if (owns_csc) csc.deinit(self.allocator);
        var output_index: usize = 0;
        for (self.column_chains.items, 0..) |chain, column| {
            csc.col_starts[column] = output_index;
            var node_index = chain.head;
            while (node_index != no_term) {
                try self.control.tick();
                const node = self.column_terms.items[@intCast(node_index)];
                if (@abs(node.value) > options.zero_tolerance) {
                    csc.row_indices[output_index] = node.row;
                    csc.values[output_index] = node.value;
                    output_index += 1;
                }
                node_index = node.next;
            }
        }
        csc.col_starts[self.columns.items.len] = output_index;
        std.debug.assert(output_index == nonzeros);

        self.column_terms.deinit(self.allocator);
        self.column_terms = .empty;
        self.column_chains.deinit(self.allocator);
        self.column_chains = .empty;
        owns_csc = false;
        return self.finishWithMatrix(options, csc);
    }

    /// Fast path for MPS COLUMNS records. When coordinates arrive ordered by
    /// `(column,row)`, canonical CSC is created in one exact-size allocation
    /// without copying triplets into a second builder or running a global sort.
    /// Unordered input transparently falls back to the general deterministic
    /// path, preserving format compatibility.
    pub fn finishColumnOrdered(self: *Builder, options: types.ReadOptions) types.IoError!ModelData {
        if (self.column_term_mode) return self.finishColumnTerms(options);
        if (!try self.termsAreOrdered()) return self.finish(options);
        return self.finishTerms(options, true);
    }

    /// Freeze the generic triplet stream, optionally trusting column/row order.
    fn finishTerms(self: *Builder, options: types.ReadOptions, already_ordered: bool) types.IoError!ModelData {
        if (!std.math.isFinite(options.zero_tolerance) or options.zero_tolerance < 0.0) return error.InvalidTolerance;
        try self.control.checkNow();
        if (!already_ordered) {
            std.sort.block(Term, self.terms.items, {}, termLessThan);
            try self.control.checkNow();
        }
        try self.mergeTermsInPlace(options.zero_tolerance);

        var csc = matrix.CscMatrix.initPackedUninitialized(self.allocator, self.rows.items.len, self.columns.items.len, self.terms.items.len) catch |err| return switch (err) {
            error.OutOfMemory => error.OutOfMemory,
            else => error.InvalidDimensions,
        };
        var owns_csc = true;
        errdefer if (owns_csc) csc.deinit(self.allocator);
        @memset(csc.col_starts, 0);
        for (self.terms.items, 0..) |term, output_index| {
            try self.control.tick();
            csc.row_indices[output_index] = term.row;
            csc.values[output_index] = term.value;
            csc.col_starts[term.col.toUsize() + 1] += 1;
        }
        for (0..self.columns.items.len) |column| {
            try self.control.tick();
            csc.col_starts[column + 1] += csc.col_starts[column];
        }

        // The final CSC no longer borrows terms. Release the largest semantic
        // builder buffer before allocating published names and bound arrays.
        self.terms.deinit(self.allocator);
        self.terms = .empty;
        owns_csc = false;
        return self.finishWithMatrix(options, csc);
    }

    /// Merge adjacent duplicate coordinates and discard final numerical zeros.
    fn mergeTermsInPlace(self: *Builder, zero_tolerance: f64) types.IoError!void {
        var read: usize = 0;
        var write: usize = 0;
        while (read < self.terms.items.len) {
            try self.control.tick();
            const coordinate = self.terms.items[read];
            var sum = coordinate.value;
            read += 1;
            while (read < self.terms.items.len and sameCoordinate(self.terms.items[read], coordinate)) : (read += 1) {
                try self.control.tick();
                sum += self.terms.items[read].value;
            }
            if (!std.math.isFinite(sum)) return error.NonFiniteValue;
            if (@abs(sum) <= zero_tolerance) continue;
            self.terms.items[write] = .{ .row = coordinate.row, .col = coordinate.col, .value = sum };
            write += 1;
        }
        self.terms.shrinkRetainingCapacity(write);
    }

    /// Check nondecreasing column/row order while polling cancellation.
    fn termsAreOrdered(self: *Builder) types.IoError!bool {
        if (self.terms.items.len < 2) return true;
        for (self.terms.items[1..], self.terms.items[0 .. self.terms.items.len - 1]) |current, previous| {
            try self.control.tick();
            const current_col = current.col.toUsize();
            const previous_col = previous.col.toUsize();
            if (current_col < previous_col) return false;
            if (current_col == previous_col and current.row.toUsize() < previous.row.toUsize()) return false;
        }
        return true;
    }

    /// Transfer canonical CSC ownership and copy semantic/name arrays to ModelData.
    fn finishWithMatrix(self: *Builder, options: types.ReadOptions, csc: matrix.CscMatrix) types.IoError!ModelData {
        var owned_csc = csc;
        var owns_csc = true;
        errdefer if (owns_csc) owned_csc.deinit(self.allocator);
        var additional_name_bytes: usize = 0;
        if (options.keep_names) {
            for (self.columns.items) |column| {
                try self.control.tick();
                additional_name_bytes = std.math.add(usize, additional_name_bytes, column.name.len) catch return error.InvalidDimensions;
            }
            for (self.rows.items) |row| if (row.name) |name| {
                try self.control.tick();
                additional_name_bytes = std.math.add(usize, additional_name_bytes, name.len) catch return error.InvalidDimensions;
            };
        }

        owns_csc = false;
        var result = try ModelData.initPacked(
            self.allocator,
            owned_csc,
            self.name,
            additional_name_bytes,
            options.keep_names,
            self.objective_sense,
            self.objective_offset,
        );
        errdefer result.deinit();
        var name_arena = StringArena.fromUsed(result.name_storage, result.name.len) catch return error.InvalidDimensions;
        for (self.columns.items, 0..) |column, index| {
            try self.control.tick();
            result.col_cost[index] = column.cost;
            result.col_lower[index] = column.lower;
            result.col_upper[index] = column.upper;
            result.col_type[index] = column.kind;
            if (options.keep_names) result.col_names[index] = name_arena.append(column.name) catch return error.InvalidDimensions;
        }
        for (self.rows.items, 0..) |row, index| {
            try self.control.tick();
            result.row_lower[index] = row.lower;
            result.row_upper[index] = row.upper;
            if (options.keep_names) {
                if (row.name) |name| result.row_names[index] = name_arena.append(name) catch return error.InvalidDimensions;
            }
        }
        std.debug.assert(name_arena.used() == self.name.len + additional_name_bytes);

        // Attribute values and retained names are now owned by `result`.
        // Release semantic metadata immediately instead of waiting for the
        // format parser's deferred Builder.deinit after parse returns.
        self.columns.deinit(self.allocator);
        self.columns = .empty;
        self.rows.deinit(self.allocator);
        self.rows = .empty;
        return result;
    }
};

/// Order temporary coordinates by column then row.
fn termLessThan(_: void, lhs: Builder.Term, rhs: Builder.Term) bool {
    const lhs_col = lhs.col.toUsize();
    const rhs_col = rhs.col.toUsize();
    if (lhs_col != rhs_col) return lhs_col < rhs_col;
    return lhs.row.toUsize() < rhs.row.toUsize();
}

/// Return whether two temporary terms address the same matrix coordinate.
fn sameCoordinate(lhs: Builder.Term, rhs: Builder.Term) bool {
    return lhs.col == rhs.col and lhs.row == rhs.row;
}

test "finish sorts terms in place merges duplicates and applies tolerance" {
    var builder = Builder.init(std.testing.allocator);
    defer builder.deinit();
    _ = try builder.addColumn(.{ .name = "x" });
    _ = try builder.addColumn(.{ .name = "y" });
    _ = try builder.addColumn(.{ .name = "z" });
    _ = try builder.addRow(.{ .name = "r0" });
    _ = try builder.addRow(.{ .name = "r1" });
    try builder.addTerm(1, 2, 4.0);
    try builder.addTerm(0, 0, 2.0);
    try builder.addTerm(1, 2, -1.0);
    try builder.addTerm(1, 1, 1e-10);

    var model = try builder.finish(.{ .zero_tolerance = 1e-9 });
    defer model.deinit();
    try std.testing.expectEqual(@as(usize, 0), builder.terms.capacity);
    try model.matrix.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 1, 1, 2 }, model.matrix.col_starts);
    try std.testing.expectEqualSlices(f64, &.{ 2.0, 3.0 }, model.matrix.values);
    try std.testing.expectEqual(@as(usize, 0), model.matrix.row_indices[0].toUsize());
    try std.testing.expectEqual(@as(usize, 1), model.matrix.row_indices[1].toUsize());
}

test "finish preserves insertion order while summing duplicate terms" {
    var builder = Builder.init(std.testing.allocator);
    defer builder.deinit();
    _ = try builder.addColumn(.{ .name = "x" });
    _ = try builder.addRow(.{ .name = "r" });
    try builder.addTerm(0, 0, 1e16);
    try builder.addTerm(0, 0, 1.0);
    try builder.addTerm(0, 0, -1e16);

    var model = try builder.finish(.{});
    defer model.deinit();
    try std.testing.expectEqual(@as(usize, 0), model.matrix.nnz());
}

test "ordered MPS freeze releases semantic terms before model publication" {
    var builder = Builder.init(std.testing.allocator);
    defer builder.deinit();
    _ = try builder.addColumn(.{ .name = "x" });
    _ = try builder.addColumn(.{ .name = "y" });
    _ = try builder.addRow(.{ .name = "r0" });
    _ = try builder.addRow(.{ .name = "r1" });
    try builder.addTerm(0, 0, 2.0);
    try builder.addTerm(0, 0, 3.0);
    try builder.addTerm(1, 0, 4.0);
    try builder.addTerm(0, 1, -1.0);

    var model = try builder.finishColumnOrdered(.{});
    defer model.deinit();
    try model.matrix.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 2, 3 }, model.matrix.col_starts);
    try std.testing.expectEqualSlices(f64, &.{ 5.0, 4.0, -1.0 }, model.matrix.values);
    try std.testing.expectEqual(@as(usize, 0), builder.terms.items.len);
    try std.testing.expectEqual(@as(usize, 0), builder.terms.capacity);
}

test "finish rejects invalid zero tolerance" {
    var builder = Builder.init(std.testing.allocator);
    defer builder.deinit();
    try std.testing.expectError(error.InvalidTolerance, builder.finish(.{ .zero_tolerance = -1.0 }));
    try std.testing.expectError(error.InvalidTolerance, builder.finish(.{ .zero_tolerance = std.math.nan(f64) }));
}

test "term storage uses compact strong IDs" {
    if (@sizeOf(foundation.HUInt) == 4) try std.testing.expectEqual(@as(usize, 16), @sizeOf(Builder.Term));
}

test "builder enforces semantic limits before growing storage" {
    var columns = Builder.init(std.testing.allocator);
    defer columns.deinit();
    columns.configureLimits(.{ .max_columns = 0 });
    try std.testing.expectError(error.ResourceLimitExceeded, columns.addColumn(.{ .name = "x" }));
    try std.testing.expectEqual(@as(usize, 0), columns.columns.capacity);

    var rows = Builder.init(std.testing.allocator);
    defer rows.deinit();
    rows.configureLimits(.{ .max_rows = 0 });
    try std.testing.expectError(error.ResourceLimitExceeded, rows.addRow(.{ .name = "r" }));
    try std.testing.expectEqual(@as(usize, 0), rows.rows.capacity);

    var names = Builder.init(std.testing.allocator);
    defer names.deinit();
    names.configureLimits(.{ .max_name_bytes = 1 });
    try std.testing.expectError(error.ResourceLimitExceeded, names.addColumn(.{ .name = "long" }));
    try std.testing.expectEqual(@as(usize, 0), names.name_bytes);

    var terms = Builder.init(std.testing.allocator);
    defer terms.deinit();
    terms.configureLimits(.{ .max_matrix_terms = 0 });
    _ = try terms.addColumn(.{ .name = "x" });
    _ = try terms.addRow(.{ .name = "r" });
    try std.testing.expectError(error.ResourceLimitExceeded, terms.addTerm(0, 0, 1.0));
    try std.testing.expectEqual(@as(usize, 0), terms.terms.capacity);
}

test "builder finalization cooperatively cancels without leaking" {
    var interrupted = std.atomic.Value(bool).init(false);
    var builder = Builder.init(std.testing.allocator);
    defer builder.deinit();
    builder.configureLimits(.{ .interrupt_flag = &interrupted, .interrupt_check_interval = 1 });
    _ = try builder.addColumn(.{ .name = "x" });
    _ = try builder.addRow(.{ .name = "r" });
    try builder.addTerm(0, 0, 1.0);
    interrupted.store(true, .release);
    try std.testing.expectError(error.Cancelled, builder.finish(.{}));
}

test "LP column terms merge at append and freeze directly to CSC" {
    var builder = Builder.init(std.testing.allocator);
    defer builder.deinit();
    builder.enableColumnTermStorage();
    _ = try builder.addColumn(.{ .name = "x" });
    _ = try builder.addColumn(.{ .name = "y" });
    _ = try builder.addRow(.{ .name = "r0" });
    _ = try builder.addRow(.{ .name = "r1" });
    _ = try builder.addRow(.{ .name = "r2" });

    try builder.addColumnTerm(0, 0, 3.0);
    try builder.addColumnTerm(0, 1, 4.0);
    try builder.addColumnTerm(0, 0, -3.0);
    try builder.addColumnTerm(1, 1, 2.0);
    try builder.addColumnTerm(1, 0, 1.0);
    try builder.addColumnTerm(1, 0, 2.0);
    try builder.addColumnTerm(2, 1, 1e-10);
    try std.testing.expectEqual(@as(usize, 0), builder.terms.items.len);

    var model = try builder.finishColumnTerms(.{ .zero_tolerance = 1e-9 });
    defer model.deinit();
    try model.matrix.validate();
    try std.testing.expectEqualSlices(usize, &.{ 0, 1, 3 }, model.matrix.col_starts);
    try std.testing.expectEqualSlices(f64, &.{ 3.0, 4.0, 2.0 }, model.matrix.values);
    try std.testing.expectEqual(@as(usize, 1), model.matrix.row_indices[0].toUsize());
    try std.testing.expectEqual(@as(usize, 0), model.matrix.row_indices[1].toUsize());
    try std.testing.expectEqual(@as(usize, 1), model.matrix.row_indices[2].toUsize());
    try std.testing.expectEqual(@as(usize, 0), builder.column_terms.capacity);
    try std.testing.expectEqual(@as(usize, 0), builder.column_chains.capacity);
}

test "LP column terms enforce row monotonicity and reject mixed storage" {
    var builder = Builder.init(std.testing.allocator);
    defer builder.deinit();
    builder.enableColumnTermStorage();
    _ = try builder.addColumn(.{ .name = "x" });
    _ = try builder.addRow(.{ .name = "r0" });
    _ = try builder.addRow(.{ .name = "r1" });
    try builder.addColumnTerm(1, 0, 1.0);
    try std.testing.expectError(error.InvalidDimensions, builder.addColumnTerm(0, 0, 2.0));
    try std.testing.expectError(error.InvalidDimensions, builder.addTerm(0, 0, 2.0));
    if (@sizeOf(foundation.HUInt) == 4) try std.testing.expectEqual(@as(usize, 16), @sizeOf(Builder.ColumnTerm));
}

test "packed model storage pools attributes and names" {
    var builder = Builder.init(std.testing.allocator);
    defer builder.deinit();
    builder.name = "pooled-model";
    _ = try builder.addColumn(.{ .name = "alpha", .cost = 2.0 });
    _ = try builder.addColumn(.{ .name = "beta", .lower = -1.0 });
    _ = try builder.addRow(.{ .name = "capacity", .upper = 4.0 });
    _ = try builder.addRow(.{ .name = null, .lower = 0.0 });

    var model = try builder.finish(.{});
    defer model.deinit();
    try std.testing.expectEqualStrings("pooled-model", model.name);
    try std.testing.expectEqualStrings("alpha", model.col_names[0].?);
    try std.testing.expectEqualStrings("beta", model.col_names[1].?);
    try std.testing.expectEqualStrings("capacity", model.row_names[0].?);
    try std.testing.expect(model.row_names[1] == null);
    try std.testing.expectEqual(@as(usize, 0), @intFromPtr(model.attribute_storage.ptr) % 64);

    const pool_begin = @intFromPtr(model.name_storage.ptr);
    const pool_end = pool_begin + model.name_storage.len;
    for ([_][]const u8{ model.name, model.col_names[0].?, model.col_names[1].?, model.row_names[0].? }) |name| {
        const address = @intFromPtr(name.ptr);
        try std.testing.expect(address >= pool_begin and address + name.len <= pool_end);
    }
}

test "keep_names false retains only pooled model name" {
    var builder = Builder.init(std.testing.allocator);
    defer builder.deinit();
    builder.name = "anonymous-columns";
    _ = try builder.addColumn(.{ .name = "x" });
    _ = try builder.addRow(.{ .name = "row" });

    var model = try builder.finish(.{ .keep_names = false });
    defer model.deinit();
    try std.testing.expectEqualStrings("anonymous-columns", model.name);
    try std.testing.expectEqual(@as(usize, 0), model.col_names.len);
    try std.testing.expectEqual(@as(usize, 0), model.row_names.len);
    try std.testing.expectEqual(model.name.len, model.name_storage.len);
    try std.testing.expectEqual(@as(usize, 0), builder.columns.capacity);
    try std.testing.expectEqual(@as(usize, 0), builder.rows.capacity);
}

test "packed storage supports an empty model" {
    var builder = Builder.init(std.testing.allocator);
    defer builder.deinit();
    var model = try builder.finish(.{});
    defer model.deinit();
    try std.testing.expectEqual(@as(usize, 0), model.col_cost.len);
    try std.testing.expectEqual(@as(usize, 0), model.row_lower.len);
    try std.testing.expectEqual(@as(usize, 0), model.matrix.nnz());
}

/// Exercise generic-term finalization under the testing failing allocator.
fn buildPackedForAllocationFailureTest(allocator: std.mem.Allocator) !void {
    var builder = Builder.init(allocator);
    defer builder.deinit();
    builder.name = "failure-test";
    _ = try builder.addColumn(.{ .name = "x" });
    _ = try builder.addColumn(.{ .name = "y" });
    _ = try builder.addRow(.{ .name = "row" });
    try builder.addTerm(0, 0, 1.0);
    try builder.addTerm(0, 1, 2.0);
    var model = try builder.finish(.{});
    defer model.deinit();
}

test "packed storage cleans up every allocation failure point" {
    try std.testing.checkAllAllocationFailures(std.testing.allocator, buildPackedForAllocationFailureTest, .{});
}

/// Exercise LP column-chain finalization under the testing failing allocator.
fn buildColumnTermsForAllocationFailureTest(allocator: std.mem.Allocator) !void {
    var builder = Builder.init(allocator);
    defer builder.deinit();
    builder.enableColumnTermStorage();
    _ = try builder.addColumn(.{ .name = "x" });
    _ = try builder.addColumn(.{ .name = "y" });
    _ = try builder.addRow(.{ .name = "row" });
    try builder.addColumnTerm(0, 0, 1.0);
    try builder.addColumnTerm(0, 1, 2.0);
    var model = try builder.finishColumnTerms(.{});
    defer model.deinit();
}

test "LP column term storage cleans up every allocation failure point" {
    try std.testing.checkAllAllocationFailures(std.testing.allocator, buildColumnTermsForAllocationFailureTest, .{});
}
