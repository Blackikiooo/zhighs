//! Data-oriented sparse basis-matrix assembly for revised simplex.
//!
//! A simplex basis contains exactly one column per row, selected from the
//! structural matrix plus logical and temporary artificial identity columns.
//! This module materializes that basis as compact CSC in caller-retained SoA
//! buffers.  It deliberately owns no pivoting policy: sparse LU, dense-LU
//! oracles, and factorization benchmarks can consume the same representation.

const std = @import("std");
const foundation = @import("foundation");
const csc = @import("csc.zig");
const target = @import("target_policy.zig");

const Offset = foundation.HUInt;
const RowId = foundation.RowId;

pub const SparseBasisError = error{
    DimensionMismatch,
    InvalidBasisColumn,
    CapacityOverflow,
    NonFiniteValue,
    OutOfMemory,
};

/// Borrowed compact CSC view valid until its `SparseBasisBuffers` owner grows
/// or is deinitialized.  Offsets use the configured HiGHS integer width, which
/// halves offset traffic in the default w32 build compared with `usize`.
pub const SparseBasisView = struct {
    dimension: usize,
    starts: []const Offset,
    rows: []const RowId,
    values: []const f64,

    pub inline fn nnz(self: SparseBasisView) usize {
        return self.values.len;
    }
};

/// Retaining SoA workspace for repeated basis reinversion.
///
/// Starts, row IDs and values are independent aligned streams because the
/// symbolic phase consumes starts/rows while numerical kernels consume
/// rows/values.  Keeping these streams separate avoids loading unused fields
/// and lets the allocator preserve capacity across simplex reinversions.
pub const SparseBasisBuffers = struct {
    allocator: std.mem.Allocator,
    starts: []align(64) Offset = &.{},
    rows: []align(64) RowId = &.{},
    values: []align(64) f64 = &.{},
    starts_capacity: usize = 0,
    entry_capacity: usize = 0,
    dimension: usize = 0,
    nonzeros: usize = 0,

    pub fn init(allocator: std.mem.Allocator) SparseBasisBuffers {
        return .{ .allocator = allocator };
    }

    pub fn deinit(self: *SparseBasisBuffers) void {
        if (self.starts_capacity != 0) self.allocator.free(self.starts.ptr[0..self.starts_capacity]);
        if (self.entry_capacity != 0) {
            self.allocator.free(self.rows.ptr[0..self.entry_capacity]);
            self.allocator.free(self.values.ptr[0..self.entry_capacity]);
        }
        self.* = .{ .allocator = self.allocator };
    }

    /// Assemble a basis without retaining any input slices.
    ///
    /// Global column numbering matches the simplex engine:
    /// structural `[0,A.num_cols)`, logical identity columns immediately
    /// afterwards, then one signed artificial identity column per row.
    pub fn assemble(
        self: *SparseBasisBuffers,
        matrix: csc.CscView,
        basic_index: []const u32,
        row_scale: []const f64,
        column_scale: []const f64,
        artificial_sign: []const f64,
    ) SparseBasisError!SparseBasisView {
        const n = matrix.num_rows;
        if (basic_index.len != n or row_scale.len != n or column_scale.len != matrix.num_cols or artificial_sign.len != n)
            return error.DimensionMismatch;

        var required_entries: usize = 0;
        for (basic_index) |global_column_u32| {
            const global_column: usize = global_column_u32;
            if (global_column < matrix.num_cols) {
                required_entries = std.math.add(
                    usize,
                    required_entries,
                    matrix.col_starts[global_column + 1] - matrix.col_starts[global_column],
                ) catch return error.CapacityOverflow;
            } else if (global_column >= matrix.num_cols + 2 * n) {
                return error.InvalidBasisColumn;
            } else {
                required_entries = std.math.add(usize, required_entries, 1) catch return error.CapacityOverflow;
            }
        }
        if (required_entries > std.math.maxInt(Offset)) return error.CapacityOverflow;
        try self.ensureCapacity(n + 1, required_entries);

        var output: usize = 0;
        self.starts[0] = 0;
        for (basic_index, 0..) |global_column_u32, basis_column| {
            // Prefetch only structural columns. Logical/artificial columns are
            // synthesized locally and touching arbitrary A storage would add
            // cache pollution rather than hide useful latency.
            const future = basis_column + target.sparse_column_prefetch_distance;
            if (future < basic_index.len) {
                const future_column: usize = basic_index[future];
                if (future_column < matrix.num_cols) {
                    const future_begin = matrix.col_starts[future_column];
                    if (future_begin < matrix.values.len) {
                        @prefetch(&matrix.row_indices[future_begin], .{ .locality = 2 });
                        @prefetch(&matrix.values[future_begin], .{ .locality = 2 });
                    }
                }
            }

            const global_column: usize = global_column_u32;
            if (global_column < matrix.num_cols) {
                const begin = matrix.col_starts[global_column];
                const end = matrix.col_starts[global_column + 1];
                for (matrix.row_indices[begin..end], matrix.values[begin..end]) |row, coefficient| {
                    if (!target.retainsModelCoefficient(coefficient)) continue;
                    const scaled = coefficient * row_scale[row.toUsize()] * column_scale[global_column];
                    if (!std.math.isFinite(scaled)) return error.NonFiniteValue;
                    // Free rows have a zero normalization scale. Do not publish
                    // explicit zeros: symbolic counts and Markowitz merit must
                    // describe the numerical basis that LU will actually see.
                    if (scaled == 0.0) continue;
                    self.rows[output] = row;
                    self.values[output] = scaled;
                    output += 1;
                }
            } else {
                const internal = global_column - matrix.num_cols;
                const row = if (internal < n) internal else internal - n;
                const value: f64 = if (internal < n) 1.0 else artificial_sign[row];
                if (!std.math.isFinite(value) or value == 0.0) return error.NonFiniteValue;
                self.rows[output] = RowId.fromUsizeAssumeValid(row);
                self.values[output] = value;
                output += 1;
            }
            self.starts[basis_column + 1] = @intCast(output);
        }
        self.dimension = n;
        self.nonzeros = output;
        return self.view();
    }

    pub inline fn view(self: *const SparseBasisBuffers) SparseBasisView {
        return .{
            .dimension = self.dimension,
            .starts = self.starts[0 .. self.dimension + 1],
            .rows = self.rows[0..self.nonzeros],
            .values = self.values[0..self.nonzeros],
        };
    }

    fn ensureCapacity(self: *SparseBasisBuffers, starts_needed: usize, entries_needed: usize) SparseBasisError!void {
        if (starts_needed > self.starts_capacity) {
            const capacity = growCapacity(self.starts_capacity, starts_needed) catch return error.CapacityOverflow;
            const replacement = self.allocator.alignedAlloc(Offset, .@"64", capacity) catch return error.OutOfMemory;
            if (self.starts_capacity != 0) self.allocator.free(self.starts.ptr[0..self.starts_capacity]);
            self.starts = replacement;
            self.starts_capacity = capacity;
        }
        if (entries_needed > self.entry_capacity) {
            const capacity = growCapacity(self.entry_capacity, entries_needed) catch return error.CapacityOverflow;
            const next_rows = self.allocator.alignedAlloc(RowId, .@"64", capacity) catch return error.OutOfMemory;
            errdefer self.allocator.free(next_rows);
            const next_values = self.allocator.alignedAlloc(f64, .@"64", capacity) catch return error.OutOfMemory;
            if (self.entry_capacity != 0) {
                self.allocator.free(self.rows.ptr[0..self.entry_capacity]);
                self.allocator.free(self.values.ptr[0..self.entry_capacity]);
            }
            self.rows = next_rows;
            self.values = next_values;
            self.entry_capacity = capacity;
        }
    }
};

fn growCapacity(current: usize, required: usize) error{Overflow}!usize {
    var capacity = @max(current, 8);
    while (capacity < required) capacity = std.math.add(usize, capacity, capacity / 2 + 8) catch return error.Overflow;
    return capacity;
}

test "basis assembly mixes structural logical and artificial columns" {
    const rows = [_]RowId{
        RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(2),
        RowId.fromUsizeAssumeValid(1), RowId.fromUsizeAssumeValid(2),
    };
    const matrix = csc.CscView.initAssumeValid(
        3,
        2,
        &[_]usize{ 0, 2, 4 },
        &rows,
        &[_]f64{ 2, 4, 3, 5 },
    );
    var buffers = SparseBasisBuffers.init(std.testing.allocator);
    defer buffers.deinit();
    const basis = try buffers.assemble(matrix, &[_]u32{ 0, 3, 6 }, &[_]f64{ 1, -1, 0.5 }, &[_]f64{ 1, 1, 1 }, &[_]f64{ 0, 0, -1 });
    try std.testing.expectEqualSlices(Offset, &[_]Offset{ 0, 2, 3, 4 }, basis.starts);
    try std.testing.expectEqualSlices(RowId, &[_]RowId{ rows[0], rows[1], RowId.fromUsizeAssumeValid(1), RowId.fromUsizeAssumeValid(2) }, basis.rows);
    try std.testing.expectEqualSlices(f64, &[_]f64{ 2, 2, 1, -1 }, basis.values);
}

test "basis assembly retains capacity across reinversions" {
    const row = [_]RowId{RowId.fromUsizeAssumeValid(0)};
    const matrix = csc.CscView.initAssumeValid(1, 1, &[_]usize{ 0, 1 }, &row, &[_]f64{2});
    var buffers = SparseBasisBuffers.init(std.testing.allocator);
    defer buffers.deinit();
    _ = try buffers.assemble(matrix, &[_]u32{0}, &[_]f64{1}, &[_]f64{1}, &[_]f64{0});
    const starts_pointer = buffers.starts.ptr;
    const values_pointer = buffers.values.ptr;
    _ = try buffers.assemble(matrix, &[_]u32{1}, &[_]f64{1}, &[_]f64{1}, &[_]f64{0});
    try std.testing.expectEqual(starts_pointer, buffers.starts.ptr);
    try std.testing.expectEqual(values_pointer, buffers.values.ptr);
}
