//! Owning model interchange representation returned by parsers.
//!
//! Numeric arrays and the canonical CSC matrix can be moved into a solver
//! model without reparsing. Attribute slices share one aligned packed
//! allocation, while the model, column, and row names share one string pool.
//! Public slices retain the same API and never expose the storage layout.
//! When names are discarded, `col_names` and `row_names` are empty slices, so
//! no per-object nullable-slice table is allocated.

const std = @import("std");
const matrix = @import("matrix");
const types = @import("types.zig");
const StringArena = @import("string_arena.zig").StringArena;

/// Fully owned canonical model produced by a format parser.
pub const ModelData = struct {
    /// Allocator that owns every allocation reachable from this value.
    allocator: std.mem.Allocator,
    /// Single 64-byte-aligned allocation backing all numeric/attribute slices.
    attribute_storage: []align(64) u8,
    /// Contiguous byte pool backing the model, column and row names.
    name_storage: []u8,
    /// Model name borrowed from `name_storage`.
    name: []u8,
    /// Original objective direction.
    objective_sense: types.ObjectiveSense = .minimize,
    /// Constant objective term.
    objective_offset: f64 = 0.0,
    /// Structural objective coefficients; length equals `matrix.num_cols`.
    col_cost: []f64,
    /// Structural column lower bounds.
    col_lower: []f64,
    /// Structural column upper bounds.
    col_upper: []f64,
    /// Structural column domains.
    col_type: []types.VariableType,
    /// Optional names borrowing from `name_storage`, or an empty table.
    col_names: []?[]u8,
    /// Row activity lower bounds; length equals `matrix.num_rows`.
    row_lower: []f64,
    /// Row activity upper bounds.
    row_upper: []f64,
    /// Optional names borrowing from `name_storage`, or an empty table.
    row_names: []?[]u8,
    /// Owned canonical constraint matrix.
    matrix: matrix.CscMatrix,

    /// Allocate packed attribute arrays and a contiguous name pool. Ownership
    /// of `owned_matrix` transfers to this function even when allocation fails.
    pub fn initPacked(
        allocator: std.mem.Allocator,
        owned_matrix: matrix.CscMatrix,
        model_name: []const u8,
        additional_name_bytes: usize,
        retain_object_names: bool,
        objective_sense: types.ObjectiveSense,
        objective_offset: f64,
    ) types.IoError!ModelData {
        var owned_csc = owned_matrix;
        errdefer owned_csc.deinit(allocator);
        const num_cols = owned_csc.num_cols;
        const num_rows = owned_csc.num_rows;

        const sizes = [_]usize{
            try arrayBytes(f64, num_cols),
            try arrayBytes(f64, num_cols),
            try arrayBytes(f64, num_cols),
            try arrayBytes(types.VariableType, num_cols),
            try arrayBytes(?[]u8, if (retain_object_names) num_cols else 0),
            try arrayBytes(f64, num_rows),
            try arrayBytes(f64, num_rows),
            try arrayBytes(?[]u8, if (retain_object_names) num_rows else 0),
        };
        const layout = packedLayout(sizes, .{
            64,
            64,
            64,
            @alignOf(types.VariableType),
            @alignOf(?[]u8),
            64,
            64,
            @alignOf(?[]u8),
        }) catch return error.InvalidDimensions;
        const attribute_storage = allocator.alignedAlloc(u8, .@"64", @max(layout.total, 1)) catch return error.OutOfMemory;
        errdefer allocator.free(attribute_storage);

        const total_name_bytes = std.math.add(usize, model_name.len, additional_name_bytes) catch return error.InvalidDimensions;
        const name_storage = allocator.alloc(u8, @max(total_name_bytes, 1)) catch return error.OutOfMemory;
        errdefer allocator.free(name_storage);
        var name_arena = StringArena.init(name_storage);
        const name = name_arena.append(model_name) catch return error.InvalidDimensions;

        const col_cost = sliceAt(f64, attribute_storage, layout.offsets[0], num_cols);
        const col_lower = sliceAt(f64, attribute_storage, layout.offsets[1], num_cols);
        const col_upper = sliceAt(f64, attribute_storage, layout.offsets[2], num_cols);
        const col_type = sliceAt(types.VariableType, attribute_storage, layout.offsets[3], num_cols);
        const col_names = sliceAt(?[]u8, attribute_storage, layout.offsets[4], if (retain_object_names) num_cols else 0);
        const row_lower = sliceAt(f64, attribute_storage, layout.offsets[5], num_rows);
        const row_upper = sliceAt(f64, attribute_storage, layout.offsets[6], num_rows);
        const row_names = sliceAt(?[]u8, attribute_storage, layout.offsets[7], if (retain_object_names) num_rows else 0);
        @memset(col_names, null);
        @memset(row_names, null);

        return .{
            .allocator = allocator,
            .attribute_storage = attribute_storage,
            .name_storage = name_storage,
            .name = name,
            .objective_sense = objective_sense,
            .objective_offset = objective_offset,
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

    /// Release packed attributes, names and CSC storage.
    pub fn deinit(self: *ModelData) void {
        const allocator = self.allocator;
        allocator.free(self.attribute_storage);
        allocator.free(self.name_storage);
        self.matrix.deinit(allocator);
        self.* = undefined;
    }

    /// Create a non-owning writer view valid until this model is mutated or freed.
    pub fn view(self: *const ModelData) types.ModelView {
        return .{
            .name = self.name,
            .objective_sense = self.objective_sense,
            .objective_offset = self.objective_offset,
            .col_cost = self.col_cost,
            .col_lower = self.col_lower,
            .col_upper = self.col_upper,
            .col_type = self.col_type,
            .col_names = @ptrCast(self.col_names),
            .row_lower = self.row_lower,
            .row_upper = self.row_upper,
            .row_names = @ptrCast(self.row_names),
            .matrix = self.matrix.view(),
        };
    }
};

/// Compute the byte count of a `T[len]` array with overflow checking.
fn arrayBytes(comptime T: type, len: usize) types.IoError!usize {
    return std.math.mul(usize, @sizeOf(T), len) catch error.InvalidDimensions;
}

/// Lay out eight arrays consecutively while respecting each requested alignment.
fn packedLayout(sizes: [8]usize, alignments: [8]usize) error{Overflow}!struct { offsets: [8]usize, total: usize } {
    var offsets: [8]usize = undefined;
    var cursor: usize = 0;
    for (sizes, alignments, 0..) |size, alignment, index| {
        if (alignment > 1 and cursor > std.math.maxInt(usize) - (alignment - 1)) return error.Overflow;
        cursor = std.mem.alignForward(usize, cursor, alignment);
        offsets[index] = cursor;
        cursor = std.math.add(usize, cursor, size) catch return error.Overflow;
    }
    return .{ .offsets = offsets, .total = cursor };
}

/// Reinterpret an aligned subrange of packed attribute storage as `[]T`.
///
/// Callers must obtain `offset` and `len` from `packedLayout`; this helper
/// deliberately performs no runtime bounds check in the construction hot path.
fn sliceAt(comptime T: type, storage: []align(64) u8, offset: usize, len: usize) []T {
    const pointer: [*]T = @ptrCast(@alignCast(storage.ptr + offset));
    return pointer[0..len];
}
