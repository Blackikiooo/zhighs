//! Owning authoritative CSC matrix with a rebuildable CSR cache.
//!
//! Structural mutations must go through MatrixStore methods so revision and
//! derived caches remain coherent. Solver algorithms receive borrowed views;
//! they do not own or mutate the authoritative CSC storage.

const std = @import("std");
const foundation = @import("foundation");
const csc_module = @import("csc.zig");
const csr_module = @import("csr_view.zig");

pub const MatrixStoreError = csc_module.MatrixError || error{
    RevisionOverflow,
};

pub const MatrixStore = struct {
    /// Authoritative matrix. Treat this field as private even though Zig struct
    /// fields are visible: replacing it directly bypasses cache invalidation.
    matrix_storage: csc_module.CscMatrix,
    matrix_revision: u64 = 0,
    /// High-water-mark output and cursor storage retained across revisions.
    /// `csr_cache_revision == null` means contents must be rebuilt before use.
    csr_buffers: ?csr_module.CsrBuffers = null,
    csr_cache_revision: ?u64 = null,

    const Self = @This();

    /// Takes ownership of a matrix after validating its canonical invariants.
    /// On validation failure ownership remains with the caller.
    pub fn init(matrix: csc_module.CscMatrix) MatrixStoreError!Self {
        try matrix.validate();
        return initAssumeValid(matrix);
    }

    /// Takes ownership of a canonical matrix without another O(nnz) scan.
    pub inline fn initAssumeValid(matrix: csc_module.CscMatrix) Self {
        return .{ .matrix_storage = matrix };
    }

    pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
        if (self.csr_buffers) |*buffers| buffers.deinit(allocator);
        self.matrix_storage.deinit(allocator);
        self.* = undefined;
    }

    /// Read-only-by-contract access to the authoritative CSC matrix.
    pub inline fn csc(self: *const Self) *const csc_module.CscMatrix {
        return &self.matrix_storage;
    }

    pub inline fn matrixRevision(self: Self) u64 {
        return self.matrix_revision;
    }

    /// Returns a cached CSR view, building it only when absent or stale.
    ///
    /// A returned view remains valid until the next mutable MatrixStore
    /// operation or deinit. The cache-hit path performs one revision comparison
    /// and no allocation.
    pub fn csr(self: *Self, allocator: std.mem.Allocator) (std.mem.Allocator.Error || MatrixStoreError)!csr_module.CsrView {
        if (self.csr_cache_revision) |revision| {
            if (revision == self.matrix_revision)
                return self.csrViewAssumeCurrent();
        }

        const rows = self.matrix_storage.num_rows;
        const nnz = self.matrix_storage.nnz();
        if (nnz > std.math.maxInt(foundation.HUInt)) return error.DimensionTooLarge;
        const has_capacity = if (self.csr_buffers) |buffers|
            buffers.row_starts.len >= rows + 1 and
                buffers.col_indices.len >= nnz and
                buffers.values.len >= nnz and
                buffers.cursor.len >= rows
        else
            false;

        if (!has_capacity) {
            // Allocate and fill the replacement before releasing the old
            // high-water buffers, preserving strong failure safety.
            var replacement = try csr_module.CsrBuffers.init(allocator, rows, nnz);
            errdefer replacement.deinit(allocator);
            try csr_module.fillFromCscAssumeValid(
                self.matrix_storage,
                replacement.row_starts,
                replacement.col_indices,
                replacement.values,
                replacement.cursor,
            );
            if (self.csr_buffers) |*old| old.deinit(allocator);
            self.csr_buffers = replacement;
        } else {
            const buffers = &self.csr_buffers.?;
            try csr_module.fillFromCscAssumeValid(
                self.matrix_storage,
                buffers.row_starts[0 .. rows + 1],
                buffers.col_indices[0..nnz],
                buffers.values[0..nnz],
                buffers.cursor[0..rows],
            );
        }
        self.csr_cache_revision = self.matrix_revision;
        return self.csrViewAssumeCurrent();
    }

    fn csrViewAssumeCurrent(self: *const Self) csr_module.CsrView {
        const buffers = self.csr_buffers.?;
        const rows = self.matrix_storage.num_rows;
        const nnz = self.matrix_storage.nnz();
        return csr_module.CsrView.initAssumeValid(
            rows,
            self.matrix_storage.num_cols,
            buffers.row_starts[0 .. rows + 1],
            buffers.col_indices[0..nnz],
            buffers.values[0..nnz],
        );
    }

    /// Replaces the authoritative matrix after validation and invalidates CSR.
    /// Ownership transfers only after validation succeeds.
    pub fn replaceMatrix(self: *Self, allocator: std.mem.Allocator, replacement: csc_module.CscMatrix) MatrixStoreError!void {
        try replacement.validate();
        try self.replaceMatrixAssumeValid(allocator, replacement);
    }

    /// Trusted replacement path for a matrix emitted by MatrixBuilder.freeze.
    pub fn replaceMatrixAssumeValid(self: *Self, allocator: std.mem.Allocator, replacement: csc_module.CscMatrix) MatrixStoreError!void {
        if (self.matrix_revision == std.math.maxInt(u64)) return error.RevisionOverflow;

        self.csr_cache_revision = null;
        self.matrix_storage.deinit(allocator);
        self.matrix_storage = replacement;
        self.matrix_revision += 1;
    }

    /// Updates existing matrix values without copying structural CSC streams.
    /// The caller guarantees matching slice lengths, valid positions, finite
    /// nonzero values, and unique positions. All fallible work happens before
    /// cache invalidation, so revision overflow leaves the store unchanged.
    pub fn updateValuesAtPositionsAssumeValid(self: *Self, allocator: std.mem.Allocator, positions: []const usize, values: []const f64) MatrixStoreError!void {
        _ = allocator; // Retained for API symmetry with structural mutations.
        std.debug.assert(positions.len == values.len);
        if (positions.len == 0) return;
        if (self.matrix_revision == std.math.maxInt(u64)) return error.RevisionOverflow;
        for (positions, values) |position, value| {
            std.debug.assert(position < self.matrix_storage.values.len);
            std.debug.assert(std.math.isFinite(value) and value != 0.0);
        }

        self.csr_cache_revision = null;
        for (positions, values) |position, value|
            self.matrix_storage.values[position] = value;
        self.matrix_revision += 1;
    }
};

test "matrix store builds CSR lazily and reuses the cache" {
    var csc = try csc_module.CscMatrix.initZero(std.testing.allocator, 2, 3);
    var model = MatrixStore.initAssumeValid(csc);
    // Ownership moved to model; prevent accidental use of the local value.
    csc = undefined;
    defer model.deinit(std.testing.allocator);

    try std.testing.expect(model.csr_cache_revision == null);
    const first = try model.csr(std.testing.allocator);
    const second = try model.csr(std.testing.allocator);
    try std.testing.expectEqual(first.row_starts.ptr, second.row_starts.ptr);
    try std.testing.expectEqual(@as(u64, 0), model.matrixRevision());
}

test "matrix replacement increments revision and rebuilds CSR" {
    var initial = try csc_module.CscMatrix.initZero(std.testing.allocator, 1, 1);
    var model = MatrixStore.initAssumeValid(initial);
    initial = undefined;
    defer model.deinit(std.testing.allocator);

    _ = try model.csr(std.testing.allocator);
    var replacement = try csc_module.CscMatrix.initZero(std.testing.allocator, 3, 2);
    try model.replaceMatrixAssumeValid(std.testing.allocator, replacement);
    replacement = undefined;

    try std.testing.expectEqual(@as(u64, 1), model.matrixRevision());
    try std.testing.expect(model.csr_cache_revision == null);
    const rebuilt = try model.csr(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 3), rebuilt.num_rows);
    try std.testing.expectEqual(@as(usize, 2), rebuilt.num_cols);
    try std.testing.expectEqual(@as(?u64, 1), model.csr_cache_revision);
    try std.testing.expectEqual(@as(usize, 3), model.csr_buffers.?.cursor.len);

    const storage_ptr = model.csr_buffers.?.storage.ptr;
    var smaller = try csc_module.CscMatrix.initZero(std.testing.allocator, 1, 1);
    try model.replaceMatrixAssumeValid(std.testing.allocator, smaller);
    smaller = undefined;
    _ = try model.csr(std.testing.allocator);
    try std.testing.expectEqual(storage_ptr, model.csr_buffers.?.storage.ptr);
    try std.testing.expectEqual(@as(usize, 3), model.csr_buffers.?.cursor.len);
}

test "revision overflow rejects replacement without changing ownership" {
    var initial = try csc_module.CscMatrix.initZero(std.testing.allocator, 1, 1);
    var model = MatrixStore.initAssumeValid(initial);
    initial = undefined;
    defer model.deinit(std.testing.allocator);
    model.matrix_revision = std.math.maxInt(u64);

    var replacement = try csc_module.CscMatrix.initZero(std.testing.allocator, 2, 2);
    defer replacement.deinit(std.testing.allocator);
    try std.testing.expectError(error.RevisionOverflow, model.replaceMatrixAssumeValid(std.testing.allocator, replacement));
    try std.testing.expectEqual(@as(usize, 1), model.csc().num_rows);
}
