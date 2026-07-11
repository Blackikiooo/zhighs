//! Owning authoritative CSC matrix with a rebuildable CSR cache.
//!
//! Structural mutations must go through MatrixStore methods so revision and
//! derived caches remain coherent. Solver algorithms receive borrowed views;
//! they do not own or mutate the authoritative CSC storage.

const std = @import("std");
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
    csr_cache: ?csr_module.CsrCache = null,

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
        if (self.csr_cache) |*cache| cache.deinit(allocator);
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
        if (self.csr_cache) |cache| {
            if (cache.isCurrent(self.matrix_revision))
                return cache.viewAssumeCurrent();
        }

        // Build first, then release the stale cache. If allocation fails, the
        // authoritative matrix and any old cache remain fully valid.
        var replacement = try csr_module.CsrCache.buildAssumeValid(
            allocator,
            self.matrix_storage,
            self.matrix_revision,
        );
        errdefer replacement.deinit(allocator);

        if (self.csr_cache) |*old_cache| old_cache.deinit(allocator);
        self.csr_cache = replacement;
        return self.csr_cache.?.viewAssumeCurrent();
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

        if (self.csr_cache) |*cache| cache.deinit(allocator);
        self.csr_cache = null;
        self.matrix_storage.deinit(allocator);
        self.matrix_storage = replacement;
        self.matrix_revision += 1;
    }
};

test "matrix store builds CSR lazily and reuses the cache" {
    var csc = try csc_module.CscMatrix.initZero(std.testing.allocator, 2, 3);
    var model = MatrixStore.initAssumeValid(csc);
    // Ownership moved to model; prevent accidental use of the local value.
    csc = undefined;
    defer model.deinit(std.testing.allocator);

    try std.testing.expect(model.csr_cache == null);
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
    try std.testing.expect(model.csr_cache == null);
    const rebuilt = try model.csr(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 3), rebuilt.num_rows);
    try std.testing.expectEqual(@as(usize, 2), rebuilt.num_cols);
    try std.testing.expect(model.csr_cache.?.isCurrent(1));
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
