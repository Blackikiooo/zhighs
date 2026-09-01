//! Vector-native numerical linear algebra boundary for revised simplex.
//!
//! `SimplexNla` is concrete so Zig can inline backend dispatch. It owns the
//! factorization and consumes the same `IndexedVector(f64)` instances that
//! PRICE, CHUZR and UPDATE exchange. A backend that cannot publish exact
//! support marks the result dense instead of fabricating indices by scanning.

const std = @import("std");
const factorization_module = @import("factorization.zig");
const indexed_vector = @import("indexed_vector.zig");

pub const Vector = indexed_vector.IndexedVector(f64);
pub const NlaError = factorization_module.FactorizationError || indexed_vector.IndexedVectorError;

pub const SimplexNla = struct {
    factorization: factorization_module.Factorization,

    pub fn init(allocator: std.mem.Allocator) SimplexNla {
        return .{ .factorization = factorization_module.Factorization.init(allocator) };
    }

    pub fn deinit(self: *SimplexNla) void {
        self.factorization.deinit();
    }

    pub fn factorize(self: *SimplexNla, dimension: usize, matrix_data: []const f64) NlaError!void {
        try self.factorization.factorize(dimension, matrix_data);
    }

    pub fn factorizeIdentity(self: *SimplexNla, dimension: usize) NlaError!void {
        try self.factorization.factorizeIdentity(dimension);
    }

    /// Compute FTRAN while retaining producer support only when the selected
    /// sparse backend actually supplies it.
    pub fn ftran(self: *SimplexNla, vector: *Vector) NlaError!void {
        if (vector.state == .indexed) {
            if (vector.index_count == 0) return;
            if (!self.factorization.usesHyperSparseRhs(vector.index_count))
                vector.materializeZeroComplement();
            try self.factorization.solveSparse(vector.values, vector.indices[0..vector.index_count]);
            const output = self.factorization.lastSparseFtranIndices();
            if (output.len != 0) {
                @memcpy(vector.indices[0..output.len], output);
                vector.index_count = output.len;
                return;
            }
        } else {
            try self.factorization.solve(vector.values);
        }
        vector.state = .dense;
        vector.index_count = 0;
    }

    /// Compute BTRAN with the same explicit indexed-to-dense state contract.
    pub fn btran(self: *SimplexNla, vector: *Vector) NlaError!void {
        if (vector.state == .indexed) {
            if (vector.index_count == 0) return;
            if (!self.factorization.usesHyperSparseRhs(vector.index_count))
                vector.materializeZeroComplement();
            try self.factorization.solveTransposeSparse(vector.values, vector.indices[0..vector.index_count]);
            const output = self.factorization.lastSparseSolveIndices();
            if (output.len != 0) {
                @memcpy(vector.indices[0..output.len], output);
                vector.index_count = output.len;
                return;
            }
        } else {
            try self.factorization.solveTranspose(vector.values);
        }
        vector.state = .dense;
        vector.index_count = 0;
    }

    /// Entering-column FTRAN preserves the captured Aq required by the next
    /// Forrest--Tomlin update and retains the mutable-U producer order.
    pub fn ftranForUpdate(self: *SimplexNla, vector: *Vector) NlaError!void {
        if (vector.state == .indexed and self.factorization.supportsSparseUpdateRhs()) {
            try self.factorization.solveForUpdateSparse(vector.values, vector.indices[0..vector.index_count]);
            const output = self.factorization.lastUpdateFtranIndices();
            @memcpy(vector.indices[0..output.len], output);
            vector.index_count = output.len;
            return;
        } else {
            vector.materializeZeroComplement();
            try self.factorization.solveForUpdate(vector.values);
        }
        vector.state = .dense;
        vector.index_count = 0;
    }

    /// Install the basis update from the same Aq vector consumed by FTRAN.
    pub fn update(
        self: *SimplexNla,
        col_aq: *const Vector,
        leaving_row: u32,
        entering_col: u32,
        column_scale: f64,
    ) factorization_module.FactorizationError!void {
        const direction_indices: ?[]const u32 = if (col_aq.values.len >= 128 and
            col_aq.state == .indexed and
            col_aq.index_count * 10 <= col_aq.values.len)
            col_aq.activeIndices()
        else
            null;
        try self.factorization.update(.{
            .leaving_row = leaving_row,
            .entering_col = entering_col,
            .direction = col_aq.values,
            .direction_indices = direction_indices,
            .direction_is_finite = true,
            .column_scale = column_scale,
        });
    }

    pub fn retainedBytes(self: *const SimplexNla) usize {
        return self.factorization.requestedBytes();
    }
};

test "simplex NLA keeps empty support and marks dense backend results" {
    var nla = SimplexNla.init(std.testing.allocator);
    defer nla.deinit();
    try nla.factorize(2, &.{ 2, 0, 0, 4 });
    var vector = try Vector.init(std.testing.allocator, 2);
    defer vector.deinit();

    try nla.ftran(&vector);
    try std.testing.expectEqual(indexed_vector.IndexState.indexed, vector.state);
    try std.testing.expectEqual(@as(usize, 0), vector.index_count);

    vector.values[0] = 6;
    vector.values[1] = 99; // Dense LU must not read outside indexed support.
    try vector.appendIndex(0);
    try nla.ftran(&vector);
    try std.testing.expectEqual(indexed_vector.IndexState.dense, vector.state);
    try std.testing.expectApproxEqAbs(@as(f64, 3), vector.values[0], 1e-12);
    try std.testing.expectEqual(@as(f64, 0), vector.values[1]);
}

test "simplex NLA retains sparse producer support for identity FTRAN and BTRAN" {
    const n = 64;
    var nla = SimplexNla.init(std.testing.allocator);
    defer nla.deinit();
    try nla.factorizeIdentity(n);
    var vector = try Vector.init(std.testing.allocator, n);
    defer vector.deinit();

    vector.values[7] = 3;
    try vector.appendIndex(7);
    try nla.ftran(&vector);
    try std.testing.expectEqual(indexed_vector.IndexState.indexed, vector.state);
    try std.testing.expectEqualSlices(u32, &.{7}, vector.activeIndices().?);
    try std.testing.expectEqual(@as(f64, 3), vector.values[7]);

    vector.clear();
    vector.values[19] = -5;
    try vector.appendIndex(19);
    try nla.btran(&vector);
    try std.testing.expectEqual(indexed_vector.IndexState.indexed, vector.state);
    try std.testing.expectEqualSlices(u32, &.{19}, vector.activeIndices().?);
    try std.testing.expectEqual(@as(f64, -5), vector.values[19]);
}

test "simplex NLA update consumes the captured Aq vector" {
    var nla = SimplexNla.init(std.testing.allocator);
    defer nla.deinit();
    try nla.factorizeIdentity(2);
    var col_aq = try Vector.init(std.testing.allocator, 2);
    defer col_aq.deinit();
    col_aq.values[0] = 2;
    col_aq.values[1] = 1;
    try col_aq.appendIndex(0);
    try col_aq.appendIndex(1);
    try nla.ftranForUpdate(&col_aq);
    try nla.update(&col_aq, 0, 0, 1);
    try std.testing.expectEqual(@as(usize, 1), nla.factorization.update_count);
}

test "sparse update FTRAN publishes support introduced by FT corrections" {
    const n = 64;
    var nla = SimplexNla.init(std.testing.allocator);
    defer nla.deinit();
    try nla.factorizeIdentity(n);

    // Replace identity column zero by e_0 + 0.25 e_5. The first update
    // captures this Aq and creates a Forrest--Tomlin correction.
    var entering = try Vector.init(std.testing.allocator, n);
    defer entering.deinit();
    entering.values[0] = 1;
    entering.values[5] = 0.25;
    try entering.appendIndex(0);
    try entering.appendIndex(5);
    try nla.ftranForUpdate(&entering);
    try nla.update(&entering, 0, 0, 1);

    // Solving the next unit RHS must scatter through that correction. The
    // published support is consumed directly by sparse DSE and clearing, so
    // every numerical nonzero must be represented even when absent in input.
    var corrected = try Vector.init(std.testing.allocator, n);
    defer corrected.deinit();
    corrected.values[0] = 1;
    try corrected.appendIndex(0);
    try nla.ftranForUpdate(&corrected);
    const support = corrected.activeIndices().?;
    var marked = [_]bool{false} ** n;
    for (support) |row| marked[row] = true;
    for (corrected.values, 0..) |value, row| {
        if (value != 0.0) try std.testing.expect(marked[row]);
    }
    try std.testing.expectApproxEqAbs(@as(f64, 1), corrected.values[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, -0.25), corrected.values[5], 1e-12);
}
