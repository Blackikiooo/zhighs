//! Global valid-cut storage with canonical sparse rows and aging metadata.

const std = @import("std");
const matrix = @import("matrix");
const foundation = @import("foundation");

pub const CutId = enum(usize) { _ };

const Metadata = struct {
    rhs: f64,
    hash: u64,
    efficacy: f64,
    age: u16,
    active: bool,
};

const MetadataList = std.MultiArrayList(Metadata);

pub const AddResult = union(enum) {
    added: CutId,
    existing: CutId,
};

pub const CutPool = struct {
    rows: matrix.DynamicRowMatrix,
    metadata: MetadataList = .empty,
    revision: usize = 0,

    pub fn init(allocator: std.mem.Allocator, num_columns: usize) !CutPool {
        return .{ .rows = try matrix.DynamicRowMatrix.init(allocator, num_columns) };
    }

    pub fn deinit(self: *CutPool, allocator: std.mem.Allocator) void {
        self.rows.deinit(allocator);
        self.metadata.deinit(allocator);
        self.* = undefined;
    }

    pub fn len(self: *const CutPool) usize {
        return self.metadata.len;
    }

    pub fn activeCount(self: *const CutPool) usize {
        var count: usize = 0;
        for (self.metadata.items(.active)) |active| count += @intFromBool(active);
        return count;
    }

    pub fn rhs(self: *const CutPool, id: CutId) f64 {
        return self.metadata.items(.rhs)[@intFromEnum(id)];
    }

    pub fn row(self: *const CutPool, id: CutId) matrix.SparseVectorView(foundation.ColId) {
        return self.rows.rowAssumeValid(@intFromEnum(id));
    }

    pub fn add(
        self: *CutPool,
        allocator: std.mem.Allocator,
        cut: matrix.SparseVectorView(foundation.ColId),
        rhs_value: f64,
        efficacy: f64,
    ) !AddResult {
        if (!std.math.isFinite(rhs_value) or !std.math.isFinite(efficacy) or efficacy < 0.0) return error.NonFiniteValue;
        try cut.validate();
        if (cut.dimension != self.rows.num_cols or cut.nnz() == 0) return error.DimensionMismatch;
        const hash = hashCut(cut, rhs_value);
        const hashes = self.metadata.items(.hash);
        for (hashes, 0..) |existing_hash, index| {
            if (existing_hash != hash or self.metadata.items(.rhs)[index] != rhs_value) continue;
            const existing = self.rows.rowAssumeValid(index);
            if (!std.mem.eql(foundation.ColId, existing.indices, cut.indices) or !std.mem.eql(f64, existing.values, cut.values)) continue;
            const was_active = self.metadata.items(.active)[index];
            self.metadata.items(.active)[index] = true;
            self.metadata.items(.age)[index] = 0;
            self.metadata.items(.efficacy)[index] = @max(self.metadata.items(.efficacy)[index], efficacy);
            if (!was_active) self.revision +%= 1;
            return .{ .existing = @enumFromInt(index) };
        }

        try self.rows.reserve(allocator, 1, cut.nnz());
        try self.metadata.ensureUnusedCapacity(allocator, 1);
        const id: CutId = @enumFromInt(self.metadata.len);
        _ = self.rows.appendRowPreReserved(cut);
        self.metadata.appendAssumeCapacity(.{ .rhs = rhs_value, .hash = hash, .efficacy = efficacy, .age = 0, .active = true });
        self.revision +%= 1;
        return .{ .added = id };
    }

    /// Ages inactive-at-the-current-point cuts and returns newly expired rows.
    pub fn ageAt(self: *CutPool, primal: []const f64, violation_tolerance: f64, max_age: u16) usize {
        std.debug.assert(primal.len == self.rows.num_cols);
        var expired: usize = 0;
        for (0..self.len()) |index| {
            if (!self.metadata.items(.active)[index]) continue;
            const row_view = self.rows.rowAssumeValid(index);
            var activity: f64 = 0.0;
            for (row_view.indices, row_view.values) |column, value| activity = @mulAdd(f64, value, primal[column.toUsize()], activity);
            if (activity > self.metadata.items(.rhs)[index] + violation_tolerance) {
                self.metadata.items(.age)[index] = 0;
            } else if (self.metadata.items(.age)[index] >= max_age) {
                self.metadata.items(.active)[index] = false;
                self.revision +%= 1;
                expired += 1;
            } else {
                self.metadata.items(.age)[index] += 1;
            }
        }
        return expired;
    }

    /// Reclaims inactive sparse rows at an explicit synchronization point.
    pub fn compact(self: *CutPool, allocator: std.mem.Allocator) !void {
        if (self.activeCount() == self.len()) return;
        const next_revision = self.revision +% 1;
        var replacement = try CutPool.init(allocator, self.rows.num_cols);
        errdefer replacement.deinit(allocator);
        for (0..self.len()) |index| {
            if (!self.metadata.items(.active)[index]) continue;
            const row_view = self.rows.rowAssumeValid(index);
            _ = try replacement.add(
                allocator,
                row_view,
                self.metadata.items(.rhs)[index],
                self.metadata.items(.efficacy)[index],
            );
            replacement.metadata.items(.age)[replacement.len() - 1] = self.metadata.items(.age)[index];
        }
        self.deinit(allocator);
        self.* = replacement;
        self.revision = next_revision;
    }

    pub fn appendActiveToCsc(self: *const CutPool, allocator: std.mem.Allocator, base: matrix.CscMatrix) !matrix.CscMatrix {
        var active_rows = try matrix.DynamicRowMatrix.init(allocator, self.rows.num_cols);
        defer active_rows.deinit(allocator);
        try active_rows.reserve(allocator, self.activeCount(), self.rows.nnz());
        for (0..self.len()) |index| {
            if (self.metadata.items(.active)[index]) _ = active_rows.appendRowPreReserved(self.rows.rowAssumeValid(index));
        }
        return active_rows.appendToCsc(allocator, base);
    }

    pub fn fillActiveUpper(self: *const CutPool, output: []f64) void {
        std.debug.assert(output.len == self.activeCount());
        var position: usize = 0;
        for (self.metadata.items(.rhs), self.metadata.items(.active)) |rhs_value, active| {
            if (!active) continue;
            output[position] = rhs_value;
            position += 1;
        }
    }
};

fn hashCut(cut: matrix.SparseVectorView(foundation.ColId), rhs: f64) u64 {
    var hasher = std.hash.Wyhash.init(0x6375_745f_706f_6f6c);
    hasher.update(std.mem.asBytes(&rhs));
    for (cut.indices, cut.values) |column, value| {
        hasher.update(std.mem.asBytes(&column));
        hasher.update(std.mem.asBytes(&value));
    }
    return hasher.final();
}

test "cut pool deduplicates reactivates ages and compacts rows" {
    var pool = try CutPool.init(std.testing.allocator, 3);
    defer pool.deinit(std.testing.allocator);
    const columns = [_]foundation.ColId{ try foundation.ColId.init(0), try foundation.ColId.init(2) };
    const values = [_]f64{ 1.0, 1.0 };
    const cut = matrix.SparseVectorView(foundation.ColId).initAssumeValid(3, &columns, &values);
    try std.testing.expectEqual(std.meta.Tag(AddResult).added, std.meta.activeTag(try pool.add(std.testing.allocator, cut, 1.0, 0.5)));
    try std.testing.expectEqual(std.meta.Tag(AddResult).existing, std.meta.activeTag(try pool.add(std.testing.allocator, cut, 1.0, 0.7)));
    try std.testing.expectEqual(@as(usize, 1), pool.len());
    try std.testing.expectEqual(@as(usize, 0), pool.ageAt(&.{ 0.0, 0.0, 0.0 }, 1e-6, 1));
    try std.testing.expectEqual(@as(usize, 1), pool.ageAt(&.{ 0.0, 0.0, 0.0 }, 1e-6, 1));
    try std.testing.expectEqual(@as(usize, 0), pool.activeCount());
    try pool.compact(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 0), pool.len());
}
