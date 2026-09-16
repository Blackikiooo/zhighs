//! Row-profile Cholesky. Symbolic storage is computed once from A's CSC
//! pattern; elimination without pivoting cannot fill outside this envelope.
//! Both triangular solves traverse contiguous rows, including L-transpose.
const std = @import("std");
const StandardForm = @import("standard_form.zig").StandardForm;

pub const ProfileCholesky = struct {
    allocator: std.mem.Allocator,
    first: []usize,
    offsets: []usize,
    data: []f64,

    pub fn init(allocator: std.mem.Allocator, sf: *const StandardForm) !ProfileCholesky {
        const first = try allocator.alloc(usize, sf.num_rows);
        errdefer allocator.free(first);
        for (first, 0..) |*value, i| value.* = i;
        for (0..sf.num_cols) |j| {
            const rows = sf.row_indices[sf.col_starts[j]..sf.col_starts[j + 1]];
            var minimum = sf.num_rows;
            for (rows) |row| minimum = @min(minimum, row);
            for (rows) |row| first[row] = @min(first[row], minimum);
        }
        const offsets = try allocator.alloc(usize, sf.num_rows + 1);
        errdefer allocator.free(offsets);
        offsets[0] = 0;
        for (first, 0..) |begin, i| offsets[i + 1] = try std.math.add(usize, offsets[i], i - begin + 1);
        const data = try allocator.alloc(f64, offsets[sf.num_rows]);
        return .{ .allocator = allocator, .first = first, .offsets = offsets, .data = data };
    }

    pub fn deinit(self: *ProfileCholesky) void {
        self.allocator.free(self.first);
        self.allocator.free(self.offsets);
        self.allocator.free(self.data);
    }

    pub fn index(self: *const ProfileCholesky, row: usize, col: usize) usize {
        std.debug.assert(col >= self.first[row] and col <= row);
        return self.offsets[row] + col - self.first[row];
    }

    pub fn factor(self: *ProfileCholesky, values: []const f64) error{NotPositiveDefinite}!void {
        @memcpy(self.data, values);
        for (self.first, 0..) |begin, i| {
            const row = self.data[self.offsets[i]..self.offsets[i + 1]];
            for (begin..i) |j| {
                var value = row[j - begin];
                const other_begin = self.first[j];
                const other = self.data[self.offsets[j]..self.offsets[j + 1]];
                for (@max(begin, other_begin)..j) |k| value -= row[k - begin] * other[k - other_begin];
                row[j - begin] = value / other[j - other_begin];
            }
            var pivot = row[i - begin];
            for (row[0 .. i - begin]) |value| pivot -= value * value;
            if (!std.math.isFinite(pivot) or pivot <= 0) return error.NotPositiveDefinite;
            row[i - begin] = @sqrt(pivot);
        }
    }

    pub fn solve(self: *const ProfileCholesky, rhs: []f64) void {
        std.debug.assert(rhs.len == self.first.len);
        for (self.first, 0..) |begin, i| {
            const row = self.data[self.offsets[i]..self.offsets[i + 1]];
            var value = rhs[i];
            for (begin..i) |j| value -= row[j - begin] * rhs[j];
            rhs[i] = value / row[i - begin];
        }
        var i = rhs.len;
        while (i > 0) {
            i -= 1;
            const begin = self.first[i];
            const row = self.data[self.offsets[i]..self.offsets[i + 1]];
            rhs[i] /= row[i - begin];
            for (begin..i) |j| rhs[j] -= row[j - begin] * rhs[i];
        }
    }
};
