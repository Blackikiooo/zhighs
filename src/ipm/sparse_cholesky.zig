//! Row-oriented sparse Cholesky with a reusable symbolic fill pattern.
//!
//! The symbolic phase works on the shared adjacency graph (`graph.zig`):
//! the elimination tree and the exact fill pattern are computed with path
//! compression and per-column reaches — O(fill · α(n)) time and O(fill)
//! memory, replacing the previous dense bit-graph elimination game.
const std = @import("std");
const StandardForm = @import("standard_form.zig").StandardForm;
const graph_mod = @import("graph.zig");

pub const SparseCholesky = struct {
    allocator: std.mem.Allocator,
    row_starts: []usize,
    col_indices: []usize,
    dot_starts: []usize,
    dot_pairs: []u64,
    tail_pairs: []u64,
    tail_lengths: []u32,
    product_count: usize,
    indirect_count: usize,
    data: []f64,

    pub fn init(allocator: std.mem.Allocator, sf: *const StandardForm, g: *const graph_mod.Graph, order: *const graph_mod.Order) !SparseCholesky {
        const n = sf.num_rows;
        const parent = try graph_mod.eliminationTree(allocator, g, order);
        defer allocator.free(parent);
        var pattern = try graph_mod.symbolicPattern(allocator, g, order, parent);
        defer pattern.deinit();

        // The pattern is already the CSR lower triangle with ascending
        // columns and the diagonal last; widen the u32 columns to usize.
        const row_starts = try allocator.dupe(usize, pattern.starts);
        errdefer allocator.free(row_starts);
        const col_indices = try allocator.alloc(usize, row_starts[n]);
        errdefer allocator.free(col_indices);
        for (pattern.cols, 0..) |col, i| col_indices[i] = col;

        var pair_capacity: usize = 0;
        for (0..n) |row| {
            const begin = row_starts[row];
            const diagonal = row_starts[row + 1] - 1;
            for (begin..diagonal) |entry| {
                const col = col_indices[entry];
                const other_begin = row_starts[col];
                const other_diagonal = row_starts[col + 1] - 1;
                pair_capacity = try std.math.add(usize, pair_capacity, @min(entry - begin, other_diagonal - other_begin));
            }
        }

        const dot_starts = try allocator.alloc(usize, row_starts[n] + 1);
        errdefer allocator.free(dot_starts);
        const tail_pairs = try allocator.alloc(u64, row_starts[n]);
        errdefer allocator.free(tail_pairs);
        const tail_lengths = try allocator.alloc(u32, row_starts[n]);
        errdefer allocator.free(tail_lengths);
        @memset(tail_pairs, 0);
        @memset(tail_lengths, 0);
        if (row_starts[n] > std.math.maxInt(u32)) return error.OutOfMemory;
        const dot_pairs = try allocator.alloc(u64, pair_capacity);
        errdefer allocator.free(dot_pairs);
        var pair: usize = 0;
        var product_count: usize = 0;
        for (0..n) |row| {
            const begin = row_starts[row];
            const diagonal = row_starts[row + 1] - 1;
            for (begin..row_starts[row + 1]) |entry| {
                dot_starts[entry] = pair;
                if (entry == diagonal) continue;
                const col = col_indices[entry];
                var left = begin;
                var right = row_starts[col];
                const other_diagonal = row_starts[col + 1] - 1;
                var run_left: usize = 0;
                var run_right: usize = 0;
                var run_length: usize = 0;
                while (left < entry and right < other_diagonal) {
                    if (col_indices[left] == col_indices[right]) {
                        dot_pairs[pair] = (@as(u64, @intCast(left)) << 32) | @as(u32, @intCast(right));
                        pair += 1;
                        product_count = try std.math.add(usize, product_count, 1);
                        if (run_length != 0 and left == run_left + run_length and right == run_right + run_length) {
                            run_length += 1;
                        } else {
                            run_left = left;
                            run_right = right;
                            run_length = 1;
                        }
                        left += 1;
                        right += 1;
                    } else if (col_indices[left] < col_indices[right]) {
                        left += 1;
                    } else {
                        right += 1;
                    }
                }
                if (run_length >= 4) {
                    pair -= run_length;
                    tail_pairs[entry] = (@as(u64, @intCast(run_left)) << 32) | @as(u32, @intCast(run_right));
                    tail_lengths[entry] = @intCast(run_length);
                }
            }
        }
        dot_starts[row_starts[n]] = pair;
        const data = try allocator.alloc(f64, row_starts[n]);
        return .{ .allocator = allocator, .row_starts = row_starts, .col_indices = col_indices, .dot_starts = dot_starts, .dot_pairs = dot_pairs, .tail_pairs = tail_pairs, .tail_lengths = tail_lengths, .product_count = product_count, .indirect_count = pair, .data = data };
    }

    pub fn deinit(self: *SparseCholesky) void {
        self.allocator.free(self.row_starts);
        self.allocator.free(self.col_indices);
        self.allocator.free(self.dot_starts);
        self.allocator.free(self.dot_pairs);
        self.allocator.free(self.tail_pairs);
        self.allocator.free(self.tail_lengths);
        self.allocator.free(self.data);
    }

    pub fn productCount(self: *const SparseCholesky) usize {
        return self.product_count;
    }

    pub fn index(self: *const SparseCholesky, row: usize, col: usize) usize {
        var low = self.row_starts[row];
        var high = self.row_starts[row + 1];
        while (low < high) {
            const mid = low + (high - low) / 2;
            if (self.col_indices[mid] < col) low = mid + 1 else high = mid;
        }
        std.debug.assert(low < self.row_starts[row + 1] and self.col_indices[low] == col);
        return low;
    }

    /// Factor the SPD matrix whose lower triangle is given in the symbolic
    /// layout. Pivots below `pivot_floor` are lifted to the floor, which is
    /// equivalent to a diagonal shift applied only where needed (per-pivot
    /// regularization); pass 0 to reject non-positive pivots instead.
    /// Non-finite pivots always fail.
    pub fn factor(self: *SparseCholesky, values: []const f64, pivot_floor: f64) error{NotPositiveDefinite}!void {
        @memcpy(self.data, values);
        for (0..self.row_starts.len - 1) |row| {
            const begin = self.row_starts[row];
            const diagonal = self.row_starts[row + 1] - 1;
            for (begin..diagonal) |entry| {
                const col = self.col_indices[entry];
                const other_diagonal = self.row_starts[col + 1] - 1;
                var value = self.data[entry];
                for (self.dot_starts[entry]..self.dot_starts[entry + 1]) |pair| {
                    const indices = self.dot_pairs[pair];
                    value -= self.data[@intCast(indices >> 32)] * self.data[@intCast(indices & 0xffff_ffff)];
                }
                const tail_length: usize = self.tail_lengths[entry];
                if (tail_length != 0) {
                    const indices = self.tail_pairs[entry];
                    const left: usize = @intCast(indices >> 32);
                    const right: usize = @intCast(indices & 0xffff_ffff);
                    value -= dotProduct4(self.data[left..][0..tail_length], self.data[right..][0..tail_length]);
                }
                self.data[entry] = value / self.data[other_diagonal];
            }
            var pivot = self.data[diagonal];
            for (self.data[begin..diagonal]) |value| pivot -= value * value;
            if (!std.math.isFinite(pivot) or (pivot_floor <= 0.0 and pivot <= 0.0)) return error.NotPositiveDefinite;
            self.data[diagonal] = @sqrt(@max(pivot, pivot_floor));
        }
    }

    pub fn solve(self: *const SparseCholesky, rhs: []f64) void {
        for (0..self.row_starts.len - 1) |row| {
            const begin = self.row_starts[row];
            const diagonal = self.row_starts[row + 1] - 1;
            var value = rhs[row];
            for (begin..diagonal) |entry| value -= self.data[entry] * rhs[self.col_indices[entry]];
            rhs[row] = value / self.data[diagonal];
        }
        var row = self.row_starts.len - 1;
        while (row > 0) {
            row -= 1;
            const begin = self.row_starts[row];
            const diagonal = self.row_starts[row + 1] - 1;
            rhs[row] /= self.data[diagonal];
            for (begin..diagonal) |entry| rhs[self.col_indices[entry]] -= self.data[entry] * rhs[row];
        }
    }
};

inline fn dotProduct4(left: []const f64, right: []const f64) f64 {
    var sum: f64 = 0.0;
    var index: usize = 0;
    while (index + 16 <= left.len) : (index += 16) {
        const l0: @Vector(4, f64) = left[index..][0..4].*;
        const l1: @Vector(4, f64) = left[index + 4 ..][0..4].*;
        const l2: @Vector(4, f64) = left[index + 8 ..][0..4].*;
        const l3: @Vector(4, f64) = left[index + 12 ..][0..4].*;
        const r0: @Vector(4, f64) = right[index..][0..4].*;
        const r1: @Vector(4, f64) = right[index + 4 ..][0..4].*;
        const r2: @Vector(4, f64) = right[index + 8 ..][0..4].*;
        const r3: @Vector(4, f64) = right[index + 12 ..][0..4].*;
        sum += @reduce(.Add, (l0 * r0 + l1 * r1) + (l2 * r2 + l3 * r3));
    }
    while (index + 4 <= left.len) : (index += 4) {
        const left_vector: @Vector(4, f64) = left[index..][0..4].*;
        const right_vector: @Vector(4, f64) = right[index..][0..4].*;
        sum += @reduce(.Add, left_vector * right_vector);
    }
    while (index < left.len) : (index += 1) sum += left[index] * right[index];
    return sum;
}