//! Shared adjacency of the normal-equations sparsity graph G(A·Aᵀ), built
//! once per solve and reused by the row ordering, the density heuristic, and
//! the symbolic Cholesky factorization.
//!
//! Previously each consumer expanded the column cliques into its own dense
//! bit graph (O(n²/8) bytes each, three expansions per solve). This module
//! builds one CSR adjacency (O(edges) bytes) and provides elimination-tree
//! based symbolic analysis (Liu's algorithm with path compression, plus the
//! per-column reach that yields the exact Cholesky fill pattern) so no dense
//! graph is ever materialized.

const std = @import("std");
const Allocator = std.mem.Allocator;

pub const Graph = struct {
    allocator: Allocator,
    n: usize,
    /// CSR adjacency, sorted ascending per node: neighbors of node i in
    /// neighbors[starts[i]..starts[i+1]].
    starts: []usize,
    neighbors: []u32,

    pub fn deinit(self: *Graph) void {
        self.allocator.free(self.starts);
        self.allocator.free(self.neighbors);
    }

    pub fn edgeCount(self: *const Graph) usize {
        return self.starts[self.n] / 2;
    }

    /// Fraction of present off-diagonal pairs, in [0, 1].
    pub fn density(self: *const Graph) f64 {
        if (self.n < 2) return 0.0;
        const present: f64 = @floatFromInt(self.starts[self.n]);
        const possible: f64 = @floatFromInt(self.n *| (self.n - 1));
        return @min(present / possible, 1.0);
    }
};

/// Build the sparsity graph of A·Aᵀ from the CSC matrix: two rows are
/// adjacent iff they share a column. Cost is proportional to the clique
/// expansion Σ_j nnz_j², paid exactly once per solve.
pub fn build(
    allocator: Allocator,
    num_rows: usize,
    col_starts: []const usize,
    row_indices: []const usize,
) !Graph {
    const n = num_rows;
    const marker = try allocator.alloc(usize, n);
    defer allocator.free(marker);
    const degree = try allocator.alloc(usize, n);
    defer allocator.free(degree);
    @memset(degree, 0);
    @memset(marker, std.math.maxInt(usize));

    const num_cols = col_starts.len - 1;
    // Pass 1: count neighbors per node. marker[q] == r deduplicates q in
    // r's list across all columns containing both.
    for (0..num_cols) |j| {
        const rows = row_indices[col_starts[j]..col_starts[j + 1]];
        for (rows) |r| {
            for (rows) |q| {
                if (q != r and marker[q] != r) {
                    marker[q] = r;
                    degree[r] += 1;
                }
            }
        }
    }

    const starts = try allocator.alloc(usize, n + 1);
    errdefer allocator.free(starts);
    starts[0] = 0;
    for (degree, 0..) |d, i| starts[i + 1] = try std.math.add(usize, starts[i], d);
    const neighbors = try allocator.alloc(u32, starts[n]);
    errdefer allocator.free(neighbors);

    // Pass 2: fill adjacency, then sort each list for deterministic walks.
    @memset(marker, std.math.maxInt(usize));
    const cursor = try allocator.dupe(usize, starts[0..n]);
    defer allocator.free(cursor);
    for (0..num_cols) |j| {
        const rows = row_indices[col_starts[j]..col_starts[j + 1]];
        for (rows) |r| {
            for (rows) |q| {
                if (q != r and marker[q] != r) {
                    marker[q] = r;
                    neighbors[cursor[r]] = @intCast(q);
                    cursor[r] += 1;
                }
            }
        }
    }
    for (0..n) |i| std.mem.sort(u32, neighbors[starts[i]..starts[i + 1]], {}, std.sort.asc(u32));

    return .{ .allocator = allocator, .n = n, .starts = starts, .neighbors = neighbors };
}

/// Elimination order: `perm[node] = position`, `iperm[position] = node`.
/// Wraps the natural order or an ordering returned by `row_ordering`.
pub const Order = struct {
    allocator: Allocator,
    perm: []u32,
    iperm: []u32,

    pub fn init(allocator: Allocator, n: usize, old_to_new: ?[]const usize) !Order {
        const perm = try allocator.alloc(u32, n);
        errdefer allocator.free(perm);
        const iperm = try allocator.alloc(u32, n);
        errdefer allocator.free(iperm);
        for (0..n) |i| {
            const pos = if (old_to_new) |map| map[i] else i;
            perm[i] = @intCast(pos);
            iperm[pos] = @intCast(i);
        }
        return .{ .allocator = allocator, .perm = perm, .iperm = iperm };
    }

    pub fn deinit(self: *Order) void {
        self.allocator.free(self.perm);
        self.allocator.free(self.iperm);
    }
};

/// Elimination tree under the given order. `parent[p]` is the parent of
/// position p — the first later position whose elimination touches p — or
/// n for component roots. Liu's algorithm with path compression; runs in
/// O(edges · α(n)).
pub fn eliminationTree(allocator: Allocator, g: *const Graph, order: *const Order) ![]usize {
    const n = g.n;
    const parent = try allocator.alloc(usize, n);
    errdefer allocator.free(parent);
    const ancestor = try allocator.alloc(usize, n);
    defer allocator.free(ancestor);
    @memset(ancestor, n);
    for (0..n) |k| {
        parent[k] = n;
        const u = order.iperm[k];
        for (g.neighbors[g.starts[u]..g.starts[u + 1]]) |v| {
            const i: usize = order.perm[v];
            if (i >= k) continue;
            // Walk from i to the root of its partial tree, compressing the
            // path onto k; attach the root under k.
            var r = i;
            while (ancestor[r] != n and ancestor[r] != k) {
                const next = ancestor[r];
                ancestor[r] = k;
                r = next;
            }
            if (ancestor[r] == n) {
                parent[r] = k;
                ancestor[r] = k;
            }
        }
    }
    return parent;
}

/// Pattern of row k of L (positions ≤ k, unsorted, including the diagonal)
/// written to scratch; returns the length. `mark` is a zero-initialized
/// timestamp array the caller reuses across rows. Each lower-position
/// neighbor of k seeds an etree walk that stops at already-reached nodes;
/// the visited nodes are exactly the filled-graph neighbors of k.
fn rowReach(
    g: *const Graph,
    order: *const Order,
    parent: []const usize,
    mark: []usize,
    k: usize,
    scratch: []u32,
) usize {
    const tag = k + 1;
    var len: usize = 0;
    mark[k] = tag;
    scratch[len] = @intCast(k);
    len += 1;
    const u = order.iperm[k];
    for (g.neighbors[g.starts[u]..g.starts[u + 1]]) |v| {
        var w: usize = order.perm[v];
        if (w >= k) continue;
        while (mark[w] != tag) {
            mark[w] = tag;
            scratch[len] = @intCast(w);
            len += 1;
            w = parent[w];
        }
    }
    return len;
}

/// Exact Cholesky fill (nonzeros of L including the diagonal) under the
/// given elimination order. Used to rank ordering candidates; cost is
/// O(fill · α(n)) per candidate.
pub fn fillCount(allocator: Allocator, g: *const Graph, old_to_new: ?[]const usize) !usize {
    const n = g.n;
    var order = try Order.init(allocator, n, old_to_new);
    defer order.deinit();
    const parent = try eliminationTree(allocator, g, &order);
    defer allocator.free(parent);
    const mark = try allocator.alloc(usize, n);
    defer allocator.free(mark);
    @memset(mark, 0);
    const scratch = try allocator.alloc(u32, n);
    defer allocator.free(scratch);
    var fill: usize = 0;
    for (0..n) |k| fill = try std.math.add(usize, fill, rowReach(g, &order, parent, mark, k, scratch));
    return fill;
}

/// CSR row patterns of the Cholesky factor L under the given order:
/// cols[starts[r]..starts[r+1]] are the column positions of row r, sorted
/// ascending — the diagonal is the last entry of each row.
pub const Pattern = struct {
    allocator: Allocator,
    starts: []usize,
    cols: []u32,

    pub fn deinit(self: *Pattern) void {
        self.allocator.free(self.starts);
        self.allocator.free(self.cols);
    }
};

pub fn symbolicPattern(
    allocator: Allocator,
    g: *const Graph,
    order: *const Order,
    parent: []const usize,
) !Pattern {
    const n = g.n;
    const mark = try allocator.alloc(usize, n);
    defer allocator.free(mark);
    @memset(mark, 0);
    const scratch = try allocator.alloc(u32, n);
    defer allocator.free(scratch);

    const starts = try allocator.alloc(usize, n + 1);
    errdefer allocator.free(starts);
    starts[0] = 0;
    for (0..n) |k| starts[k + 1] = try std.math.add(usize, starts[k], rowReach(g, order, parent, mark, k, scratch));

    const cols = try allocator.alloc(u32, starts[n]);
    errdefer allocator.free(cols);
    @memset(mark, 0);
    for (0..n) |k| {
        const len = rowReach(g, order, parent, mark, k, scratch);
        std.debug.assert(len == starts[k + 1] - starts[k]);
        @memcpy(cols[starts[k]..starts[k + 1]], scratch[0..len]);
        std.mem.sort(u32, cols[starts[k]..starts[k + 1]], {}, std.sort.asc(u32));
    }
    return .{ .allocator = allocator, .starts = starts, .cols = cols };
}