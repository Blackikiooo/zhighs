//! Structural row orderings for the sparse normal-equation factorization.
//!
//! Both orderings operate on the shared CSR adjacency of the normal
//! equations graph (see `graph.zig`) and return `old_to_new` row maps.
//! Reverse Cuthill-McKee reduces the factor envelope; minimum degree
//! greedily minimizes fill. The caller ranks candidates by exact symbolic
//! fill before applying one.

const std = @import("std");
const graph_mod = @import("graph.zig");
const Graph = graph_mod.Graph;

/// Reverse Cuthill-McKee: BFS from the minimum-degree node of each
/// component, visiting neighbors by increasing degree; each component's
/// visit order is reversed. O(nnz log Δ).
pub fn reverseCuthillMcKee(allocator: std.mem.Allocator, g: *const Graph) ![]usize {
    const n = g.n;
    const visited = try allocator.alloc(bool, n);
    defer allocator.free(visited);
    @memset(visited, false);
    const queue = try allocator.alloc(u32, n);
    defer allocator.free(queue);
    const neighbors = try allocator.alloc(u32, n);
    defer allocator.free(neighbors);
    const order = try allocator.alloc(u32, n);
    defer allocator.free(order);

    var order_len: usize = 0;
    while (order_len < n) {
        var seed: ?u32 = null;
        for (0..n) |row| {
            if (visited[row]) continue;
            if (seed == null or degreeLess(g, row, seed.?)) seed = @intCast(row);
        }
        const component_begin = order_len;
        var head: usize = 0;
        var tail: usize = 1;
        queue[0] = seed.?;
        visited[seed.?] = true;
        while (head < tail) : (head += 1) {
            const row = queue[head];
            order[order_len] = row;
            order_len += 1;

            var neighbor_len: usize = 0;
            for (g.neighbors[g.starts[row]..g.starts[row + 1]]) |candidate| {
                if (visited[candidate]) continue;
                visited[candidate] = true;
                neighbors[neighbor_len] = candidate;
                neighbor_len += 1;
            }
            std.mem.sort(u32, neighbors[0..neighbor_len], g, byDegreeThenId);
            for (neighbors[0..neighbor_len]) |candidate| {
                queue[tail] = candidate;
                tail += 1;
            }
        }
        std.mem.reverse(u32, order[component_begin..order_len]);
    }

    const old_to_new = try allocator.alloc(usize, n);
    for (order, 0..) |old, new| old_to_new[old] = new;
    return old_to_new;
}

fn degreeLess(g: *const Graph, left: usize, right: usize) bool {
    const dl = g.starts[left + 1] - g.starts[left];
    const dr = g.starts[right + 1] - g.starts[right];
    return dl < dr or (dl == dr and left < right);
}

fn byDegreeThenId(g: *const Graph, left: u32, right: u32) bool {
    return degreeLess(g, left, right);
}

/// Exact minimum-degree ordering on a quotient graph: pivots are drawn from
/// degree buckets and external degrees are recomputed exactly after every
/// elimination. Best fill quality, highest ordering cost.
pub fn minimumDegree(allocator: std.mem.Allocator, g: *const Graph) ![]usize {
    return minimumDegreeImpl(true, allocator, g);
}

/// Approximate minimum degree (AMD): same quotient-graph elimination, but
/// degrees are bounded estimates (Amestoy–Davis–Duff), avoiding the exact
/// union-marking recomputation. Near-minimum fill at a fraction of the cost.
pub fn approximateMinimumDegree(allocator: std.mem.Allocator, g: *const Graph) ![]usize {
    return minimumDegreeImpl(false, allocator, g);
}

/// Minimum-degree elimination on a quotient graph. Eliminated variables
/// become elements whose clique is stored once, so fill edges are implicit;
/// a variable's adjacency is its surviving original neighbors plus the
/// element cliques containing it. Elements adjacent to the pivot are always
/// fully absorbed, so live element cliques never list eliminated variables.
///
/// This replaces the previous dense bit-graph elimination game, which cost
/// O(n³/64) time and O(n²/8) bytes; here memory is O(nnz + fill) and time
/// is proportional to the clique sizes actually touched.
fn minimumDegreeImpl(comptime exact: bool, allocator: std.mem.Allocator, g: *const Graph) ![]usize {
    const n = g.n;
    if (n == 0) return allocator.alloc(usize, 0);

    const alive = try allocator.alloc(bool, n);
    defer allocator.free(alive);
    @memset(alive, true);
    // Element cliques, pooled: element e spans pool[estart[e]..][0..elen[e]].
    var pool: std.ArrayListUnmanaged(u32) = .empty;
    defer pool.deinit(allocator);
    const estart = try allocator.alloc(usize, n);
    defer allocator.free(estart);
    const elen = try allocator.alloc(usize, n);
    defer allocator.free(elen);
    const element_dead = try allocator.alloc(bool, n);
    defer allocator.free(element_dead);
    @memset(element_dead, false);
    // Element membership lists per variable, pooled per append.
    var elists = try allocator.alloc(std.ArrayListUnmanaged(u32), n);
    defer allocator.free(elists);
    for (elists) |*list| list.* = .empty;
    defer for (elists) |*list| list.deinit(allocator);

    const degree = try allocator.alloc(usize, n);
    defer allocator.free(degree);
    // Degree buckets: head[d] first variable with degree d, linked by next/prev.
    const head = try allocator.alloc(i64, n);
    defer allocator.free(head);
    @memset(head, -1);
    const bnext = try allocator.alloc(i64, n);
    defer allocator.free(bnext);
    const bprev = try allocator.alloc(i64, n);
    defer allocator.free(bprev);
    for (0..n) |i| {
        degree[i] = g.starts[i + 1] - g.starts[i];
        bucketInsert(head, bnext, bprev, i, degree[i]);
    }

    const mark = try allocator.alloc(usize, n);
    defer allocator.free(mark);
    @memset(mark, 0);
    var tag: usize = 0;
    var clique: std.ArrayListUnmanaged(u32) = .empty;
    defer clique.deinit(allocator);
    // AMD-only: per-element overlap with the current pivot clique, rebuilt
    // each elimination and reset via the touched list.
    const AmdState = struct { overlap: []usize, touched: std.ArrayListUnmanaged(u32) };
    var amd: if (exact) void else AmdState = if (exact) {} else .{
        .overlap = try allocator.alloc(usize, n),
        .touched = .empty,
    };
    defer if (!exact) {
        allocator.free(amd.overlap);
        amd.touched.deinit(allocator);
    };
    if (!exact) @memset(amd.overlap, 0);

    const old_to_new = try allocator.alloc(usize, n);
    errdefer allocator.free(old_to_new);
    var current_min: usize = 0;
    var dbg_clique: usize = 0;
    var dbg_deg: usize = 0;
    defer std.debug.print("amd n={} clique_walk={} deg_walk={}\n", .{ n, dbg_clique, dbg_deg });

    for (0..n) |step| {
        while (current_min < n and head[current_min] < 0) current_min += 1;
        const pivot: usize = @intCast(head[current_min]);
        bucketRemove(head, bnext, bprev, pivot, degree[pivot]);
        alive[pivot] = false;
        old_to_new[pivot] = step;

        // Clique of the new element: surviving original neighbors of the
        // pivot plus the members of its live element cliques.
        tag += 1;
        clique.clearRetainingCapacity();
        for (g.neighbors[g.starts[pivot]..g.starts[pivot + 1]]) |v| {
            if (alive[v] and mark[v] != tag) {
                mark[v] = tag;
                try clique.append(allocator, v);
            }
        }
        for (elists[pivot].items) |e| {
            if (element_dead[e]) continue;
            dbg_clique += elen[e];
            for (pool.items[estart[e]..][0..elen[e]]) |v| {
                if (alive[v] and mark[v] != tag) {
                    mark[v] = tag;
                    try clique.append(allocator, v);
                }
            }
        }

        // Every element adjacent to the pivot is covered by the new clique
        // and dies; live clique lists therefore never contain dead members.
        for (elists[pivot].items) |e| element_dead[e] = true;

        // The pivot becomes element `pivot`.
        estart[pivot] = pool.items.len;
        elen[pivot] = clique.items.len;
        try pool.appendSlice(allocator, clique.items);

        if (!exact) {
            // Overlap of each surviving element clique with L_pivot, used by
            // the degree estimates below. Element lists are compacted on the
            // way in: dead entries are dropped in place, so every later scan
            // of the same list touches live elements only. Each element is
            // appended once and removed at most once, keeping the total
            // compaction work proportional to the fill actually created —
            // without this the lists grow monotonically and the degree loop
            // degenerates into rescanning long-dead cliques.
            amd.touched.clearRetainingCapacity();
            for (clique.items) |v| {
                const list = &elists[v];
                var live: usize = 0;
                for (list.items) |e| {
                    if (element_dead[e]) continue;
                    list.items[live] = e;
                    live += 1;
                    if (amd.overlap[e] == 0) try amd.touched.append(allocator, e);
                    amd.overlap[e] += 1;
                }
                list.shrinkRetainingCapacity(live);
            }
        }

        // Join the neighbors to the new element and refresh degrees; the
        // minimum can only move to updated variables.
        var updated_min: usize = n;
        for (clique.items) |v| {
            if (exact) {
                // Same lazy compaction as the AMD pre-pass, done here so the
                // external-degree scan below walks live elements only.
                const list = &elists[v];
                var live: usize = 0;
                for (list.items) |e| {
                    if (!element_dead[e]) {
                        list.items[live] = e;
                        live += 1;
                    }
                }
                list.shrinkRetainingCapacity(live);
            }
            try elists[v].append(allocator, @intCast(pivot));
            const new_degree = if (exact) blk: {
                // Exact external degree: mark the union of surviving
                // original neighbors and live element cliques.
                tag += 1;
                var external: usize = 0;
                for (g.neighbors[g.starts[v]..g.starts[v + 1]]) |u| {
                    if (alive[u] and u != v and mark[u] != tag) {
                        mark[u] = tag;
                        external += 1;
                    }
                }
                for (elists[v].items) |e| {
                    for (pool.items[estart[e]..][0..elen[e]]) |u| {
                        if (alive[u] and u != v and mark[u] != tag) {
                            mark[u] = tag;
                            external += 1;
                        }
                    }
                }
                break :blk external;
            } else blk: {
                // AMD bound: original neighbors outside the pivot clique,
                // plus the rest of the pivot clique, plus the non-overlapping
                // remainder of every other live element clique containing v.
                var outside: usize = 0;
                for (g.neighbors[g.starts[v]..g.starts[v + 1]]) |u| {
                    if (alive[u] and mark[u] != tag) outside += 1;
                }
                var remainder: usize = 0;
                dbg_deg += elists[v].items.len;
                for (elists[v].items) |e| {
                    if (e == pivot) continue;
                    remainder += elen[e] - amd.overlap[e];
                }
                const remaining = n - step - 1;
                break :blk @min(outside + clique.items.len - 1 + remainder, remaining - 1);
            };
            bucketRemove(head, bnext, bprev, v, degree[v]);
            degree[v] = new_degree;
            bucketInsert(head, bnext, bprev, v, new_degree);
            updated_min = @min(updated_min, new_degree);
        }
        if (!exact) for (amd.touched.items) |e| {
            amd.overlap[e] = 0;
        };
        current_min = @min(current_min, updated_min);
    }
    return old_to_new;
}

fn bucketInsert(head: []i64, bnext: []i64, bprev: []i64, v: usize, d: usize) void {
    bprev[v] = -1;
    bnext[v] = head[d];
    if (head[d] >= 0) bprev[@intCast(head[d])] = @intCast(v);
    head[d] = @intCast(v);
}

fn bucketRemove(head: []i64, bnext: []i64, bprev: []i64, v: usize, d: usize) void {
    if (bprev[v] >= 0) bnext[@intCast(bprev[v])] = bnext[v] else head[d] = bnext[v];
    if (bnext[v] >= 0) bprev[@intCast(bnext[v])] = bprev[v];
}

test "minimum-degree orderings eliminate star fill and are valid permutations" {
    // Star normal-equations graph: column j holds rows {0, j}. Eliminating
    // the hub first would fill the whole graph; degree-based orderings must
    // take the leaves first and produce zero fill beyond the star edges.
    const a = std.testing.allocator;
    const n = 16;
    var col_starts: [n + 1]usize = undefined;
    var row_indices: [2 * (n - 1)]usize = undefined;
    for (0..n) |j| {
        col_starts[j] = if (j == 0) 0 else 2 * (j - 1);
        if (j > 0) {
            row_indices[2 * (j - 1)] = 0;
            row_indices[2 * (j - 1) + 1] = j;
        }
    }
    col_starts[n] = 2 * (n - 1);
    var g = try graph_mod.build(a, n, &col_starts, &row_indices);
    defer g.deinit();

    const natural_fill = try graph_mod.fillCount(a, &g, null);
    try std.testing.expect(natural_fill > n * n / 4);
    for ([_]*const fn (std.mem.Allocator, *const Graph) anyerror![]usize{
        minimumDegree,
        approximateMinimumDegree,
        reverseCuthillMcKee,
    }) |ordering| {
        const map = try ordering(a, &g);
        defer a.free(map);
        var seen = [_]bool{false} ** n;
        for (map) |new| {
            try std.testing.expect(new < n and !seen[new]);
            seen[new] = true;
        }
        // Hub last: leaves have degree 1, so the hub is eliminated last and
        // the factor keeps exactly one diagonal plus one off-diagonal per row.
        try std.testing.expectEqual(n * 2 - 1, try graph_mod.fillCount(a, &g, map));
    }
}