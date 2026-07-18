//! Reusable initial symbolic pivot analysis for sparse LU.
//!
//! The planner consumes compact basis CSC and builds only a row-to-entry
//! companion index. Numerical values remain single-copy in CSC. Active row and
//! column counts are maintained incrementally; columns live in intrusive count
//! buckets and singleton rows enter a monotonic queue. This is the DOD state
//! needed by a later packed numerical L/U kernel without coupling pivot policy
//! to simplex orchestration. The planner intentionally stops after selecting
//! the first non-singleton pivot: fill-in and updated numerical values must be
//! committed before another threshold Markowitz choice is valid.

const std = @import("std");
const foundation = @import("foundation");
const sparse_basis = @import("sparse_basis.zig");

const Offset = foundation.HUInt;
const none = std.math.maxInt(u32);

pub const SymbolicError = error{
    DimensionTooLarge,
    InvalidBasis,
    Singular,
    OutOfMemory,
};

pub const PivotKind = enum(u8) { row_singleton, column_singleton, markowitz };

/// Borrowed pivot order. Slices remain valid until the workspace is reused.
pub const SymbolicPlanView = struct {
    pivot_rows: []const u32,
    pivot_columns: []const u32,
    pivot_kinds: []const PivotKind,
    singleton_pivots: usize,
    markowitz_pivots: usize,
    remaining_dimension: usize,
    active_rows: []const bool,
    active_columns: []const bool,
    row_counts: []const u32,
    column_counts: []const u32,
    row_starts: []const Offset,
    row_entries: []const u32,
    entry_columns: []const u32,

    pub inline fn dimension(self: SymbolicPlanView) usize {
        return self.pivot_rows.len;
    }
};

/// Retaining SoA workspace for singleton elimination and threshold Markowitz.
pub const SymbolicWorkspace = struct {
    allocator: std.mem.Allocator,
    dimension_capacity: usize = 0,
    entry_capacity: usize = 0,

    row_starts: []Offset = &.{},
    row_entries: []u32 = &.{},
    entry_columns: []u32 = &.{},
    row_counts: []u32 = &.{},
    column_counts: []u32 = &.{},
    row_active: []bool = &.{},
    column_active: []bool = &.{},

    bucket_first: []u32 = &.{},
    bucket_next: []u32 = &.{},
    bucket_previous: []u32 = &.{},
    singleton_rows: []u32 = &.{},

    pivot_rows: []u32 = &.{},
    pivot_columns: []u32 = &.{},
    pivot_kinds: []PivotKind = &.{},
    // Borrowed only during `plan`; fields avoid threading a wide view through
    // the hottest elimination helper. They are assigned before the pivot loop.
    current_starts: []const Offset = &.{},
    current_rows: []const foundation.RowId = &.{},

    pub fn init(allocator: std.mem.Allocator) SymbolicWorkspace {
        return .{ .allocator = allocator };
    }

    pub fn deinit(self: *SymbolicWorkspace) void {
        self.freeDimensionArrays();
        self.allocator.free(self.row_entries);
        self.allocator.free(self.entry_columns);
        self.* = .{ .allocator = self.allocator };
    }

    /// Requested bytes retained for the next reinversion. Allocator metadata
    /// and alignment padding outside the returned slices are intentionally not
    /// guessed; process RSS is measured separately by dataset runners.
    pub fn retainedBytes(self: *const SymbolicWorkspace) usize {
        const n = self.dimension_capacity;
        return (n + 1) * @sizeOf(Offset) +
            self.entry_capacity * (@sizeOf(u32) + @sizeOf(u32)) +
            n * (@sizeOf(u32) * 7 + @sizeOf(bool) * 2 + @sizeOf(PivotKind)) +
            (n + 1) * @sizeOf(u32);
    }

    /// Eliminate the fill-free singleton prefix and select at most one kernel
    /// pivot by threshold Markowitz. `pivot_threshold` is relative to the
    /// largest active entry in each candidate column and must be in `(0, 1]`.
    /// The numerical factorizer must apply that kernel pivot and its fill-in
    /// before asking for another choice.
    pub fn plan(self: *SymbolicWorkspace, basis: sparse_basis.SparseBasisView, pivot_threshold: f64) SymbolicError!SymbolicPlanView {
        return self.planImpl(basis, pivot_threshold, false);
    }

    /// Return only the fill-free singleton prefix. The active masks and row
    /// companion describe the reduced kernel without materializing it.
    pub fn planSingletonPrefix(self: *SymbolicWorkspace, basis: sparse_basis.SparseBasisView) SymbolicError!SymbolicPlanView {
        return self.planImpl(basis, 0.1, true);
    }

    fn planImpl(self: *SymbolicWorkspace, basis: sparse_basis.SparseBasisView, pivot_threshold: f64, comptime singleton_only: bool) SymbolicError!SymbolicPlanView {
        const n = basis.dimension;
        if (basis.starts.len != n + 1 or basis.rows.len != basis.values.len or
            !std.math.isFinite(pivot_threshold) or pivot_threshold <= 0.0 or pivot_threshold > 1.0 or
            n > std.math.maxInt(u32) or basis.nnz() > std.math.maxInt(u32) or basis.nnz() > std.math.maxInt(Offset))
            return error.InvalidBasis;
        try self.ensureCapacity(n, basis.nnz());
        self.current_starts = basis.starts;
        self.current_rows = basis.rows;
        try self.initialize(basis);

        var queue_head: usize = 0;
        var queue_tail: usize = 0;
        for (self.row_counts[0..n], 0..) |count, row| {
            if (count == 0) return error.Singular;
            if (count == 1) {
                self.singleton_rows[queue_tail] = @intCast(row);
                queue_tail += 1;
            }
        }

        var pivot_count: usize = 0;
        var singleton_count: usize = 0;
        var markowitz_count: usize = 0;
        while (pivot_count < n) {
            var selected_row: ?u32 = null;
            var selected_column: ?u32 = null;
            var selected_kind: PivotKind = .markowitz;

            while (queue_head < queue_tail) {
                const row = self.singleton_rows[queue_head];
                queue_head += 1;
                if (!self.row_active[row] or self.row_counts[row] != 1) continue;
                const column = self.onlyActiveColumn(row) orelse return error.Singular;
                selected_row = row;
                selected_column = column;
                selected_kind = .row_singleton;
                break;
            }

            if (selected_row == null) {
                if (singleton_only) {
                    const column = self.bucket_first[1];
                    if (column == none) break;
                    const begin: usize = @intCast(basis.starts[column]);
                    const end: usize = @intCast(basis.starts[column + 1]);
                    for (basis.rows[begin..end]) |row_id| if (self.row_active[row_id.toUsize()]) {
                        selected_row = @intCast(row_id.toUsize());
                        selected_column = column;
                        selected_kind = .column_singleton;
                        break;
                    };
                    if (selected_row == null) return error.Singular;
                } else {
                    const choice = self.chooseMarkowitz(basis, pivot_threshold) orelse return error.Singular;
                    selected_row = choice.row;
                    selected_column = choice.column;
                    selected_kind = if (self.column_counts[choice.column] == 1) .column_singleton else .markowitz;
                }
            }

            const row = selected_row.?;
            const column = selected_column.?;
            self.pivot_rows[pivot_count] = row;
            self.pivot_columns[pivot_count] = column;
            self.pivot_kinds[pivot_count] = selected_kind;
            if (selected_kind == .markowitz) markowitz_count += 1 else singleton_count += 1;
            pivot_count += 1;
            if (selected_kind == .markowitz) break;
            self.eliminate(row, column, &queue_tail);
        }

        return .{
            .pivot_rows = self.pivot_rows[0..pivot_count],
            .pivot_columns = self.pivot_columns[0..pivot_count],
            .pivot_kinds = self.pivot_kinds[0..pivot_count],
            .singleton_pivots = singleton_count,
            .markowitz_pivots = markowitz_count,
            .remaining_dimension = n - pivot_count,
            .active_rows = self.row_active[0..n],
            .active_columns = self.column_active[0..n],
            .row_counts = self.row_counts[0..n],
            .column_counts = self.column_counts[0..n],
            .row_starts = self.row_starts[0 .. n + 1],
            .row_entries = self.row_entries[0..basis.nnz()],
            .entry_columns = self.entry_columns[0..basis.nnz()],
        };
    }

    fn initialize(self: *SymbolicWorkspace, basis: sparse_basis.SparseBasisView) SymbolicError!void {
        const n = basis.dimension;
        @memset(self.row_counts[0..n], 0);
        @memset(self.column_counts[0..n], 0);
        @memset(self.row_active[0..n], true);
        @memset(self.column_active[0..n], true);
        @memset(self.bucket_first[0 .. n + 1], none);

        for (0..n) |column| {
            const begin: usize = @intCast(basis.starts[column]);
            const end: usize = @intCast(basis.starts[column + 1]);
            if (begin > end or end > basis.nnz() or end - begin > n) return error.InvalidBasis;
            self.column_counts[column] = @intCast(end - begin);
            var previous_row: ?usize = null;
            for (begin..end) |entry| {
                const row = basis.rows[entry].toUsize();
                if (row >= n or !std.math.isFinite(basis.values[entry]) or basis.values[entry] == 0.0)
                    return error.InvalidBasis;
                if (previous_row) |previous| if (row <= previous) return error.InvalidBasis;
                previous_row = row;
                self.row_counts[row] += 1;
                self.entry_columns[entry] = @intCast(column);
            }
        }

        var running: usize = 0;
        for (self.row_counts[0..n], 0..) |count, row| {
            self.row_starts[row] = @intCast(running);
            running += count;
        }
        self.row_starts[n] = @intCast(running);
        // Reuse row_counts as insertion cursors, then restore the counts from
        // adjacent row starts. No fourth nnz-sized workspace is needed.
        @memset(self.row_counts[0..n], 0);
        for (basis.rows, 0..) |row_id, entry| {
            const row = row_id.toUsize();
            const output: usize = @intCast(self.row_starts[row] + self.row_counts[row]);
            self.row_entries[output] = @intCast(entry);
            self.row_counts[row] += 1;
        }
        for (0..n) |column| {
            const count = self.column_counts[column];
            if (count == 0) return error.Singular;
            self.bucketInsert(@intCast(column), count);
        }
    }

    const Choice = struct { row: u32, column: u32 };

    fn chooseMarkowitz(self: *SymbolicWorkspace, basis: sparse_basis.SparseBasisView, threshold: f64) ?Choice {
        const n = basis.dimension;
        var best: ?Choice = null;
        var best_merit: u64 = std.math.maxInt(u64);
        var column_count: usize = 1;
        while (column_count <= n) : (column_count += 1) {
            // Every non-singleton active row contributes at least one to the
            // row factor. Once the column factor cannot improve the best merit,
            // denser buckets cannot improve it either.
            if (best != null and column_count - 1 > best_merit) break;
            var column = self.bucket_first[column_count];
            while (column != none) : (column = self.bucket_next[column]) {
                var maximum: f64 = 0.0;
                const begin: usize = @intCast(basis.starts[column]);
                const end: usize = @intCast(basis.starts[column + 1]);
                for (begin..end) |entry| {
                    const row = basis.rows[entry].toUsize();
                    if (self.row_active[row]) maximum = @max(maximum, @abs(basis.values[entry]));
                }
                if (maximum == 0.0 or !std.math.isFinite(maximum)) continue;
                const minimum_pivot = threshold * maximum;
                for (begin..end) |entry| {
                    const row = basis.rows[entry].toUsize();
                    if (!self.row_active[row] or @abs(basis.values[entry]) < minimum_pivot) continue;
                    const merit = @as(u64, self.row_counts[row] - 1) * @as(u64, self.column_counts[column] - 1);
                    if (merit < best_merit or (merit == best_merit and isEarlier(row, column, best))) {
                        best_merit = merit;
                        best = .{ .row = @intCast(row), .column = column };
                    }
                }
            }
        }
        return best;
    }

    fn isEarlier(row: usize, column: u32, current: ?Choice) bool {
        const old = current orelse return true;
        return column < old.column or (column == old.column and row < old.row);
    }

    fn onlyActiveColumn(self: *const SymbolicWorkspace, row: u32) ?u32 {
        const begin: usize = @intCast(self.row_starts[row]);
        const end: usize = @intCast(self.row_starts[row + 1]);
        for (self.row_entries[begin..end]) |entry| {
            const column = self.entry_columns[entry];
            if (self.column_active[column]) return column;
        }
        return null;
    }

    fn eliminate(self: *SymbolicWorkspace, pivot_row: u32, pivot_column: u32, queue_tail: *usize) void {
        self.row_active[pivot_row] = false;
        self.column_active[pivot_column] = false;
        self.bucketRemove(pivot_column, self.column_counts[pivot_column]);

        // Removing the pivot row shortens every other active column touching
        // it. Move each column between count buckets in O(1).
        const row_begin: usize = @intCast(self.row_starts[pivot_row]);
        const row_end: usize = @intCast(self.row_starts[pivot_row + 1]);
        for (self.row_entries[row_begin..row_end]) |entry| {
            const column = self.entry_columns[entry];
            if (!self.column_active[column]) continue;
            const old_count = self.column_counts[column];
            self.bucketRemove(column, old_count);
            self.column_counts[column] = old_count - 1;
            self.bucketInsert(column, old_count - 1);
        }

        // Removing the pivot column shortens every other active row. Each row
        // can become singleton only once, so the fixed n-entry queue suffices.
        const column_begin: usize = @intCast(self.current_starts[pivot_column]);
        const column_end: usize = @intCast(self.current_starts[pivot_column + 1]);
        for (self.current_rows[column_begin..column_end]) |row_id| {
            const row: u32 = @intCast(row_id.toUsize());
            if (!self.row_active[row]) continue;
            self.row_counts[row] -= 1;
            if (self.row_counts[row] == 1) {
                self.singleton_rows[queue_tail.*] = row;
                queue_tail.* += 1;
            }
        }
    }

    fn bucketInsert(self: *SymbolicWorkspace, column: u32, count: u32) void {
        const first = self.bucket_first[count];
        self.bucket_previous[column] = none;
        self.bucket_next[column] = first;
        if (first != none) self.bucket_previous[first] = column;
        self.bucket_first[count] = column;
    }

    fn bucketRemove(self: *SymbolicWorkspace, column: u32, count: u32) void {
        const previous = self.bucket_previous[column];
        const next = self.bucket_next[column];
        if (previous == none) self.bucket_first[count] = next else self.bucket_next[previous] = next;
        if (next != none) self.bucket_previous[next] = previous;
        self.bucket_previous[column] = none;
        self.bucket_next[column] = none;
    }

    fn ensureCapacity(self: *SymbolicWorkspace, n: usize, nnz: usize) SymbolicError!void {
        if (n > self.dimension_capacity) {
            const capacity = growCapacity(self.dimension_capacity, n) catch return error.DimensionTooLarge;
            // Sequential realloc is failure-safe: every successful growth is
            // immediately owned by the workspace, while a failing realloc
            // leaves that particular old slice valid for `deinit` or retry.
            self.row_starts = self.allocator.realloc(self.row_starts, capacity + 1) catch return error.OutOfMemory;
            self.row_counts = self.allocator.realloc(self.row_counts, capacity) catch return error.OutOfMemory;
            self.column_counts = self.allocator.realloc(self.column_counts, capacity) catch return error.OutOfMemory;
            self.row_active = self.allocator.realloc(self.row_active, capacity) catch return error.OutOfMemory;
            self.column_active = self.allocator.realloc(self.column_active, capacity) catch return error.OutOfMemory;
            self.bucket_first = self.allocator.realloc(self.bucket_first, capacity + 1) catch return error.OutOfMemory;
            self.bucket_next = self.allocator.realloc(self.bucket_next, capacity) catch return error.OutOfMemory;
            self.bucket_previous = self.allocator.realloc(self.bucket_previous, capacity) catch return error.OutOfMemory;
            self.singleton_rows = self.allocator.realloc(self.singleton_rows, capacity) catch return error.OutOfMemory;
            self.pivot_rows = self.allocator.realloc(self.pivot_rows, capacity) catch return error.OutOfMemory;
            self.pivot_columns = self.allocator.realloc(self.pivot_columns, capacity) catch return error.OutOfMemory;
            self.pivot_kinds = self.allocator.realloc(self.pivot_kinds, capacity) catch return error.OutOfMemory;
            self.dimension_capacity = capacity;
        }
        if (nnz > self.entry_capacity) {
            const capacity = growCapacity(self.entry_capacity, nnz) catch return error.DimensionTooLarge;
            self.row_entries = self.allocator.realloc(self.row_entries, capacity) catch return error.OutOfMemory;
            self.entry_columns = self.allocator.realloc(self.entry_columns, capacity) catch return error.OutOfMemory;
            self.entry_capacity = capacity;
        }
    }

    fn freeDimensionArrays(self: *SymbolicWorkspace) void {
        self.allocator.free(self.row_starts);
        self.allocator.free(self.row_counts);
        self.allocator.free(self.column_counts);
        self.allocator.free(self.row_active);
        self.allocator.free(self.column_active);
        self.allocator.free(self.bucket_first);
        self.allocator.free(self.bucket_next);
        self.allocator.free(self.bucket_previous);
        self.allocator.free(self.singleton_rows);
        self.allocator.free(self.pivot_rows);
        self.allocator.free(self.pivot_columns);
        self.allocator.free(self.pivot_kinds);
        self.row_starts = &.{};
        self.row_counts = &.{};
        self.column_counts = &.{};
        self.row_active = &.{};
        self.column_active = &.{};
        self.bucket_first = &.{};
        self.bucket_next = &.{};
        self.bucket_previous = &.{};
        self.singleton_rows = &.{};
        self.pivot_rows = &.{};
        self.pivot_columns = &.{};
        self.pivot_kinds = &.{};
        self.dimension_capacity = 0;
    }
};

fn growCapacity(current: usize, required: usize) error{Overflow}!usize {
    var capacity = @max(current, 8);
    while (capacity < required) capacity = std.math.add(usize, capacity, capacity / 2 + 8) catch return error.Overflow;
    return capacity;
}

test "symbolic planner eliminates triangular basis as singletons" {
    const RowId = foundation.RowId;
    const basis = sparse_basis.SparseBasisView{
        .dimension = 3,
        .starts = &[_]Offset{ 0, 2, 4, 5 },
        .rows = &[_]RowId{
            RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(1),
            RowId.fromUsizeAssumeValid(1), RowId.fromUsizeAssumeValid(2),
            RowId.fromUsizeAssumeValid(2),
        },
        .values = &[_]f64{ 4, 1, 3, 1, 2 },
    };
    var workspace = SymbolicWorkspace.init(std.testing.allocator);
    defer workspace.deinit();
    const plan_view = try workspace.plan(basis, 0.1);
    try std.testing.expectEqual(@as(usize, 3), plan_view.singleton_pivots);
    try std.testing.expectEqual(@as(usize, 0), plan_view.markowitz_pivots);
    try std.testing.expectEqual(@as(usize, 0), plan_view.remaining_dimension);
}

test "threshold Markowitz rejects a tiny low-merit pivot" {
    const RowId = foundation.RowId;
    const basis = sparse_basis.SparseBasisView{
        .dimension = 3,
        .starts = &[_]Offset{ 0, 2, 4, 7 },
        .rows = &[_]RowId{
            RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(1),
            RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(2),
            RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(1),
            RowId.fromUsizeAssumeValid(2),
        },
        .values = &[_]f64{ 1e-12, 1, 2, 1, 1, 1, 3 },
    };
    var workspace = SymbolicWorkspace.init(std.testing.allocator);
    defer workspace.deinit();
    const plan_view = try workspace.plan(basis, 0.1);
    try std.testing.expect(!(plan_view.pivot_columns[0] == 0 and plan_view.pivot_rows[0] == 0));
    try std.testing.expectEqual(@as(usize, 1), plan_view.markowitz_pivots);
    try std.testing.expectEqual(@as(usize, 2), plan_view.remaining_dimension);
}

test "symbolic workspace retains storage and produces deterministic ties" {
    const RowId = foundation.RowId;
    const basis = sparse_basis.SparseBasisView{
        .dimension = 2,
        .starts = &[_]Offset{ 0, 2, 4 },
        .rows = &[_]RowId{
            RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(1),
            RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(1),
        },
        .values = &[_]f64{ 2, 1, 1, 2 },
    };
    var workspace = SymbolicWorkspace.init(std.testing.allocator);
    defer workspace.deinit();
    const first = try workspace.plan(basis, 0.1);
    const first_column = first.pivot_columns[0];
    const retained = workspace.retainedBytes();
    const row_pointer = workspace.row_entries.ptr;
    const second = try workspace.plan(basis, 0.1);
    try std.testing.expectEqual(first_column, second.pivot_columns[0]);
    try std.testing.expectEqual(row_pointer, workspace.row_entries.ptr);
    try std.testing.expectEqual(retained, workspace.retainedBytes());
}

test "symbolic planner reports an empty structural row as singular" {
    const RowId = foundation.RowId;
    const basis = sparse_basis.SparseBasisView{
        .dimension = 2,
        .starts = &[_]Offset{ 0, 1, 2 },
        .rows = &[_]RowId{ RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(0) },
        .values = &[_]f64{ 1, 2 },
    };
    var workspace = SymbolicWorkspace.init(std.testing.allocator);
    defer workspace.deinit();
    try std.testing.expectError(error.Singular, workspace.plan(basis, 0.1));
}
