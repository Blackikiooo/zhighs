//! Mutable data-oriented sparse elimination kernel.
//!
//! Entries live in reusable SoA pools and participate in intrusive row and
//! column lists. Removing an entry is O(1); fill lookup scans one active row,
//! which is preferable to a hash table for the small fronts produced after
//! singleton elimination. Freed slots are recycled before the pool grows.

const std = @import("std");
const foundation = @import("foundation");
const sparse_basis = @import("sparse_basis.zig");

const none = std.math.maxInt(u32);

pub const KernelError = error{ InvalidBasis, Singular, OutOfMemory, CapacityOverflow };

pub const PivotChoice = struct { row: u32, column: u32, value: f64, merit: u64 };

/// Borrowed elimination data valid until the next kernel mutation.
pub const PivotView = struct {
    pivot_row: u32,
    pivot_column: u32,
    pivot_value: f64,
    l_rows: []const u32,
    l_values: []const f64,
    u_columns: []const u32,
    u_values: []const f64,
    inserted_fill: usize,
    removed_entries: usize,
};

pub const MutableSparseKernel = struct {
    allocator: std.mem.Allocator,
    dimension_capacity: usize = 0,
    entry_capacity: usize = 0,
    dimension: usize = 0,
    entry_high_water: usize = 0,
    free_head: u32 = none,
    buckets_ready: bool = false,

    row_head: []u32 = &.{},
    column_head: []u32 = &.{},
    row_count: []u32 = &.{},
    column_count: []u32 = &.{},
    row_active: []bool = &.{},
    column_active: []bool = &.{},
    row_bucket_first: []u32 = &.{},
    column_bucket_first: []u32 = &.{},
    row_bucket_next: []u32 = &.{},
    row_bucket_previous: []u32 = &.{},
    column_bucket_next: []u32 = &.{},
    column_bucket_previous: []u32 = &.{},

    entry_row: []u32 = &.{},
    entry_column: []u32 = &.{},
    entry_value: []f64 = &.{},
    row_next: []u32 = &.{},
    row_previous: []u32 = &.{},
    column_next: []u32 = &.{},
    column_previous: []u32 = &.{},
    free_next: []u32 = &.{},

    scratch_rows: []u32 = &.{},
    scratch_columns: []u32 = &.{},
    scratch_l: []f64 = &.{},
    scratch_u: []f64 = &.{},

    pub fn init(allocator: std.mem.Allocator) MutableSparseKernel {
        return .{ .allocator = allocator };
    }

    pub fn deinit(self: *MutableSparseKernel) void {
        inline for (.{
            "row_head",         "column_head",         "row_count",       "column_count",        "row_active",         "column_active",
            "row_bucket_first", "column_bucket_first", "row_bucket_next", "row_bucket_previous", "column_bucket_next", "column_bucket_previous",
            "entry_row",        "entry_column",        "entry_value",     "row_next",            "row_previous",       "column_next",
            "column_previous",  "free_next",           "scratch_rows",    "scratch_columns",     "scratch_l",          "scratch_u",
        }) |field_name| self.allocator.free(@field(self, field_name));
        self.* = .{ .allocator = self.allocator };
    }

    /// Bytes requested from the allocator by the retained SoA capacities.
    /// This intentionally excludes allocator metadata and resident pages; it
    /// is the stable, cross-allocator number used beside process peak RSS.
    pub fn requestedBytes(self: *const MutableSparseKernel) usize {
        var total: usize = 0;
        inline for (.{
            "row_head", "column_head", "row_count", "column_count",
            "row_bucket_first", "column_bucket_first", "row_bucket_next", "row_bucket_previous",
            "column_bucket_next", "column_bucket_previous", "entry_row", "entry_column",
            "row_next", "row_previous", "column_next", "column_previous", "free_next",
            "scratch_rows", "scratch_columns",
        }) |name| total += @sizeOf(std.meta.Elem(@TypeOf(@field(self, name)))) * @field(self, name).len;
        inline for (.{ "row_active", "column_active" }) |name|
            total += @sizeOf(std.meta.Elem(@TypeOf(@field(self, name)))) * @field(self, name).len;
        inline for (.{ "entry_value", "scratch_l", "scratch_u" }) |name|
            total += @sizeOf(std.meta.Elem(@TypeOf(@field(self, name)))) * @field(self, name).len;
        return total;
    }

    pub fn load(self: *MutableSparseKernel, basis: sparse_basis.SparseBasisView) KernelError!void {
        return self.loadImpl(basis, true);
    }

    /// Trusted hot-path load for basis views emitted by `SparseBasisBuffers`.
    /// Dimensions, sorted rows and finite nonzero values are not rescanned.
    pub fn loadAssumeValid(self: *MutableSparseKernel, basis: sparse_basis.SparseBasisView) KernelError!void {
        return self.loadImpl(basis, false);
    }

    fn loadImpl(self: *MutableSparseKernel, basis: sparse_basis.SparseBasisView, comptime validate: bool) KernelError!void {
        const n = basis.dimension;
        if (n > std.math.maxInt(u32) or basis.nnz() > std.math.maxInt(u32)) return error.InvalidBasis;
        if (validate and (basis.starts.len != n + 1 or basis.rows.len != basis.values.len)) return error.InvalidBasis;
        try self.ensureDimension(n);
        try self.ensureEntries(@max(basis.nnz(), n));
        self.dimension = n;
        self.entry_high_water = basis.nnz();
        self.free_head = none;
        self.buckets_ready = false;
        @memset(self.row_head[0..n], none);
        @memset(self.column_head[0..n], none);
        @memset(self.row_count[0..n], 0);
        @memset(self.column_count[0..n], 0);
        @memset(self.row_active[0..n], true);
        @memset(self.column_active[0..n], true);
        @memset(self.row_bucket_first[0 .. n + 1], none);
        @memset(self.column_bucket_first[0 .. n + 1], none);

        for (0..n) |column| {
            const begin: usize = @intCast(basis.starts[column]);
            const end: usize = @intCast(basis.starts[column + 1]);
            if (validate and (begin > end or end > basis.nnz())) return error.InvalidBasis;
            var previous_row: ?usize = null;
            for (begin..end) |position| {
                const row = basis.rows[position].toUsize();
                const value = basis.values[position];
                if (validate and (row >= n or !std.math.isFinite(value) or value == 0.0)) return error.InvalidBasis;
                if (validate) if (previous_row) |previous| if (row <= previous) return error.InvalidBasis;
                previous_row = row;
                // Initial CSC positions are already unique dense pool IDs.
                // Materialize both intrusive views directly: routing the cold
                // load through the dynamic fill allocator would add free-list,
                // capacity and bucket branches to every original nonzero.
                const entry: u32 = @intCast(position);
                self.entry_row[entry] = @intCast(row);
                self.entry_column[entry] = @intCast(column);
                self.entry_value[entry] = value;
                self.row_previous[entry] = none;
                self.row_next[entry] = self.row_head[row];
                if (self.row_head[row] != none) self.row_previous[self.row_head[row]] = entry;
                self.row_head[row] = entry;
                self.column_previous[entry] = none;
                self.column_next[entry] = self.column_head[column];
                if (self.column_head[column] != none) self.column_previous[self.column_head[column]] = entry;
                self.column_head[column] = entry;
                self.row_count[row] += 1;
                self.column_count[column] += 1;
            }
        }
        for (self.row_count[0..n], self.column_count[0..n]) |row_entries, column_entries|
            if (row_entries == 0 or column_entries == 0) return error.Singular;
        for (0..n) |index| {
            self.rowBucketInsert(@intCast(index), self.row_count[index]);
            self.columnBucketInsert(@intCast(index), self.column_count[index]);
        }
        self.buckets_ready = true;
    }

    /// Threshold Markowitz over the current numerical kernel. Column maxima
    /// are recomputed after every fill update, so pivot eligibility never uses
    /// stale pre-elimination values.
    pub fn choosePivot(self: *const MutableSparseKernel, threshold: f64) ?PivotChoice {
        if (!std.math.isFinite(threshold) or threshold <= 0.0 or threshold > 1.0) return null;
        // Singleton pivots create no fill and their sole entry is necessarily
        // the column maximum. Bypass both floating-point threshold scans and
        // the general Markowitz loop; deterministic ID selection keeps runs
        // and factor layouts reproducible despite intrusive-list order.
        if (self.lowestBucketMember(self.row_bucket_first[1], self.row_bucket_next)) |row| {
            const entry = self.row_head[row];
            if (entry != none) return .{
                .row = row,
                .column = self.entry_column[entry],
                .value = self.entry_value[entry],
                .merit = 0,
            };
        }
        if (self.lowestBucketMember(self.column_bucket_first[1], self.column_bucket_next)) |column| {
            const entry = self.column_head[column];
            if (entry != none) return .{
                .row = self.entry_row[entry],
                .column = column,
                .value = self.entry_value[entry],
                .merit = 0,
            };
        }
        var best: ?PivotChoice = null;
        var minimum_row_count: usize = 2;
        while (minimum_row_count <= self.dimension and self.row_bucket_first[minimum_row_count] == none) : (minimum_row_count += 1) {}
        if (minimum_row_count > self.dimension) return null;
        var count: usize = 1;
        while (count <= self.dimension) : (count += 1) {
            if (best != null and count - 1 > best.?.merit) break;
            var bucket_column = self.column_bucket_first[count];
            while (bucket_column != none) : (bucket_column = self.column_bucket_next[bucket_column]) {
                const column: usize = bucket_column;
                var maximum: f64 = 0.0;
                var entry = self.column_head[column];
                while (entry != none) : (entry = self.column_next[entry])
                    maximum = @max(maximum, @abs(self.entry_value[entry]));
                if (maximum == 0.0 or !std.math.isFinite(maximum)) continue;
                const minimum = threshold * maximum;
                entry = self.column_head[column];
                while (entry != none) : (entry = self.column_next[entry]) {
                    const value = self.entry_value[entry];
                    if (@abs(value) < minimum) continue;
                    const row = self.entry_row[entry];
                    const merit = @as(u64, self.row_count[row] - 1) * @as(u64, self.column_count[column] - 1);
                    if (best == null or merit < best.?.merit or
                        (merit == best.?.merit and (column < best.?.column or (column == best.?.column and row < best.?.row))))
                        best = .{ .row = row, .column = @intCast(column), .value = value, .merit = merit };
                    const lower_bound = @as(u64, @intCast(count - 1)) * @as(u64, @intCast(minimum_row_count - 1));
                    if (best.?.merit == lower_bound) return best;
                }
            }
        }
        return best;
    }

    fn lowestBucketMember(_: *const MutableSparseKernel, first: u32, next: []const u32) ?u32 {
        if (first == none) return null;
        var lowest = first;
        var current = next[first];
        while (current != none) : (current = next[current]) lowest = @min(lowest, current);
        return lowest;
    }

    /// Apply one Gaussian pivot, materialize L multipliers/U row, insert fill,
    /// drop numerical zeros, and remove the completed pivot row and column.
    pub fn applyPivot(self: *MutableSparseKernel, choice: PivotChoice, zero_tolerance: f64) KernelError!PivotView {
        if (choice.row >= self.dimension or choice.column >= self.dimension or
            !self.row_active[choice.row] or !self.column_active[choice.column] or
            !std.math.isFinite(choice.value) or choice.value == 0.0 or zero_tolerance < 0.0)
            return error.Singular;
        const pivot_entry = self.find(choice.row, choice.column) orelse return error.Singular;
        const pivot_value = self.entry_value[pivot_entry];

        var l_count: usize = 0;
        var entry = self.column_head[choice.column];
        while (entry != none) : (entry = self.column_next[entry]) {
            const row = self.entry_row[entry];
            if (row != choice.row) {
                self.scratch_rows[l_count] = row;
                self.scratch_l[l_count] = self.entry_value[entry] / pivot_value;
                l_count += 1;
            }
        }
        var u_count: usize = 0;
        entry = self.row_head[choice.row];
        while (entry != none) : (entry = self.row_next[entry]) {
            const column = self.entry_column[entry];
            if (column != choice.column) {
                self.scratch_columns[u_count] = column;
                self.scratch_u[u_count] = self.entry_value[entry];
                u_count += 1;
            }
        }

        var inserted_fill: usize = 0;
        for (self.scratch_rows[0..l_count], self.scratch_l[0..l_count]) |row, multiplier| {
            for (self.scratch_columns[0..u_count], self.scratch_u[0..u_count]) |column, upper| {
                const delta = multiplier * upper;
                if (self.find(row, column)) |existing| {
                    const updated = self.entry_value[existing] - delta;
                    if (!std.math.isFinite(updated)) return error.Singular;
                    if (@abs(updated) <= zero_tolerance) self.remove(existing) else self.entry_value[existing] = updated;
                } else if (@abs(delta) > zero_tolerance) {
                    _ = try self.insert(row, column, -delta);
                    inserted_fill += 1;
                }
            }
        }

        var removed_entries: usize = 0;
        // Retire the pivot dimensions before unlinking their entries. Their
        // own counts are dead after this pivot, so moving them through every
        // intermediate count bucket is wasted work. Neighboring active rows
        // and columns still relocate exactly once per removed entry.
        self.rowBucketRemove(choice.row, self.row_count[choice.row]);
        self.columnBucketRemove(choice.column, self.column_count[choice.column]);
        self.row_active[choice.row] = false;
        self.column_active[choice.column] = false;
        entry = self.row_head[choice.row];
        while (entry != none) {
            const next = self.row_next[entry];
            self.remove(entry);
            removed_entries += 1;
            entry = next;
        }
        entry = self.column_head[choice.column];
        while (entry != none) {
            const next = self.column_next[entry];
            self.remove(entry);
            removed_entries += 1;
            entry = next;
        }
        return .{
            .pivot_row = choice.row,
            .pivot_column = choice.column,
            .pivot_value = pivot_value,
            .l_rows = self.scratch_rows[0..l_count],
            .l_values = self.scratch_l[0..l_count],
            .u_columns = self.scratch_columns[0..u_count],
            .u_values = self.scratch_u[0..u_count],
            .inserted_fill = inserted_fill,
            .removed_entries = removed_entries,
        };
    }

    pub fn activeEntries(self: *const MutableSparseKernel) usize {
        var total: usize = 0;
        for (self.row_count[0..self.dimension]) |count| total += count;
        return total;
    }

    fn find(self: *const MutableSparseKernel, row: u32, column: u32) ?u32 {
        var entry = self.row_head[row];
        while (entry != none) : (entry = self.row_next[entry])
            if (self.entry_column[entry] == column) return entry;
        return null;
    }

    fn insert(self: *MutableSparseKernel, row: u32, column: u32, value: f64) KernelError!u32 {
        if (self.entry_high_water == self.entry_capacity and self.free_head == none)
            try self.ensureEntries(self.entry_capacity + self.entry_capacity / 2 + 8);
        const entry: u32 = if (self.free_head != none) blk: {
            const reused = self.free_head;
            self.free_head = self.free_next[reused];
            break :blk reused;
        } else blk: {
            const fresh: u32 = @intCast(self.entry_high_water);
            self.entry_high_water += 1;
            break :blk fresh;
        };
        self.entry_row[entry] = row;
        self.entry_column[entry] = column;
        self.entry_value[entry] = value;
        self.row_previous[entry] = none;
        self.row_next[entry] = self.row_head[row];
        if (self.row_head[row] != none) self.row_previous[self.row_head[row]] = entry;
        self.row_head[row] = entry;
        self.column_previous[entry] = none;
        self.column_next[entry] = self.column_head[column];
        if (self.column_head[column] != none) self.column_previous[self.column_head[column]] = entry;
        self.column_head[column] = entry;
        if (self.buckets_ready) {
            self.rowBucketRemove(row, self.row_count[row]);
            self.columnBucketRemove(column, self.column_count[column]);
        }
        self.row_count[row] += 1;
        self.column_count[column] += 1;
        if (self.buckets_ready) {
            self.rowBucketInsert(row, self.row_count[row]);
            self.columnBucketInsert(column, self.column_count[column]);
        }
        return entry;
    }

    fn remove(self: *MutableSparseKernel, entry: u32) void {
        const row = self.entry_row[entry];
        const column = self.entry_column[entry];
        const row_prev = self.row_previous[entry];
        const row_after = self.row_next[entry];
        if (row_prev == none) self.row_head[row] = row_after else self.row_next[row_prev] = row_after;
        if (row_after != none) self.row_previous[row_after] = row_prev;
        const col_prev = self.column_previous[entry];
        const col_after = self.column_next[entry];
        if (col_prev == none) self.column_head[column] = col_after else self.column_next[col_prev] = col_after;
        if (col_after != none) self.column_previous[col_after] = col_prev;
        if (self.buckets_ready and self.row_active[row]) {
            self.rowBucketRemove(row, self.row_count[row]);
        }
        if (self.buckets_ready and self.column_active[column]) {
            self.columnBucketRemove(column, self.column_count[column]);
        }
        self.row_count[row] -= 1;
        self.column_count[column] -= 1;
        if (self.buckets_ready and self.row_active[row]) {
            self.rowBucketInsert(row, self.row_count[row]);
        }
        if (self.buckets_ready and self.column_active[column]) {
            self.columnBucketInsert(column, self.column_count[column]);
        }
        self.free_next[entry] = self.free_head;
        self.free_head = entry;
    }

    fn ensureDimension(self: *MutableSparseKernel, required: usize) KernelError!void {
        if (required <= self.dimension_capacity) return;
        const capacity = grow(self.dimension_capacity, required) catch return error.CapacityOverflow;
        inline for (.{ "row_head", "column_head", "row_count", "column_count", "row_active", "column_active", "row_bucket_next", "row_bucket_previous", "column_bucket_next", "column_bucket_previous", "scratch_rows", "scratch_columns", "scratch_l", "scratch_u" }) |field_name|
            @field(self, field_name) = self.allocator.realloc(@field(self, field_name), capacity) catch return error.OutOfMemory;
        self.row_bucket_first = self.allocator.realloc(self.row_bucket_first, capacity + 1) catch return error.OutOfMemory;
        self.column_bucket_first = self.allocator.realloc(self.column_bucket_first, capacity + 1) catch return error.OutOfMemory;
        self.dimension_capacity = capacity;
    }

    fn ensureEntries(self: *MutableSparseKernel, required: usize) KernelError!void {
        if (required <= self.entry_capacity) return;
        const capacity = grow(self.entry_capacity, required) catch return error.CapacityOverflow;
        inline for (.{ "entry_row", "entry_column", "entry_value", "row_next", "row_previous", "column_next", "column_previous", "free_next" }) |field_name|
            @field(self, field_name) = self.allocator.realloc(@field(self, field_name), capacity) catch return error.OutOfMemory;
        self.entry_capacity = capacity;
    }

    fn rowBucketInsert(self: *MutableSparseKernel, row: u32, count: u32) void {
        const first = self.row_bucket_first[count];
        self.row_bucket_previous[row] = none;
        self.row_bucket_next[row] = first;
        if (first != none) self.row_bucket_previous[first] = row;
        self.row_bucket_first[count] = row;
    }

    fn rowBucketRemove(self: *MutableSparseKernel, row: u32, count: u32) void {
        const previous = self.row_bucket_previous[row];
        const next = self.row_bucket_next[row];
        if (previous == none) self.row_bucket_first[count] = next else self.row_bucket_next[previous] = next;
        if (next != none) self.row_bucket_previous[next] = previous;
        self.row_bucket_previous[row] = none;
        self.row_bucket_next[row] = none;
    }

    fn columnBucketInsert(self: *MutableSparseKernel, column: u32, count: u32) void {
        const first = self.column_bucket_first[count];
        self.column_bucket_previous[column] = none;
        self.column_bucket_next[column] = first;
        if (first != none) self.column_bucket_previous[first] = column;
        self.column_bucket_first[count] = column;
    }

    fn columnBucketRemove(self: *MutableSparseKernel, column: u32, count: u32) void {
        const previous = self.column_bucket_previous[column];
        const next = self.column_bucket_next[column];
        if (previous == none) self.column_bucket_first[count] = next else self.column_bucket_next[previous] = next;
        if (next != none) self.column_bucket_previous[next] = previous;
        self.column_bucket_previous[column] = none;
        self.column_bucket_next[column] = none;
    }
};

fn grow(current: usize, required: usize) error{Overflow}!usize {
    var capacity = @max(current, 8);
    while (capacity < required) capacity = std.math.add(usize, capacity, capacity / 2 + 8) catch return error.Overflow;
    if (capacity > std.math.maxInt(u32)) return error.Overflow;
    return capacity;
}

test "mutable kernel inserts fill and updates row column counts" {
    const RowId = foundation.RowId;
    const Offset = foundation.HUInt;
    const basis = sparse_basis.SparseBasisView{
        .dimension = 3,
        .starts = &[_]Offset{ 0, 2, 4, 6 },
        .rows = &[_]RowId{
            RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(1),
            RowId.fromUsizeAssumeValid(0), RowId.fromUsizeAssumeValid(2),
            RowId.fromUsizeAssumeValid(1), RowId.fromUsizeAssumeValid(2),
        },
        .values = &[_]f64{ 2, 1, 1, 2, 3, 4 },
    };
    var kernel = MutableSparseKernel.init(std.testing.allocator);
    defer kernel.deinit();
    try kernel.load(basis);
    const pivot = kernel.choosePivot(0.1).?;
    const result = try kernel.applyPivot(pivot, 1e-14);
    try std.testing.expectEqual(@as(usize, 1), result.inserted_fill);
    try std.testing.expectEqual(@as(usize, 4), kernel.activeEntries());
    try std.testing.expect(kernel.choosePivot(0.1) != null);
}
