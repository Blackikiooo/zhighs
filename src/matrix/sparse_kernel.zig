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
pub const KernelShape = struct { dimension: usize, nonzeros: usize, maximum_row_count: u32, maximum_column_count: u32 };

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
    cache_column_maxima: bool = false,
    cache_local_candidates: bool = false,
    /// Maximum non-singleton columns inspected after an eligible pivot is
    /// known. Low-count buckets remain first, preserving the Markowitz bias.
    markowitz_candidate_limit: usize = 32,
    /// Candidate budget for the bounded row/column Markowitz search. The
    /// budget is adapted only at fixed windows so one irregular pivot cannot
    /// make ordering policy oscillate.
    markowitz_search_limit: usize = 8,
    // Eight is the validated HiGHS baseline. Adaptation may spend more work
    // when that demonstrably improves merit, but never undercut the baseline.
    markowitz_search_minimum: usize = 8,
    markowitz_search_maximum: usize = 64,
    /// Below this kernel size, bounded-search overhead dominates the possible
    /// reduction in later fill; retain the proven eight-candidate baseline.
    adaptive_markowitz_minimum_dimension: usize = 512,
    last_markowitz_searches: usize = 0,
    last_markowitz_budget_exhausted: bool = false,
    last_markowitz_extension_improved: bool = false,
    markowitz_window_pivots: usize = 0,
    markowitz_window_exhaustions: usize = 0,
    markowitz_window_improvements: usize = 0,

    row_head: []u32 = &.{},
    column_head: []u32 = &.{},
    row_count: []u32 = &.{},
    column_count: []u32 = &.{},
    column_maximum: []f64 = &.{},
    column_maximum_dirty: []bool = &.{},
    local_candidate_dirty: []bool = &.{},
    local_candidate_valid: []bool = &.{},
    local_candidate_row: []u32 = &.{},
    local_candidate_value: []f64 = &.{},
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
    scratch_lookup: []u32 = &.{},
    scratch_lookup_generation: []u32 = &.{},
    lookup_generation: u32 = 0,

    pub fn init(allocator: std.mem.Allocator) MutableSparseKernel {
        return .{ .allocator = allocator };
    }

    pub fn deinit(self: *MutableSparseKernel) void {
        self.allocator.free(self.row_head);
        self.allocator.free(self.column_head);
        self.allocator.free(self.row_count);
        self.allocator.free(self.column_count);
        self.allocator.free(self.column_maximum);
        self.allocator.free(self.column_maximum_dirty);
        self.allocator.free(self.local_candidate_dirty);
        self.allocator.free(self.local_candidate_valid);
        self.allocator.free(self.local_candidate_row);
        self.allocator.free(self.local_candidate_value);
        self.allocator.free(self.row_active);
        self.allocator.free(self.column_active);
        self.allocator.free(self.row_bucket_first);
        self.allocator.free(self.column_bucket_first);
        self.allocator.free(self.row_bucket_next);
        self.allocator.free(self.row_bucket_previous);
        self.allocator.free(self.column_bucket_next);
        self.allocator.free(self.column_bucket_previous);
        self.allocator.free(self.entry_row);
        self.allocator.free(self.entry_column);
        self.allocator.free(self.entry_value);
        self.allocator.free(self.row_next);
        self.allocator.free(self.row_previous);
        self.allocator.free(self.column_next);
        self.allocator.free(self.column_previous);
        self.allocator.free(self.free_next);
        self.allocator.free(self.scratch_rows);
        self.allocator.free(self.scratch_columns);
        self.allocator.free(self.scratch_l);
        self.allocator.free(self.scratch_u);
        self.allocator.free(self.scratch_lookup);
        self.allocator.free(self.scratch_lookup_generation);
        self.* = .{ .allocator = self.allocator };
    }

    /// Bytes requested from the allocator by the retained SoA capacities.
    /// This intentionally excludes allocator metadata and resident pages; it
    /// is the stable, cross-allocator number used beside process peak RSS.
    pub fn requestedBytes(self: *const MutableSparseKernel) usize {
        var total: usize = 0;
        total += std.mem.sliceAsBytes(self.row_head).len;
        total += std.mem.sliceAsBytes(self.column_head).len;
        total += std.mem.sliceAsBytes(self.row_count).len;
        total += std.mem.sliceAsBytes(self.column_count).len;
        total += std.mem.sliceAsBytes(self.column_maximum).len;
        total += std.mem.sliceAsBytes(self.column_maximum_dirty).len;
        total += std.mem.sliceAsBytes(self.local_candidate_dirty).len;
        total += std.mem.sliceAsBytes(self.local_candidate_valid).len;
        total += std.mem.sliceAsBytes(self.local_candidate_row).len;
        total += std.mem.sliceAsBytes(self.local_candidate_value).len;
        total += std.mem.sliceAsBytes(self.row_active).len;
        total += std.mem.sliceAsBytes(self.column_active).len;
        total += std.mem.sliceAsBytes(self.row_bucket_first).len;
        total += std.mem.sliceAsBytes(self.column_bucket_first).len;
        total += std.mem.sliceAsBytes(self.row_bucket_next).len;
        total += std.mem.sliceAsBytes(self.row_bucket_previous).len;
        total += std.mem.sliceAsBytes(self.column_bucket_next).len;
        total += std.mem.sliceAsBytes(self.column_bucket_previous).len;
        total += std.mem.sliceAsBytes(self.entry_row).len;
        total += std.mem.sliceAsBytes(self.entry_column).len;
        total += std.mem.sliceAsBytes(self.entry_value).len;
        total += std.mem.sliceAsBytes(self.row_next).len;
        total += std.mem.sliceAsBytes(self.row_previous).len;
        total += std.mem.sliceAsBytes(self.column_next).len;
        total += std.mem.sliceAsBytes(self.column_previous).len;
        total += std.mem.sliceAsBytes(self.free_next).len;
        total += std.mem.sliceAsBytes(self.scratch_rows).len;
        total += std.mem.sliceAsBytes(self.scratch_columns).len;
        total += std.mem.sliceAsBytes(self.scratch_l).len;
        total += std.mem.sliceAsBytes(self.scratch_u).len;
        total += std.mem.sliceAsBytes(self.scratch_lookup).len;
        total += std.mem.sliceAsBytes(self.scratch_lookup_generation).len;
        return total;
    }

    pub fn load(self: *MutableSparseKernel, basis: sparse_basis.SparseBasisView) KernelError!void {
        return self.loadImpl(basis, true, null, null, false, &.{}, &.{});
    }

    /// Trusted hot-path load for basis views emitted by `SparseBasisBuffers`.
    /// Dimensions, sorted rows and finite nonzero values are not rescanned.
    pub fn loadAssumeValid(self: *MutableSparseKernel, basis: sparse_basis.SparseBasisView) KernelError!void {
        return self.loadImpl(basis, false, null, null, false, &.{}, &.{});
    }

    /// Load only the active reduced kernel after symbolic singleton peeling.
    pub fn loadReducedAssumeValid(self: *MutableSparseKernel, basis: sparse_basis.SparseBasisView, active_rows: []const bool, active_columns: []const bool) KernelError!void {
        if (active_rows.len != basis.dimension or active_columns.len != basis.dimension) return error.InvalidBasis;
        return self.loadImpl(basis, false, active_rows, active_columns, false, &.{}, &.{});
    }

    /// Load a singleton-reduced CSC while reusing the exact active degrees
    /// already maintained by the symbolic planner. CSC materialization remains
    /// column-local, but avoids two count writes per retained entry.
    pub fn loadReducedSymbolicAssumeValid(
        self: *MutableSparseKernel,
        basis: sparse_basis.SparseBasisView,
        active_rows: []const bool,
        active_columns: []const bool,
        row_counts: []const u32,
        column_counts: []const u32,
    ) KernelError!void {
        const n = basis.dimension;
        if (active_rows.len != n or active_columns.len != n or row_counts.len != n or column_counts.len != n)
            return error.InvalidBasis;
        return self.loadImpl(basis, false, active_rows, active_columns, true, row_counts, column_counts);
    }

    fn loadImpl(
        self: *MutableSparseKernel,
        basis: sparse_basis.SparseBasisView,
        comptime validate: bool,
        reduced_rows: ?[]const bool,
        reduced_columns: ?[]const bool,
        comptime reuse_symbolic_counts: bool,
        symbolic_row_counts: []const u32,
        symbolic_column_counts: []const u32,
    ) KernelError!void {
        const n = basis.dimension;
        if (n > std.math.maxInt(u32) or basis.nnz() > std.math.maxInt(u32)) return error.InvalidBasis;
        if (validate and (basis.starts.len != n + 1 or basis.rows.len != basis.values.len)) return error.InvalidBasis;
        try self.ensureDimension(n);
        try self.ensureEntries(@max(basis.nnz(), n));
        if (self.dimension != n) self.markowitz_search_limit = 8;
        self.dimension = n;
        self.entry_high_water = 0;
        self.free_head = none;
        self.buckets_ready = false;
        self.cache_column_maxima = n != 0 and basis.nnz() / n > 4;
        self.cache_local_candidates = self.cache_column_maxima;
        self.last_markowitz_searches = 0;
        self.last_markowitz_budget_exhausted = false;
        self.last_markowitz_extension_improved = false;
        self.markowitz_window_pivots = 0;
        self.markowitz_window_exhaustions = 0;
        self.markowitz_window_improvements = 0;
        @memset(self.row_head[0..n], none);
        @memset(self.column_head[0..n], none);
        if (reuse_symbolic_counts) {
            @memcpy(self.row_count[0..n], symbolic_row_counts);
            @memcpy(self.column_count[0..n], symbolic_column_counts);
        } else {
            @memset(self.row_count[0..n], 0);
            @memset(self.column_count[0..n], 0);
        }
        @memset(self.column_maximum[0..n], 0.0);
        @memset(self.column_maximum_dirty[0..n], false);
        @memset(self.local_candidate_dirty[0..n], true);
        @memset(self.local_candidate_valid[0..n], false);
        @memset(self.scratch_lookup_generation[0..n], 0);
        self.lookup_generation = 0;
        if (reduced_rows) |mask| @memcpy(self.row_active[0..n], mask) else @memset(self.row_active[0..n], true);
        if (reduced_columns) |mask| @memcpy(self.column_active[0..n], mask) else @memset(self.column_active[0..n], true);
        @memset(self.row_bucket_first[0 .. n + 1], none);
        @memset(self.column_bucket_first[0 .. n + 1], none);

        for (0..n) |column| {
            if (!self.column_active[column]) continue;
            const begin: usize = @intCast(basis.starts[column]);
            const end: usize = @intCast(basis.starts[column + 1]);
            if (validate and (begin > end or end > basis.nnz())) return error.InvalidBasis;
            var previous_row: ?usize = null;
            for (begin..end) |position| {
                const row = basis.rows[position].toUsize();
                const value = basis.values[position];
                if (!self.row_active[row]) continue;
                if (validate and (row >= n or !std.math.isFinite(value) or value == 0.0)) return error.InvalidBasis;
                if (validate) if (previous_row) |previous| if (row <= previous) return error.InvalidBasis;
                previous_row = row;
                // Initial CSC positions are already unique dense pool IDs.
                // Materialize both intrusive views directly: routing the cold
                // load through the dynamic fill allocator would add free-list,
                // capacity and bucket branches to every original nonzero.
                const entry: u32 = @intCast(self.entry_high_water);
                self.entry_high_water += 1;
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
                if (!reuse_symbolic_counts) {
                    self.row_count[row] += 1;
                    self.column_count[column] += 1;
                }
                self.column_maximum[column] = @max(self.column_maximum[column], @abs(value));
            }
        }
        for (self.row_count[0..n], self.column_count[0..n], self.row_active[0..n], self.column_active[0..n]) |row_entries, column_entries, row_is_active, column_is_active|
            if ((row_is_active and row_entries == 0) or (column_is_active and column_entries == 0)) return error.Singular;
        for (0..n) |index| {
            if (self.row_active[index]) self.rowBucketInsert(@intCast(index), self.row_count[index]);
            if (self.column_active[index]) self.columnBucketInsert(@intCast(index), self.column_count[index]);
        }
        self.buckets_ready = true;
    }

    /// Threshold Markowitz over the current numerical kernel. Column maxima
    /// are recomputed after every fill update, so pivot eligibility never uses
    /// stale pre-elimination values.
    pub fn choosePivot(self: *MutableSparseKernel, threshold: f64) ?PivotChoice {
        if (!std.math.isFinite(threshold) or threshold <= 0.0 or threshold > 1.0) return null;
        // Singleton pivots create no fill and their sole entry is necessarily
        // the column maximum. Bypass both floating-point threshold scans and
        // the general Markowitz loop. Small kernels retain lowest-ID tie
        // breaking; larger kernels take the deterministic bucket head and
        // avoid a repeated scan that becomes quadratic across many singletons.
        if (self.singletonMember(self.row_bucket_first[1], self.row_bucket_next)) |row| {
            const entry = self.row_head[row];
            if (entry != none) return .{
                .row = row,
                .column = self.entry_column[entry],
                .value = self.entry_value[entry],
                .merit = 0,
            };
        }
        if (self.singletonMember(self.column_bucket_first[1], self.column_bucket_next)) |column| {
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
        var columns_examined: usize = 0;
        while (count <= self.dimension) : (count += 1) {
            if (best != null and count - 1 > best.?.merit) break;
            var bucket_column = self.column_bucket_first[count];
            while (bucket_column != none) : (bucket_column = self.column_bucket_next[bucket_column]) {
                const column: usize = bucket_column;
                columns_examined += 1;
                if (!self.cache_local_candidates) {
                    const maximum = self.columnMaximum(column);
                    if (maximum == 0.0 or !std.math.isFinite(maximum)) continue;
                    const minimum = threshold * maximum;
                    var entry = self.column_head[column];
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
                } else if (self.localColumnCandidate(column, threshold)) |candidate| {
                    if (best == null or candidate.merit < best.?.merit or
                        (candidate.merit == best.?.merit and (column < best.?.column or (column == best.?.column and candidate.row < best.?.row))))
                        best = candidate;
                    const lower_bound = @as(u64, @intCast(count - 1)) * @as(u64, @intCast(minimum_row_count - 1));
                    if (best.?.merit == lower_bound) return best;
                }
                if (best != null and columns_examined >= self.markowitz_candidate_limit) return best;
            }
        }
        return best;
    }

    /// Return a unit/column-singleton or row-singleton pivot without entering
    /// either general Markowitz backend.
    pub fn chooseSingleton(self: *const MutableSparseKernel) ?PivotChoice {
        if (self.column_bucket_first[1] != none) {
            const column = self.column_bucket_first[1];
            const entry = self.column_head[column];
            if (entry != none) return .{ .row = self.entry_row[entry], .column = column, .value = self.entry_value[entry], .merit = 0 };
        }
        if (self.row_bucket_first[1] != none) {
            const row = self.row_bucket_first[1];
            const entry = self.row_head[row];
            if (entry != none) return .{ .row = row, .column = self.entry_column[entry], .value = self.entry_value[entry], .merit = 0 };
        }
        return null;
    }

    /// HiGHS-style bounded kernel search: visit low-count column and row
    /// buckets alternately, accept an opposing side with a smaller degree,
    /// and stop after eight bucket members once a stable candidate exists.
    pub fn choosePivotHighs(self: *MutableSparseKernel, threshold: f64) ?PivotChoice {
        self.last_markowitz_searches = 0;
        self.last_markowitz_budget_exhausted = false;
        self.last_markowitz_extension_improved = false;
        if (!std.math.isFinite(threshold) or threshold <= 0.0 or threshold > 1.0) return null;
        if (self.column_bucket_first[1] != none) {
            const column = self.column_bucket_first[1];
            const entry = self.column_head[column];
            if (entry != none) return .{ .row = self.entry_row[entry], .column = column, .value = self.entry_value[entry], .merit = 0 };
        }
        if (self.row_bucket_first[1] != none) {
            const row = self.row_bucket_first[1];
            const entry = self.row_head[row];
            if (entry != none) return .{ .row = row, .column = self.entry_column[entry], .value = self.entry_value[entry], .merit = 0 };
        }
        var best: ?PivotChoice = null;
        var searched: usize = 0;
        var budget_exhausted = false;
        var half_budget_merit: u64 = std.math.maxInt(u64);
        var half_budget_recorded = false;
        defer {
            self.last_markowitz_searches = searched;
            self.last_markowitz_budget_exhausted = budget_exhausted;
            self.last_markowitz_extension_improved = half_budget_recorded and best != null and
                best.?.merit < half_budget_merit and
                best.?.merit <= half_budget_merit - half_budget_merit / 4;
        }
        const search_limit = self.markowitz_search_limit;
        const half_search_limit = @max(search_limit / 2, 1);
        var count: usize = 2;
        while (count <= self.dimension) : (count += 1) {
            var column = self.column_bucket_first[count];
            while (column != none) : (column = self.column_bucket_next[column]) {
                searched += 1;
                const maximum = self.columnMaximum(column);
                const minimum = threshold * maximum;
                var entry = self.column_head[column];
                while (entry != none) : (entry = self.column_next[entry]) {
                    const value = self.entry_value[entry];
                    if (@abs(value) < minimum) continue;
                    const row = self.entry_row[entry];
                    const merit = @as(u64, self.row_count[row] - 1) * @as(u64, self.column_count[column] - 1);
                    if (best == null or merit < best.?.merit)
                        best = .{ .row = row, .column = column, .value = value, .merit = merit };
                    if (self.row_count[row] < self.column_count[column]) return best;
                }
                if (!half_budget_recorded and searched >= half_search_limit) {
                    half_budget_recorded = true;
                    if (best) |candidate| half_budget_merit = candidate.merit;
                }
                if (searched >= search_limit and best != null) {
                    budget_exhausted = true;
                    return best;
                }
            }
            var row = self.row_bucket_first[count];
            while (row != none) : (row = self.row_bucket_next[row]) {
                searched += 1;
                var entry = self.row_head[row];
                while (entry != none) : (entry = self.row_next[entry]) {
                    const candidate_column = self.entry_column[entry];
                    const merit = @as(u64, self.row_count[row] - 1) * @as(u64, self.column_count[candidate_column] - 1);
                    if (best != null and merit >= best.?.merit) continue;
                    const maximum = self.columnMaximum(candidate_column);
                    const value = self.entry_value[entry];
                    if (@abs(value) < threshold * maximum) continue;
                    best = .{ .row = row, .column = candidate_column, .value = value, .merit = merit };
                    if (self.column_count[candidate_column] <= self.row_count[row]) return best;
                }
                if (!half_budget_recorded and searched >= half_search_limit) {
                    half_budget_recorded = true;
                    if (best) |candidate| half_budget_merit = candidate.merit;
                }
                if (searched >= search_limit and best != null) {
                    budget_exhausted = true;
                    return best;
                }
            }
        }
        return best;
    }

    /// Set the number of low-degree row/column frontier members inspected by
    /// the bounded search. Shape-aware callers may reduce this for compact,
    /// already-peeled kernels where extra candidates cost more than the fill
    /// they avoid; numerical threshold checks remain unchanged.
    pub fn setMarkowitzSearchBudget(self: *MutableSparseKernel, limit: usize) void {
        self.markowitz_search_limit = @max(limit, 1);
    }

    /// Feed the measured search result back into the bounded Markowitz budget.
    /// A wider search is justified only when candidates found in the second
    /// half of the current budget repeatedly improve merit. The fixed
    /// eight-pivot window provides hysteresis and deterministic results.
    pub fn observeMarkowitzPivot(self: *MutableSparseKernel) void {
        if (self.dimension < self.adaptive_markowitz_minimum_dimension) return;
        if (self.last_markowitz_searches == 0) return;
        self.markowitz_window_pivots += 1;
        self.markowitz_window_exhaustions += @intFromBool(self.last_markowitz_budget_exhausted);
        self.markowitz_window_improvements += @intFromBool(self.last_markowitz_extension_improved);
        if (self.markowitz_window_pivots < 8) return;

        if (self.markowitz_window_exhaustions == self.markowitz_window_pivots and self.markowitz_window_improvements >= 4) {
            self.markowitz_search_limit = @min(self.markowitz_search_maximum, self.markowitz_search_limit * 2);
        } else if (self.markowitz_window_improvements == 0) {
            self.markowitz_search_limit = @max(self.markowitz_search_minimum, self.markowitz_search_limit / 2);
        }
        self.markowitz_window_pivots = 0;
        self.markowitz_window_exhaustions = 0;
        self.markowitz_window_improvements = 0;
    }

    fn localColumnCandidate(self: *MutableSparseKernel, column: usize, threshold: f64) ?PivotChoice {
        if (self.local_candidate_dirty[column]) {
            self.local_candidate_dirty[column] = false;
            self.local_candidate_valid[column] = false;
            const maximum = self.columnMaximum(column);
            if (maximum == 0.0 or !std.math.isFinite(maximum)) return null;
            const minimum = threshold * maximum;
            var best_merit: u64 = std.math.maxInt(u64);
            var entry = self.column_head[column];
            while (entry != none) : (entry = self.column_next[entry]) {
                const value = self.entry_value[entry];
                if (@abs(value) < minimum) continue;
                const row = self.entry_row[entry];
                const merit = @as(u64, self.row_count[row] - 1) * @as(u64, self.column_count[column] - 1);
                if (!self.local_candidate_valid[column] or merit < best_merit or
                    (merit == best_merit and row < self.local_candidate_row[column]))
                {
                    self.local_candidate_valid[column] = true;
                    self.local_candidate_row[column] = row;
                    self.local_candidate_value[column] = value;
                    best_merit = merit;
                }
            }
        }
        if (!self.local_candidate_valid[column]) return null;
        const row = self.local_candidate_row[column];
        const merit = @as(u64, self.row_count[row] - 1) * @as(u64, self.column_count[column] - 1);
        return .{ .row = row, .column = @intCast(column), .value = self.local_candidate_value[column], .merit = merit };
    }

    /// Resolve a recorded pivot against the current numerical kernel. Used by
    /// fixed-trace benchmarks to separate ordering cost from update cost.
    pub fn chooseRecordedPivot(self: *const MutableSparseKernel, row: u32, column: u32) ?PivotChoice {
        if (row >= self.dimension or column >= self.dimension or !self.row_active[row] or !self.column_active[column]) return null;
        const entry = self.find(row, column) orelse return null;
        const value = self.entry_value[entry];
        if (value == 0.0 or !std.math.isFinite(value)) return null;
        return .{ .row = row, .column = column, .value = value, .merit = 0 };
    }

    /// Resolve a recorded pivot and verify the same numerical threshold used
    /// by a fresh Markowitz search. A trace from a neighbouring simplex basis
    /// is reusable only while its pivot remains large enough relative to the
    /// current active column; existence alone is not a stability guarantee.
    pub fn chooseRecordedPivotThreshold(self: *MutableSparseKernel, row: u32, column: u32, threshold: f64) ?PivotChoice {
        if (!std.math.isFinite(threshold) or threshold <= 0.0 or threshold > 1.0) return null;
        const choice = self.chooseRecordedPivot(row, column) orelse return null;
        const maximum = self.columnMaximum(column);
        if (maximum == 0.0 or !std.math.isFinite(maximum) or @abs(choice.value) < threshold * maximum) return null;
        return choice;
    }

    fn singletonMember(self: *const MutableSparseKernel, first: u32, next: []const u32) ?u32 {
        if (first == none) return null;
        if (self.dimension >= 128) return first;
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

        // A pivot is an atomic mutation from the perspective of Markowitz
        // selection. Detach every affected dimension once, update counts
        // without exposing intermediate states, then publish final buckets.
        // This replaces potentially many remove/insert bucket pairs per row
        // and column with one pair for the complete Schur update.
        self.rowBucketRemove(choice.row, self.row_count[choice.row]);
        for (self.scratch_rows[0..l_count]) |row|
            self.rowBucketRemove(row, self.row_count[row]);
        self.columnBucketRemove(choice.column, self.column_count[choice.column]);
        for (self.scratch_columns[0..u_count]) |column|
            self.columnBucketRemove(column, self.column_count[column]);
        self.buckets_ready = false;

        var inserted_fill: usize = 0;
        for (self.scratch_rows[0..l_count], self.scratch_l[0..l_count]) |row, multiplier| {
            // Two or more U probes amortize one traversal of the target row.
            // Generation marks avoid clearing the dense lookup between rows.
            const use_lookup = u_count >= 2;
            var generation: u32 = 0;
            if (use_lookup) {
                generation = self.nextLookupGeneration();
                entry = self.row_head[row];
                while (entry != none) : (entry = self.row_next[entry]) {
                    self.scratch_lookup[self.entry_column[entry]] = entry;
                    self.scratch_lookup_generation[self.entry_column[entry]] = generation;
                }
            }
            for (self.scratch_columns[0..u_count], self.scratch_u[0..u_count]) |column, upper| {
                const delta = multiplier * upper;
                const existing_entry = if (use_lookup)
                    (if (self.scratch_lookup_generation[column] == generation) self.scratch_lookup[column] else null)
                else
                    self.find(row, column);
                if (existing_entry) |existing| {
                    const updated = self.entry_value[existing] - delta;
                    if (!std.math.isFinite(updated)) return error.Singular;
                    if (@abs(updated) <= zero_tolerance) self.removeUnbucketed(existing) else {
                        self.entry_value[existing] = updated;
                        if (self.cache_column_maxima) self.column_maximum_dirty[column] = true;
                        if (self.cache_local_candidates) self.local_candidate_dirty[column] = true;
                    }
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
        self.row_active[choice.row] = false;
        self.column_active[choice.column] = false;
        entry = self.row_head[choice.row];
        while (entry != none) {
            const next = self.row_next[entry];
            self.removePivotRowEntry(entry);
            removed_entries += 1;
            entry = next;
        }
        self.row_head[choice.row] = none;
        self.row_count[choice.row] = 0;
        entry = self.column_head[choice.column];
        while (entry != none) {
            const next = self.column_next[entry];
            self.removePivotColumnEntry(entry);
            removed_entries += 1;
            entry = next;
        }
        self.column_head[choice.column] = none;
        self.column_count[choice.column] = 0;
        for (self.scratch_rows[0..l_count]) |row|
            self.rowBucketInsert(row, self.row_count[row]);
        for (self.scratch_columns[0..u_count]) |column|
            self.columnBucketInsert(column, self.column_count[column]);
        self.buckets_ready = true;
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

    pub fn shape(self: *const MutableSparseKernel) KernelShape {
        var nonzeros: usize = 0;
        var maximum_row_count: u32 = 0;
        var maximum_column_count: u32 = 0;
        var active_dimension: usize = 0;
        for (0..self.dimension) |index| {
            if (self.row_active[index]) {
                active_dimension += 1;
                nonzeros += self.row_count[index];
                maximum_row_count = @max(maximum_row_count, self.row_count[index]);
            }
            if (self.column_active[index]) maximum_column_count = @max(maximum_column_count, self.column_count[index]);
        }
        return .{ .dimension = active_dimension, .nonzeros = nonzeros, .maximum_row_count = maximum_row_count, .maximum_column_count = maximum_column_count };
    }

    fn find(self: *const MutableSparseKernel, row: u32, column: u32) ?u32 {
        if (self.row_count[row] <= self.column_count[column]) {
            var entry = self.row_head[row];
            while (entry != none) : (entry = self.row_next[entry])
                if (self.entry_column[entry] == column) return entry;
        } else {
            var entry = self.column_head[column];
            while (entry != none) : (entry = self.column_next[entry])
                if (self.entry_row[entry] == row) return entry;
        }
        return null;
    }

    /// Begin a new logical sparse-lookup row without clearing its dense index
    /// array. A full O(n) clear occurs only after the u32 generation wraps.
    fn nextLookupGeneration(self: *MutableSparseKernel) u32 {
        self.lookup_generation +%= 1;
        if (self.lookup_generation == 0) {
            @memset(self.scratch_lookup_generation[0..self.dimension], 0);
            self.lookup_generation = 1;
        }
        return self.lookup_generation;
    }

    fn columnMaximum(self: *MutableSparseKernel, column: usize) f64 {
        if (self.cache_column_maxima and !self.column_maximum_dirty[column]) return self.column_maximum[column];
        var maximum: f64 = 0.0;
        var entry = self.column_head[column];
        while (entry != none) : (entry = self.column_next[entry])
            maximum = @max(maximum, @abs(self.entry_value[entry]));
        self.column_maximum[column] = maximum;
        if (self.cache_column_maxima) self.column_maximum_dirty[column] = false;
        return maximum;
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
        if (self.cache_column_maxima) self.column_maximum_dirty[column] = true;
        if (self.cache_local_candidates) self.local_candidate_dirty[column] = true;
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

    /// Remove a Schur entry while all affected dimensions are detached from
    /// buckets. Keeping this hot path branch-free avoids repeated active-state
    /// and bucket checks inside the numerical update.
    fn removeUnbucketed(self: *MutableSparseKernel, entry: u32) void {
        const row = self.entry_row[entry];
        const column = self.entry_column[entry];
        if (self.cache_column_maxima) self.column_maximum_dirty[column] = true;
        if (self.cache_local_candidates) self.local_candidate_dirty[column] = true;
        const row_prev = self.row_previous[entry];
        const row_after = self.row_next[entry];
        if (row_prev == none) self.row_head[row] = row_after else self.row_next[row_prev] = row_after;
        if (row_after != none) self.row_previous[row_after] = row_prev;
        const col_prev = self.column_previous[entry];
        const col_after = self.column_next[entry];
        if (col_prev == none) self.column_head[column] = col_after else self.column_next[col_prev] = col_after;
        if (col_after != none) self.column_previous[col_after] = col_prev;
        self.row_count[row] -= 1;
        self.column_count[column] -= 1;
        self.releaseEntry(entry);
    }

    /// Retire one entry from the completed pivot row. The entire row list is
    /// discarded, so only its column links and live column count are touched.
    fn removePivotRowEntry(self: *MutableSparseKernel, entry: u32) void {
        const column = self.entry_column[entry];
        if (self.column_active[column]) {
            if (self.cache_column_maxima) self.column_maximum_dirty[column] = true;
            if (self.cache_local_candidates) self.local_candidate_dirty[column] = true;
            self.column_count[column] -= 1;
        }
        const previous = self.column_previous[entry];
        const next = self.column_next[entry];
        if (previous == none) self.column_head[column] = next else self.column_next[previous] = next;
        if (next != none) self.column_previous[next] = previous;
        self.releaseEntry(entry);
    }

    /// Retire one entry from the completed pivot column. Its column links are
    /// already dead, so only the live row list and row count are updated.
    fn removePivotColumnEntry(self: *MutableSparseKernel, entry: u32) void {
        const row = self.entry_row[entry];
        const previous = self.row_previous[entry];
        const next = self.row_next[entry];
        if (previous == none) self.row_head[row] = next else self.row_next[previous] = next;
        if (next != none) self.row_previous[next] = previous;
        self.row_count[row] -= 1;
        self.releaseEntry(entry);
    }

    inline fn releaseEntry(self: *MutableSparseKernel, entry: u32) void {
        self.free_next[entry] = self.free_head;
        self.free_head = entry;
    }

    fn ensureDimension(self: *MutableSparseKernel, required: usize) KernelError!void {
        if (required <= self.dimension_capacity) return;
        const capacity = grow(self.dimension_capacity, required) catch return error.CapacityOverflow;
        try self.resizeRetained(&self.row_head, capacity);
        try self.resizeRetained(&self.column_head, capacity);
        try self.resizeRetained(&self.row_count, capacity);
        try self.resizeRetained(&self.column_count, capacity);
        try self.resizeRetained(&self.column_maximum, capacity);
        try self.resizeRetained(&self.column_maximum_dirty, capacity);
        try self.resizeRetained(&self.local_candidate_dirty, capacity);
        try self.resizeRetained(&self.local_candidate_valid, capacity);
        try self.resizeRetained(&self.local_candidate_row, capacity);
        try self.resizeRetained(&self.local_candidate_value, capacity);
        try self.resizeRetained(&self.row_active, capacity);
        try self.resizeRetained(&self.column_active, capacity);
        try self.resizeRetained(&self.row_bucket_next, capacity);
        try self.resizeRetained(&self.row_bucket_previous, capacity);
        try self.resizeRetained(&self.column_bucket_next, capacity);
        try self.resizeRetained(&self.column_bucket_previous, capacity);
        try self.resizeRetained(&self.scratch_rows, capacity);
        try self.resizeRetained(&self.scratch_columns, capacity);
        try self.resizeRetained(&self.scratch_l, capacity);
        try self.resizeRetained(&self.scratch_u, capacity);
        try self.resizeRetained(&self.scratch_lookup, capacity);
        try self.resizeRetained(&self.scratch_lookup_generation, capacity);
        self.row_bucket_first = self.allocator.realloc(self.row_bucket_first, capacity + 1) catch return error.OutOfMemory;
        self.column_bucket_first = self.allocator.realloc(self.column_bucket_first, capacity + 1) catch return error.OutOfMemory;
        self.dimension_capacity = capacity;
    }

    fn ensureEntries(self: *MutableSparseKernel, required: usize) KernelError!void {
        if (required <= self.entry_capacity) return;
        const capacity = grow(self.entry_capacity, required) catch return error.CapacityOverflow;
        try self.resizeRetained(&self.entry_row, capacity);
        try self.resizeRetained(&self.entry_column, capacity);
        try self.resizeRetained(&self.entry_value, capacity);
        try self.resizeRetained(&self.row_next, capacity);
        try self.resizeRetained(&self.row_previous, capacity);
        try self.resizeRetained(&self.column_next, capacity);
        try self.resizeRetained(&self.column_previous, capacity);
        try self.resizeRetained(&self.free_next, capacity);
        self.entry_capacity = capacity;
    }

    fn resizeRetained(self: *MutableSparseKernel, slice: anytype, capacity: usize) KernelError!void {
        slice.* = self.allocator.realloc(slice.*, capacity) catch return error.OutOfMemory;
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
    kernel.dimension = 512;
    try kernel.load(basis);
    const pivot = kernel.choosePivot(0.1).?;
    const result = try kernel.applyPivot(pivot, 1e-14);
    try std.testing.expectEqual(@as(usize, 1), result.inserted_fill);
    try std.testing.expectEqual(@as(usize, 4), kernel.activeEntries());
    try std.testing.expect(kernel.choosePivot(0.1) != null);
}

test "symbolic degree reuse matches ordinary reduced load" {
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
    const active = [_]bool{ false, true, true };
    const counts = [_]u32{ 0, 1, 2 };
    var ordinary = MutableSparseKernel.init(std.testing.allocator);
    defer ordinary.deinit();
    var fused = MutableSparseKernel.init(std.testing.allocator);
    defer fused.deinit();
    try ordinary.loadReducedAssumeValid(basis, &active, &active);
    try fused.loadReducedSymbolicAssumeValid(basis, &active, &active, &counts, &counts);

    try std.testing.expectEqual(ordinary.shape(), fused.shape());
    try std.testing.expectEqual(ordinary.choosePivot(0.1).?, fused.choosePivot(0.1).?);
}

test "Markowitz search budget adapts only after stable windows" {
    var kernel = MutableSparseKernel.init(std.testing.allocator);
    defer kernel.deinit();
    kernel.markowitz_search_limit = 8;
    for (0..7) |_| {
        kernel.last_markowitz_searches = 8;
        kernel.last_markowitz_budget_exhausted = true;
        kernel.last_markowitz_extension_improved = true;
        kernel.observeMarkowitzPivot();
    }
    try std.testing.expectEqual(@as(usize, 8), kernel.markowitz_search_limit);
    kernel.last_markowitz_searches = 8;
    kernel.last_markowitz_budget_exhausted = true;
    kernel.last_markowitz_extension_improved = true;
    kernel.observeMarkowitzPivot();
    try std.testing.expectEqual(@as(usize, 16), kernel.markowitz_search_limit);

    kernel.markowitz_search_minimum = 4;
    for (0..8) |_| {
        kernel.last_markowitz_searches = 2;
        kernel.last_markowitz_budget_exhausted = false;
        kernel.last_markowitz_extension_improved = false;
        kernel.observeMarkowitzPivot();
    }
    try std.testing.expectEqual(@as(usize, 8), kernel.markowitz_search_limit);
}
