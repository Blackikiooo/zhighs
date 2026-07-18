//! Mutable Forrest--Tomlin upper factor and row-correction workspace.
//!
//! Pivot identifiers are simplex basis-column IDs. U is stored simultaneously
//! by column and by row (UR) through one intrusive SoA entry pool, allowing a
//! pivotal row and column to be removed without rebuilding unrelated segments.

const std = @import("std");

const none = std.math.maxInt(u32);

pub const FtError = error{ DimensionMismatch, Singular, NumericalFailure, OutOfMemory, CapacityOverflow };

pub const MutableUpperView = struct {
    pivot_ids: []const u32,
    pivot_values: []f64,
    pivot_lookup: []const u32,
    column_heads: []const u32,
    row_heads: []const u32,
    entry_rows: []const u32,
    entry_columns: []const u32,
    entry_values: []f64,
};

pub const SparseForrestTomlin = struct {
    allocator: std.mem.Allocator,
    dimension: usize = 0,
    logical_count: usize = 0,
    logical_capacity: usize = 0,
    entry_high_water: usize = 0,
    entry_capacity: usize = 0,
    free_head: u32 = none,
    update_count: usize = 0,

    pivot_ids: []u32 = &.{},
    pivot_values: []f64 = &.{},
    pivot_lookup: []u32 = &.{},
    column_head: []u32 = &.{},
    row_head: []u32 = &.{},

    entry_row: []u32 = &.{},
    entry_column: []u32 = &.{},
    entry_value: []f64 = &.{},
    column_next: []u32 = &.{},
    column_previous: []u32 = &.{},
    row_next: []u32 = &.{},
    row_previous: []u32 = &.{},
    free_next: []u32 = &.{},

    correction_starts: []usize = &.{},
    correction_pivots: []u32 = &.{},
    correction_indices: []u32 = &.{},
    correction_values: []f64 = &.{},
    correction_count: usize = 0,
    correction_capacity: usize = 0,

    work: []f64 = &.{},
    captured_aq: []f64 = &.{},
    aq_ready: bool = false,

    pub fn init(allocator: std.mem.Allocator) SparseForrestTomlin {
        return .{ .allocator = allocator };
    }

    pub fn deinit(self: *SparseForrestTomlin) void {
        self.allocator.free(self.pivot_ids);
        self.allocator.free(self.pivot_values);
        self.allocator.free(self.pivot_lookup);
        self.allocator.free(self.column_head);
        self.allocator.free(self.row_head);
        self.allocator.free(self.entry_row);
        self.allocator.free(self.entry_column);
        self.allocator.free(self.entry_value);
        self.allocator.free(self.column_next);
        self.allocator.free(self.column_previous);
        self.allocator.free(self.row_next);
        self.allocator.free(self.row_previous);
        self.allocator.free(self.free_next);
        self.allocator.free(self.correction_starts);
        self.allocator.free(self.correction_pivots);
        self.allocator.free(self.correction_indices);
        self.allocator.free(self.correction_values);
        self.allocator.free(self.work);
        self.allocator.free(self.captured_aq);
        self.* = .{ .allocator = self.allocator };
    }

    /// Initialize mutable U/UR from packed U columns. `column_rows` contains
    /// logical pivot positions and is remapped to stable basis-column IDs.
    pub fn reset(
        self: *SparseForrestTomlin,
        pivot_ids: []const u32,
        pivot_values: []const f64,
        column_starts: []const usize,
        column_rows: []const u32,
        column_values: []const f64,
    ) FtError!void {
        const n = pivot_ids.len;
        if (pivot_values.len != n or column_starts.len != n + 1 or column_rows.len != column_values.len or column_starts[n] != column_rows.len)
            return error.DimensionMismatch;
        try self.ensureDimension(n);
        try self.ensureLogical(n + 64);
        try self.ensureEntries(@max(column_rows.len + n * 4, n));
        try self.ensureCorrections(@max(n * 4, n));
        try self.ensureCorrectionStarts(66);
        self.dimension = n;
        self.logical_count = n;
        self.entry_high_water = 0;
        self.free_head = none;
        self.update_count = 0;
        self.correction_count = 0;
        self.aq_ready = false;
        @memset(self.row_head[0..n], none);
        @memset(self.column_head[0..self.logical_capacity], none);
        @memset(self.pivot_lookup[0..n], none);
        @memcpy(self.pivot_ids[0..n], pivot_ids);
        @memcpy(self.pivot_values[0..n], pivot_values);
        for (pivot_ids, 0..) |pivot_id, logical| {
            if (pivot_id >= n or self.pivot_lookup[pivot_id] != none) return error.DimensionMismatch;
            self.pivot_lookup[pivot_id] = @intCast(logical);
        }
        for (0..n) |column| {
            for (column_starts[column]..column_starts[column + 1]) |entry| {
                const row_position = column_rows[entry];
                if (row_position >= n) return error.DimensionMismatch;
                try self.insert(pivot_ids[row_position], @intCast(column), column_values[entry]);
            }
        }
        self.correction_starts[0] = 0;
    }

    pub fn mutableUpperView(self: *SparseForrestTomlin) MutableUpperView {
        return .{
            .pivot_ids = self.pivot_ids[0..self.logical_count],
            .pivot_values = self.pivot_values[0..self.logical_count],
            .pivot_lookup = self.pivot_lookup[0..self.dimension],
            .column_heads = self.column_head[0..self.logical_count],
            .row_heads = self.row_head[0..self.dimension],
            .entry_rows = self.entry_row[0..self.entry_high_water],
            .entry_columns = self.entry_column[0..self.entry_high_water],
            .entry_values = self.entry_value[0..self.entry_high_water],
        };
    }

    pub fn hasUpdates(self: *const SparseForrestTomlin) bool {
        return self.update_count != 0;
    }

    pub fn retainedBytes(self: *const SparseForrestTomlin) usize {
        var total: usize = 0;
        total += std.mem.sliceAsBytes(self.pivot_ids).len;
        total += std.mem.sliceAsBytes(self.pivot_values).len;
        total += std.mem.sliceAsBytes(self.pivot_lookup).len;
        total += std.mem.sliceAsBytes(self.column_head).len;
        total += std.mem.sliceAsBytes(self.row_head).len;
        total += std.mem.sliceAsBytes(self.entry_row).len;
        total += std.mem.sliceAsBytes(self.entry_column).len;
        total += std.mem.sliceAsBytes(self.entry_value).len;
        total += std.mem.sliceAsBytes(self.column_next).len;
        total += std.mem.sliceAsBytes(self.column_previous).len;
        total += std.mem.sliceAsBytes(self.row_next).len;
        total += std.mem.sliceAsBytes(self.row_previous).len;
        total += std.mem.sliceAsBytes(self.free_next).len;
        total += std.mem.sliceAsBytes(self.correction_starts).len;
        total += std.mem.sliceAsBytes(self.correction_pivots).len;
        total += std.mem.sliceAsBytes(self.correction_indices).len;
        total += std.mem.sliceAsBytes(self.correction_values).len;
        total += std.mem.sliceAsBytes(self.work).len;
        total += std.mem.sliceAsBytes(self.captured_aq).len;
        return total;
    }

    /// Apply stored FT row corrections and capture the spike immediately
    /// before the mutable upper solve.
    pub fn prepareFtran(self: *SparseForrestTomlin, values: []f64, capture: bool) FtError!void {
        if (values.len != self.dimension) return error.DimensionMismatch;
        for (0..self.correction_count) |correction| {
            const pivot = self.correction_pivots[correction];
            var value = values[pivot];
            for (self.correction_starts[correction]..self.correction_starts[correction + 1]) |entry|
                value -= values[self.correction_indices[entry]] * self.correction_values[entry];
            values[pivot] = value;
        }
        if (capture) {
            @memcpy(self.captured_aq[0..self.dimension], values);
            self.aq_ready = true;
        }
    }

    pub fn solveUpper(self: *SparseForrestTomlin, values: []f64) FtError!void {
        if (values.len != self.dimension) return error.DimensionMismatch;
        var logical = self.logical_count;
        while (logical > 0) {
            logical -= 1;
            const pivot = self.pivot_ids[logical];
            if (pivot == none) continue;
            const diagonal = self.pivot_values[logical];
            if (diagonal == 0.0 or !std.math.isFinite(diagonal)) return error.Singular;
            const solved = values[pivot] / diagonal;
            if (!std.math.isFinite(solved)) return error.NumericalFailure;
            values[pivot] = solved;
            var entry = self.column_head[logical];
            while (entry != none) : (entry = self.column_next[entry])
                values[self.entry_row[entry]] -= self.entry_value[entry] * solved;
        }
    }

    /// U^-T followed by the reverse FT row-correction product. `captureEp`
    /// runs this directly in retained work storage so the next update can
    /// consume the partial BTRAN without another dense copy.
    pub fn solveUpperTranspose(self: *SparseForrestTomlin, values: []f64) FtError!void {
        if (values.len != self.dimension) return error.DimensionMismatch;
        for (0..self.logical_count) |logical| {
            const pivot = self.pivot_ids[logical];
            if (pivot == none) continue;
            var value = values[pivot];
            var entry = self.column_head[logical];
            while (entry != none) : (entry = self.column_next[entry])
                value -= values[self.entry_row[entry]] * self.entry_value[entry];
            value /= self.pivot_values[logical];
            if (!std.math.isFinite(value)) return error.NumericalFailure;
            values[pivot] = value;
        }
        var correction = self.correction_count;
        while (correction > 0) {
            correction -= 1;
            const multiplier = values[self.correction_pivots[correction]];
            for (self.correction_starts[correction]..self.correction_starts[correction + 1]) |entry|
                values[self.correction_indices[entry]] -= multiplier * self.correction_values[entry];
        }
    }

    pub fn captureEp(self: *SparseForrestTomlin, leaving_id: u32) FtError!void {
        if (leaving_id >= self.dimension) return error.DimensionMismatch;
        @memset(self.work[0..self.dimension], 0.0);
        self.work[leaving_id] = 1.0;
        try self.solveUpperTranspose(self.work[0..self.dimension]);
    }

    /// Delete the pivotal U row/column, append the captured FTRAN spike, and
    /// append the FT row correction formed from the captured partial BTRAN.
    pub fn update(self: *SparseForrestTomlin, leaving_id: u32, alpha: f64, zero_tolerance: f64) FtError!void {
        if (leaving_id >= self.dimension or !self.aq_ready or !std.math.isFinite(alpha)) return error.DimensionMismatch;
        const old_logical = self.pivot_lookup[leaving_id];
        if (old_logical == none) return error.Singular;
        const old_pivot = self.pivot_values[old_logical];
        const new_pivot = old_pivot * alpha;
        if (!std.math.isFinite(new_pivot) or @abs(new_pivot) <= zero_tolerance) return error.Singular;

        // Reserve every stream before the first destructive unlink. An OOM
        // therefore leaves the published U/UR and correction chain intact.
        try self.ensureLogical(self.logical_count + 1);
        try self.ensureEntries(self.entry_high_water + self.dimension);
        try self.ensureCorrectionStarts(self.correction_count + 2);
        try self.ensureCorrections(self.correction_starts[self.correction_count] + self.dimension);

        var entry = self.row_head[leaving_id];
        while (entry != none) {
            const next = self.row_next[entry];
            self.remove(entry);
            entry = next;
        }
        entry = self.column_head[old_logical];
        while (entry != none) {
            const next = self.column_next[entry];
            self.remove(entry);
            entry = next;
        }
        self.pivot_ids[old_logical] = none;

        const new_logical: u32 = @intCast(self.logical_count);
        self.logical_count += 1;
        self.pivot_ids[new_logical] = leaving_id;
        self.pivot_values[new_logical] = new_pivot;
        self.column_head[new_logical] = none;
        self.pivot_lookup[leaving_id] = new_logical;
        for (self.captured_aq[0..self.dimension], 0..) |value, row_id| {
            if (row_id == leaving_id or @abs(value) <= zero_tolerance) continue;
            try self.insert(@intCast(row_id), new_logical, value);
        }

        self.correction_pivots[self.correction_count] = leaving_id;
        var correction_output = self.correction_starts[self.correction_count];
        for (self.work[0..self.dimension], 0..) |value, index| {
            if (index == leaving_id or @abs(value) <= zero_tolerance) continue;
            self.correction_indices[correction_output] = @intCast(index);
            self.correction_values[correction_output] = -value * old_pivot;
            correction_output += 1;
        }
        self.correction_count += 1;
        self.correction_starts[self.correction_count] = correction_output;
        self.update_count += 1;
        self.aq_ready = false;
    }

    fn insert(self: *SparseForrestTomlin, row: u32, column: u32, value: f64) FtError!void {
        if (!std.math.isFinite(value)) return error.NumericalFailure;
        if (self.free_head == none and self.entry_high_water == self.entry_capacity) try self.ensureEntries(self.entry_capacity + self.entry_capacity / 2 + 8);
        const entry: u32 = if (self.free_head != none) blk: {
            const result = self.free_head;
            self.free_head = self.free_next[result];
            break :blk result;
        } else blk: {
            const result: u32 = @intCast(self.entry_high_water);
            self.entry_high_water += 1;
            break :blk result;
        };
        self.entry_row[entry] = row;
        self.entry_column[entry] = column;
        self.entry_value[entry] = value;
        self.column_previous[entry] = none;
        self.column_next[entry] = self.column_head[column];
        if (self.column_head[column] != none) self.column_previous[self.column_head[column]] = entry;
        self.column_head[column] = entry;
        self.row_previous[entry] = none;
        self.row_next[entry] = self.row_head[row];
        if (self.row_head[row] != none) self.row_previous[self.row_head[row]] = entry;
        self.row_head[row] = entry;
    }

    fn remove(self: *SparseForrestTomlin, entry: u32) void {
        const row = self.entry_row[entry];
        const column = self.entry_column[entry];
        const cp = self.column_previous[entry];
        const cn = self.column_next[entry];
        if (cp == none) self.column_head[column] = cn else self.column_next[cp] = cn;
        if (cn != none) self.column_previous[cn] = cp;
        const rp = self.row_previous[entry];
        const rn = self.row_next[entry];
        if (rp == none) self.row_head[row] = rn else self.row_next[rp] = rn;
        if (rn != none) self.row_previous[rn] = rp;
        self.free_next[entry] = self.free_head;
        self.free_head = entry;
    }

    fn ensureDimension(self: *SparseForrestTomlin, required: usize) FtError!void {
        if (required <= self.work.len) return;
        const capacity = grow(self.work.len, required) catch return error.CapacityOverflow;
        try self.resizeRetained(&self.pivot_lookup, capacity);
        try self.resizeRetained(&self.row_head, capacity);
        try self.resizeRetained(&self.work, capacity);
        try self.resizeRetained(&self.captured_aq, capacity);
    }
    fn ensureLogical(self: *SparseForrestTomlin, required: usize) FtError!void {
        if (required <= self.logical_capacity) return;
        const capacity = grow(self.logical_capacity, required) catch return error.CapacityOverflow;
        try self.resizeRetained(&self.pivot_ids, capacity);
        try self.resizeRetained(&self.pivot_values, capacity);
        try self.resizeRetained(&self.column_head, capacity);
        self.logical_capacity = capacity;
    }
    fn ensureEntries(self: *SparseForrestTomlin, required: usize) FtError!void {
        if (required <= self.entry_capacity) return;
        const capacity = grow(self.entry_capacity, required) catch return error.CapacityOverflow;
        try self.resizeRetained(&self.entry_row, capacity);
        try self.resizeRetained(&self.entry_column, capacity);
        try self.resizeRetained(&self.entry_value, capacity);
        try self.resizeRetained(&self.column_next, capacity);
        try self.resizeRetained(&self.column_previous, capacity);
        try self.resizeRetained(&self.row_next, capacity);
        try self.resizeRetained(&self.row_previous, capacity);
        try self.resizeRetained(&self.free_next, capacity);
        self.entry_capacity = capacity;
    }

    fn resizeRetained(self: *SparseForrestTomlin, slice: anytype, capacity: usize) FtError!void {
        slice.* = self.allocator.realloc(slice.*, capacity) catch return error.OutOfMemory;
    }
    fn ensureCorrections(self: *SparseForrestTomlin, required: usize) FtError!void {
        if (required <= self.correction_capacity) return;
        const capacity = grow(self.correction_capacity, required) catch return error.CapacityOverflow;
        self.correction_indices = self.allocator.realloc(self.correction_indices, capacity) catch return error.OutOfMemory;
        self.correction_values = self.allocator.realloc(self.correction_values, capacity) catch return error.OutOfMemory;
        self.correction_capacity = capacity;
    }
    fn ensureCorrectionStarts(self: *SparseForrestTomlin, required: usize) FtError!void {
        if (required <= self.correction_starts.len) return;
        const capacity = grow(self.correction_starts.len, required) catch return error.CapacityOverflow;
        self.correction_starts = self.allocator.realloc(self.correction_starts, capacity) catch return error.OutOfMemory;
        self.correction_pivots = self.allocator.realloc(self.correction_pivots, capacity) catch return error.OutOfMemory;
    }
};

fn grow(current: usize, required: usize) error{Overflow}!usize {
    var capacity = @max(current, 8);
    while (capacity < required) capacity = std.math.add(usize, capacity, capacity / 2 + 8) catch return error.Overflow;
    if (capacity > std.math.maxInt(u32)) return error.Overflow;
    return capacity;
}
