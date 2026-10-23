//! Packed mutable Forrest--Tomlin upper factor.
//!
//! U columns are append-only contiguous segments. Reciprocal row membership
//! lives in growable contiguous slabs and is used only by UPDATE deletion;
//! both FTRAN and BTRAN retain the legacy column traversal and floating-point
//! order. This mirrors HiGHS' U/UR ownership without per-entry linked lists.

const std = @import("std");
const legacy = @import("sparse_ft.zig");

const none = std.math.maxInt(u32);
const highs_tiny: f64 = 1e-14;

pub const FtError = legacy.FtError;
pub const MutableUpperView = legacy.MutableUpperView;

pub const SparseForrestTomlin = struct {
    allocator: std.mem.Allocator,
    dimension: usize = 0,
    logical_count: usize = 0,
    logical_capacity: usize = 0,
    update_count: usize = 0,

    pivot_ids: []u32 = &.{},
    pivot_values: []f64 = &.{},
    pivot_lookup: []u32 = &.{},

    column_start: []usize = &.{},
    column_end: []usize = &.{},
    column_rows: []u32 = &.{},
    column_values: []f64 = &.{},
    column_high_water: usize = 0,
    column_live_count: usize = 0,
    column_capacity: usize = 0,

    row_start: []usize = &.{},
    row_end: []usize = &.{},
    row_space: []usize = &.{},
    row_columns: []u32 = &.{},
    row_high_water: usize = 0,
    row_capacity: usize = 0,

    correction_starts: []usize = &.{},
    correction_pivots: []u32 = &.{},
    correction_indices: []u32 = &.{},
    correction_values: []f64 = &.{},
    correction_count: usize = 0,
    correction_capacity: usize = 0,

    work: []f64 = &.{},
    ep_touched: []u32 = &.{},
    ep_touched_count: usize = 0,
    captured_aq: []f64 = &.{},
    aq_ready: bool = false,
    ftran_output: []u32 = &.{},
    ftran_output_count: usize = 0,

    pub fn init(allocator: std.mem.Allocator) SparseForrestTomlin {
        return .{ .allocator = allocator };
    }

    pub fn deinit(self: *SparseForrestTomlin) void {
        self.allocator.free(self.pivot_ids);
        self.allocator.free(self.pivot_values);
        self.allocator.free(self.pivot_lookup);
        self.allocator.free(self.column_start);
        self.allocator.free(self.column_end);
        self.allocator.free(self.column_rows);
        self.allocator.free(self.column_values);
        self.allocator.free(self.row_start);
        self.allocator.free(self.row_end);
        self.allocator.free(self.row_space);
        self.allocator.free(self.row_columns);
        self.allocator.free(self.correction_starts);
        self.allocator.free(self.correction_pivots);
        self.allocator.free(self.correction_indices);
        self.allocator.free(self.correction_values);
        self.allocator.free(self.work);
        self.allocator.free(self.ep_touched);
        self.allocator.free(self.captured_aq);
        self.allocator.free(self.ftran_output);
        self.* = init(self.allocator);
    }

    pub fn reset(
        self: *SparseForrestTomlin,
        pivot_ids: []const u32,
        pivot_values: []const f64,
        column_starts: []const usize,
        column_rows: []const u32,
        column_values: []const f64,
    ) FtError!void {
        const n = pivot_ids.len;
        if (pivot_values.len != n or column_starts.len != n + 1 or
            column_rows.len != column_values.len or column_starts[n] != column_rows.len)
            return error.DimensionMismatch;
        try self.ensureDimension(n);
        try self.ensureLogical(n + 64);
        try self.ensureColumnEntries(@max(column_rows.len + n * 4, n));
        try self.ensureRowEntries(@max(column_rows.len + n * 5, n));
        try self.ensureCorrections(@max(n * 4, n));
        try self.ensureCorrectionStarts(66);

        self.dimension = n;
        self.logical_count = n;
        self.update_count = 0;
        self.column_high_water = column_rows.len;
        self.column_live_count = column_rows.len;
        self.row_high_water = 0;
        self.correction_count = 0;
        self.ep_touched_count = 0;
        self.aq_ready = false;
        self.ftran_output_count = 0;
        @memset(self.work[0..n], 0.0);
        @memset(self.pivot_lookup[0..n], none);
        @memcpy(self.pivot_ids[0..n], pivot_ids);
        @memcpy(self.pivot_values[0..n], pivot_values);

        for (pivot_ids, pivot_values, 0..) |pivot_id, pivot_value, logical| {
            if (pivot_id >= n or self.pivot_lookup[pivot_id] != none) return error.DimensionMismatch;
            if (pivot_value == 0.0 or !std.math.isFinite(pivot_value)) return error.Singular;
            self.pivot_lookup[pivot_id] = @intCast(logical);
            self.column_start[logical] = column_starts[logical];
            self.column_end[logical] = column_starts[logical + 1];
        }
        for (column_rows, column_values, 0..) |row_position, value, entry| {
            if (row_position >= n or !std.math.isFinite(value)) return error.DimensionMismatch;
            self.column_rows[entry] = pivot_ids[row_position];
            self.column_values[entry] = value;
        }

        @memset(self.row_space[0..n], 0);
        for (column_rows) |row_position| self.row_space[pivot_ids[row_position]] += 1;
        for (0..n) |row_position| {
            const pivot_id = pivot_ids[row_position];
            const count = self.row_space[pivot_id];
            self.row_start[pivot_id] = self.row_high_water;
            self.row_end[pivot_id] = self.row_high_water;
            self.row_space[pivot_id] = count + 5;
            self.row_high_water += count + 5;
        }
        if (self.row_high_water > self.row_capacity) return error.CapacityOverflow;
        for (0..n) |logical| {
            for (column_starts[logical]..column_starts[logical + 1]) |entry| {
                const row = self.column_rows[entry];
                const destination = self.row_end[row];
                self.row_columns[destination] = @intCast(logical);
                self.row_end[row] += 1;
                self.row_space[row] -= 1;
            }
        }
        self.correction_starts[0] = 0;
    }

    pub fn mutableUpperView(self: *SparseForrestTomlin) MutableUpperView {
        return .{
            .pivot_ids = self.pivot_ids[0..self.logical_count],
            .pivot_values = self.pivot_values[0..self.logical_count],
            .pivot_lookup = self.pivot_lookup[0..self.dimension],
            .column_heads = &.{},
            .row_heads = &.{},
            .entry_rows = self.column_rows[0..self.column_high_water],
            .entry_columns = &.{},
            .entry_values = self.column_values[0..self.column_high_water],
        };
    }

    pub fn hasUpdates(self: *const SparseForrestTomlin) bool {
        return self.update_count != 0;
    }

    pub fn retainedBytes(self: *const SparseForrestTomlin) usize {
        var total: usize = 0;
        inline for (.{
            self.pivot_ids,
            self.pivot_values,
            self.pivot_lookup,
            self.column_start,
            self.column_end,
            self.column_rows,
            self.column_values,
            self.row_start,
            self.row_end,
            self.row_space,
            self.row_columns,
            self.correction_starts,
            self.correction_pivots,
            self.correction_indices,
            self.correction_values,
            self.work,
            self.ep_touched,
            self.captured_aq,
            self.ftran_output,
        }) |slice| total += std.mem.sliceAsBytes(slice).len;
        return total;
    }

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
        self.ftran_output_count = 0;
        var logical = self.logical_count;
        while (logical > 0) {
            logical -= 1;
            const pivot = self.pivot_ids[logical];
            if (pivot == none) continue;
            const pivot_multiplier = values[pivot];
            if (@abs(pivot_multiplier) <= highs_tiny) {
                values[pivot] = 0.0;
                continue;
            }
            const solved = pivot_multiplier / self.pivot_values[logical];
            if (!std.math.isFinite(solved)) return error.NumericalFailure;
            self.ftran_output[self.ftran_output_count] = pivot;
            self.ftran_output_count += 1;
            values[pivot] = solved;
            var entry = self.column_end[logical];
            while (entry > self.column_start[logical]) {
                entry -= 1;
                values[self.column_rows[entry]] -= self.column_values[entry] * solved;
            }
        }
    }

    pub fn lastFtranOutputIndices(self: *const SparseForrestTomlin) []const u32 {
        return self.ftran_output[0..self.ftran_output_count];
    }

    pub fn solveUpperTranspose(self: *SparseForrestTomlin, values: []f64) FtError!void {
        if (values.len != self.dimension) return error.DimensionMismatch;
        try self.solveMutableUpperTransposeFrom(values, 0);
        self.applyCorrectionsTranspose(values);
    }

    fn solveMutableUpperTransposeFrom(self: *SparseForrestTomlin, values: []f64, first_logical: usize) FtError!void {
        for (first_logical..self.logical_count) |logical| {
            const pivot = self.pivot_ids[logical];
            if (pivot == none) continue;
            var value = values[pivot];
            var entry = self.column_end[logical];
            while (entry > self.column_start[logical]) {
                entry -= 1;
                value -= values[self.column_rows[entry]] * self.column_values[entry];
            }
            value /= self.pivot_values[logical];
            if (!std.math.isFinite(value)) return error.NumericalFailure;
            values[pivot] = value;
        }
    }

    fn applyCorrectionsTranspose(self: *const SparseForrestTomlin, values: []f64) void {
        var correction = self.correction_count;
        while (correction > 0) {
            correction -= 1;
            const multiplier = values[self.correction_pivots[correction]];
            if (multiplier == 0.0) continue;
            for (self.correction_starts[correction]..self.correction_starts[correction + 1]) |entry|
                values[self.correction_indices[entry]] -= multiplier * self.correction_values[entry];
        }
    }

    pub fn captureEp(self: *SparseForrestTomlin, leaving_id: u32) FtError!void {
        if (leaving_id >= self.dimension) return error.DimensionMismatch;
        const first_logical = self.pivot_lookup[leaving_id];
        if (first_logical == none) return error.Singular;
        for (self.ep_touched[0..self.ep_touched_count]) |pivot| self.work[pivot] = 0.0;
        self.ep_touched_count = 0;
        self.work[leaving_id] = 1.0;
        for (first_logical..self.logical_count) |logical| {
            const pivot = self.pivot_ids[logical];
            if (pivot == none) continue;
            var value = self.work[pivot];
            var entry = self.column_end[logical];
            while (entry > self.column_start[logical]) {
                entry -= 1;
                value -= self.work[self.column_rows[entry]] * self.column_values[entry];
            }
            value /= self.pivot_values[logical];
            if (!std.math.isFinite(value)) return error.NumericalFailure;
            self.work[pivot] = value;
            self.ep_touched[self.ep_touched_count] = pivot;
            self.ep_touched_count += 1;
        }
    }

    pub fn update(self: *SparseForrestTomlin, leaving_id: u32, alpha: f64, zero_tolerance: f64) FtError!void {
        if (leaving_id >= self.dimension or !self.aq_ready or !std.math.isFinite(alpha)) return error.DimensionMismatch;
        const old_logical = self.pivot_lookup[leaving_id];
        if (old_logical == none) return error.Singular;
        const old_pivot = self.pivot_values[old_logical];
        const new_pivot = old_pivot * alpha;
        if (!std.math.isFinite(new_pivot) or @abs(new_pivot) <= zero_tolerance) return error.Singular;

        try self.ensureLogical(self.logical_count + 1);
        try self.ensureColumnEntries(std.math.add(usize, self.column_high_water, self.dimension) catch return error.CapacityOverflow);
        try self.ensureCorrectionStarts(self.correction_count + 2);
        try self.ensureCorrections(self.correction_starts[self.correction_count] + self.dimension);

        // Build the unpublished replacement column and determine every UR
        // relocation before mutating the live U. This keeps UPDATE
        // transactional while scanning the dense spike only once.
        const new_column_start = self.column_high_water;
        var new_column_end = new_column_start;
        var relocation_entries: usize = 0;
        for (self.captured_aq[0..self.dimension], 0..) |value, row_usize| {
            if (row_usize == leaving_id or @abs(value) <= zero_tolerance) continue;
            const row: u32 = @intCast(row_usize);
            self.column_rows[new_column_end] = row;
            self.column_values[new_column_end] = value;
            new_column_end += 1;
            if (self.row_space[row] == 0) {
                const count = self.row_end[row] - self.row_start[row];
                relocation_entries = std.math.add(usize, relocation_entries, count + count / 10 + 5) catch return error.CapacityOverflow;
            }
        }
        try self.ensureRowEntries(std.math.add(usize, self.row_high_water, relocation_entries) catch return error.CapacityOverflow);

        // Remove the pivotal stable row from every U column that contains it.
        const leaving_row_start = self.row_start[leaving_id];
        const leaving_row_end = self.row_end[leaving_id];
        for (self.row_columns[leaving_row_start..leaving_row_end]) |logical|
            self.removeColumnRow(logical, leaving_id);
        self.row_end[leaving_id] = leaving_row_start;
        self.row_space[leaving_id] += leaving_row_end - leaving_row_start;

        // Remove the old pivotal logical column from every reciprocal row.
        for (self.column_rows[self.column_start[old_logical]..self.column_end[old_logical]]) |row|
            self.removeRowColumn(row, old_logical);
        self.column_live_count -= self.column_end[old_logical] - self.column_start[old_logical];
        self.column_end[old_logical] = self.column_start[old_logical];
        self.pivot_ids[old_logical] = none;

        const new_logical: u32 = @intCast(self.logical_count);
        self.logical_count += 1;
        self.pivot_ids[new_logical] = leaving_id;
        self.pivot_values[new_logical] = new_pivot;
        self.pivot_lookup[leaving_id] = new_logical;
        self.column_start[new_logical] = new_column_start;
        self.column_end[new_logical] = new_column_end;
        self.column_high_water = new_column_end;
        self.column_live_count += new_column_end - new_column_start;
        for (self.column_rows[new_column_start..new_column_end]) |row| self.appendRowColumn(row, new_logical);

        self.correction_pivots[self.correction_count] = leaving_id;
        var correction_output = self.correction_starts[self.correction_count];
        for (self.ep_touched[0..self.ep_touched_count]) |index| {
            if (index == leaving_id) continue;
            const value = self.work[index];
            if (@abs(value) <= zero_tolerance) continue;
            self.correction_indices[correction_output] = index;
            self.correction_values[correction_output] = -value * old_pivot;
            correction_output += 1;
        }
        self.correction_count += 1;
        self.correction_starts[self.correction_count] = correction_output;
        self.update_count += 1;
        self.aq_ready = false;
        if (self.dimension >= 128 and self.column_high_water > self.column_live_count + self.column_live_count / 2) self.compactColumns();
    }

    fn removeColumnRow(self: *SparseForrestTomlin, logical: u32, row: u32) void {
        const start = self.column_start[logical];
        var end = self.column_end[logical];
        var found = start;
        while (found < end and self.column_rows[found] != row) : (found += 1) {}
        std.debug.assert(found < end);
        end -= 1;
        var position = found;
        while (position < end) : (position += 1) {
            self.column_rows[position] = self.column_rows[position + 1];
            self.column_values[position] = self.column_values[position + 1];
        }
        self.column_end[logical] = end;
        self.column_live_count -= 1;
    }

    /// Reclaim abandoned column segments once they materially enlarge the U
    /// working set. Logical columns are created in append order, so every
    /// live source segment is at or above its compacted destination and a
    /// forward in-place copy is safe.
    noinline fn compactColumns(self: *SparseForrestTomlin) void {
        var output: usize = 0;
        for (0..self.logical_count) |logical| {
            if (self.pivot_ids[logical] == none) continue;
            const start = self.column_start[logical];
            const count = self.column_end[logical] - start;
            for (0..count) |offset| {
                self.column_rows[output + offset] = self.column_rows[start + offset];
                self.column_values[output + offset] = self.column_values[start + offset];
            }
            self.column_start[logical] = output;
            self.column_end[logical] = output + count;
            output += count;
        }
        std.debug.assert(output == self.column_live_count);
        self.column_high_water = output;
    }

    fn removeRowColumn(self: *SparseForrestTomlin, row: u32, logical: u32) void {
        const start = self.row_start[row];
        var end = self.row_end[row];
        var found = start;
        while (found < end and self.row_columns[found] != logical) : (found += 1) {}
        std.debug.assert(found < end);
        end -= 1;
        var position = found;
        while (position < end) : (position += 1)
            self.row_columns[position] = self.row_columns[position + 1];
        self.row_end[row] = end;
        self.row_space[row] += 1;
    }

    fn appendRowColumn(self: *SparseForrestTomlin, row: u32, logical: u32) void {
        if (self.row_space[row] == 0) self.relocateRow(row);
        self.row_columns[self.row_end[row]] = logical;
        self.row_end[row] += 1;
        self.row_space[row] -= 1;
    }

    fn relocateRow(self: *SparseForrestTomlin, row: u32) void {
        const old_start = self.row_start[row];
        const count = self.row_end[row] - old_start;
        const spare = count / 10 + 5;
        const new_start = self.row_high_water;
        std.debug.assert(new_start + count + spare <= self.row_capacity);
        @memcpy(self.row_columns[new_start .. new_start + count], self.row_columns[old_start .. old_start + count]);
        self.row_start[row] = new_start;
        self.row_end[row] = new_start + count;
        self.row_space[row] = spare;
        self.row_high_water += count + spare;
    }

    pub fn validate(self: *const SparseForrestTomlin) bool {
        if (self.dimension == 0) return true;
        var live: usize = 0;
        for (self.pivot_ids[0..self.logical_count], 0..) |pivot, logical| {
            if (pivot == none) continue;
            if (pivot >= self.dimension or self.pivot_lookup[pivot] != logical) return false;
            if (self.column_start[logical] > self.column_end[logical] or self.column_end[logical] > self.column_high_water) return false;
            live += 1;
            for (self.column_rows[self.column_start[logical]..self.column_end[logical]]) |row| {
                if (row >= self.dimension or self.row_start[row] > self.row_end[row] or self.row_end[row] > self.row_high_water) return false;
                var reciprocal_count: usize = 0;
                for (self.row_columns[self.row_start[row]..self.row_end[row]]) |column| {
                    if (column == logical) reciprocal_count += 1;
                }
                if (reciprocal_count != 1) return false;
            }
        }
        if (live != self.dimension) return false;
        var column_live_count: usize = 0;
        for (0..self.logical_count) |logical|
            column_live_count += self.column_end[logical] - self.column_start[logical];
        if (column_live_count != self.column_live_count) return false;
        for (0..self.dimension) |row| {
            if (self.row_start[row] > self.row_end[row] or self.row_end[row] > self.row_high_water) return false;
            for (self.row_columns[self.row_start[row]..self.row_end[row]]) |logical| {
                if (logical >= self.logical_count or self.pivot_ids[logical] == none) return false;
                var column_count: usize = 0;
                for (self.column_rows[self.column_start[logical]..self.column_end[logical]]) |column_row| {
                    if (column_row == row) column_count += 1;
                }
                if (column_count != 1) return false;
            }
        }
        return true;
    }

    fn ensureDimension(self: *SparseForrestTomlin, required: usize) FtError!void {
        if (required <= self.work.len) return;
        const capacity = grow(self.work.len, required) catch return error.CapacityOverflow;
        try self.resizeRetained(&self.pivot_lookup, capacity);
        try self.resizeRetained(&self.row_start, capacity);
        try self.resizeRetained(&self.row_end, capacity);
        try self.resizeRetained(&self.row_space, capacity);
        try self.resizeRetained(&self.work, capacity);
        try self.resizeRetained(&self.ep_touched, capacity);
        try self.resizeRetained(&self.captured_aq, capacity);
        try self.resizeRetained(&self.ftran_output, capacity);
    }

    fn ensureLogical(self: *SparseForrestTomlin, required: usize) FtError!void {
        if (required <= self.logical_capacity) return;
        const capacity = grow(self.logical_capacity, required) catch return error.CapacityOverflow;
        try self.resizeRetained(&self.pivot_ids, capacity);
        try self.resizeRetained(&self.pivot_values, capacity);
        try self.resizeRetained(&self.column_start, capacity);
        try self.resizeRetained(&self.column_end, capacity);
        self.logical_capacity = capacity;
    }

    fn ensureColumnEntries(self: *SparseForrestTomlin, required: usize) FtError!void {
        if (required <= self.column_capacity) return;
        const capacity = grow(self.column_capacity, required) catch return error.CapacityOverflow;
        try self.resizeRetained(&self.column_rows, capacity);
        try self.resizeRetained(&self.column_values, capacity);
        self.column_capacity = capacity;
    }

    fn ensureRowEntries(self: *SparseForrestTomlin, required: usize) FtError!void {
        if (required <= self.row_capacity) return;
        const capacity = grow(self.row_capacity, required) catch return error.CapacityOverflow;
        try self.resizeRetained(&self.row_columns, capacity);
        self.row_capacity = capacity;
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

    fn resizeRetained(self: *SparseForrestTomlin, slice: anytype, capacity: usize) FtError!void {
        slice.* = self.allocator.realloc(slice.*, capacity) catch return error.OutOfMemory;
    }
};

fn grow(current: usize, required: usize) error{Overflow}!usize {
    var capacity = @max(current, 8);
    while (capacity < required) capacity = std.math.add(usize, capacity, capacity / 2 + 8) catch return error.Overflow;
    if (capacity > std.math.maxInt(u32)) return error.Overflow;
    return capacity;
}

test "packed FT stays byte-equivalent to intrusive FT across updates" {
    const pivots = [_]u32{ 0, 1, 2, 3 };
    const diagonal = [_]f64{ 2, 3, 4, 5 };
    const starts = [_]usize{ 0, 0, 1, 3, 5 };
    const rows = [_]u32{ 0, 0, 1, 1, 2 };
    const values = [_]f64{ 0.25, -0.5, 0.125, 0.375, -0.25 };

    var packed_ft = SparseForrestTomlin.init(std.testing.allocator);
    defer packed_ft.deinit();
    var linked = legacy.SparseForrestTomlin.init(std.testing.allocator);
    defer linked.deinit();
    try packed_ft.reset(&pivots, &diagonal, &starts, &rows, &values);
    try linked.reset(&pivots, &diagonal, &starts, &rows, &values);
    try std.testing.expect(packed_ft.validate());

    const leaving = [_]u32{ 1, 3, 0, 2 };
    const raw = [_][4]f64{
        .{ 1, 4, -1, 2 },
        .{ 2, -1, 3, 5 },
        .{ 6, 1, -2, 1 },
        .{ -1, 2, 7, 3 },
    };
    for (leaving, raw) |leaving_id, input| {
        var packed_aq = input;
        var linked_aq = input;
        try packed_ft.prepareFtran(&packed_aq, true);
        try linked.prepareFtran(&linked_aq, true);
        try packed_ft.solveUpper(&packed_aq);
        try linked.solveUpper(&linked_aq);
        try std.testing.expectEqualSlices(f64, &linked_aq, &packed_aq);
        try std.testing.expectEqualSlices(u32, linked.lastFtranOutputIndices(), packed_ft.lastFtranOutputIndices());
        try packed_ft.captureEp(leaving_id);
        try linked.captureEp(leaving_id);
        try packed_ft.update(leaving_id, packed_aq[leaving_id], highs_tiny);
        try linked.update(leaving_id, linked_aq[leaving_id], highs_tiny);
        try std.testing.expect(packed_ft.validate());

        var packed_rhs = [_]f64{ 3, -2, 5, 1 };
        var linked_rhs = packed_rhs;
        try packed_ft.prepareFtran(&packed_rhs, false);
        try linked.prepareFtran(&linked_rhs, false);
        try std.testing.expectEqualSlices(f64, &linked_rhs, &packed_rhs);
        try packed_ft.solveUpper(&packed_rhs);
        try linked.solveUpper(&linked_rhs);
        try std.testing.expectEqualSlices(f64, &linked_rhs, &packed_rhs);
        packed_rhs = .{ -1, 4, 2, 6 };
        linked_rhs = packed_rhs;
        try packed_ft.solveUpperTranspose(&packed_rhs);
        try linked.solveUpperTranspose(&linked_rhs);
        try std.testing.expectEqualSlices(f64, &linked_rhs, &packed_rhs);
    }
}

test "packed FT randomized differential preserves solves and reciprocal structure" {
    const n = 12;
    var prng = std.Random.DefaultPrng.init(0x7061_636b_6564_4654);
    const random = prng.random();

    var pivots: [n]u32 = undefined;
    var diagonal: [n]f64 = undefined;
    var starts: [n + 1]usize = undefined;
    var rows: [n * (n - 1) / 2]u32 = undefined;
    var coefficients: [rows.len]f64 = undefined;
    for (0..n) |index| {
        pivots[index] = @intCast(index);
        diagonal[index] = @floatFromInt(3 + index % 5);
    }
    var entry_count: usize = 0;
    for (0..n) |column| {
        starts[column] = entry_count;
        for (0..column) |row| {
            if (random.intRangeLessThan(u8, 0, 4) == 0) continue;
            rows[entry_count] = @intCast(row);
            const integer = random.intRangeAtMost(i8, -4, 4);
            coefficients[entry_count] = if (integer == 0) 0.125 else @as(f64, @floatFromInt(integer)) * 0.125;
            entry_count += 1;
        }
    }
    starts[n] = entry_count;

    var packed_ft = SparseForrestTomlin.init(std.testing.allocator);
    defer packed_ft.deinit();
    var linked = legacy.SparseForrestTomlin.init(std.testing.allocator);
    defer linked.deinit();
    try packed_ft.reset(&pivots, &diagonal, &starts, rows[0..entry_count], coefficients[0..entry_count]);
    try linked.reset(&pivots, &diagonal, &starts, rows[0..entry_count], coefficients[0..entry_count]);

    for (0..32) |update_index| {
        const leaving_id: u32 = @intCast((update_index * 7 + 3) % n);
        var packed_aq: [n]f64 = undefined;
        var linked_aq: [n]f64 = undefined;
        var usable = false;
        for (0..16) |_| {
            for (&packed_aq) |*value| {
                const integer = random.intRangeAtMost(i8, -8, 8);
                value.* = @floatFromInt(integer);
            }
            linked_aq = packed_aq;
            try packed_ft.prepareFtran(&packed_aq, true);
            try linked.prepareFtran(&linked_aq, true);
            try std.testing.expectEqualSlices(f64, &linked_aq, &packed_aq);
            try packed_ft.solveUpper(&packed_aq);
            try linked.solveUpper(&linked_aq);
            try std.testing.expectEqualSlices(f64, &linked_aq, &packed_aq);
            try std.testing.expectEqualSlices(u32, linked.lastFtranOutputIndices(), packed_ft.lastFtranOutputIndices());
            const alpha = packed_aq[leaving_id];
            if (std.math.isFinite(alpha) and @abs(alpha) >= 0.01 and @abs(alpha) <= 100.0) {
                usable = true;
                break;
            }
        }
        try std.testing.expect(usable);
        try packed_ft.captureEp(leaving_id);
        try linked.captureEp(leaving_id);
        try packed_ft.update(leaving_id, packed_aq[leaving_id], highs_tiny);
        try linked.update(leaving_id, linked_aq[leaving_id], highs_tiny);
        try std.testing.expect(packed_ft.validate());

        var packed_rhs: [n]f64 = undefined;
        for (&packed_rhs) |*value| {
            const integer = random.intRangeAtMost(i8, -8, 8);
            value.* = @floatFromInt(integer);
        }
        var linked_rhs = packed_rhs;
        try packed_ft.prepareFtran(&packed_rhs, false);
        try linked.prepareFtran(&linked_rhs, false);
        try std.testing.expectEqualSlices(f64, &linked_rhs, &packed_rhs);
        try packed_ft.solveUpper(&packed_rhs);
        try linked.solveUpper(&linked_rhs);
        try std.testing.expectEqualSlices(f64, &linked_rhs, &packed_rhs);

        for (&packed_rhs) |*value| {
            const integer = random.intRangeAtMost(i8, -8, 8);
            value.* = @floatFromInt(integer);
        }
        linked_rhs = packed_rhs;
        try packed_ft.solveUpperTranspose(&packed_rhs);
        try linked.solveUpperTranspose(&linked_rhs);
        try std.testing.expectEqualSlices(f64, &linked_rhs, &packed_rhs);
    }
}
