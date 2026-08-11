//! Mutable Forrest--Tomlin upper factor and row-correction workspace.
//!
//! Pivot identifiers are simplex basis-column IDs. U is stored simultaneously
//! by column and by row (UR) through one intrusive SoA entry pool, allowing a
//! pivotal row and column to be removed without rebuilding unrelated segments.

const std = @import("std");

const none = std.math.maxInt(u32);
const highs_tiny: f64 = 1e-14;

pub const FtError = error{ DimensionMismatch, Singular, NumericalFailure, OutOfMemory, CapacityOverflow };

pub const MutableUpperView = struct {
    /// Stable basis-column ID at each logical U position.
    pivot_ids: []const u32,
    /// Mutable U diagonal parallel to `pivot_ids`.
    pivot_values: []f64,
    /// Map from stable pivot ID to current logical position.
    pivot_lookup: []const u32,
    /// Intrusive entry-list head for each logical U column.
    column_heads: []const u32,
    /// Intrusive entry-list head for each stable pivot row.
    row_heads: []const u32,
    /// Stable row ID for every allocated entry slot.
    entry_rows: []const u32,
    /// Logical column ID for every allocated entry slot.
    entry_columns: []const u32,
    /// Mutable coefficient for every allocated entry slot.
    entry_values: []f64,
};

pub const SparseForrestTomlin = struct {
    /// Owner of mutable U, correction, and solve-work storage.
    allocator: std.mem.Allocator,
    /// Fixed number of stable basis pivot IDs.
    dimension: usize = 0,
    /// Current logical U column count, including update spikes.
    logical_count: usize = 0,
    /// Allocated capacity for logical U columns.
    logical_capacity: usize = 0,
    /// First never-used intrusive entry slot.
    entry_high_water: usize = 0,
    /// Allocated intrusive-entry capacity.
    entry_capacity: usize = 0,
    /// Head of the recycled-entry free list.
    free_head: u32 = none,
    /// Number of installed Forrest--Tomlin updates.
    update_count: usize = 0,

    /// Stable pivot ID at each logical column position.
    pivot_ids: []u32 = &.{},
    /// Mutable diagonal value at each logical column position.
    pivot_values: []f64 = &.{},
    /// Stable pivot ID to logical-position lookup.
    pivot_lookup: []u32 = &.{},
    /// Intrusive list heads by logical column.
    column_head: []u32 = &.{},
    /// Intrusive list heads by stable row.
    row_head: []u32 = &.{},

    /// Stable row ID per entry slot.
    entry_row: []u32 = &.{},
    /// Logical column ID per entry slot.
    entry_column: []u32 = &.{},
    /// Numerical coefficient per entry slot.
    entry_value: []f64 = &.{},
    /// Next entry in the same column.
    column_next: []u32 = &.{},
    /// Previous entry in the same column.
    column_previous: []u32 = &.{},
    /// Next entry in the same row.
    row_next: []u32 = &.{},
    /// Previous entry in the same row.
    row_previous: []u32 = &.{},
    /// Next recycled entry slot.
    free_next: []u32 = &.{},

    /// Prefix offsets for stored row-correction entries.
    correction_starts: []usize = &.{},
    /// Stable pivotal row for each correction.
    correction_pivots: []u32 = &.{},
    /// Source row IDs for correction coefficients.
    correction_indices: []u32 = &.{},
    /// Correction coefficients parallel to `correction_indices`.
    correction_values: []f64 = &.{},
    /// Number of installed row corrections.
    correction_count: usize = 0,
    /// Allocated correction-entry capacity.
    correction_capacity: usize = 0,

    /// Reusable dense solve workspace.
    work: []f64 = &.{},
    /// Stable IDs written by the previous `captureEp` solve. Tracking the
    /// producer-owned support permits sparse clearing without relying on a
    /// caller-maintained zero complement.
    ep_touched: []u32 = &.{},
    ep_touched_count: usize = 0,
    /// Captured pre-upper-solve spike used to construct an update.
    captured_aq: []f64 = &.{},
    /// Whether `captured_aq` contains a valid pending spike.
    aq_ready: bool = false,
    /// Stable pivot IDs emitted by the latest forward upper solve, in the
    /// reverse logical traversal order consumed by simplex HVector clients.
    ftran_output: []u32 = &.{},
    ftran_output_count: usize = 0,

    /// Construct an empty lazily allocated FT update store.
    pub fn init(allocator: std.mem.Allocator) SparseForrestTomlin {
        return .{ .allocator = allocator };
    }

    /// Release all mutable-U, correction, and solve buffers.
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
        self.allocator.free(self.ep_touched);
        self.allocator.free(self.captured_aq);
        self.allocator.free(self.ftran_output);
        self.* = init(self.allocator);
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
        self.ep_touched_count = 0;
        self.aq_ready = false;
        self.ftran_output_count = 0;
        @memset(self.work[0..n], 0.0);
        @memset(self.row_head[0..n], none);
        @memset(self.column_head[0..self.logical_capacity], none);
        @memset(self.pivot_lookup[0..n], none);
        @memcpy(self.pivot_ids[0..n], pivot_ids);
        @memcpy(self.pivot_values[0..n], pivot_values);
        for (pivot_ids, pivot_values, 0..) |pivot_id, pivot_value, logical| {
            if (pivot_id >= n or self.pivot_lookup[pivot_id] != none) return error.DimensionMismatch;
            if (pivot_value == 0.0 or !std.math.isFinite(pivot_value)) return error.Singular;
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

    /// Borrow mutable U/UR arrays for sparse LU solve integration.
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

    /// Return whether the base factorization has a nonempty FT update chain.
    pub fn hasUpdates(self: *const SparseForrestTomlin) bool {
        return self.update_count != 0;
    }

    /// Sum bytes represented by all retained slices.
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
        total += std.mem.sliceAsBytes(self.ep_touched).len;
        total += std.mem.sliceAsBytes(self.captured_aq).len;
        total += std.mem.sliceAsBytes(self.ftran_output).len;
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

    /// Solve the current mutable upper triangular system in place.
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
            var entry = self.column_head[logical];
            while (entry != none) : (entry = self.column_next[entry])
                values[self.entry_row[entry]] -= self.entry_value[entry] * solved;
        }
    }

    /// Producer-owned HVector support from the latest forward upper solve.
    pub fn lastFtranOutputIndices(self: *const SparseForrestTomlin) []const u32 {
        return self.ftran_output[0..self.ftran_output_count];
    }

    /// U^-T followed by the reverse FT row-correction product. `captureEp`
    /// runs this directly in retained work storage so the next update can
    /// consume the partial BTRAN without another dense copy.
    pub fn solveUpperTranspose(self: *SparseForrestTomlin, values: []f64) FtError!void {
        if (values.len != self.dimension) return error.DimensionMismatch;
        try self.solveMutableUpperTransposeFrom(values, 0);
        self.applyCorrectionsTranspose(values);
    }

    /// Solve only the active suffix of mutable U^-T. FT update construction
    /// must not apply historical row corrections here: those corrections are
    /// already represented to the left of the current mutable U.
    fn solveMutableUpperTransposeFrom(self: *SparseForrestTomlin, values: []f64, first_logical: usize) FtError!void {
        for (first_logical..self.logical_count) |logical| {
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
    }

    /// Apply stored row corrections in reverse transpose-solve order.
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

    /// Capture the pivotal BTRAN row needed to form the next FT correction.
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
            var entry = self.column_head[logical];
            while (entry != none) : (entry = self.column_next[entry])
                value -= self.work[self.entry_row[entry]] * self.entry_value[entry];
            value /= self.pivot_values[logical];
            if (!std.math.isFinite(value)) return error.NumericalFailure;
            self.work[pivot] = value;
            self.ep_touched[self.ep_touched_count] = pivot;
            self.ep_touched_count += 1;
        }
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

    /// Allocate and link one coefficient into both U column and row lists.
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

    /// Unlink one coefficient from both orientations and recycle its slot.
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

    /// Grow arrays addressed by stable basis pivot ID.
    fn ensureDimension(self: *SparseForrestTomlin, required: usize) FtError!void {
        if (required <= self.work.len) return;
        const capacity = grow(self.work.len, required) catch return error.CapacityOverflow;
        try self.resizeRetained(&self.pivot_lookup, capacity);
        try self.resizeRetained(&self.row_head, capacity);
        try self.resizeRetained(&self.work, capacity);
        try self.resizeRetained(&self.ep_touched, capacity);
        try self.resizeRetained(&self.captured_aq, capacity);
        try self.resizeRetained(&self.ftran_output, capacity);
    }
    /// Grow arrays addressed by mutable logical U column.
    fn ensureLogical(self: *SparseForrestTomlin, required: usize) FtError!void {
        if (required <= self.logical_capacity) return;
        const capacity = grow(self.logical_capacity, required) catch return error.CapacityOverflow;
        try self.resizeRetained(&self.pivot_ids, capacity);
        try self.resizeRetained(&self.pivot_values, capacity);
        try self.resizeRetained(&self.column_head, capacity);
        self.logical_capacity = capacity;
    }
    /// Grow every parallel intrusive-entry stream.
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

    /// Reallocate one retained slice through the workspace allocator.
    fn resizeRetained(self: *SparseForrestTomlin, slice: anytype, capacity: usize) FtError!void {
        slice.* = self.allocator.realloc(slice.*, capacity) catch return error.OutOfMemory;
    }
    /// Grow parallel correction-index and coefficient streams.
    fn ensureCorrections(self: *SparseForrestTomlin, required: usize) FtError!void {
        if (required <= self.correction_capacity) return;
        const capacity = grow(self.correction_capacity, required) catch return error.CapacityOverflow;
        self.correction_indices = self.allocator.realloc(self.correction_indices, capacity) catch return error.OutOfMemory;
        self.correction_values = self.allocator.realloc(self.correction_values, capacity) catch return error.OutOfMemory;
        self.correction_capacity = capacity;
    }
    /// Grow correction prefix offsets and pivotal-row metadata.
    fn ensureCorrectionStarts(self: *SparseForrestTomlin, required: usize) FtError!void {
        if (required <= self.correction_starts.len) return;
        const capacity = grow(self.correction_starts.len, required) catch return error.CapacityOverflow;
        self.correction_starts = self.allocator.realloc(self.correction_starts, capacity) catch return error.OutOfMemory;
        self.correction_pivots = self.allocator.realloc(self.correction_pivots, capacity) catch return error.OutOfMemory;
    }
};

/// Compute a geometric retained-capacity target.
fn grow(current: usize, required: usize) error{Overflow}!usize {
    var capacity = @max(current, 8);
    while (capacity < required) capacity = std.math.add(usize, capacity, capacity / 2 + 8) catch return error.Overflow;
    if (capacity > std.math.maxInt(u32)) return error.Overflow;
    return capacity;
}

test "captureEp clears retained workspace by its producer-owned support" {
    var ft = SparseForrestTomlin.init(std.testing.allocator);
    defer ft.deinit();
    const pivots = [_]u32{ 0, 1, 2 };
    const diagonal = [_]f64{ 2, 3, 4 };
    const starts = [_]usize{ 0, 0, 0, 0 };
    try ft.reset(&pivots, &diagonal, &starts, &.{}, &.{});

    // A retained allocation may contain values from an older factorization.
    @memset(ft.work[0..3], 99.0);
    try ft.reset(&pivots, &diagonal, &starts, &.{}, &.{});
    try ft.captureEp(1);
    try std.testing.expectEqualSlices(f64, &[_]f64{ 0, 1.0 / 3.0, 0 }, ft.work[0..3]);

    // The next capture clears exactly the IDs written by the previous suffix.
    try ft.captureEp(0);
    try std.testing.expectEqualSlices(f64, &[_]f64{ 0.5, 0, 0 }, ft.work[0..3]);
}

test "forward solve publishes mutable U traversal order after an update" {
    var ft = SparseForrestTomlin.init(std.testing.allocator);
    defer ft.deinit();
    const pivots = [_]u32{ 0, 1, 2 };
    const diagonal = [_]f64{ 2, 3, 4 };
    const starts = [_]usize{ 0, 0, 0, 0 };
    try ft.reset(&pivots, &diagonal, &starts, &.{}, &.{});

    var entering = [_]f64{ 1, 1, 1 };
    try ft.prepareFtran(&entering, true);
    try ft.solveUpper(&entering);
    try ft.captureEp(1);
    try ft.update(1, entering[1], highs_tiny);

    var rhs = [_]f64{ 1, 1, 1 };
    try ft.prepareFtran(&rhs, false);
    try ft.solveUpper(&rhs);
    // Row 1 was removed from its original logical position and appended as
    // the newest FT pivot, so reverse traversal must publish it first.
    try std.testing.expectEqualSlices(u32, &.{ 1, 2, 0 }, ft.lastFtranOutputIndices());

    rhs = .{ 1, 1e-15, 1 };
    try ft.prepareFtran(&rhs, false);
    try ft.solveUpper(&rhs);
    try std.testing.expectEqual(@as(f64, 0), rhs[1]);
    try std.testing.expectEqualSlices(u32, &.{ 2, 0 }, ft.lastFtranOutputIndices());
}
