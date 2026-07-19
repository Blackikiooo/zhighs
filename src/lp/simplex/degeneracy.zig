//! Reusable state for deterministic degeneracy policies.
//!
//! Perturbations are deterministic bounded lexicographic margins. Ratio tests
//! may take the corresponding sub-tolerance primal step to leave a degenerate
//! face; terminal cleanup always rebuilds in original model coordinates.

const std = @import("std");

/// Record of a single basis change, stored in the taboo circular buffer.
pub const TabooChange = struct {
    entering: u32 = 0,
    leaving: u32 = 0,
    direction: i8 = 0,
    generation: u32 = 0,
};

pub const Workspace = struct {
    allocator: std.mem.Allocator,
    // Perturbation magnitudes for rows and columns, deterministically generated per basis_epoch
    row_rank: []f64 = &.{},
    column_rank: []f64 = &.{},
    // Tracks which generation each row/column last participated in perturbation (for cache invalidation)
    row_generation: []u32 = &.{},
    column_generation: []u32 = &.{},
    // Column index -> taboo expiration (iteration count), 0 means not tabooed
    taboo_until: []usize = &.{},
    // Reduced cost snapshot for detecting progress during perturbation
    reduced_cost_snapshot: []f64 = &.{},
    // Circular buffer: records the last 64 basis changes (entering/leaving/direction)
    taboo: [64]TabooChange = @splat(.{}),
    // Current perturbation generation (incrementing invalidates old records)
    generation: u32 = 0,
    // Basis change counter (one of the inputs to deterministicUnit)
    basis_epoch: u64 = 0,
    // Iteration count since activation (used to decide when to escalate)
    active_age: usize = 0,
    // Perturbation intensity level 0-3 (escalation multiplier)
    level: u8 = 0,
    // Whether perturbation is currently active
    active: bool = false,
    // Whether perturbation has ever been activated (for diagnostics)
    ever_active: bool = false,
    // Write cursor for the taboo circular buffer
    taboo_cursor: usize = 0,

    /// Construct an empty workspace; arrays are allocated lazily by `ensureCapacity`.
    pub fn init(allocator: std.mem.Allocator) Workspace {
        return .{ .allocator = allocator };
    }

    /// Free all allocated buffers and reset to the empty state.
    pub fn deinit(self: *Workspace) void {
        self.allocator.free(self.row_rank);
        self.allocator.free(self.column_rank);
        self.allocator.free(self.row_generation);
        self.allocator.free(self.column_generation);
        self.allocator.free(self.taboo_until);
        self.allocator.free(self.reduced_cost_snapshot);
        self.* = .{ .allocator = self.allocator };
    }

    /// Grow each buffer as needed to hold `rows` rows and `columns` columns.
    pub fn ensureCapacity(self: *Workspace, rows: usize, columns: usize) !void {
        if (self.row_rank.len < rows) self.row_rank = try self.allocator.realloc(self.row_rank, rows);
        if (self.column_rank.len < columns) self.column_rank = try self.allocator.realloc(self.column_rank, columns);
        if (self.row_generation.len < rows) self.row_generation = try self.allocator.realloc(self.row_generation, rows);
        if (self.column_generation.len < columns) self.column_generation = try self.allocator.realloc(self.column_generation, columns);
        if (self.taboo_until.len < columns) self.taboo_until = try self.allocator.realloc(self.taboo_until, columns);
        if (self.reduced_cost_snapshot.len < columns)
            self.reduced_cost_snapshot = try self.allocator.realloc(self.reduced_cost_snapshot, columns);
    }

    /// Start a new deterministic perturbation epoch. Magnitudes are bounded by
    /// the solver feasibility tolerances and never accumulate across epochs.
    pub fn activate(
        self: *Workspace,
        row_scale: []const f64,
        column_scale: []const f64,
        primal_tolerance: f64,
        dual_tolerance: f64,
        objective_scale: f64,
    ) !void {
        const rows = row_scale.len;
        const columns = column_scale.len;
        try self.ensureCapacity(rows, columns);
        self.generation +%= 1;
        if (self.generation == 0) {
            @memset(self.row_generation, 0);
            @memset(self.column_generation, 0);
            self.generation = 1;
        }
        self.basis_epoch +%= 1;
        self.active_age = 0;
        self.active = true;
        self.ever_active = true;
        const escalation = std.math.scalbn(@as(f64, 1.0), @as(i32, @intCast(@min(self.level, 3))));
        var minimum_row_scale: f64 = 1.0;
        for (row_scale) |scale| minimum_row_scale = @min(minimum_row_scale, @abs(scale));
        var minimum_column_scale: f64 = 1.0;
        for (column_scale) |scale| minimum_column_scale = @min(minimum_column_scale, @abs(scale));
        const row_bound = primal_tolerance * minimum_row_scale;
        const column_bound = dual_tolerance * objective_scale * minimum_column_scale;
        const row_magnitude = @min(row_bound * 0.25 * escalation, row_bound * 2.0);
        const column_magnitude = @min(column_bound * 0.25 * escalation, column_bound * 2.0);
        // A common epoch magnitude preserves deterministic id ordering. Using
        // a per-entry magnitude would silently reorder ties solely because of
        // equilibration and proved much less stable on high-fill bases.
        for (self.row_rank[0..rows], self.row_generation[0..rows], 0..) |*rank, *generation, row| {
            rank.* = deterministicUnit(row, self.basis_epoch) * row_magnitude;
            generation.* = self.generation;
        }
        for (self.column_rank[0..columns], self.column_generation[0..columns], 0..) |*rank, *generation, column| {
            rank.* = deterministicUnit(column + rows, self.basis_epoch) * column_magnitude;
            generation.* = self.generation;
        }
    }

    /// Advance the perturbation age and escalate to the next level if it exceeds `maximum_age`.
    /// Returns `true` when an escalation has been triggered.
    pub fn advance(self: *Workspace, maximum_age: usize) bool {
        if (!self.active) return false;
        self.active_age += 1;
        if (self.active_age < maximum_age) return false;
        self.active = false;
        self.level = @min(self.level + 1, 3);
        return true;
    }

    /// Escalate to a stronger perturbation level immediately after detecting a repeat pattern.
    pub fn escalateAfterRepeat(self: *Workspace) void {
        self.active = false;
        self.active_age = 0;
        self.level = @min(self.level + 1, 3);
    }

    /// Reset perturbation state after observing real progress; clears taboo state too.
    pub fn clearAfterProgress(self: *Workspace) void {
        self.active = false;
        self.active_age = 0;
        self.level = 0;
        self.invalidateTaboo();
    }

    /// Full reset between solves; clears all perturbation history and taboo state.
    pub fn resetSolve(self: *Workspace) void {
        self.active = false;
        self.ever_active = false;
        self.active_age = 0;
        self.level = 0;
        self.invalidateTaboo();
    }

    /// Record a forbidden basis change: store it in the circular buffer and
    /// mark the entering column taboo until `iteration + lifetime`.
    pub fn recordTaboo(
        self: *Workspace,
        entering: u32,
        leaving: u32,
        direction: f64,
        iteration: usize,
        lifetime: usize,
    ) void {
        self.taboo[self.taboo_cursor] = .{
            .entering = entering,
            .leaving = leaving,
            .direction = if (direction < 0.0) -1 else 1,
            .generation = self.generation,
        };
        self.taboo_cursor = (self.taboo_cursor + 1) % self.taboo.len;
        const column: usize = @intCast(entering);
        if (column < self.taboo_until.len)
            self.taboo_until[column] = std.math.add(usize, iteration, lifetime) catch std.math.maxInt(usize);
    }

    /// Invalidate all taboo entries by bumping the generation counter and clearing expiration times.
    pub fn invalidateTaboo(self: *Workspace) void {
        self.generation +%= 1;
        if (self.generation == 0) self.generation = 1;
        self.taboo_cursor = 0;
        @memset(self.taboo_until, 0);
    }

    /// Total bytes currently held by all dynamic buffers (used for memory budgeting).
    pub fn requestedBytes(self: *const Workspace) usize {
        return std.mem.sliceAsBytes(self.row_rank).len +
            std.mem.sliceAsBytes(self.column_rank).len +
            std.mem.sliceAsBytes(self.row_generation).len +
            std.mem.sliceAsBytes(self.column_generation).len +
            std.mem.sliceAsBytes(self.taboo_until).len +
            std.mem.sliceAsBytes(self.reduced_cost_snapshot).len +
            @sizeOf(@TypeOf(self.taboo));
    }

    /// Generates a deterministic pseudo-random value in [0, 1) from an index and epoch.
    ///
    /// Uses splitmix64-style linear congruential generator with no mutable state.
    /// The mapping from (index, epoch) to output value is pure and injective-ish,
    /// making it suitable for lexicographic perturbation tiebreaking where results
    /// must be reproducible across solver restarts.
    ///
    /// - `index`: unique identifier for the row/column being perturbed
    /// - `epoch`: monotonically increasing basis generation counter
    ///
    /// The 64-bit state is derived by hashing `index` with `epoch * GOLDEN_RATIO`,
    /// then passed through three xorshift-multiply stages for bit mixing.
    /// Final output takes the upper 53 bits of the state, scaled to [0, 1) in
    /// IEEE 754 double precision (exactly representable as 2^-53 steps).
    fn deterministicUnit(index: usize, epoch: u64) f64 {
        // Stage 1: seed the 64-bit state using the golden ratio constant.
        // Adding epoch first ensures different epochs produce disjoint value sequences.
        var value = @as(u64, @intCast(index)) +% epoch *% 0x9e3779b97f4a7c15;

        // Stage 2-4: three rounds of xorshift-multiply for full state diffusion.
        // Each multiply uses a different odd constant to avoid cycles.
        value = (value ^ (value >> 30)) *% 0xbf58476d1ce4e5b9;
        value = (value ^ (value >> 27)) *% 0x94d049bb133111eb;
        value ^= value >> 31;

        // Take upper 53 bits and normalize to [0, 1).
        // Using >> 11 keeps the mantissa bits while discarding the 11 lowest bits
        // (including the sign bit which is always 0 after xorshift).
        // Adding 1.0 ensures the result is in (1, 2) before scaling.
        return (@as(f64, @floatFromInt(value >> 11)) + 1.0) * 0x1.0p-53;
    }
};

test "perturbation ranks are deterministic and bounded" {
    var workspace = Workspace.init(std.testing.allocator);
    defer workspace.deinit();
    try workspace.activate(&.{ 1.0, 0.5, 2.0 }, &.{ 1.0, 0.25, 2.0, 1.0 }, 1e-7, 2e-7, 1.0);
    const first = workspace.row_rank[0];
    try std.testing.expect(first >= 0.0 and first <= 2.5e-8);
    workspace.basis_epoch -%= 1;
    try workspace.activate(&.{ 1.0, 0.5, 2.0 }, &.{ 1.0, 0.25, 2.0, 1.0 }, 1e-7, 2e-7, 1.0);
    try std.testing.expectEqual(first, workspace.row_rank[0]);
}
