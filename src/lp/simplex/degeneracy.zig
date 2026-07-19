//! Reusable state for deterministic degeneracy policies.
//!
//! Perturbations are virtual lexicographic ranks: they alter only tie order,
//! while ratio tests still return the exact unperturbed step. This preserves
//! basis equations and keeps cleanup in the original model coordinates.

const std = @import("std");

pub const TabooChange = struct {
    entering: u32 = 0,
    leaving: u32 = 0,
    direction: i8 = 0,
    generation: u32 = 0,
};

pub const Workspace = struct {
    allocator: std.mem.Allocator,
    row_rank: []f64 = &.{},
    column_rank: []f64 = &.{},
    row_generation: []u32 = &.{},
    column_generation: []u32 = &.{},
    taboo_until: []usize = &.{},
    reduced_cost_snapshot: []f64 = &.{},
    taboo: [64]TabooChange = @splat(.{}),
    generation: u32 = 0,
    basis_epoch: u64 = 0,
    active_age: usize = 0,
    level: u8 = 0,
    active: bool = false,
    ever_active: bool = false,
    taboo_cursor: usize = 0,

    pub fn init(allocator: std.mem.Allocator) Workspace {
        return .{ .allocator = allocator };
    }

    pub fn deinit(self: *Workspace) void {
        self.allocator.free(self.row_rank);
        self.allocator.free(self.column_rank);
        self.allocator.free(self.row_generation);
        self.allocator.free(self.column_generation);
        self.allocator.free(self.taboo_until);
        self.allocator.free(self.reduced_cost_snapshot);
        self.* = .{ .allocator = self.allocator };
    }

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

    pub fn advance(self: *Workspace, maximum_age: usize) bool {
        if (!self.active) return false;
        self.active_age += 1;
        if (self.active_age < maximum_age) return false;
        self.active = false;
        self.level = @min(self.level + 1, 3);
        return true;
    }

    pub fn escalateAfterRepeat(self: *Workspace) void {
        self.active = false;
        self.active_age = 0;
        self.level = @min(self.level + 1, 3);
    }

    pub fn clearAfterProgress(self: *Workspace) void {
        self.active = false;
        self.active_age = 0;
        self.level = 0;
        self.invalidateTaboo();
    }

    pub fn resetSolve(self: *Workspace) void {
        self.active = false;
        self.ever_active = false;
        self.active_age = 0;
        self.level = 0;
        self.invalidateTaboo();
    }

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

    pub fn invalidateTaboo(self: *Workspace) void {
        self.generation +%= 1;
        if (self.generation == 0) self.generation = 1;
        self.taboo_cursor = 0;
        @memset(self.taboo_until, 0);
    }

    pub fn requestedBytes(self: *const Workspace) usize {
        return std.mem.sliceAsBytes(self.row_rank).len +
            std.mem.sliceAsBytes(self.column_rank).len +
            std.mem.sliceAsBytes(self.row_generation).len +
            std.mem.sliceAsBytes(self.column_generation).len +
            std.mem.sliceAsBytes(self.taboo_until).len +
            std.mem.sliceAsBytes(self.reduced_cost_snapshot).len +
            @sizeOf(@TypeOf(self.taboo));
    }

    fn deterministicUnit(index: usize, epoch: u64) f64 {
        var value = @as(u64, @intCast(index)) +% epoch *% 0x9e3779b97f4a7c15;
        value = (value ^ (value >> 30)) *% 0xbf58476d1ce4e5b9;
        value = (value ^ (value >> 27)) *% 0x94d049bb133111eb;
        value ^= value >> 31;
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
