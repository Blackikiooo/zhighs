//! Solve-local coordination of cut storage, separation, aging, and LP sync.

const std = @import("std");
const solver = @import("solver");
const plugin = @import("plugin");
const CutPool = @import("cut_pool.zig").CutPool;
const CutRelaxation = @import("cut_relaxation.zig").CutRelaxation;
const CoverSeparator = @import("separation.zig").CoverSeparator;

pub const RoundResult = struct {
    rows_considered: usize = 0,
    cuts_added: usize = 0,
    duplicates: usize = 0,
    structure_changed: bool = false,
};

pub const CutGeneration = struct {
    allocator: std.mem.Allocator,
    pool: CutPool,
    relaxation: CutRelaxation,
    separator: CoverSeparator,

    pub fn init(allocator: std.mem.Allocator, problem: solver.LpProblemView) std.mem.Allocator.Error!CutGeneration {
        var pool = CutPool.init(allocator, problem.num_cols) catch |err| switch (err) {
            error.OutOfMemory => return error.OutOfMemory,
            else => unreachable,
        };
        errdefer pool.deinit(allocator);
        var relaxation = CutRelaxation.init(allocator, problem) catch |err| switch (err) {
            error.OutOfMemory => return error.OutOfMemory,
            else => unreachable,
        };
        errdefer relaxation.deinit();
        const separator = try CoverSeparator.init(allocator, problem.num_cols);
        return .{ .allocator = allocator, .pool = pool, .relaxation = relaxation, .separator = separator };
    }

    pub fn deinit(self: *CutGeneration) void {
        self.separator.deinit();
        self.relaxation.deinit();
        self.pool.deinit(self.allocator);
        self.* = undefined;
    }

    pub fn view(self: *const CutGeneration, lower: []const f64, upper: []const f64) solver.LpProblemView {
        return self.relaxation.view(lower, upper);
    }

    pub fn structureRevision(self: *const CutGeneration) usize {
        return self.relaxation.structure_revision;
    }

    pub fn synchronize(self: *CutGeneration) std.mem.Allocator.Error!bool {
        return self.relaxation.synchronize(&self.pool) catch |err| switch (err) {
            error.OutOfMemory => return error.OutOfMemory,
            else => unreachable,
        };
    }

    pub fn separateRoot(
        self: *CutGeneration,
        integrality: []const plugin.Integrality,
        primal: []const f64,
        tolerance: f64,
        max_cuts: usize,
    ) std.mem.Allocator.Error!RoundResult {
        const separated = self.separator.separate(
            self.relaxation.base,
            self.relaxation.originalRows(),
            integrality,
            primal,
            &self.pool,
            tolerance,
            max_cuts,
        ) catch |err| switch (err) {
            error.OutOfMemory => return error.OutOfMemory,
            else => unreachable,
        };
        const changed = if (separated.cuts_added != 0) try self.synchronize() else false;
        return .{
            .rows_considered = separated.rows_considered,
            .cuts_added = separated.cuts_added,
            .duplicates = separated.duplicates,
            .structure_changed = changed,
        };
    }

    pub fn age(self: *CutGeneration, primal: []const f64, tolerance: f64, max_age: u16) usize {
        return self.pool.ageAt(primal, tolerance, max_age);
    }
};

test "cut generation synchronizes a separated cover into the LP view" {
    const matrix = @import("matrix");
    const foundation = @import("foundation");
    const problem = solver.LpProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 3.0, 2.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ 1.0, 1.0 },
        .row_lower = &.{-std.math.inf(f64)},
        .row_upper = &.{4.0},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &.{ 0, 1, 2 }, &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) }, &.{ 2.0, 3.0 }),
        .objective_sense = .maximize,
        .objective_offset = 0.0,
    };
    var generation = try CutGeneration.init(std.testing.allocator, problem);
    defer generation.deinit();
    const round = try generation.separateRoot(&.{ .binary, .binary }, &.{ 1.0, 2.0 / 3.0 }, 1e-6, 8);
    try std.testing.expectEqual(@as(usize, 1), round.cuts_added);
    try std.testing.expect(round.structure_changed);
    try std.testing.expectEqual(@as(usize, 2), generation.view(problem.col_lower, problem.col_upper).num_rows);
}
