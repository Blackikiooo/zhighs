//! MIP-specific ownership boundary around the reusable LP solve session.
//!
//! The root may use the solver's reversible presolve pipeline. Tree nodes use
//! the stable original LP mapping, allowing bound-only basis reoptimization.

const std = @import("std");
const solver = @import("solver");
const matrix = @import("matrix");

pub const Statistics = struct {
    cold_solves: usize,
    reoptimizations: usize,
    strong_branching_probes: usize,
    strong_branching_iterations: usize,
    strong_branching_work: u64,
    cut_rounds: usize,
    cuts_added: usize,
    cuts_expired: usize,
    clique_edges: usize,
    clique_tightenings: usize,
    reduced_cost_tightenings: usize,
    root_reduced_cost_tightenings: usize,
    conflicts_learned: usize,
    conflict_tightenings: usize,
    conflicts_expired: usize,
    conflicts_compacted: usize,
    conflict_checks: usize,
};

pub const Relaxation = struct {
    session: solver.LpSolveSession,
    strong_branching_probes: usize = 0,
    strong_branching_iterations: usize = 0,
    strong_branching_work: u64 = 0,
    cut_rounds: usize = 0,
    cuts_added: usize = 0,
    cuts_expired: usize = 0,
    clique_edges: usize = 0,
    clique_tightenings: usize = 0,
    reduced_cost_tightenings: usize = 0,
    root_reduced_cost_tightenings: usize = 0,
    conflicts_learned: usize = 0,
    conflict_tightenings: usize = 0,
    conflicts_expired: usize = 0,
    conflicts_compacted: usize = 0,
    conflict_checks: usize = 0,

    pub fn init(allocator: std.mem.Allocator) Relaxation {
        return .{ .session = solver.LpSolveSession.init(allocator) };
    }

    pub fn deinit(self: *Relaxation) void {
        self.session.deinit();
        self.* = undefined;
    }

    pub fn solve(
        self: *Relaxation,
        problem: solver.LpProblemView,
        control: solver.LpSolveControl,
        node_number: usize,
        revision: solver.LpRevisionStamp,
    ) solver.LpSolveStatus {
        return if (node_number == 1)
            self.session.solve(problem, control, revision)
        else
            self.session.solveRelaxation(problem, control, revision);
    }

    pub fn resultView(self: *const Relaxation, problem: solver.LpProblemView, status: solver.LpSolveStatus) ?solver.LpSolveResultView {
        return self.session.resultView(problem, status);
    }

    pub fn lastIterations(self: *const Relaxation) usize {
        return self.session.engine.iterations;
    }

    pub fn lastWork(self: *const Relaxation) u64 {
        return self.session.engine.work_used;
    }

    pub fn recordStrongBranching(self: *Relaxation, probes: usize, iterations: usize, work: u64) void {
        self.strong_branching_probes = std.math.add(usize, self.strong_branching_probes, probes) catch std.math.maxInt(usize);
        self.strong_branching_iterations = std.math.add(usize, self.strong_branching_iterations, iterations) catch std.math.maxInt(usize);
        self.strong_branching_work = std.math.add(u64, self.strong_branching_work, work) catch std.math.maxInt(u64);
    }

    pub fn recordCutRound(self: *Relaxation, cuts_added: usize) void {
        self.cut_rounds = std.math.add(usize, self.cut_rounds, 1) catch std.math.maxInt(usize);
        self.cuts_added = std.math.add(usize, self.cuts_added, cuts_added) catch std.math.maxInt(usize);
    }

    pub fn recordExpiredCuts(self: *Relaxation, expired: usize) void {
        self.cuts_expired = std.math.add(usize, self.cuts_expired, expired) catch std.math.maxInt(usize);
    }

    pub fn recordCliqueTable(self: *Relaxation, edges: usize) void {
        self.clique_edges = edges;
    }

    pub fn recordCliquePropagation(self: *Relaxation, tightenings: usize) void {
        self.clique_tightenings = std.math.add(usize, self.clique_tightenings, tightenings) catch std.math.maxInt(usize);
    }

    pub fn recordReducedCostFixing(self: *Relaxation, tightenings: usize) void {
        self.reduced_cost_tightenings = std.math.add(usize, self.reduced_cost_tightenings, tightenings) catch std.math.maxInt(usize);
    }

    pub fn recordRootReducedCostFixing(self: *Relaxation, tightenings: usize) void {
        self.root_reduced_cost_tightenings = std.math.add(usize, self.root_reduced_cost_tightenings, tightenings) catch std.math.maxInt(usize);
    }

    pub fn recordLearnedConflict(self: *Relaxation) void {
        self.conflicts_learned = std.math.add(usize, self.conflicts_learned, 1) catch std.math.maxInt(usize);
    }

    pub fn recordConflictPropagation(self: *Relaxation, checks: usize, tightenings: usize, expired: usize, compacted: usize) void {
        self.conflict_checks = std.math.add(usize, self.conflict_checks, checks) catch std.math.maxInt(usize);
        self.conflict_tightenings = std.math.add(usize, self.conflict_tightenings, tightenings) catch std.math.maxInt(usize);
        self.conflicts_expired = std.math.add(usize, self.conflicts_expired, expired) catch std.math.maxInt(usize);
        self.conflicts_compacted = std.math.add(usize, self.conflicts_compacted, compacted) catch std.math.maxInt(usize);
    }

    pub fn statistics(self: *const Relaxation) Statistics {
        return .{
            .cold_solves = self.session.cold_solves,
            .reoptimizations = self.session.reoptimizations,
            .strong_branching_probes = self.strong_branching_probes,
            .strong_branching_iterations = self.strong_branching_iterations,
            .strong_branching_work = self.strong_branching_work,
            .cut_rounds = self.cut_rounds,
            .cuts_added = self.cuts_added,
            .cuts_expired = self.cuts_expired,
            .clique_edges = self.clique_edges,
            .clique_tightenings = self.clique_tightenings,
            .reduced_cost_tightenings = self.reduced_cost_tightenings,
            .root_reduced_cost_tightenings = self.root_reduced_cost_tightenings,
            .conflicts_learned = self.conflicts_learned,
            .conflict_tightenings = self.conflict_tightenings,
            .conflicts_expired = self.conflicts_expired,
            .conflicts_compacted = self.conflicts_compacted,
            .conflict_checks = self.conflict_checks,
        };
    }
};

test "tree relaxations leave root presolve once then reoptimize stable coordinates" {
    var lower = [_]f64{0.0};
    var upper = [_]f64{0.0};
    const problem = solver.LpProblemView{
        .num_rows = 0,
        .num_cols = 1,
        .col_cost = &.{1.0},
        .col_lower = &lower,
        .col_upper = &upper,
        .row_lower = &.{},
        .row_upper = &.{},
        .matrix = matrix.CscView.initAssumeValid(0, 1, &.{ 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var relaxation = Relaxation.init(std.testing.allocator);
    defer relaxation.deinit();

    try std.testing.expectEqual(solver.LpSolveStatus.optimal, relaxation.solve(problem, .{}, 1, .{ .structure = 0, .matrix_values = 0, .bounds = 1, .objective = 0 }));
    try std.testing.expect(relaxation.session.presolved_problem != null);

    upper[0] = 1.0;
    try std.testing.expectEqual(solver.LpSolveStatus.optimal, relaxation.solve(problem, .{}, 2, .{ .structure = 0, .matrix_values = 0, .bounds = 2, .objective = 0 }));
    try std.testing.expect(relaxation.session.presolved_problem == null);
    upper[0] = 2.0;
    try std.testing.expectEqual(solver.LpSolveStatus.optimal, relaxation.solve(problem, .{}, 3, .{ .structure = 0, .matrix_values = 0, .bounds = 3, .objective = 0 }));

    const stats = relaxation.statistics();
    try std.testing.expectEqual(@as(usize, 2), stats.cold_solves);
    try std.testing.expectEqual(@as(usize, 1), stats.reoptimizations);
}
