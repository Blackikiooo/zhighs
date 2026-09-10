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
};

pub const Relaxation = struct {
    session: solver.LpSolveSession,

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

    pub fn statistics(self: *const Relaxation) Statistics {
        return .{
            .cold_solves = self.session.cold_solves,
            .reoptimizations = self.session.reoptimizations,
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
