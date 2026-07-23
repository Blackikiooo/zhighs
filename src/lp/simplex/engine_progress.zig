//! Solve control, statistics/timing, and density observers for
//! `SimplexEngine`.
//!
//! ## Responsibility
//!
//! Owns work/time/callback stop conditions, per-iteration progress events,
//! deterministic clocks, `SimplexStats` accounting, and adaptive pricing
//! density observation.

const std = @import("std");
const basis_module = @import("basis.zig");
const basis_snapshot_module = @import("basis_snapshot.zig");
const factorization_module = @import("factorization.zig");
const pricing_module = @import("pricing.zig");
const ratio_module = @import("ratio_test.zig");
const numerical_module = @import("numerical.zig");
const dual_phase_one_module = @import("dual_phase_one.zig");
const crash_module = @import("crash.zig");
const degeneracy_module = @import("degeneracy.zig");
const pricing_workspace_module = @import("pricing_workspace.zig");
const problem_module = @import("problem.zig");
const solution_module = @import("solution.zig");
const foundation = @import("foundation");
const matrix = @import("matrix");
const SimplexEngine = @import("engine.zig").SimplexEngine;
const SolvePhase = @import("engine.zig").SolvePhase;
const CallbackAction = @import("engine.zig").CallbackAction;
const ProgressEventView = @import("engine.zig").ProgressEventView;
const SolveStatus = @import("engine.zig").SolveStatus;
const SolveControl = @import("engine.zig").SolveControl;
const SimplexStats = @import("engine.zig").SimplexStats;

/// Check interruption, deterministic work and wall-clock limits.
///
/// Returns null when the solve may continue; this function does not charge work.
pub fn controlledStop(self: *SimplexEngine, control: SolveControl) ?SolveStatus {
    if (control.interrupt_flag) |flag| {
        if (flag.load(.acquire)) return .interrupted;
    }
    if (control.work_limit) |limit| {
        if (self.work_used >= limit) return .work_limit;
    }
    if (control.time_limit_ns) |limit| {
        if (limit == 0) return .time_limit;
        if (self.solve_start_ns) |start| {
            const io = self.solve_clock_io orelse return null;
            const now = std.Io.Clock.awake.now(io).nanoseconds;
            if (now >= start and @as(u128, @intCast(now - start)) >= limit) return .time_limit;
        }
    }
    return null;
}

/// Publish due callbacks/log events, then charge one attempted work unit.
pub fn beginIteration(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl, phase: SolvePhase) ?SolveStatus {
    self.current_phase = phase;
    if (self.controlledStop(control)) |status| return status;
    const callback_due = control.iteration_callback != null and
        self.work_used % @max(control.callback_interval_work, 1) == 0;
    const log_due = control.log_level == .iterations and control.log_callback != null and
        self.work_used % @max(control.log_interval_work, 1) == 0;
    if (callback_due or log_due) {
        const event = self.progressEvent(problem, phase);
        if (callback_due) {
            if (control.iteration_callback.?(event, control.callback_user_data) == .stop) return .interrupted;
        }
        if (log_due) control.log_callback.?(event, control.log_user_data);
    }
    self.iteration_counters.attempted_iterations = std.math.add(
        usize,
        self.iteration_counters.attempted_iterations,
        1,
    ) catch std.math.maxInt(usize);
    self.work_used = std.math.add(u64, self.work_used, 1) catch std.math.maxInt(u64);
    return null;
}

/// Build an allocation-free scalar progress snapshot from current engine state.
pub fn progressEvent(self: *const SimplexEngine, problem: problem_module.ProblemView, phase: SolvePhase) ProgressEventView {
    var primal_infeasibility: f64 = 0.0;
    var dual_infeasibility: f64 = 0.0;
    var current_objective = problem.objective_offset;
    if (self.basis) |*basis| {
        for (problem.col_cost, basis.primal[0..problem.num_cols], basis.column_scale[0..problem.num_cols]) |cost, value, scale|
            current_objective += cost * value * scale;
        for (basis.basic_value, basis.basic_lower, basis.basic_upper) |value, lower, upper| {
            primal_infeasibility = @max(primal_infeasibility, @max(lower - value, value - upper));
        }
        for (basis.reduced_cost, basis.col_status) |reduced, status| {
            const violation: f64 = switch (status) {
                .at_lower => -reduced,
                .at_upper => reduced,
                .free, .superbasic => @abs(reduced),
                .basic, .fixed => 0.0,
            };
            dual_infeasibility = @max(dual_infeasibility, violation);
        }
    }
    return .{
        .phase = phase,
        .algorithm = self.algorithm,
        .iterations = self.iterations,
        .work_used = self.work_used,
        .objective_value = current_objective,
        .primal_infeasibility = @max(primal_infeasibility, 0.0),
        .dual_infeasibility = @max(dual_infeasibility, 0.0),
    };
}

/// Bind caller trace buffers and initialize the optional wall-clock deadline.
pub fn startSolveClock(self: *SimplexEngine, control: SolveControl) void {
    self.active_pivot_trace = control.pivot_trace;
    self.pivot_trace_count = 0;
    self.active_degeneracy_trace = control.degeneracy_trace;
    self.degeneracy_trace_count = 0;
    self.degeneracy_basis_fingerprints = @splat(0);
    self.degeneracy_basis_fingerprint_count = 0;
    self.degeneracy_basis_fingerprint_cursor = 0;
    self.dual_phase_one_failure = null;
    self.active_dual_phase_one_candidate_trace = control.dual_phase_one_candidate_trace;
    self.dual_phase_one_candidate_trace_count = 0;
    self.active_dual_phase_one_ep_trace = control.dual_phase_one_ep_trace;
    self.dual_phase_one_ep_trace_count = 0;
    if (control.time_limit_ns == null) {
        self.solve_start_ns = null;
        self.solve_clock_io = null;
        return;
    }
    const io = control.clock_io orelse std.Io.Threaded.global_single_threaded.io();
    self.solve_clock_io = io;
    self.solve_start_ns = std.Io.Clock.awake.now(io).nanoseconds;
}

/// Reset per-solve statistics while preserving counters across recursive restarts.
pub fn resetStatistics(self: *SimplexEngine, control: SolveControl) void {
    const saved_cold_restarts = .{
        self.stats.cold_restart_solves,
        self.stats.cold_restart_phase_one,
    };
    self.stats = .{};
    // Preserve cold-restart counters across recursive solveProblem calls
    // so that restartSolveWithoutPerturbation and restartPhaseOneWithout-
    // Perturbation counts survive the nested resetStatistics.
    self.stats.cold_restart_solves = saved_cold_restarts[0];
    self.stats.cold_restart_phase_one = saved_cold_restarts[1];
    self.statistics_io = if (control.collect_statistics)
        control.clock_io orelse std.Io.Threaded.global_single_threaded.io()
    else
        null;
    self.factorization.resetStatistics(self.statistics_io);
}

/// Read the optional statistics clock; null means timing collection is disabled.
pub fn statisticsTimestamp(self: *const SimplexEngine) ?i96 {
    const io = self.statistics_io orelse return null;
    return std.Io.Clock.awake.now(io).nanoseconds;
}

/// Return saturated nonnegative nanoseconds since `started`, or zero when disabled.
pub fn elapsedSince(self: *const SimplexEngine, started: ?i96) u64 {
    const begin = started orelse return 0;
    const io = self.statistics_io orelse return 0;
    const end = std.Io.Clock.awake.now(io).nanoseconds;
    if (end <= begin) return 0;
    return @intCast(end - begin);
}

/// Accumulate elapsed time in the counter belonging to `phase`.
pub fn recordPhaseElapsed(self: *SimplexEngine, phase: SolvePhase, started: ?i96) void {
    const elapsed = self.elapsedSince(started);
    const target = switch (phase) {
        .phase_one => &self.stats.phase_one_ns,
        .dual_phase_one => &self.stats.dual_phase_one_ns,
        .dual_feasibility_repair => &self.stats.dual_repair_ns,
        .phase_two => &self.stats.phase_two_ns,
    };
    target.* = std.math.add(u64, target.*, elapsed) catch std.math.maxInt(u64);
}

/// Accumulate attempted iterations performed since phase entry.
pub fn recordPhaseIterations(self: *SimplexEngine, phase: SolvePhase, started: usize) void {
    const count = self.iterations -| started;
    const target = switch (phase) {
        .phase_one => &self.stats.phase_one_iterations,
        .dual_phase_one => &self.stats.dual_phase_one_iterations,
        .dual_feasibility_repair => &self.stats.dual_repair_iterations,
        .phase_two => &self.stats.phase_two_iterations,
    };
    target.* = std.math.add(usize, target.*, count) catch std.math.maxInt(usize);
}

/// Record one column-oriented pricing dispatch and its elapsed time.
pub fn recordPricingElapsed(self: *SimplexEngine, started: ?i96) void {
    self.stats.pricing_calls += 1;
    self.stats.column_pricing_dispatches += 1;
    self.stats.pricing_ns = std.math.add(u64, self.stats.pricing_ns, self.elapsedSince(started)) catch std.math.maxInt(u64);
}

/// Record one row-oriented pricing dispatch and its elapsed time.
pub fn recordRowPricingElapsed(self: *SimplexEngine, started: ?i96) void {
    self.stats.pricing_calls += 1;
    self.stats.row_pricing_dispatches += 1;
    self.stats.pricing_ns = std.math.add(u64, self.stats.pricing_ns, self.elapsedSince(started)) catch std.math.maxInt(u64);
}

/// Choose the matrix orientation once for the complete pricing operation.
/// Automatic mode estimates actual CSR work from the active row support;
/// coefficient loops remain branch-free with respect to representation.
pub fn selectRowPricing(self: *const SimplexEngine, vector: []const f64) bool {
    switch (self.active_pricing_kernel) {
        .column => return false,
        .row => return self.pricing_row_view.num_rows == vector.len,
        .automatic => {},
    }
    if (self.pricing_row_view.num_rows != vector.len or self.pricing_row_view.nonzeros == 0) return false;
    var touched_entries: usize = 0;
    for (vector, 0..) |value, row| {
        if (@abs(value) <= self.numerical.zero_tolerance) continue;
        touched_entries = std.math.add(usize, touched_entries, self.pricing_row_view.rowDegree(row)) catch return false;
    }
    // Row traversal wins only with a useful work margin, absorbing scatter
    // stores and the support scan that CSC dot products do not pay.
    return touched_entries * 4 <= self.pricing_row_view.nonzeros * 3;
}

/// Sample the density of an entering FTRAN result when statistics are enabled.
pub fn observeAqDensity(self: *SimplexEngine, vector: []const f64) void {
    if (self.statistics_io == null) return;
    self.stats.aq_samples += 1;
    for (vector) |value| if (@abs(value) > self.numerical.zero_tolerance) {
        self.stats.aq_nonzeros += 1;
    };
}

/// Sample the density of a pivotal BTRAN row when statistics are enabled.
pub fn observeEpDensity(self: *SimplexEngine, vector: []const f64) void {
    if (self.statistics_io == null) return;
    self.stats.ep_samples += 1;
    for (vector) |value| if (@abs(value) > self.numerical.zero_tolerance) {
        self.stats.ep_nonzeros += 1;
    };
}

/// Sample the density of a pricing vector when statistics are enabled.
pub fn observePricingDensity(self: *SimplexEngine, vector: []const f64) void {
    if (self.statistics_io == null) return;
    self.stats.pricing_samples += 1;
    self.stats.pricing_entries += vector.len;
    for (vector) |value| if (@abs(value) > self.numerical.dual_tolerance) {
        self.stats.pricing_nonzeros += 1;
    };
}

test "deterministic work limit spans Phase I and Phase II" {
    const rows = [_]foundation.RowId{foundation.RowId.fromUsizeAssumeValid(0)};
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{1},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{std.math.inf(f64)},
        .row_lower = &[_]f64{2},
        .row_upper = &[_]f64{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &[_]usize{ 0, 1 }, &rows, &[_]f64{1}),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.work_limit, engine.solveProblem(problem, .{ .work_limit = 1 }));
    try std.testing.expectEqual(@as(u64, 1), engine.work_used);
}

test "iteration callback can stop without allocating callback state" {
    const Context = struct {
        calls: usize = 0,
        last_work: u64 = 0,

        /// Test callback recording its invocation before requesting interruption.
        fn callback(event: ProgressEventView, context_ptr: ?*anyopaque) CallbackAction {
            const self: *@This() = @ptrCast(@alignCast(context_ptr.?));
            self.calls += 1;
            self.last_work = event.work_used;
            return .stop;
        }
    };
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{-1},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{1},
        .row_lower = &[_]f64{-std.math.inf(f64)},
        .row_upper = &[_]f64{1},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &[_]usize{ 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var context = Context{};
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.interrupted, engine.solveProblem(problem, .{
        .iteration_callback = Context.callback,
        .callback_user_data = &context,
    }));
    try std.testing.expectEqual(@as(usize, 1), context.calls);
    try std.testing.expectEqual(@as(u64, 0), context.last_work);
    try std.testing.expectEqual(@as(u64, 0), engine.work_used);
}

test "structured iteration logging obeys its deterministic interval" {
    const Context = struct {
        calls: usize = 0,
        /// Test logger counting deterministic interval publications.
        fn log(_: ProgressEventView, context_ptr: ?*anyopaque) void {
            const self: *@This() = @ptrCast(@alignCast(context_ptr.?));
            self.calls += 1;
        }
    };
    const problem = problem_module.ProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &[_]f64{-1},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{1},
        .row_lower = &[_]f64{-std.math.inf(f64)},
        .row_upper = &[_]f64{1},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &[_]usize{ 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var context = Context{};
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, engine.solveProblem(problem, .{
        .log_level = .iterations,
        .log_callback = Context.log,
        .log_user_data = &context,
        .log_interval_work = 1,
    }));
    try std.testing.expect(context.calls >= 1);
}

test "engine honors an immediate time limit" {
    const problem = problem_module.ProblemView{
        .num_rows = 0,
        .num_cols = 1,
        .col_cost = &[_]f64{0},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{1},
        .row_lower = &.{},
        .row_upper = &.{},
        .matrix = matrix.CscView.initAssumeValid(0, 1, &[_]usize{ 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.time_limit, engine.solveProblem(problem, .{ .time_limit_ns = 0 }));
}

test "engine honors a caller-owned atomic interrupt flag" {
    const problem = problem_module.ProblemView{
        .num_rows = 0,
        .num_cols = 1,
        .col_cost = &[_]f64{0},
        .col_lower = &[_]f64{0},
        .col_upper = &[_]f64{1},
        .row_lower = &.{},
        .row_upper = &.{},
        .matrix = matrix.CscView.initAssumeValid(0, 1, &[_]usize{ 0, 0 }, &.{}, &.{}),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var interrupted = std.atomic.Value(bool).init(true);
    var engine = SimplexEngine.init(std.testing.allocator);
    defer engine.deinit();
    try std.testing.expectEqual(SolveStatus.interrupted, engine.solveProblem(problem, .{ .interrupt_flag = &interrupted }));
}
