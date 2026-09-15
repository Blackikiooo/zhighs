//! Plugin-driven serial branch-and-bound search over the LP solver boundary.
//!
//! The engine owns node management, LP reoptimization, incumbent bookkeeping,
//! and limit enforcement. Every algorithmic decision is delegated to
//! components registered through a `framework.PluginRegistry`:
//!
//! * branching rules pick the disjunction that splits a fractional node,
//! * primal heuristics propose incumbent candidates (validated here),
//! * a node selector picks the next open node to process.
//!
//! The default component set reproduces the historical behaviour exactly:
//! most-fractional branching, LP rounding, and depth-first search.

const std = @import("std");
const solver = @import("solver");
const plugin = @import("plugin");
const framework = @import("framework");
const builtin_plugins = @import("plugins_builtin");
const matrix = @import("matrix");
const foundation = @import("foundation");
const types = @import("types.zig");
const domain = @import("domain.zig");
const incumbent_candidate = @import("incumbent.zig");
const objective_bounds = @import("bounds.zig");
const propagation = @import("propagation.zig");
const CliqueTable = @import("clique_table.zig").CliqueTable;
const reduced_cost_fixing = @import("reduced_cost_fixing.zig");
const RootReducedCosts = @import("root_reduced_costs.zig").RootReducedCosts;
const ConflictPool = @import("conflict_pool.zig").ConflictPool;
const Pseudocost = @import("pseudocost.zig").Pseudocost;
const StrongBranching = @import("strong_branching.zig").StrongBranching;
const CutGeneration = @import("cut_generation.zig").CutGeneration;
const NodeQueue = @import("node_queue.zig").NodeQueue;
const Relaxation = @import("relaxation.zig").Relaxation;

pub const Integrality = types.Integrality;
pub const ProblemView = types.ProblemView;
pub const SolveStatus = types.SolveStatus;
pub const SolveComponents = types.SolveComponents;
pub const SolveControl = types.SolveControl;
pub const SolveResult = types.SolveResult;

pub fn solve(allocator: std.mem.Allocator, problem: ProblemView, control: SolveControl) std.mem.Allocator.Error!SolveResult {
    problem.validate() catch return .{ .allocator = allocator, .status = .invalid_model };
    if (!types.validControl(control)) return .{ .allocator = allocator, .status = .invalid_model };

    var registry = framework.PluginRegistry.init(allocator);
    defer registry.deinit();

    var default_branching = builtin_plugins.branching.Pseudocost{};
    var default_heuristic = builtin_plugins.heuristics.Rounding{};
    var default_diving = builtin_plugins.heuristics.Diving{};
    // The diving heuristic retains an LP session across invocations; release
    // it before the search-local model slices go away.
    defer default_diving.deinit();
    if (control.components.branching_rule) |rule| {
        try registry.addBranchingRule(rule, 0);
    } else {
        try registry.addBranchingRule(plugin.bindBranching(builtin_plugins.branching.Pseudocost, &default_branching), 0);
    }
    if (control.components.node_selector) |selector| try registry.addNodeSelector(selector, 0);
    if (control.components.heuristics.len != 0) {
        for (control.components.heuristics) |heuristic| try registry.addHeuristic(heuristic, 0);
    } else if (control.rounding_heuristic) {
        try registry.addHeuristic(plugin.bindHeuristic(builtin_plugins.heuristics.Rounding, &default_heuristic), 0);
        // Diving re-solves LPs and is therefore the more expensive fallback;
        // it runs after the cheap rounding pass (lower priority).
        try registry.addHeuristic(plugin.bindHeuristic(builtin_plugins.heuristics.Diving, &default_diving), -1);
    }
    const node_selector = registry.primaryNodeSelector();

    var current_domain = (try domain.Domain.init(
        allocator,
        problem.integrality,
        problem.relaxation.col_lower,
        problem.relaxation.col_upper,
        control.integer_tolerance,
    )) orelse return .{ .allocator = allocator, .status = .infeasible };
    defer current_domain.deinit(allocator);
    var change_arena = domain.DomainChangeArena{};
    defer change_arena.deinit(allocator);

    var open_nodes = NodeQueue{};
    defer open_nodes.deinit(allocator);
    try open_nodes.append(allocator, .{});

    var incumbent: []f64 = &.{};
    errdefer allocator.free(incumbent);
    var incumbent_objective: f64 = 0.0;
    var have_incumbent = false;
    var nodes: usize = 0;
    var iterations: usize = 0;
    var work: u64 = 0;
    var solutions: usize = 0;
    var lp_relaxation = Relaxation.init(allocator);
    defer lp_relaxation.deinit();
    var pseudocost = try Pseudocost.init(allocator, problem.relaxation.num_cols, control.pseudocost_min_reliable);
    defer pseudocost.deinit(allocator);
    var root_reduced_costs = try RootReducedCosts.init(allocator, problem.relaxation.num_cols);
    defer root_reduced_costs.deinit();
    var conflict_pool = try ConflictPool.init(allocator, problem.relaxation.num_cols, control.max_conflict_entries, control.max_conflict_length);
    defer conflict_pool.deinit();
    var strong_branching: ?StrongBranching = if (control.reliability_branching and control.components.branching_rule == null)
        try StrongBranching.init(allocator, problem.relaxation.num_cols)
    else
        null;
    defer if (strong_branching) |*engine| engine.deinit();
    var cut_generation: ?CutGeneration = if (control.cut_separation)
        try CutGeneration.init(allocator, problem.relaxation)
    else
        null;
    defer if (cut_generation) |*generation| generation.deinit();
    var propagation_workspace: ?propagation.Workspace = if (control.propagation)
        try propagation.Workspace.init(allocator, problem.relaxation.num_rows, problem.relaxation.num_cols)
    else
        null;
    defer if (propagation_workspace) |*workspace| workspace.deinit(allocator);
    var clique_table: ?CliqueTable = if (control.clique_propagation)
        try CliqueTable.init(
            allocator,
            problem.relaxation,
            problem.integrality,
            current_domain.root_lower,
            current_domain.root_upper,
            control.integer_tolerance,
            control.max_clique_row_entries,
            control.max_clique_edges,
        )
    else
        null;
    defer if (clique_table) |*table| table.deinit();
    if (clique_table) |*table| lp_relaxation.recordCliqueTable(table.edgeCount());
    const clock_io = if (control.lp.time_limit_ns != null)
        control.lp.clock_io orelse std.Io.Threaded.global_single_threaded.io()
    else
        null;
    const start_ns = if (clock_io) |io| std.Io.Clock.awake.now(io).nanoseconds else 0;

    if (control.initial_solution) |start| {
        if (try incumbent_candidate.copyIfFeasible(allocator, problem, start, control.integer_tolerance)) |candidate| {
            incumbent = candidate.primal;
            incumbent_objective = candidate.objective;
            have_incumbent = true;
            solutions = 1;
            if (solutions >= control.solution_limit)
                return finish(allocator, .solution_limit, incumbent, incumbent_objective, objective_bounds.unknown(problem.relaxation.objective_sense), nodes, iterations, work, solutions, &lp_relaxation, pseudocost.observations);
        }
    }

    while (open_nodes.len() != 0) {
        if (nodes >= control.node_limit)
            return finish(allocator, .node_limit, incumbent, incumbent_objective, objective_bounds.bestOpen(problem.relaxation.objective_sense, open_nodes.items(), incumbent_objective, have_incumbent), nodes, iterations, work, solutions, &lp_relaxation, pseudocost.observations);
        if (iterations >= control.lp.max_iterations)
            return finish(allocator, .iteration_limit, incumbent, incumbent_objective, objective_bounds.bestOpen(problem.relaxation.objective_sense, open_nodes.items(), incumbent_objective, have_incumbent), nodes, iterations, work, solutions, &lp_relaxation, pseudocost.observations);
        if (control.lp.work_limit) |limit| if (work >= limit)
            return finish(allocator, .work_limit, incumbent, incumbent_objective, objective_bounds.bestOpen(problem.relaxation.objective_sense, open_nodes.items(), incumbent_objective, have_incumbent), nodes, iterations, work, solutions, &lp_relaxation, pseudocost.observations);
        if (remainingTime(control.lp.time_limit_ns, clock_io, start_ns)) |remaining| if (remaining == 0)
            return finish(allocator, .time_limit, incumbent, incumbent_objective, objective_bounds.bestOpen(problem.relaxation.objective_sense, open_nodes.items(), incumbent_objective, have_incumbent), nodes, iterations, work, solutions, &lp_relaxation, pseudocost.observations);

        var node = (try open_nodes.take(allocator, problem.relaxation.objective_sense, node_selector)) orelse break;
        defer node.deinit(allocator);
        if (!current_domain.load(&change_arena, node.change_head)) continue;
        nodes += 1;
        if (control.propagation) {
            const propagated = try propagation.propagate(
                &propagation_workspace.?,
                allocator,
                problem.relaxation,
                problem.integrality,
                current_domain.lower,
                current_domain.upper,
                &change_arena,
                node.change_head,
                control.integer_tolerance,
                control.max_propagation_rounds,
            );
            if (!propagated.feasible) continue;
            node.change_head = propagated.change_head;
        }
        if (clique_table) |*table| {
            const propagated = try table.propagate(
                allocator,
                current_domain.lower,
                current_domain.upper,
                &change_arena,
                node.change_head,
                control.integer_tolerance,
            );
            lp_relaxation.recordCliquePropagation(propagated.tightenings);
            if (!propagated.feasible) continue;
            node.change_head = propagated.change_head;
        }
        if (control.conflict_learning and conflict_pool.activeCount() != 0) {
            const propagated = try conflict_pool.propagate(
                allocator,
                current_domain.lower,
                current_domain.upper,
                &change_arena,
                node.change_head,
                current_domain.root_revision,
                control.integer_tolerance,
                control.max_conflict_propagation_rounds,
            );
            const expired = conflict_pool.age(control.conflict_max_age);
            const compacted = if (expired != 0 and conflict_pool.shouldCompact()) try conflict_pool.compact() else 0;
            lp_relaxation.recordConflictPropagation(propagated.checks, propagated.tightenings, expired, compacted);
            if (!propagated.feasible) continue;
            node.change_head = propagated.change_head;
        }
        if (cut_generation) |*generation| _ = try generation.synchronize();
        var relaxation = if (cut_generation) |*generation|
            generation.view(current_domain.lower, current_domain.upper)
        else block: {
            var base = problem.relaxation;
            base.col_lower = current_domain.lower;
            base.col_upper = current_domain.upper;
            break :block base;
        };
        var lp_status: solver.LpSolveStatus = undefined;
        var solved_view: ?solver.LpSolveResultView = null;
        var separation_rounds: usize = 0;
        var node_lp_solves: usize = 0;
        while (true) {
            var node_control = control.lp;
            node_control.max_iterations = control.lp.max_iterations - iterations;
            node_control.work_limit = if (control.lp.work_limit) |limit| limit - work else null;
            node_control.initial_basis = if (nodes == 1 and node_lp_solves == 0) control.lp.initial_basis else null;
            node_control.time_limit_ns = remainingTime(control.lp.time_limit_ns, clock_io, start_ns);
            const structure_revision = if (cut_generation) |*generation| generation.structureRevision() else 0;
            lp_status = lp_relaxation.solve(relaxation, node_control, nodes, .{
                .structure = structure_revision,
                .matrix_values = structure_revision,
                .bounds = nodes,
                .objective = 0,
            });
            node_lp_solves += 1;
            iterations = std.math.add(usize, iterations, lp_relaxation.lastIterations()) catch std.math.maxInt(usize);
            work = std.math.add(u64, work, lp_relaxation.lastWork()) catch std.math.maxInt(u64);
            if (lp_status != .optimal) break;
            const current_view = lp_relaxation.resultView(relaxation, lp_status) orelse
                return finish(allocator, .numerical_failure, incumbent, incumbent_objective, 0.0, nodes, iterations, work, solutions, &lp_relaxation, pseudocost.observations);
            if (nodes == 1 and separation_rounds < control.max_cut_rounds) {
                if (cut_generation) |*generation| {
                    const round = try generation.separateRoot(
                        problem.integrality,
                        current_view.primal,
                        control.integer_tolerance,
                        control.max_cuts_per_round,
                    );
                    lp_relaxation.recordCutRound(round.cuts_added);
                    separation_rounds += 1;
                    if (round.structure_changed) {
                        relaxation = generation.view(current_domain.lower, current_domain.upper);
                        continue;
                    }
                }
            }
            solved_view = current_view;
            break;
        }

        switch (lp_status) {
            .infeasible => {
                if (control.conflict_learning) switch (try conflict_pool.learnPath(
                    &change_arena,
                    node.change_head,
                    problem.integrality,
                    current_domain.root_lower,
                    current_domain.root_upper,
                    control.integer_tolerance,
                )) {
                    .added => lp_relaxation.recordLearnedConflict(),
                    else => {},
                };
                continue;
            },
            .unbounded => return finish(allocator, if (domain.hasDiscreteColumns(problem.integrality)) .unbounded_or_infeasible else .unbounded, incumbent, incumbent_objective, 0.0, nodes, iterations, work, solutions, &lp_relaxation, pseudocost.observations),
            .iteration_limit => return finish(allocator, .iteration_limit, incumbent, incumbent_objective, node.bound orelse 0.0, nodes, iterations, work, solutions, &lp_relaxation, pseudocost.observations),
            .work_limit => return finish(allocator, .work_limit, incumbent, incumbent_objective, node.bound orelse 0.0, nodes, iterations, work, solutions, &lp_relaxation, pseudocost.observations),
            .time_limit => return finish(allocator, .time_limit, incumbent, incumbent_objective, node.bound orelse 0.0, nodes, iterations, work, solutions, &lp_relaxation, pseudocost.observations),
            .interrupted => return finish(allocator, .interrupted, incumbent, incumbent_objective, node.bound orelse 0.0, nodes, iterations, work, solutions, &lp_relaxation, pseudocost.observations),
            .invalid_model => return finish(allocator, .invalid_model, incumbent, incumbent_objective, 0.0, nodes, iterations, work, solutions, &lp_relaxation, pseudocost.observations),
            .numerical_failure, .not_implemented => return finish(allocator, .numerical_failure, incumbent, incumbent_objective, 0.0, nodes, iterations, work, solutions, &lp_relaxation, pseudocost.observations),
            .optimal => {},
        }
        const view = solved_view orelse lp_relaxation.resultView(relaxation, lp_status) orelse
            return finish(allocator, .numerical_failure, incumbent, incumbent_objective, 0.0, nodes, iterations, work, solutions, &lp_relaxation, pseudocost.observations);
        if (nodes == 1 and view.dual_valid) root_reduced_costs.capture(current_domain.lower, current_domain.upper, view.reduced_cost, view.objective_value);
        if (node.branch_origin) |origin| pseudocost.observe(problem.relaxation.objective_sense, origin, view.objective_value);
        if (nodes > 1) if (cut_generation) |*generation| {
            const expired = generation.age(view.primal, control.integer_tolerance, control.cut_max_age);
            lp_relaxation.recordExpiredCuts(expired);
        };
        if (have_incumbent and objective_bounds.cannotImprove(problem.relaxation.objective_sense, view.objective_value, incumbent_objective)) continue;
        if (have_incumbent and control.reduced_cost_fixing) {
            const global = root_reduced_costs.apply(problem.relaxation.objective_sense, problem.integrality, &current_domain, incumbent_objective, control.integer_tolerance);
            lp_relaxation.recordRootReducedCostFixing(global.tightenings);
            if (!global.feasible) continue;
        }
        if (have_incumbent and control.reduced_cost_fixing and view.dual_valid) {
            const fixed = try reduced_cost_fixing.tighten(
                allocator,
                problem.relaxation.objective_sense,
                problem.integrality,
                current_domain.lower,
                current_domain.upper,
                view.reduced_cost,
                view.objective_value,
                incumbent_objective,
                &change_arena,
                node.change_head,
                control.integer_tolerance,
            );
            lp_relaxation.recordReducedCostFixing(fixed.tightenings);
            if (!fixed.feasible) continue;
            node.change_head = fixed.change_head;
            if (fixed.tightenings != 0) if (clique_table) |*table| {
                const implied = try table.propagate(
                    allocator,
                    current_domain.lower,
                    current_domain.upper,
                    &change_arena,
                    node.change_head,
                    control.integer_tolerance,
                );
                lp_relaxation.recordCliquePropagation(implied.tightenings);
                if (!implied.feasible) continue;
                node.change_head = implied.change_head;
            };
        }

        if (strong_branching) |*engine| {
            var probe_control = control.lp;
            probe_control.initial_basis = null;
            probe_control.max_iterations = control.lp.max_iterations - iterations;
            probe_control.work_limit = if (control.lp.work_limit) |limit| limit - work else null;
            probe_control.time_limit_ns = remainingTime(control.lp.time_limit_ns, clock_io, start_ns);
            const refined = engine.refine(
                relaxation,
                problem.integrality,
                view.primal,
                view.objective_value,
                &pseudocost,
                probe_control,
                .{
                    .max_candidates = control.strong_branching_max_candidates,
                    .max_iterations_per_probe = control.strong_branching_max_iterations,
                    .max_work_per_node = control.strong_branching_max_work,
                },
                control.integer_tolerance,
            );
            iterations = std.math.add(usize, iterations, refined.statistics.iterations) catch std.math.maxInt(usize);
            work = std.math.add(u64, work, refined.statistics.work) catch std.math.maxInt(u64);
            lp_relaxation.recordStrongBranching(refined.statistics.probes, refined.statistics.iterations, refined.statistics.work);
            if (refined.limit_status) |status| switch (status) {
                .time_limit => return finish(allocator, .time_limit, incumbent, incumbent_objective, node.bound orelse view.objective_value, nodes, iterations, work, solutions, &lp_relaxation, pseudocost.observations),
                .interrupted => return finish(allocator, .interrupted, incumbent, incumbent_objective, node.bound orelse view.objective_value, nodes, iterations, work, solutions, &lp_relaxation, pseudocost.observations),
                else => {},
            };
            if (iterations >= control.lp.max_iterations)
                return finish(allocator, .iteration_limit, incumbent, incumbent_objective, node.bound orelse view.objective_value, nodes, iterations, work, solutions, &lp_relaxation, pseudocost.observations);
            if (control.lp.work_limit) |limit| if (work >= limit)
                return finish(allocator, .work_limit, incumbent, incumbent_objective, node.bound orelse view.objective_value, nodes, iterations, work, solutions, &lp_relaxation, pseudocost.observations);
        }

        const pseudocost_scores = pseudocost.prepareScores(view.primal);
        const branch_context = plugin.BranchingContext{
            .problem = problem,
            .node_lower = current_domain.lower,
            .node_upper = current_domain.upper,
            .lp_primal = view.primal,
            .active_lower = current_domain.active_lower,
            .active_upper = current_domain.active_upper,
            .tolerance = control.integer_tolerance,
            .pseudocost_scores = pseudocost_scores,
        };
        var violation: ?plugin.BranchingCandidate = null;
        for (registry.branchingRuleEntries()) |entry| {
            violation = entry.rule.select(branch_context);
            if (violation != null) break;
        }

        if (violation != null) {
            const heuristic_context = plugin.HeuristicContext{
                .allocator = allocator,
                .problem = problem,
                .node_lower = current_domain.lower,
                .node_upper = current_domain.upper,
                .lp_primal = view.primal,
                .lp_basis = if (view.basis_valid and view.basic_index.len == problem.relaxation.num_rows) .{
                    .structural_status = view.structural_status,
                    .logical_status = view.logical_status,
                    .basic_index = view.basic_index,
                } else null,
                .active_lower = current_domain.active_lower,
                .active_upper = current_domain.active_upper,
                .tolerance = control.integer_tolerance,
            };
            for (registry.heuristicEntries()) |entry| {
                const proposed = entry.heuristic.trySolution(heuristic_context) catch return error.OutOfMemory;
                const raw = proposed orelse continue;
                defer allocator.free(raw);
                if (try incumbent_candidate.copyIfFeasible(allocator, problem, raw, control.integer_tolerance)) |candidate| {
                    if (!have_incumbent or objective_bounds.betterObjective(problem.relaxation.objective_sense, candidate.objective, incumbent_objective)) {
                        allocator.free(incumbent);
                        incumbent = candidate.primal;
                        incumbent_objective = candidate.objective;
                        have_incumbent = true;
                        solutions += 1;
                        if (solutions >= control.solution_limit) {
                            const pending_bound = objective_bounds.combine(problem.relaxation.objective_sense, objective_bounds.bestOpen(problem.relaxation.objective_sense, open_nodes.items(), view.objective_value, true), view.objective_value);
                            return finish(allocator, .solution_limit, incumbent, incumbent_objective, pending_bound, nodes, iterations, work, solutions, &lp_relaxation, pseudocost.observations);
                        }
                    } else {
                        allocator.free(candidate.primal);
                    }
                }
            }
        }

        if (violation == null) {
            const candidate = try allocator.dupe(f64, view.primal);
            allocator.free(incumbent);
            incumbent = candidate;
            incumbent_objective = view.objective_value;
            have_incumbent = true;
            solutions += 1;
            if (solutions >= control.solution_limit)
                return finish(allocator, .solution_limit, incumbent, incumbent_objective, objective_bounds.bestOpen(problem.relaxation.objective_sense, open_nodes.items(), incumbent_objective, true), nodes, iterations, work, solutions, &lp_relaxation, pseudocost.observations);
            const open_bound = objective_bounds.bestOpen(problem.relaxation.objective_sense, open_nodes.items(), incumbent_objective, true);
            if (objective_bounds.gapSatisfied(open_bound, incumbent_objective, control))
                return finish(allocator, .optimal, incumbent, incumbent_objective, open_bound, nodes, iterations, work, solutions, &lp_relaxation, pseudocost.observations);
            continue;
        }

        const column = violation.?.column;
        const value = view.primal[column];
        switch (violation.?.kind) {
            .integer => {
                const down_probe = if (strong_branching) |*engine| engine.directionResult(column, .down) else null;
                const up_probe = if (strong_branching) |*engine| engine.directionResult(column, .up) else null;
                if ((down_probe == null or !down_probe.?.infeasible) and @floor(value) >= current_domain.lower[column] - control.integer_tolerance) {
                    var child = try node.childWith(&change_arena, allocator, if (down_probe) |probe| probe.bound orelse view.objective_value else view.objective_value, .{
                        .column = column,
                        .upper = @floor(value),
                    });
                    child.branch_origin = .{ .column = column, .direction = .down, .distance = value - @floor(value), .parent_objective = view.objective_value };
                    errdefer child.deinit(allocator);
                    try open_nodes.append(allocator, child);
                }
                if ((up_probe == null or !up_probe.?.infeasible) and @ceil(value) <= current_domain.upper[column] + control.integer_tolerance) {
                    var child = try node.childWith(&change_arena, allocator, if (up_probe) |probe| probe.bound orelse view.objective_value else view.objective_value, .{
                        .column = column,
                        .lower = @ceil(value),
                    });
                    child.branch_origin = .{ .column = column, .direction = .up, .distance = @ceil(value) - value, .parent_objective = view.objective_value };
                    errdefer child.deinit(allocator);
                    try open_nodes.append(allocator, child);
                }
            },
            .semi_activation => {
                if (current_domain.lower[column] <= 0 and current_domain.upper[column] >= 0) {
                    var zero = try node.childWith(&change_arena, allocator, view.objective_value, .{
                        .column = column,
                        .lower = 0,
                        .upper = 0,
                    });
                    errdefer zero.deinit(allocator);
                    try open_nodes.append(allocator, zero);
                }
                const active_lower = @max(current_domain.lower[column], current_domain.active_lower[column]);
                const active_upper = @min(current_domain.upper[column], current_domain.active_upper[column]);
                if (active_lower <= active_upper) {
                    var active = try node.childWith(&change_arena, allocator, view.objective_value, .{
                        .column = column,
                        .lower = active_lower,
                        .upper = active_upper,
                    });
                    errdefer active.deinit(allocator);
                    try open_nodes.append(allocator, active);
                }
            },
        }
    }
    if (!have_incumbent) {
        const statistics = lp_relaxation.statistics();
        return .{
            .allocator = allocator,
            .status = .infeasible,
            .nodes = nodes,
            .lp_iterations = iterations,
            .lp_work = work,
            .lp_cold_solves = statistics.cold_solves,
            .lp_reoptimizations = statistics.reoptimizations,
            .pseudocost_observations = pseudocost.observations,
            .strong_branching_probes = statistics.strong_branching_probes,
            .strong_branching_iterations = statistics.strong_branching_iterations,
            .strong_branching_work = statistics.strong_branching_work,
            .cut_rounds = statistics.cut_rounds,
            .cuts_added = statistics.cuts_added,
            .cuts_expired = statistics.cuts_expired,
            .clique_edges = statistics.clique_edges,
            .clique_tightenings = statistics.clique_tightenings,
            .reduced_cost_tightenings = statistics.reduced_cost_tightenings,
            .root_reduced_cost_tightenings = statistics.root_reduced_cost_tightenings,
            .conflicts_learned = statistics.conflicts_learned,
            .conflict_tightenings = statistics.conflict_tightenings,
            .conflicts_expired = statistics.conflicts_expired,
            .conflicts_compacted = statistics.conflicts_compacted,
            .conflict_checks = statistics.conflict_checks,
        };
    }
    return finish(allocator, .optimal, incumbent, incumbent_objective, incumbent_objective, nodes, iterations, work, solutions, &lp_relaxation, pseudocost.observations);
}

fn remainingTime(limit: ?u64, clock_io: ?std.Io, start_ns: i96) ?u64 {
    const total = limit orelse return null;
    if (total == 0) return 0;
    const io = clock_io orelse return total;
    const now = std.Io.Clock.awake.now(io).nanoseconds;
    if (now <= start_ns) return total;
    const elapsed: u128 = @intCast(now - start_ns);
    if (elapsed >= total) return 0;
    return total - @as(u64, @intCast(elapsed));
}

fn finish(allocator: std.mem.Allocator, status: SolveStatus, primal: []f64, objective: f64, bound: f64, nodes: usize, iterations: usize, work: u64, solutions: usize, relaxation: *const Relaxation, pseudocost_observations: usize) SolveResult {
    const statistics = relaxation.statistics();
    return .{
        .allocator = allocator,
        .status = status,
        .primal = primal,
        .objective_value = objective,
        .objective_bound = bound,
        .nodes = nodes,
        .lp_iterations = iterations,
        .lp_work = work,
        .lp_cold_solves = statistics.cold_solves,
        .lp_reoptimizations = statistics.reoptimizations,
        .pseudocost_observations = pseudocost_observations,
        .strong_branching_probes = statistics.strong_branching_probes,
        .strong_branching_iterations = statistics.strong_branching_iterations,
        .strong_branching_work = statistics.strong_branching_work,
        .cut_rounds = statistics.cut_rounds,
        .cuts_added = statistics.cuts_added,
        .cuts_expired = statistics.cuts_expired,
        .clique_edges = statistics.clique_edges,
        .clique_tightenings = statistics.clique_tightenings,
        .reduced_cost_tightenings = statistics.reduced_cost_tightenings,
        .root_reduced_cost_tightenings = statistics.root_reduced_cost_tightenings,
        .conflicts_learned = statistics.conflicts_learned,
        .conflict_tightenings = statistics.conflict_tightenings,
        .conflicts_expired = statistics.conflicts_expired,
        .conflicts_compacted = statistics.conflicts_compacted,
        .conflict_checks = statistics.conflict_checks,
        .solution_count = solutions,
    };
}

fn fractionalBinaryRelaxation() solver.LpProblemView {
    return .{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 1.0, 1.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ 1.0, 1.0 },
        .row_lower = &.{1.5},
        .row_upper = &.{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &.{ 0, 1, 2 }, &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) }, &.{ 1.0, 1.0 }),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
}

test "branch-and-bound solves a fractional binary relaxation" {
    const relaxation = fractionalBinaryRelaxation();
    var result = try solve(std.testing.allocator, .{ .relaxation = relaxation, .integrality = &.{ .integer, .integer } }, .{ .propagation = false });
    defer result.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, result.status);
    try std.testing.expectApproxEqAbs(@as(f64, 2), result.objective_value, 1e-9);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0 }, result.primal);
    try std.testing.expect(result.nodes >= 2);
    try std.testing.expectEqual(result.nodes, result.lp_cold_solves + result.lp_reoptimizations);
    try std.testing.expectEqual(@as(usize, 1), result.lp_cold_solves);
    try std.testing.expectEqual(result.nodes - 1, result.lp_reoptimizations);
    try std.testing.expect(result.pseudocost_observations > 0);

    var limited = try solve(std.testing.allocator, .{ .relaxation = relaxation, .integrality = &.{ .integer, .integer } }, .{ .propagation = false, .lp = .{ .max_iterations = 1 } });
    defer limited.deinit();
    try std.testing.expectEqual(SolveStatus.iteration_limit, limited.status);
    try std.testing.expect(limited.lp_iterations <= 1);

    var work_limited = try solve(std.testing.allocator, .{ .relaxation = relaxation, .integrality = &.{ .integer, .integer } }, .{ .lp = .{ .work_limit = 0 } });
    defer work_limited.deinit();
    try std.testing.expectEqual(SolveStatus.work_limit, work_limited.status);
    try std.testing.expectEqual(@as(u64, 0), work_limited.lp_work);
    try std.testing.expectEqual(@as(usize, 0), work_limited.nodes);

    var time_limited = try solve(std.testing.allocator, .{ .relaxation = relaxation, .integrality = &.{ .integer, .integer } }, .{ .lp = .{ .time_limit_ns = 0 } });
    defer time_limited.deinit();
    try std.testing.expectEqual(SolveStatus.time_limit, time_limited.status);

    var heuristic = try solve(std.testing.allocator, .{ .relaxation = relaxation, .integrality = &.{ .integer, .integer } }, .{ .solution_limit = 1 });
    defer heuristic.deinit();
    try std.testing.expectEqual(SolveStatus.solution_limit, heuristic.status);
    try std.testing.expectEqual(@as(usize, 1), heuristic.nodes);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0 }, heuristic.primal);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), heuristic.objective_value, 1e-9);
}

test "row propagation can close the root domain before branching" {
    const relaxation = fractionalBinaryRelaxation();
    var result = try solve(std.testing.allocator, .{ .relaxation = relaxation, .integrality = &.{ .binary, .binary } }, .{});
    defer result.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, result.status);
    try std.testing.expectEqual(@as(usize, 1), result.nodes);
    try std.testing.expectEqual(@as(usize, 1), result.lp_cold_solves);
    try std.testing.expectEqual(@as(usize, 0), result.lp_reoptimizations);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0 }, result.primal);
}

test "clique propagation rejects a conflicting child before its LP solve" {
    const relaxation = solver.LpProblemView{
        .num_rows = 3,
        .num_cols = 2,
        .col_cost = &.{ 0.0, 0.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ 1.0, 1.0 },
        .row_lower = &.{ -std.math.inf(f64), 0.5, 0.5 },
        .row_upper = &.{ 1.0, std.math.inf(f64), std.math.inf(f64) },
        .matrix = matrix.CscView.initAssumeValid(3, 2, &.{ 0, 2, 4 }, &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(2) }, &.{ 1.0, 1.0, 1.0, 1.0 }),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    const problem = ProblemView{ .relaxation = relaxation, .integrality = &.{ .binary, .binary } };
    var with_cliques = try solve(std.testing.allocator, problem, .{
        .propagation = false,
        .cut_separation = false,
        .reliability_branching = false,
    });
    defer with_cliques.deinit();
    var without_cliques = try solve(std.testing.allocator, problem, .{
        .propagation = false,
        .clique_propagation = false,
        .cut_separation = false,
        .reliability_branching = false,
    });
    defer without_cliques.deinit();
    try std.testing.expectEqual(SolveStatus.infeasible, with_cliques.status);
    try std.testing.expectEqual(with_cliques.status, without_cliques.status);
    try std.testing.expectEqual(@as(usize, 1), with_cliques.clique_edges);
    try std.testing.expect(with_cliques.clique_tightenings > 0);
    try std.testing.expectEqual(@as(usize, 0), without_cliques.clique_edges);
    try std.testing.expectEqual(@as(usize, 0), without_cliques.clique_tightenings);
    try std.testing.expect(without_cliques.conflicts_learned > 0);
}

test "an incumbent activates reduced-cost domain tightening" {
    const relaxation = solver.LpProblemView{
        .num_rows = 1,
        .num_cols = 2,
        .col_cost = &.{ 2.0, 1.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ 10.0, 1.0 },
        .row_lower = &.{0.5},
        .row_upper = &.{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(1, 2, &.{ 0, 0, 1 }, &.{foundation.RowId.fromUsizeAssumeValid(0)}, &.{1.0}),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    const problem = ProblemView{ .relaxation = relaxation, .integrality = &.{ .integer, .integer } };
    var with_fixing = try solve(std.testing.allocator, problem, .{ .initial_solution = &.{ 2.0, 1.0 } });
    defer with_fixing.deinit();
    var without_fixing = try solve(std.testing.allocator, problem, .{
        .initial_solution = &.{ 2.0, 1.0 },
        .reduced_cost_fixing = false,
    });
    defer without_fixing.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, with_fixing.status);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), with_fixing.objective_value, 1e-9);
    try std.testing.expectApproxEqAbs(with_fixing.objective_value, without_fixing.objective_value, 1e-9);
    try std.testing.expect(with_fixing.reduced_cost_tightenings + with_fixing.root_reduced_cost_tightenings > 0);
    try std.testing.expectEqual(@as(usize, 0), without_fixing.reduced_cost_tightenings);
    try std.testing.expectEqual(@as(usize, 0), without_fixing.root_reduced_cost_tightenings);
}

test "root cover separation closes a binary knapsack without branching" {
    const relaxation = solver.LpProblemView{
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
    const problem = ProblemView{ .relaxation = relaxation, .integrality = &.{ .binary, .binary } };
    var with_cuts = try solve(std.testing.allocator, problem, .{ .propagation = false, .reliability_branching = false });
    defer with_cuts.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, with_cuts.status);
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), with_cuts.objective_value, 1e-9);
    try std.testing.expectEqual(@as(usize, 1), with_cuts.nodes);
    try std.testing.expectEqual(@as(usize, 1), with_cuts.cuts_added);
    try std.testing.expectEqual(@as(usize, 2), with_cuts.cut_rounds);
    try std.testing.expectEqual(@as(usize, 2), with_cuts.lp_cold_solves);

    var without_cuts = try solve(std.testing.allocator, problem, .{ .propagation = false, .reliability_branching = false, .cut_separation = false });
    defer without_cuts.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, without_cuts.status);
    try std.testing.expectApproxEqAbs(with_cuts.objective_value, without_cuts.objective_value, 1e-9);
    try std.testing.expect(without_cuts.nodes > with_cuts.nodes);
    try std.testing.expectEqual(@as(usize, 0), without_cuts.cuts_added);
}

test "reliability branching probes multiple fractional candidates" {
    const relaxation = solver.LpProblemView{
        .num_rows = 2,
        .num_cols = 2,
        .col_cost = &.{ 1.0, 1.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ 2.0, 2.0 },
        .row_lower = &.{ 0.5, 0.5 },
        .row_upper = &.{ std.math.inf(f64), std.math.inf(f64) },
        .matrix = matrix.CscView.initAssumeValid(2, 2, &.{ 0, 1, 2 }, &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1) }, &.{ 1.0, 1.0 }),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    var result = try solve(std.testing.allocator, .{ .relaxation = relaxation, .integrality = &.{ .integer, .integer } }, .{
        .propagation = false,
        .strong_branching_max_candidates = 1,
    });
    defer result.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, result.status);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), result.objective_value, 1e-9);
    try std.testing.expectEqual(@as(usize, 2), result.strong_branching_probes);
    try std.testing.expect(result.strong_branching_iterations <= result.lp_iterations);
    try std.testing.expect(result.strong_branching_work <= result.lp_work);
}

test "branch-and-bound reoptimizes stable node LP mappings" {
    const relaxation = solver.LpProblemView{
        .num_rows = 2,
        .num_cols = 2,
        .col_cost = &.{ 1.0, 0.0 },
        .col_lower = &.{ 0.0, 0.0 },
        .col_upper = &.{ 10.0, 10.0 },
        .row_lower = &.{ -std.math.inf(f64), -std.math.inf(f64) },
        .row_upper = &.{ 4.5, 9.0 },
        .matrix = matrix.CscView.initAssumeValid(2, 2, &.{ 0, 2, 4 }, &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1) }, &.{ 1.0, 1.0, 1.0, -1.0 }),
        .objective_sense = .maximize,
        .objective_offset = 0,
    };
    var result = try solve(std.testing.allocator, .{ .relaxation = relaxation, .integrality = &.{ .integer, .continuous } }, .{ .propagation = false });
    defer result.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, result.status);
    try std.testing.expectApproxEqAbs(@as(f64, 4), result.objective_value, 1e-9);
    try std.testing.expect(result.lp_reoptimizations >= 1);
    try std.testing.expectEqual(result.nodes, result.lp_cold_solves + result.lp_reoptimizations);
}

test "branch-and-bound proves integer infeasibility" {
    const relaxation = solver.LpProblemView{
        .num_rows = 1,
        .num_cols = 1,
        .col_cost = &.{0.0},
        .col_lower = &.{0.0},
        .col_upper = &.{1.0},
        .row_lower = &.{0.5},
        .row_upper = &.{0.5},
        .matrix = matrix.CscView.initAssumeValid(1, 1, &.{ 0, 1 }, &.{foundation.RowId.fromUsizeAssumeValid(0)}, &.{1.0}),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var result = try solve(std.testing.allocator, .{ .relaxation = relaxation, .integrality = &.{.integer} }, .{});
    defer result.deinit();
    try std.testing.expectEqual(SolveStatus.infeasible, result.status);
}

test "plugin components can override node selection and branching" {
    const FirstFractional = struct {
        pub const name = "first-fractional-test";
        pub fn select(_: *const @This(), ctx: plugin.BranchingContext) ?plugin.BranchingCandidate {
            for (ctx.problem.integrality, ctx.lp_primal, 0..) |kind, value, column| {
                if (kind == .continuous or kind == .semi_continuous) continue;
                if (@abs(value - @round(value)) > ctx.tolerance) return .{ .column = column, .kind = .integer };
            }
            return null;
        }
    };
    const Fifo = struct {
        pub const name = "fifo-test";
        pub fn select(_: *const @This(), ctx: plugin.NodeSelectionContext) ?usize {
            return if (ctx.open.len == 0) null else 0;
        }
    };
    var branching = FirstFractional{};
    var selector = Fifo{};
    var result = try solve(std.testing.allocator, .{
        .relaxation = fractionalBinaryRelaxation(),
        .integrality = &.{ .integer, .integer },
    }, .{ .components = .{
        .branching_rule = plugin.bindBranching(FirstFractional, &branching),
        .node_selector = plugin.bindNodeSelector(Fifo, &selector),
    } });
    defer result.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, result.status);
    try std.testing.expectApproxEqAbs(@as(f64, 2), result.objective_value, 1e-9);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0 }, result.primal);
}

test "pseudocost branching stays competitive with most-fractional branching" {
    const relaxation = solver.LpProblemView{
        .num_rows = 2,
        .num_cols = 4,
        .col_cost = &.{ 1.0, 2.0, 3.0, 4.0 },
        .col_lower = &.{ 0.0, 0.0, 0.0, 0.0 },
        .col_upper = &.{ 1.0, 1.0, 1.0, 1.0 },
        .row_lower = &.{ 1.5, 1.5 },
        .row_upper = &.{ std.math.inf(f64), std.math.inf(f64) },
        .matrix = matrix.CscView.initAssumeValid(2, 4, &.{ 0, 1, 2, 3, 4 }, &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(1), foundation.RowId.fromUsizeAssumeValid(1) }, &.{ 1.0, 1.0, 1.0, 1.0 }),
        .objective_sense = .minimize,
        .objective_offset = 0.0,
    };
    const problem = ProblemView{ .relaxation = relaxation, .integrality = &.{ .binary, .binary, .binary, .binary } };
    var most_fractional = builtin_plugins.branching.MostFractional{};
    var reference = try solve(std.testing.allocator, problem, .{
        .propagation = false,
        .reliability_branching = false,
        .cut_separation = false,
        .components = .{ .branching_rule = plugin.bindBranching(builtin_plugins.branching.MostFractional, &most_fractional) },
    });
    defer reference.deinit();
    var pseudocost = try solve(std.testing.allocator, problem, .{
        .propagation = false,
        .reliability_branching = false,
        .cut_separation = false,
    });
    defer pseudocost.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, reference.status);
    try std.testing.expectEqual(SolveStatus.optimal, pseudocost.status);
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), reference.objective_value, 1e-9);
    try std.testing.expectApproxEqAbs(reference.objective_value, pseudocost.objective_value, 1e-9);
    try std.testing.expect(pseudocost.pseudocost_observations > 0);
    try std.testing.expect(pseudocost.nodes <= reference.nodes);
}

test "best-bound node selection preserves optimal objective" {
    var selector = builtin_plugins.node_selectors.BestBound{};
    var result = try solve(std.testing.allocator, .{
        .relaxation = fractionalBinaryRelaxation(),
        .integrality = &.{ .integer, .integer },
    }, .{ .components = .{
        .node_selector = plugin.bindNodeSelector(builtin_plugins.node_selectors.BestBound, &selector),
    } });
    defer result.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, result.status);
    try std.testing.expectApproxEqAbs(@as(f64, 2), result.objective_value, 1e-9);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0 }, result.primal);
}

test "disabling the rounding heuristic keeps the search correct" {
    var result = try solve(std.testing.allocator, .{
        .relaxation = fractionalBinaryRelaxation(),
        .integrality = &.{ .integer, .integer },
    }, .{ .rounding_heuristic = false });
    defer result.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, result.status);
    try std.testing.expectApproxEqAbs(@as(f64, 2), result.objective_value, 1e-9);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0 }, result.primal);
}

test "the default diving heuristic finds an incumbent and stays optimal" {
    // min x1 + x2 + x3 subject to x1 + x2 + x3 >= 2.5, all binary. The root
    // LP point is (5/6, 5/6, 5/6); diving repeatedly fixes the rounded column
    // and re-solves to the incumbent (1, 1, 1), which is optimal with
    // objective 3.
    const relaxation = solver.LpProblemView{
        .num_rows = 1,
        .num_cols = 3,
        .col_cost = &.{ 1.0, 1.0, 1.0 },
        .col_lower = &.{ 0.0, 0.0, 0.0 },
        .col_upper = &.{ 1.0, 1.0, 1.0 },
        .row_lower = &.{2.5},
        .row_upper = &.{std.math.inf(f64)},
        .matrix = matrix.CscView.initAssumeValid(1, 3, &.{ 0, 1, 2, 3 }, &.{ foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0), foundation.RowId.fromUsizeAssumeValid(0) }, &.{ 1.0, 1.0, 1.0 }),
        .objective_sense = .minimize,
        .objective_offset = 0,
    };
    var result = try solve(std.testing.allocator, .{
        .relaxation = relaxation,
        .integrality = &.{ .binary, .binary, .binary },
    }, .{});
    defer result.deinit();
    try std.testing.expectEqual(SolveStatus.optimal, result.status);
    try std.testing.expectApproxEqAbs(@as(f64, 3), result.objective_value, 1e-9);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 1.0, 1.0 }, result.primal);
}