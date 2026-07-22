//! End-to-end MPS parser + revised-simplex acceptance runner.
//!
//! Output is a single TSV row so corpus scripts can compare status, objective,
//! iteration count, residuals, factorization lifecycle, and elapsed solve time.

const std = @import("std");
const zhighs = @import("zhighs");

fn nowNs() i128 {
    var timestamp: std.posix.timespec = undefined;
    switch (std.posix.errno(std.posix.system.clock_gettime(.MONOTONIC, &timestamp))) {
        .SUCCESS => {},
        else => unreachable,
    }
    return @as(i128, timestamp.sec) * std.time.ns_per_s + timestamp.nsec;
}

fn violation(value: f64, lower: f64, upper: f64) f64 {
    return @max(@max(lower - value, value - upper), 0.0);
}

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;
    const io_context = init.io;
    var args = std.process.Args.iterate(init.minimal.args);
    _ = args.next();
    const path = args.next() orelse return error.MissingModelPath;
    const max_iterations = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 1_000_000;
    const max_updates = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 100;
    const trace_enabled = if (args.next()) |text| std.mem.eql(u8, text, "trace") else false;
    const degenerate_limit = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 8;
    const refinement_steps = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 2;
    const sparse_threshold = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 64;
    const collect_statistics = if (args.next()) |text| std.mem.eql(u8, text, "stats") else false;
    const fresh_recovery_pivots = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 32;
    const phase_one_strategy: zhighs.lp.simplex.engine.PhaseOneStrategy = if (args.next()) |text|
        if (std.mem.eql(u8, text, "dual")) .dual else if (std.mem.eql(u8, text, "auto")) .automatic else if (std.mem.eql(u8, text, "primal")) .primal else return error.InvalidArguments
    else
        .primal;
    const crash_strategy: zhighs.lp.simplex.engine.CrashStrategy = if (args.next()) |text|
        if (std.mem.eql(u8, text, "ltssf")) .ltssf else if (std.mem.eql(u8, text, "bixby")) .bixby else if (std.mem.eql(u8, text, "auto")) .automatic else if (std.mem.eql(u8, text, "logical")) .logical else return error.InvalidArguments
    else
        .logical;
    const crash_max_columns = if (args.next()) |text| blk: {
        const value = try std.fmt.parseUnsigned(usize, text, 10);
        break :blk if (value == 0) null else value;
    } else null;
    const degeneracy_strategy: zhighs.lp.simplex.engine.DegeneracyStrategy = if (args.next()) |text|
        if (std.mem.eql(u8, text, "perturb")) .perturbation else if (std.mem.eql(u8, text, "taboo")) .perturbation_taboo else if (std.mem.eql(u8, text, "auto")) .automatic else if (std.mem.eql(u8, text, "baseline")) .baseline else return error.InvalidArguments
    else
        .automatic;
    const phase_one_pricing: zhighs.lp.simplex.engine.PhaseOnePricingStrategy = if (args.next()) |text|
        if (std.mem.eql(u8, text, "dantzig")) .dantzig else if (std.mem.eql(u8, text, "devex")) .devex else if (std.mem.eql(u8, text, "steepest")) .steepest_edge else if (std.mem.eql(u8, text, "inherit")) .inherit else return error.InvalidArguments
    else
        .inherit;
    const adaptive_reprice = if (args.next()) |text|
        if (std.mem.eql(u8, text, "adaptive")) true else if (std.mem.eql(u8, text, "fixed")) false else return error.InvalidArguments
    else
        false;
    const pricing_kernel: zhighs.lp.simplex.engine.PricingKernel = if (args.next()) |text|
        if (std.mem.eql(u8, text, "row")) .row else if (std.mem.eql(u8, text, "auto")) .automatic else if (std.mem.eql(u8, text, "column")) .column else return error.InvalidArguments
    else
        .column;
    const devex_strategy: zhighs.lp.simplex.engine.DevexStrategy = if (args.next()) |text|
        if (std.mem.eql(u8, text, "framework")) .framework else if (std.mem.eql(u8, text, "legacy")) .legacy else return error.InvalidArguments
    else
        .legacy;
    const dual_edge_weight_strategy: zhighs.lp.simplex.engine.DualEdgeWeightStrategy = if (args.next()) |text|
        if (std.mem.eql(u8, text, "steepest-devex")) .steepest_devex else if (std.mem.eql(u8, text, "inherit")) .inherit else return error.InvalidArguments
    else
        .steepest_devex;
    const dual_dse_update_budget = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 64;
    const primal_pricing_strategy: zhighs.lp.simplex.engine.PrimalPricingStrategy = if (args.next()) |text|
        if (std.mem.eql(u8, text, "partial")) .partial else if (std.mem.eql(u8, text, "inherit")) .inherit else return error.InvalidArguments
    else
        .inherit;
    if (args.next() != null) return error.InvalidArguments;

    const started = nowNs();
    var model = try zhighs.io.readFile(io_context, allocator, path, .{ .keep_names = false });
    defer model.deinit();
    const parsed_ns: u64 = @intCast(nowNs() - started);
    const problem = zhighs.lp.ProblemView{
        .num_rows = model.row_lower.len,
        .num_cols = model.col_cost.len,
        .col_cost = model.col_cost,
        .col_lower = model.col_lower,
        .col_upper = model.col_upper,
        .row_lower = model.row_lower,
        .row_upper = model.row_upper,
        .matrix = model.matrix.view(),
        .objective_sense = if (model.objective_sense == .minimize) .minimize else .maximize,
        .objective_offset = model.objective_offset,
    };
    if (trace_enabled) {
        var nonzero_costs: usize = 0;
        var nonzero_row_bounds: usize = 0;
        for (model.col_cost) |value| if (value != 0.0) {
            nonzero_costs += 1;
        };
        for (model.row_lower, model.row_upper) |lower, upper| if ((std.math.isFinite(lower) and lower != 0.0) or
            (std.math.isFinite(upper) and upper != 0.0))
        {
            nonzero_row_bounds += 1;
        };
        const meta = try std.fmt.allocPrint(allocator, "model\t{s}\tcost_nz={d}\tbound_nz={d}\n", .{
            @tagName(model.objective_sense), nonzero_costs, nonzero_row_bounds,
        });
        defer allocator.free(meta);
        try std.Io.File.stderr().writeStreamingAll(io_context, meta);
    }
    var engine = zhighs.lp.simplex.engine.SimplexEngine.init(allocator);
    defer engine.deinit();
    engine.numerical.max_update_count = max_updates;
    engine.numerical.degenerate_pivot_limit = degenerate_limit;
    engine.numerical.max_refinement_steps = refinement_steps;
    engine.numerical.fresh_factorization_recovery_pivots = fresh_recovery_pivots;
    engine.factorization.sparse_dimension_threshold = sparse_threshold;
    const trace: []zhighs.lp.simplex.engine.PivotTraceEvent = if (trace_enabled)
        try allocator.alloc(zhighs.lp.simplex.engine.PivotTraceEvent, @min(max_iterations, 100_000))
    else
        &.{};
    defer if (trace_enabled) allocator.free(trace);
    const degeneracy_trace: []zhighs.lp.simplex.engine.DegeneracyTraceEvent = if (trace_enabled)
        try allocator.alloc(zhighs.lp.simplex.engine.DegeneracyTraceEvent, @min(max_iterations, 100_000))
    else
        &.{};
    defer if (trace_enabled) allocator.free(degeneracy_trace);
    const dual_phase_one_candidates: []zhighs.lp.simplex.engine.DualPhaseOneCandidateTraceEvent = if (trace_enabled)
        try allocator.alloc(zhighs.lp.simplex.engine.DualPhaseOneCandidateTraceEvent, problem.num_cols + problem.num_rows)
    else
        &.{};
    defer if (trace_enabled) allocator.free(dual_phase_one_candidates);
    const dual_phase_one_ep: []zhighs.lp.simplex.engine.DualPhaseOneEpTraceEvent = if (trace_enabled)
        try allocator.alloc(zhighs.lp.simplex.engine.DualPhaseOneEpTraceEvent, problem.num_rows)
    else
        &.{};
    defer if (trace_enabled) allocator.free(dual_phase_one_ep);
    const solve_started = nowNs();
    const status = engine.solveProblem(problem, .{
        .max_iterations = max_iterations,
        .pivot_trace = trace,
        .degeneracy_trace = degeneracy_trace,
        .dual_phase_one_candidate_trace = dual_phase_one_candidates,
        .dual_phase_one_ep_trace = dual_phase_one_ep,
        .collect_statistics = collect_statistics,
        .phase_one_strategy = phase_one_strategy,
        .crash_strategy = crash_strategy,
        .crash_max_columns = crash_max_columns,
        .degeneracy_strategy = degeneracy_strategy,
        .phase_one_pricing = phase_one_pricing,
        .adaptive_reprice = adaptive_reprice,
        .pricing_kernel = pricing_kernel,
        .devex_strategy = devex_strategy,
        .dual_edge_weight_strategy = dual_edge_weight_strategy,
        .dual_dse_update_budget = dual_dse_update_budget,
        .primal_pricing_strategy = primal_pricing_strategy,
    });
    const solve_ns: u64 = @intCast(nowNs() - solve_started);
    if (trace_enabled) for (trace[0..engine.pivot_trace_count]) |event| {
        const trace_line = try std.fmt.allocPrint(allocator, "pivot\t{s}\t{d}\t{d}\t{d}\t{d}\t{d:.17}\t{d:.17}\t{d}\t{e:.6}\t{e:.6}\t{d}\n", .{
            @tagName(event.phase),  event.iteration, event.entering_column, event.leaving_column,          event.leaving_row,
            event.pivot,            event.step,      event.update_count,    event.ftran_relative_residual, event.condition_estimate,
            event.bound_flip_count,
        });
        defer allocator.free(trace_line);
        try std.Io.File.stderr().writeStreamingAll(io_context, trace_line);
    };
    if (trace_enabled) if (engine.basis) |basis| {
        for (basis.row_scale, 0..) |scale, row| {
            const line = try std.fmt.allocPrint(allocator, "scale_row\t{d}\t{e:.17}\n", .{ row, scale });
            defer allocator.free(line);
            try std.Io.File.stderr().writeStreamingAll(io_context, line);
        }
        for (basis.column_scale[0..problem.num_cols], 0..) |scale, column| {
            const line = try std.fmt.allocPrint(allocator, "scale_column\t{d}\t{e:.17}\n", .{ column, scale });
            defer allocator.free(line);
            try std.Io.File.stderr().writeStreamingAll(io_context, line);
        }
    };
    if (trace_enabled) for (degeneracy_trace[0..engine.degeneracy_trace_count]) |event| {
        const trace_line = try std.fmt.allocPrint(allocator, "degenerate\t{s}\t{d}\t{s}\t{d}\t{d}\t{e:.6}\t{e:.6}\t{x}\n", .{
            @tagName(event.phase),
            event.iteration,
            @tagName(event.reason),
            event.entering_column,
            event.leaving_column,
            event.step,
            event.objective_change,
            event.basis_fingerprint,
        });
        defer allocator.free(trace_line);
        try std.Io.File.stderr().writeStreamingAll(io_context, trace_line);
    };
    if (trace_enabled) if (engine.dual_phase_one_failure) |failure| {
        const failure_line = try std.fmt.allocPrint(allocator, "dual_phase1_failure\titeration={d}\tleaving_row={d}\tleaving_column={d}\tleaving_bound={s}\tviolation={e:.6}\tvalue={e:.17}\tworking_lower={e:.17}\tworking_upper={e:.17}\toriginal_lower={e:.17}\toriginal_upper={e:.17}\tep_nonzeros={d}\tep_max_abs={e:.6}\tsmall_tableau={d}\tbasic_or_fixed={d}\twrong_pivot_sign={d}\taccepted_bound_flips={d}\teligible_unselected={d}\n", .{
            failure.iteration,
            failure.leaving_row,
            failure.leaving_column,
            @tagName(failure.leaving_bound),
            failure.violation,
            failure.leaving_value,
            failure.working_lower,
            failure.working_upper,
            failure.original_lower,
            failure.original_upper,
            failure.ep_nonzeros,
            failure.ep_max_abs,
            failure.small_tableau,
            failure.basic_or_fixed,
            failure.wrong_pivot_sign,
            failure.accepted_bound_flips,
            failure.eligible_unselected,
        });
        defer allocator.free(failure_line);
        try std.Io.File.stderr().writeStreamingAll(io_context, failure_line);
        for (dual_phase_one_ep[0..engine.dual_phase_one_ep_trace_count]) |event| {
            const line = try std.fmt.allocPrint(allocator, "dual_phase1_ep\t{d}\t{e:.17}\n", .{ event.row, event.value });
            defer allocator.free(line);
            try std.Io.File.stderr().writeStreamingAll(io_context, line);
        }
        for (dual_phase_one_candidates[0..engine.dual_phase_one_candidate_trace_count]) |event| {
            const line = try std.fmt.allocPrint(allocator, "dual_phase1_candidate\t{d}\t{s}\t{d}\t{e:.17}\t{e:.17}\t{e:.17}\t{e:.17}\t{e:.17}\t{e:.17}\t{e:.17}\t{e:.17}\t{s}\n", .{
                event.column,
                @tagName(event.status),
                event.explicit_move,
                event.tableau,
                event.direction,
                event.signed_pivot,
                event.reduced_cost,
                event.lower,
                event.upper,
                event.primal,
                event.flip_capacity,
                @tagName(event.reason),
            });
            defer allocator.free(line);
            try std.Io.File.stderr().writeStreamingAll(io_context, line);
        }
    };

    var primal_residual: f64 = std.math.inf(f64);
    var dual_residual: f64 = std.math.inf(f64);
    var ray_residual: f64 = std.math.inf(f64);
    var ray_objective: f64 = 0.0;
    if (engine.solutionView(problem, status)) |solution| {
        primal_residual = 0.0;
        dual_residual = 0.0;
        for (solution.primal, model.col_lower, model.col_upper) |value, lower, upper|
            primal_residual = @max(primal_residual, violation(value, lower, upper));
        const activity = try allocator.alloc(f64, problem.num_rows);
        defer allocator.free(activity);
        @memset(activity, 0.0);
        for (0..problem.num_cols) |column| {
            const begin = problem.matrix.col_starts[column];
            const end = problem.matrix.col_starts[column + 1];
            for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient|
                activity[row.toUsize()] += coefficient * solution.primal[column];
        }
        for (activity, model.row_lower, model.row_upper) |value, lower, upper|
            primal_residual = @max(primal_residual, violation(value, lower, upper));
        if (engine.basis) |basis| for (solution.reduced_cost, basis.col_status[0..problem.num_cols]) |cost, column_status| {
            const current = switch (column_status) {
                .at_lower => @max(-cost, 0.0),
                .at_upper => @max(cost, 0.0),
                .basic, .free, .superbasic => @abs(cost),
                .fixed => 0.0,
            };
            dual_residual = @max(dual_residual, current);
        };
        if (solution.unbounded_ray.len == problem.num_cols) {
            ray_residual = 0.0;
            @memset(activity, 0.0);
            for (0..problem.num_cols) |column| {
                const direction = solution.unbounded_ray[column];
                ray_objective += problem.col_cost[column] * direction;
                if (std.math.isFinite(problem.col_upper[column])) ray_residual = @max(ray_residual, direction);
                if (std.math.isFinite(problem.col_lower[column])) ray_residual = @max(ray_residual, -direction);
                const begin = problem.matrix.col_starts[column];
                const end = problem.matrix.col_starts[column + 1];
                for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
                    if (!zhighs.matrix.MatrixTargetPolicy.retainsModelCoefficient(coefficient)) continue;
                    activity[row.toUsize()] += coefficient * direction;
                }
            }
            for (activity, model.row_lower, model.row_upper) |direction, lower, upper| {
                if (std.math.isFinite(upper)) ray_residual = @max(ray_residual, direction);
                if (std.math.isFinite(lower)) ray_residual = @max(ray_residual, -direction);
            }
            ray_residual = @max(ray_residual, 0.0);
        }
    }

    const stats = engine.factorization.stats;
    const reinversions = stats.update_limit_reinversions + stats.update_growth_reinversions + stats.solve_residual_reinversions;
    const total_ns = parsed_ns + solve_ns;
    const line = try std.fmt.allocPrint(
        allocator,
        "zhighs\t{s}\t{s}\t{s}\t{d:.17}\t{d}\t{e:.6}\t{e:.6}\t{e:.6}\t{e:.6}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\n",
        .{ path, @tagName(status), @tagName(engine.failure_site), engine.objective_value, engine.iteration_counters.committed_pivots, primal_residual, dual_residual, ray_residual, ray_objective, stats.factorizations, reinversions, stats.update_limit_reinversions, stats.update_growth_reinversions, stats.ft_updates, engine.factorization.update_count, parsed_ns, solve_ns, total_ns },
    );
    defer allocator.free(line);
    try std.Io.File.stdout().writeStreamingAll(io_context, line);
    if (!collect_statistics) return;

    const simplex_stats = engine.stats;
    const factor_stats = engine.factorization.stats;
    const aq_density = if (simplex_stats.aq_samples == 0 or problem.num_rows == 0)
        0.0
    else
        @as(f64, @floatFromInt(simplex_stats.aq_nonzeros)) /
            @as(f64, @floatFromInt(simplex_stats.aq_samples * problem.num_rows));
    const ep_density = if (simplex_stats.ep_samples == 0 or problem.num_rows == 0)
        0.0
    else
        @as(f64, @floatFromInt(simplex_stats.ep_nonzeros)) /
            @as(f64, @floatFromInt(simplex_stats.ep_samples * problem.num_rows));
    const pricing_density = if (simplex_stats.pricing_entries == 0)
        0.0
    else
        @as(f64, @floatFromInt(simplex_stats.pricing_nonzeros)) /
            @as(f64, @floatFromInt(simplex_stats.pricing_entries));
    const ftran_rhs_density = if (factor_stats.ftran_rhs_samples == 0 or problem.num_rows == 0)
        0.0
    else
        @as(f64, @floatFromInt(factor_stats.ftran_rhs_nonzeros)) /
            @as(f64, @floatFromInt(factor_stats.ftran_rhs_samples * problem.num_rows));
    const btran_rhs_density = if (factor_stats.btran_rhs_samples == 0 or problem.num_rows == 0)
        0.0
    else
        @as(f64, @floatFromInt(factor_stats.btran_rhs_nonzeros)) /
            @as(f64, @floatFromInt(factor_stats.btran_rhs_samples * problem.num_rows));
    const requested_bytes = engine.requestedBytes();
    const peak_rss_kb: u64 = @intCast(std.posix.getrusage(std.posix.rusage.SELF).maxrss);
    const stats_line = try std.fmt.allocPrint(
        allocator,
        "stats\t{s}\tphase1_iterations={d}\tdual_phase1_iterations={d}\tdual_phase1_fallbacks={d}\tdual_repair_iterations={d}\tphase2_iterations={d}\tphase1_ns={d}\tdual_phase1_ns={d}\tdual_repair_ns={d}\tphase2_ns={d}\tcleanup_ns={d}\tcrash_attempts={d}\tcrash_fallbacks={d}\tcrash_planned_columns={d}\tcrash_structural_columns={d}\tcrash_basis_nonzeros={d}\tcrash_condition={e:.6}\trebuild_calls={d}\trebuild_ns={d}\tinvert_calls={d}\tinvert_ns={d}\tftran_calls={d}\tftran_ns={d}\tbtran_calls={d}\tbtran_ns={d}\tupdate_calls={d}\tupdate_ns={d}\tpricing_calls={d}\tpricing_ns={d}\tdegenerate_pivots={d}\tanti_cycling_activations={d}\tbound_flips={d}\n",
        .{
            path,
            simplex_stats.phase_one_iterations,
            simplex_stats.dual_phase_one_iterations,
            simplex_stats.dual_phase_one_fallbacks,
            simplex_stats.dual_repair_iterations,
            simplex_stats.phase_two_iterations,
            simplex_stats.phase_one_ns,
            simplex_stats.dual_phase_one_ns,
            simplex_stats.dual_repair_ns,
            simplex_stats.phase_two_ns,
            simplex_stats.cleanup_ns,
            simplex_stats.crash_attempts,
            simplex_stats.crash_fallbacks,
            simplex_stats.crash_planned_columns,
            simplex_stats.crash_structural_columns,
            simplex_stats.crash_basis_nonzeros,
            simplex_stats.crash_condition_estimate,
            simplex_stats.rebuild_calls,
            simplex_stats.rebuild_ns,
            factor_stats.factorizations,
            factor_stats.invert_ns,
            factor_stats.ftran_calls,
            factor_stats.ftran_ns,
            factor_stats.btran_calls,
            factor_stats.btran_ns,
            factor_stats.eta_updates + factor_stats.ft_updates,
            factor_stats.update_ns,
            simplex_stats.pricing_calls,
            simplex_stats.pricing_ns,
            engine.numerical.degenerate_pivot_count,
            engine.numerical.anti_cycling_activations,
            simplex_stats.bound_flips,
        },
    );
    defer allocator.free(stats_line);
    try std.Io.File.stdout().writeStreamingAll(io_context, stats_line);

    const iteration_counters = engine.iteration_counters;
    const canonical_iterations_line = try std.fmt.allocPrint(
        allocator,
        "stats\t{s}\tattempted_iterations={d}\tcommitted_pivots={d}\tbound_moves={d}\tshifted_dual_pivots={d}\tdual_phase1_pivots={d}\tdual_phase2_pivots={d}\tprimal_phase1_pivots={d}\tprimal_phase2_pivots={d}\tdual_repair_pivots={d}\tcleanup_pivots={d}\tclassified_pivots={d}\n",
        .{
            path,
            iteration_counters.attempted_iterations,
            iteration_counters.committed_pivots,
            iteration_counters.bound_moves,
            iteration_counters.shifted_dual_pivots,
            iteration_counters.dual_phase_one_pivots,
            iteration_counters.dual_phase_two_pivots,
            iteration_counters.primal_phase_one_pivots,
            iteration_counters.primal_phase_two_pivots,
            iteration_counters.dual_repair_pivots,
            iteration_counters.cleanup_pivots,
            iteration_counters.classifiedPivots(),
        },
    );
    defer allocator.free(canonical_iterations_line);
    try std.Io.File.stdout().writeStreamingAll(io_context, canonical_iterations_line);

    const dual_snapshot_stats_line = try std.fmt.allocPrint(
        allocator,
        "stats\t{s}\tdual_phase1_snapshot_retries={d}\n",
        .{ path, simplex_stats.dual_phase_one_snapshot_retries },
    );
    defer allocator.free(dual_snapshot_stats_line);
    try std.Io.File.stdout().writeStreamingAll(io_context, dual_snapshot_stats_line);

    const infeasibility_certificate_line = try std.fmt.allocPrint(
        allocator,
        "stats\t{s}\tinfeasibility_ray_valid={}\tinfeasibility_certificate_gap={e:.17}\tfailure={s}\tinfinite_row_mass={e:.17}\tinfinite_column_mass={e:.17}\n",
        .{
            path,
            engine.infeasibility_ray_valid,
            engine.infeasibility_certificate_gap,
            @tagName(engine.infeasibility_certificate_failure),
            engine.infeasibility_certificate_infinite_row_mass,
            engine.infeasibility_certificate_infinite_column_mass,
        },
    );
    defer allocator.free(infeasibility_certificate_line);
    try std.Io.File.stdout().writeStreamingAll(io_context, infeasibility_certificate_line);

    const bound_flip_stats_line = try std.fmt.allocPrint(
        allocator,
        "stats\t{s}\tbound_flip_batches={d}\tbound_flip_ftran_savings={d}\n",
        .{ path, simplex_stats.bound_flip_batches, simplex_stats.bound_flip_ftran_savings },
    );
    defer allocator.free(bound_flip_stats_line);
    try std.Io.File.stdout().writeStreamingAll(io_context, bound_flip_stats_line);

    const dual_reprice_stats_line = try std.fmt.allocPrint(
        allocator,
        "stats\t{s}\tdual_reduced_cost_updates={d}\tdual_exact_reprices={d}\tdual_dse_updates={d}\tdual_devex_updates={d}\tdual_dse_invalid_fallbacks={d}\tdual_dse_budget_fallbacks={d}\n",
        .{ path, simplex_stats.dual_reduced_cost_updates, simplex_stats.dual_exact_reprices, simplex_stats.dual_dse_updates, simplex_stats.dual_devex_updates, simplex_stats.dual_dse_invalid_fallbacks, simplex_stats.dual_dse_budget_fallbacks },
    );
    defer allocator.free(dual_reprice_stats_line);
    try std.Io.File.stdout().writeStreamingAll(io_context, dual_reprice_stats_line);

    const shifted_dual_stats_line = try std.fmt.allocPrint(
        allocator,
        "stats\t{s}\tshifted_dual_exit={s}\tshifted_dual_failure_site={s}\tshifted_dual_iterations={d}\tshifted_cleanup_iterations={d}\n",
        .{ path, @tagName(engine.shifted_dual_exit), @tagName(engine.shifted_dual_failure_site), simplex_stats.shifted_dual_iterations, simplex_stats.shifted_cleanup_iterations },
    );
    defer allocator.free(shifted_dual_stats_line);
    try std.Io.File.stdout().writeStreamingAll(io_context, shifted_dual_stats_line);

    const devex_stats_line = try std.fmt.allocPrint(
        allocator,
        "stats\t{s}\tdevex_frameworks={d}\tdevex_framework_updates={d}\tdevex_bad_weights={d}\n",
        .{ path, simplex_stats.devex_frameworks, simplex_stats.devex_framework_updates, simplex_stats.devex_bad_weights },
    );
    defer allocator.free(devex_stats_line);
    try std.Io.File.stdout().writeStreamingAll(io_context, devex_stats_line);

    const partial_pricing_stats_line = try std.fmt.allocPrint(
        allocator,
        "stats\t{s}\tpartial_searches={d}\tpartial_scanned_entries={d}\tpartial_full_scans={d}\n",
        .{ path, engine.pricing.partial_searches, engine.pricing.partial_scanned_entries, engine.pricing.partial_full_scans },
    );
    defer allocator.free(partial_pricing_stats_line);
    try std.Io.File.stdout().writeStreamingAll(io_context, partial_pricing_stats_line);

    const degeneracy_stats_line = try std.fmt.allocPrint(
        allocator,
        "stats\t{s}\tdegenerate_classified={d}\tdegenerate_bound_tie={d}\tdegenerate_ratio_tie={d}\tdegenerate_zero_step={d}\tdegenerate_phase1_stall={d}\tdegenerate_repeated_basis={d}\tdegenerate_small_pivot={d}\tdegenerate_bound_flip={d}\tperturbation_activations={d}\tperturbation_expirations={d}\tperturbation_cleanups={d}\tcold_restart_solves={d}\tcold_restart_phase_one={d}\ttaboo_records={d}\ttaboo_retries={d}\texact_reprices={d}\tmax_reduced_cost_drift={e:.6}\n",
        .{
            path,
            simplex_stats.classifiedDegeneratePivots(),
            simplex_stats.degeneracy_bound_ties,
            simplex_stats.degeneracy_ratio_ties,
            simplex_stats.degeneracy_zero_primal_steps,
            simplex_stats.degeneracy_phase_one_objective_stalls,
            simplex_stats.degeneracy_repeated_bases,
            simplex_stats.degeneracy_small_pivot_retries,
            simplex_stats.degeneracy_bound_flips,
            simplex_stats.perturbation_activations,
            simplex_stats.perturbation_expirations,
            simplex_stats.perturbation_cleanups,
            simplex_stats.cold_restart_solves,
            simplex_stats.cold_restart_phase_one,
            simplex_stats.taboo_records,
            simplex_stats.taboo_retries,
            engine.exact_reprices,
            engine.maximum_reduced_cost_drift,
        },
    );
    defer allocator.free(degeneracy_stats_line);
    try std.Io.File.stdout().writeStreamingAll(io_context, degeneracy_stats_line);

    const kernel_stats_line = try std.fmt.allocPrint(
        allocator,
        "stats\t{s}\tupdate_limit_reinversions={d}\tupdate_growth_reinversions={d}\tsolve_residual_reinversions={d}\tsmall_pivot_reinversions={d}\tmax_ft_chain={d}\tmax_update_growth={e:.6}\tmax_ftran_residual={e:.6}\tftran_rhs_density={e:.6}\tbtran_rhs_density={e:.6}\taq_density={e:.6}\tep_density={e:.6}\tpricing_density={e:.6}\tdense_ftran={d}\thyper_ftran={d}\tdense_btran={d}\thyper_btran={d}\tdense_pricing={d}\thyper_pricing={d}\trow_pricing={d}\tcolumn_pricing={d}\trequested_bytes={d}\tpeak_rss_kb={d}\n",
        .{
            path,
            factor_stats.update_limit_reinversions,
            factor_stats.update_growth_reinversions,
            factor_stats.solve_residual_reinversions,
            factor_stats.small_pivot_reinversions,
            factor_stats.maximum_update_count,
            factor_stats.maximum_update_growth,
            engine.numerical.max_ftran_relative_residual,
            ftran_rhs_density,
            btran_rhs_density,
            aq_density,
            ep_density,
            pricing_density,
            factor_stats.dense_ftran_dispatches,
            factor_stats.hyper_ftran_dispatches,
            factor_stats.dense_btran_dispatches,
            factor_stats.hyper_btran_dispatches,
            simplex_stats.dense_pricing_dispatches,
            simplex_stats.hyper_pricing_dispatches,
            simplex_stats.row_pricing_dispatches,
            simplex_stats.column_pricing_dispatches,
            requested_bytes,
            peak_rss_kb,
        },
    );
    defer allocator.free(kernel_stats_line);
    try std.Io.File.stdout().writeStreamingAll(io_context, kernel_stats_line);

    const first_update_failure = if (simplex_stats.first_update_failure_kind) |kind| @tagName(kind) else "none";
    const rebuild_stats_line = try std.fmt.allocPrint(
        allocator,
        "stats\t{s}\trebuild_phase_one_setup={d}\trebuild_solve_residual={d}\trebuild_small_pivot={d}\trebuild_update_limit={d}\trebuild_update_growth={d}\trebuild_direction_refinement={d}\trebuild_fresh_mode={d}\trebuild_update_rejected={d}\trebuild_cleanup={d}\trebuild_edge_weight_reset={d}\trebuild_numerical_policy={d}\tupdate_dimension_failures={d}\tupdate_unsupported_failures={d}\tupdate_singular_failures={d}\tupdate_numerical_failures={d}\tupdate_out_of_memory_failures={d}\tdense_update_failures={d}\tsparse_update_failures={d}\tfirst_update_failure={s}\tfirst_update_failure_iteration={d}\tfirst_update_failure_entering={d}\tfirst_update_failure_leaving_row={d}\n",
        .{
            path,
            simplex_stats.rebuild_phase_one_setup,
            simplex_stats.rebuild_solve_residual,
            simplex_stats.rebuild_small_pivot,
            simplex_stats.rebuild_update_limit,
            simplex_stats.rebuild_update_growth,
            simplex_stats.rebuild_direction_refinement,
            simplex_stats.rebuild_fresh_mode,
            simplex_stats.rebuild_update_rejected,
            simplex_stats.rebuild_cleanup,
            simplex_stats.rebuild_edge_weight_reset,
            simplex_stats.rebuild_numerical_policy,
            factor_stats.update_dimension_failures,
            factor_stats.update_unsupported_failures,
            factor_stats.update_singular_failures,
            factor_stats.update_numerical_failures,
            factor_stats.update_out_of_memory_failures,
            factor_stats.dense_update_failures,
            factor_stats.sparse_update_failures,
            first_update_failure,
            simplex_stats.first_update_failure_iteration,
            simplex_stats.first_update_failure_entering,
            simplex_stats.first_update_failure_leaving_row,
        },
    );
    defer allocator.free(rebuild_stats_line);
    try std.Io.File.stdout().writeStreamingAll(io_context, rebuild_stats_line);
}
