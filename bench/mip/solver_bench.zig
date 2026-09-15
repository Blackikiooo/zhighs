//! Repeatable end-to-end MIP benchmark and profiling entry point.
//!
//! Parses once, then reports aggregate branch-and-bound timings and structural
//! counters as TSV. Build with ReleaseFast and run the installed binary under
//! samply, perf, VTune, or another sampling profiler.

const std = @import("std");
const zhighs = @import("zhighs");

fn nowNs(io: std.Io) i128 {
    return std.Io.Clock.now(.awake, io).nanoseconds;
}

/// Windows in-process sampling profiler: a helper thread suspends the main
/// thread, reads its instruction pointer, and resumes it. After the run the
/// raw addresses are resolved against the PDB via std.debug.SelfInfo.
const SampleProfiler = struct {
    const windows = std.os.windows;
    const max_samples = 1 << 16;
    const sample_interval_100ns: i64 = -10000; // 1ms (rounded up to timer resolution)

    extern "kernel32" fn GetCurrentProcess() callconv(.c) windows.HANDLE;
    extern "kernel32" fn GetCurrentThread() callconv(.c) windows.HANDLE;
    extern "kernel32" fn DuplicateHandle(
        src_proc: windows.HANDLE,
        src: windows.HANDLE,
        dst_proc: windows.HANDLE,
        out: *windows.HANDLE,
        access: windows.DWORD,
        inherit: windows.BOOL,
        options: windows.DWORD,
    ) callconv(.c) windows.BOOL;
    extern "kernel32" fn SuspendThread(handle: windows.HANDLE) callconv(.c) windows.DWORD;
    extern "kernel32" fn ResumeThread(handle: windows.HANDLE) callconv(.c) windows.DWORD;
    extern "kernel32" fn GetThreadContext(handle: windows.HANDLE, ctx: *windows.CONTEXT) callconv(.c) windows.BOOL;
    extern "ntdll" fn NtDelayExecution(alertable: windows.BOOLEAN, delay: *const windows.LARGE_INTEGER) callconv(.c) windows.NTSTATUS;
    extern "winmm" fn timeBeginPeriod(period: windows.DWORD) callconv(.c) windows.DWORD;
    extern "winmm" fn timeEndPeriod(period: windows.DWORD) callconv(.c) windows.DWORD;

    stop: std.atomic.Value(bool) = .init(false),
    target: windows.HANDLE,
    samples: [max_samples]usize = undefined,
    count: std.atomic.Value(usize) = .init(0),

    fn currentThreadHandle() windows.HANDLE {
        var out: windows.HANDLE = undefined;
        if (DuplicateHandle(GetCurrentProcess(), GetCurrentThread(), GetCurrentProcess(), &out, 0, .FALSE, 2) == .FALSE)
            @panic("DuplicateHandle failed");
        return out;
    }

    fn run(self: *SampleProfiler) void {
        _ = timeBeginPeriod(1);
        defer _ = timeEndPeriod(1);
        while (!self.stop.load(.acquire)) {
            if (SuspendThread(self.target) != std.math.maxInt(windows.DWORD)) {
                var ctx: windows.CONTEXT = undefined;
                ctx.ContextFlags = 0x100003; // CONTEXT_AMD64 | CONTROL | INTEGER
                if (GetThreadContext(self.target, &ctx) != .FALSE) {
                    const rip: usize = @intCast(ctx.Rip);
                    const index = self.count.fetchAdd(1, .monotonic);
                    if (index < max_samples) self.samples[index] = rip;
                }
                _ = ResumeThread(self.target);
            }
            const interval = sample_interval_100ns;
            _ = NtDelayExecution(.FALSE, &interval);
        }
    }
};

fn printSampleReport(
    io_context: std.Io,
    allocator: std.mem.Allocator,
    profiler: *SampleProfiler,
    label: []const u8,
) !void {
    const total = @min(profiler.count.load(.acquire), SampleProfiler.max_samples);
    if (total == 0) return;
    var buckets = std.AutoHashMap(usize, usize).init(allocator);
    defer buckets.deinit();
    for (profiler.samples[0..total]) |rip| {
        const key = rip & ~@as(usize, 0x3F);
        const entry = try buckets.getOrPut(key);
        if (entry.found_existing) entry.value_ptr.* += 1 else entry.value_ptr.* = 1;
    }
    var info = std.debug.SelfInfo.init;
    defer info.deinit(io_context);
    var text_arena = std.heap.ArenaAllocator.init(allocator);
    defer text_arena.deinit();
    var symbol_totals = std.StringHashMap(usize).init(allocator);
    defer symbol_totals.deinit();
    var unresolved: usize = 0;
    var iterator = buckets.iterator();
    while (iterator.next()) |entry| {
        var symbols: std.ArrayList(std.debug.Symbol) = .empty;
        defer symbols.deinit(allocator);
        info.getSymbols(io_context, allocator, text_arena.allocator(), entry.key_ptr.*, false, &symbols) catch {
            unresolved += entry.value_ptr.*;
            continue;
        };
        if (symbols.items.len == 0) {
            unresolved += entry.value_ptr.*;
            continue;
        }
        const name = symbols.items[symbols.items.len - 1].name orelse "?";
        const symbol_entry = try symbol_totals.getOrPut(name);
        if (symbol_entry.found_existing) symbol_entry.value_ptr.* += entry.value_ptr.* else symbol_entry.value_ptr.* = entry.value_ptr.*;
    }
    const SortEntry = struct { name: []const u8, count: usize };
    var ordered = std.ArrayList(SortEntry).empty;
    defer ordered.deinit(allocator);
    var symbol_iterator = symbol_totals.iterator();
    while (symbol_iterator.next()) |entry| try ordered.append(allocator, .{ .name = entry.key_ptr.*, .count = entry.value_ptr.* });
    std.mem.sort(SortEntry, ordered.items, {}, comptime struct {
        fn lessThan(_: void, a: SortEntry, b: SortEntry) bool {
            return a.count > b.count;
        }
    }.lessThan);
    const header = try std.fmt.allocPrint(allocator, "zhighs-sample\t{s}\tsamples={d}\tunresolved={d}\tsymbols={d}\n", .{ label, total, unresolved, ordered.items.len });
    defer allocator.free(header);
    try std.Io.File.stdout().writeStreamingAll(io_context, header);
    for (ordered.items, 0..) |item, rank| {
        if (rank >= 40) break;
        const line = try std.fmt.allocPrint(allocator, "zhighs-sample\t{s}\trank={d}\tpct={d:.2}\tcount={d}\t{s}\n", .{ label, rank, 100.0 * @as(f64, @floatFromInt(item.count)) / @as(f64, @floatFromInt(total)), item.count, item.name });
        defer allocator.free(line);
        try std.Io.File.stdout().writeStreamingAll(io_context, line);
    }
}

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;
    const io_context = init.io;
    var args = try init.minimal.args.iterateAllocator(allocator);
    defer args.deinit();
    _ = args.next();
    const path = args.next() orelse return error.MissingModelPath;
    const repetitions = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 7;
    const profile = args.next() orelse "default";
    const sample_runs = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 0;
    if (repetitions == 0 or args.next() != null) return error.InvalidArguments;
    var mip_control = if (std.mem.eql(u8, profile, "default"))
        zhighs.mip.SolveControl{}
    else if (std.mem.eql(u8, profile, "no-cuts"))
        zhighs.mip.SolveControl{ .cut_separation = false }
    else if (std.mem.eql(u8, profile, "no-cliques"))
        zhighs.mip.SolveControl{ .clique_propagation = false }
    else if (std.mem.eql(u8, profile, "no-conflicts"))
        zhighs.mip.SolveControl{ .conflict_learning = false }
    else if (std.mem.eql(u8, profile, "no-heuristics"))
        zhighs.mip.SolveControl{ .rounding_heuristic = false }
    else if (std.mem.eql(u8, profile, "rounding-only"))
        zhighs.mip.SolveControl{}
    else if (std.mem.eql(u8, profile, "branching"))
        zhighs.mip.SolveControl{ .propagation = false }
    else if (std.mem.eql(u8, profile, "branching-no-probes"))
        zhighs.mip.SolveControl{ .propagation = false, .reliability_branching = false }
    else if (std.mem.eql(u8, profile, "most-fractional"))
        zhighs.mip.SolveControl{ .propagation = false, .reliability_branching = false }
    else
        return error.InvalidProfile;
    var most_fractional = zhighs.plugins_builtin.branching.MostFractional{};
    var rounding = zhighs.plugins_builtin.heuristics.Rounding{};
    var rounding_plugins = [_]zhighs.plugin.PrimalHeuristic{
        zhighs.plugin.bindHeuristic(zhighs.plugins_builtin.heuristics.Rounding, &rounding),
    };
    if (std.mem.eql(u8, profile, "rounding-only")) mip_control.components.heuristics = &rounding_plugins;
    if (std.mem.eql(u8, profile, "most-fractional"))
        mip_control.components.branching_rule = zhighs.plugin.bindBranching(zhighs.plugins_builtin.branching.MostFractional, &most_fractional);

    var model = try zhighs.io.readFile(io_context, allocator, path, .{ .keep_names = false });
    defer model.deinit();
    const integrality = try allocator.alloc(zhighs.mip.Integrality, model.col_type.len);
    defer allocator.free(integrality);
    for (model.col_type, integrality) |source, *target| target.* = switch (source) {
        .continuous => .continuous,
        .binary => .binary,
        .integer => .integer,
        .semi_continuous => .semi_continuous,
        .semi_integer => .semi_integer,
    };
    const problem = zhighs.mip.ProblemView{
        .relaxation = .{
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
        },
        .integrality = integrality,
    };

    var lp_elapsed_ns: u128 = 0;
    for (0..repetitions) |_| {
        const started = nowNs(io_context);
        var result = zhighs.solver.solveLp(allocator, problem.relaxation, .{});
        lp_elapsed_ns += @intCast(nowNs(io_context) - started);
        if (result.status != .optimal) return error.UnexpectedRelaxationStatus;
        result.deinit();
    }

    // Warm reoptimization microbenchmark: one cold root solve, then N bound
    // changes each followed by a warm reoptimization through the same session.
    {
        var session = zhighs.solver.LpSolveSession.init(allocator);
        defer session.deinit();
        var bounds = try allocator.dupe(f64, problem.relaxation.col_upper);
        defer allocator.free(bounds);
        var warm_problem = problem.relaxation;
        warm_problem.col_upper = bounds;
        _ = session.solve(warm_problem, .{ .collect_statistics = true }, .{ .structure = 0, .matrix_values = 0, .bounds = 1, .objective = 0 });
        {
            const stats = session.engine.stats;
            const line = try std.fmt.allocPrint(allocator, "zhighs-root\t{s}\troot_iterations={d}\tphase1_ns={d}\tdual_phase1_ns={d}\trepair_ns={d}\tphase2_ns={d}\tcleanup_ns={d}\trebuild_ns={d}\trebuilds={d}\tpricing_ns={d}\tpricing_calls={d}\tphase1_iter={d}\tdual_phase1_iter={d}\tphase2_iter={d}\tbound_flips={d}\tdual_dse_updates={d}\tdual_devex_updates={d}\tdense_price={d}\thyper_price={d}\n", .{ path, session.engine.iterations, stats.phase_one_ns, stats.dual_phase_one_ns, stats.dual_repair_ns, stats.phase_two_ns, stats.cleanup_ns, stats.rebuild_ns, stats.rebuild_calls, stats.pricing_ns, stats.pricing_calls, stats.phase_one_iterations, stats.dual_phase_one_iterations, stats.phase_two_iterations, stats.bound_flips, stats.dual_dse_updates, stats.dual_devex_updates, stats.dense_pricing_dispatches, stats.hyper_pricing_dispatches });
            defer allocator.free(line);
            try std.Io.File.stdout().writeStreamingAll(io_context, line);
        }
        var warm_elapsed_ns: u128 = 0;
        const warm_rounds = 200;
        for (0..warm_rounds) |round| {
            bounds[round % bounds.len] = if (round % 2 == 0) 2.0 else 1.0;
            const started = nowNs(io_context);
            const status = session.solveRelaxation(warm_problem, .{}, .{ .structure = 0, .matrix_values = 0, .bounds = round + 2, .objective = 0 });
            warm_elapsed_ns += @intCast(nowNs(io_context) - started);
            if (status != .optimal) break;
        }
        const line = try std.fmt.allocPrint(allocator, "zhighs-warm\t{s}\twarm_mean_ns={d}\titerations={d}\tcold={d}\treopt={d}\n", .{ path, warm_elapsed_ns / warm_rounds, session.engine.iterations, session.cold_solves, session.reoptimizations });
        defer allocator.free(line);
        try std.Io.File.stdout().writeStreamingAll(io_context, line);
    }

    var mip_elapsed_ns: u128 = 0;
    var last: zhighs.mip.SolveResult = undefined;
    for (0..repetitions) |run| {
        const started = nowNs(io_context);
        var result = try zhighs.mip.solve(allocator, problem, mip_control);
        mip_elapsed_ns += @intCast(nowNs(io_context) - started);
        if (run + 1 == repetitions) {
            last = result;
        } else {
            result.deinit();
        }
    }
    defer last.deinit();

    if (sample_runs != 0) {
        var profiler = SampleProfiler{ .target = SampleProfiler.currentThreadHandle() };
        const sampler_thread = try std.Thread.spawn(.{}, SampleProfiler.run, .{&profiler});
        for (0..sample_runs) |_| {
            var result = try zhighs.mip.solve(allocator, problem, mip_control);
            result.deinit();
        }
        profiler.stop.store(true, .release);
        sampler_thread.join();
        try printSampleReport(io_context, allocator, &profiler, profile);
    }

    const line = try std.fmt.allocPrint(
        allocator,
        "zhighs-mip\t{s}\tprofile={s}\tstatus={s}\truns={d}\tlp_mean_ns={d}\tmip_mean_ns={d}\tobjective={d:.17}\tbound={d:.17}\tnodes={d}\tlp_iterations={d}\tlp_work={d}\tcold_solves={d}\treoptimizations={d}\tpseudocost_observations={d}\tstrong_probes={d}\tstrong_iterations={d}\tstrong_work={d}\tcut_rounds={d}\tcuts_added={d}\tcuts_expired={d}\tclique_edges={d}\tclique_tightenings={d}\treduced_cost_tightenings={d}\troot_reduced_cost_tightenings={d}\tconflicts_learned={d}\tconflict_checks={d}\tconflict_tightenings={d}\tconflicts_expired={d}\tconflicts_compacted={d}\n",
        .{ path, profile, @tagName(last.status), repetitions, lp_elapsed_ns / repetitions, mip_elapsed_ns / repetitions, last.objective_value, last.objective_bound, last.nodes, last.lp_iterations, last.lp_work, last.lp_cold_solves, last.lp_reoptimizations, last.pseudocost_observations, last.strong_branching_probes, last.strong_branching_iterations, last.strong_branching_work, last.cut_rounds, last.cuts_added, last.cuts_expired, last.clique_edges, last.clique_tightenings, last.reduced_cost_tightenings, last.root_reduced_cost_tightenings, last.conflicts_learned, last.conflict_checks, last.conflict_tightenings, last.conflicts_expired, last.conflicts_compacted },
    );
    defer allocator.free(line);
    try std.Io.File.stdout().writeStreamingAll(io_context, line);
}