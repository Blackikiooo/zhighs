//! Mehrotra predictor-corrector primal-dual interior-point engine.
//!
//! Solves the standard form min ĉᵀx̂, Âx̂ = b̂, 0 ≤ x̂ ≤ û produced by
//! `standard_form.zig` and maps the converged iterate back to model
//! coordinates. Finite upper bounds stay native: each bounded column carries
//! a slack w = û − x̂ with dual t, so the normal matrix keeps the original
//! row count. Each iteration assembles the lower triangle of A·Θ·Aᵀ
//! with Θ = (X⁻¹S + W⁻¹T)⁻¹ in reusable sparse storage (or the reference
//! dense kernel), factors it with Cholesky (degenerate pivots are lifted to
//! an escalating per-pivot floor), and combines an affine predictor with a
//! centrality-corrected second step under a fraction-to-boundary line search.

const std = @import("std");
const Allocator = std.mem.Allocator;
const lp = @import("lp");

const NormalFactor = @import("normal_factor.zig").NormalFactor;
const graph_mod = @import("graph.zig");
const options_mod = @import("options.zig");
const row_ordering = @import("row_ordering.zig");
const standard_form = @import("standard_form.zig");

pub const Options = options_mod.Options;
pub const StandardForm = standard_form.StandardForm;

const inf = std.math.inf(f64);

/// Owning solve result in original model coordinates.
pub const Result = struct {
    allocator: Allocator,
    status: lp.SolveStatus,
    /// Structural primal values, length `num_cols` (empty on failure paths).
    primal: []f64,
    /// Row multipliers in the original objective convention, length `num_rows`.
    dual: []f64,
    /// Column reduced costs in the original objective convention, length `num_cols`.
    reduced_cost: []f64,
    /// Objective including sense and offset; meaningful when status is optimal.
    objective_value: f64,
    iterations: usize,

    pub fn deinit(self: *Result) void {
        self.allocator.free(self.primal);
        self.allocator.free(self.dual);
        self.allocator.free(self.reduced_cost);
    }

    /// Borrowed view compatible with the simplex-layer solution contract.
    pub fn solutionView(self: *const Result) lp.SolutionView {
        return .{
            .status = self.status,
            .primal = self.primal,
            .dual = self.dual,
            .reduced_cost = self.reduced_cost,
            .unbounded_ray = &.{},
            .infeasibility_ray = &.{},
            .objective_value = self.objective_value,
            .iterations = self.iterations,
        };
    }
};

fn emptyResult(allocator: Allocator, status: lp.SolveStatus) Result {
    return .{
        .allocator = allocator,
        .status = status,
        .primal = &.{},
        .dual = &.{},
        .reduced_cost = &.{},
        .objective_value = 0.0,
        .iterations = 0,
    };
}

/// Solve `problem` with the given options. Invalid input is reported through
/// the returned status. Invalid options and a missing time-limit clock are
/// errors, as are allocation failures and storage-size overflow.
pub fn solve(allocator: Allocator, problem: lp.ProblemView, options: Options) !Result {
    return solveWithStats(allocator, problem, options, null, null);
}

/// Solve while accumulating per-phase timings into `stats` when `io` is set.
pub fn solveWithStats(
    allocator: Allocator,
    problem: lp.ProblemView,
    options: Options,
    io: ?std.Io,
    stats: ?*Stats,
) !Result {
    if (!options.valid()) return error.InvalidOptions;
    if (options.time_limit_ns != null and io == null) return error.ClockRequired;
    const started = clockNs(io);
    var initial_options = options;
    // Large, wide sparse models benefit from equilibration even when their
    // raw coefficient range is moderate: it reduces Cholesky growth and
    // usually lowers the factorization cost on the first attempt.
    const scale_first = !options.equilibrate and options.retry_equilibrated_on_failure and
        (coefficientRange(problem) >= 1e7 or
            (problem.num_rows >= 500 and problem.num_cols >= 1000));
    if (!options.equilibrate and options.retry_equilibrated_on_failure) {
        initial_options.max_iterations = @min(initial_options.max_iterations, 40);
        if (options.time_limit_ns) |limit| initial_options.time_limit_ns = @max(limit / 2, 1);
    }
    if (scale_first) initial_options.equilibrate = true;
    var first = try solveAttempt(allocator, problem, initial_options, io, stats);
    if (first.status == .optimal or first.status == .interrupted or first.status == .work_limit or
        first.status == .infeasible or first.status == .unbounded or
        options.max_iterations == 0 or options.equilibrate or !options.retry_equilibrated_on_failure)
        return first;

    var retry_options = options;
    retry_options.equilibrate = !scale_first;
    retry_options.retry_equilibrated_on_failure = false;
    retry_options.max_iterations = @max(retry_options.max_iterations, 300);
    if (retry_options.equilibrate) {
        retry_options.refinement_steps = @max(retry_options.refinement_steps, 1);
        retry_options.primal_correction_steps = @max(retry_options.primal_correction_steps, 20);
        retry_options.tolerances.max_regularization = @max(retry_options.tolerances.max_regularization, 1e-4);
    }
    if (options.time_limit_ns) |limit| {
        const elapsed: u64 = @intCast(@max(clockNs(io) - started, 0));
        if (elapsed >= limit) return first;
        retry_options.time_limit_ns = limit - elapsed;
    }
    first.deinit();
    if (stats) |value| value.* = .{};
    return solveAttempt(allocator, problem, retry_options, io, stats);
}

fn coefficientRange(problem: lp.ProblemView) f64 {
    var minimum = inf;
    var maximum: f64 = 0.0;
    for (problem.matrix.values) |value| {
        const magnitude = @abs(value);
        if (magnitude == 0.0 or !std.math.isFinite(magnitude)) continue;
        minimum = @min(minimum, magnitude);
        maximum = @max(maximum, magnitude);
    }
    return if (minimum == inf) 1.0 else maximum / minimum;
}

fn solveAttempt(
    allocator: Allocator,
    problem: lp.ProblemView,
    options: Options,
    io: ?std.Io,
    stats: ?*Stats,
) !Result {
    const init_started = clockNs(io);
    var sf = standard_form.build(allocator, problem) catch |err| switch (err) {
        error.InvalidModel => return emptyResult(allocator, .invalid_model),
        else => return err,
    };
    defer sf.deinit();
    if (options.equilibrate) sf.equilibrate();
    const t_sf = clockNs(io);

    // One shared adjacency of A·Aᵀ feeds the ordering heuristic, the fill
    // estimates, and the symbolic factorization.
    var graph_storage: ?graph_mod.Graph = null;
    defer if (graph_storage) |*g| g.deinit();
    var permutation: ?[]usize = null;
    defer if (permutation) |map| allocator.free(map);
    if (!options.dense_normal) {
        graph_storage = try graph_mod.build(allocator, sf.num_rows, sf.col_starts, sf.row_indices);
        if (io != null) {
            var max_col: usize = 0;
            var dense_cols: usize = 0;
            var total_len: usize = 0;
            for (0..sf.num_cols) |j| {
                const len = sf.col_starts[j + 1] - sf.col_starts[j];
                max_col = @max(max_col, len);
                total_len += len;
                if (len >= 32) dense_cols += 1;
            }
            // Top column lengths, descending, to guide the dense-column cut.
            var lens = try allocator.alloc(usize, sf.num_cols);
            defer allocator.free(lens);
            for (0..sf.num_cols) |j| lens[j] = sf.col_starts[j + 1] - sf.col_starts[j];
            std.mem.sort(usize, lens, {}, std.sort.desc(usize));
            const g0 = &graph_storage.?;
            std.debug.print("  graph_build={} us nnz={} density={d:.3} max_col={} dense32={} avg={d:.1} m={} toplens=", .{ @divTrunc(clockNs(io) - t_sf, 1000), g0.starts[g0.n], g0.density(), max_col, dense_cols, @as(f64, @floatFromInt(total_len)) / @as(f64, @floatFromInt(sf.num_cols)), sf.num_rows });
            for (lens[0..@min(24, lens.len)]) |l| std.debug.print("{} ", .{l});
            std.debug.print("\n", .{});
        }
        const g = &graph_storage.?;
        if (sf.num_rows >= 2) {
            const selected: options_mod.RowOrdering = switch (options.row_ordering) {
                .automatic => if (g.density() >= 0.3 and
                    (sf.num_rows >= 300 or sf.num_cols >= 3 * sf.num_rows))
                    .reverse_cuthill_mckee
                else
                    .amd,
                else => options.row_ordering,
            };
            const candidate: ?[]usize = switch (selected) {
                .automatic, .natural => null,
                .reverse_cuthill_mckee => try row_ordering.reverseCuthillMcKee(allocator, g),
                .minimum_degree => try row_ordering.minimumDegree(allocator, g),
                .amd => try row_ordering.approximateMinimumDegree(allocator, g),
            };
            if (candidate) |map| {
                // Rank the candidate against the natural order by exact
                // symbolic fill; the loser is discarded.
                const t_ord = clockNs(io);
                const reordered_fill = try graph_mod.fillCount(allocator, g, map);
                const natural_fill = try graph_mod.fillCount(allocator, g, null);
                if (io != null) std.debug.print("  ord={} fill2={} us\n", .{ @divTrunc(t_ord - t_sf, 1000), @divTrunc(clockNs(io) - t_ord, 1000) });
                if (reordered_fill < natural_fill) {
                    try sf.permuteRows(map);
                    permutation = map;
                } else {
                    allocator.free(map);
                }
            }
        }
    }
    const t_graph = clockNs(io);
    var order = try graph_mod.Order.init(allocator, sf.num_rows, permutation);
    defer order.deinit();

    var engine = try Engine.init(allocator, &sf, options, if (graph_storage) |*g| g else null, &order);
    defer engine.deinit();
    if (io != null) std.debug.print("init sf={} graph+ord={} engine={} us\n", .{ @divTrunc(t_sf - init_started, 1000), @divTrunc(t_graph - t_sf, 1000), @divTrunc(clockNs(io) - t_graph, 1000) });
    engine.stats_io = io;
    engine.stats = stats;
    if (stats) |value| addElapsed(io, init_started, &value.init_ns);
    const outcome = engine.run();

    return mapResult(allocator, problem, &sf, &engine, outcome);
}

const Outcome = struct {
    status: lp.SolveStatus,
    iterations: usize,
};

/// Slice `len` elements off `block` at `cursor` and advance the cursor.
fn carve(block: []f64, cursor: *usize, len: usize) []f64 {
    const slice = block[cursor.*..][0..len];
    cursor.* += len;
    return slice;
}

/// Per-phase wall-clock accumulation for one solve, populated only when
/// `stats_io` is set on the engine (benchmarking hook, zero cost otherwise).
pub const Stats = struct {
    normal_entries: usize = 0,
    factor_products: usize = 0,
    factor_indirect_products: usize = 0,
    /// Standard-form build, ordering, and symbolic factorization.
    init_ns: u64 = 0,
    assembly_ns: u64 = 0,
    factor_ns: u64 = 0,
    direction_ns: u64 = 0,
    residual_ns: u64 = 0,
    total_ns: u64 = 0,
    final_primal_residual: f64 = 0,
    final_dual_residual: f64 = 0,
    final_gap: f64 = 0,
    final_mu: f64 = 0,
    final_max_x: f64 = 0,
    /// 1=nonfinite mu, 2=mu growth, 3=primal escape, 4=initial factor,
    /// 5=iteration factor, 6/8=primal infeasible, 7/9=unbounded,
    /// 10=dual escape.
    failure_code: u8 = 0,
};

fn clockNs(io: ?std.Io) i96 {
    const clock_io = io orelse return 0;
    return std.Io.Clock.awake.now(clock_io).nanoseconds;
}

fn addElapsed(io: ?std.Io, start: i96, counter: *u64) void {
    if (io == null) return;
    counter.* += @intCast(@max(clockNs(io) - start, 0));
}

const Engine = struct {
    allocator: Allocator,
    sf: *const StandardForm,
    options: Options,
    stats_io: ?std.Io = null,
    stats: ?*Stats = null,

    x: []f64,
    /// Lower-bound dual slack (reduced-cost part), one per column.
    s: []f64,
    /// Upper slack w = û − x; inf for columns without a finite upper bound.
    w: []f64,
    /// Upper-bound dual slack; exactly zero for unbounded columns.
    t: []f64,
    y: []f64,
    dx: []f64,
    ds: []f64,
    dt: []f64,
    dy: []f64,
    dx_aff: []f64,
    ds_aff: []f64,
    dt_aff: []f64,
    dy_aff: []f64,
    rp: []f64,
    rd: []f64,
    rxs: []f64,
    rwt: []f64,
    work_n: []f64,
    rhs: []f64,
    refinement: []f64,
    normal: []f64,
    normal_positions: []usize,
    chol: NormalFactor,
    /// Backing allocation for every f64 slice above (single block).
    block: []f64,

    /// Complementarity pairs: n lower bounds plus one per finite upper bound.
    num_pairs: usize,
    mu0: f64,

    fn init(allocator: Allocator, sf: *const StandardForm, options: Options, g: ?*const graph_mod.Graph, order: *const graph_mod.Order) !Engine {
        const n = sf.num_cols;
        const m = sf.num_rows;
        var chol = try NormalFactor.init(allocator, sf, options.dense_normal, g, order);
        errdefer chol.deinit();
        // All iterate and workspace vectors plus the normal-matrix values are
        // carved out of one contiguous block: a single allocation instead of
        // two dozen, with the iterate vectors sharing cache lines.
        var block_len = try std.math.mul(usize, 14, n);
        block_len = try std.math.add(usize, block_len, try std.math.mul(usize, 6, m));
        block_len = try std.math.add(usize, block_len, chol.storageLen());
        const block = try allocator.alloc(f64, block_len);
        errdefer allocator.free(block);
        var cursor: usize = 0;
        const x = carve(block, &cursor, n);
        const s = carve(block, &cursor, n);
        const w = carve(block, &cursor, n);
        const t = carve(block, &cursor, n);
        const dx = carve(block, &cursor, n);
        const ds = carve(block, &cursor, n);
        const dt = carve(block, &cursor, n);
        const dx_aff = carve(block, &cursor, n);
        const ds_aff = carve(block, &cursor, n);
        const dt_aff = carve(block, &cursor, n);
        const rd = carve(block, &cursor, n);
        const rxs = carve(block, &cursor, n);
        const rwt = carve(block, &cursor, n);
        const work_n = carve(block, &cursor, n);
        const y = carve(block, &cursor, m);
        const dy = carve(block, &cursor, m);
        const dy_aff = carve(block, &cursor, m);
        const rp = carve(block, &cursor, m);
        const rhs = carve(block, &cursor, m);
        const refinement = carve(block, &cursor, m);
        const normal = carve(block, &cursor, chol.storageLen());
        std.debug.assert(cursor == block.len);
        var position_count: usize = 0;
        for (0..sf.num_cols) |j| {
            const begin = sf.col_starts[j];
            const end = sf.col_starts[j + 1];
            for (begin..end) |p| for (begin..end) |q| {
                if (sf.row_indices[q] <= sf.row_indices[p]) position_count += 1;
            };
        }
        const normal_positions = try allocator.alloc(usize, position_count);
        errdefer allocator.free(normal_positions);
        var position: usize = 0;
        for (0..sf.num_cols) |j| {
            const begin = sf.col_starts[j];
            const end = sf.col_starts[j + 1];
            for (begin..end) |p| for (begin..end) |q| {
                const row_p = sf.row_indices[p];
                const row_q = sf.row_indices[q];
                if (row_q <= row_p) {
                    normal_positions[position] = chol.index(row_p, row_q);
                    position += 1;
                }
            };
        }

        var num_pairs: usize = n;
        for (sf.upper) |u| {
            if (u != inf) num_pairs += 1;
        }

        return .{
            .allocator = allocator,
            .sf = sf,
            .options = options,
            .x = x,
            .s = s,
            .w = w,
            .t = t,
            .y = y,
            .dx = dx,
            .ds = ds,
            .dt = dt,
            .dy = dy,
            .dx_aff = dx_aff,
            .ds_aff = ds_aff,
            .dt_aff = dt_aff,
            .dy_aff = dy_aff,
            .rp = rp,
            .rd = rd,
            .rxs = rxs,
            .rwt = rwt,
            .work_n = work_n,
            .rhs = rhs,
            .refinement = refinement,
            .normal = normal,
            .normal_positions = normal_positions,
            .chol = chol,
            .block = block,
            .num_pairs = num_pairs,
            .mu0 = 0.0,
        };
    }

    fn deinit(self: *Engine) void {
        const a = self.allocator;
        a.free(self.block);
        a.free(self.normal_positions);
        self.chol.deinit();
    }

    /// Diagonal weight applied to each column during normal-matrix assembly.
    const NormalWeights = union(enum) {
        /// Θ = (X⁻¹S + W⁻¹T)⁻¹ from the current iterate.
        theta,
        /// All-ones diagonal: assembles A·Aᵀ for the starting-point projection.
        unit,
        /// Explicit per-column weights (e.g. a 0/1 mask).
        mask: []const f64,
    };

    /// Lower-triangle assembly of A·W·Aᵀ with diagonal W given by `weights`.
    fn assembleNormal(self: *Engine, weights: NormalWeights) void {
        const started = clockNs(self.stats_io);
        defer if (self.stats) |stats| addElapsed(self.stats_io, started, &stats.assembly_ns);
        const sf = self.sf;
        @memset(self.normal, 0.0);
        var position: usize = 0;
        for (0..sf.num_cols) |j| {
            const weight = switch (weights) {
                .theta => self.theta(j),
                .unit => 1.0,
                .mask => |mask| mask[j],
            };
            const begin = sf.col_starts[j];
            const end = sf.col_starts[j + 1];
            for (begin..end) |p| {
                const row_p = sf.row_indices[p];
                const vp = weight * sf.values[p];
                for (begin..end) |q| {
                    const row_q = sf.row_indices[q];
                    if (row_q <= row_p) {
                        self.normal[self.normal_positions[position]] += vp * sf.values[q];
                        position += 1;
                    }
                }
            }
        }
    }

    /// Diagonal of the normal-equation scaling Θ = (X⁻¹S + W⁻¹T)⁻¹.
    inline fn theta(self: *const Engine, j: usize) f64 {
        var diagonal = self.s[j] / self.x[j];
        if (self.sf.upper[j] != inf) diagonal += self.t[j] / self.w[j];
        return 1.0 / diagonal;
    }

    /// Factor the normal matrix, escalating the per-pivot floor on failure.
    /// The factorization lifts any pivot below the floor to the floor — a
    /// diagonal shift applied only where degenerate — so `self.normal` keeps
    /// the true matrix and iterative refinement sees unshifted residuals.
    /// Returns false when no floor recovers positive definiteness.
    fn factorWithRegularization(self: *Engine) bool {
        const started = clockNs(self.stats_io);
        defer if (self.stats) |stats| addElapsed(self.stats_io, started, &stats.factor_ns);
        const tol = self.options.tolerances;
        var floor = tol.pivot_regularization;
        var attempts: usize = 0;
        while (attempts < 16) : (attempts += 1) {
            if (self.chol.factor(self.normal, floor)) return true else |_| {}
            if (floor > tol.max_regularization) return false;
            floor *= tol.regularization_growth;
        }
        return false;
    }

    /// y := A·u for a dense standard-form vector u.
    fn matvec(self: *const Engine, u: []const f64, out: []f64) void {
        const sf = self.sf;
        @memset(out, 0.0);
        for (0..sf.num_cols) |j| {
            const uj = u[j];
            if (uj == 0.0) continue;
            for (sf.col_starts[j]..sf.col_starts[j + 1]) |k| {
                out[sf.row_indices[k]] += sf.values[k] * uj;
            }
        }
    }

    /// out := Aᵀ·y for a dense row-space vector y.
    fn matvecTranspose(self: *const Engine, y: []const f64, out: []f64) void {
        const sf = self.sf;
        for (0..sf.num_cols) |j| {
            var acc: f64 = 0.0;
            for (sf.col_starts[j]..sf.col_starts[j + 1]) |k| {
                acc += sf.values[k] * y[sf.row_indices[k]];
            }
            out[j] = acc;
        }
    }

    /// Solve one normal-equations right-hand side and expand the full
    /// (Δx, Δy, Δs, Δt) direction from the stored residuals and the
    /// complementarity targets `rxs` (lower bounds) and `rwt` (upper bounds).
    fn direction(self: *Engine, rxs: []const f64, rwt: []const f64, dx: []f64, dy: []f64, ds: []f64, dt: []f64) void {
        const started = clockNs(self.stats_io);
        defer if (self.stats) |stats| addElapsed(self.stats_io, started, &stats.direction_ns);
        const sf = self.sf;
        for (0..sf.num_cols) |j| {
            var target = self.rd[j] + rxs[j] / self.x[j];
            if (sf.upper[j] != inf) target -= rwt[j] / self.w[j];
            self.work_n[j] = self.theta(j) * target;
        }
        self.matvec(self.work_n, self.rhs);
        for (0..sf.num_rows) |i| self.rhs[i] += self.rp[i];
        @memcpy(dy, self.rhs);
        self.chol.solve(dy);
        for (0..self.options.refinement_steps) |_| {
            self.chol.multiply(self.normal, dy, self.refinement);
            for (self.refinement, self.rhs) |*residual, original| residual.* = original - residual.*;
            self.chol.solve(self.refinement);
            for (dy, self.refinement) |*value, correction| value.* += correction;
        }
        self.matvecTranspose(dy, dx);
        for (0..sf.num_cols) |j| {
            var target = dx[j] - self.rd[j] - rxs[j] / self.x[j];
            if (sf.upper[j] != inf) target += rwt[j] / self.w[j];
            dx[j] = self.theta(j) * target;
        }
        for (0..sf.num_cols) |j| {
            ds[j] = -(rxs[j] + self.s[j] * dx[j]) / self.x[j];
            dt[j] = if (sf.upper[j] != inf) (self.t[j] * dx[j] - rwt[j]) / self.w[j] else 0.0;
        }
    }

    /// Mehrotra starting-point heuristic: project onto the constraints with
    /// Θ = I, split the dual slack into lower/upper parts, then lift all
    /// complementarity pairs to strict positivity while keeping clearance
    /// from finite upper bounds.
    fn startingPoint(self: *Engine) void {
        const sf = self.sf;
        // x̃ = Aᵀ(AAᵀ)⁻¹b.
        @memcpy(self.rhs, sf.rhs);
        self.chol.solve(self.rhs);
        self.matvecTranspose(self.rhs, self.x);
        // ỹ = (AAᵀ)⁻¹Ac, c − Aᵀỹ split into s (lower) and t (upper) parts.
        self.matvec(sf.cost, self.rhs);
        self.chol.solve(self.rhs);
        @memcpy(self.y, self.rhs);
        self.matvecTranspose(self.y, self.s);
        for (0..sf.num_cols) |j| {
            const slack = sf.cost[j] - self.s[j];
            if (sf.upper[j] == inf) {
                self.s[j] = slack;
                self.t[j] = 0.0;
            } else {
                self.s[j] = @max(slack, 0.0);
                self.t[j] = @max(-slack, 0.0);
            }
        }

        const margin = self.options.start_margin;
        var min_unbounded_x: f64 = 0.0;
        for (self.x, sf.upper) |x, u| {
            if (u == inf) min_unbounded_x = @min(min_unbounded_x, x);
        }
        const lift_x = @max(-1.5 * min_unbounded_x, 0.0) + margin;
        var min_dual: f64 = 0.0;
        for (0..sf.num_cols) |j| {
            min_dual = @min(min_dual, self.s[j]);
            if (sf.upper[j] != inf) min_dual = @min(min_dual, self.t[j]);
        }
        const lift_s = @max(-1.5 * min_dual, 0.0) + margin;
        for (0..sf.num_cols) |j| {
            if (sf.upper[j] == inf) self.x[j] += lift_x;
        }
        self.clampBelowUpper(margin);
        for (0..sf.num_cols) |j| {
            self.s[j] += lift_s;
            if (sf.upper[j] != inf) self.t[j] += lift_s;
        }
        var dot_val: f64 = 0.0;
        var sum_primal: f64 = 0.0;
        var sum_dual: f64 = 0.0;
        for (0..sf.num_cols) |j| {
            sum_primal += self.x[j];
            sum_dual += self.s[j];
            dot_val += self.x[j] * self.s[j];
            if (sf.upper[j] != inf) {
                sum_primal += self.w[j];
                sum_dual += self.t[j];
                dot_val += self.w[j] * self.t[j];
            }
        }
        const extra_x = 0.5 * dot_val / sum_dual;
        const extra_s = 0.5 * dot_val / sum_primal;
        for (0..sf.num_cols) |j| {
            if (sf.upper[j] == inf) self.x[j] += extra_x;
        }
        self.clampBelowUpper(margin);
        for (0..sf.num_cols) |j| {
            self.s[j] += extra_s;
            if (sf.upper[j] != inf) self.t[j] += extra_s;
        }
    }

    /// Enforce 0 < x < û on bounded columns and refresh w = û − x.
    fn clampBelowUpper(self: *Engine, margin: f64) void {
        const sf = self.sf;
        for (0..sf.num_cols) |j| {
            const u = sf.upper[j];
            if (u == inf) {
                self.w[j] = inf;
                continue;
            }
            const clearance = @min(0.1 * u, @max(margin, 1e-8 * u));
            if (self.x[j] >= u or self.x[j] <= 0.0)
                self.x[j] = clearance
            else
                self.x[j] = std.math.clamp(self.x[j], clearance, u - clearance);
            self.w[j] = u - self.x[j];
        }
    }

    /// w = û − x on bounded columns; called after every primal update.
    fn refreshUpperSlack(self: *Engine) void {
        const sf = self.sf;
        for (0..sf.num_cols) |j| {
            if (sf.upper[j] != inf) self.w[j] = sf.upper[j] - self.x[j];
        }
    }

    fn computeResiduals(self: *Engine) void {
        const started = clockNs(self.stats_io);
        defer if (self.stats) |stats| addElapsed(self.stats_io, started, &stats.residual_ns);
        const sf = self.sf;
        self.matvec(self.x, self.rp);
        for (0..sf.num_rows) |i| self.rp[i] = sf.rhs[i] - self.rp[i];
        self.matvecTranspose(self.y, self.rd);
        for (0..sf.num_cols) |j| self.rd[j] = sf.cost[j] - self.rd[j] - self.s[j] + self.t[j];
    }

    /// Complementarity measure over lower and upper bound pairs.
    fn complementarityMu(self: *const Engine) f64 {
        var total = dot(self.x, self.s);
        for (0..self.sf.num_cols) |j| {
            if (self.sf.upper[j] != inf) total += self.w[j] * self.t[j];
        }
        return total / @as(f64, @floatFromInt(@max(self.num_pairs, 1)));
    }

    /// Dual objective bᵀy − ûᵀt in standard-form (minimization) convention.
    fn dualObjective(self: *const Engine) f64 {
        var value = dot(self.sf.rhs, self.y);
        for (0..self.sf.num_cols) |j| {
            if (self.sf.upper[j] != inf) value -= self.sf.upper[j] * self.t[j];
        }
        return value;
    }

    fn run(self: *Engine) Outcome {
        const started = clockNs(self.stats_io);
        defer if (self.stats) |stats| addElapsed(self.stats_io, started, &stats.total_ns);
        const tol = self.options.tolerances;
        const sf = self.sf;
        const n = sf.num_cols;

        if (self.limitStatus(started, 0)) |status| return .{ .status = status, .iterations = 0 };
        if (self.stats) |stats| stats.normal_entries = self.normal.len;
        if (self.stats) |stats| stats.factor_products = self.chol.productCount();
        if (self.stats) |stats| stats.factor_indirect_products = self.chol.indirectProductCount();
        @memset(self.y, 0);
        if (n == 0) return .{
            .status = if (maxAbs(sf.rhs) <= tol.primal_feasibility) .optimal else .infeasible,
            .iterations = 0,
        };
        if (sf.num_rows == 0) {
            // min cᵀx over a box: negative costs run to their upper bound.
            for (0..n) |j| {
                const u = sf.upper[j];
                if (sf.cost[j] < 0.0) {
                    if (u == inf) return .{ .status = .unbounded, .iterations = 0 };
                    self.x[j] = u;
                    self.s[j] = 0.0;
                    self.t[j] = -sf.cost[j];
                } else {
                    self.x[j] = 0.0;
                    self.s[j] = sf.cost[j];
                    self.t[j] = 0.0;
                }
            }
            return .{ .status = .optimal, .iterations = 0 };
        }

        self.assembleNormal(.unit);
        if (!self.factorWithRegularization()) {
            if (self.stats) |stats| stats.failure_code = 4;
            return .{ .status = .numerical_failure, .iterations = 0 };
        }
        self.startingPoint();
        self.mu0 = @max(self.complementarityMu(), 1.0);
        self.computeResiduals();

        var iteration: usize = 0;
        var best_pres: f64 = inf;
        var best_dres: f64 = inf;
        // Residuals captured before the most recent step. Divergence
        // certificates classify against these so one escaping step cannot
        // mask which side was still feasible.
        var prev_pres = scaledPrimalResidual(self.rp, sf.rhs, sf.row_scale);
        var prev_dres = scaledMaxResidual(self.rd, sf.cost);
        var stall_primal: usize = 0;
        var stall_dual: usize = 0;
        while (iteration < self.options.max_iterations) : (iteration += 1) {
            if (self.options.iteration_callback) |callback| {
                if (callback(iteration, self.options.callback_context)) return .{ .status = .interrupted, .iterations = iteration };
            }
            if (self.limitStatus(started, iteration)) |status| return .{ .status = status, .iterations = iteration };
            const mu = self.complementarityMu();
            const pres = scaledPrimalResidual(self.rp, sf.rhs, sf.row_scale);
            const dres = scaledMaxResidual(self.rd, sf.cost);
            const primal_obj = dot(sf.cost, self.x);
            const dual_obj = self.dualObjective();
            const gap = @abs(primal_obj - dual_obj) / (1.0 + @abs(primal_obj));
            if (self.stats) |stats| {
                stats.final_primal_residual = pres;
                stats.final_dual_residual = dres;
                stats.final_gap = gap;
                stats.final_mu = mu;
                stats.final_max_x = maxAbs(self.x);
            }

            if (pres <= tol.primal_feasibility and dres <= tol.dual_feasibility and gap <= tol.relative_gap) {
                return .{ .status = .optimal, .iterations = iteration };
            }
            // No-progress certificates: one residual halved at most over ten
            // iterations while the other side has already converged.
            if (pres < 0.5 * best_pres) {
                best_pres = pres;
                stall_primal = 0;
            } else {
                stall_primal += 1;
            }
            if (dres < 0.5 * best_dres) {
                best_dres = dres;
                stall_dual = 0;
            } else {
                stall_dual += 1;
            }

            if (stall_primal >= 10 and
                pres > tol.primal_feasibility * tol.heuristic_residual_relaxation and
                dres <= tol.dual_feasibility)
            {
                if (self.stats) |stats| stats.failure_code = 8;
                return .{ .status = .infeasible, .iterations = iteration };
            }
            if (stall_dual >= 10 and
                dres > tol.dual_feasibility * tol.heuristic_residual_relaxation and
                pres <= tol.primal_feasibility)
            {
                if (self.stats) |stats| stats.failure_code = 9;
                return .{ .status = .unbounded, .iterations = iteration };
            }
            // A global norm can hide a small-rhs row behind an unrelated large
            // rhs. Before taking another central-path step, project such a
            // near-optimal candidate back onto Ax=b with the current Θ metric.
            if (globalScaledResidual(self.rp, sf.rhs) <= tol.primal_feasibility and
                globalScaledResidual(self.rd, sf.cost) <= tol.dual_feasibility and
                gap <= tol.relative_gap)
            {
                if (!self.correctPrimalFeasibility())
                    return .{ .status = .numerical_failure, .iterations = iteration };
                const corrected_pres = scaledPrimalResidual(self.rp, sf.rhs, sf.row_scale);
                const corrected_dres = scaledMaxResidual(self.rd, sf.cost);
                const corrected_primal_obj = dot(sf.cost, self.x);
                const corrected_dual_obj = self.dualObjective();
                const corrected_gap = @abs(corrected_primal_obj - corrected_dual_obj) /
                    (1.0 + @abs(corrected_primal_obj));
                if (corrected_pres <= tol.primal_feasibility and
                    corrected_dres <= tol.dual_feasibility and corrected_gap <= tol.relative_gap)
                    return .{ .status = .optimal, .iterations = iteration };
                continue;
            }
            // Either side escaping to infinity certifies the opposite problem:
            // x -> inf while Ax=b keeps holding exhibits a primal ray (dual
            // infeasible); y/s/t -> inf while the dual residual keeps holding
            // exhibits a dual ray (primal infeasible).
            const x_escaped = maxAbs(self.x) > tol.escape_magnitude;
            const dual_escaped = maxAbs(self.y) > tol.escape_magnitude or
                maxAbs(self.s) > tol.escape_magnitude or
                maxAbs(self.t) > tol.escape_magnitude;
            if (!std.math.isFinite(mu) or mu > tol.divergence_growth * self.mu0 or x_escaped or dual_escaped) {
                if (self.stats) |stats| stats.failure_code = if (!std.math.isFinite(mu)) 1 else if (mu > tol.divergence_growth * self.mu0) 2 else if (x_escaped) 3 else 10;
                // Classify against the pre-step residuals (min convention): a
                // dual that stayed feasible while the primal residual stalled
                // certifies primal infeasibility; a primal that stayed
                // feasible while the dual residual stalled certifies
                // unboundedness. Both stuck or both feasible is a genuine
                // numerical failure.
                const primal_stuck = prev_pres > tol.primal_feasibility * tol.heuristic_residual_relaxation;
                const dual_stuck = prev_dres > tol.dual_feasibility * tol.heuristic_residual_relaxation;
                if (primal_stuck and !dual_stuck) {
                    if (self.stats) |stats| stats.failure_code = 6;
                    return .{ .status = .infeasible, .iterations = iteration };
                }
                if (dual_stuck and !primal_stuck) {
                    if (self.stats) |stats| stats.failure_code = 7;
                    return .{ .status = .unbounded, .iterations = iteration };
                }
                return .{ .status = .numerical_failure, .iterations = iteration };
            }

            self.assembleNormal(.theta);
            if (!self.factorWithRegularization()) {
                if (self.stats) |stats| stats.failure_code = 5;
                return .{ .status = .numerical_failure, .iterations = iteration };
            }

            // Affine predictor: complementarity targets x∘s and w∘t.
            for (0..n) |j| {
                self.rxs[j] = self.x[j] * self.s[j];
                self.rwt[j] = if (sf.upper[j] != inf) self.w[j] * self.t[j] else 0.0;
            }
            self.direction(self.rxs, self.rwt, self.dx_aff, self.dy_aff, self.ds_aff, self.dt_aff);
            const alpha_p_aff = fractionToBoundaryUpper(self.x, self.w, sf.upper, self.dx_aff, 1.0);
            const alpha_d_aff = @min(fractionToBoundary(self.s, self.ds_aff, 1.0), fractionToBoundary(self.t, self.dt_aff, 1.0));
            var mu_aff: f64 = 0.0;
            for (0..n) |j| {
                mu_aff += (self.x[j] + alpha_p_aff * self.dx_aff[j]) * (self.s[j] + alpha_d_aff * self.ds_aff[j]);
                if (sf.upper[j] != inf) {
                    mu_aff += (self.w[j] - alpha_p_aff * self.dx_aff[j]) * (self.t[j] + alpha_d_aff * self.dt_aff[j]);
                }
            }
            mu_aff /= @floatFromInt(@max(self.num_pairs, 1));
            const ratio = if (mu > 0.0) mu_aff / mu else 0.0;
            const sigma = @min(ratio * ratio * ratio, 1.0);

            // Corrector with Mehrotra centrality and second-order terms; the
            // upper-pair second-order term uses Δw = −Δx.
            for (0..n) |j| {
                self.rxs[j] = self.x[j] * self.s[j] - sigma * mu + self.dx_aff[j] * self.ds_aff[j];
                self.rwt[j] = if (sf.upper[j] != inf)
                    self.w[j] * self.t[j] - sigma * mu - self.dx_aff[j] * self.dt_aff[j]
                else
                    0.0;
            }
            self.direction(self.rxs, self.rwt, self.dx, self.dy, self.ds, self.dt);

            // Gondzio centrality corrections: reuse the factorization to pull
            // outlying trial complementarity products into [0.1 μ, 10 μ].
            // The corrections accumulate on the Mehrotra right-hand side, and
            // a corrected direction is kept only when it enlarges the step.
            var alpha_p_trial = fractionToBoundaryUpper(self.x, self.w, sf.upper, self.dx, 1.0);
            var alpha_d_trial = @min(fractionToBoundary(self.s, self.ds, 1.0), fractionToBoundary(self.t, self.dt, 1.0));
            var best_step = @min(alpha_p_trial, alpha_d_trial);
            for (0..self.options.centrality_corrections) |_| {
                const lo = 0.1 * mu;
                const hi = 10.0 * mu;
                var outlier = false;
                for (0..n) |j| {
                    const lower_prod = (self.x[j] + alpha_p_trial * self.dx[j]) * (self.s[j] + alpha_d_trial * self.ds[j]);
                    var adjust: f64 = 0.0;
                    if (lower_prod < lo) {
                        adjust = lo - lower_prod;
                    } else if (lower_prod > hi) {
                        adjust = hi - lower_prod;
                    }
                    if (adjust != 0.0) outlier = true;
                    self.rxs[j] -= adjust;
                    if (sf.upper[j] != inf) {
                        const upper_prod = (self.w[j] - alpha_p_trial * self.dx[j]) * (self.t[j] + alpha_d_trial * self.dt[j]);
                        var adjust_upper: f64 = 0.0;
                        if (upper_prod < lo) {
                            adjust_upper = lo - upper_prod;
                        } else if (upper_prod > hi) {
                            adjust_upper = hi - upper_prod;
                        }
                        if (adjust_upper != 0.0) outlier = true;
                        self.rwt[j] -= adjust_upper;
                    }
                }
                if (!outlier) break;
                self.direction(self.rxs, self.rwt, self.dx_aff, self.dy_aff, self.ds_aff, self.dt_aff);
                const candidate_p = fractionToBoundaryUpper(self.x, self.w, sf.upper, self.dx_aff, 1.0);
                const candidate_d = @min(fractionToBoundary(self.s, self.ds_aff, 1.0), fractionToBoundary(self.t, self.dt_aff, 1.0));
                if (@min(candidate_p, candidate_d) <= best_step * 1.0001) break;
                @memcpy(self.dx, self.dx_aff);
                @memcpy(self.dy, self.dy_aff);
                @memcpy(self.ds, self.ds_aff);
                @memcpy(self.dt, self.dt_aff);
                alpha_p_trial = candidate_p;
                alpha_d_trial = candidate_d;
                best_step = @min(candidate_p, candidate_d);
            }

            const step = self.options.step_fraction;
            const alpha_p = @min(1.0, step * fractionToBoundaryUpper(self.x, self.w, sf.upper, self.dx, 1.0));
            const alpha_d = @min(1.0, step * @min(fractionToBoundary(self.s, self.ds, 1.0), fractionToBoundary(self.t, self.dt, 1.0)));
            for (0..n) |j| {
                self.x[j] += alpha_p * self.dx[j];
                self.s[j] += alpha_d * self.ds[j];
                self.t[j] += alpha_d * self.dt[j];
            }
            for (0..sf.num_rows) |i| self.y[i] += alpha_d * self.dy[i];
            self.refreshUpperSlack();
            self.computeResiduals();
            prev_pres = pres;
            prev_dres = dres;
        }
        // Report optimality if the final iterate already meets tolerances.
        const pres = scaledPrimalResidual(self.rp, sf.rhs, sf.row_scale);
        const dres = scaledMaxResidual(self.rd, sf.cost);
        const primal_obj = dot(sf.cost, self.x);
        const dual_obj = self.dualObjective();
        const gap = @abs(primal_obj - dual_obj) / (1.0 + @abs(primal_obj));
        if (pres <= tol.primal_feasibility and dres <= tol.dual_feasibility and gap <= tol.relative_gap) {
            return .{ .status = .optimal, .iterations = iteration };
        }
        return .{ .status = .iteration_limit, .iterations = iteration };
    }

    /// Correct Ax=b while retaining 0 <= x <= û. After using explicit slack
    /// columns, apply a minimum-Θ inverse projection. Repeated
    /// fraction-to-boundary steps handle a correction that would otherwise
    /// cross a bound. The projection metric is assembled and factored once
    /// and reused for every correction step: Θ drifts only with x, and a
    /// stale metric still defines a valid projection.
    fn correctPrimalFeasibility(self: *Engine) bool {
        if (self.options.equilibrate) return self.correctPrimalFeasibilityActive();
        var factored = false;
        var unit_metric = false;
        for (0..self.options.primal_correction_steps) |_| {
            if (scaledPrimalResidual(self.rp, self.sf.rhs, self.sf.row_scale) <= self.options.tolerances.primal_feasibility) return true;
            // Prefer zero-cost singleton columns (ordinary inequality slacks):
            // they repair one row exactly without perturbing any other row or
            // the objective. `refinement` marks rows already assigned.
            @memset(self.dx, 0.0);
            @memset(self.refinement, 0.0);
            for (0..self.sf.num_cols) |j| {
                const begin = self.sf.col_starts[j];
                const end = self.sf.col_starts[j + 1];
                if (end - begin != 1 or self.sf.cost[j] != 0.0) continue;
                const row = self.sf.row_indices[begin];
                if (self.refinement[row] != 0.0) continue;
                const correction = self.rp[row] / self.sf.values[begin];
                if (self.x[j] + correction < 0.0) continue;
                if (self.sf.upper[j] != inf and self.x[j] + correction > self.sf.upper[j]) continue;
                self.dx[j] = correction;
                self.refinement[row] = 1.0;
            }
            for (self.x, self.dx) |*x, correction| x.* += correction;
            self.refreshUpperSlack();
            self.computeResiduals();
            if (scaledPrimalResidual(self.rp, self.sf.rhs, self.sf.row_scale) <= self.options.tolerances.primal_feasibility) return true;

            if (!factored) {
                self.assembleNormal(.theta);
                if (!self.factorWithRegularization()) {
                    // Late Θ may be too ill-conditioned to factor even though a
                    // small feasibility projection remains well-defined.
                    self.assembleNormal(.unit);
                    if (!self.factorWithRegularization()) return false;
                    unit_metric = true;
                }
                factored = true;
            }
            @memcpy(self.rhs, self.rp);
            self.chol.solve(self.rhs);
            self.matvecTranspose(self.rhs, self.work_n);
            for (self.dx, self.work_n, 0..) |*correction, aty, j| {
                correction.* = if (unit_metric) aty else self.theta(j) * aty;
            }
            const alpha = @min(1.0, self.options.step_fraction * fractionToBoundaryUpper(self.x, self.w, self.sf.upper, self.dx, 1.0));
            for (self.x, self.dx) |*x, correction| x.* += alpha * correction;
            self.refreshUpperSlack();
            self.computeResiduals();
        }
        return true;
    }

    fn correctPrimalFeasibilityActive(self: *Engine) bool {
        @memset(self.ds, 1.0);
        // The masked unit metric changes only when a column is excluded, so
        // the factorization is rebuilt solely after an exclusion.
        var factored = false;
        for (0..self.options.primal_correction_steps) |_| {
            if (scaledPrimalResidual(self.rp, self.sf.rhs, self.sf.row_scale) <= self.options.tolerances.primal_feasibility) return true;
            if (!factored) {
                self.assembleNormal(.{ .mask = self.ds });
                if (!self.factorWithRegularization()) return false;
                factored = true;
            }
            @memcpy(self.rhs, self.rp);
            self.chol.solve(self.rhs);
            self.matvecTranspose(self.rhs, self.work_n);
            var excluded = false;
            for (self.dx, self.work_n, self.x, self.ds, 0..) |*correction, aty, x, *allowed, j| {
                correction.* = allowed.* * aty;
                if (allowed.* != 0.0 and (x + correction.* < 0.0 or
                    (self.sf.upper[j] != inf and x + correction.* > self.sf.upper[j])))
                {
                    allowed.* = 0.0;
                    excluded = true;
                }
            }
            if (excluded) {
                factored = false;
                continue;
            }
            for (self.x, self.dx) |*x, correction| x.* += correction;
            self.refreshUpperSlack();
            self.computeResiduals();
        }
        return true;
    }

    fn limitStatus(self: *const Engine, started: i96, iteration: usize) ?lp.SolveStatus {
        if (self.options.interrupt_flag) |flag| if (flag.load(.acquire)) return .interrupted;
        if (self.options.time_limit_ns) |limit| if (clockNs(self.stats_io) - started >= limit) return .time_limit;
        if (self.options.work_limit) |limit| if (iteration >= limit) return .work_limit;
        return null;
    }
};

fn mapResult(
    allocator: Allocator,
    problem: lp.ProblemView,
    sf: *const StandardForm,
    engine: *Engine,
    outcome: Outcome,
) !Result {
    var result = emptyResult(allocator, outcome.status);
    result.iterations = outcome.iterations;
    if (outcome.status != .optimal) return result;

    errdefer result.deinit();
    // A split free variable may accumulate a large common component in its
    // positive and negative parts. Removing it preserves A*x and c*x exactly
    // in real arithmetic and avoids catastrophic cancellation when mapping
    // back to the original free variable.
    for (sf.col_records) |record| {
        if (record.kind != .free_split) continue;
        const negative = record.negative.?;
        const common = @min(engine.x[record.primary], engine.x[negative]);
        engine.x[record.primary] -= common;
        engine.x[negative] -= common;
    }
    result.primal = try allocator.alloc(f64, problem.num_cols);
    result.dual = try allocator.alloc(f64, problem.num_rows);
    result.reduced_cost = try allocator.alloc(f64, problem.num_cols);

    const sense: f64 = @floatFromInt(@intFromEnum(problem.objective_sense));

    // Primal: reassemble each column from its standard-form parts.
    for (0..problem.num_cols) |j| {
        const record = sf.col_records[j];
        result.primal[j] = switch (record.kind) {
            .fixed => record.base,
            .shifted => record.base + sf.col_scale[record.primary] * engine.x[record.primary],
            .flipped => record.base - sf.col_scale[record.primary] * engine.x[record.primary],
            .free_split => sf.col_scale[record.primary] * engine.x[record.primary] -
                sf.col_scale[record.negative.?] * engine.x[record.negative.?],
        };
    }

    // Row multipliers in min convention, then flipped into the model sense.
    const dual_min = try allocator.alloc(f64, problem.num_rows);
    defer allocator.free(dual_min);
    for (0..problem.num_rows) |i| {
        dual_min[i] = if (sf.row_to_std[i]) |r| sf.row_scale[r] * engine.y[r] else 0.0;
        result.dual[i] = sense * dual_min[i];
    }

    // Reduced costs follow the simplex-layer convention c − Aᵀy, computed in
    // min convention from the row multipliers and flipped into the model sense.
    for (0..problem.num_cols) |j| {
        result.reduced_cost[j] = sense * sf.columnReducedCostMin(problem, j, dual_min);
    }

    result.objective_value = sf.originalObjective(problem, dot(sf.cost, engine.x));
    return result;
}

fn maxAbs(v: []const f64) f64 {
    var acc: f64 = 0.0;
    for (v) |value| {
        if (!std.math.isFinite(value)) return inf;
        acc = @max(acc, @abs(value));
    }
    return acc;
}

fn scaledMaxResidual(residual: []const f64, reference: []const f64) f64 {
    var result: f64 = 0.0;
    for (residual, reference) |value, scale| {
        if (!std.math.isFinite(value) or !std.math.isFinite(scale)) return inf;
        result = @max(result, @abs(value) / (1.0 + @abs(scale)));
    }
    return result;
}

fn scaledPrimalResidual(residual: []const f64, rhs: []const f64, row_scale: []const f64) f64 {
    var result: f64 = 0.0;
    for (residual, rhs, row_scale) |value, reference, scale| {
        if (!std.math.isFinite(value) or !std.math.isFinite(reference) or scale <= 0.0) return inf;
        // rhs and residual are both multiplied by row_scale. This is the
        // relative residual in the unscaled standard-form row.
        result = @max(result, @abs(value) / (scale + @abs(reference)));
    }
    return result;
}

fn globalScaledResidual(residual: []const f64, reference: []const f64) f64 {
    return maxAbs(residual) / (1.0 + maxAbs(reference));
}

fn minValue(v: []const f64) f64 {
    var acc: f64 = 0.0;
    for (v) |value| acc = @min(acc, value);
    return acc;
}

fn sum(v: []const f64) f64 {
    var acc: f64 = 0.0;
    for (v) |value| acc += value;
    return acc;
}

fn dot(a: []const f64, b: []const f64) f64 {
    var acc: f64 = 0.0;
    for (a, b) |x, y| acc += x * y;
    return acc;
}

/// Largest α ≤ cap such that v + α·d stays non-negative.
fn fractionToBoundary(v: []const f64, d: []const f64, cap: f64) f64 {
    var alpha = cap;
    for (v, d) |value, step| {
        if (step < 0.0) alpha = @min(alpha, -value / step);
    }
    return alpha;
}

/// Largest α ≤ cap keeping 0 ≤ x + α·d ≤ û: lower boundary from x, upper
/// boundary from the stored slack w = û − x on columns with finite û.
fn fractionToBoundaryUpper(x: []const f64, w: []const f64, upper: []const f64, d: []const f64, cap: f64) f64 {
    var alpha = cap;
    for (x, d, 0..) |value, step, j| {
        if (step < 0.0) {
            alpha = @min(alpha, -value / step);
        } else if (step > 0.0 and upper[j] != inf) {
            alpha = @min(alpha, w[j] / step);
        }
    }
    return alpha;
}