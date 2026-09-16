//! Mehrotra predictor-corrector primal-dual interior-point engine.
//!
//! Solves the standard form min ĉᵀx̂, Âx̂ = b̂, x̂ ≥ 0 produced by
//! `standard_form.zig` and maps the converged iterate back to model
//! coordinates. Each iteration assembles the lower triangle of A·Θ·Aᵀ
//! with Θ = X·S⁻¹ in reusable profile storage (or the reference dense kernel),
//! factors it with Cholesky (escalating diagonal
//! regularization on failure), and combines an affine predictor with a
//! centrality-corrected second step under a fraction-to-boundary line search.

const std = @import("std");
const Allocator = std.mem.Allocator;
const lp = @import("lp");

const NormalFactor = @import("normal_factor.zig").NormalFactor;
const options_mod = @import("options.zig");
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
    var sf = standard_form.build(allocator, problem) catch |err| switch (err) {
        error.InvalidModel => return emptyResult(allocator, .invalid_model),
        else => return err,
    };
    defer sf.deinit();

    var engine = try Engine.init(allocator, &sf, options);
    defer engine.deinit();
    engine.stats_io = io;
    engine.stats = stats;
    const outcome = engine.run();

    return mapResult(allocator, problem, &sf, &engine, outcome);
}

const Outcome = struct {
    status: lp.SolveStatus,
    iterations: usize,
};

/// Per-phase wall-clock accumulation for one solve, populated only when
/// `stats_io` is set on the engine (benchmarking hook, zero cost otherwise).
pub const Stats = struct {
    normal_entries: usize = 0,
    assembly_ns: u64 = 0,
    factor_ns: u64 = 0,
    direction_ns: u64 = 0,
    residual_ns: u64 = 0,
    total_ns: u64 = 0,
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
    s: []f64,
    y: []f64,
    dx: []f64,
    ds: []f64,
    dy: []f64,
    dx_aff: []f64,
    ds_aff: []f64,
    dy_aff: []f64,
    rp: []f64,
    rd: []f64,
    rxs: []f64,
    work_n: []f64,
    rhs: []f64,
    normal: []f64,
    chol: NormalFactor,

    norm_b: f64,
    norm_c: f64,
    mu0: f64,

    fn init(allocator: Allocator, sf: *const StandardForm, options: Options) !Engine {
        const n = sf.num_cols;
        const m = sf.num_rows;
        const x = try allocator.alloc(f64, n);
        errdefer allocator.free(x);
        const s = try allocator.alloc(f64, n);
        errdefer allocator.free(s);
        const dx = try allocator.alloc(f64, n);
        errdefer allocator.free(dx);
        const ds = try allocator.alloc(f64, n);
        errdefer allocator.free(ds);
        const dx_aff = try allocator.alloc(f64, n);
        errdefer allocator.free(dx_aff);
        const ds_aff = try allocator.alloc(f64, n);
        errdefer allocator.free(ds_aff);
        const rd = try allocator.alloc(f64, n);
        errdefer allocator.free(rd);
        const rxs = try allocator.alloc(f64, n);
        errdefer allocator.free(rxs);
        const work_n = try allocator.alloc(f64, n);
        errdefer allocator.free(work_n);
        const y = try allocator.alloc(f64, m);
        errdefer allocator.free(y);
        const dy = try allocator.alloc(f64, m);
        errdefer allocator.free(dy);
        const dy_aff = try allocator.alloc(f64, m);
        errdefer allocator.free(dy_aff);
        const rp = try allocator.alloc(f64, m);
        errdefer allocator.free(rp);
        const rhs = try allocator.alloc(f64, m);
        errdefer allocator.free(rhs);
        var chol = try NormalFactor.init(allocator, sf, options.dense_normal);
        errdefer chol.deinit();
        const normal = try allocator.alloc(f64, chol.storageLen());
        errdefer allocator.free(normal);

        return .{
            .allocator = allocator,
            .sf = sf,
            .options = options,
            .x = x,
            .s = s,
            .y = y,
            .dx = dx,
            .ds = ds,
            .dy = dy,
            .dx_aff = dx_aff,
            .ds_aff = ds_aff,
            .dy_aff = dy_aff,
            .rp = rp,
            .rd = rd,
            .rxs = rxs,
            .work_n = work_n,
            .rhs = rhs,
            .normal = normal,
            .chol = chol,
            .norm_b = 1.0 + maxAbs(sf.rhs),
            .norm_c = 1.0 + maxAbs(sf.cost),
            .mu0 = 0.0,
        };
    }

    fn deinit(self: *Engine) void {
        const a = self.allocator;
        a.free(self.x);
        a.free(self.s);
        a.free(self.dx);
        a.free(self.ds);
        a.free(self.dx_aff);
        a.free(self.ds_aff);
        a.free(self.rd);
        a.free(self.rxs);
        a.free(self.work_n);
        a.free(self.y);
        a.free(self.dy);
        a.free(self.dy_aff);
        a.free(self.rp);
        a.free(self.rhs);
        a.free(self.normal);
        self.chol.deinit();
    }

    /// Lower-triangle assembly of A·Θ·Aᵀ with Θ_i = x_i / s_i.
    fn assembleNormal(self: *Engine) void {
        const started = clockNs(self.stats_io);
        defer if (self.stats) |stats| addElapsed(self.stats_io, started, &stats.assembly_ns);
        const sf = self.sf;
        @memset(self.normal, 0.0);
        for (0..sf.num_cols) |j| {
            const weight = self.x[j] / self.s[j];
            if (weight == 0.0) continue;
            const begin = sf.col_starts[j];
            const end = sf.col_starts[j + 1];
            for (begin..end) |p| {
                const row_p = sf.row_indices[p];
                const vp = weight * sf.values[p];
                for (begin..end) |q| {
                    const row_q = sf.row_indices[q];
                    if (row_q <= row_p) self.normal[self.chol.index(row_p, row_q)] += vp * sf.values[q];
                }
            }
        }
    }

    /// Factor the normal matrix, escalating the diagonal regularization shift
    /// on failure. Returns false when no shift recovers positive definiteness.
    fn factorWithRegularization(self: *Engine) bool {
        const started = clockNs(self.stats_io);
        defer if (self.stats) |stats| addElapsed(self.stats_io, started, &stats.factor_ns);
        const tol = self.options.tolerances;
        const m = self.sf.num_rows;
        var applied: f64 = 0.0;
        var target = tol.pivot_regularization;
        var attempts: usize = 0;
        while (attempts < 16) : (attempts += 1) {
            if (self.chol.factor(self.normal)) return true else |_| {}
            if (target > tol.max_regularization) return false;
            // The buffer still holds the unshifted matrix; only the
            // incremental shift is added before each retry.
            for (0..m) |i| self.normal[self.chol.index(i, i)] += target - applied;
            applied = target;
            target *= tol.regularization_growth;
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
    /// (Δx, Δy, Δs) direction from the stored residuals and complementarity
    /// target `rxs`.
    fn direction(self: *Engine, rxs: []const f64, dx: []f64, dy: []f64, ds: []f64) void {
        const started = clockNs(self.stats_io);
        defer if (self.stats) |stats| addElapsed(self.stats_io, started, &stats.direction_ns);
        const sf = self.sf;
        for (0..sf.num_cols) |j| {
            self.work_n[j] = (self.x[j] * self.rd[j] + rxs[j]) / self.s[j];
        }
        self.matvec(self.work_n, self.rhs);
        for (0..sf.num_rows) |i| self.rhs[i] += self.rp[i];
        @memcpy(dy, self.rhs);
        self.chol.solve(dy);
        self.matvecTranspose(dy, dx);
        for (0..sf.num_cols) |j| {
            ds[j] = self.rd[j] - dx[j];
            dx[j] = (self.x[j] * (dx[j] - self.rd[j]) - rxs[j]) / self.s[j];
        }
    }

    /// Mehrotra starting-point heuristic: project onto the constraints with
    /// Θ = I, then lift both non-negative blocks to strict positivity.
    fn startingPoint(self: *Engine) void {
        const sf = self.sf;
        // x̃ = Aᵀ(AAᵀ)⁻¹b.
        @memcpy(self.rhs, sf.rhs);
        self.chol.solve(self.rhs);
        self.matvecTranspose(self.rhs, self.x);
        // ỹ = (AAᵀ)⁻¹Ac, s̃ = c − Aᵀỹ.
        self.matvec(sf.cost, self.rhs);
        self.chol.solve(self.rhs);
        @memcpy(self.y, self.rhs);
        self.matvecTranspose(self.y, self.s);
        for (0..sf.num_cols) |j| self.s[j] = sf.cost[j] - self.s[j];

        const lift_x = @max(-1.5 * minValue(self.x), 0.0) + self.options.start_margin;
        const lift_s = @max(-1.5 * minValue(self.s), 0.0) + self.options.start_margin;
        for (0..sf.num_cols) |j| self.x[j] += lift_x;
        for (0..sf.num_cols) |j| self.s[j] += lift_s;
        const dot_val = dot(self.x, self.s);
        const extra_x = 0.5 * dot_val / sum(self.s);
        const extra_s = 0.5 * dot_val / sum(self.x);
        for (0..sf.num_cols) |j| self.x[j] += extra_x;
        for (0..sf.num_cols) |j| self.s[j] += extra_s;
    }

    fn computeResiduals(self: *Engine) void {
        const started = clockNs(self.stats_io);
        defer if (self.stats) |stats| addElapsed(self.stats_io, started, &stats.residual_ns);
        const sf = self.sf;
        self.matvec(self.x, self.rp);
        for (0..sf.num_rows) |i| self.rp[i] = sf.rhs[i] - self.rp[i];
        self.matvecTranspose(self.y, self.rd);
        for (0..sf.num_cols) |j| self.rd[j] = sf.cost[j] - self.rd[j] - self.s[j];
    }

    fn run(self: *Engine) Outcome {
        const started = clockNs(self.stats_io);
        defer if (self.stats) |stats| addElapsed(self.stats_io, started, &stats.total_ns);
        const tol = self.options.tolerances;
        const sf = self.sf;
        const n = sf.num_cols;

        if (self.limitStatus(started, 0)) |status| return .{ .status = status, .iterations = 0 };
        if (self.stats) |stats| stats.normal_entries = self.normal.len;
        @memset(self.y, 0);
        if (n == 0) return .{
            .status = if (maxAbs(sf.rhs) <= tol.primal_feasibility) .optimal else .infeasible,
            .iterations = 0,
        };
        if (sf.num_rows == 0) {
            @memset(self.x, 0);
            for (sf.cost) |cost| if (cost < 0) return .{ .status = .unbounded, .iterations = 0 };
            return .{ .status = .optimal, .iterations = 0 };
        }

        self.assembleUnitNormal();
        if (!self.factorWithRegularization()) return .{ .status = .numerical_failure, .iterations = 0 };
        self.startingPoint();
        self.mu0 = @max(dot(self.x, self.s) / @as(f64, @floatFromInt(@max(n, 1))), 1.0);
        self.computeResiduals();

        var iteration: usize = 0;
        while (iteration < self.options.max_iterations) : (iteration += 1) {
            if (self.options.iteration_callback) |callback| {
                if (callback(iteration, self.options.callback_context)) return .{ .status = .interrupted, .iterations = iteration };
            }
            if (self.limitStatus(started, iteration)) |status| return .{ .status = status, .iterations = iteration };
            const mu = dot(self.x, self.s) / @as(f64, @floatFromInt(@max(n, 1)));
            const pres = maxAbs(self.rp) / self.norm_b;
            const dres = maxAbs(self.rd) / self.norm_c;
            const primal_obj = dot(sf.cost, self.x);
            const dual_obj = dot(sf.rhs, self.y);
            const gap = @abs(primal_obj - dual_obj) / (1.0 + @abs(primal_obj));

            if (pres <= tol.primal_feasibility and dres <= tol.dual_feasibility and gap <= tol.relative_gap) {
                return .{ .status = .optimal, .iterations = iteration };
            }
            if (!std.math.isFinite(mu) or mu > tol.divergence_growth * self.mu0 or
                maxAbs(self.x) > tol.escape_magnitude)
            {
                // Divergence alone is not a primal or dual infeasibility certificate.
                return .{ .status = .numerical_failure, .iterations = iteration };
            }

            self.assembleNormal();
            if (!self.factorWithRegularization()) return .{ .status = .numerical_failure, .iterations = iteration };

            // Affine predictor: complementarity target r_xs = x∘s.
            for (0..n) |j| self.rxs[j] = self.x[j] * self.s[j];
            self.direction(self.rxs, self.dx_aff, self.dy_aff, self.ds_aff);
            const alpha_p_aff = fractionToBoundary(self.x, self.dx_aff, 1.0);
            const alpha_d_aff = fractionToBoundary(self.s, self.ds_aff, 1.0);
            var mu_aff: f64 = 0.0;
            for (0..n) |j| {
                mu_aff += (self.x[j] + alpha_p_aff * self.dx_aff[j]) * (self.s[j] + alpha_d_aff * self.ds_aff[j]);
            }
            mu_aff /= @floatFromInt(@max(n, 1));
            const ratio = if (mu > 0.0) mu_aff / mu else 0.0;
            const sigma = @min(ratio * ratio * ratio, 1.0);

            // Corrector with Mehrotra centrality and second-order term.
            for (0..n) |j| {
                self.rxs[j] = self.x[j] * self.s[j] - sigma * mu + self.dx_aff[j] * self.ds_aff[j];
            }
            self.direction(self.rxs, self.dx, self.dy, self.ds);

            const step = self.options.step_fraction;
            const alpha_p = @min(1.0, step * fractionToBoundary(self.x, self.dx, 1.0));
            const alpha_d = @min(1.0, step * fractionToBoundary(self.s, self.ds, 1.0));
            for (0..n) |j| {
                self.x[j] += alpha_p * self.dx[j];
                self.s[j] += alpha_d * self.ds[j];
            }
            for (0..sf.num_rows) |i| self.y[i] += alpha_d * self.dy[i];
            self.computeResiduals();
        }
        // Report optimality if the final iterate already meets tolerances.
        const pres = maxAbs(self.rp) / self.norm_b;
        const dres = maxAbs(self.rd) / self.norm_c;
        const primal_obj = dot(sf.cost, self.x);
        const dual_obj = dot(sf.rhs, self.y);
        const gap = @abs(primal_obj - dual_obj) / (1.0 + @abs(primal_obj));
        if (pres <= tol.primal_feasibility and dres <= tol.dual_feasibility and gap <= tol.relative_gap) {
            return .{ .status = .optimal, .iterations = iteration };
        }
        return .{ .status = .iteration_limit, .iterations = iteration };
    }

    /// Assemble A·Aᵀ (Θ = I) for the starting-point projection.
    fn assembleUnitNormal(self: *Engine) void {
        const started = clockNs(self.stats_io);
        defer if (self.stats) |stats| addElapsed(self.stats_io, started, &stats.assembly_ns);
        const sf = self.sf;
        @memset(self.normal, 0.0);
        for (0..sf.num_cols) |j| {
            const begin = sf.col_starts[j];
            const end = sf.col_starts[j + 1];
            for (begin..end) |p| {
                const row_p = sf.row_indices[p];
                const vp = sf.values[p];
                for (begin..end) |q| {
                    const row_q = sf.row_indices[q];
                    if (row_q <= row_p) self.normal[self.chol.index(row_p, row_q)] += vp * sf.values[q];
                }
            }
        }
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
    result.primal = try allocator.alloc(f64, problem.num_cols);
    result.dual = try allocator.alloc(f64, problem.num_rows);
    result.reduced_cost = try allocator.alloc(f64, problem.num_cols);

    const sense: f64 = @floatFromInt(@intFromEnum(problem.objective_sense));

    // Primal: reassemble each column from its standard-form parts.
    for (0..problem.num_cols) |j| {
        const record = sf.col_records[j];
        result.primal[j] = switch (record.kind) {
            .fixed => record.base,
            .shifted => record.base + engine.x[record.primary],
            .flipped => record.base - engine.x[record.primary],
            .free_split => engine.x[record.primary] - engine.x[record.negative.?],
        };
    }

    // Row multipliers in min convention, then flipped into the model sense.
    const dual_min = try allocator.alloc(f64, problem.num_rows);
    defer allocator.free(dual_min);
    for (0..problem.num_rows) |i| {
        dual_min[i] = if (sf.row_to_std[i]) |r| engine.y[r] else 0.0;
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
