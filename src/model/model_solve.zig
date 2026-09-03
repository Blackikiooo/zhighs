//! Solve, reset, copy, and MIP-query methods for `Model`.
//!
//! ## Responsibility
//!
//! Owns the main optimization lifecycle: flushing updates, compiling solver IR,
//! dispatching by problem class, resetting solution state, copying committed
//! models, and classifying MIP models.  Presolve-only and other auxiliary
//! workflows belong to `model_advanced.zig`.

const std = @import("std");
const types = @import("types.zig");
const Model = @import("model.zig").Model;
const compile_model_module = @import("compile_model.zig");
const compiled_model_view_module = @import("compiled_model_view.zig");
const matrix = @import("matrix");
const solver = @import("solver");
const mip = @import("mip");

const ModelError = types.ModelError;
const VarType = types.VarType;
const Sense = types.Sense;
const MatrixStore = matrix.MatrixStore;

const SimplexCallbackBridge = struct {
    model: *Model,

    fn iteration(_: solver.LpProgressEventView, context_ptr: ?*anyopaque) solver.LpCallbackAction {
        const self: *SimplexCallbackBridge = @ptrCast(@alignCast(context_ptr orelse return .continue_solve));
        if (self.model.env.callback) |callback| callback(.simplex, self.model.env.usrstate);
        return if (self.model.interrupted.load(.acquire)) .stop else .continue_solve;
    }

    fn log(_: solver.LpProgressEventView, context_ptr: ?*anyopaque) void {
        const self: *SimplexCallbackBridge = @ptrCast(@alignCast(context_ptr orelse return));
        if (self.model.env.log_callback) |callback| callback(.message, self.model.env.log_usrstate);
    }
};

/// Optimise the model.
///
/// Flushes pending changes, compiles the model into solver‑internal IR, and
/// dispatches to the appropriate solver engine based on the problem class.
/// On return, solution attributes are populated and `status` reflects the
/// result.
pub fn optimize(self: *Model) ModelError!void {
    self.interrupted.store(false, .release);
    // Flush any queued modifications.
    try self.updateModel();

    // A genuinely empty 0-by-0 model keeps the historical API error, but a
    // constraint-only model is a valid constant-feasibility LP. The simplex
    // boundary supports zero structural columns and can publish either the
    // constant objective or an empty-row Farkas certificate.
    if (self.num_vars == 0 and self.num_constrs == 0) return error.EmptyModel;

    // Continuous LPs use the zero-copy fast path: committed column SoA and
    // CSC buffers are borrowed directly, while derived row bounds are cached
    // by model revision.
    if (isContinuousLp(self)) {
        const csc = self.matrix.csc();
        const linear_view = self.compiled_view_cache.compileLinearView(self.allocator, .{
            .revision = self.revision,
            .objective_sense = self.sense,
            .objective_offset = self.obj_con,
            .col_cost = self.var_obj,
            .col_lower = self.var_lb,
            .col_upper = self.var_ub,
            .row_sense = self.constr_sense,
            .row_rhs = self.constr_rhs,
            .matrix = matrix.CscView.initAssumeValid(csc.num_rows, csc.num_cols, csc.col_starts, csc.row_indices, csc.values),
        }) catch |err| return switch (err) {
            error.OutOfMemory => error.OutOfMemory,
            error.DimensionMismatch => error.InvalidArgument,
        };
        return solveLinearView(self, linear_view);
    }

    // Unsupported and future problem classes retain the owning compilation
    // path so their transformations can evolve independently of LP views.
    var compiled = compile_model_module.compileModel(self.allocator, self) catch |err| switch (err) {
        error.FeatureNotAvailable => return error.FeatureNotAvailable,
        error.OutOfMemory => return error.OutOfMemory,
        error.ColumnOutOfRange => return error.InvalidArgument,
        error.IndexOutOfRange => return error.InvalidArgument,
    };
    defer compiled.deinit();

    // Dispatch based on problem class.
    switch (compiled.problemClass()) {
        .lp => {
            const linear = compiled.linearData();
            const csc = linear.matrix.csc();
            return solveLinearView(self, .{
                .source_revision = self.revision,
                .objective_sense = linear.objective_sense,
                .objective_offset = linear.objective_offset,
                .num_rows = linear.num_rows,
                .num_cols = linear.num_cols,
                .col_cost = linear.col_cost,
                .col_lower = linear.col_lower,
                .col_upper = linear.col_upper,
                .row_lower = linear.row_lower,
                .row_upper = linear.row_upper,
                .matrix = matrix.CscView.initAssumeValid(csc.num_rows, csc.num_cols, csc.col_starts, csc.row_indices, csc.values),
            });
        },
        .milp => {
            const linear = compiled.linearData();
            return solveMixedIntegerLinear(self, linear);
        },
        .qp,
        .miqp,
        .qcp,
        .miqcp,
        .nlp,
        .minlp,
        => return error.FeatureNotAvailable,
    }
}

fn solveMixedIntegerLinear(self: *Model, linear: *const @import("linear_model.zig").LinearModel) ModelError!void {
    const integrality = self.allocator.alloc(mip.Integrality, linear.num_cols) catch return error.OutOfMemory;
    defer self.allocator.free(integrality);
    for (self.var_type, integrality) |source, *target| target.* = switch (source) {
        .continuous => .continuous,
        .binary => .binary,
        .integer => .integer,
        .semicont => .semi_continuous,
        .semiint => .semi_integer,
    };

    const csc = linear.matrix.csc();
    const relaxation = solver.LpProblemView{
        .num_rows = linear.num_rows,
        .num_cols = linear.num_cols,
        .col_cost = linear.col_cost,
        .col_lower = linear.col_lower,
        .col_upper = linear.col_upper,
        .row_lower = linear.row_lower,
        .row_upper = linear.row_upper,
        .matrix = matrix.CscView.initAssumeValid(csc.num_rows, csc.num_cols, csc.col_starts, csc.row_indices, csc.values),
        .objective_sense = switch (linear.objective_sense) {
            .minimize => .minimize,
            .maximize => .maximize,
        },
        .objective_offset = linear.objective_offset,
    };
    const node_limit_value = self.env.getIntParam("NodeLimit") catch std.math.maxInt(i64);
    const solution_limit_value = self.env.getIntParam("SolutionLimit") catch std.math.maxInt(i64);
    const iteration_limit_value = self.env.getIntParam("IterationLimit") catch std.math.maxInt(i64);
    const time_limit_seconds = self.env.getDblParam("TimeLimit") catch std.math.inf(f64);
    const work_limit_value = self.env.getDblParam("WorkLimit") catch std.math.inf(f64);
    var callback_bridge = SimplexCallbackBridge{ .model = self };
    var result = mip.solve(self.allocator, .{ .relaxation = relaxation, .integrality = integrality }, .{
        .node_limit = nonnegativeLimit(node_limit_value),
        .solution_limit = nonnegativeLimit(solution_limit_value),
        .integer_tolerance = self.env.getDblParam("FeasibilityTol") catch 1e-6,
        .absolute_gap = self.env.getDblParam("MIPGapAbs") catch 1e-10,
        .relative_gap = self.env.getDblParam("MIPGap") catch 1e-4,
        .initial_solution = if (self.mip_start.len == linear.num_cols) self.mip_start else null,
        .lp = .{
            .max_iterations = nonnegativeLimit(iteration_limit_value),
            .time_limit_ns = secondsToNanoseconds(time_limit_seconds),
            .work_limit = workLimit(work_limit_value),
            .interrupt_flag = &self.interrupted,
            .iteration_callback = if (self.env.callback != null) SimplexCallbackBridge.iteration else null,
            .callback_user_data = if (self.env.callback != null) &callback_bridge else null,
        },
    }) catch return error.OutOfMemory;
    defer result.deinit();

    self.clearPublishedResult();
    self.status = mapMipStatus(result.status);
    self.node_count = std.math.cast(i64, result.nodes) orelse std.math.maxInt(i64);
    self.iter_count = std.math.cast(i64, result.lp_iterations) orelse std.math.maxInt(i64);
    if (result.primal.len == linear.num_cols and result.solution_count != 0) {
        const slack = self.allocator.alloc(f64, linear.num_rows) catch return error.OutOfMemory;
        errdefer self.allocator.free(slack);
        @memset(slack, 0.0);
        for (0..linear.num_cols) |column| {
            const begin = csc.col_starts[column];
            const end = csc.col_starts[column + 1];
            for (csc.row_indices[begin..end], csc.values[begin..end]) |row, coefficient|
                slack[row.toUsize()] += coefficient * result.primal[column];
        }
        for (slack, self.constr_sense, self.constr_rhs) |*value, sense, rhs| value.* = switch (sense) {
            .less_equal, .equal => rhs - value.*,
            .greater_equal => value.* - rhs,
        };
        self.solution = result.primal;
        result.primal = &.{};
        self.slack = slack;
        self.obj_val = result.objective_value;
        self.obj_bound = result.objective_bound;
        self.sol_count = std.math.cast(i64, result.solution_count) orelse std.math.maxInt(i64);
    }
}

fn nonnegativeLimit(value: i64) usize {
    if (value <= 0) return 0;
    return std.math.cast(usize, value) orelse std.math.maxInt(usize);
}

fn secondsToNanoseconds(seconds: f64) ?u64 {
    if (!std.math.isFinite(seconds)) return null;
    if (seconds <= 0) return 0;
    const value = seconds * std.time.ns_per_s;
    const maximum: f64 = @floatFromInt(std.math.maxInt(u64));
    return if (value >= maximum) std.math.maxInt(u64) else @intFromFloat(value);
}

fn workLimit(value: f64) ?u64 {
    if (!std.math.isFinite(value)) return null;
    if (value <= 0) return 0;
    const maximum: f64 = @floatFromInt(std.math.maxInt(u64));
    return if (value >= maximum) std.math.maxInt(u64) else @intFromFloat(value);
}

fn mapMipStatus(status: mip.SolveStatus) types.Status {
    return switch (status) {
        .optimal => .optimal,
        .infeasible => .infeasible,
        .unbounded => .unbounded,
        .node_limit => .node_limit,
        .solution_limit => .solution_limit,
        .iteration_limit => .iteration_limit,
        .work_limit => .work_limit,
        .time_limit => .time_limit,
        .interrupted => .interrupted,
        .unbounded_or_infeasible => .inf_or_unbd,
        .invalid_model => .invalid_model,
        .numerical_failure => .numeric,
    };
}

fn isContinuousLp(self: *const Model) bool {
    if (self.q_nz != 0 or self.qconstr_count != 0 or self.sos_count != 0 or
        self.genconstr_count != 0 or self.pwlobj_count != 0)
        return false;
    for (self.var_type) |var_type| {
        if (var_type != .continuous) return false;
    }
    return true;
}

fn solveLinearView(self: *Model, linear: compiled_model_view_module.CompiledLinearModelView) ModelError!void {
    const problem = solver.LpProblemView{
        .num_rows = linear.num_rows,
        .num_cols = linear.num_cols,
        .col_cost = linear.col_cost,
        .col_lower = linear.col_lower,
        .col_upper = linear.col_upper,
        .row_lower = linear.row_lower,
        .row_upper = linear.row_upper,
        .matrix = linear.matrix,
        .objective_sense = switch (linear.objective_sense) {
            .minimize => .minimize,
            .maximize => .maximize,
        },
        .objective_offset = linear.objective_offset,
    };
    const iteration_limit = self.env.getIntParam("IterationLimit") catch std.math.maxInt(i64);
    const max_iterations = if (iteration_limit <= 0)
        0
    else
        std.math.cast(usize, iteration_limit) orelse std.math.maxInt(usize);
    const time_limit_seconds = self.env.getDblParam("TimeLimit") catch std.math.inf(f64);
    const time_limit_ns: ?u64 = if (!std.math.isFinite(time_limit_seconds))
        null
    else if (time_limit_seconds <= 0.0)
        0
    else blk: {
        const nanoseconds = time_limit_seconds * std.time.ns_per_s;
        const max_ns: f64 = @floatFromInt(std.math.maxInt(u64));
        break :blk if (nanoseconds >= max_ns) std.math.maxInt(u64) else @intFromFloat(nanoseconds);
    };
    const work_limit_value = self.env.getDblParam("WorkLimit") catch std.math.inf(f64);
    const work_limit: ?u64 = if (!std.math.isFinite(work_limit_value))
        null
    else if (work_limit_value <= 0.0)
        0
    else blk: {
        const max_work: f64 = @floatFromInt(std.math.maxInt(u64));
        break :blk if (work_limit_value >= max_work) std.math.maxInt(u64) else @intFromFloat(work_limit_value);
    };
    const output_enabled = (self.env.getIntParam("OutputFlag") catch 1) != 0;
    const log_interval_value = self.env.getIntParam("SimplexLogInterval") catch 100;
    const log_interval: u64 = if (log_interval_value <= 0) 1 else @intCast(log_interval_value);
    const pricing_value = self.env.getIntParam("SimplexPricing") catch -1;
    self.lp_session.engine.pricing.rule = switch (pricing_value) {
        0 => .dantzig,
        2 => .steepest_edge,
        3 => .partial,
        4 => .hyper_sparse,
        else => .devex,
    };
    var callback_bridge = SimplexCallbackBridge{ .model = self };
    const warm_basis: ?solver.LpBasisView = if (self.basis_available and
        self.vbasis.len == problem.num_cols and self.cbasis.len == problem.num_rows and self.basis_head.len == problem.num_rows)
        .{
            .structural_status = basisStatusView(self.vbasis),
            .logical_status = basisStatusView(self.cbasis),
            .basic_index = self.basis_head,
        }
    else
        null;
    const control = solver.LpSolveControl{
        .max_iterations = max_iterations,
        .time_limit_ns = time_limit_ns,
        .interrupt_flag = &self.interrupted,
        .initial_basis = warm_basis,
        .work_limit = work_limit,
        .iteration_callback = if (self.env.callback != null) SimplexCallbackBridge.iteration else null,
        .callback_user_data = if (self.env.callback != null) &callback_bridge else null,
        .log_level = if (output_enabled and self.env.log_callback != null) .iterations else .off,
        .log_callback = if (output_enabled and self.env.log_callback != null) SimplexCallbackBridge.log else null,
        .log_user_data = if (output_enabled and self.env.log_callback != null) &callback_bridge else null,
        .log_interval_work = log_interval,
    };
    const solve_status = self.lp_session.solve(problem, control, .{
        .structure = self.revisions.structure,
        .matrix_values = self.revisions.matrix_values,
        .bounds = self.revisions.bounds,
        .objective = self.revisions.objective,
    });

    self.status = mapLpStatus(solve_status);
    self.iter_count = std.math.cast(i64, self.lp_session.engine.iterations) orelse std.math.maxInt(i64);
    switch (solve_status) {
        .optimal => {
            const result_view = self.lp_session.resultView(problem, solve_status) orelse {
                self.status = .numeric;
                self.clearPublishedResult();
                return;
            };
            if (!result_view.primal_valid or !result_view.dual_valid or !result_view.basis_valid) {
                self.status = .numeric;
                self.clearPublishedResult();
                return;
            }
            publishLpSolution(self, problem, result_view) catch |err| {
                self.status = .numeric;
                return err;
            };
        },
        .infeasible, .unbounded => {
            const result_view = self.lp_session.resultView(problem, solve_status) orelse {
                self.status = .numeric;
                self.clearPublishedResult();
                return;
            };
            const certificate = if (solve_status == .unbounded)
                result_view.unbounded_ray
            else
                result_view.infeasibility_ray;
            const expected_len = if (solve_status == .unbounded) problem.num_cols else problem.num_rows;
            if (certificate.len != expected_len) {
                self.status = .numeric;
                self.clearPublishedResult();
                return;
            }
            publishLpCertificate(self, solve_status, certificate) catch |err| {
                self.status = .numeric;
                return err;
            };
        },
        else => self.clearPublishedResult(),
    }
}

fn mapLpStatus(status: solver.LpSolveStatus) types.Status {
    return switch (status) {
        .optimal => .optimal,
        .infeasible => .infeasible,
        .unbounded => .unbounded,
        .iteration_limit => .iteration_limit,
        .work_limit => .work_limit,
        .time_limit => .time_limit,
        .interrupted => .interrupted,
        .invalid_model => .invalid_model,
        .numerical_failure, .not_implemented => .numeric,
    };
}

test "LP invalid model status remains distinct from numerical failure" {
    try std.testing.expectEqual(types.Status.invalid_model, mapLpStatus(.invalid_model));
    try std.testing.expectEqualStrings("INVALID_MODEL", types.Status.invalid_model.label());
}

fn mapBasisStatus(status: solver.LpBasisStatus) types.BasisStatus {
    return switch (status) {
        .basic => .basic,
        .at_lower, .fixed => .non_basic_lower,
        .at_upper => .non_basic_upper,
        .superbasic, .free => .super_non_basic,
    };
}

/// Copy one terminal LP certificate into stable public model storage. The
/// simplex Farkas convention is positive for row lower bounds and negative for
/// row upper bounds; `FarkasDual` uses the public/Gurobi convention, so its
/// sign is reversed. An unbounded primal direction needs no transformation.
fn publishLpCertificate(
    self: *Model,
    status: solver.LpSolveStatus,
    certificate: []const f64,
) ModelError!void {
    const owned = self.allocator.dupe(f64, certificate) catch return error.OutOfMemory;
    errdefer self.allocator.free(owned);
    if (status == .infeasible) {
        for (owned) |*value| value.* = -value.*;
    }
    self.clearPublishedResult();
    switch (status) {
        .unbounded => self.unbounded_ray = owned,
        .infeasible => self.farkas_dual = owned,
        else => unreachable,
    }
}

/// Copy the transient engine result exactly once into stable public model
/// attributes. All allocations complete before the previous solution is
/// released, so an allocation failure cannot partially publish a solve.
fn publishLpSolution(self: *Model, problem: solver.LpProblemView, result: solver.LpSolveResultView) ModelError!void {
    const allocator = self.allocator;
    const primal = allocator.alloc(f64, problem.num_cols) catch return error.OutOfMemory;
    errdefer allocator.free(primal);
    const reduced_cost = allocator.alloc(f64, problem.num_cols) catch return error.OutOfMemory;
    errdefer allocator.free(reduced_cost);
    const slack = allocator.alloc(f64, problem.num_rows) catch return error.OutOfMemory;
    errdefer allocator.free(slack);
    const pi = allocator.alloc(f64, problem.num_rows) catch return error.OutOfMemory;
    errdefer allocator.free(pi);
    const vbasis = allocator.alloc(types.BasisStatus, problem.num_cols) catch return error.OutOfMemory;
    errdefer allocator.free(vbasis);
    const cbasis = allocator.alloc(types.BasisStatus, problem.num_rows) catch return error.OutOfMemory;
    errdefer allocator.free(cbasis);
    const basis_head = allocator.dupe(u32, result.basic_index) catch return error.OutOfMemory;
    errdefer allocator.free(basis_head);

    @memcpy(primal, result.primal);
    const objective_sign: f64 = switch (problem.objective_sense) {
        .minimize => 1.0,
        .maximize => -1.0,
    };
    for (reduced_cost, result.reduced_cost) |*published, internal| published.* = objective_sign * internal;
    for (pi, result.dual, result.row_scale) |*published, internal, scale| published.* = objective_sign * scale * internal;
    for (vbasis, result.structural_status) |*published, internal| published.* = mapBasisStatus(internal);
    for (cbasis, result.logical_status) |*published, internal| published.* = mapBasisStatus(internal);

    @memset(slack, 0.0);
    for (0..problem.num_cols) |column| {
        const begin = problem.matrix.col_starts[column];
        const end = problem.matrix.col_starts[column + 1];
        for (problem.matrix.row_indices[begin..end], problem.matrix.values[begin..end]) |row, coefficient| {
            slack[row.toUsize()] += coefficient * primal[column];
        }
    }
    for (slack, self.constr_sense, self.constr_rhs) |*value, sense, rhs| {
        value.* = switch (sense) {
            .less_equal, .equal => rhs - value.*,
            .greater_equal => value.* - rhs,
        };
    }

    self.clearPublishedResult();
    self.solution = primal;
    self.reduced_cost = reduced_cost;
    self.slack = slack;
    self.pi = pi;
    self.vbasis = vbasis;
    self.cbasis = cbasis;
    self.basis_head = basis_head;
    self.basis_available = true;
    self.obj_val = result.objective_value;
    // For a proven-optimal continuous LP the incumbent objective and global
    // bound coincide exactly at the public solver boundary.
    self.obj_bound = result.objective_value;
    self.sol_count = 1;
}

fn basisStatusView(status: []const types.BasisStatus) []const solver.LpBasisStatus {
    comptime {
        std.debug.assert(@sizeOf(types.BasisStatus) == @sizeOf(solver.LpBasisStatus));
        std.debug.assert(@intFromEnum(types.BasisStatus.basic) == @intFromEnum(solver.LpBasisStatus.basic));
        std.debug.assert(@intFromEnum(types.BasisStatus.non_basic_lower) == @intFromEnum(solver.LpBasisStatus.at_lower));
        std.debug.assert(@intFromEnum(types.BasisStatus.non_basic_upper) == @intFromEnum(solver.LpBasisStatus.at_upper));
        std.debug.assert(@intFromEnum(types.BasisStatus.super_non_basic) == @intFromEnum(solver.LpBasisStatus.superbasic));
    }
    const pointer: [*]const solver.LpBasisStatus = @ptrCast(status.ptr);
    return pointer[0..status.len];
}

/// Discard the solution and reset the status to LOADED
pub fn reset(self: *Model, clear_all: bool) void {
    _ = clear_all;
    self.clearPublishedResult();
    self.status = .loaded;
    self.obj_bound = 0.0;
    self.iter_count = 0;
    self.node_count = 0;
    self.bar_iter_count = 0;
    self.interrupted.store(false, .release);
}

/// Create an independent copy of the model
pub fn copy(self: *Model, new_name: []const u8) ModelError!Model {
    var new = try Model.init(self.allocator, self.env, new_name);
    errdefer new.deinit();

    // Flush pending changes first.
    try self.updateModel();

    // Clone variable data.
    if (self.num_vars > 0) {
        new.var_lb = try self.allocator.dupe(f64, self.var_lb);
        new.var_ub = try self.allocator.dupe(f64, self.var_ub);
        new.var_obj = try self.allocator.dupe(f64, self.var_obj);
        new.var_type = try self.allocator.dupe(VarType, self.var_type);
        new.mip_start = try self.allocator.dupe(f64, self.mip_start);
        new.p_start = try self.allocator.dupe(f64, self.p_start);
        new.d_start = try self.allocator.dupe(f64, self.d_start);
        new.var_names = try self.allocator.alloc(?[]const u8, self.num_vars);
        for (self.var_names, 0..) |n, i| {
            new.var_names[i] = if (n) |s| try self.allocator.dupe(u8, s) else null;
        }
        new.num_vars = self.num_vars;
    }

    // Clone constraint data.
    if (self.num_constrs > 0) {
        new.constr_sense = try self.allocator.dupe(Sense, self.constr_sense);
        new.constr_rhs = try self.allocator.dupe(f64, self.constr_rhs);
        new.constr_names = try self.allocator.alloc(?[]const u8, self.num_constrs);
        for (self.constr_names, 0..) |n, i| {
            new.constr_names[i] = if (n) |s| try self.allocator.dupe(u8, s) else null;
        }
        new.num_constrs = self.num_constrs;
    }

    // Clone the matrix.
    const csc_src = self.matrix.csc();
    const csc_copy = try csc_src.clone(self.allocator);
    // Free the initial zero matrix from init().
    new.matrix.deinit(self.allocator);
    new.matrix = MatrixStore.initAssumeValid(csc_copy);

    new.sense = self.sense;
    // Transient solve state belongs to the source model's persistent session;
    // a structural copy starts unsolved and owns no stale solution/certificate.
    new.status = .loaded;
    new.revision = self.revision;
    return new;
}

/// Return whether the model contains any non-continuous variable.
pub fn isMip(self: Model) bool {
    for (self.var_type) |vt| {
        if (vt != .continuous) return true;
    }
    return false;
}
