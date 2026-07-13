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
const foundation = @import("foundation");
const matrix = @import("matrix");
const solver = @import("solver");

const ModelError = types.ModelError;
const VarType = types.VarType;
const Sense = types.Sense;
const RowId = foundation.RowId;
const CscMatrix = matrix.CscMatrix;
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

    if (self.num_vars == 0) return error.EmptyModel;

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
            .matrix = .{
                .num_rows = csc.num_rows,
                .num_cols = csc.num_cols,
                .col_starts = csc.col_starts,
                .row_indices = csc.row_indices,
                .values = csc.values,
            },
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
                .matrix = .{
                    .num_rows = csc.num_rows,
                    .num_cols = csc.num_cols,
                    .col_starts = csc.col_starts,
                    .row_indices = csc.row_indices,
                    .values = csc.values,
                },
            });
        },
        .milp,
        .qp,
        .miqp,
        .qcp,
        .miqcp,
        .nlp,
        .minlp,
        => return error.FeatureNotAvailable,
    }
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
    if (solve_status == .optimal) {
        const result_view = self.lp_session.resultView(problem, solve_status) orelse {
            self.status = .numeric;
            clearPublishedSolution(self);
            return;
        };
        publishLpSolution(self, problem, result_view) catch |err| {
            self.status = .numeric;
            return err;
        };
    } else {
        clearPublishedSolution(self);
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
        .numerical_failure, .not_implemented => .numeric,
    };
}

fn mapBasisStatus(status: solver.LpBasisStatus) types.BasisStatus {
    return switch (status) {
        .basic => .basic,
        .at_lower, .fixed => .non_basic_lower,
        .at_upper => .non_basic_upper,
        .superbasic, .free => .super_non_basic,
    };
}

fn clearPublishedSolution(self: *Model) void {
    self.allocator.free(self.solution);
    self.allocator.free(self.reduced_cost);
    self.allocator.free(self.slack);
    self.allocator.free(self.pi);
    self.allocator.free(self.vbasis);
    self.allocator.free(self.cbasis);
    self.allocator.free(self.basis_head);
    self.solution = &.{};
    self.reduced_cost = &.{};
    self.slack = &.{};
    self.pi = &.{};
    self.vbasis = &.{};
    self.cbasis = &.{};
    self.basis_head = &.{};
    self.basis_available = false;
    self.obj_val = 0.0;
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

    clearPublishedSolution(self);
    self.solution = primal;
    self.reduced_cost = reduced_cost;
    self.slack = slack;
    self.pi = pi;
    self.vbasis = vbasis;
    self.cbasis = cbasis;
    self.basis_head = basis_head;
    self.basis_available = true;
    self.obj_val = result.objective_value;
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
    self.status = .loaded;
    self.obj_val = 0.0;
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
    const csc_copy = CscMatrix{
        .num_rows = csc_src.num_rows,
        .num_cols = csc_src.num_cols,
        .col_starts = try self.allocator.dupe(usize, csc_src.col_starts),
        .row_indices = try self.allocator.dupe(RowId, csc_src.row_indices),
        .values = try self.allocator.dupe(f64, csc_src.values),
    };
    // Free the initial zero matrix from init().
    new.matrix.deinit(self.allocator);
    new.matrix = MatrixStore.initAssumeValid(csc_copy);

    new.sense = self.sense;
    new.status = self.status;
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
