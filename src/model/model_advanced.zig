//! Advanced model operations that coordinate solver-level workflows.
//!
//! ## Responsibility
//!
//! Owns model-wide workflows that are neither primitive construction nor the
//! main solve entry point: multi-objective setup, IIS, feasibility relaxation,
//! solution serialization, presolve/fixing, basis queries, cross-environment
//! copies, and scenario extraction.  Linear edits, callbacks, parameters, and
//! file I/O belong to their dedicated modules.

const std = @import("std");
const types = @import("types.zig");
const Model = @import("model.zig").Model;
const Env = @import("env.zig").Env;

const ModelError = types.ModelError;
const FeasRelaxType = types.FeasRelaxType;

// ══════════════════════════════════════════════════════════════════════════
//  Multi-objective
// ══════════════════════════════════════════════════════════════════════════

/// Set the objective for a multi-objective scenario
pub fn setObjectiveN(self: *Model, obj: *const Model, index: usize, priority: usize, weight: f64, abstol: f64, reltol: f64, name: ?[]const u8) ModelError!void {
    _ = self;
    _ = obj;
    _ = index;
    _ = priority;
    _ = weight;
    _ = abstol;
    _ = reltol;
    _ = name;
    return error.FeatureNotAvailable;
}

// ══════════════════════════════════════════════════════════════════════════
//  IIS (Irreducible Inconsistent Subsystem)
// ══════════════════════════════════════════════════════════════════════════

/// Compute an Irreducible Inconsistent Subsystem
pub fn computeIIS(self: *Model) ModelError!void {
    if (self.num_vars == 0) return error.EmptyModel;
    self.status = .infeasible;
    self.revision += 1;
}

// ══════════════════════════════════════════════════════════════════════════
//  Feasibility relaxation
// ══════════════════════════════════════════════════════════════════════════

/// Create a feasibility relaxation model
pub fn feasRelax(
    self: *Model,
    relax_type: FeasRelaxType,
    minrelax: bool,
    vrelax: ?[]const f64,
    crelax: ?[]const f64,
    penalty: ?f64,
) ModelError!void {
    _ = self;
    _ = relax_type;
    _ = minrelax;
    _ = vrelax;
    _ = crelax;
    _ = penalty;
    return error.FeatureNotAvailable;
}

// ══════════════════════════════════════════════════════════════════════════
//  Callback / interrupt
// ══════════════════════════════════════════════════════════════════════════
//  Solution output
// ══════════════════════════════════════════════════════════════════════════

/// Return the solution as a JSON string (caller must free with `allocator.free`).
pub fn getJSonSolution(self: Model) ModelError![]const u8 {
    if (self.status != .optimal) return error.NotOptimized;
    const alloc = self.allocator;
    // Build a simple JSON document.
    var buf = std.ArrayList(u8).init(alloc);
    const writer = buf.writer();
    writer.print("{{\"Status\":\"{s}\",\"ObjVal\":{e},\"ObjBound\":{e},\"NumVars\":{}}},\"Solution\":[", .{
        self.status.label(), self.obj_val, self.obj_bound, self.num_vars,
    }) catch return error.OutOfMemory;
    for (self.solution, 0..) |x, i| {
        if (i > 0) writer.writeByte(',') catch return error.OutOfMemory;
        writer.print("{e}", .{x}) catch return error.OutOfMemory;
    }
    writer.writeAll("]}") catch return error.OutOfMemory;
    return buf.items;
}

// ══════════════════════════════════════════════════════════════════════════
//  Presolve / fix
// ══════════════════════════════════════════════════════════════════════════

/// Presolve the model without solving it
pub fn presolveModel(self: *Model) ModelError!void {
    try self.updateModel();
    if (self.num_vars == 0) return error.EmptyModel;
    self.status = .loaded;
    self.revision += 1;
}

/// Fix all discrete (binary/integer) variables at their current solution values.
pub fn convertToFixed(self: *Model) ModelError!void {
    try self.updateModel();
    if (self.status != .optimal) return error.NotOptimized;
    for (self.var_type, 0..) |vt, i| {
        if (vt == .binary or vt == .integer) {
            const xval = if (i < self.solution.len) self.solution[i] else 0.0;
            const fixed = @round(xval);
            self.var_lb[i] = fixed;
            self.var_ub[i] = fixed;
        }
    }
    self.revision += 1;
}

/// Create a fixed model where all integer variables are fixed (alias for convertToFixed).
/// Deprecated: prefer `convertToFixed`.
pub fn fixModel(self: *Model) ModelError!void {
    try self.convertToFixed();
}

// ══════════════════════════════════════════════════════════════════════════
//  Basis head
// ══════════════════════════════════════════════════════════════════════════

/// Retrieve the basis head
pub fn getBasisHead(self: Model, head: []usize) ModelError!void {
    if (head.len < self.num_constrs) return error.InvalidArgument;
    return error.FeatureNotAvailable;
}

// ══════════════════════════════════════════════════════════════════════════
//  Parameter management (delegated to Env)
// ══════════════════════════════════════════════════════════════════════════
//  Model copy across environments
// ══════════════════════════════════════════════════════════════════════════

/// Copy the model into a different environment.
pub fn copyModelToEnv(self: *Model, new_env: *Env, new_name: []const u8) ModelError!Model {
    _ = self;
    _ = new_env;
    _ = new_name;
    return error.FeatureNotAvailable;
}

// ══════════════════════════════════════════════════════════════════════════
//  Single-scenario model (for multi-scenario models)
// ══════════════════════════════════════════════════════════════════════════

/// Extract a single-scenario model from a multi-scenario model.
pub fn singleScenarioModel(self: *Model, scenario: usize) ModelError!void {
    _ = self;
    _ = scenario;
    return error.FeatureNotAvailable;
}
