//! Model‑to‑CompiledModel compilation path.
//!
//! Translates a fully‑committed Gurobi‑style `Model` into the solver‑internal
//! [`CompiledModel`] IR by deep‑copying all data and applying discriminant
//! mappings (VarType → Integrality, Sense → row bounds, etc.).
//!
//! ## Preconditions
//!
//! The source `Model` must be fully committed — call `model.updateModel()`
//! before `compileModel`.  The `optimize` method in `model_solve.zig` already
//! guarantees this; standalone callers must ensure it themselves.
//!
//! ## Unsupported features
//!
//! Models containing SOS constraints, general constraints, or piecewise‑linear
//! objective data cause `CompileError.FeatureNotAvailable`.  These features
//! may be supported in a later phase.
//!
//! Experimental API.

const std = @import("std");
const types = @import("types.zig");
const Model = @import("model.zig").Model;
const problem_class = @import("problem_class.zig");
const linear_model_module = @import("linear_model.zig");
const linear_model_builder_module = @import("linear_model_builder.zig");
const hessian_module = @import("hessian.zig");
const quadratic_model_module = @import("quadratic_model.zig");
const compiled_model_module = @import("compiled_model.zig");
const matrix = @import("matrix");
const foundation = @import("foundation");

const ObjectiveSense = types.ObjectiveSense;
const VarType = types.VarType;
const Sense = types.Sense;
const INFINITY = types.INFINITY;
const Integrality = linear_model_module.Integrality;
const LinearModel = linear_model_module.LinearModel;
const ProblemClass = problem_class.ProblemClass;
const Hessian = hessian_module.Hessian;
const HessianFormat = hessian_module.HessianFormat;
const QuadraticConstraint = quadratic_model_module.QuadraticConstraint;
const QuadraticModel = quadratic_model_module.QuadraticModel;
const CompiledModel = compiled_model_module.CompiledModel;
const CscMatrix = matrix.CscMatrix;
const MatrixStore = matrix.MatrixStore;
const RowId = foundation.RowId;
const ColId = foundation.ColId;

// ── CompileError ─────────────────────────────────────────────────────────────

/// Errors that can occur during model compilation.
pub const CompileError = error{
    /// The model contains features (SOS, general constraints, PWL, multi‑obj)
    /// that cannot yet be compiled into the solver‑internal IR.
    FeatureNotAvailable,

    /// A memory allocation failed.
    OutOfMemory,

    /// A column index in the quadratic data is out of range.
    ColumnOutOfRange,

    /// Row index out of range.
    IndexOutOfRange,
};

// ── compileModel ─────────────────────────────────────────────────────────────

/// Compile a fully‑committed Gurobi‑style `Model` into a `CompiledModel`.
///
/// The source model is unchanged (all data is deep‑copied).
pub fn compileModel(allocator: std.mem.Allocator, model: *const Model) (CompileError || std.mem.Allocator.Error)!CompiledModel {
    // ── Check for unsupported features ──────────────────────────────────────
    if (model.sos_count > 0) return error.FeatureNotAvailable;
    if (model.genconstr_count > 0) return error.FeatureNotAvailable;
    if (model.pwlobj_count > 0) return error.FeatureNotAvailable;

    const num_vars = model.num_vars;
    const num_constrs = model.num_constrs;

    // ── Build the LinearModel in a block ──────────────────────────────────
    //
    // By constructing the LinearModel inside a `blk`, all intermediate
    // errdefers (for the individual array allocations, the CscMatrix, and
    // the MatrixStore) are scoped to that block.  When the block returns
    // the LinearModel, those errdefers are consumed — they will NOT fire on
    // a subsequent error outside the block.  This avoids double‑free when
    // later steps (Hessian, QuadraticConstraint) fail.

    var linear = blk: {
        // Objective sense + offset
        const sense = model.sense;
        const offset = model.obj_con;

        // Column data (deep‑copy)
        const col_cost = try allocator.dupe(f64, model.var_obj[0..num_vars]);
        errdefer allocator.free(col_cost);

        const col_lower = try allocator.dupe(f64, model.var_lb[0..num_vars]);
        errdefer allocator.free(col_lower);

        const col_upper = try allocator.dupe(f64, model.var_ub[0..num_vars]);
        errdefer allocator.free(col_upper);

        // Row bounds (Sense + RHS → lower/upper)
        const row_lower = try allocator.alloc(f64, num_constrs);
        errdefer allocator.free(row_lower);
        const row_upper = try allocator.alloc(f64, num_constrs);
        errdefer allocator.free(row_upper);

        for (0..num_constrs) |i| {
            const s = model.constr_sense[i];
            const rhs = model.constr_rhs[i];
            switch (s) {
                .less_equal => {
                    row_lower[i] = -INFINITY;
                    row_upper[i] = rhs;
                },
                .equal => {
                    row_lower[i] = rhs;
                    row_upper[i] = rhs;
                },
                .greater_equal => {
                    row_lower[i] = rhs;
                    row_upper[i] = INFINITY;
                },
            }
        }

        // Integrality (VarType → Integrality mapping)
        const has_integer = hasNonContinuousVars(model.var_type[0..num_vars]);
        const integrality: ?[]Integrality = if (has_integer) int_blk: {
            const int = try allocator.alloc(Integrality, num_vars);
            errdefer allocator.free(int);
            for (0..num_vars) |i| {
                int[i] = switch (model.var_type[i]) {
                    .continuous => Integrality.continuous,
                    .binary => Integrality.integer,
                    .integer => Integrality.integer,
                    .semicont => Integrality.semi_continuous,
                    .semiint => Integrality.semi_integer,
                };
            }
            break :int_blk int;
        } else null;

        // ── Copy the constraint matrix ─────────────────────────────────────
        const src_csc = model.matrix.csc();

        const csc_col_starts = try allocator.dupe(usize, src_csc.col_starts);
        errdefer allocator.free(csc_col_starts);
        const csc_row_indices = try allocator.dupe(RowId, src_csc.row_indices);
        errdefer allocator.free(csc_row_indices);
        const csc_values = try allocator.dupe(f64, src_csc.values);
        errdefer allocator.free(csc_values);

        // Build CscMatrix + MatrixStore.
        var csc_buf = CscMatrix{
            .num_rows = src_csc.num_rows,
            .num_cols = src_csc.num_cols,
            .col_starts = csc_col_starts,
            .row_indices = csc_row_indices,
            .values = csc_values,
        };
        _ = &csc_buf;
        errdefer csc_buf.deinit(allocator);

        var ms_buf = MatrixStore.initAssumeValid(csc_buf);
        _ = &ms_buf;
        errdefer ms_buf.deinit(allocator);

        // Assemble and return the LinearModel from the block.
        break :blk LinearModel{
            .allocator = allocator,
            .objective_sense = sense,
            .objective_offset = offset,
            .num_rows = num_constrs,
            .num_cols = num_vars,
            .col_cost = col_cost,
            .col_lower = col_lower,
            .col_upper = col_upper,
            .row_lower = row_lower,
            .row_upper = row_upper,
            .integrality = integrality,
            .matrix = ms_buf,
            .revision = 1,
        };
    };
    // ── Block ends — errdefers inside are consumed. linear owns all data. ──

    // ── Check for quadratic objective (Hessian) ────────────────────────────
    const has_qp = model.q_nz > 0;
    const has_qconstr = model.qconstr_count > 0;

    if (!has_qp and !has_qconstr) {
        // Pure LP / MILP.
        return CompiledModel{ .linear = linear };
    }

    // Build objective Hessian.  On failure, clean up the LinearModel manually.
    var obj_hessian: ?Hessian = if (has_qp) buildHessianFromTriples(
        allocator,
        num_vars,
        model.q_row[0..model.q_nz],
        model.q_col[0..model.q_nz],
        model.q_val[0..model.q_nz],
    ) catch |err| {
        linear.deinit();
        return err;
    } else null;
    _ = &obj_hessian;

    // ── Build QuadraticConstraint list ─────────────────────────────────────
    var qconstraints: []QuadraticConstraint = &.{};
    if (has_qconstr) {
        qconstraints = allocator.alloc(QuadraticConstraint, model.qconstr_count) catch |err| {
            if (obj_hessian) |*h| h.deinit();
            linear.deinit();
            return err;
        };
        var built: usize = 0;
        errdefer {
            for (qconstraints[0..built]) |*qc| qc.deinit();
            allocator.free(qconstraints);
            if (obj_hessian) |*h| h.deinit();
            linear.deinit();
        }
        for (0..model.qconstr_count) |i| {
            const qb = model.qconstr_qbegin[i];
            const qe = model.qconstr_qbegin[i + 1];
            const lb = model.qconstr_lbegin[i];
            const le = model.qconstr_lbegin[i + 1];
            qconstraints[i] = buildQuadraticConstraint(
                allocator,
                num_vars,
                model.qconstr_qrow[qb..qe],
                model.qconstr_qcol[qb..qe],
                model.qconstr_qval[qb..qe],
                model.qconstr_lind[lb..le],
                model.qconstr_lval[lb..le],
                model.qconstr_sense[i],
                model.qconstr_rhs[i],
            ) catch |err| return err;
            built += 1;
        }
    }

    // ── Assemble the QuadraticModel ────────────────────────────────────────
    const quad = QuadraticModel{
        .allocator = allocator,
        .linear = linear,
        .objective_hessian = obj_hessian,
        .quadratic_constraints = qconstraints,
    };

    return CompiledModel{ .quadratic = quad };
}

// ── Helpers ──────────────────────────────────────────────────────────────────

/// Returns `true` when any variable type is non‑continuous.
fn hasNonContinuousVars(var_type: []const VarType) bool {
    for (var_type) |vt| {
        if (vt != .continuous) return true;
    }
    return false;
}

/// Build a lower‑triangle CSC Hessian from COO triplets.
///
/// `rows` and `cols` are `i32` index arrays (as stored in the Model's
/// `q_row` / `q_col` fields).  All triples must be in the lower triangle
/// of the `dimension × dimension` matrix.
fn buildHessianFromTriples(
    allocator: std.mem.Allocator,
    dimension: usize,
    rows: []const i32,
    cols: []const i32,
    values: []const f64,
) (CompileError || std.mem.Allocator.Error)!Hessian {
    if (dimension == 0) {
        return Hessian.initEmpty(allocator);
    }

    if (rows.len != values.len or cols.len != values.len) {
        // All three must have the same length; if not, we can still proceed
        // using the minimum length.
    }

    const nnz_in = @min(rows.len, @min(cols.len, values.len));

    // Count non‑zeros per column (skip explicit zeros).
    var per_col = try allocator.alloc(usize, dimension);
    defer allocator.free(per_col);
    @memset(per_col, 0);

    for (0..nnz_in) |i| {
        const v = values[i];
        if (v == 0.0) continue;
        const col = @as(usize, @intCast(cols[i]));
        if (col >= dimension) return error.ColumnOutOfRange;
        per_col[col] += 1;
    }

    // Build starts (prefix sum).
    const starts = try allocator.alloc(usize, dimension + 1);
    errdefer allocator.free(starts);
    starts[0] = 0;
    for (0..dimension) |j| {
        starts[j + 1] = starts[j] + per_col[j];
    }

    const total_nnz = starts[dimension];
    const indices = try allocator.alloc(ColId, total_nnz);
    errdefer allocator.free(indices);
    const vals = try allocator.alloc(f64, total_nnz);
    errdefer allocator.free(vals);

    // Fill using a cursor per column (copy of starts).
    var cursor = try allocator.dupe(usize, starts[0..dimension]);
    defer allocator.free(cursor);

    for (0..nnz_in) |i| {
        const v = values[i];
        if (v == 0.0) continue;
        const col = @as(usize, @intCast(cols[i]));
        const row = @as(usize, @intCast(rows[i]));
        if (row >= dimension) return error.ColumnOutOfRange;

        const pos = cursor[col];
        indices[pos] = ColId.fromUsizeAssumeValid(row);
        vals[pos] = v;
        cursor[col] = pos + 1;
    }

    return Hessian{
        .allocator = allocator,
        .dimension = dimension,
        .starts = starts,
        .indices = indices,
        .values = vals,
        .format = .triangular,
    };
}

/// Build a single `QuadraticConstraint` from Model's packed‑data fields.
fn buildQuadraticConstraint(
    allocator: std.mem.Allocator,
    num_cols: usize,
    qrows: []const i32,
    qcols: []const i32,
    qvals: []const f64,
    lind: []const usize,
    lval: []const f64,
    sense: Sense,
    rhs: f64,
) (CompileError || std.mem.Allocator.Error)!QuadraticConstraint {
    // Build the quadratic part (Hessian).
    var hessian = try buildHessianFromTriples(allocator, num_cols, qrows, qcols, qvals);
    _ = &hessian;
    errdefer hessian.deinit();

    // Build the linear part.
    const linear_indices = try allocator.alloc(ColId, lind.len);
    errdefer allocator.free(linear_indices);
    const linear_values = try allocator.dupe(f64, lval);
    errdefer allocator.free(linear_values);

    for (lind, 0..) |idx, j| {
        if (idx >= num_cols) {
            allocator.free(linear_indices);
            allocator.free(linear_values);
            hessian.deinit();
            return error.ColumnOutOfRange;
        }
        linear_indices[j] = ColId.fromUsizeAssumeValid(idx);
    }

    // Sense + RHS → lower/upper bounds.
    const lower: f64 = switch (sense) {
        .less_equal => -INFINITY,
        .equal => rhs,
        .greater_equal => rhs,
    };
    const upper: f64 = switch (sense) {
        .less_equal => rhs,
        .equal => rhs,
        .greater_equal => INFINITY,
    };

    return QuadraticConstraint{
        .allocator = allocator,
        .linear_indices = linear_indices,
        .linear_values = linear_values,
        .hessian = hessian,
        .lower = lower,
        .upper = upper,
    };
}

// ── Tests ────────────────────────────────────────────────────────────────────

const testing = std.testing;
const Env = @import("env.zig").Env;
const LinearModelBuilder = linear_model_builder_module.LinearModelBuilder;

/// Helper: build a simple LP model for testing compilation.
fn buildTestLpModel(allocator: std.mem.Allocator) !Model {
    var env = try Env.initSimple(allocator);
    errdefer env.deinit();

    var model = try Model.init(allocator, &env, "test_lp");
    errdefer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, "c0");
    try model.addConstr(0, &.{}, &.{}, .equal, 5.0, "c1");
    try model.addVar(2, &[_]usize{ 0, 1 }, &[_]f64{ 1.0, 2.0 }, 1.0, 0.0, 10.0, .continuous, "x1");
    try model.addVar(2, &[_]usize{ 0, 1 }, &[_]f64{ 3.0, 4.0 }, 2.0, -5.0, 5.0, .continuous, "x2");
    try model.addVar(0, &.{}, &.{}, 3.0, 0.0, INFINITY, .integer, "x3");
    try model.updateModel();
    return model;
}

test "compileModel empty model compiles to LP" {
    var env = try Env.initSimple(testing.allocator);
    defer env.deinit();
    var model = try Model.init(testing.allocator, &env, "empty");
    defer model.deinit();

    var cm = try compileModel(testing.allocator, &model);
    defer cm.deinit();

    try testing.expectEqual(ProblemClass.lp, cm.problemClass());
    try testing.expectEqual(@as(usize, 0), cm.linearData().num_cols);
    try testing.expectEqual(@as(usize, 0), cm.linearData().num_rows);
}

test "compileModel simple LP" {
    var env = try Env.initSimple(testing.allocator);
    defer env.deinit();
    var model = try Model.init(testing.allocator, &env, "lp");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 100.0, "c0");
    try model.addVar(1, &[_]usize{0}, &[_]f64{5.0}, 1.0, 0.0, 10.0, .continuous, null);
    try model.addVar(1, &[_]usize{0}, &[_]f64{3.0}, 2.0, 0.0, 20.0, .continuous, null);
    try model.setDblAttr(.obj_con, 1.5);
    try model.updateModel();

    var cm = try compileModel(testing.allocator, &model);
    defer cm.deinit();

    try testing.expectEqual(ProblemClass.lp, cm.problemClass());
    const lm = cm.linearData();
    try testing.expectEqual(@as(usize, 1), lm.num_rows);
    try testing.expectEqual(@as(usize, 2), lm.num_cols);
    try testing.expectEqual(@as(f64, 1.0), lm.col_cost[0]);
    try testing.expectEqual(@as(f64, 2.0), lm.col_cost[1]);
    try testing.expectEqual(@as(f64, 1.5), lm.objective_offset);
    try testing.expect(lm.integrality == null);
    try testing.expect(!lm.isMixedInteger());
}

test "compileModel MILP with binary and integer" {
    var env = try Env.initSimple(testing.allocator);
    defer env.deinit();
    var model = try Model.init(testing.allocator, &env, "milp");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, "c0");
    try model.addVar(1, &[_]usize{0}, &[_]f64{1.0}, 1.0, 0.0, 1.0, .binary, "x1");
    try model.addVar(1, &[_]usize{0}, &[_]f64{2.0}, 2.0, -5.0, 10.0, .integer, "x2");
    try model.addVar(1, &[_]usize{0}, &[_]f64{3.0}, 3.0, 0.0, 5.0, .semicont, "x3");
    try model.updateModel();

    var cm = try compileModel(testing.allocator, &model);
    defer cm.deinit();

    try testing.expectEqual(ProblemClass.milp, cm.problemClass());
    const lm = cm.linearData();
    try testing.expect(lm.isMixedInteger());
    try testing.expect(lm.integrality != null);

    const int = lm.integrality.?;
    try testing.expectEqual(Integrality.integer, int[0]); // binary → integer
    try testing.expectEqual(Integrality.integer, int[1]); // integer → integer
    try testing.expectEqual(Integrality.semi_continuous, int[2]); // semicont → semi_continuous

    // Bounds should be preserved.
    try testing.expectEqual(@as(f64, 0.0), lm.col_lower[0]);
    try testing.expectEqual(@as(f64, 1.0), lm.col_upper[0]); // binary bounds stay [0,1]
}

test "compileModel Sense-to-bounds mapping" {
    var env = try Env.initSimple(testing.allocator);
    defer env.deinit();
    var model = try Model.init(testing.allocator, &env, "sense_test");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, "le");
    try model.addConstr(0, &.{}, &.{}, .equal, 20.0, "eq");
    try model.addConstr(0, &.{}, &.{}, .greater_equal, 30.0, "ge");
    try model.addVar(0, &.{}, &.{}, 1.0, 0.0, 1.0, .continuous, null);
    try model.updateModel();

    var cm = try compileModel(testing.allocator, &model);
    defer cm.deinit();

    const lm = cm.linearData();
    try testing.expectEqual(@as(usize, 3), lm.num_rows);
    // less_equal → lower=-inf, upper=10
    try testing.expectEqual(-INFINITY, lm.row_lower[0]);
    try testing.expectEqual(@as(f64, 10.0), lm.row_upper[0]);
    // equal → lower=20, upper=20
    try testing.expectEqual(@as(f64, 20.0), lm.row_lower[1]);
    try testing.expectEqual(@as(f64, 20.0), lm.row_upper[1]);
    // greater_equal → lower=30, upper=inf
    try testing.expectEqual(@as(f64, 30.0), lm.row_lower[2]);
    try testing.expectEqual(INFINITY, lm.row_upper[2]);
}

test "compileModel with SOS returns FeatureNotAvailable" {
    var env = try Env.initSimple(testing.allocator);
    defer env.deinit();
    var model = try Model.init(testing.allocator, &env, "sos_test");
    defer model.deinit();

    try model.addVar(0, &.{}, &.{}, 1.0, 0.0, 1.0, .continuous, null);
    try model.addSOS(.sos1, 1, &[_]usize{0}, null, "sos1");

    try testing.expectError(error.FeatureNotAvailable, compileModel(testing.allocator, &model));
}

test "compileModel with general constraint returns FeatureNotAvailable" {
    var env = try Env.initSimple(testing.allocator);
    defer env.deinit();
    var model = try Model.init(testing.allocator, &env, "gencon_test");
    defer model.deinit();

    try model.addVar(0, &.{}, &.{}, 1.0, 0.0, 1.0, .continuous, null);
    try model.addGenConstrMax(0, 1, &[_]usize{0}, 0.0, "max1");

    try testing.expectError(error.FeatureNotAvailable, compileModel(testing.allocator, &model));
}

test "compileModel source Model unchanged after compile" {
    var env = try Env.initSimple(testing.allocator);
    defer env.deinit();
    var model = try Model.init(testing.allocator, &env, "check");
    defer model.deinit();

    try model.addConstr(1, &[_]usize{0}, &[_]f64{2.0}, .less_equal, 10.0, "c0");
    try model.addVar(1, &[_]usize{0}, &[_]f64{2.0}, 1.0, 0.0, 5.0, .continuous, null);
    try model.updateModel();

    const orig_nv = model.num_vars;
    const orig_nc = model.num_constrs;
    const orig_cost = model.var_obj[0];

    var cm = try compileModel(testing.allocator, &model);
    defer cm.deinit();

    try testing.expectEqual(orig_nv, model.num_vars);
    try testing.expectEqual(orig_nc, model.num_constrs);
    try testing.expectEqual(orig_cost, model.var_obj[0]);
}

test "compileModel QP with Hessian" {
    var env = try Env.initSimple(testing.allocator);
    defer env.deinit();
    var model = try Model.init(testing.allocator, &env, "qp");
    defer model.deinit();

    // 2 variables, with a 2x2 Q matrix: Q = [[2, 0], [0, 3]]
    // Lower triangle: (0,0)=2, (1,1)=3
    try model.addQPterms(&[_]i32{ 0, 1 }, &[_]i32{ 0, 1 }, &[_]f64{ 2.0, 3.0 });
    try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, "c0");
    try model.addVar(0, &.{}, &.{}, 1.0, 0.0, 10.0, .continuous, null);
    try model.addVar(0, &.{}, &.{}, 2.0, 0.0, 10.0, .continuous, null);
    try model.updateModel();

    var cm = try compileModel(testing.allocator, &model);
    defer cm.deinit();

    try testing.expectEqual(ProblemClass.qp, cm.problemClass());
}

test "compileModel QCP with one quadratic constraint" {
    var env = try Env.initSimple(testing.allocator);
    defer env.deinit();
    var model = try Model.init(testing.allocator, &env, "qcp");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, "c0");
    try model.addVar(0, &.{}, &.{}, 1.0, 0.0, 5.0, .continuous, null);
    try model.addVar(0, &.{}, &.{}, 2.0, 0.0, 5.0, .continuous, null);
    try model.addQConstr(
        2, // qnz
        &[_]i32{ 0, 1 }, // qrow
        &[_]i32{ 0, 1 }, // qcol
        &[_]f64{ 1.0, 2.0 }, // qval
        1, // lnz
        &[_]usize{0}, // lind
        &[_]f64{1.0}, // lval
        .less_equal,
        5.0,
        "qc0",
    );
    try model.updateModel();

    var cm = try compileModel(testing.allocator, &model);
    defer cm.deinit();

    try testing.expectEqual(ProblemClass.qcp, cm.problemClass());
}

test "compileModel matrix data is copied correctly" {
    var env = try Env.initSimple(testing.allocator);
    defer env.deinit();
    var model = try Model.init(testing.allocator, &env, "matrix_cp");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 100.0, "c0");
    try model.addConstr(0, &.{}, &.{}, .equal, 50.0, "c1");
    try model.addVar(2, &[_]usize{ 0, 1 }, &[_]f64{ 1.0, 2.0 }, 1.0, 0.0, 10.0, .continuous, null);
    try model.addVar(2, &[_]usize{ 0, 1 }, &[_]f64{ 3.0, 4.0 }, 2.0, 0.0, 20.0, .continuous, null);
    try model.updateModel();

    var cm = try compileModel(testing.allocator, &model);
    defer cm.deinit();

    const src = model.matrix.csc();
    const dst = cm.linearData().matrix.csc();

    try testing.expectEqual(src.num_rows, dst.num_rows);
    try testing.expectEqual(src.num_cols, dst.num_cols);
    try testing.expectEqual(src.nnz(), dst.nnz());

    // Check values match.
    for (0..src.nnz()) |i| {
        try testing.expectEqual(src.values[i], dst.values[i]);
    }

    // Check they are independent copies.
    try testing.expect(src.col_starts.ptr != dst.col_starts.ptr);
}

test "compileModel with all continuous vars yields null integrality" {
    var env = try Env.initSimple(testing.allocator);
    defer env.deinit();
    var model = try Model.init(testing.allocator, &env, "continuous");
    defer model.deinit();

    try model.addVar(0, &.{}, &.{}, 1.0, 0.0, 1.0, .continuous, null);
    try model.addVar(0, &.{}, &.{}, 2.0, 0.0, 1.0, .continuous, null);
    try model.addConstr(0, &.{}, &.{}, .less_equal, 5.0, "c0");
    try model.updateModel();

    var cm = try compileModel(testing.allocator, &model);
    defer cm.deinit();

    try testing.expect(cm.linearData().integrality == null);
}
