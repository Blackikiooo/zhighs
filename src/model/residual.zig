//! Primal–dual residual computation and KKT verification.
//!
//! These functions check how well a candidate [`Solution`] satisfies the
//! optimality conditions for a given [`LinearModel`].  They are intended for
//! solver debugging, solution validation, and convergence monitoring — not
//! for performance‑sensitive inner loops.
//!
//! ## KKT conditions (LP)
//!
//! 1. **Primal feasibility**
//!    - `col_lower ≤ x ≤ col_upper`
//!    - `row_lower ≤ Ax ≤ row_upper`
//!
//! 2. **Dual feasibility**
//!    - `rc = c − Aᵀπ`  (definition)
//!    - Reduced cost sign conditions per bound type (minimise convention)
//!
//! 3. **Complementarity**
//!    - `(x − lb) · rc ≥ −tol`
//!    - `(ub − x) · (−rc) ≥ −tol`
//!
//! Experimental API.

const std = @import("std");
const linear_model_module = @import("linear_model.zig");
const solution_module = @import("solution.zig");
const matrix = @import("matrix");

const LinearModel = linear_model_module.LinearModel;
const Solution = solution_module.Solution;
const INFINITY = @import("types.zig").INFINITY;
const EPSILON = @import("types.zig").EPSILON;

// ── Tolerance ────────────────────────────────────────────────────────────────

/// Default tolerance for KKT checks.
pub const DEFAULT_TOL: f64 = 1e-8;

// ── Residual ────────────────────────────────────────────────────────────────

/// Summary of primal‑dual infeasibility measures.
pub const Residual = struct {
    /// Maximum constraint violation: max(0, lb − Ax, Ax − ub).
    max_primal_infeasibility: f64,
    /// Maximum bound violation: max(0, lb − x, x − ub).
    max_bound_violation: f64,
    /// Maximum |rc − (c − Aᵀπ)|.
    max_dual_infeasibility: f64,
    /// Maximum complementarity product violation.
    max_complementarity: f64,
    /// Relative objective gap |primal − dual| / max(1, |dual|).
    objective_gap: f64,
};

// ── KKTStatus ────────────────────────────────────────────────────────────────

/// Outcome of a full KKT check.
pub const KKTStatus = enum {
    optimal,
    primal_infeasible,
    dual_infeasible,
    not_optimal,
    undetermined,
};

// ── computePrimalResidual ────────────────────────────────────────────────────

/// Compute primal infeasibility measures for a candidate solution.
///
/// `work` must have length ≥ `model.num_rows` (used as a scratch buffer for
/// the matrix‑vector product).  Pass `null` to skip the row‑constraint check.
pub fn computePrimalResidual(
    model: *const LinearModel,
    sol: *const Solution,
    work: ?[]f64,
) Residual {
    const num_rows = model.num_rows;
    const num_cols = model.num_cols;
    const csc = model.matrix.csc();

    // ── Bound violations ───────────────────────────────────────────────
    var max_bound_viol: f64 = 0.0;
    const safe_ncols = @min(sol.primal.len, num_cols);
    for (0..safe_ncols) |j| {
        const x = sol.primal[j];
        const lb = model.col_lower[j];
        const ub = model.col_upper[j];
        if (lb != -INFINITY and x < lb - DEFAULT_TOL) {
            max_bound_viol = @max(max_bound_viol, lb - x);
        }
        if (ub != INFINITY and x > ub + DEFAULT_TOL) {
            max_bound_viol = @max(max_bound_viol, x - ub);
        }
    }

    // ── Row constraint violations ──────────────────────────────────────
    var max_row_viol: f64 = 0.0;
    if (work) |ax| {
        if (num_rows > 0 and ax.len >= num_rows) {
            @memset(ax[0..num_rows], 0.0);
            if (num_cols > 0 and sol.primal.len > 0) {
                matrix.addProduct(csc.*, 1.0, sol.primal, ax[0..num_rows]) catch {};
            }
            for (0..num_rows) |i| {
                const val = ax[i];
                const lb = model.row_lower[i];
                const ub = model.row_upper[i];
                if (lb != -INFINITY and val < lb - DEFAULT_TOL) {
                    max_row_viol = @max(max_row_viol, lb - val);
                }
                if (ub != INFINITY and val > ub + DEFAULT_TOL) {
                    max_row_viol = @max(max_row_viol, val - ub);
                }
            }
        }
    }

    return Residual{
        .max_primal_infeasibility = @max(max_row_viol, 0.0),
        .max_bound_violation = @max(max_bound_viol, 0.0),
        .max_dual_infeasibility = 0.0,
        .max_complementarity = 0.0,
        .objective_gap = 0.0,
    };
}

// ── computeDualResidual ──────────────────────────────────────────────────────

/// Compute dual infeasibility and complementarity measures.
///
/// `work` must have length ≥ `model.num_cols` (used as a scratch buffer for
/// reduced cost computation).  Pass `null` to skip.
pub fn computeDualResidual(
    model: *const LinearModel,
    sol: *const Solution,
    work: ?[]f64,
) Residual {
    const num_rows = model.num_rows;
    const num_cols = model.num_cols;
    const csc = model.matrix.csc();

    var max_dual_inf: f64 = 0.0;
    var max_comp: f64 = 0.0;

    if (work) |rc| {
        if (rc.len >= num_cols) {
            // rc = c − Aᵀπ
            @memcpy(rc[0..num_cols], model.col_cost);
            if (num_rows > 0 and num_cols > 0 and sol.dual.len > 0) {
                matrix.addTransposeProduct(csc.*, -1.0, sol.dual, rc[0..num_cols]) catch {};
            }

            // Dual infeasibility: max |c − Aᵀπ − rc|
            const safe_ncols = @min(sol.reduced_cost.len, num_cols);
            for (0..safe_ncols) |j| {
                const dev = @abs(rc[j] - sol.reduced_cost[j]);
                if (dev > max_dual_inf) max_dual_inf = dev;
            }

            // Complementarity (minimise convention).
            const primal_end = @min(sol.primal.len, num_cols);
            for (0..primal_end) |j| {
                const x = sol.primal[j];
                const lb = model.col_lower[j];
                const ub = model.col_upper[j];
                const rc_j = rc[j];

                // (x − lb)·rc ≥ −tol
                if (lb != -INFINITY) {
                    const prod = (x - lb) * rc_j;
                    if (prod < -DEFAULT_TOL) {
                        max_comp = @max(max_comp, -prod);
                    }
                }
                // (ub − x)·(−rc) ≥ −tol
                if (ub != INFINITY) {
                    const prod = (ub - x) * (-rc_j);
                    if (prod < -DEFAULT_TOL) {
                        max_comp = @max(max_comp, -prod);
                    }
                }
            }
        }
    }

    // Objective gap (minimise: dual = bᵀπ for standard form).
    var obj_gap: f64 = 0.0;
    {
        const primal_obj = if (sol.primal.len >= num_cols) blk: {
            var v = model.objective_offset;
            for (0..num_cols) |j| v += model.col_cost[j] * sol.primal[j];
            break :blk v;
        } else 0.0;
        // Simplified dual objective using row_upper as b (minimise).
        const dual_obj = if (sol.dual.len >= num_rows) blk: {
            var v = model.objective_offset;
            for (0..num_rows) |i| {
                v += sol.dual[i] * model.row_upper[i];
            }
            break :blk v;
        } else 0.0;
        const denom = @max(1.0, @abs(dual_obj));
        obj_gap = @abs(primal_obj - dual_obj) / denom;
    }

    return Residual{
        .max_primal_infeasibility = 0.0,
        .max_bound_violation = 0.0,
        .max_dual_infeasibility = max_dual_inf,
        .max_complementarity = max_comp,
        .objective_gap = obj_gap,
    };
}

// ── checkKKT ─────────────────────────────────────────────────────────────────

/// Perform a full KKT check.
///
/// `work` must be at least `max(num_rows, num_cols)` elements long for
/// scratch space, or `null` to skip matrix‑vector products.
pub fn checkKKT(
    model: *const LinearModel,
    sol: *const Solution,
    tol: f64,
    work: ?[]f64,
) KKTStatus {
    const pres = computePrimalResidual(model, sol, work);
    if (pres.max_primal_infeasibility > tol or pres.max_bound_violation > tol) {
        return .primal_infeasible;
    }

    const dres = computeDualResidual(model, sol, work);
    if (dres.max_dual_infeasibility > tol or dres.max_complementarity > tol) {
        return .dual_infeasible;
    }

    if (pres.max_primal_infeasibility <= tol and
        dres.max_dual_infeasibility <= tol and
        dres.max_complementarity <= tol and
        dres.objective_gap <= tol)
    {
        return .optimal;
    }

    return .not_optimal;
}

// ── Tests ────────────────────────────────────────────────────────────────────

const testing = std.testing;

test "computePrimalResidual zero on zero model" {
    var lm = try LinearModel.initEmpty(testing.allocator);
    defer lm.deinit();
    var sol = try Solution.init(testing.allocator, 0, 0);
    defer sol.deinit();

    const res = computePrimalResidual(&lm, &sol, null);
    try testing.expectEqual(@as(f64, 0.0), res.max_primal_infeasibility);
    try testing.expectEqual(@as(f64, 0.0), res.max_bound_violation);
}

test "computePrimalResidual detects bound violations" {
    const Builder = @import("linear_model_builder.zig").LinearModelBuilder;
    var bld = try Builder.init(testing.allocator, 1, 1);
    defer bld.deinit();
    try bld.setColBounds(0, 0.0, 1.0);
    try bld.appendCoefficient(0, 0, 1.0);
    var lm = try bld.freeze();
    defer lm.deinit();

    var sol = try Solution.init(testing.allocator, 1, 1);
    defer sol.deinit();

    // No violation.
    sol.primal[0] = 0.5;
    var res = computePrimalResidual(&lm, &sol, null);
    try testing.expectEqual(@as(f64, 0.0), res.max_bound_violation);

    // Upper bound violation.
    sol.primal[0] = 2.0;
    res = computePrimalResidual(&lm, &sol, null);
    try testing.expect(res.max_bound_violation > 0.0);
}

test "computePrimalResidual detects row constraint violation" {
    const Builder = @import("linear_model_builder.zig").LinearModelBuilder;
    var bld = try Builder.init(testing.allocator, 1, 1);
    defer bld.deinit();
    try bld.setRowBounds(0, 0.0, 5.0);
    try bld.appendCoefficient(0, 0, 1.0);
    var lm = try bld.freeze();
    defer lm.deinit();

    var sol = try Solution.init(testing.allocator, 1, 1);
    defer sol.deinit();

    var work: [1]f64 = undefined;
    sol.primal[0] = 10.0; // Ax = 10, violates row_upper=5
    const res = computePrimalResidual(&lm, &sol, &work);
    try testing.expect(res.max_primal_infeasibility > 0.0);
}

test "computeDualResidual zero on empty model" {
    var lm = try LinearModel.initEmpty(testing.allocator);
    defer lm.deinit();
    var sol = try Solution.init(testing.allocator, 0, 0);
    defer sol.deinit();

    const res = computeDualResidual(&lm, &sol, null);
    try testing.expectEqual(@as(f64, 0.0), res.max_dual_infeasibility);
    try testing.expectEqual(@as(f64, 0.0), res.max_complementarity);
}

test "checkKKT returns optimal for empty feasible model" {
    var lm = try LinearModel.initEmpty(testing.allocator);
    defer lm.deinit();
    var sol = try Solution.init(testing.allocator, 0, 0);
    defer sol.deinit();

    const kkt = checkKKT(&lm, &sol, DEFAULT_TOL, null);
    try testing.expectEqual(KKTStatus.optimal, kkt);
}

test "checkKKT returns primal_infeasible for bound violation" {
    const Builder = @import("linear_model_builder.zig").LinearModelBuilder;
    var bld = try Builder.init(testing.allocator, 1, 1);
    defer bld.deinit();
    try bld.setColBounds(0, 0.0, 1.0);
    try bld.appendCoefficient(0, 0, 1.0);
    var lm = try bld.freeze();
    defer lm.deinit();

    var sol = try Solution.init(testing.allocator, 1, 1);
    defer sol.deinit();
    sol.primal[0] = 5.0;

    const kkt = checkKKT(&lm, &sol, DEFAULT_TOL, null);
    try testing.expectEqual(KKTStatus.primal_infeasible, kkt);
}
