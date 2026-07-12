//! Temporary constraint — expression paired with sense and RHS.

const std = @import("std");
const Sense = @import("../types.zig").Sense;
const LinExpr = @import("lin_expr.zig").LinExpr;
const QuadExpr = @import("quad_expr.zig").QuadExpr;

/// A temporary constraint produced by comparing an expression to a value.
/// Designed for use with `Model.addConstr`:
///
/// ```zig
/// const tc = expr <= 10.0;
/// try model.addConstr(tc, "my_constr");
/// ```
pub const TempConstr = struct {
    /// The linear part of the constraint.
    lin_expr: ?LinExpr = null,
    /// The quadratic part (for quadratic constraints).
    quad_expr: ?QuadExpr = null,
    /// Constraint sense.
    sense: Sense,
    /// Right-hand side value.
    rhs: f64,

    const Self = @This();

    /// Create a temporary constraint from a linear term.
    pub fn initLin(expr: LinExpr, sense: Sense, rhs: f64) Self {
        return .{ .lin_expr = expr, .sense = sense, .rhs = rhs };
    }

    /// Create a temporary constraint from a quadratic term.
    pub fn initQuad(expr: QuadExpr, sense: Sense, rhs: f64) Self {
        return .{ .quad_expr = expr, .sense = sense, .rhs = rhs };
    }
};

// ── Operator overloads ─────────────────────────────────────────────────────

/// Compare a `LinExpr` with a scalar: `expr <= rhs`.
pub fn leExpr(expr: *LinExpr, rhs: f64) TempConstr {
    return TempConstr.initLin(expr.*, .less_equal, rhs);
}

/// Compare a `LinExpr` with a scalar: `expr == rhs`.
pub fn eqExpr(expr: *LinExpr, rhs: f64) TempConstr {
    return TempConstr.initLin(expr.*, .equal, rhs);
}

/// Compare a `LinExpr` with a scalar: `expr >= rhs`.
pub fn geExpr(expr: *LinExpr, rhs: f64) TempConstr {
    return TempConstr.initLin(expr.*, .greater_equal, rhs);
}
