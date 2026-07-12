//! Optimisation model module.
//!
//! This module provides a lazy-update, attribute-based model layer.
//! The central types are:
//!
//! - [`types`] — enums and constants (`VarType`, `Sense`, `Status`, …)
//! - [`Env`] — environment (parameter container, logging)
//! - [`Model`] — problem data, lazy updates, attribute-based access
//! - [`attrs`] — attribute name constants and metadata
//! - Entity sub-modules (`var`, `constraint`, `sos`, `qconstr`, `genconstr`)

const std = @import("std");

pub const types = @import("types.zig");
pub const Env = @import("env.zig").Env;
pub const Model = @import("model.zig").Model;
pub const attrs = @import("attrs.zig");
pub const expr = @import("expr/root.zig");

// Entity-type sub-modules.
pub const var_ = @import("var/root.zig");
pub const constraint = @import("constraint/root.zig");
pub const sos = @import("sos/root.zig");
pub const qconstr = @import("qconstr/root.zig");
pub const genconstr = @import("genconstr/root.zig");

// Re-export the public type aliases for convenience.
pub const VarType = types.VarType;
pub const Sense = types.Sense;
pub const ObjectiveSense = types.ObjectiveSense;
pub const Status = types.Status;
pub const BasisStatus = types.BasisStatus;
pub const SosType = types.SosType;
pub const GenConstrType = types.GenConstrType;
pub const FeasRelaxType = types.FeasRelaxType;
pub const CallbackFunc = types.CallbackFunc;
pub const CallbackWhere = types.CallbackWhere;
pub const ModelError = types.ModelError;
pub const INFINITY = types.INFINITY;
pub const EPSILON = types.EPSILON;

// Re-export object wrapper types (from entity sub-modules).
pub const Var = var_.Var;
pub const VarData = var_.VarData;
pub const Constr = constraint.Constr;
pub const ConstrData = constraint.ConstrData;
pub const QConstr = qconstr.QConstr;
pub const SOS = sos.SOS;
pub const GenConstr = genconstr.GenConstr;

// Re-export expression / column types.
pub const Column = expr.Column;
pub const LinExpr = expr.LinExpr;
pub const QuadExpr = expr.QuadExpr;
pub const TempConstr = expr.TempConstr;

test {
    std.testing.refAllDecls(@This());
}
