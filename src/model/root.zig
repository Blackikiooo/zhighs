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
pub const entity_handle = @import("entity_handle.zig");
pub const VarId = entity_handle.VarId;
pub const ConstrId = entity_handle.ConstrId;
pub const QConstrId = entity_handle.QConstrId;
pub const SosId = entity_handle.SosId;
pub const GenConstrId = entity_handle.GenConstrId;
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
pub const Constr = constraint.Constr;
pub const QConstr = qconstr.QConstr;
pub const SOS = sos.SOS;
pub const GenConstr = genconstr.GenConstr;

// Re-export expression / column types.
pub const Column = expr.Column;
pub const LinExpr = expr.LinExpr;
pub const QuadExpr = expr.QuadExpr;
pub const TempConstr = expr.TempConstr;

// ── Solver-internal core model types (Experimental API) ──────────────────
//
// These types form the solver-internal mathematical IR layer.  They sit
// between the Gurobi-style user `Model` and the presolve / simplex layers.
//
// All APIs here are Experimental until explicitly marked Stable API.

pub const problem_class = @import("problem_class.zig");
pub const linear_model = @import("linear_model.zig");
pub const linear_model_builder = @import("linear_model_builder.zig");
pub const hessian = @import("hessian.zig");
pub const quadratic_model = @import("quadratic_model.zig");
pub const expression_graph = @import("expression_graph.zig");
pub const nonlinear_model = @import("nonlinear_model.zig");
pub const compiled_model = @import("compiled_model.zig");
pub const compiled_model_view = @import("compiled_model_view.zig");
pub const revisions = @import("revisions.zig");
pub const validate = @import("validate.zig");
pub const solution = @import("solution.zig");
pub const compile_model_module = @import("compile_model.zig");
pub const residual = @import("residual.zig");

// Re-export core types at the module top level for discovery.
pub const Integrality = linear_model.Integrality;
pub const LinearModel = linear_model.LinearModel;
pub const LinearModelBuilder = linear_model_builder.LinearModelBuilder;
pub const ProblemClass = problem_class.ProblemClass;
pub const DomainClass = problem_class.DomainClass;
pub const ObjectiveClass = problem_class.ObjectiveClass;
pub const ConstraintClass = problem_class.ConstraintClass;
pub const SolverCapability = problem_class.SolverCapability;
pub const classify = problem_class.classify;
pub const Hessian = hessian.Hessian;
pub const HessianFormat = hessian.HessianFormat;
pub const Curvature = hessian.Curvature;
pub const QuadraticConstraint = quadratic_model.QuadraticConstraint;
pub const QuadraticModel = quadratic_model.QuadraticModel;
pub const CompiledModel = compiled_model.CompiledModel;
pub const CompiledLinearModelView = compiled_model_view.CompiledLinearModelView;
pub const ExpressionGraph = expression_graph.ExpressionGraph;
pub const ExpressionGraphBuilder = expression_graph.ExpressionGraphBuilder;
pub const NodeId = expression_graph.NodeId;
pub const Opcode = expression_graph.Opcode;
pub const evaluate = expression_graph.evaluate;
pub const NonlinearConstraint = nonlinear_model.NonlinearConstraint;
pub const NonlinearModel = nonlinear_model.NonlinearModel;
pub const ValidationError = validate.ValidationError;
pub const validateLinearModel = validate.validateLinearModel;
pub const validateHessian = validate.validateHessian;
pub const validateQuadraticConstraint = validate.validateQuadraticConstraint;
pub const validateQuadraticModel = validate.validateQuadraticModel;
pub const validateExpressionGraph = validate.validateExpressionGraph;
pub const validateNonlinearModel = validate.validateNonlinearModel;
pub const validateCompiledModel = validate.validateCompiledModel;

// Solution / Basis / SolveInfo
pub const Solution = solution.Solution;
pub const Basis = solution.Basis;
pub const SolveInfo = solution.SolveInfo;

// Compile model
pub const CompileError = compile_model_module.CompileError;
pub const compileModel = compile_model_module.compileModel;

// Residual / KKT
pub const Residual = residual.Residual;
pub const KKTStatus = residual.KKTStatus;
pub const computePrimalResidual = residual.computePrimalResidual;
pub const computeDualResidual = residual.computeDualResidual;
pub const checkKKT = residual.checkKKT;

test {
    std.testing.refAllDecls(@This());
}
