//! Stable Zig-facing API for the zhighs optimisation library.
//!
//! This module exposes solver construction, model mutation, options, status,
//! callbacks, and result access without leaking internal details.
//!
//! Usage pattern:
//!
//! ```zig
//! const zhighs = @import("zhighs");
//! const api = zhighs.api;
//!
//! // 1. Create environment
//! var env = try api.Env.init(allocator, "solver.log");
//! defer env.deinit();
//!
//! // 2. Create an empty model
//! var model = try api.Model.init(allocator, &env, "my_problem");
//! defer model.deinit();
//!
//! // 3. Add variables (CSC column-format)
//! try model.addVars(2, &[_]usize{0, 0}, &[_]usize{}, &[_]f64{},
//!                   &[_]f64{1.0, 2.0}, &[_]f64{0.0, 0.0},
//!                   &[_]f64{INFINITY, INFINITY}, null, null);
//!
//! // 4. Add constraints (CSR row-format)
//! try model.addConstr(2, &[_]usize{0, 1}, &[_]f64{1.0, 1.0},
//!                     .less_equal, 10.0, "c0");
//!
//! // 5. Optimize
//! try model.optimize();
//!
//! // 6. Query results via attributes
//! const obj_val = try model.getDblAttr("ObjVal");
//! const x0 = try model.getDblAttrElement("X", 0);
//! ```

const std = @import("std");
const model_root = @import("model");

pub const Env = model_root.Env;
pub const Model = model_root.Model;

pub const VarType = model_root.VarType;
pub const Sense = model_root.Sense;
pub const ObjectiveSense = model_root.ObjectiveSense;
pub const Status = model_root.Status;
pub const BasisStatus = model_root.BasisStatus;
pub const SosType = model_root.SosType;
pub const GenConstrType = model_root.GenConstrType;
pub const FeasRelaxType = model_root.FeasRelaxType;
pub const CallbackFunc = model_root.CallbackFunc;
pub const CallbackWhere = model_root.CallbackWhere;
pub const ModelError = model_root.ModelError;
pub const INFINITY = model_root.INFINITY;
pub const EPSILON = model_root.EPSILON;

// Object wrapper types.
pub const Var = model_root.Var;
pub const Constr = model_root.Constr;
pub const QConstr = model_root.QConstr;
pub const SOS = model_root.SOS;
pub const GenConstr = model_root.GenConstr;

// Expression / column types.
pub const Column = model_root.Column;
pub const LinExpr = model_root.LinExpr;
pub const QuadExpr = model_root.QuadExpr;
pub const TempConstr = model_root.TempConstr;

// Attribute name constants for use with get*Attr / set*Attr.
pub const attrs = model_root.attrs;

// Parameter name constants for use with env.set*Param / env.get*Param.
pub const params = struct {
    pub const METHOD = "Method";
    pub const THREADS = "Threads";
    pub const PRESOLVE = "Presolve";
    pub const TIMELIMIT = "TimeLimit";
    pub const FEASIBILITY_TOL = "FeasibilityTol";
    pub const OPTIMALITY_TOL = "OptimalityTol";
    pub const OUTPUT_FLAG = "OutputFlag";
    pub const MIP_GAP = "MIPGap";
    pub const MIP_GAP_ABS = "MIPGapAbs";
    pub const ITERATION_LIMIT = "IterationLimit";
    pub const NODE_LIMIT = "NodeLimit";
    pub const SOLUTION_LIMIT = "SolutionLimit";
    pub const BAR_ITER_LIMIT = "BarIterLimit";
    pub const CROSSOVER = "Crossover";
    pub const BAR_HOMOGENEOUS = "BarHomogeneous";
    pub const DUAL_REDUCTIONS = "DualReductions";
    pub const INF_UNBD_INFO = "InfUnbdInfo";
    pub const NORM_ADJUST = "NormAdjust";
    pub const LP_WARM_START = "LPWarmStart";
    pub const OBJ_SCALE = "ObjScale";
    pub const BAR_QCP_CONV_TOL = "BarQCPConvTol";
    pub const MARKOWITZ_TOL = "MarkowitzTol";
    pub const PSD_TOL = "PSDTol";
    pub const DISPLAY_INTERVAL = "DisplayInterval";
    pub const LOG_TO_CONSOLE = "LogToConsole";
    pub const QCP_DUAL = "QCPDual";
};

test {
    std.testing.refAllDecls(@This());
}
