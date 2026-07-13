//! Main optimisation model.
//!
//! A `Model` owns the problem data — variables, linear constraints, quadratic
//! objective, SOS and general constraints — along with solution results and
//! status.  The design follows a lazy-update, attribute-based convention:
//!
//! 1. **Lazy updates** – modifications are queued and applied in batch by
//!    `updateModel`, `optimize`, or `writeModel`.
//! 2. **Attribute-based access** – problem and solution data are read/written
//!    through uniform `get*Attr` / `set*Attr` methods keyed by string names
//!    (e.g. `"LB"`, `"Obj"`, `"X"`, `"Status"`).
//! 3. **Index-based references** – variables and constraints are identified by
//!    contiguous `usize` indices.
//! 4. **Incremental construction** – call `addVar`/`addConstr` one at a time
//!    or use the batch `addVars`/`addConstrs` counterparts.
//!
//! Methods are split across sub-files for maintainability and re-exported
//! via `pub const` so they remain callable as `model.method(...)`.
//!
//! ## Responsibility
//!
//! This file defines the `Model` data layout, lifetime management, basic size
//! queries, and the public method surface.  Domain behaviour belongs in the
//! corresponding `model_*.zig` implementation file; new operations should not
//! be implemented here unless they directly manage the model's core lifetime
//! or representation.

const std = @import("std");
const types = @import("types.zig");
const env_module = @import("env.zig");
const Attr = @import("attrs.zig").Attr;

// Foundation types (strongly-typed indices, etc.)
const foundation = @import("foundation");
const RowId = foundation.RowId;

// Matrix infrastructure
const matrix = @import("matrix");
const MatrixStore = matrix.MatrixStore;
const CscMatrix = matrix.CscMatrix;
const MatrixBuilder = matrix.MatrixBuilder;

const Env = env_module.Env;
const VarType = types.VarType;
const Sense = types.Sense;
const ObjectiveSense = types.ObjectiveSense;
const Status = types.Status;
const BasisStatus = types.BasisStatus;
const SosType = types.SosType;
const GenConstrType = types.GenConstrType;
const FeasRelaxType = types.FeasRelaxType;
const CallbackWhere = types.CallbackWhere;
const CallbackFunc = types.CallbackFunc;
const ModelError = types.ModelError;
const INFINITY = types.INFINITY;

// ── Model ───────────────────────────────────────────────────────────────

pub const Model = struct {
    allocator: std.mem.Allocator,
    env: *Env,

    // ── Identity ───────────────────────────────────────────────────────
    name: []const u8,

    // ── Variable data (committed) ──────────────────────────────────────
    num_vars: usize = 0,
    var_lb: []f64 = &.{},
    var_ub: []f64 = &.{},
    var_obj: []f64 = &.{},
    var_type: []VarType = &.{},
    var_names: [](?[]const u8) = &.{},

    // ── Constraint data (committed) ────────────────────────────────────
    num_constrs: usize = 0,
    constr_sense: []Sense = &.{},
    constr_rhs: []f64 = &.{},
    constr_names: [](?[]const u8) = &.{},

    // Constraint matrix in canonical CSC storage.
    matrix: MatrixStore = undefined,

    // ── Quadratic objective terms (lower-triangle triples) ─────────────
    q_nz: usize = 0,
    q_row: []i32 = &.{},
    q_col: []i32 = &.{},
    q_val: []f64 = &.{},

    // ── Model state ────────────────────────────────────────────────────
    sense: ObjectiveSense = .minimize,
    status: Status = .loaded,
    revision: u64 = 0,

    // ── Solutions (populated after `optimize`) ─────────────────────────
    solution: []f64 = &.{}, // X
    reduced_cost: []f64 = &.{}, // RC
    slack: []f64 = &.{}, // Slack
    pi: []f64 = &.{}, // Pi (dual)

    obj_val: f64 = 0.0,
    obj_bound: f64 = 0.0,
    obj_con: f64 = 0.0,
    iter_count: i64 = 0,
    node_count: i64 = 0,
    bar_iter_count: i64 = 0,

    // ── Basis ──────────────────────────────────────────────────────────
    vbasis: []BasisStatus = &.{},
    cbasis: []BasisStatus = &.{},

    // ── Start vectors (warm-start / MIP start) ─────────────────────────
    mip_start: []f64 = &.{}, // MIP start (Start attribute)
    p_start: []f64 = &.{}, // primal start (PStart)
    d_start: []f64 = &.{}, // dual start (DStart)

    // ── Pending changes (lazy update) ──────────────────────────────────
    pending: std.ArrayListUnmanaged(PendingChange) = .empty,
    has_pending: bool = false,

    // ── SOS constraints ────────────────────────────────────────────────
    sos_count: usize = 0,
    sos_types: []SosType = &.{},
    sos_begin: []usize = &.{}, // start index into sos_indices / sos_weights
    sos_indices: []usize = &.{},
    sos_weights: []f64 = &.{},
    sos_names: [](?[]const u8) = &.{},

    // ── Quadratic constraints ──────────────────────────────────────────
    qconstr_count: usize = 0,
    qconstr_qrow: []i32 = &.{},
    qconstr_qcol: []i32 = &.{},
    qconstr_qval: []f64 = &.{},
    qconstr_lind: []usize = &.{},
    qconstr_lval: []f64 = &.{},
    qconstr_sense: []Sense = &.{},
    qconstr_rhs: []f64 = &.{},
    qconstr_names: [](?[]const u8) = &.{},

    // ── General constraints ────────────────────────────────────────────
    genconstr_count: usize = 0,
    genconstr_types: []GenConstrType = &.{},
    genconstr_resvar: []usize = &.{},
    genconstr_nvars: []usize = &.{}, // number of vars per constraint
    genconstr_indices: []usize = &.{}, // packed variable indices
    genconstr_names: [](?[]const u8) = &.{},

    // ── Piecewise-linear objective data ──────────────────────────────────
    pwlobj_count: usize = 0,
    pwlobj_var: []usize = &.{}, // variable index for each PWL objective
    pwlobj_npts: []usize = &.{}, // number of points per entry
    pwlobj_xdata: []f64 = &.{}, // packed x values
    pwlobj_ydata: []f64 = &.{}, // packed y values

    // ── Callback / interrupt state ─────────────────────────────────────
    interrupted: bool = false,

    const Self = @This();

    // ── Construction / destruction ─────────────────────────────────────

    /// Create a new empty model within the given environment.
    pub fn init(allocator: std.mem.Allocator, env: *Env, name: []const u8) (ModelError)!Self {
        var self = Self{
            .allocator = allocator,
            .env = env,
            .name = "",
            .matrix = undefined,
        };
        self.name = try allocator.dupe(u8, name);
        const empty = CscMatrix.initZero(allocator, 0, 0) catch return error.OutOfMemory;
        self.matrix = MatrixStore.initAssumeValid(empty);
        return self;
    }

    /// Release all dynamically allocated memory.
    pub fn deinit(self: *Self) void {
        const alloc = self.allocator;

        // Flush (but don't hold onto) any pending changes.
        self.discardPending();

        alloc.free(self.name);
        alloc.free(self.var_lb);
        alloc.free(self.var_ub);
        alloc.free(self.var_obj);
        alloc.free(self.var_type);
        for (self.var_names) |n| if (n) |s| alloc.free(s);
        alloc.free(self.var_names);
        alloc.free(self.constr_sense);
        alloc.free(self.constr_rhs);
        for (self.constr_names) |n| if (n) |s| alloc.free(s);
        alloc.free(self.constr_names);
        alloc.free(self.q_row);
        alloc.free(self.q_col);
        alloc.free(self.q_val);
        alloc.free(self.solution);
        alloc.free(self.reduced_cost);
        alloc.free(self.slack);
        alloc.free(self.pi);
        alloc.free(self.vbasis);
        alloc.free(self.cbasis);
        alloc.free(self.mip_start);
        alloc.free(self.p_start);
        alloc.free(self.d_start);
        // SOS
        alloc.free(self.sos_types);
        alloc.free(self.sos_begin);
        alloc.free(self.sos_indices);
        alloc.free(self.sos_weights);
        for (self.sos_names) |n| if (n) |s| alloc.free(s);
        alloc.free(self.sos_names);
        // Quadratic constraints
        alloc.free(self.qconstr_qrow);
        alloc.free(self.qconstr_qcol);
        alloc.free(self.qconstr_qval);
        alloc.free(self.qconstr_lind);
        alloc.free(self.qconstr_lval);
        alloc.free(self.qconstr_sense);
        alloc.free(self.qconstr_rhs);
        for (self.qconstr_names) |n| if (n) |s| alloc.free(s);
        alloc.free(self.qconstr_names);
        // General constraints
        alloc.free(self.genconstr_types);
        alloc.free(self.genconstr_resvar);
        alloc.free(self.genconstr_nvars);
        alloc.free(self.genconstr_indices);
        for (self.genconstr_names) |n| if (n) |s| alloc.free(s);
        alloc.free(self.genconstr_names);
        // PWL objective
        alloc.free(self.pwlobj_var);
        alloc.free(self.pwlobj_npts);
        alloc.free(self.pwlobj_xdata);
        alloc.free(self.pwlobj_ydata);
        self.matrix.deinit(alloc);
        self.* = undefined;
    }

    // ── Model size queries ─────────────────────────────────────────────

    pub fn numVars(self: Self) usize {
        return self.num_vars + if (self.has_pending) self.countPendingAddVar() else 0;
    }

    pub fn numConstrs(self: Self) usize {
        return self.num_constrs + if (self.has_pending) self.countPendingAddConstr() else 0;
    }

    pub fn numNz(self: Self) usize {
        return self.matrix.csc().nnz();
    }

    // ── Re-exported sub-module methods ─────────────────────────────────

    pub const PendingChange = @import("model_pending.zig").PendingChange;
    pub const enqueue = @import("model_pending.zig").enqueue;
    pub const discardPending = @import("model_pending.zig").discardPending;
    pub const countPendingAddVar = @import("model_pending.zig").countPendingAddVar;
    pub const countPendingAddConstr = @import("model_pending.zig").countPendingAddConstr;

    pub const updateModel = @import("model_update.zig").updateModel;
    pub const applyPending = @import("model_update.zig").applyPending;

    pub const optimize = @import("model_solve.zig").optimize;
    pub const reset = @import("model_solve.zig").reset;
    pub const copy = @import("model_solve.zig").copy;
    pub const isMip = @import("model_solve.zig").isMip;

    pub const getIntAttr = @import("model_attr.zig").getIntAttr;
    pub const setIntAttr = @import("model_attr.zig").setIntAttr;
    pub const getDblAttr = @import("model_attr.zig").getDblAttr;
    pub const setDblAttr = @import("model_attr.zig").setDblAttr;
    pub const getStrAttr = @import("model_attr.zig").getStrAttr;
    pub const setStrAttr = @import("model_attr.zig").setStrAttr;
    pub const getIntAttrElement = @import("model_attr.zig").getIntAttrElement;
    pub const setIntAttrElement = @import("model_attr.zig").setIntAttrElement;
    pub const getDblAttrElement = @import("model_attr.zig").getDblAttrElement;
    pub const setDblAttrElement = @import("model_attr.zig").setDblAttrElement;
    pub const getCharAttrElement = @import("model_attr.zig").getCharAttrElement;
    pub const setCharAttrElement = @import("model_attr.zig").setCharAttrElement;
    pub const getStrAttrElement = @import("model_attr.zig").getStrAttrElement;
    pub const setStrAttrElement = @import("model_attr.zig").setStrAttrElement;
    pub const setDblAttrArray = @import("model_attr.zig").setDblAttrArray;
    pub const setIntAttrArray = @import("model_attr.zig").setIntAttrArray;
    pub const setStrAttrArray = @import("model_attr.zig").setStrAttrArray;
    pub const setCharAttrArray = @import("model_attr.zig").setCharAttrArray;
    pub const getStrAttrArray = @import("model_attr.zig").getStrAttrArray;
    pub const getDblAttrArray = @import("model_attr.zig").getDblAttrArray;
    pub const getIntAttrArray = @import("model_attr.zig").getIntAttrArray;
    pub const getCharAttrArray = @import("model_attr.zig").getCharAttrArray;
    pub const setDblAttrList = @import("model_attr.zig").setDblAttrList;
    pub const setIntAttrList = @import("model_attr.zig").setIntAttrList;
    pub const setStrAttrList = @import("model_attr.zig").setStrAttrList;
    pub const setCharAttrList = @import("model_attr.zig").setCharAttrList;
    pub const getDblAttrList = @import("model_attr.zig").getDblAttrList;
    pub const getIntAttrList = @import("model_attr.zig").getIntAttrList;
    pub const getStrAttrList = @import("model_attr.zig").getStrAttrList;
    pub const getCharAttrList = @import("model_attr.zig").getCharAttrList;

    pub const addVar = @import("model_linear.zig").addVar;
    pub const addVars = @import("model_linear.zig").addVars;
    pub const addConstr = @import("model_linear.zig").addConstr;
    pub const addConstrs = @import("model_linear.zig").addConstrs;
    pub const addSOS = @import("model_constraints.zig").addSOS;
    pub const addQConstr = @import("model_constraints.zig").addQConstr;
    pub const addGenConstrMax = @import("model_genconstr.zig").addGenConstrMax;
    pub const addGenConstrMin = @import("model_genconstr.zig").addGenConstrMin;
    pub const addGenConstrAbs = @import("model_genconstr.zig").addGenConstrAbs;
    pub const addGenConstrAnd = @import("model_genconstr.zig").addGenConstrAnd;
    pub const addGenConstrOr = @import("model_genconstr.zig").addGenConstrOr;
    pub const addGenConstrIndicator = @import("model_genconstr.zig").addGenConstrIndicator;
    pub const addGenConstrPWL = @import("model_genconstr.zig").addGenConstrPWL;
    pub const addGenConstrPoly = @import("model_genconstr.zig").addGenConstrPoly;
    pub const addGenConstrExp = @import("model_genconstr.zig").addGenConstrExp;
    pub const addGenConstrExpA = @import("model_genconstr.zig").addGenConstrExpA;
    pub const addGenConstrLog = @import("model_genconstr.zig").addGenConstrLog;
    pub const addGenConstrLogA = @import("model_genconstr.zig").addGenConstrLogA;
    pub const addGenConstrPow = @import("model_genconstr.zig").addGenConstrPow;
    pub const addGenConstrSin = @import("model_genconstr.zig").addGenConstrSin;
    pub const addGenConstrCos = @import("model_genconstr.zig").addGenConstrCos;
    pub const addGenConstrTan = @import("model_genconstr.zig").addGenConstrTan;
    pub const addGenConstrLogistic = @import("model_genconstr.zig").addGenConstrLogistic;
    pub const addGenConstrNorm = @import("model_genconstr.zig").addGenConstrNorm;
    pub const addGenConstrNL = @import("model_genconstr.zig").addGenConstrNL;
    pub const getGenConstrMax = @import("model_genconstr.zig").getGenConstrMax;
    pub const getGenConstrMin = @import("model_genconstr.zig").getGenConstrMin;
    pub const getGenConstrAbs = @import("model_genconstr.zig").getGenConstrAbs;
    pub const getGenConstrAnd = @import("model_genconstr.zig").getGenConstrAnd;
    pub const getGenConstrOr = @import("model_genconstr.zig").getGenConstrOr;
    pub const getGenConstrIndicator = @import("model_genconstr.zig").getGenConstrIndicator;
    pub const getGenConstrPWL = @import("model_genconstr.zig").getGenConstrPWL;
    pub const getGenConstrPoly = @import("model_genconstr.zig").getGenConstrPoly;
    pub const getGenConstrExp = @import("model_genconstr.zig").getGenConstrExp;
    pub const getGenConstrExpA = @import("model_genconstr.zig").getGenConstrExpA;
    pub const getGenConstrLog = @import("model_genconstr.zig").getGenConstrLog;
    pub const getGenConstrLogA = @import("model_genconstr.zig").getGenConstrLogA;
    pub const getGenConstrPow = @import("model_genconstr.zig").getGenConstrPow;
    pub const getGenConstrSin = @import("model_genconstr.zig").getGenConstrSin;
    pub const getGenConstrCos = @import("model_genconstr.zig").getGenConstrCos;
    pub const getGenConstrTan = @import("model_genconstr.zig").getGenConstrTan;
    pub const getGenConstrLogistic = @import("model_genconstr.zig").getGenConstrLogistic;
    pub const getGenConstrNorm = @import("model_genconstr.zig").getGenConstrNorm;
    pub const getGenConstrNL = @import("model_genconstr.zig").getGenConstrNL;
    pub const delGenConstrs = @import("model_genconstr.zig").delGenConstrs;
    pub const getQConstrByName = @import("model_constraints.zig").getQConstrByName;
    pub const setPWLObj = @import("model_objective.zig").setPWLObj;
    pub const getPWLObj = @import("model_objective.zig").getPWLObj;
    pub const getJSonSolution = @import("model_advanced.zig").getJSonSolution;
    pub const convertToFixed = @import("model_advanced.zig").convertToFixed;
    pub const chgCoeff = @import("model_linear.zig").chgCoeff;
    pub const chgCoeffs = @import("model_linear.zig").chgCoeffs;
    pub const chgBounds = @import("model_linear.zig").chgBounds;
    pub const chgObj = @import("model_linear.zig").chgObj;
    pub const chgRHS = @import("model_linear.zig").chgRHS;
    pub const chgSense = @import("model_linear.zig").chgSense;
    pub const chgVarType = @import("model_linear.zig").chgVarType;
    pub const delVars = @import("model_linear.zig").delVars;
    pub const delConstrs = @import("model_linear.zig").delConstrs;
    pub const addRangeConstr = @import("model_linear.zig").addRangeConstr;
    pub const addRangeConstrs = @import("model_linear.zig").addRangeConstrs;
    pub const addQPterms = @import("model_constraints.zig").addQPterms;
    pub const delQ = @import("model_constraints.zig").delQ;
    pub const getQ = @import("model_constraints.zig").getQ;
    pub const delQConstrs = @import("model_constraints.zig").delQConstrs;
    pub const getQConstr = @import("model_constraints.zig").getQConstr;
    pub const delSOS = @import("model_constraints.zig").delSOS;
    pub const getSOS = @import("model_constraints.zig").getSOS;
    pub const getVarByName = @import("model_linear.zig").getVarByName;
    pub const getConstrByName = @import("model_linear.zig").getConstrByName;
    pub const getCoeff = @import("model_linear.zig").getCoeff;
    pub const getVars = @import("model_linear.zig").getVars;
    pub const getConstrs = @import("model_linear.zig").getConstrs;
    pub const setObjectiveN = @import("model_advanced.zig").setObjectiveN;
    pub const computeIIS = @import("model_advanced.zig").computeIIS;
    pub const feasRelax = @import("model_advanced.zig").feasRelax;
    pub const getCallbackFunc = @import("model_callback.zig").getCallbackFunc;
    pub const setCallbackFunc = @import("model_callback.zig").setCallbackFunc;
    pub const terminate = @import("model_callback.zig").terminate;
    pub const presolveModel = @import("model_advanced.zig").presolveModel;
    pub const fixModel = @import("model_advanced.zig").fixModel;
    pub const getBasisHead = @import("model_advanced.zig").getBasisHead;
    pub const getEnv = @import("model_params.zig").getEnv;
    pub const getAttrInfo = @import("model_params.zig").getAttrInfo;
    pub const msg = @import("model_params.zig").msg;
    pub const setLogCallbackFunc = @import("model_params.zig").setLogCallbackFunc;
    pub const getLogCallbackFunc = @import("model_params.zig").getLogCallbackFunc;
    pub const cbGet = @import("model_callback.zig").cbGet;
    pub const cbCut = @import("model_callback.zig").cbCut;
    pub const cbLazy = @import("model_callback.zig").cbLazy;
    pub const cbSolution = @import("model_callback.zig").cbSolution;
    pub const cbSetDblParam = @import("model_callback.zig").cbSetDblParam;
    pub const cbSetIntParam = @import("model_callback.zig").cbSetIntParam;
    pub const cbSetStrParam = @import("model_callback.zig").cbSetStrParam;
    pub const cbSetParam = @import("model_callback.zig").cbSetParam;
    pub const cbProceed = @import("model_callback.zig").cbProceed;
    pub const cbStopOneMultiObj = @import("model_callback.zig").cbStopOneMultiObj;
    pub const setCallbackFuncAdv = @import("model_callback.zig").setCallbackFuncAdv;
    pub const tuneModel = @import("model_params.zig").tuneModel;
    pub const getTuneResult = @import("model_params.zig").getTuneResult;
    pub const getDblParamInfo = @import("model_params.zig").getDblParamInfo;
    pub const getIntParamInfo = @import("model_params.zig").getIntParamInfo;
    pub const getStrParamInfo = @import("model_params.zig").getStrParamInfo;
    pub const copyModelToEnv = @import("model_advanced.zig").copyModelToEnv;
    pub const version = @import("model_params.zig").version;
    pub const singleScenarioModel = @import("model_advanced.zig").singleScenarioModel;
    pub const getErrormsg = @import("model_params.zig").getErrormsg;
    pub const setIntParam = @import("model_params.zig").setIntParam;
    pub const getIntParam = @import("model_params.zig").getIntParam;
    pub const setDblParam = @import("model_params.zig").setDblParam;
    pub const getDblParam = @import("model_params.zig").getDblParam;
    pub const setStrParam = @import("model_params.zig").setStrParam;
    pub const getStrParam = @import("model_params.zig").getStrParam;
    pub const setParam = @import("model_params.zig").setParam;
    pub const writeParams = @import("model_params.zig").writeParams;
    pub const readParams = @import("model_params.zig").readParams;
    pub const resetParams = @import("model_params.zig").resetParams;

    pub const writeModel = @import("model_io.zig").writeModel;
    pub const readModel = @import("model_io.zig").readModel;
    pub const read = @import("model_io.zig").read;
};

// ── Tests ───────────────────────────────────────────────────────────────

test "Model.init creates an empty model" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();

    var model = try Model.init(std.testing.allocator, &env, "empty");
    defer model.deinit();

    try std.testing.expectEqual(@as(usize, 0), model.numVars());
    try std.testing.expectEqual(@as(usize, 0), model.numConstrs());
    try std.testing.expectEqualStrings("empty", model.name);
}

test "Model.addVar queues a variable and applyPending commits it" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();

    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    // Add a constraint first so the variable's column entry is valid.
    try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, "c0");

    const vind = [_]usize{0};
    const vval = [_]f64{2.0};
    try model.addVar(1, &vind, &vval, 1.0, 0.0, 1.0, .continuous, "x1");
    try std.testing.expect(model.has_pending);

    try model.updateModel();
    try std.testing.expect(!model.has_pending);
    try std.testing.expectEqual(@as(usize, 1), model.num_vars);
    try std.testing.expectEqual(@as(usize, 1), model.num_constrs);
    try std.testing.expectEqual(@as(f64, 0.0), model.var_lb[0]);
    try std.testing.expectEqual(@as(f64, 1.0), model.var_ub[0]);
    try std.testing.expectEqual(@as(f64, 1.0), model.var_obj[0]);
}

test "Model.addVar with too-large lb returns error" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();

    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try std.testing.expectError(error.InvalidArgument, model.addVar(0, &.{}, &.{}, 0.0, 5.0, 1.0, .continuous, null));
}

test "Model attribute get/set round-trips" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();

    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.setIntAttr(.model_sense, -1);
    try std.testing.expectEqual(@as(i64, -1), try model.getIntAttr(.model_sense));
    try std.testing.expectEqual(ObjectiveSense.maximize, model.sense);

    try model.setStrAttr(.model_name, "renamed");
    try std.testing.expectEqualStrings("renamed", try model.getStrAttr(.model_name));

    try std.testing.expectEqual(@as(i64, 0), try model.getIntAttr(.status));
}

test "Model attribute get on unknown returns error" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();

    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try std.testing.expectError(error.InvalidAttribute, model.getIntAttr(.model_name));
    try std.testing.expectError(error.InvalidAttribute, model.getDblAttr(.model_name));
}

test "Model.optimize on empty model returns error" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();

    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try std.testing.expectError(error.EmptyModel, model.optimize());
}

test "Model.reset clears solution state" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();

    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    model.obj_val = 42.0;
    model.iter_count = 100;
    model.status = .optimal;
    model.reset(false);

    try std.testing.expectEqual(@as(f64, 0.0), model.obj_val);
    try std.testing.expectEqual(@as(i64, 0), model.iter_count);
    try std.testing.expectEqual(Status.loaded, model.status);
}

// ═══════════════════════════════════════════════════════════════════════════
//  New API method tests
// ═══════════════════════════════════════════════════════════════════════════

test "Model.getVarByName finds a named variable" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, "c0");
    try model.addVar(0, &.{}, &.{}, 1.0, 0.0, 1.0, .continuous, "x1");
    try model.updateModel();

    const idx = try model.getVarByName("x1");
    try std.testing.expectEqual(@as(usize, 0), idx);
    try std.testing.expectError(error.NotInModel, model.getVarByName("nonexistent"));
}

test "Model.getConstrByName finds a named constraint" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addConstr(1, &[_]usize{0}, &[_]f64{1.0}, .less_equal, 10.0, "c0");
    try model.addVar(1, &[_]usize{0}, &[_]f64{1.0}, 1.0, 0.0, 1.0, .continuous, null);
    try model.updateModel();

    const idx = try model.getConstrByName("c0");
    try std.testing.expectEqual(@as(usize, 0), idx);
}

test "Model.getCoeff returns stored coefficient" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, "c0");
    try model.addVar(1, &[_]usize{0}, &[_]f64{2.5}, 1.0, 0.0, 1.0, .continuous, null);
    try model.updateModel();

    const val = try model.getCoeff(0, 0);
    try std.testing.expectApproxEqAbs(@as(f64, 2.5), val, 1e-12);
}

test "Model.addRangeConstr creates two constraints for finite lower and upper" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addVar(0, &.{}, &.{}, 1.0, 0.0, 1.0, .continuous, null);
    try model.addRangeConstr(1, &[_]usize{0}, &[_]f64{1.0}, 5.0, 10.0, "range");
    try model.updateModel();

    try std.testing.expectEqual(@as(usize, 2), model.num_constrs);
}

test "Model.addQPterms appends quadratic objective terms" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addQPterms(&[_]i32{ 0, 0 }, &[_]i32{ 0, 1 }, &[_]f64{ 2.0, 1.0 });
    try std.testing.expectEqual(@as(usize, 2), model.q_nz);
    try std.testing.expectEqual(@as(f64, 2.0), model.q_val[0]);
}

test "Model.addQPterms delQ removes all terms" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addQPterms(&[_]i32{0}, &[_]i32{0}, &[_]f64{1.0});
    try std.testing.expectEqual(@as(usize, 1), model.q_nz);
    model.delQ();
    try std.testing.expectEqual(@as(usize, 0), model.q_nz);
}

test "Model.addQConstr respects explicit nonzero counts" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addQConstr(
        1,
        &[_]i32{ 0, 99 },
        &[_]i32{ 0, 99 },
        &[_]f64{ 2.0, 99.0 },
        1,
        &[_]usize{ 0, 99 },
        &[_]f64{ 1.0, 99.0 },
        .less_equal,
        3.0,
        "qconstr",
    );

    try std.testing.expectEqual(@as(usize, 1), model.qconstr_qrow.len);
    try std.testing.expectEqual(@as(usize, 1), model.qconstr_lind.len);
    try std.testing.expectEqual(@as(f64, 2.0), model.qconstr_qval[0]);
    try std.testing.expectEqual(@as(f64, 1.0), model.qconstr_lval[0]);
}

test "Model.addSOS stores SOS constraint" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addVar(0, &.{}, &.{}, 1.0, 0.0, 1.0, .continuous, null);
    try model.addVar(0, &.{}, &.{}, 1.0, 0.0, 1.0, .continuous, null);
    try model.addSOS(.sos1, 2, &[_]usize{ 0, 1 }, null, "sos1");
    try std.testing.expectEqual(@as(usize, 1), model.sos_count);

    var st: SosType = .sos1;
    var nm: usize = 0;
    var indices: [2]usize = undefined;
    var weights: [2]f64 = undefined;
    try model.getSOS(0, &st, &nm, &indices, &weights);
    try std.testing.expectEqual(SosType.sos1, st);
    try std.testing.expectEqual(@as(usize, 2), nm);
    try std.testing.expectEqual(@as(f64, 1.0), weights[0]);
    try std.testing.expectEqual(@as(f64, 2.0), weights[1]);
}

test "Model.addGenConstrMax stores general constraint" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addVar(0, &.{}, &.{}, 1.0, 0.0, 1.0, .continuous, null);
    try model.addVar(0, &.{}, &.{}, 1.0, 0.0, 1.0, .continuous, null);
    try model.addGenConstrMax(0, 2, &[_]usize{ 0, 1 }, 0.0, "max1");
    try std.testing.expectEqual(@as(usize, 1), model.genconstr_count);
    try std.testing.expectEqual(GenConstrType.max, model.genconstr_types[0]);
}

test "packed general constraints preserve following constraint offsets" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addGenConstrPWL(0, 1, 2, &[_]f64{ 0.0, 1.0 }, &[_]f64{ 0.0, 2.0 }, "pwl");
    try model.addGenConstrMax(2, 2, &[_]usize{ 3, 4 }, 0.0, "max");

    var resvar: usize = undefined;
    var num_vars: usize = undefined;
    var vars: [2]usize = undefined;
    try model.getGenConstrMax(1, &resvar, &num_vars, &vars);
    try std.testing.expectEqual(@as(usize, 2), resvar);
    try std.testing.expectEqual(@as(usize, 2), num_vars);
    try std.testing.expectEqualSlices(usize, &[_]usize{ 3, 4 }, &vars);
}

test "Model.terminate sets interrupt flag" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try std.testing.expect(!model.interrupted);
    model.terminate();
    try std.testing.expect(model.interrupted);
}

test "Model.chgCoeffs queues multiple coefficient changes" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, "c0");
    try model.addVar(1, &[_]usize{0}, &[_]f64{1.0}, 1.0, 0.0, 1.0, .continuous, null);
    try model.chgCoeffs(1, &[_]usize{0}, &[_]usize{0}, &[_]f64{5.0});
    try model.updateModel();

    try std.testing.expectApproxEqAbs(@as(f64, 5.0), try model.getCoeff(0, 0), 1e-12);
}

test "Model.setDblAttrArray sets multiple LB values" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, "c0");
    try model.addVar(0, &.{}, &.{}, 1.0, 0.0, 1.0, .continuous, null);
    try model.addVar(0, &.{}, &.{}, 1.0, 0.0, 1.0, .continuous, null);
    try model.updateModel();

    try model.setDblAttrArray(.lb, 0, &[_]f64{ 1.0, 2.0 });
    try std.testing.expectEqual(@as(f64, 1.0), try model.getDblAttrElement(.lb, 0));
    try std.testing.expectEqual(@as(f64, 2.0), try model.getDblAttrElement(.lb, 1));
}

test "Model.setDblAttrList sets selected LB values" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, "c0");
    try model.addVar(0, &.{}, &.{}, 1.0, 0.0, 1.0, .continuous, null);
    try model.addVar(0, &.{}, &.{}, 1.0, 0.0, 1.0, .continuous, null);
    try model.addVar(0, &.{}, &.{}, 1.0, 0.0, 1.0, .continuous, null);
    try model.updateModel();

    try model.setDblAttrList(.lb, &[_]usize{ 0, 2 }, &[_]f64{ 5.0, 10.0 });
    try std.testing.expectEqual(@as(f64, 5.0), try model.getDblAttrElement(.lb, 0));
    try std.testing.expectEqual(@as(f64, 0.0), try model.getDblAttrElement(.lb, 1)); // unchanged
    try std.testing.expectEqual(@as(f64, 10.0), try model.getDblAttrElement(.lb, 2));
}

test "Model.getDblAttrList retrieves selected attribute values" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, "c0");
    try model.addVar(0, &.{}, &.{}, 5.0, 1.0, 10.0, .continuous, null);
    try model.addVar(0, &.{}, &.{}, 8.0, 2.0, 20.0, .continuous, null);
    try model.updateModel();

    var obj_vals: [2]f64 = undefined;
    try model.getDblAttrList(.obj, &[_]usize{ 0, 1 }, &obj_vals);
    try std.testing.expectEqual(@as(f64, 5.0), obj_vals[0]);
    try std.testing.expectEqual(@as(f64, 8.0), obj_vals[1]);
}

test "Model.copy creates independent copy" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "original");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, "c0");
    try model.addVar(1, &[_]usize{0}, &[_]f64{1.0}, 1.0, 0.0, 1.0, .continuous, "x1");
    try model.addVar(1, &[_]usize{0}, &[_]f64{2.0}, 2.0, 0.0, 2.0, .continuous, "x2");
    try model.updateModel();

    var copy = try model.copy("copy");
    defer copy.deinit();

    try std.testing.expectEqual(@as(usize, 2), copy.numVars());
    try std.testing.expectEqual(@as(usize, 1), copy.numConstrs());
    try std.testing.expectEqual(@as(f64, 0.0), try copy.getDblAttrElement(.lb, 0));
    try std.testing.expectEqual(@as(f64, 2.0), try copy.getDblAttrElement(.ub, 1));
    try std.testing.expectEqualStrings("copy", copy.name);
}
