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
//! 3. **Stable handles over dense storage** – public variable/constraint
//!    handles carry a generation-checked identity, while the numerical core
//!    remains packed in contiguous SoA arrays for cache-friendly traversal.
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
const compiled_model_view = @import("compiled_model_view.zig");
const revision_module = @import("revisions.zig");
const solver = @import("solver");

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
const VarId = @import("entity_handle.zig").VarId;
const VarHandleTable = @import("entity_handle.zig").HandleTableFor(VarId);
const ConstrId = @import("entity_handle.zig").ConstrId;
const ConstrHandleTable = @import("entity_handle.zig").HandleTableFor(ConstrId);
const QConstrId = @import("entity_handle.zig").QConstrId;
const QConstrHandleTable = @import("entity_handle.zig").HandleTableFor(QConstrId);
const SosId = @import("entity_handle.zig").SosId;
const SosHandleTable = @import("entity_handle.zig").HandleTableFor(SosId);
const GenConstrId = @import("entity_handle.zig").GenConstrId;
const GenConstrHandleTable = @import("entity_handle.zig").HandleTableFor(GenConstrId);

// ── Model ───────────────────────────────────────────────────────────────

pub const Model = struct {
    allocator: std.mem.Allocator,
    env: *Env,

    // ── Identity ───────────────────────────────────────────────────────
    name: []const u8,

    // ── Variable data (committed, dense SoA) ──────────────────────────
    // These arrays are deliberately mutable and compacted after deletion;
    // they are not append-only logs.  Stable VarId generations decouple API
    // references from the moving dense slots used by numerical kernels.
    num_vars: usize = 0,
    var_lb: []f64 = &.{},
    var_ub: []f64 = &.{},
    var_obj: []f64 = &.{},
    var_type: []VarType = &.{},
    var_names: [](?[]const u8) = &.{},
    // Stable API identities mapped to the dense variable arrays above.
    var_handles: VarHandleTable = .{},

    // ── Constraint data (committed) ────────────────────────────────────
    num_constrs: usize = 0,
    constr_sense: []Sense = &.{},
    constr_rhs: []f64 = &.{},
    constr_names: [](?[]const u8) = &.{},
    // Stable API identities mapped to the dense constraint arrays above.
    constr_handles: ConstrHandleTable = .{},

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
    revisions: revision_module.RevisionSet = .{},
    // Derived canonical row bounds reused by zero-copy LP solve views.
    compiled_view_cache: compiled_model_view.CompiledModelViewCache = .{},
    lp_session: solver.LpSolveSession,

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
    /// Global structural/logical column occupying every basis row.
    basis_head: []u32 = &.{},
    basis_available: bool = false,

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
    sos_handles: SosHandleTable = .{},

    // ── Quadratic constraints ──────────────────────────────────────────
    qconstr_count: usize = 0,
    qconstr_qrow: []i32 = &.{},
    qconstr_qcol: []i32 = &.{},
    qconstr_qval: []f64 = &.{},
    qconstr_qbegin: []usize = &.{},
    qconstr_lind: []usize = &.{},
    qconstr_lval: []f64 = &.{},
    qconstr_lbegin: []usize = &.{},
    qconstr_sense: []Sense = &.{},
    qconstr_rhs: []f64 = &.{},
    qconstr_names: [](?[]const u8) = &.{},
    qconstr_handles: QConstrHandleTable = .{},

    // ── General constraints ────────────────────────────────────────────
    genconstr_count: usize = 0,
    genconstr_types: []GenConstrType = &.{},
    genconstr_resvar: []usize = &.{},
    genconstr_nvars: []usize = &.{}, // number of vars per constraint
    genconstr_indices: []usize = &.{}, // packed variable indices
    genconstr_begin: []usize = &.{}, // packed-data offsets, length count + 1
    genconstr_names: [](?[]const u8) = &.{},
    genconstr_handles: GenConstrHandleTable = .{},

    // ── Piecewise-linear objective data ──────────────────────────────────
    pwlobj_count: usize = 0,
    pwlobj_var: []usize = &.{}, // variable index for each PWL objective
    pwlobj_npts: []usize = &.{}, // number of points per entry
    pwlobj_xdata: []f64 = &.{}, // packed x values
    pwlobj_ydata: []f64 = &.{}, // packed y values

    // ── Callback / interrupt state ─────────────────────────────────────
    interrupted: std.atomic.Value(bool) = std.atomic.Value(bool).init(false),

    const Self = @This();

    pub fn markRevision(self: *Self, kind: revision_module.RevisionKind) ModelError!void {
        self.revisions.bump(kind) catch return error.RevisionOverflow;
    }

    // ── Construction / destruction ─────────────────────────────────────

    /// Create a new empty model within the given environment.
    pub fn init(allocator: std.mem.Allocator, env: *Env, name: []const u8) (ModelError)!Self {
        var self = Self{
            .allocator = allocator,
            .env = env,
            .name = "",
            .matrix = undefined,
            .lp_session = solver.LpSolveSession.init(allocator),
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
        self.var_handles.deinit(alloc);
        alloc.free(self.constr_sense);
        alloc.free(self.constr_rhs);
        for (self.constr_names) |n| if (n) |s| alloc.free(s);
        alloc.free(self.constr_names);
        self.constr_handles.deinit(alloc);
        alloc.free(self.q_row);
        alloc.free(self.q_col);
        alloc.free(self.q_val);
        self.compiled_view_cache.deinit(alloc);
        self.lp_session.deinit();
        alloc.free(self.solution);
        alloc.free(self.reduced_cost);
        alloc.free(self.slack);
        alloc.free(self.pi);
        alloc.free(self.vbasis);
        alloc.free(self.cbasis);
        alloc.free(self.basis_head);
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
        self.sos_handles.deinit(alloc);
        // Quadratic constraints
        alloc.free(self.qconstr_qrow);
        alloc.free(self.qconstr_qcol);
        alloc.free(self.qconstr_qval);
        alloc.free(self.qconstr_qbegin);
        alloc.free(self.qconstr_lind);
        alloc.free(self.qconstr_lval);
        alloc.free(self.qconstr_lbegin);
        alloc.free(self.qconstr_sense);
        alloc.free(self.qconstr_rhs);
        for (self.qconstr_names) |n| if (n) |s| alloc.free(s);
        alloc.free(self.qconstr_names);
        self.qconstr_handles.deinit(alloc);
        // General constraints
        alloc.free(self.genconstr_types);
        alloc.free(self.genconstr_resvar);
        alloc.free(self.genconstr_nvars);
        alloc.free(self.genconstr_indices);
        alloc.free(self.genconstr_begin);
        for (self.genconstr_names) |n| if (n) |s| alloc.free(s);
        alloc.free(self.genconstr_names);
        self.genconstr_handles.deinit(alloc);
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

    /// Resolve a stable variable identity to its current dense position.
    pub fn resolveVarId(self: Self, id: VarId) ModelError!usize {
        const dense = self.var_handles.resolve(id) catch |err| return switch (err) {
            error.InvalidHandle => error.NotInModel,
            error.HandleExhausted => error.OutOfMemory,
            error.DenseIndexOutOfRange => error.InvalidArgument,
        };
        return @intCast(dense);
    }

    /// Return the stable identity at a committed dense variable position.
    pub fn varIdAt(self: *Self, index: usize) ModelError!VarId {
        if (index >= self.num_vars) return error.IndexOutOfRange;
        while (self.var_handles.liveLen() < self.num_vars) {
            const id = self.var_handles.allocate(self.allocator) catch return error.OutOfMemory;
            self.var_handles.bindDenseWithAllocator(self.allocator, id, @intCast(self.var_handles.liveLen())) catch return error.OutOfMemory;
        }
        return self.var_handles.idAtDense(index) catch return error.NotInModel;
    }

    pub fn resolveConstrId(self: Self, id: ConstrId) ModelError!usize {
        const dense = self.constr_handles.resolve(id) catch |err| return switch (err) {
            error.InvalidHandle => error.NotInModel,
            error.HandleExhausted => error.OutOfMemory,
            error.DenseIndexOutOfRange => error.InvalidArgument,
        };
        return @intCast(dense);
    }

    pub fn constrIdAt(self: *Self, index: usize) ModelError!ConstrId {
        if (index >= self.num_constrs) return error.IndexOutOfRange;
        while (self.constr_handles.liveLen() < self.num_constrs) {
            const id = self.constr_handles.allocate(self.allocator) catch return error.OutOfMemory;
            self.constr_handles.bindDenseWithAllocator(self.allocator, id, @intCast(self.constr_handles.liveLen())) catch return error.OutOfMemory;
        }
        return self.constr_handles.idAtDense(index) catch return error.NotInModel;
    }

    pub fn resolveQConstrId(self: Self, id: QConstrId) ModelError!usize {
        return @intCast(self.qconstr_handles.resolve(id) catch return error.NotInModel);
    }

    pub fn qconstrIdAt(self: *Self, index: usize) ModelError!QConstrId {
        if (index >= self.qconstr_count) return error.IndexOutOfRange;
        while (self.qconstr_handles.liveLen() < self.qconstr_count) {
            const id = self.qconstr_handles.allocate(self.allocator) catch return error.OutOfMemory;
            self.qconstr_handles.bindDenseWithAllocator(self.allocator, id, @intCast(self.qconstr_handles.liveLen())) catch return error.OutOfMemory;
        }
        return self.qconstr_handles.idAtDense(index) catch return error.NotInModel;
    }

    pub fn resolveSosId(self: Self, id: SosId) ModelError!usize {
        return @intCast(self.sos_handles.resolve(id) catch return error.NotInModel);
    }

    pub fn sosIdAt(self: *Self, index: usize) ModelError!SosId {
        if (index >= self.sos_count) return error.IndexOutOfRange;
        while (self.sos_handles.liveLen() < self.sos_count) {
            const id = self.sos_handles.allocate(self.allocator) catch return error.OutOfMemory;
            self.sos_handles.bindDenseWithAllocator(self.allocator, id, @intCast(self.sos_handles.liveLen())) catch return error.OutOfMemory;
        }
        return self.sos_handles.idAtDense(index) catch return error.NotInModel;
    }

    pub fn resolveGenConstrId(self: Self, id: GenConstrId) ModelError!usize {
        return @intCast(self.genconstr_handles.resolve(id) catch return error.NotInModel);
    }

    pub fn genconstrIdAt(self: *Self, index: usize) ModelError!GenConstrId {
        if (index >= self.genconstr_count) return error.IndexOutOfRange;
        while (self.genconstr_handles.liveLen() < self.genconstr_count) {
            const id = self.genconstr_handles.allocate(self.allocator) catch return error.OutOfMemory;
            self.genconstr_handles.bindDenseWithAllocator(self.allocator, id, @intCast(self.genconstr_handles.liveLen())) catch return error.OutOfMemory;
        }
        return self.genconstr_handles.idAtDense(index) catch return error.NotInModel;
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
    pub const addVarColumn = @import("model_linear.zig").addVarColumn;
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
    pub const addQPtermsExpr = @import("model_constraints.zig").addQPtermsExpr;
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
    pub const exportBasisSnapshot = @import("model_advanced.zig").exportBasisSnapshot;
    pub const importBasisView = @import("model_advanced.zig").importBasisView;
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
    pub const write = @import("model_io.zig").write;
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

test "Model.optimize publishes an LP solution and basis attributes" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "lp");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 5.0, "capacity");
    try model.addVar(1, &[_]usize{0}, &[_]f64{1.0}, 1.0, 0.0, std.math.inf(f64), .continuous, "x");
    try model.setIntAttr(.model_sense, -1);
    try model.optimize();

    try std.testing.expectEqual(Status.optimal, model.status);
    try std.testing.expectApproxEqAbs(@as(f64, 5.0), model.obj_val, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 5.0), model.solution[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), model.slack[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), model.pi[0], 1e-12);
    try std.testing.expectEqual(BasisStatus.basic, model.vbasis[0]);
    try std.testing.expectEqual(@as(usize, 1), model.lp_session.cold_solves);
    const factor_buffer = model.lp_session.engine.factorization.dense_lu.lu.ptr;

    const cached_revision = model.compiled_view_cache.source_revision.?;
    const cached_row_upper_ptr = model.compiled_view_cache.row_upper.ptr;
    try model.optimize();
    try std.testing.expectEqual(cached_revision, model.compiled_view_cache.source_revision.?);
    try std.testing.expectEqual(cached_row_upper_ptr, model.compiled_view_cache.row_upper.ptr);
    try std.testing.expectEqual(@as(usize, 1), model.lp_session.reoptimizations);
    try std.testing.expectEqual(factor_buffer, model.lp_session.engine.factorization.dense_lu.lu.ptr);

    try model.setDblAttrElement(.rhs, 0, 4.0);
    try model.optimize();
    try std.testing.expect(model.compiled_view_cache.source_revision.? > cached_revision);
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), model.compiled_view_cache.row_upper[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), model.solution[0], 1e-12);
    try std.testing.expectEqual(@as(usize, 2), model.lp_session.reoptimizations);
}

test "matrix-value revisions invalidate a persistent LP session" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "matrix_revision");
    defer model.deinit();
    try model.addConstr(0, &.{}, &.{}, .less_equal, 4.0, null);
    try model.addVar(1, &[_]usize{0}, &[_]f64{1.0}, -1.0, 0.0, std.math.inf(f64), .continuous, null);
    try model.optimize();
    try std.testing.expectEqual(@as(usize, 1), model.lp_session.cold_solves);

    try model.chgCoeff(0, 0, 2.0);
    try model.optimize();
    try std.testing.expectEqual(@as(usize, 2), model.lp_session.cold_solves);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), model.solution[0], 1e-12);
}

test "Model.optimize publishes infeasible status without stale solution arrays" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "infeasible_lp");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 1.0, null);
    try model.addVar(1, &[_]usize{0}, &[_]f64{1.0}, 0.0, 2.0, std.math.inf(f64), .continuous, null);
    try model.optimize();

    try std.testing.expectEqual(Status.infeasible, model.status);
    try std.testing.expectEqual(@as(usize, 0), model.solution.len);
    try std.testing.expectEqual(@as(usize, 0), model.pi.len);
}

test "Model.optimize maps an immediate TimeLimit" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    try env.setDblParam("TimeLimit", 0.0);
    var model = try Model.init(std.testing.allocator, &env, "limited_lp");
    defer model.deinit();

    try model.addVar(0, &.{}, &.{}, 0.0, 0.0, 1.0, .continuous, null);
    try model.optimize();

    try std.testing.expectEqual(Status.time_limit, model.status);
    try std.testing.expectEqual(@as(usize, 0), model.solution.len);
}

test "Model.optimize maps an immediate deterministic WorkLimit" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    try env.setDblParam("WorkLimit", 0.0);
    var model = try Model.init(std.testing.allocator, &env, "work_limited_lp");
    defer model.deinit();

    try model.addVar(0, &.{}, &.{}, -1.0, 0.0, 1.0, .continuous, null);
    try model.optimize();

    try std.testing.expectEqual(Status.work_limit, model.status);
    try std.testing.expectEqual(@as(u64, 0), model.lp_session.engine.work_used);
    try std.testing.expectEqual(@as(usize, 0), model.solution.len);
}

test "Model callback bridge reports simplex iterations without allocation" {
    const CallbackState = struct {
        calls: usize = 0,
        last_where: CallbackWhere = .polling,

        fn callback(where: CallbackWhere, user_data: ?*anyopaque) void {
            const self: *@This() = @ptrCast(@alignCast(user_data.?));
            self.calls += 1;
            self.last_where = where;
        }
    };
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var callback_state = CallbackState{};
    env.setCallbackFunc(CallbackState.callback, &callback_state);
    var model = try Model.init(std.testing.allocator, &env, "callback_lp");
    defer model.deinit();

    try model.addVar(0, &.{}, &.{}, -1.0, 0.0, 1.0, .continuous, null);
    try model.optimize();

    try std.testing.expectEqual(Status.optimal, model.status);
    try std.testing.expect(callback_state.calls >= 1);
    try std.testing.expectEqual(CallbackWhere.simplex, callback_state.last_where);
}

test "Model.optimize warm-starts RHS changes through dual simplex" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "dual_reopt");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .greater_equal, 2.0, null);
    try model.addVar(1, &[_]usize{0}, &[_]f64{1.0}, 1.0, 0.0, std.math.inf(f64), .continuous, null);
    try model.optimize();
    try std.testing.expectApproxEqAbs(@as(f64, 2), model.solution[0], 1e-12);
    try std.testing.expectEqual(@as(usize, 1), model.basis_head.len);

    try model.setDblAttrElement(.rhs, 0, -1.0);
    try model.optimize();
    try std.testing.expectEqual(Status.optimal, model.status);
    try std.testing.expectApproxEqAbs(@as(f64, 0), model.solution[0], 1e-12);
    try std.testing.expect(model.iter_count <= 1);
}

test "Model.optimize warm-starts bound changes through dual simplex" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "bound_reopt");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 5.0, null);
    try model.addVar(1, &[_]usize{0}, &[_]f64{1.0}, 1.0, 0.0, std.math.inf(f64), .continuous, null);
    try model.setIntAttr(.model_sense, -1);
    try model.optimize();
    try std.testing.expectApproxEqAbs(@as(f64, 5), model.solution[0], 1e-12);

    try model.setDblAttrElement(.ub, 0, 3.0);
    try model.optimize();
    try std.testing.expectEqual(Status.optimal, model.status);
    try std.testing.expectApproxEqAbs(@as(f64, 3), model.solution[0], 1e-12);
    try std.testing.expect(model.iter_count <= 1);
}

test "Model.optimize warm-starts objective changes through primal simplex" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "objective_reopt");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 5.0, null);
    try model.addVar(1, &[_]usize{0}, &[_]f64{1.0}, -1.0, 0.0, std.math.inf(f64), .continuous, null);
    try model.optimize();
    try std.testing.expectApproxEqAbs(@as(f64, 5), model.solution[0], 1e-12);

    try model.setDblAttrElement(.obj, 0, 1.0);
    try model.optimize();
    try std.testing.expectEqual(Status.optimal, model.status);
    try std.testing.expectApproxEqAbs(@as(f64, 0), model.solution[0], 1e-12);
    try std.testing.expect(model.iter_count <= 1);
}

test "dual BFRT flips a boxed variable before the entering pivot" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "bfrt_reopt");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .greater_equal, 1.0, null);
    try model.addVar(1, &[_]usize{0}, &[_]f64{1.0}, 2.0, 0.0, std.math.inf(f64), .continuous, "basic");
    try model.addVar(1, &[_]usize{0}, &[_]f64{1.0}, 1.0, 0.0, 0.5, .continuous, "boxed");
    try model.optimize();
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), model.solution[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), model.solution[1], 1e-12);

    try model.setDblAttrElement(.rhs, 0, -1.0);
    try model.optimize();
    try std.testing.expectEqual(Status.optimal, model.status);
    try std.testing.expectApproxEqAbs(@as(f64, 0), model.solution[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0), model.solution[1], 1e-12);
    try std.testing.expect(model.iter_count <= 1);
}

test "dual warm reoptimization proves infeasibility after a bound change" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "dual_infeasible");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .greater_equal, 2.0, null);
    try model.addVar(1, &[_]usize{0}, &[_]f64{1.0}, 1.0, 0.0, std.math.inf(f64), .continuous, null);
    try model.optimize();
    try std.testing.expectEqual(Status.optimal, model.status);

    try model.setDblAttrElement(.ub, 0, 1.0);
    try model.optimize();
    try std.testing.expectEqual(Status.infeasible, model.status);
    try std.testing.expectEqual(@as(usize, 0), model.solution.len);
}

test "Model exports and imports an owning basis snapshot" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var source = try Model.init(std.testing.allocator, &env, "basis_source");
    defer source.deinit();
    try source.addConstr(0, &.{}, &.{}, .less_equal, 5.0, null);
    try source.addVar(1, &[_]usize{0}, &[_]f64{1.0}, -1.0, 0.0, std.math.inf(f64), .continuous, null);
    try source.optimize();
    var snapshot = try source.exportBasisSnapshot(std.testing.allocator);
    defer snapshot.deinit();

    var target = try Model.init(std.testing.allocator, &env, "basis_target");
    defer target.deinit();
    try target.addConstr(0, &.{}, &.{}, .less_equal, 5.0, null);
    try target.addVar(1, &[_]usize{0}, &[_]f64{1.0}, -1.0, 0.0, std.math.inf(f64), .continuous, null);
    try target.updateModel();
    try target.importBasisView(snapshot.view());
    try target.optimize();
    try std.testing.expectEqual(Status.optimal, target.status);
    try std.testing.expectApproxEqAbs(@as(f64, 5), target.solution[0], 1e-12);
    try std.testing.expectEqual(@as(i64, 0), target.iter_count);
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

test "Var.at resolves through a stable id" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addVar(0, &.{}, &.{}, 0.0, 2.0, 4.0, .continuous, "x");
    try model.updateModel();

    var variable = try @import("var/index.zig").Var.at(&model, 0);
    try std.testing.expectEqual(@as(f64, 2.0), try variable.getLB());
    try std.testing.expect(variable.id != null);

    const stable_id = variable.id.?;
    try model.addVar(0, &.{}, &.{}, 0.0, 3.0, 6.0, .continuous, "y");
    try model.updateModel();
    variable.id = stable_id;
    try std.testing.expectEqual(@as(f64, 2.0), try variable.getLB());
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

test "Constr.at resolves through a stable id" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 5.0, "c0");
    try model.updateModel();

    var constr = try @import("constraint/index.zig").Constr.at(&model, 0);
    try std.testing.expectEqual(@as(f64, 5.0), try constr.getRHS());
    try std.testing.expect(constr.id != null);
}

test "linear deletion compacts dense data and preserves surviving Var handles" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addVar(0, &.{}, &.{}, 0.0, 1.0, 2.0, .continuous, "x");
    try model.addVar(0, &.{}, &.{}, 0.0, 3.0, 4.0, .continuous, "y");
    try model.updateModel();
    var survivor = try @import("var/index.zig").Var.at(&model, 1);

    try model.delVars(&[_]usize{0});
    try model.updateModel();
    try std.testing.expectEqual(@as(usize, 1), model.num_vars);
    survivor.id = survivor.id.?;
    try std.testing.expectEqual(@as(f64, 3.0), try survivor.getLB());
    try std.testing.expectError(error.IndexOutOfRange, @import("var/index.zig").Var.at(&model, 1));
}

test "linear deletion rebuilds canonical matrix in packed storage" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "packed_delete");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, null);
    try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, null);
    try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, null);
    try model.addVar(2, &[_]usize{ 0, 1 }, &[_]f64{ 1.0, 2.0 }, 0.0, 0.0, 1.0, .continuous, null);
    try model.addVar(2, &[_]usize{ 0, 2 }, &[_]f64{ 3.0, 4.0 }, 0.0, 0.0, 1.0, .continuous, null);
    try model.addVar(2, &[_]usize{ 1, 2 }, &[_]f64{ 5.0, 6.0 }, 0.0, 0.0, 1.0, .continuous, null);
    try model.updateModel();

    try model.delVars(&[_]usize{1});
    try model.delConstrs(&[_]usize{1});
    try model.updateModel();

    const rebuilt = model.matrix.csc();
    try rebuilt.validate();
    try std.testing.expect(rebuilt.storage != null);
    try std.testing.expectEqual(@as(usize, 2), rebuilt.num_rows);
    try std.testing.expectEqual(@as(usize, 2), rebuilt.num_cols);
    try std.testing.expectEqualSlices(usize, &.{ 0, 1, 2 }, rebuilt.col_starts);
    try std.testing.expectEqual(@as(usize, 0), rebuilt.row_indices[0].toUsize());
    try std.testing.expectEqual(@as(usize, 1), rebuilt.row_indices[1].toUsize());
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 6.0 }, rebuilt.values);
}

test "addVarColumn resolves stable constraint ids" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, "c0");
    try model.updateModel();
    const constr = try @import("constraint/index.zig").Constr.at(&model, 0);
    var column = @import("expr/column.zig").Column.init();
    defer column.deinit(std.testing.allocator);
    try column.addTerm(std.testing.allocator, constr, 2.5);

    try model.addVarColumn(column, 0.0, 0.0, 1.0, .continuous, "x");
    try model.updateModel();
    try std.testing.expectEqual(@as(usize, 1), model.num_vars);
    try std.testing.expectEqual(@as(usize, 1), model.matrix.csc().nnz());
    try std.testing.expectEqual(@as(usize, 0), model.matrix.csc().row_indices[0].toUsize());
}

test "addQPtermsExpr resolves stable variable ids" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addVar(0, &.{}, &.{}, 0.0, 0.0, 1.0, .continuous, "x");
    try model.addVar(0, &.{}, &.{}, 0.0, 0.0, 1.0, .continuous, "y");
    try model.updateModel();
    const x = try @import("var/index.zig").Var.at(&model, 0);
    const y = try @import("var/index.zig").Var.at(&model, 1);
    var expr = @import("expr/quad_expr.zig").QuadExpr.init(std.testing.allocator);
    defer expr.deinit();
    try expr.addQTerm(x, y, 3.0);

    try model.addQPtermsExpr(expr);
    try std.testing.expectEqual(@as(usize, 1), model.q_nz);
    try std.testing.expectEqual(@as(f64, 3.0), model.q_val[0]);
}

test "linear variable deletion remaps quadratic objective terms" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addVar(0, &.{}, &.{}, 0.0, 0.0, 1.0, .continuous, "x");
    try model.addVar(0, &.{}, &.{}, 0.0, 0.0, 1.0, .continuous, "y");
    try model.updateModel();
    try model.addQPterms(&[_]i32{ 0, 1 }, &[_]i32{ 0, 1 }, &[_]f64{ 1.0, 2.0 });
    try model.delVars(&[_]usize{0});
    try model.updateModel();
    try std.testing.expectEqual(@as(usize, 1), model.q_nz);
    try std.testing.expectEqual(@as(i32, 0), model.q_row[0]);
    try std.testing.expectEqual(@as(i32, 0), model.q_col[0]);
    try std.testing.expectEqual(@as(f64, 2.0), model.q_val[0]);
}

test "linear variable deletion remaps piecewise objective storage" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addVar(0, &.{}, &.{}, 0.0, 0.0, 1.0, .continuous, "x");
    try model.addVar(0, &.{}, &.{}, 0.0, 0.0, 1.0, .continuous, "y");
    try model.updateModel();
    try model.setPWLObj(0, 2, &[_]f64{ 0.0, 1.0 }, &[_]f64{ 0.0, 2.0 });
    try model.setPWLObj(1, 2, &[_]f64{ 0.0, 1.0 }, &[_]f64{ 0.0, 3.0 });
    try model.delVars(&[_]usize{0});
    try model.updateModel();
    try std.testing.expectEqual(@as(usize, 1), model.pwlobj_count);
    try std.testing.expectEqual(@as(usize, 0), model.pwlobj_var[0]);
    try std.testing.expectEqual(@as(f64, 3.0), model.pwlobj_ydata[1]);
}

test "linear variable deletion remaps SOS members" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addVar(0, &.{}, &.{}, 0.0, 0.0, 1.0, .continuous, "x");
    try model.addVar(0, &.{}, &.{}, 0.0, 0.0, 1.0, .continuous, "y");
    try model.updateModel();
    try model.addSOS(.sos1, 2, &[_]usize{ 0, 1 }, &[_]f64{ 2.0, 4.0 }, "s");
    try model.delVars(&[_]usize{0});
    try model.updateModel();
    try std.testing.expectEqual(@as(usize, 1), model.sos_indices.len);
    try std.testing.expectEqual(@as(usize, 0), model.sos_indices[0]);
    try std.testing.expectEqual(@as(f64, 4.0), model.sos_weights[0]);
    try std.testing.expectEqual(@as(usize, 0), model.sos_begin[0]);
    try std.testing.expectEqual(@as(usize, 1), model.sos_begin[1]);
}

test "linear variable deletion remaps quadratic constraint terms" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addVar(0, &.{}, &.{}, 0.0, 0.0, 1.0, .continuous, "x");
    try model.addVar(0, &.{}, &.{}, 0.0, 0.0, 1.0, .continuous, "y");
    try model.updateModel();
    try model.addQConstr(2, &[_]i32{ 0, 1 }, &[_]i32{ 0, 1 }, &[_]f64{ 1.0, 2.0 }, 1, &[_]usize{1}, &[_]f64{3.0}, .less_equal, 4.0, "q");
    try model.delVars(&[_]usize{0});
    try model.updateModel();
    try std.testing.expectEqual(@as(usize, 1), model.qconstr_qrow.len);
    try std.testing.expectEqual(@as(i32, 0), model.qconstr_qrow[0]);
    try std.testing.expectEqual(@as(f64, 2.0), model.qconstr_qval[0]);
    try std.testing.expectEqual(@as(usize, 1), model.qconstr_lind.len);
    try std.testing.expectEqual(@as(usize, 0), model.qconstr_lind[0]);
}

test "linear variable deletion remaps and filters general constraints" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();

    try model.addVar(0, &.{}, &.{}, 0.0, 0.0, 1.0, .continuous, "x");
    try model.addVar(0, &.{}, &.{}, 0.0, 0.0, 1.0, .continuous, "y");
    try model.addVar(0, &.{}, &.{}, 0.0, 0.0, 1.0, .continuous, "z");
    try model.updateModel();
    try model.addGenConstrMax(2, 2, &[_]usize{ 0, 1 }, 0.0, "drop");
    try model.addGenConstrExp(2, 1, "keep");
    try model.delVars(&[_]usize{0});
    try model.updateModel();
    try std.testing.expectEqual(@as(usize, 2), model.genconstr_count);
    try std.testing.expectEqual(@as(usize, 1), model.genconstr_resvar[0]);
    try std.testing.expectEqual(@as(usize, 0), model.genconstr_indices[0]);
    try std.testing.expectEqual(@as(usize, 1), model.genconstr_nvars[0]);
    try std.testing.expectEqual(@as(usize, 0), model.genconstr_indices[1]);
}

test "delQConstrs compacts quadratic constraint storage" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();
    try model.addQConstr(1, &[_]i32{0}, &[_]i32{0}, &[_]f64{1.0}, 0, &.{}, &.{}, .less_equal, 1.0, "q0");
    try model.addQConstr(1, &[_]i32{0}, &[_]i32{0}, &[_]f64{2.0}, 0, &.{}, &.{}, .less_equal, 2.0, "q1");
    var qnz: usize = 0;
    var lnz: usize = 0;
    var qrow: [1]i32 = undefined;
    var qcol: [1]i32 = undefined;
    var qval: [1]f64 = undefined;
    try model.getQConstr(1, &qnz, &qrow, &qcol, &qval, &lnz, &.{}, &.{});
    try std.testing.expectEqual(@as(usize, 1), qnz);
    try std.testing.expectEqual(@as(f64, 2.0), qval[0]);
    try model.delQConstrs(&[_]usize{0});
    try std.testing.expectEqual(@as(usize, 1), model.qconstr_count);
    try std.testing.expectEqual(@as(f64, 2.0), model.qconstr_qval[0]);
    try std.testing.expectEqualStrings("q1", model.qconstr_names[0].?);
}

test "delSOS compacts SOS storage" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "test");
    defer model.deinit();
    try model.addSOS(.sos1, 1, &[_]usize{0}, &[_]f64{1.0}, "s0");
    try model.addSOS(.sos2, 2, &[_]usize{ 0, 1 }, &[_]f64{ 2.0, 3.0 }, "s1");
    try model.delSOS(&[_]usize{0});
    try std.testing.expectEqual(@as(usize, 1), model.sos_count);
    try std.testing.expectEqual(@as(usize, 2), model.sos_indices.len);
    try std.testing.expectEqual(@as(usize, 0), model.sos_begin[0]);
    try std.testing.expectEqual(@as(usize, 2), model.sos_begin[1]);
    try std.testing.expectEqualStrings("s1", model.sos_names[0].?);
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

    try std.testing.expect(!model.interrupted.load(.acquire));
    model.terminate();
    try std.testing.expect(model.interrupted.load(.acquire));
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

test "Model batches coefficient replacement insertion and deletion into one revision" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "coefficient_batch");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, null);
    try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, null);
    try model.addVar(1, &[_]usize{0}, &[_]f64{1.0}, 0.0, 0.0, 1.0, .continuous, null);
    try model.addVar(1, &[_]usize{1}, &[_]f64{2.0}, 0.0, 0.0, 1.0, .continuous, null);
    try model.updateModel();
    _ = try model.matrix.csr(std.testing.allocator);
    const old_revision = model.matrix.matrixRevision();
    const old_starts_ptr = model.matrix.csc().col_starts.ptr;
    const old_rows_ptr = model.matrix.csc().row_indices.ptr;
    const old_values_ptr = model.matrix.csc().values.ptr;

    // Existing-value-only batches mutate the values stream in place. The
    // second write to (0,0) wins without replacing any CSC allocation.
    try model.chgCoeff(0, 0, 3.0);
    try model.chgCoeff(0, 0, 4.0);
    try model.updateModel();
    try std.testing.expectEqual(old_revision + 1, model.matrix.matrixRevision());
    try std.testing.expectEqual(old_starts_ptr, model.matrix.csc().col_starts.ptr);
    try std.testing.expectEqual(old_rows_ptr, model.matrix.csc().row_indices.ptr);
    try std.testing.expectEqual(old_values_ptr, model.matrix.csc().values.ptr);
    try std.testing.expect(model.matrix.csr_cache_revision == null);
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), try model.getCoeff(0, 0), 1e-12);

    // A structural delta batch is still materialized and committed only once.
    _ = try model.matrix.csr(std.testing.allocator);
    const structural_revision = model.matrix.matrixRevision();
    try model.chgCoeff(1, 0, 5.0);
    try model.chgCoeff(1, 1, 0.0);
    try model.updateModel();

    try std.testing.expectEqual(structural_revision + 1, model.matrix.matrixRevision());
    try std.testing.expect(model.matrix.csr_cache_revision == null);
    try std.testing.expectApproxEqAbs(@as(f64, 5.0), try model.getCoeff(1, 0), 1e-12);
    try std.testing.expectEqual(@as(f64, 0.0), try model.getCoeff(1, 1));
    try std.testing.expectEqual(@as(usize, 2), model.matrix.csc().nnz());
    try model.matrix.csc().validate();
}

test "Model coefficient batches match a dense last-write-wins oracle" {
    const num_rows: usize = 13;
    const num_cols: usize = 11;
    const batch_count: usize = 40;
    const changes_per_batch: usize = 80;

    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "coefficient_batch_differential");
    defer model.deinit();

    for (0..num_rows) |_|
        try model.addConstr(0, &.{}, &.{}, .less_equal, 10.0, null);
    for (0..num_cols) |col| {
        const row = col % num_rows;
        try model.addVar(1, &.{row}, &.{1.0}, 0.0, 0.0, 1.0, .continuous, null);
    }
    try model.updateModel();

    var expected = [_]f64{0.0} ** (num_rows * num_cols);
    for (0..num_cols) |col| expected[col * num_rows + col % num_rows] = 1.0;

    var prng = std.Random.DefaultPrng.init(0xc0ef_f1c1_e17d_2026);
    const random = prng.random();
    for (0..batch_count) |_| {
        const revision_before = model.matrix.matrixRevision();
        for (0..changes_per_batch) |_| {
            const row = random.intRangeLessThan(usize, 0, num_rows);
            const col = random.intRangeLessThan(usize, 0, num_cols);
            // Zero exercises deletion and absent-zero no-ops; revisiting a
            // coordinate in the same batch exercises last-write-wins.
            const value: f64 = @floatFromInt(random.intRangeAtMost(i8, -3, 3));
            try model.chgCoeff(row, col, value);
            expected[col * num_rows + row] = value;
        }
        try model.updateModel();

        // A batch may be a semantic no-op, but it must never commit the
        // authoritative matrix more than once.
        try std.testing.expect(model.matrix.matrixRevision() <= revision_before + 1);
        try model.matrix.csc().validate();
        for (0..num_cols) |col| {
            for (0..num_rows) |row| {
                try std.testing.expectEqual(expected[col * num_rows + row], try model.getCoeff(row, col));
            }
        }
    }
}

test "Model edit plan coalesces repeated scalar targets into one committed value" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "scalar_edit_plan");
    defer model.deinit();

    try model.addConstr(0, &.{}, &.{}, .less_equal, 1.0, null);
    try model.addVar(0, &.{}, &.{}, 0.0, 0.0, 10.0, .continuous, null);
    try model.updateModel();
    const revision_before = model.revision;
    const bounds_revision_before = model.revisions.bounds;
    const objective_revision_before = model.revisions.objective;

    try model.chgBounds(0, 1.0, 9.0);
    try model.chgBounds(0, 2.0, 8.0);
    try model.chgObj(0, 3.0);
    try model.chgObj(0, 4.0);
    try model.chgRHS(0, 5.0);
    try model.chgRHS(0, 6.0);
    try model.updateModel();

    try std.testing.expectEqual(@as(f64, 2.0), model.var_lb[0]);
    try std.testing.expectEqual(@as(f64, 8.0), model.var_ub[0]);
    try std.testing.expectEqual(@as(f64, 4.0), model.var_obj[0]);
    try std.testing.expectEqual(@as(f64, 6.0), model.constr_rhs[0]);
    try std.testing.expectEqual(revision_before + 1, model.revision);
    try std.testing.expectEqual(bounds_revision_before + 1, model.revisions.bounds);
    try std.testing.expectEqual(objective_revision_before + 1, model.revisions.objective);
}

test "Model appends nonempty constraint row payload into existing CSC columns" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();
    var model = try Model.init(std.testing.allocator, &env, "append_constraint_row");
    defer model.deinit();

    try model.addVar(0, &.{}, &.{}, 0.0, 0.0, 10.0, .continuous, null);
    try model.addVar(0, &.{}, &.{}, 0.0, 0.0, 10.0, .continuous, null);
    try model.updateModel();
    try model.addConstr(2, &.{ 0, 1 }, &.{ 2.0, -3.0 }, .less_equal, 4.0, null);
    try model.updateModel();

    try std.testing.expectEqual(@as(f64, 2.0), try model.getCoeff(0, 0));
    try std.testing.expectEqual(@as(f64, -3.0), try model.getCoeff(0, 1));
    try std.testing.expectEqualSlices(usize, &.{ 0, 1, 2 }, model.matrix.csc().col_starts);
    try model.matrix.csc().validate();
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
