//! Core type definitions for the model layer.
//!
//! Uses named enums for variable types, constraint senses, optimisation status
//! codes, and basis status values.  The attribute and parameter systems in the
//! sibling modules reference these types by name.

const std = @import("std");

// ── Variable types ──────────────────────────────────────────────────────

/// Variable type: continuous, binary, integer, semi-continuous, semi-integer.
/// Single-character discriminant values follow the convention `'C'`, `'B'`,
/// `'I'`, `'S'`, `'N'`.
pub const VarType = enum(u8) {
    continuous = 'C',
    binary = 'B',
    integer = 'I',
    semicont = 'S',
    semiint = 'N',

    pub const DEFAULT: VarType = .continuous;

    /// Parse from a single-character code.
    pub fn fromCode(code: u8) !VarType {
        return switch (code) {
            'C', 'c' => .continuous,
            'B', 'b' => .binary,
            'I', 'i' => .integer,
            'S', 's' => .semicont,
            'N', 'n' => .semiint,
            else => error.InvalidVarType,
        };
    }
};

// ── Constraint sense ────────────────────────────────────────────────────

/// Constraint sense: less-or-equal, equal, greater-or-equal.
/// Stored as the corresponding ASCII character.
pub const Sense = enum(u8) {
    less_equal = '<',
    equal = '=',
    greater_equal = '>',

    pub const DEFAULT: Sense = .less_equal;

    pub fn fromCode(code: u8) !Sense {
        return switch (code) {
            '<', 'L', 'l' => .less_equal,
            '=', 'E', 'e' => .equal,
            '>', 'G', 'g' => .greater_equal,
            else => error.InvalidSense,
        };
    }
};

// ── Objective sense ─────────────────────────────────────────────────────

/// Objective sense: minimise (1) or maximise (-1).
///
/// The numeric encoding (1 = minimise, -1 = maximise) matches the convention
/// used by several mature solver APIs, making interop straightforward.
pub const ObjectiveSense = enum(i8) {
    minimize = 1,
    maximize = -1,

    pub const DEFAULT: ObjectiveSense = .minimize;

    /// Return the numeric value: 1 for minimise, -1 for maximise.
    pub inline fn toModelSenseValue(self: ObjectiveSense) i8 {
        return @intFromEnum(self);
    }

    pub inline fn fromModelSenseValue(v: i8) !ObjectiveSense {
        return switch (v) {
            1 => .minimize,
            -1 => .maximize,
            else => error.InvalidObjectiveSense,
        };
    }
};

// ── Optimisation status ─────────────────────────────────────────────────

/// Optimisation status codes following the standard convention:
/// LOADED (0), OPTIMAL (1), INFEASIBLE (2), UNBOUNDED (3), etc.
pub const Status = enum(u8) {
    loaded = 0,
    optimal = 1,
    infeasible = 2,
    unbounded = 3,
    inf_or_unbd = 4,
    suboptimal = 5,
    iteration_limit = 6,
    time_limit = 7,
    node_limit = 8,
    solution_limit = 9,
    interrupted = 10,
    numeric = 11,
    in_progress = 12,
    user_obj_limit = 13,

    pub const DEFAULT: Status = .loaded;

    /// Human-readable label for the status.
    pub fn label(self: Status) []const u8 {
        return switch (self) {
            .loaded => "LOADED",
            .optimal => "OPTIMAL",
            .infeasible => "INFEASIBLE",
            .unbounded => "UNBOUNDED",
            .inf_or_unbd => "INF_OR_UNBD",
            .suboptimal => "SUBOPTIMAL",
            .iteration_limit => "ITERATION_LIMIT",
            .time_limit => "TIME_LIMIT",
            .node_limit => "NODE_LIMIT",
            .solution_limit => "SOLUTION_LIMIT",
            .interrupted => "INTERRUPTED",
            .numeric => "NUMERIC",
            .in_progress => "IN_PROGRESS",
            .user_obj_limit => "USER_OBJ_LIMIT",
        };
    }
};

// ── Basis status ────────────────────────────────────────────────────────

/// Simplex basis status for a variable or constraint: basic (0),
/// non-basic at lower bound (-1), non-basic at upper bound (1),
/// super-non-basic (2).
pub const BasisStatus = enum(i8) {
    basic = 0,
    non_basic_lower = -1,
    non_basic_upper = 1,
    super_non_basic = 2,
};

// ── SOS type ────────────────────────────────────────────────────────────

pub const SosType = enum(u8) {
    sos1 = 1,
    sos2 = 2,
};

// ── Error codes (return-value convention) ──────────────────────────────

/// Model-layer errors.  The function-level API uses a return-value convention
/// (`ModelError!void`) rather than propagating Zig error unions on every call.
pub const ModelError = error{
    OutOfMemory,
    InvalidArgument,
    InvalidAttribute,
    InvalidParameter,
    InvalidVarType,
    InvalidSense,
    InvalidObjectiveSense,
    IndexOutOfRange,
    DuplicateName,
    EmptyModel,
    OptimizationInProgress,
    NotOptimized,
    NotInModel,
    FeatureNotAvailable,
    RevisionOverflow,
    IoError,
};

// ── Numeric constants ───────────────────────────────────────────────────

pub const INFINITY: f64 = std.math.inf(f64);
pub const EPSILON: f64 = 1e-9;

// ── Attribute / parameter type tags ─────────────────────────────────────

/// Discriminates the storage type of a model attribute or a parameter.
pub const ValueType = enum(u8) {
    int,
    double,
    string,
    char,
};

// ── Attribute classification ───────────────────────────────────────────

/// Where the attribute lives.
pub const AttributeScope = enum(u8) {
    model,
    variable,
    constraint,
    quadratic,
    sos,
    general,
};

// ── General constraint types ─────────────────────────────────────────────

/// Type of a general (non-linear) constraint.
pub const GenConstrType = enum(u8) {
    abs = 0,
    and_ = 1,
    or_ = 2,
    min = 3,
    max = 4,
    indicator = 5,
    pwl = 6,
    poly = 7,
    exp = 8,
    expa = 9,
    log = 10,
    loga = 11,
    pow = 12,
    sin = 13,
    cos = 14,
    tan = 15,
    logistic = 16,
    norm = 17,
    nl = 18,
};

// ── Feasibility relaxation type ──────────────────────────────────────────

/// How to relax an infeasible model.
pub const FeasRelaxType = enum(u8) {
    /// Relax bounds on continuous variables and objective (LP/QP).
    relax_linear = 0,
    /// Relax bounds, objective, and constraint RHS (LP/QP).
    relax_quad = 1,
    /// For MIP models: also relax variable types.
    relax_mip = 2,
};

// ── Callback context ─────────────────────────────────────────────────────

/// Where the callback is being invoked.
pub const CallbackWhere = enum(u32) {
    polling = 0,
    presolved = 1,
    simplex = 2,
    mip = 3,
    mipsol = 4,
    mipnode = 5,
    message = 6,
    barrier = 7,
    multiobj = 8,
};

/// Signature for user-provided callback functions.
/// Takes the callback location and an opaque user-data pointer.
pub const CallbackFunc = *const fn (cb_where: CallbackWhere, usrstate: ?*anyopaque) void;

// ── Version information ──────────────────────────────────────────────────

/// Version triple (major, minor, technical).
pub const Version = struct {
    major: c_uint,
    minor: c_uint,
    technical: c_uint,
};
