//! Problem classification for routing, logging, and user queries.
//!
//! The final `ProblemClass` is derived from three orthogonal dimensions and is
//! never set directly — it always reflects the model's actual contents.
//!
//! Classification priority (highest to lowest):
//!   1. Nonlinear (constraint or objective) → NLP / MINLP
//!   2. Quadratic constraint              → QCP / MIQCP
//!   3. Quadratic objective               → QP / MIQP
//!   4. Linear                            → LP / MILP
//!
//! Experimental API. Classification rules are stable but the enum representation
//! may acquire new members as additional problem types are added.

const std = @import("std");

/// Whether the model contains discrete (integer / semi-) variables.
pub const DomainClass = enum {
    continuous,
    mixed_integer,
};

/// Objective function type.
pub const ObjectiveClass = enum {
    linear,
    quadratic,
    nonlinear,
};

/// Constraint type.
pub const ConstraintClass = enum {
    linear,
    quadratic,
    nonlinear,
};

/// Final problem class for dispatcher routing and user queries.
///
/// Derived from model content — never set directly.
pub const ProblemClass = enum {
    lp,
    milp,
    qp,
    miqp,
    qcp,
    miqcp,
    nlp,
    minlp,

    /// Short uppercase string suitable for logging (e.g. "LP", "MILP", "QCP").
    pub fn label(self: ProblemClass) []const u8 {
        return switch (self) {
            .lp => "LP",
            .milp => "MILP",
            .qp => "QP",
            .miqp => "MIQP",
            .qcp => "QCP",
            .miqcp => "MIQCP",
            .nlp => "NLP",
            .minlp => "MINLP",
        };
    }
};

// ── Derived classification ────────────────────────────────────────────────

/// Derive the final problem class from the three orthogonal dimensions.
///
/// Priority: nonlinear > quadratic constraint > quadratic objective > linear.
pub fn classify(
    domain: DomainClass,
    objective: ObjectiveClass,
    constraints: ConstraintClass,
) ProblemClass {
    // Nonlinear dominates everything.
    if (objective == .nonlinear or constraints == .nonlinear) {
        return switch (domain) {
            .continuous => .nlp,
            .mixed_integer => .minlp,
        };
    }
    // Quadratic constraint dominates over quadratic objective.
    if (constraints == .quadratic) {
        return switch (domain) {
            .continuous => .qcp,
            .mixed_integer => .miqcp,
        };
    }
    // Quadratic objective with linear constraints.
    if (objective == .quadratic) {
        return switch (domain) {
            .continuous => .qp,
            .mixed_integer => .miqp,
        };
    }
    // Fully linear.
    return switch (domain) {
        .continuous => .lp,
        .mixed_integer => .milp,
    };
}

// ── Solver capability ─────────────────────────────────────────────────────

/// Declares which problem classes the current solver configuration supports.
///
/// This only indicates solver *readiness*, not model *expressiveness*. A model
/// layer may represent classes the solver cannot yet handle; the dispatcher
/// returns `error.UnsupportedProblemClass` in that case.
///
/// Default: all bits cleared. Solvers set the bits they can handle.
pub const SolverCapability = packed struct(u8) {
    lp: bool = false,
    milp: bool = false,
    qp: bool = false,
    miqp: bool = false,
    qcp: bool = false,
    miqcp: bool = false,
    nlp: bool = false,
    minlp: bool = false,
};

// ── Tests ──────────────────────────────────────────────────────────────────

test "classify all eight combinations" {
    try std.testing.expectEqual(ProblemClass.lp, classify(.continuous, .linear, .linear));
    try std.testing.expectEqual(ProblemClass.milp, classify(.mixed_integer, .linear, .linear));
    try std.testing.expectEqual(ProblemClass.qp, classify(.continuous, .quadratic, .linear));
    try std.testing.expectEqual(ProblemClass.miqp, classify(.mixed_integer, .quadratic, .linear));
    try std.testing.expectEqual(ProblemClass.qcp, classify(.continuous, .quadratic, .quadratic));
    try std.testing.expectEqual(ProblemClass.miqcp, classify(.mixed_integer, .quadratic, .quadratic));
    try std.testing.expectEqual(ProblemClass.nlp, classify(.continuous, .nonlinear, .nonlinear));
    try std.testing.expectEqual(ProblemClass.minlp, classify(.mixed_integer, .nonlinear, .nonlinear));
}

test "classify priority: nonlinear > quadratic constraint > quadratic objective" {
    // Quadratic objective + quadratic constraint → QCP / MIQCP (constraint dominates)
    try std.testing.expectEqual(ProblemClass.qcp, classify(.continuous, .quadratic, .quadratic));
    try std.testing.expectEqual(ProblemClass.miqcp, classify(.mixed_integer, .quadratic, .quadratic));

    // Nonlinear objective → NLP / MINLP (nonlinear beats quadratic)
    try std.testing.expectEqual(ProblemClass.nlp, classify(.continuous, .nonlinear, .linear));
    try std.testing.expectEqual(ProblemClass.minlp, classify(.mixed_integer, .nonlinear, .linear));

    // Nonlinear constraint → NLP / MINLP
    try std.testing.expectEqual(ProblemClass.nlp, classify(.continuous, .linear, .nonlinear));
    try std.testing.expectEqual(ProblemClass.minlp, classify(.mixed_integer, .linear, .nonlinear));

    // Nonlinear beats both quadratic
    try std.testing.expectEqual(ProblemClass.nlp, classify(.continuous, .quadratic, .nonlinear));
    try std.testing.expectEqual(ProblemClass.minlp, classify(.mixed_integer, .quadratic, .nonlinear));
    try std.testing.expectEqual(ProblemClass.nlp, classify(.continuous, .nonlinear, .quadratic));
    try std.testing.expectEqual(ProblemClass.minlp, classify(.mixed_integer, .nonlinear, .quadratic));
}

test "classify linear domain combinations" {
    try std.testing.expectEqual(ProblemClass.lp, classify(.continuous, .linear, .linear));
    try std.testing.expectEqual(ProblemClass.milp, classify(.mixed_integer, .linear, .linear));
}

test "ProblemClass label is uppercase and short" {
    try std.testing.expectEqualStrings("LP", ProblemClass.lp.label());
    try std.testing.expectEqualStrings("MILP", ProblemClass.milp.label());
    try std.testing.expectEqualStrings("QP", ProblemClass.qp.label());
    try std.testing.expectEqualStrings("MIQP", ProblemClass.miqp.label());
    try std.testing.expectEqualStrings("QCP", ProblemClass.qcp.label());
    try std.testing.expectEqualStrings("MIQCP", ProblemClass.miqcp.label());
    try std.testing.expectEqualStrings("NLP", ProblemClass.nlp.label());
    try std.testing.expectEqualStrings("MINLP", ProblemClass.minlp.label());
}

test "SolverCapability defaults to all unsupported" {
    const caps = SolverCapability{};
    try std.testing.expect(!caps.lp);
    try std.testing.expect(!caps.milp);
    try std.testing.expect(!caps.qp);
    try std.testing.expect(!caps.minlp);
}

test "SolverCapability can set individual bits" {
    const caps = SolverCapability{ .lp = true, .milp = true };
    try std.testing.expect(caps.lp);
    try std.testing.expect(caps.milp);
    try std.testing.expect(!caps.qp);
    try std.testing.expect(!caps.minlp);
}
