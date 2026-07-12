//! Attribute name enum and metadata.
//!
//! Problem and solution data are accessed through the uniform `get*Attr` /
//! `set*Attr` interface keyed by string names.  This module provides an
//! `Attr` enum with name resolution, metadata queries, and a runtime
//! `lookup` helper for callers that receive attribute names dynamically.

const std = @import("std");
const types = @import("types.zig");

const ValueType = types.ValueType;
const AttributeScope = types.AttributeScope;

/// Describes one attribute: its name, Zig type, scope, and whether it is
/// user-settable or read-only.
pub const AttributeInfo = struct {
    name: []const u8,
    value_type: ValueType,
    scope: AttributeScope,
    settable: bool,
};

// ══════════════════════════════════════════════════════════════════════════
//  Attr — strongly typed attribute name
// ══════════════════════════════════════════════════════════════════════════

/// Strongly typed attribute identifier.
///
/// Convert to a string name via `.name()`, look up metadata via `.info()`,
/// or parse a runtime string via `Attr.fromName(...)`.
///
/// ```zig
/// // Compile-time lookups
/// try model.getIntAttr(Attr.num_vars);
/// try model.getDblAttrElement("LB", i);
///
/// // Runtime lookup
/// const a = Attr.fromName(user_input) orelse return error.InvalidAttribute;
/// ```
pub const Attr = enum {
    // ── Model scalar attributes ─────────────────────────────────────────

    /// Number of variables in the model.
    num_vars,
    /// Number of linear constraints.
    num_constrs,
    /// Number of non‑zero coefficients in the constraint matrix.
    num_nz,
    /// Optimisation status code (`Status` enum).
    status,
    /// Objective value of the current solution.
    obj_val,
    /// Best available objective bound.
    obj_bound,
    /// Objective constant (offset) term `c₀`.
    obj_con,
    /// Model objective sense: 1 = minimise, -1 = maximise.
    model_sense,
    /// Whether the model contains any integer variables (`0`/`1`).
    is_mip,
    /// Number of simplex iterations performed.
    iter_count,
    /// Number of branch‑and‑bound nodes explored.
    node_count,
    /// Number of barrier iterations performed.
    bar_iter_count,
    /// Number of non‑zero entries in the quadratic objective matrix `Q`.
    num_qnz,
    /// Number of SOS constraints.
    num_sos,
    /// Number of solutions stored in the solution pool.
    sol_count,
    /// Model name string.
    model_name,
    /// Human‑readable label for the current optimisation status.
    status_label,

    // ── Variable attributes (per‑element) ──────────────────────────────

    /// Variable lower bound.
    lb,
    /// Variable upper bound.
    ub,
    /// Variable objective (cost) coefficient.
    obj,
    /// Variable type (`VarType` enum: C/B/I/S/N).
    v_type,
    /// Primal solution vector `x`.
    x,
    /// Reduced cost vector `c − Aᵀπ`.
    rc,
    /// Variable name string.
    var_name,
    /// MIP start vector (Start attribute).
    start,
    /// Primal start vector (PStart, warm‑start hint).
    p_start,
    /// Dual start vector (DStart, warm‑start hint).
    d_start,
    /// Variable basis status (`BasisStatus` enum).
    v_basis,

    // ── Constraint attributes (per‑element) ────────────────────────────

    /// Constraint sense (`Sense` enum: `<`, `=`, `>`).
    sense,
    /// Constraint right‑hand side value.
    rhs,
    /// Constraint name string.
    constr_name,
    /// Dual solution (pi) vector `π`.
    pi,
    /// Constraint slack vector `Ax − b`.
    slack,
    /// Constraint basis status (`BasisStatus` enum).
    c_basis,

    // ── IIS attributes ─────────────────────────────────────────────────

    /// Whether the computed IIS is guaranteed minimal (`0`/`1`).
    iis_minimal,
    /// Whether the variable lower bound participates in the IIS.
    iis_lb,
    /// Whether the variable upper bound participates in the IIS.
    iis_ub,
    /// Whether the constraint sense participates in the IIS.
    iis_sense,
    /// Whether the constraint RHS participates in the IIS.
    iis_rhs,
    /// Whether a quadratic constraint participates in the IIS.
    iis_qconstr,
    /// Whether a general constraint participates in the IIS.
    iis_genconstr,
    /// Whether an SOS constraint participates in the IIS.
    iis_sos,

    // ── MIP pool attributes ─────────────────────────────────────────────

    /// Objective value of a solution in the MIP solution pool.
    pool_obj_val,
    /// Objective bound for the MIP solution pool.
    pool_obj_bound,
    /// Number of solutions stored in the MIP pool.
    pool_solutions,
    /// MIP search mode: 1 = traditional, 2 = solution pool, 3 = sos.
    pool_search_mode,

    // ── MIP gap / quality ───────────────────────────────────────────────

    /// Relative MIP optimality gap `|best − bound| / |best|`.
    mip_gap,
    /// Absolute MIP optimality gap `|best − bound|`.
    mip_gap_abs,

    // ── Sensitivity ─────────────────────────────────────────────────────

    /// Lowest objective coefficient for which the basis remains optimal.
    sa_obj_low,
    /// Highest objective coefficient for which the basis remains optimal.
    sa_obj_up,
    /// Lowest RHS value for which the basis remains optimal.
    sa_rhs_low,
    /// Highest RHS value for which the basis remains optimal.
    sa_rhs_up,

    // ── Multi‑objective ─────────────────────────────────────────────────

    /// Objective value for the `n`‑th sub‑objective.
    obj_n_val,
    /// Weight assigned to the `n`‑th sub‑objective.
    obj_n_weight,
    /// Priority level of the `n`‑th sub‑objective (lower = higher priority).
    obj_n_priority,
    /// Relative tolerance for the `n`‑th sub‑objective.
    obj_n_rel_tol,
    /// Absolute tolerance for the `n`‑th sub‑objective.
    obj_n_abs_tol,
    /// Name for the `n`‑th sub‑objective.
    obj_n_name,

    // ── Count / misc ────────────────────────────────────────────────────

    /// Number of binary variables (`vtype = B`).
    num_bin_vars,
    /// Number of integer variables (`vtype = I`, excluding binary).
    num_int_vars,
    /// Number of general constraints.
    num_gen_constrs,
    /// Number of quadratic constraints.
    num_q_constrs,
    /// Number of piecewise‑linear objective functions.
    num_pwl_obj,

    // ──────────────────────────────────────────────────────────────────────
    //  Methods
    // ──────────────────────────────────────────────────────────────────────

    /// Return the wire-format attribute name.
    ///
    /// Derived automatically from the snake_case variant name:
    /// segments that are known abbreviations or single letters are
    /// uppercased; all other segments are PascalCased.
    pub fn name(self: Attr) []const u8 {
        return switch (self) {
            inline else => |attr| comptime blk: {
                const tag = @tagName(attr);
                const out = tagToAttrName(tag);
                break :blk out;
            },
        };
    }

    /// Convert a snake_case tag to the canonical attribute name.
    /// Known all-caps abbreviations: lb, ub, rc, sa, nz, rhs, iis, mip, qnz, sos, pwl.
    /// Single-letter segments are always uppercased.  Everything else is PascalCased.
    fn tagToAttrName(comptime tag: []const u8) []const u8 {
        @setEvalBranchQuota(5000);
        const segs = comptime calcSegs: {
            var n: usize = 1;
            for (tag) |c| {
                if (c == '_') n += 1;
            }
            break :calcSegs n;
        };
        const out_len = tag.len - (segs - 1);
        const result: [out_len]u8 = comptime calc: {
            var arr: [out_len]u8 = undefined;
            var pos: usize = 0;
            var i: usize = 0;
            while (i < tag.len) {
                const start = i;
                while (i < tag.len and tag[i] != '_') : (i += 1) {}
                const seg = tag[start..i];
                const caps = seg.len == 1 or isAbbreviation(seg);
                for (seg, 0..) |c, j| {
                    const uc = if (c >= 'a' and c <= 'z') c - 32 else c;
                    const lc = if (c >= 'A' and c <= 'Z') c + 32 else c;
                    arr[pos] = if (caps) uc else if (j == 0) uc else lc;
                    pos += 1;
                }
                if (i < tag.len) i += 1;
            }
            break :calc arr;
        };
        return &result;
    }

    /// Return true if `seg` is a known all-caps abbreviation.
    fn isAbbreviation(comptime seg: []const u8) bool {
        return comptime switch (seg.len) {
            2 => switch (seg[0]) {
                'l' => seg[1] == 'b',
                'u' => seg[1] == 'b',
                'r' => seg[1] == 'c',
                's' => seg[1] == 'a',
                'n' => seg[1] == 'z',
                else => false,
            },
            3 => switch (seg[0]) {
                'r' => seg[1] == 'h' and seg[2] == 's',
                'i' => seg[1] == 'i' and seg[2] == 's',
                'm' => seg[1] == 'i' and seg[2] == 'p',
                'q' => seg[1] == 'n' and seg[2] == 'z',
                's' => seg[1] == 'o' and seg[2] == 's',
                'p' => seg[1] == 'w' and seg[2] == 'l',
                else => false,
            },
            else => false,
        };
    }

    /// Return metadata for this attribute.
    pub fn info(self: Attr) AttributeInfo {
        return switch (self) {
            // Model scalars
            .num_vars => makeInfo(.int, .model, false),
            .num_constrs => makeInfo(.int, .model, false),
            .num_nz => makeInfo(.int, .model, false),
            .status => makeInfo(.int, .model, false),
            .obj_val => makeInfo(.double, .model, false),
            .obj_bound => makeInfo(.double, .model, false),
            .obj_con => makeInfo(.double, .model, true),
            .model_sense => makeInfo(.int, .model, true),
            .is_mip => makeInfo(.int, .model, false),
            .iter_count => makeInfo(.int, .model, false),
            .node_count => makeInfo(.int, .model, false),
            .bar_iter_count => makeInfo(.int, .model, false),
            .num_qnz => makeInfo(.int, .model, false),
            .num_sos => makeInfo(.int, .model, false),
            .sol_count => makeInfo(.int, .model, false),
            .model_name => makeInfo(.string, .model, true),
            .status_label => makeInfo(.string, .model, false),
            // Variable
            .lb => makeInfo(.double, .variable, true),
            .ub => makeInfo(.double, .variable, true),
            .obj => makeInfo(.double, .variable, true),
            .v_type => makeInfo(.char, .variable, true),
            .x => makeInfo(.double, .variable, false),
            .rc => makeInfo(.double, .variable, false),
            .var_name => makeInfo(.string, .variable, true),
            .v_basis => makeInfo(.int, .variable, true),
            .start => makeInfo(.double, .variable, true),
            .p_start => makeInfo(.double, .variable, true),
            .d_start => makeInfo(.double, .variable, true),
            // Constraint
            .sense => makeInfo(.char, .constraint, true),
            .rhs => makeInfo(.double, .constraint, true),
            .constr_name => makeInfo(.string, .constraint, true),
            .pi => makeInfo(.double, .constraint, false),
            .slack => makeInfo(.double, .constraint, false),
            .c_basis => makeInfo(.int, .constraint, true),
            // IIS
            .iis_minimal => makeInfo(.int, .model, false),
            .iis_lb => makeInfo(.int, .variable, false),
            .iis_ub => makeInfo(.int, .variable, false),
            .iis_sense => makeInfo(.int, .constraint, false),
            .iis_rhs => makeInfo(.int, .constraint, false),
            .iis_qconstr => makeInfo(.int, .quadratic, false),
            .iis_genconstr => makeInfo(.int, .general, false),
            .iis_sos => makeInfo(.int, .sos, false),
            // MIP pool
            .pool_obj_val => makeInfo(.double, .model, false),
            .pool_obj_bound => makeInfo(.double, .model, false),
            .pool_solutions => makeInfo(.int, .model, true),
            .pool_search_mode => makeInfo(.int, .model, true),
            // MIP gap
            .mip_gap => makeInfo(.double, .model, false),
            .mip_gap_abs => makeInfo(.double, .model, false),
            // Sensitivity
            .sa_obj_low => makeInfo(.double, .variable, false),
            .sa_obj_up => makeInfo(.double, .variable, false),
            .sa_rhs_low => makeInfo(.double, .constraint, false),
            .sa_rhs_up => makeInfo(.double, .constraint, false),
            // Multi-objective
            .obj_n_val => makeInfo(.double, .model, false),
            .obj_n_weight => makeInfo(.double, .model, true),
            .obj_n_priority => makeInfo(.int, .model, true),
            .obj_n_rel_tol => makeInfo(.double, .model, true),
            .obj_n_abs_tol => makeInfo(.double, .model, true),
            .obj_n_name => makeInfo(.string, .model, true),
            // Count
            .num_bin_vars => makeInfo(.int, .model, false),
            .num_int_vars => makeInfo(.int, .model, false),
            .num_gen_constrs => makeInfo(.int, .model, false),
            .num_q_constrs => makeInfo(.int, .model, false),
            .num_pwl_obj => makeInfo(.int, .model, false),
        };
    }

    /// Parse an `Attr` from a wire-format string name.
    /// Returns `null` if the name is not recognised.
    pub inline fn fromName(s: []const u8) ?Attr {
        return name_to_attr.get(s);
    }
};

/// Comptime-constructed map from wire-format name → Attr.
const name_to_attr = blk: {
    const variants = std.meta.fields(Attr);
    var entries: [variants.len]struct { []const u8, Attr } = undefined;
    for (variants, 0..) |f, i| {
        entries[i] = .{ @field(Attr, f.name).name(), @field(Attr, f.name) };
    }
    break :blk std.StaticStringMap(Attr).initComptime(entries);
};

// ── Helpers ─────────────────────────────────────────────────────────────

/// Shorthand to build an `AttributeInfo` at comptime.
inline fn makeInfo(value_type: ValueType, scope: AttributeScope, settable: bool) AttributeInfo {
    return .{ .name = "", .value_type = value_type, .scope = scope, .settable = settable };
}

// ── Runtime lookup (backward-compat) ────────────────────────────────────

/// Return the `AttributeInfo` for a known attribute name, or `null`.
///
/// This is the runtime-friendly entry point for code that receives attribute
/// names as strings.  For compile-time-known names prefer `Attr.fromName(...)`
/// or using enum variants directly.
pub fn lookup(name: []const u8) ?AttributeInfo {
    const a = Attr.fromName(name) orelse return null;
    var info_rec = a.info();
    info_rec.name = a.name();
    return info_rec;
}

// ── Tests ───────────────────────────────────────────────────────────────

test "Attr.name returns expected string" {
    try std.testing.expectEqualStrings("NumVars", Attr.num_vars.name());
    try std.testing.expectEqualStrings("LB", Attr.lb.name());
    try std.testing.expectEqualStrings("ObjVal", Attr.obj_val.name());
    try std.testing.expectEqualStrings("VarName", Attr.var_name.name());
    try std.testing.expectEqualStrings("ConstrName", Attr.constr_name.name());
    try std.testing.expectEqualStrings("ModelSense", Attr.model_sense.name());
    try std.testing.expectEqualStrings("IISMinimal", Attr.iis_minimal.name());
    try std.testing.expectEqualStrings("PoolObjVal", Attr.pool_obj_val.name());
    try std.testing.expectEqualStrings("MIPGap", Attr.mip_gap.name());
    try std.testing.expectEqualStrings("NumBinVars", Attr.num_bin_vars.name());
}

test "Attr.info returns correct metadata" {
    {
        const md = Attr.lb.info();
        try std.testing.expect(md.value_type == .double);
        try std.testing.expect(md.scope == .variable);
        try std.testing.expect(md.settable);
    }
    {
        const md = Attr.x.info();
        try std.testing.expect(md.value_type == .double);
        try std.testing.expect(md.scope == .variable);
        try std.testing.expect(!md.settable);
    }
    {
        const md = Attr.num_vars.info();
        try std.testing.expect(md.value_type == .int);
        try std.testing.expect(md.scope == .model);
        try std.testing.expect(!md.settable);
    }
}

test "Attr.fromName round-trips known names" {
    try std.testing.expectEqual(Attr.num_vars, Attr.fromName("NumVars").?);
    try std.testing.expectEqual(Attr.lb, Attr.fromName("LB").?);
    try std.testing.expectEqual(Attr.obj_val, Attr.fromName("ObjVal").?);
    try std.testing.expectEqual(Attr.constr_name, Attr.fromName("ConstrName").?);
    try std.testing.expectEqual(Attr.iis_minimal, Attr.fromName("IISMinimal").?);
}

test "Attr.fromName returns null for unknown names" {
    try std.testing.expect(Attr.fromName("NonExistent") == null);
}

test "lookup returns info for known attributes" {
    try std.testing.expect(lookup("NumVars") != null);
    try std.testing.expect(lookup("LB") != null);
    try std.testing.expect(lookup("ObjVal") != null);
    try std.testing.expect(lookup("VarName") != null);
    try std.testing.expect(lookup("ConstrName") != null);
}

test "lookup returns null for unknown attributes" {
    try std.testing.expect(lookup("NonExistent") == null);
}

test "Attr.info is correct for LB" {
    const md = Attr.lb.info();
    try std.testing.expect(md.value_type == .double);
    try std.testing.expect(md.scope == .variable);
    try std.testing.expect(md.settable);
}

test "Attr.info is correct for X" {
    const md = Attr.x.info();
    try std.testing.expect(md.value_type == .double);
    try std.testing.expect(md.scope == .variable);
    try std.testing.expect(!md.settable);
}
