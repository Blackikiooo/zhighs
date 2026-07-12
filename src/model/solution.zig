//! Solver‑internal solution, basis, and solve‑info data structures.
//!
//! These types sit between the solver dispatcher and the optimisation
//! engines.  A [`Solution`] holds the primal–dual vectors, [`Basis`]
//! holds simplex basis statuses, and [`SolveInfo`] is a lightweight
//! POD capturing the outcome metadata.
//!
//! All arrays are heap‑allocated and owned by the struct.  There is no
//! coupling to the Gurobi‑style `Model` — this is pure solver‑internal IR.
//!
//! Experimental API.

const std = @import("std");
const types = @import("types.zig");
const linear_model_module = @import("linear_model.zig");

const BasisStatus = types.BasisStatus;
const Status = types.Status;
const LinearModel = linear_model_module.LinearModel;

// ── Solution ─────────────────────────────────────────────────────────────────

/// Experimental API: owning primal–dual solution vector set.
///
/// | Field        | Length     | Meaning                    |
/// |--------------|------------|----------------------------|
/// | `primal`     | `num_cols` | Primal variable values x   |
/// | `dual`       | `num_rows` | Dual variables π           |
/// | `reduced_cost` | `num_cols` | Reduced costs c − Aᵀπ   |
/// | `slack`      | `num_rows` | Constraint slack Ax − b    |
pub const Solution = struct {
    allocator: std.mem.Allocator,
    primal: []f64,
    dual: []f64,
    reduced_cost: []f64,
    slack: []f64,

    const Self = @This();

    /// Allocate zero‑initialised vectors.
    pub fn init(allocator: std.mem.Allocator, num_cols: usize, num_rows: usize) !Self {
        const primal = try allocator.alloc(f64, num_cols);
        errdefer allocator.free(primal);
        const dual = try allocator.alloc(f64, num_rows);
        errdefer allocator.free(dual);
        const reduced_cost = try allocator.alloc(f64, num_cols);
        errdefer allocator.free(reduced_cost);
        const slack = try allocator.alloc(f64, num_rows);
        errdefer allocator.free(slack);

        @memset(primal, 0.0);
        @memset(dual, 0.0);
        @memset(reduced_cost, 0.0);
        @memset(slack, 0.0);

        return Self{
            .allocator = allocator,
            .primal = primal,
            .dual = dual,
            .reduced_cost = reduced_cost,
            .slack = slack,
        };
    }

    /// Release all owned arrays.
    pub fn deinit(self: *Self) void {
        const a = self.allocator;
        if (self.primal.len > 0) a.free(self.primal);
        if (self.dual.len > 0) a.free(self.dual);
        if (self.reduced_cost.len > 0) a.free(self.reduced_cost);
        if (self.slack.len > 0) a.free(self.slack);
        self.* = undefined;
    }

    /// Compute the linear objective value cᵀx + offset.
    ///
    /// Returns `NaN` when the primal vector is shorter than `model.num_cols`.
    pub fn objectiveValue(self: *const Self, model: *const LinearModel) f64 {
        if (self.primal.len < model.num_cols) return std.math.nan(f64);
        var val = model.objective_offset;
        for (0..model.num_cols) |col| {
            val += model.col_cost[col] * self.primal[col];
        }
        return val;
    }
};

// ── Basis ────────────────────────────────────────────────────────────────────

/// Experimental API: owning simplex basis statuses.
pub const Basis = struct {
    allocator: std.mem.Allocator,
    /// Column (variable) basis statuses.  Length = num_cols.
    col_status: []BasisStatus,
    /// Row (constraint) basis statuses.  Length = num_rows.
    row_status: []BasisStatus,

    const Self = @This();

    /// Allocate and initialise to a default logical basis
    /// (slacks basic, variables non‑basic at lower bound).
    pub fn init(allocator: std.mem.Allocator, num_cols: usize, num_rows: usize) !Self {
        const cs = try allocator.alloc(BasisStatus, num_cols);
        errdefer allocator.free(cs);
        const rs = try allocator.alloc(BasisStatus, num_rows);
        errdefer allocator.free(rs);
        @memset(cs, .non_basic_lower);
        @memset(rs, .basic);
        return Self{ .allocator = allocator, .col_status = cs, .row_status = rs };
    }

    /// Release all owned arrays.
    pub fn deinit(self: *Self) void {
        const a = self.allocator;
        if (self.col_status.len > 0) a.free(self.col_status);
        if (self.row_status.len > 0) a.free(self.row_status);
        self.* = undefined;
    }

    /// Count the number of basic variables.
    pub fn numBasic(self: *const Self) usize {
        var n: usize = 0;
        for (self.col_status) |s| {
            if (s == .basic) n += 1;
        }
        for (self.row_status) |s| {
            if (s == .basic) n += 1;
        }
        return n;
    }
};

// ── SolveInfo ────────────────────────────────────────────────────────────────

/// Experimental API: lightweight solve outcome metadata.
///
/// This is a POD value type — no allocator, no `deinit`.
pub const SolveInfo = struct {
    status: Status,
    obj_val: f64,
    obj_bound: f64,
    iter_count: i64,
    node_count: i64,
};

// ── Tests ────────────────────────────────────────────────────────────────────

test "Solution.init and deinit does not leak" {
    var sol = try Solution.init(std.testing.allocator, 3, 2);
    defer sol.deinit();
    try std.testing.expectEqual(@as(usize, 3), sol.primal.len);
    try std.testing.expectEqual(@as(usize, 2), sol.dual.len);
    try std.testing.expectEqual(@as(usize, 3), sol.reduced_cost.len);
    try std.testing.expectEqual(@as(usize, 2), sol.slack.len);
}

test "Solution.objectiveValue on zero model" {
    var lm = try LinearModel.initEmpty(std.testing.allocator);
    defer lm.deinit();
    var sol = try Solution.init(std.testing.allocator, 0, 0);
    defer sol.deinit();
    try std.testing.expectEqual(@as(f64, 0.0), sol.objectiveValue(&lm));
}

test "Solution.objectiveValue with offset and costs" {
    const Builder = @import("linear_model_builder.zig").LinearModelBuilder;
    var bld = try Builder.init(std.testing.allocator, 1, 2);
    defer bld.deinit();
    try bld.setColCost(0, 3.0);
    try bld.setColCost(1, 5.0);
    bld.setObjectiveOffset(1.0);
    var lm = try bld.freeze();
    defer lm.deinit();

    var sol = try Solution.init(std.testing.allocator, 2, 1);
    defer sol.deinit();
    sol.primal[0] = 2.0;
    sol.primal[1] = 4.0;
    // 3*2 + 5*4 + 1 = 6 + 20 + 1 = 27
    try std.testing.expectApproxEqAbs(@as(f64, 27.0), sol.objectiveValue(&lm), 1e-12);
}

test "Solution.objectiveValue returns NaN on dimension mismatch" {
    var lm = try LinearModel.initEmpty(std.testing.allocator);
    defer lm.deinit();
    var sol = try Solution.init(std.testing.allocator, 1, 0);
    defer sol.deinit();
    // primal.len (1) >= model.num_cols (0) passes, but if we had
    // primal.len < model.num_cols we'd get NaN.
    // For empty model, primal.len (1) > num_cols (0), so no NaN.
    // We want a case where primal.len < num_cols:
    var sol_short = try Solution.init(std.testing.allocator, 0, 0);
    defer sol_short.deinit();
    lm.num_cols = 1;
    try std.testing.expect(std.math.isNan(sol_short.objectiveValue(&lm)));
}

test "Basis.init and deinit does not leak" {
    var b = try Basis.init(std.testing.allocator, 3, 2);
    defer b.deinit();
    try std.testing.expectEqual(@as(usize, 3), b.col_status.len);
    try std.testing.expectEqual(@as(usize, 2), b.row_status.len);
}

test "Basis default initialisation" {
    var b = try Basis.init(std.testing.allocator, 3, 2);
    defer b.deinit();
    try std.testing.expectEqual(BasisStatus.non_basic_lower, b.col_status[0]);
    try std.testing.expectEqual(BasisStatus.basic, b.row_status[0]);
    try std.testing.expectEqual(@as(usize, 2), b.numBasic());
}

test "Basis.numBasic counts correctly" {
    var b = try Basis.init(std.testing.allocator, 2, 1);
    defer b.deinit();
    try std.testing.expectEqual(@as(usize, 1), b.numBasic()); // 0 col basic, 1 row basic
    b.col_status[0] = .basic;
    try std.testing.expectEqual(@as(usize, 2), b.numBasic());
}

test "SolveInfo is a POD value" {
    const info = SolveInfo{
        .status = .optimal,
        .obj_val = 42.0,
        .obj_bound = 0.0,
        .iter_count = 10,
        .node_count = 5,
    };
    try std.testing.expectEqual(Status.optimal, info.status);
    try std.testing.expectEqual(@as(f64, 42.0), info.obj_val);
}
