//! SCIP-style plugin contracts for solver components.
//!
//! Each contract is a fat pointer (instance pointer + vtable) in the style of
//! `std.mem.Allocator`, so concrete plugins can be plain structs living
//! anywhere, registered at runtime, and dispatched without inheritance or
//! allocation. Use the comptime `bind*` adapters to turn a struct that declares
//! the conventional method set into the corresponding contract value.

const std = @import("std");
const solver = @import("solver");

/// Integrality classification of a structural column in a mixed-integer model.
pub const Integrality = enum { continuous, integer, binary, semi_continuous, semi_integer };

/// A mixed-integer problem described by its LP relaxation plus column
/// integrality. Borrowed, never owned, by every plugin context.
pub const ProblemView = struct {
    relaxation: solver.LpProblemView,
    integrality: []const Integrality,

    pub fn validate(self: ProblemView) !void {
        try self.relaxation.validate();
        if (self.integrality.len != self.relaxation.num_cols) return error.DimensionMismatch;
    }
};

/// Which disjunction the engine must apply to resolve a branching candidate.
pub const BranchingKind = enum { integer, semi_activation };

/// A column the engine should branch on.
pub const BranchingCandidate = struct {
    column: usize,
    kind: BranchingKind,
};

/// Read-only state the engine hands to branching rules when the node LP is
/// fractional. `active_lower`/`active_upper` are the root-level nonzero-branch
/// bounds of semi-continuous / semi-integer columns; both are constant across
/// the search.
pub const BranchingContext = struct {
    problem: ProblemView,
    node_lower: []const f64,
    node_upper: []const f64,
    lp_primal: []const f64,
    active_lower: []const f64,
    active_upper: []const f64,
    tolerance: f64,
};

/// Branching-rule contract. Implementations inspect the fractional LP solution
/// and return the column to branch on, or null to defer to lower-priority
/// rules.
pub const BranchingRule = struct {
    ptr: *anyopaque,
    vtable: *const VTable,

    pub const VTable = struct {
        name: []const u8,
        select: *const fn (ptr: *anyopaque, ctx: BranchingContext) ?BranchingCandidate,
    };

    pub fn name(self: BranchingRule) []const u8 {
        return self.vtable.name;
    }

    pub fn select(self: BranchingRule, ctx: BranchingContext) ?BranchingCandidate {
        return self.vtable.select(self.ptr, ctx);
    }
};

/// State a primal heuristic may use to construct a candidate solution.
pub const HeuristicContext = struct {
    /// Allocator for the returned candidate; the engine frees it after
    /// validation, so heuristics never manage incumbent memory themselves.
    allocator: std.mem.Allocator,
    problem: ProblemView,
    node_lower: []const f64,
    node_upper: []const f64,
    lp_primal: []const f64,
    active_lower: []const f64,
    active_upper: []const f64,
    tolerance: f64,
};

/// Primal-heuristic contract. `trySolution` returns an owned candidate solution
/// allocated from `ctx.allocator`, or null when no candidate was found. The
/// engine independently validates feasibility and objective improvement, so a
/// heuristic may return aggressive guesses at no correctness risk.
pub const PrimalHeuristic = struct {
    ptr: *anyopaque,
    vtable: *const VTable,

    pub const VTable = struct {
        name: []const u8,
        trySolution: *const fn (ptr: *anyopaque, ctx: HeuristicContext) std.mem.Allocator.Error!?[]f64,
    };

    pub fn name(self: PrimalHeuristic) []const u8 {
        return self.vtable.name;
    }

    pub fn trySolution(self: PrimalHeuristic, ctx: HeuristicContext) std.mem.Allocator.Error!?[]f64 {
        return self.vtable.trySolution(self.ptr, ctx);
    }
};

/// One open node as seen by node selectors: its dual bound (null before the
/// parent relaxation was solved) and its depth in the search tree.
pub const OpenNode = struct {
    bound: ?f64,
    depth: u32,
};

/// Read-only snapshot of the open-node set handed to node selectors.
pub const NodeSelectionContext = struct {
    objective_sense: solver.LpObjectiveSense,
    open: []const OpenNode,
};

/// Node-selector contract. Returns the index into `ctx.open` of the node the
/// engine should process next, or null when no node is selectable.
pub const NodeSelector = struct {
    ptr: *anyopaque,
    vtable: *const VTable,

    pub const VTable = struct {
        name: []const u8,
        select: *const fn (ptr: *anyopaque, ctx: NodeSelectionContext) ?usize,
    };

    pub fn name(self: NodeSelector) []const u8 {
        return self.vtable.name;
    }

    pub fn select(self: NodeSelector, ctx: NodeSelectionContext) ?usize {
        return self.vtable.select(self.ptr, ctx);
    }
};

/// Adapt an instance into a [`BranchingRule`]. `T` must declare
/// `pub const name` and
/// `pub fn select(self: *T, ctx: BranchingContext) ?BranchingCandidate`.
/// The instance must outlive every call made through the returned rule.
pub fn bindBranching(comptime T: type, instance: *T) BranchingRule {
    const Adapter = struct {
        fn select(ptr: *anyopaque, ctx: BranchingContext) ?BranchingCandidate {
            const self: *T = @ptrCast(@alignCast(ptr));
            return self.select(ctx);
        }
        const vtable = BranchingRule.VTable{ .name = T.name, .select = select };
    };
    return .{ .ptr = instance, .vtable = &Adapter.vtable };
}

/// Adapt an instance into a [`PrimalHeuristic`]. `T` must declare
/// `pub const name` and
/// `pub fn trySolution(self: *T, ctx: HeuristicContext) std.mem.Allocator.Error!?[]f64`.
pub fn bindHeuristic(comptime T: type, instance: *T) PrimalHeuristic {
    const Adapter = struct {
        fn trySolution(ptr: *anyopaque, ctx: HeuristicContext) std.mem.Allocator.Error!?[]f64 {
            const self: *T = @ptrCast(@alignCast(ptr));
            return self.trySolution(ctx);
        }
        const vtable = PrimalHeuristic.VTable{ .name = T.name, .trySolution = trySolution };
    };
    return .{ .ptr = instance, .vtable = &Adapter.vtable };
}

/// Adapt an instance into a [`NodeSelector`]. `T` must declare
/// `pub const name` and
/// `pub fn select(self: *T, ctx: NodeSelectionContext) ?usize`.
pub fn bindNodeSelector(comptime T: type, instance: *T) NodeSelector {
    const Adapter = struct {
        fn select(ptr: *anyopaque, ctx: NodeSelectionContext) ?usize {
            const self: *T = @ptrCast(@alignCast(ptr));
            return self.select(ctx);
        }
        const vtable = NodeSelector.VTable{ .name = T.name, .select = select };
    };
    return .{ .ptr = instance, .vtable = &Adapter.vtable };
}

test "bind adapters forward calls to plugin instances" {
    const Recorder = struct {
        calls: usize = 0,
        pub const name = "recorder";
        pub fn select(self: *@This(), _: BranchingContext) ?BranchingCandidate {
            self.calls += 1;
            return .{ .column = 7, .kind = .integer };
        }
    };
    var recorder = Recorder{};
    const rule = bindBranching(Recorder, &recorder);
    try std.testing.expectEqualStrings("recorder", rule.name());
    const candidate = rule.select(.{
        .problem = .{
            .relaxation = undefined,
            .integrality = &.{},
        },
        .node_lower = &.{},
        .node_upper = &.{},
        .lp_primal = &.{},
        .active_lower = &.{},
        .active_upper = &.{},
        .tolerance = 1e-6,
    });
    try std.testing.expectEqual(@as(usize, 7), candidate.?.column);
    try std.testing.expectEqual(@as(usize, 1), recorder.calls);

    const NopHeuristic = struct {
        pub const name = "nop";
        pub fn trySolution(_: *@This(), _: HeuristicContext) std.mem.Allocator.Error!?[]f64 {
            return null;
        }
    };
    var nop = NopHeuristic{};
    const heuristic = bindHeuristic(NopHeuristic, &nop);
    try std.testing.expectEqualStrings("nop", heuristic.name());
    try std.testing.expectEqual(@as(?[]f64, null), try heuristic.trySolution(.{
        .allocator = std.testing.allocator,
        .problem = .{ .relaxation = undefined, .integrality = &.{} },
        .node_lower = &.{},
        .node_upper = &.{},
        .lp_primal = &.{},
        .active_lower = &.{},
        .active_upper = &.{},
        .tolerance = 1e-6,
    }));

    const FirstNode = struct {
        pub const name = "first";
        pub fn select(_: *@This(), ctx: NodeSelectionContext) ?usize {
            return if (ctx.open.len == 0) null else 0;
        }
    };
    var first = FirstNode{};
    const selector = bindNodeSelector(FirstNode, &first);
    try std.testing.expectEqualStrings("first", selector.name());
    try std.testing.expectEqual(@as(?usize, 0), selector.select(.{ .objective_sense = .minimize, .open = &.{.{ .bound = null, .depth = 0 }} }));
    try std.testing.expectEqual(@as(?usize, null), selector.select(.{ .objective_sense = .minimize, .open = &.{} }));
}
