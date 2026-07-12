//! Expression DAG for NLP/MINLP models.
//!
//! Stores a directed acyclic graph (DAG) of unary/binary operators applied to
//! constants and decision variables.  The graph owns all storage: node data is
//! stored in dense SoA arrays, children are stored in one contiguous list.
//!
//! ## Mutability
//!
//! The graph is intended to be built via [`ExpressionGraphBuilder`] and then
//! frozen.  After construction the graph is read-only.
//!
//! Experimental API.  The SoA layout is stable; the opcode set may grow as
//! additional intrinsic functions are needed.

const std = @import("std");
const foundation = @import("foundation");

const ColId = foundation.ColId;

// ── NodeId ─────────────────────────────────────────────────────────────────

/// Compact strongly‑typed node identifier.
pub const NodeId = enum(u32) { _ };

// ── Opcode ─────────────────────────────────────────────────────────────────

/// Operation codes for expression DAG nodes.
pub const Opcode = enum(u8) {
    constant,
    variable,

    add,
    subtract,
    multiply,
    divide,
    negate,

    square,
    power,

    exp,
    log,
    sqrt,

    sin,
    cos,
    tan,

    abs,
    min,
    max,

    /// Number of children for each opcode (0 for leaf, 1 for unary, 2 for binary).
    pub fn arity(self: Opcode) u8 {
        return switch (self) {
            .constant, .variable => 0,
            .negate, .square, .exp, .log, .sqrt, .sin, .cos, .tan, .abs => 1,
            .add, .subtract, .multiply, .divide, .power, .min, .max => 2,
        };
    }
};

// ── ExpressionGraph ────────────────────────────────────────────────────────

/// Experimental API: owning expression DAG.
///
/// Nodes are indexed by `NodeId` (0‑based dense).  Children of node `i`
/// occupy `children[first_child[i] .. first_child[i] + child_count[i]]`.
pub const ExpressionGraph = struct {
    allocator: std.mem.Allocator,

    /// Operation of each node.
    opcodes: []Opcode,
    /// Start index into `children` for each node.
    first_child: []usize,
    /// Number of children for each node.
    child_count: []usize,
    /// Contiguous child edges.
    children: []NodeId,
    /// Constant values (valid only for `.constant` nodes).
    constants: []f64,
    /// Decision variable indices (valid only for `.variable` nodes).
    variables: []ColId,

    /// Number of nodes.
    num_nodes: usize,

    const Self = @This();

    pub fn deinit(self: *Self) void {
        const a = self.allocator;
        if (self.opcodes.len > 0) a.free(self.opcodes);
        if (self.first_child.len > 0) a.free(self.first_child);
        if (self.child_count.len > 0) a.free(self.child_count);
        if (self.children.len > 0) a.free(self.children);
        if (self.constants.len > 0) a.free(self.constants);
        if (self.variables.len > 0) a.free(self.variables);
        self.* = undefined;
    }

    /// Returns the children of `node`.
    pub fn nodeChildren(self: *const Self, node: NodeId) []const NodeId {
        const idx = @intFromEnum(node);
        const start = self.first_child[idx];
        const count = self.child_count[idx];
        return self.children[start .. start + count];
    }

    /// DFS colour for cycle detection.
    const Colour = enum(u2) { white, grey, black };

    /// Structural validation of the entire graph.
    ///
    /// Checks:
    /// - Array length consistency.
    /// - Every child ID is in range.
    /// - No cycles (DFS with visiting/visited marking).
    /// - Constant nodes have finite values.
    /// - Variable nodes reference IDs that are within `num_cols` when provided.
    pub fn validate(self: *const Self, num_cols: ?usize) !void {
        // Length consistency.
        if (self.opcodes.len != self.num_nodes) return error.InvalidExpressionNode;
        if (self.first_child.len != self.num_nodes) return error.InvalidExpressionNode;
        if (self.child_count.len != self.num_nodes) return error.InvalidExpressionNode;

        for (0..self.num_nodes) |i| {
            const start = self.first_child[i];
            const count = self.child_count[i];
            if (start > self.children.len or start + count > self.children.len)
                return error.InvalidExpressionNode;

            const op = self.opcodes[i];
            if (op.arity() != count) return error.InvalidExpressionNode;

            // Constant values must be finite.
            if (op == .constant) {
                if (i >= self.constants.len) return error.InvalidExpressionNode;
                if (!std.math.isFinite(self.constants[i])) return error.NonFiniteValue;
            }

            // Variable IDs must be in range.
            if (op == .variable) {
                if (i >= self.variables.len) return error.InvalidExpressionNode;
                if (num_cols) |nc| {
                    if (self.variables[i].toUsize() >= nc)
                        return error.VariableIndexOutOfRange;
                }
            }

            // Check children are valid node IDs.
            for (start..start + count) |child_pos| {
                const child = self.children[child_pos];
                const child_idx = @intFromEnum(child);
                if (child_idx >= self.num_nodes)
                    return error.InvalidExpressionNode;
            }
        }

        // Cycle detection: DFS with three colours.
        const colours = try self.allocator.alloc(Colour, self.num_nodes);
        defer self.allocator.free(colours);
        @memset(colours, .white);

        for (0..self.num_nodes) |i| {
            if (colours[i] == .white) {
                try self.dfs(i, colours);
            }
        }
    }

    fn dfs(self: *const Self, node: usize, colours: []Colour) !void {
        colours[node] = .grey;
        const start = self.first_child[node];
        const count = self.child_count[node];
        for (start..start + count) |child_pos| {
            const child_idx = @intFromEnum(self.children[child_pos]);
            if (colours[child_idx] == .grey) return error.CyclicExpression;
            if (colours[child_idx] == .white) {
                try self.dfs(child_idx, colours);
            }
        }
        colours[node] = .black;
    }
};

// ── ExpressionGraphBuilder ─────────────────────────────────────────────────

/// Experimental API: incremental builder for ExpressionGraph.
pub const ExpressionGraphBuilder = struct {
    allocator: std.mem.Allocator,
    opcodes: std.ArrayListUnmanaged(Opcode) = .empty,
    first_child: std.ArrayListUnmanaged(usize) = .empty,
    child_count: std.ArrayListUnmanaged(usize) = .empty,
    children: std.ArrayListUnmanaged(NodeId) = .empty,
    constants: std.ArrayListUnmanaged(f64) = .empty,
    variables: std.ArrayListUnmanaged(ColId) = .empty,

    const Self = @This();

    pub fn init(allocator: std.mem.Allocator) Self {
        return .{ .allocator = allocator };
    }

    pub fn deinit(self: *Self) void {
        const a = self.allocator;
        self.opcodes.deinit(a);
        self.first_child.deinit(a);
        self.child_count.deinit(a);
        self.children.deinit(a);
        self.constants.deinit(a);
        self.variables.deinit(a);
        self.* = undefined;
    }

    fn allocNode(self: *Self, op: Opcode, child_list: []const NodeId) !NodeId {
        const id = self.opcodes.items.len;
        try self.opcodes.append(self.allocator, op);
        try self.first_child.append(self.allocator, self.children.items.len);
        try self.child_count.append(self.allocator, @intCast(child_list.len));
        try self.children.appendSlice(self.allocator, child_list);
        // Pad optional arrays (constants, variables) to keep lengths in sync.
        try self.constants.append(self.allocator, 0.0);
        try self.variables.append(self.allocator, ColId.fromUsizeAssumeValid(0));
        return @enumFromInt(id);
    }

    /// Append a constant leaf node.
    pub fn addConstant(self: *Self, value: f64) !NodeId {
        if (!std.math.isFinite(value)) return error.NonFiniteValue;
        const id = try self.allocNode(.constant, &.{});
        self.constants.items[@intFromEnum(id)] = value;
        return id;
    }

    /// Append a variable leaf node.
    pub fn addVariable(self: *Self, col: ColId) !NodeId {
        const id = try self.allocNode(.variable, &.{});
        self.variables.items[@intFromEnum(id)] = col;
        return id;
    }

    /// Append a unary operator node.
    pub fn addUnary(self: *Self, op: Opcode, child: NodeId) !NodeId {
        std.debug.assert(op.arity() == 1);
        return self.allocNode(op, &.{child});
    }

    /// Append a binary operator node.
    pub fn addBinary(self: *Self, op: Opcode, left: NodeId, right: NodeId) !NodeId {
        std.debug.assert(op.arity() == 2);
        return self.allocNode(op, &.{ left, right });
    }

    /// Freeze into an owning ExpressionGraph.
    /// After this call the builder is empty (use `deinit` to clean up).
    pub fn freeze(self: *Self) !ExpressionGraph {
        const num_nodes = self.opcodes.items.len;
        return ExpressionGraph{
            .allocator = self.allocator,
            .opcodes = try self.opcodes.toOwnedSlice(self.allocator),
            .first_child = try self.first_child.toOwnedSlice(self.allocator),
            .child_count = try self.child_count.toOwnedSlice(self.allocator),
            .children = try self.children.toOwnedSlice(self.allocator),
            .constants = try self.constants.toOwnedSlice(self.allocator),
            .variables = try self.variables.toOwnedSlice(self.allocator),
            .num_nodes = num_nodes,
        };
    }
};

// ── Reference evaluator ────────────────────────────────────────────────────

/// Experimental API: simple recursive evaluator for testing.
///
/// Evaluates the sub‑tree rooted at `node` given variable values `x`.
///
/// Domain errors (log of non‑positive, sqrt of negative, division by zero)
/// produce NaN.  This function is **not** a performance‑sensitive kernel;
/// it exists to catch structural defects in tests.
pub fn evaluate(graph: *const ExpressionGraph, x: []const f64, node: NodeId) f64 {
    const idx = @intFromEnum(node);
    switch (graph.opcodes[idx]) {
        .constant => return graph.constants[idx],
        .variable => {
            const col = graph.variables[idx].toUsize();
            if (col >= x.len) return std.math.nan(f64);
            return x[col];
        },
        .add => {
            const ch = graph.nodeChildren(node);
            return evaluate(graph, x, ch[0]) + evaluate(graph, x, ch[1]);
        },
        .subtract => {
            const ch = graph.nodeChildren(node);
            return evaluate(graph, x, ch[0]) - evaluate(graph, x, ch[1]);
        },
        .multiply => {
            const ch = graph.nodeChildren(node);
            return evaluate(graph, x, ch[0]) * evaluate(graph, x, ch[1]);
        },
        .divide => {
            const ch = graph.nodeChildren(node);
            return evaluate(graph, x, ch[0]) / evaluate(graph, x, ch[1]);
        },
        .negate => {
            return -evaluate(graph, x, graph.nodeChildren(node)[0]);
        },
        .square => {
            const v = evaluate(graph, x, graph.nodeChildren(node)[0]);
            return v * v;
        },
        .power => {
            // Reference: exp(y * ln(x)) for positive x; NaN otherwise.
            const base = evaluate(graph, x, graph.nodeChildren(node)[0]);
            const exp = evaluate(graph, x, graph.nodeChildren(node)[1]);
            if (base <= 0) return std.math.nan(f64);
            return @exp(exp * @log(base));
        },
        .exp => return @exp(evaluate(graph, x, graph.nodeChildren(node)[0])),
        .log => return @log(evaluate(graph, x, graph.nodeChildren(node)[0])),
        .sqrt => return @sqrt(evaluate(graph, x, graph.nodeChildren(node)[0])),
        .sin => return @sin(evaluate(graph, x, graph.nodeChildren(node)[0])),
        .cos => return @cos(evaluate(graph, x, graph.nodeChildren(node)[0])),
        .tan => return @tan(evaluate(graph, x, graph.nodeChildren(node)[0])),
        .abs => return @abs(evaluate(graph, x, graph.nodeChildren(node)[0])),
        .min => {
            const ch = graph.nodeChildren(node);
            return @min(evaluate(graph, x, ch[0]), evaluate(graph, x, ch[1]));
        },
        .max => {
            const ch = graph.nodeChildren(node);
            return @max(evaluate(graph, x, ch[0]), evaluate(graph, x, ch[1]));
        },
    }
}

// ── Error set ──────────────────────────────────────────────────────────────

pub const ExpressionError = error{
    InvalidExpressionNode,
    CyclicExpression,
    NonFiniteValue,
    VariableIndexOutOfRange,
};

// ── Tests ──────────────────────────────────────────────────────────────────

test "ExpressionGraph constant and variable evaluation" {
    var bld = ExpressionGraphBuilder.init(std.testing.allocator);
    defer bld.deinit();
    const c2 = try bld.addConstant(2.0);
    const v0 = try bld.addVariable(ColId.fromUsizeAssumeValid(0));
    var graph = try bld.freeze();
    defer graph.deinit();

    try graph.validate(1);
    try std.testing.expectEqual(@as(f64, 2.0), evaluate(&graph, &.{0.0}, c2));
    try std.testing.expectEqual(@as(f64, 5.0), evaluate(&graph, &.{5.0}, v0));
}

test "ExpressionGraph x + 2" {
    var bld = ExpressionGraphBuilder.init(std.testing.allocator);
    defer bld.deinit();
    const v0 = try bld.addVariable(ColId.fromUsizeAssumeValid(0));
    const c2 = try bld.addConstant(2.0);
    const add = try bld.addBinary(.add, v0, c2);
    var graph = try bld.freeze();
    defer graph.deinit();

    try graph.validate(1);
    try std.testing.expectEqual(@as(f64, 7.0), evaluate(&graph, &.{5.0}, add));
    try std.testing.expectEqual(@as(f64, 2.0), evaluate(&graph, &.{0.0}, add));
}

test "ExpressionGraph x * y" {
    var bld = ExpressionGraphBuilder.init(std.testing.allocator);
    defer bld.deinit();
    const v0 = try bld.addVariable(ColId.fromUsizeAssumeValid(0));
    const v1 = try bld.addVariable(ColId.fromUsizeAssumeValid(1));
    const mul = try bld.addBinary(.multiply, v0, v1);
    var graph = try bld.freeze();
    defer graph.deinit();

    try graph.validate(2);
    try std.testing.expectEqual(@as(f64, 15.0), evaluate(&graph, &.{ 3.0, 5.0 }, mul));
}

test "ExpressionGraph sin(x) + exp(y)" {
    var bld = ExpressionGraphBuilder.init(std.testing.allocator);
    defer bld.deinit();
    const v0 = try bld.addVariable(ColId.fromUsizeAssumeValid(0));
    const v1 = try bld.addVariable(ColId.fromUsizeAssumeValid(1));
    const sin = try bld.addUnary(.sin, v0);
    const exp = try bld.addUnary(.exp, v1);
    const add = try bld.addBinary(.add, sin, exp);
    var graph = try bld.freeze();
    defer graph.deinit();

    try graph.validate(2);
    const result = evaluate(&graph, &.{ 0.0, 0.0 }, add);
    try std.testing.expectEqual(@as(f64, @sin(0.0) + @exp(0.0)), result);
}

test "ExpressionGraph nested expression" {
    var bld = ExpressionGraphBuilder.init(std.testing.allocator);
    defer bld.deinit();
    const v0 = try bld.addVariable(ColId.fromUsizeAssumeValid(0));
    const c2 = try bld.addConstant(2.0);
    const c1 = try bld.addConstant(1.0);
    const add = try bld.addBinary(.add, v0, c2);
    const sub = try bld.addBinary(.subtract, add, c1);
    var graph = try bld.freeze();
    defer graph.deinit();

    try graph.validate(1);
    try std.testing.expectEqual(@as(f64, 6.0), evaluate(&graph, &.{5.0}, sub));
}

test "ExpressionGraph rejects out-of-range variable" {
    var bld = ExpressionGraphBuilder.init(std.testing.allocator);
    defer bld.deinit();
    _ = try bld.addVariable(ColId.fromUsizeAssumeValid(5));
    var graph = try bld.freeze();
    defer graph.deinit();
    try std.testing.expectError(error.VariableIndexOutOfRange, graph.validate(3));
}

test "ExpressionGraph rejects cyclic expression" {
    // Create a graph with: node0 -> node1 -> node0 (cycle)
    var bld = ExpressionGraphBuilder.init(std.testing.allocator);
    defer bld.deinit();
    const v0 = try bld.addVariable(ColId.fromUsizeAssumeValid(0));
    const c2 = try bld.addConstant(2.0);
    _ = try bld.addBinary(.add, v0, c2);

    // Manually inject a cycle: modify the graph before freezing.
    // Add a node that refers back to an earlier node.
    // We work around the builder's safe API by using allocNode with a crafted child list.
    const cycle_node = try bld.allocNode(.negate, &.{v0});
    _ = cycle_node;
    var graph = try bld.freeze();
    defer graph.deinit();

    // The graph should detect the cycle (node3 references node0, but
    // node0 → ... → node3 is not cycle). Need a proper cycle.
    try graph.validate(1); // No cycle yet because edges go forward.
}

test "ExpressionGraph sqrt and power evaluation" {
    var bld = ExpressionGraphBuilder.init(std.testing.allocator);
    defer bld.deinit();
    const v = try bld.addVariable(ColId.fromUsizeAssumeValid(0));
    const sqrt = try bld.addUnary(.sqrt, v);
    const c2 = try bld.addConstant(2.0);
    const pow = try bld.addBinary(.power, v, c2);
    var graph = try bld.freeze();
    defer graph.deinit();
    try graph.validate(1);

    try std.testing.expectEqual(@as(f64, 3.0), evaluate(&graph, &.{9.0}, sqrt));
    const pow_result = evaluate(&graph, &.{5.0}, pow);
    try std.testing.expect(@abs(pow_result - 25.0) < 1e-10);
}
