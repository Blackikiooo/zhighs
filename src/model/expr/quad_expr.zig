//! Quadratic expression `LinExpr + Σ Q[j,k] · x[j] · x[k]`.

const std = @import("std");
const Var = @import("../var/index.zig").Var;
const ModelError = @import("../types.zig").ModelError;
const LinExpr = @import("lin_expr.zig").LinExpr;
const VarId = @import("../entity_handle.zig").VarId;
const Model = @import("../model.zig").Model;

/// A quadratic expression:
///   `LinExpr + Σ Q[j,k] · x[j] · x[k]`
///
/// where the quadratic terms are stored as lower-triangle triples
/// `(row, col, coeff)`.
pub const QuadExpr = struct {
    pub const ResolvedQTerm = struct { row: usize, col: usize, coeff: f64 };
    allocator: std.mem.Allocator,
    linear: LinExpr,
    q_terms: std.ArrayListUnmanaged(struct { row: usize, col: usize, row_id: ?VarId = null, col_id: ?VarId = null, coeff: f64 }) = .empty,

    const Self = @This();

    /// Create an empty quadratic expression.
    pub fn init(allocator: std.mem.Allocator) Self {
        return .{
            .allocator = allocator,
            .linear = LinExpr.init(allocator),
        };
    }

    /// Release memory.
    pub fn deinit(self: *Self) void {
        self.linear.deinit();
        self.q_terms.deinit(self.allocator);
    }

    /// Add a term `coeff * var` to the linear part.
    pub fn addTerm(self: *Self, variable: Var, coeff: f64) ModelError!void {
        try self.linear.addTerm(variable, coeff);
    }

    /// Add a term `coeff * x[row] * x[col]` to the quadratic part.
    pub fn addQTerm(self: *Self, row: Var, col: Var, coeff: f64) ModelError!void {
        self.q_terms.append(self.allocator, .{
            .row = row.index,
            .col = col.index,
            .row_id = row.id,
            .col_id = col.id,
            .coeff = coeff,
        }) catch return error.OutOfMemory;
    }

    /// Add a constant to the linear part.
    pub fn addConstant(self: *Self, value: f64) void {
        self.linear.addConstant(value);
    }

    /// Return the number of quadratic terms.
    pub fn numQTerms(self: Self) usize {
        return self.q_terms.items.len;
    }

    /// Resolve stable variable handles for the current dense model layout.
    pub fn resolveQTerms(self: Self, allocator: std.mem.Allocator, model: Model) ModelError![]ResolvedQTerm {
        const resolved = allocator.alloc(ResolvedQTerm, self.q_terms.items.len) catch return error.OutOfMemory;
        errdefer allocator.free(resolved);
        for (self.q_terms.items, 0..) |term, i| {
            resolved[i] = .{
                .row = if (term.row_id) |id| try model.resolveVarId(id) else term.row,
                .col = if (term.col_id) |id| try model.resolveVarId(id) else term.col,
                .coeff = term.coeff,
            };
        }
        return resolved;
    }
};
