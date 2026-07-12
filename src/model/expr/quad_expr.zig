//! Quadratic expression `LinExpr + Σ Q[j,k] · x[j] · x[k]`.

const std = @import("std");
const Var = @import("../var/index.zig").Var;
const ModelError = @import("../types.zig").ModelError;
const LinExpr = @import("lin_expr.zig").LinExpr;

/// A quadratic expression:
///   `LinExpr + Σ Q[j,k] · x[j] · x[k]`
///
/// where the quadratic terms are stored as lower-triangle triples
/// `(row, col, coeff)`.
pub const QuadExpr = struct {
    allocator: std.mem.Allocator,
    linear: LinExpr,
    q_terms: std.ArrayListUnmanaged(struct { row: usize, col: usize, coeff: f64 }) = .empty,

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
        self.q_terms.append(self.allocator, .{ .row = row.index, .col = col.index, .coeff = coeff }) catch return error.OutOfMemory;
    }

    /// Add a constant to the linear part.
    pub fn addConstant(self: *Self, value: f64) void {
        self.linear.addConstant(value);
    }

    /// Return the number of quadratic terms.
    pub fn numQTerms(self: Self) usize {
        return self.q_terms.items.len;
    }
};
