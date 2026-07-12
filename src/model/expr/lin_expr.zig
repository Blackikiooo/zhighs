//! Linear expression `c + Σ aᵢ xᵢ`.

const std = @import("std");
const Model = @import("../model.zig").Model;
const Var = @import("../var/index.zig").Var;
const ModelError = @import("../types.zig").ModelError;

/// A linear expression `c + Σ aᵢ xᵢ` where `c` is a constant term and each
/// pair `(xᵢ, aᵢ)` is a variable-coefficient term.
///
/// ```zig
/// var expr = LinExpr.init(allocator);
/// defer expr.deinit();
/// try expr.addTerm(var_0, 2.0);
/// try expr.addTerm(var_1, -1.0);
/// try expr.addConstant(3.0);  // 2*x0 - x1 + 3
/// ```
pub const LinExpr = struct {
    allocator: std.mem.Allocator,
    terms: std.ArrayListUnmanaged(struct { var_idx: usize, coeff: f64 }) = .empty,
    constant: f64 = 0.0,

    const Self = @This();

    /// Create an empty linear expression.
    pub fn init(allocator: std.mem.Allocator) Self {
        return .{ .allocator = allocator };
    }

    /// Release memory.
    pub fn deinit(self: *Self) void {
        self.terms.deinit(self.allocator);
    }

    /// Add a term `coeff * var` to the expression.
    pub fn addTerm(self: *Self, variable: Var, coeff: f64) ModelError!void {
        self.terms.append(self.allocator, .{ .var_idx = variable.index, .coeff = coeff }) catch return error.OutOfMemory;
    }

    /// Add a term by raw variable index.
    pub fn addTermByIndex(self: *Self, var_idx: usize, coeff: f64) ModelError!void {
        self.terms.append(self.allocator, .{ .var_idx = var_idx, .coeff = coeff }) catch return error.OutOfMemory;
    }

    /// Add a constant term to the expression.
    pub fn addConstant(self: *Self, value: f64) void {
        self.constant += value;
    }

    /// Add another linear expression to this one.
    pub fn addExpr(self: *Self, other: LinExpr) ModelError!void {
        for (other.terms.items) |t| {
            try self.addTermByIndex(t.var_idx, t.coeff);
        }
        self.constant += other.constant;
    }

    /// Evaluate the expression using the model's solution vector.
    pub fn getValue(self: Self, model: *Model) ModelError!f64 {
        var val = self.constant;
        for (self.terms.items) |t| {
            const x = try model.getDblAttrElement(.x, t.var_idx);
            val += t.coeff * x;
        }
        return val;
    }

    /// Return the number of variable terms.
    pub fn numTerms(self: Self) usize {
        return self.terms.items.len;
    }

    /// Return the constant term.
    pub fn getConstant(self: Self) f64 {
        return self.constant;
    }
};
