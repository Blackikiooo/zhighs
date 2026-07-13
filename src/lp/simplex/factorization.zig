//! Basis factorization policy boundary.
//!
//! The first implementation may use dense LU for correctness tests; the
//! interface deliberately permits sparse LU, eta/Forrest--Tomlin updates, and
//! periodic reinversion without changing the simplex engine.

const std = @import("std");
const matrix = @import("matrix");

pub const FactorizationError = error{ DimensionMismatch, NotImplemented, Singular, NumericalFailure, OutOfMemory };
pub const PivotUpdate = struct { leaving_row: u32, entering_col: u32, pivot: f64 };

pub const Factorization = struct {
    allocator: std.mem.Allocator,
    update_count: usize = 0,
    dense_lu: matrix.DenseLU,

    pub fn init(allocator: std.mem.Allocator) Factorization {
        return .{ .allocator = allocator, .dense_lu = matrix.DenseLU.init(allocator) };
    }
    pub fn deinit(self: *Factorization) void {
        self.dense_lu.deinit();
    }
    pub fn factorize(self: *Factorization, n: usize, matrix_data: []const f64) FactorizationError!void {
        self.dense_lu.factorize(n, matrix_data) catch |err| return switch (err) {
            error.DimensionMismatch => error.DimensionMismatch,
            error.Singular => error.Singular,
            error.NumericalFailure => error.NumericalFailure,
            error.OutOfMemory => error.OutOfMemory,
        };
    }
    pub fn solve(self: *Factorization, rhs: []f64) FactorizationError!void {
        self.dense_lu.solve(rhs) catch |err| return switch (err) {
            error.DimensionMismatch => error.DimensionMismatch,
            error.Singular => error.Singular,
            error.NumericalFailure => error.NumericalFailure,
            error.OutOfMemory => error.OutOfMemory,
        };
    }
    pub fn solveTranspose(self: *Factorization, rhs: []f64) FactorizationError!void {
        self.dense_lu.solveTranspose(rhs) catch |err| return switch (err) {
            error.DimensionMismatch => error.DimensionMismatch,
            error.Singular => error.Singular,
            error.NumericalFailure => error.NumericalFailure,
            error.OutOfMemory => error.OutOfMemory,
        };
    }
    pub fn update(self: *Factorization, _: PivotUpdate) FactorizationError!void {
        self.update_count += 1;
        return error.NotImplemented;
    }
    pub fn needsRefactor(self: Factorization, limit: usize) bool {
        return self.update_count >= limit;
    }
};

test {
    std.testing.refAllDecls(@This());
}
