//! Revised-simplex orchestration boundary; policies remain replaceable.
const std = @import("std");
const basis_module = @import("basis.zig");
const factorization_module = @import("factorization.zig");
const pricing_module = @import("pricing.zig");
const ratio_module = @import("ratio_test.zig");
const numerical_module = @import("numerical.zig");
const problem_module = @import("problem.zig");

pub const Algorithm = enum { primal_revised, dual_revised };
pub const SolveStatus = enum { optimal, infeasible, unbounded, iteration_limit, numerical_failure, not_implemented };
pub const SolveControl = struct { max_iterations: usize = 1_000_000, time_limit_ns: u64 = 0 };

pub const SimplexEngine = struct {
    allocator: std.mem.Allocator,
    algorithm: Algorithm = .dual_revised,
    basis: ?basis_module.BasisState = null,
    factorization: factorization_module.Factorization,
    pricing: pricing_module.Pricing = .{},
    ratio_test: ratio_module.RatioTest = .{},
    numerical: numerical_module.NumericalState = .{},

    pub fn init(a: std.mem.Allocator) SimplexEngine {
        return .{ .allocator = a, .factorization = factorization_module.Factorization.init(a) };
    }
    pub fn deinit(self: *SimplexEngine) void {
        if (self.basis) |*b| b.deinit();
        self.factorization.deinit();
    }
    pub fn solve(_: *SimplexEngine, _: usize, _: usize, _: SolveControl) SolveStatus {
        return .not_implemented;
    }

    /// Solve entry point consuming a borrowed LP `ProblemView`.
    ///
    /// The view must outlive this call; the engine never takes ownership of
    /// model arrays. Basis storage is owned by the engine instance.
    pub fn solveProblem(self: *SimplexEngine, problem: problem_module.ProblemView, control: SolveControl) SolveStatus {
        if (self.basis) |*old| old.deinit();
        self.basis = basis_module.BasisState.init(self.allocator, problem.num_rows, problem.num_cols) catch return .numerical_failure;
        self.basis.?.initializeSlackBasis();
        self.numerical.markRefactorized();
        return self.solve(problem.num_rows, problem.num_cols, control);
    }
};

test {
    std.testing.refAllDecls(@This());
}
