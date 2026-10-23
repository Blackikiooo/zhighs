//! Numerical linear algebra used by optimization algorithms.
//!
//! Covers dense and sparse LU, basis factorization, FTRAN/BTRAN, update
//! strategies, singularity recovery, and reusable linear-system contracts.
//!
//! Layer boundary:
//! - `nla` owns numerical kernels, their workspaces, and factorization state.
//!   It consumes canonical sparse storage from `matrix` (for example
//!   `CscView` / `SparseBasisView`) but never owns model data.
//! - `matrix` must not import `nla`; the dependency direction is
//!   `nla -> matrix -> foundation`. Keep structural storage and numeric
//!   elimination separate so matrix benchmarks stay free of solver policy.
//!
//! API stability policy (same as `matrix`):
//! - A declaration documented as `Stable API` is intended for production use.
//! - Every other public declaration is currently an `Experimental API`. Its
//!   signature, ownership model, validation contract, or implementation may be
//!   corrected as LP/IPM/MIP layers expose new requirements.

const std = @import("std");
const dense_lu = @import("dense_lu.zig");
const sparse_basis = @import("sparse_basis.zig");
const sparse_symbolic = @import("sparse_symbolic.zig");
const sparse_kernel = @import("sparse_kernel.zig");
const sparse_lu = @import("sparse_lu.zig");
const sparse_ft = @import("sparse_ft.zig");
const sparse_ft_packed = @import("sparse_ft_packed.zig");
const sparse_ft_adaptive = @import("sparse_ft_adaptive.zig");

/// Dense LU with partial pivoting; correctness baseline for simplex.
pub const DenseLU = dense_lu.DenseLU;
pub const DenseLuError = dense_lu.DenseLuError;

/// Sparse basis assembly as compact CSC in caller-retained SoA buffers.
pub const SparseBasisError = sparse_basis.SparseBasisError;
pub const SparseBasisView = sparse_basis.SparseBasisView;
pub const SparseBasisBuffers = sparse_basis.SparseBasisBuffers;

/// Initial symbolic pivot analysis (threshold Markowitz).
pub const SymbolicError = sparse_symbolic.SymbolicError;
pub const PivotKind = sparse_symbolic.PivotKind;
pub const SymbolicPlanView = sparse_symbolic.SymbolicPlanView;
pub const SymbolicWorkspace = sparse_symbolic.SymbolicWorkspace;

/// Mutable data-oriented elimination kernel.
pub const KernelError = sparse_kernel.KernelError;
pub const SparsePivotChoice = sparse_kernel.PivotChoice;
pub const SparsePivotView = sparse_kernel.PivotView;
pub const MutableSparseKernel = sparse_kernel.MutableSparseKernel;

/// Packed sparse LU (`P * B * Q = L * U`) plus FTRAN/BTRAN and updates.
pub const SparseLuError = sparse_lu.SparseLuError;
pub const SparseLU = sparse_lu.SparseLU;

/// Forrest--Tomlin upper factor and its packed/adaptive backends.
pub const FtError = sparse_ft_adaptive.FtError;
pub const SparseForrestTomlin = sparse_ft.SparseForrestTomlin;
pub const PackedSparseForrestTomlin = sparse_ft_packed.SparseForrestTomlin;
pub const AdaptiveSparseForrestTomlin = sparse_ft_adaptive.SparseForrestTomlin;
pub const MutableUpperView = sparse_ft.MutableUpperView;

test {
    std.testing.refAllDecls(@This());
}
