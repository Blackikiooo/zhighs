//! Assembly-optimized kernels for sparse matrix hot paths.
//!
//! These routines address specific LLVM code-generation limitations
//! identified during benchmarking. Each is written in inline assembly
//! with a structured Zig fallback for other architectures.
//!
//! ## Architecture support
//!
//! | Arch   | Status | Notes                        |
//! |--------|--------|------------------------------|
//! | x86_64 | Full   | SSE2 baseline (all x86-64)   |
//! | other  | Zig    | No hand-tuned asm (fallback) |
//!
//! ## Module structure
//!
//! - [`clear`](clear.zig): memory zeroing (clearF64)
//! - [`csr`](csr.zig): CSR A^T x scatter-add (csrTransposeMultiply)
//!
//! ## Why hand-tuned assembly
//!
//! LLVM spills base pointers (yp, vs, ci) to the stack in scatter-add loops
//! because its alias analysis cannot prove that stores to y[col] don't alias
//! the pointer variables themselves. Hand-tuned assembly keeps all pointers
//! in registers, eliminating spill traffic entirely.

const builtin = @import("builtin");

pub const clear = @import("clear.zig");
pub const csr = @import("csr.zig");

test "asm module smoke test" {
    _ = clear;
    _ = csr;
}
