# Sparse simplex basis design

This document fixes the performance and ownership contract for the sparse
basis factorization work. It prevents microbenchmark improvements from
silently changing simplex semantics or being compared with a different amount
of work in HiGHS.

## Current boundary

`SparseBasisBuffers.assemble` converts the immutable model CSC plus the current
`basic_index` into compact basis CSC. Structural columns are copied with row
scaling; logical and artificial columns are synthesized as signed unit columns.
The buffers retain capacity across reinversions and own no model or basis
indices.

The representation is data-oriented SoA:

```text
starts: HUInt[dimension + 1]
rows:   RowId[nnz]
values: f64[nnz]
```

The default w32 build therefore uses four-byte starts and row IDs. All streams
are independently 64-byte aligned. Symbolic factorization can read
`starts/rows` without pulling values into cache; numerical factorization and
FTRAN/BTRAN can stream `rows/values` without AoS padding.

Target decisions use Zig's `@import("builtin")` through `target_policy.zig`.
SIMD width and sparse prefetch distance are compile-time constants, so there is
no CPUID, virtual dispatch, or target branch in the assembly loop.

## HiGHS comparison contract

The local reference is `/home/godv/codefiles/cppfiles/HiGHS`, primarily
`highs/util/HFactor.cpp` and `HFactor.h`. HiGHS `HFactor::build` performs basis
assembly, singleton elimination, Markowitz pivoting, numerical factorization,
and factor construction. Comparing it directly with assembly alone would be
invalid.

Fair gates are separated as follows:

1. **Assembly:** zhighs retaining assembly versus a small C++ oracle that only
   materializes the same selected CSC columns and signed unit columns.
2. **INVERT:** complete zhighs sparse factorization versus `HFactor::build`,
   starting from identical CSC and `basic_index` arrays.
3. **FTRAN/BTRAN:** identical factorization, RHS density, warmups, repetitions,
   checksums, and solution residuals.
4. **End-to-end simplex:** identical LP file, options, tolerances, CPU affinity,
   allocator lifetime, and status/KKT checks.

Every report records median, MAD, min/max, requested bytes, peak RSS, factor
nonzeros, fill ratio, pivots, reinversions, residuals, and checksums. A speed
claim requires interleaved process order and low-noise samples; a single best
run is never used.

## Dataset gates

- Synthetic diagonal, block triangular, arrow, network, and controlled-fill
  bases isolate symbolic and numerical behavior.
- SuiteSparse matrices already available under
  `/home/godv/codefiles/datasets/suitesparse/matrices` exercise cache and fill:
  `thermal1`, `cage12`, and `webbase-1M`.
- Netlib LP models provide end-to-end status, objective, KKT and basis checks.
- Larger Mittelmann models are added only after the Netlib gate is clean.

Matrix Market inputs are not automatically valid simplex bases. Dataset
benchmarks must construct a nonsingular basis deterministically and verify it,
rather than timing arbitrary square slices that may be rank deficient.

## Next implementation stages

1. Reusable symbolic workspace: row counts, column counts, permutations and
   count-linked lists in SoA storage.
2. Singleton row/column elimination before the Markowitz kernel.
3. Threshold Markowitz pivoting with deterministic ties and rank-deficiency
   reporting.
4. Packed sparse L/U plus row-wise companion views required by BTRAN.
5. Dense/sparse FTRAN and BTRAN selected by measured RHS density.
6. Dense-LU oracle comparison, iterative refinement and residual gates.
7. Forrest--Tomlin updates after clean reinversion and solve benchmarks.

The existing dense LU remains the small-basis fallback and correctness oracle
until the sparse backend passes every numerical and end-to-end gate.
