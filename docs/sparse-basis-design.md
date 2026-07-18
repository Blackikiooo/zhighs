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

Completed initial symbolic layer:

- Retaining CSC-to-row-entry companion construction without duplicating values.
- Incremental active row/column counts, intrusive column-count buckets and a
  monotonic singleton-row queue.
- Fill-free row/column singleton elimination.
- Threshold Markowitz kernel-pivot selection with deterministic ties.

The symbolic API deliberately stops after the first non-singleton pivot.
Selecting later pivots before numerical elimination would ignore fill-in and
updated values, producing an invalid order. The next numerical layer must
apply the selected pivot, update the Schur complement and bucket membership,
then request the next choice.

Remaining stages:

Completed numerical MVP:

- Reusable SoA entry pool with intrusive row/column lists and recycled slots.
- O(1) row/column count-bucket relocation on fill insertion and removal.
- Numerical Schur updates, configurable zero dropping, and threshold Markowitz
  re-evaluation after every pivot.
- Packed L columns and U rows with explicit `P B Q = L U` permutations.
- Allocation-free packed FTRAN and BTRAN after factorization.

Remaining stages:

1. Rank-deficiency repair and stronger pivot-growth/condition monitoring.
2. Hyper-sparse FTRAN/BTRAN selected by measured RHS density.
3. Iterative refinement and integration behind the simplex factorization API.
4. SuiteSparse/Netlib gates and complete HiGHS HFactor parity reports.
5. Forrest--Tomlin updates after clean reinversion and solve benchmarks.

The existing dense LU remains the small-basis fallback and correctness oracle
until the sparse backend passes every numerical and end-to-end gate.

## Initial numerical MVP measurement

The cyclic tridiagonal benchmark uses identical CSC values, basic indices,
warmup and retaining factor objects. zhighs uses ReleaseFast `-Dcpu=native`;
HiGHS `de09bbad9f` uses Release with Clang `-O3 -march=native -flto`. Setup and
initial allocations are outside both measured INVERT regions.

The first MVP is correct but not yet faster than HiGHS:

| dimension | zhighs INVERT | HiGHS HFactor | ratio |
|---:|---:|---:|---:|
| 128 | 27.4 us | 18.1 us | 1.51x |
| 256 | 55.9 us | 35.3 us | 1.58x |
| 512 | 125.0 us | 69.6 us | 1.80x |
| 1024 | 241.5 us | 139.6 us | 1.73x |

No performance-leading claim is made. On the 512 case, packed zhighs FTRAN
plus BTRAN is about 10 us, but the current comparison runner measures only the
fair common INVERT boundary. `perf` on the 1024 case identifies general
Markowitz scanning and intrusive entry removal as the next tuning targets.

### Non-allocator tuning follow-up

The next pass kept the factor allocator fixed while making four measured
changes: retiring pivot dimensions before unlinking entries, direct CSC pool
loading, removing the redundant entry-alive stream, and stopping Markowitz
search when the theoretical count-merit lower bound is reached. The latter
retains deterministic traversal but does not scan an entire bucket merely to
choose the smallest ID among mathematically equivalent pivots.

With `c_allocator` on the same cyclic tridiagonal workload (501 medians):

| dimension | zhighs INVERT | HiGHS HFactor | zhighs speedup |
|---:|---:|---:|---:|
| 128 | 10.2 us | 17.8 us | 1.75x |
| 256 | 20.0 us | 34.5 us | 1.73x |
| 512 | 39.1 us | 69.3 us | 1.77x |
| 1024 | 79.3 us | 138.9 us | 1.75x |

This is a narrow structured-basis result, not a general sparse-LU superiority
claim. More fill-heavy synthetic bases and SuiteSparse/Netlib-derived valid
bases remain required.

An interleaved seven-process allocator audit at dimension 1024 found warm GPA
and `c_allocator` medians both around 185 us before the Markowitz early-stop
change, with less than one percent separation. The factorization benchmark now
links libc and uses `c_allocator` for the factor workspace by default in the
HiGHS comparison script. Fixtures and samples remain outside the timed region.
The warm allocation-count test independently proves that reinversion and
FTRAN/BTRAN issue no allocator calls once retained capacity is sufficient.
