# Real sparse-matrix acceptance results (2026-07-14)

## Environment and method

- CPU: Intel Core i9-10900KF, 10 cores / 20 threads, 20 MiB shared L3.
- Toolchain: Zig 0.16.0 `ReleaseFast -Dcpu=native`; GCC 13.3 `-O3 -march=native -flto`.
- CPU affinity: core 2.
- Reference: HiGHS 1.14.0 checkout
  `de09bbad9fb7c5d39a1a464a7641bbb5531c6e9d`.
- Corpus: SuiteSparse Matrix Collection; URLs and SHA-256 values are pinned in
  `tools/matrix_gate/suitesparse-corpus.lock.tsv`.
- SpMV reports the mean of 5--50 warm iterations, selected to visit about 100
  million nonzeros. Format conversion reports the mean of seven warm,
  capacity-reusing iterations. Both implementations use identical values,
  input vectors, CPU affinity and repeat counts.

The zhighs acceptance runner additionally validates canonical CSC, CSC/CSR
semantic equality, explicit transpose, exact power-of-two scaling round-trip,
column slicing and permutation round-trip. Its peak RSS therefore describes the
complete acceptance workflow, not the memory footprint of CSC alone; it must
not be compared directly with the narrower HiGHS kernel runner RSS.

## Dataset results

| dataset | rows | columns | nnz | parse + canonical build | CSC SpMV | CSR SpMV | CSC -> CSR | transpose | acceptance peak RSS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| thermal1 | 82,654 | 82,654 | 574,458 | 317.867 ms | 0.628 ms | 0.587 ms | 1.546 ms | 1.601 ms | 74.5 MiB |
| cage12 | 130,228 | 130,228 | 2,032,536 | 1,028.162 ms | 1.882 ms | 1.767 ms | 8.814 ms | 8.924 ms | 130.7 MiB |
| webbase-1M | 1,000,005 | 1,000,005 | 3,105,536 | 762.466 ms | 4.986 ms | 4.520 ms | 13.568 ms | 14.687 ms | 204.7 MiB |

All three datasets passed every semantic check and the fail-closed large
dataset gate.

## Same-machine comparison with HiGHS

| dataset | kernel | zhighs | HiGHS | zhighs relative to HiGHS |
| --- | --- | ---: | ---: | ---: |
| thermal1 | CSC SpMV | 0.628 ms | 0.613 ms | 2.4% slower |
| cage12 | CSC SpMV | 1.882 ms | 1.783 ms | 5.6% slower |
| webbase-1M | CSC SpMV | 4.986 ms | 4.859 ms | 2.6% slower |
| thermal1 | CSR SpMV | 0.587 ms | 0.523 ms | 12.2% slower |
| cage12 | CSR SpMV | 1.767 ms | 1.781 ms | 0.8% faster |
| webbase-1M | CSR SpMV | 4.520 ms | 5.001 ms | 9.6% faster |
| thermal1 | CSC -> CSR | 1.546 ms | 1.851 ms | 16.5% faster |
| cage12 | CSC -> CSR | 8.814 ms | 8.723 ms | 1.0% slower |
| webbase-1M | CSC -> CSR | 13.568 ms | 14.365 ms | 5.5% faster |

Interpretation: current zhighs kernels are in the same single-thread CPU
performance class as HiGHS on these matrices. There is no evidence of a broad
performance lead: CSC SpMV consistently trails by 2--6%, while CSR behavior is
matrix-dependent. Reusable CSC-to-CSR conversion is competitive and usually
faster. The sample count is sufficient for a production acceptance baseline,
but not for claims covering every SuiteSparse sparsity distribution or other
hardware architectures.

## Commercial-readiness decision

The agreed matrix acceptance boundary now passes all four gates: HiGHS
differential, structural fuzz/OOM safety, real datasets, and multi-configuration
regression. The matrix storage and kernel layer is therefore acceptable as a
production candidate for its currently implemented scope.

This does **not** make the complete solver commercially ready. Commercial
deployment still needs sparse factorization/update kernels, end-to-end
presolve/simplex/MIP workloads, long-running soak and concurrency tests,
isolated memory accounting under configured limits, more CPU architectures,
ABI/versioning policy, observability, and real customer-model incident data.
The honest conclusion is: stop speculative matrix-layout work, integrate this
layer into the solver, and reopen matrix optimization only when end-to-end
profiles identify a measured bottleneck.

## 2026-07-14 perf/assembly-driven kernel update

The earlier conclusion that CSC consistently trailed HiGHS was reopened after
`perf record` showed a benchmark code-generation artifact. Zig had inlined the
kernel into a very large profiling `main`, spilling array bases to the stack,
while HiGHS called an out-of-line shared-library function. The production
change establishes stable leaf boundaries, uses explicit `@mulAdd` (matching
the FMA emitted for HiGHS), and dispatches once to the existing compact column
offsets. There is no per-nonzero format branch.

Five independent processes were run per implementation on CPU core 2. Values
below are process medians. Runs were grouped by implementation; an interleaved
real-dataset runner remains required before treating sub-5% differences as a
portable lead.

| dataset | kernel | zhighs median | HiGHS median | zhighs relative |
| --- | --- | ---: | ---: | ---: |
| thermal1 | CSC SpMV | 0.561 ms | 0.602 ms | 6.8% faster |
| cage12 | CSC SpMV | 1.774 ms | 1.842 ms | 3.7% faster |
| webbase-1M | CSC SpMV | 4.712 ms | 5.162 ms | 8.7% faster |
| thermal1 | CSR SpMV | 0.526 ms | 0.540 ms | 2.6% faster |
| cage12 | CSR SpMV | 1.648 ms | 1.845 ms | 10.7% faster |
| webbase-1M | CSR SpMV | 4.314 ms | 5.094 ms | 15.3% faster |
| thermal1 | CSC -> CSR | 1.550 ms | 1.910 ms | 18.8% faster |
| cage12 | CSC -> CSR | 8.594 ms | 8.944 ms | 3.9% faster |
| webbase-1M | CSC -> CSR | 13.675 ms | 14.862 ms | 8.0% faster |

Synthetic fixed-core `perf stat` corroborated the direction: CSC dropped from
about 363M to 282M cycles per profiling process and CSR from about 355M to
222M. The HiGHS references were about 318M and 283M cycles. Full matrix
acceptance subsequently passed all four gates.

## 2026-07-14 post-layout/SIMD/allocator rerun

This rerun supersedes the grouped five-process comparison above. Synthetic
kernels used 11 alternating-order processes on CPU 2; real datasets used seven
alternating-order processes. Zig was built with `ReleaseFast -Dcpu=native` and
the C++ harness with `-O3 -march=native -DNDEBUG -flto`, against the same HiGHS
commit `de09bbad9fb7c5d39a1a464a7641bbb5531c6e9d`. Every checksum and structural
hash matched. Raw synthetic results are in
`/tmp/zhighs-matrix-after-layout/results`, and real-process reports are in
`/tmp/zhighs-matrix-after-layout/real`.

### Interleaved real-data medians

| dataset | kernel | zhighs | HiGHS | zhighs relative |
| --- | --- | ---: | ---: | ---: |
| thermal1 | CSC SpMV | 0.559 ms | 0.603 ms | 7.3% faster |
| cage12 | CSC SpMV | 1.841 ms | 1.860 ms | 1.0% faster |
| webbase-1M | CSC SpMV | 4.604 ms | 5.332 ms | 13.7% faster |
| thermal1 | CSR SpMV | 0.529 ms | 0.533 ms | 0.8% faster |
| cage12 | CSR SpMV | 1.696 ms | 1.803 ms | 5.9% faster |
| webbase-1M | CSR SpMV | 4.273 ms | 5.128 ms | 16.7% faster |
| thermal1 | CSC -> CSR reusable | 1.571 ms | 1.879 ms | 16.4% faster |
| cage12 | CSC -> CSR reusable | 9.050 ms | 8.994 ms | 0.6% slower |
| webbase-1M | CSC -> CSR reusable | 14.002 ms | 14.612 ms | 4.2% faster |

The sub-1% thermal CSR and cage conversion differences are parity, not a
meaningful lead. The larger webbase advantages and thermal conversion advantage
were repeatable across the seven samples. The current conclusion is therefore
that zhighs SpMV is at least in the same class as HiGHS and usually faster on
this corpus, but it is not uniformly faster for every sparsity distribution.

### Stable synthetic ownership and construction results

| kernel (50k dimension, 149,998 nnz) | zhighs | C++/HiGHS harness | result |
| --- | ---: | ---: | ---: |
| CSR dense SpMV | 89.990 us | 113.018 us | zhighs 20.4% faster |
| apply full row/column scale | 107.836 us | 118.864 us | zhighs 9.3% faster |
| CSC -> CSR reusable | 272.527 us | 334.395 us | zhighs 18.5% faster |
| CSC -> CSR owning | 1,365.465 us | 1,396.231 us | zhighs 2.2% faster |
| transpose reusable/into | 291.938 us | 288.059 us | parity; zhighs 1.3% slower |
| transpose owning | 920.299 us | 337.142 us | zhighs 2.73x slower |
| builder canonical owning | 745.046 us | 256.373 us | zhighs 2.91x slower |
| builder reusable | 277.446 us | 310.629 us | zhighs 10.7% faster |
| builder general/sort | 3,009.064 us | 3,433.172 us | zhighs 12.4% faster |
| sparse accumulate | 81.419 us | 243.653 us | zhighs 2.99x faster |

The owning-transpose result improved about 7.7% versus the immediately prior
Zig implementation (compact output reused as cursor), but remains the largest
HiGHS gap. The near-equal reusable transpose shows that counting/prefix/scatter
is not the cause. The owning Zig result constructs and retains both `usize` and
`HUInt` starts, whereas the C++ result stores one `HighsInt` stream. The same
representation and output-construction cost is visible in canonical builder;
the reusable builder is already faster. These are the next justified perf and
assembly targets. `csc_ax_dense`, `csc_atx_dense`, `alpha_ax_plus_y`, and
`product_quad` had 13-45% median absolute deviation in at least one side during
this synthetic run, so their single-run rankings are not used for a verdict;
the real-data SpMV table is the authoritative product comparison.

### Fair memory interpretation

The acceptance runner peak RSS is not comparable with the narrower C++ runner:
it deliberately keeps CSC, CSR workspace, transpose workspace, scaling copies,
slice and permutation validation objects alive. For authoritative w32 CSC
streams, the comparable formulas are currently `12*(cols+1) + 12*nnz` bytes
for zhighs (wide plus compact starts) and `4*(cols+1) + 12*nnz` for HiGHS.

| dataset | zhighs authoritative CSC | HiGHS authoritative CSC | zhighs overhead |
| --- | ---: | ---: | ---: |
| thermal1 | 7.52 MiB | 6.89 MiB | 9.2% |
| cage12 | 24.75 MiB | 23.76 MiB | 4.2% |
| webbase-1M | 46.98 MiB | 39.35 MiB | 19.4% |

Thus the earlier whole-process impression that zhighs was universally crushed
on memory was not a fair CSC-only comparison. A real and material deficit does
remain for matrices with many columns: the duplicate wide/compact offset
streams cost 7.63 MiB on webbase-1M and directly affect owning construction.
