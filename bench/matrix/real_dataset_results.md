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

### Corrected fair-allocator synthetic results

The first synthetic rerun used Zig `smp_allocator` while C++ `std::vector`
used libc malloc. That was not allocator-fair for repeated multi-megabyte
owning operations: smp returned large mappings to the OS on every iteration,
while malloc retained and reused them. The corrected runner uses Zig
`c_allocator` only for timed owning allocations; fixtures and allocation-free
kernels are unchanged.

| kernel (50k dimension, 149,998 nnz) | zhighs | C++/HiGHS harness | result |
| --- | ---: | ---: | ---: |
| CSC dense SpMV | 114.519 us | 126.147 us | zhighs 9.2% faster |
| CSR dense SpMV | 89.136 us | 112.666 us | zhighs 20.9% faster |
| apply full row/column scale | 108.555 us | 119.462 us | zhighs 9.1% faster |
| CSC -> CSR reusable | 271.822 us | 339.225 us | zhighs 19.9% faster |
| CSC -> CSR owning | 1,435.375 us | 1,363.925 us | zhighs 5.2% slower |
| transpose reusable/into | 297.416 us | 296.905 us | parity; zhighs 0.2% slower |
| transpose owning | 360.382 us | 339.685 us | zhighs 6.1% slower |
| builder sorted owning | 646.968 us | 1,258.392 us | zhighs 48.6% faster |
| builder prepopulated owning | 414.441 us | 691.635 us | zhighs 40.1% faster |
| builder canonical owning | 231.101 us | 262.858 us | zhighs 12.1% faster |
| builder reusable | 278.216 us | 309.656 us | zhighs 10.2% faster |
| builder general/sort | 1,591.723 us | 3,437.349 us | zhighs 53.7% faster |
| sparse accumulate | 80.115 us | 241.703 us | zhighs 3.02x faster |

The same Zig transpose code measured 905 us with smp, 357 us with c allocator,
and 340 us in C++. `perf stat` retired about 403M instructions for both Zig
allocator variants, but smp incurred about 61,867 page faults versus 3,853 for
c and 4,038 for C++. With allocator policy equalized, Zig and C++ transpose
also had essentially equal instructions and user cycles; the remaining 6.1%
wall-time difference does not justify changing the already-parity scatter
kernel. Canonical builder similarly changed from a false 2.91x deficit to a
12.1% lead; perf measured about 1.19B Zig versus 1.86B C++ instructions.

All corrected medians had MAD below 4.4%. Explicit high-side anomalies remain
in the raw report: C++ owning transpose reached 436 us versus its 340 us median,
and Zig sorted builder reached 808 us versus its 647 us median. They were not
discarded; median/MAD conclusions remain unchanged. `alpha_ax_plus_y` is the
remaining stable compute deficit at 23.4%; CSR transpose-product is 3.8% slower
and `product_quad` 2.4% slower. These are separate compute-kernel tasks and were
not modified during owning-allocation work.

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

An RSS-only probe then ran each dataset in a fresh process, retained only the
canonical CSC, and used libc allocation on both sides. Both parsers retained
the complete Matrix Market text during construction; the C++ RSS probe used a
non-copying span stream so the input text existed exactly once. Seven runs per
dataset were invariant or varied by at most 164 KiB; there were no material
outliers.

| dataset | current RSS zhighs / HiGHS | peak build RSS zhighs / HiGHS |
| --- | ---: | ---: |
| thermal1 | 9.46 / 12.05 MiB | 32.19 / 39.22 MiB |
| cage12 | 26.70 / 28.93 MiB | 131.88 / 149.53 MiB |
| webbase-1M | 48.93 / 44.52 MiB | 187.34 / 206.56 MiB |

Current RSS is page-residency and allocator dependent, so requested bytes are
the authoritative final-storage comparison. Peak construction RSS is useful:
zhighs was 9-18% lower because its builder merges canonical entries in place,
while the C++ reference held both input and canonical triplet vectors. This
does not erase the final dual-offset overhead on webbase; it separates that
19.4% storage cost from transient construction memory.
