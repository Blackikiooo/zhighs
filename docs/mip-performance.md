# MIP performance baseline

## Benchmark contract

`bench/mip/solver_bench.zig` parses the model once and then measures two paths
with the same allocator, process, model arrays, repetition count, and clock:

- one independent solve of the root LP relaxation;
- one complete MIP solve, including node LPs and branch-and-bound overhead.

Every reported row also includes status, objective, global bound, nodes, LP
iterations/work, cold solves, and reoptimizations. Timing samples are invalid if
these correctness and structural counters differ.

## 2026-09-10 smoke baseline

Environment: Windows, Zig 0.16.0, `ReleaseFast`, 1,000 repetitions per process,
`test/model/data/tiny_mip.mps`.

| sample | LP mean | MIP mean | nodes | objective | bound |
|---:|---:|---:|---:|---:|---:|
| 1 | 204,548 ns | 692,135 ns | 5 | 3 | 3 |
| 2 | 181,341 ns | 775,024 ns | 5 | 3 | 3 |
| 3 | 225,344 ns | 689,830 ns | 5 | 3 | 3 |

All samples used 3 LP iterations, 7 LP work units, 5 cold node solves, and no
reoptimization. The complete MIP solve costs roughly 3.1--4.3 independent root
LP solves despite processing 5 nodes. On this smoke model, LP/presolve dominates
and the node/domain layer is not the primary bottleneck.

An earlier 25-run result of 36,688 ns for the complete MIP path is below the
later directly measured root-LP time and was not reproducible. It is retained
as an invalid historical observation, not a regression baseline.

## Interpretation limits

This two-variable model validates instrumentation and result parity only. It
does not exercise a large open queue, deep domain paths, cut separation, or
large sparse LP relaxations. Performance claims require a pinned public MIP
corpus and comparison with a matching HiGHS build. Until that corpus exists,
this document must not claim parity or leadership over HiGHS.

The next profiling target is node LP cold-start behavior: this model reports
five cold solves and zero reoptimizations. Improving stable presolve mappings
and bound-only LP reoptimization is likely to matter more than further tuning
the already allocation-free default node-selection path.
