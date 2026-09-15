# MIP performance baseline

## Benchmark contract

`bench/mip/solver_bench.zig` parses the model once and then measures two paths
with the same allocator, process, model arrays, repetition count, and clock:

- one independent solve of the root LP relaxation;
- one complete MIP solve, including node LPs and branch-and-bound overhead.

Every reported row also includes status, objective, global bound, nodes, LP
iterations/work, cold solves, reoptimizations, and pseudocost observations.
Strong-branch probes, iterations, and work are reported separately while also
counting toward total LP budgets. Cut rounds, additions, and expirations expose
the structural work performed by separation; clique edge and tightening counts
do the same for binary conflict propagation. Timing samples are invalid if
these correctness and structural counters differ.

Reduced-cost tightenings are also reported. They count persistent integer
domain changes justified by the incumbent cutoff and an optimal LP dual
solution; they are zero when no incumbent is available or the feature is
disabled.
`root_reduced_cost_tightenings` separately counts bounds promoted to the shared
root domain from the retained root LP certificate, so all queued nodes inherit
them without duplicating sparse node changes.

Conflict instrumentation reports learned no-goods, unit-propagated bound
tightenings, examined conflict ranges, aged-out entries, and physically
compacted ranges. Learning uses
retained scratch capacity and is skipped unless an infeasible LP path consists
entirely of safely negatable ordinary integer bounds.
The `no-conflicts` profile changes only `conflict_learning` for paired
measurement of learning, propagation, aging, and compaction overhead.

The `no-cliques` benchmark profile changes only `clique_propagation`, allowing
the extraction and propagation overhead to be isolated without changing cut,
branching, or heuristic controls.

On `tiny_mip.mps`, two 10,000-solve ReleaseFast samples measured
18,281--18,589 ns with the one-edge literal graph and 18,324--18,470 ns with
`no-cliques`. Both configurations used one node, seven LP work units, two cold
solves, and the same objective and bound. The intervals overlap, so this smoke
test shows no distinguishable regression or speedup; models with substantial
implication chains are still needed for a meaningful benefit measurement.

## 2026-09-10 pre-reoptimization baseline

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

This result identified node LP cold-start behavior as the next target: the
model reported five cold solves and zero reoptimizations.

## 2026-09-10 bound-reoptimization result

After introducing the MIP relaxation layer and an original-coordinate
bound-only solve path, the same command and model produced:

| sample | LP mean | MIP mean | nodes | objective | bound |
|---:|---:|---:|---:|---:|---:|
| 1 | 12,691 ns | 21,459 ns | 5 | 3 | 3 |
| 2 | 18,385 ns | 23,602 ns | 5 | 3 | 3 |
| 3 | 25,471 ns | 28,059 ns | 5 | 3 | 3 |

Every sample used 7 LP iterations, 12 LP work units, 1 cold solve, and 4
reoptimizations. The structural result is the important invariant: after the
first node, every node reused the LP session. The timing improvement is large
on this smoke model, but it remains an environment-local result and is not a
claim of performance parity with HiGHS.

## 2026-09-10 linear-row propagation result

The first propagation stage computes conservative row activities, rounds
deduced integer bounds, persists deductions in the domain-change arena, and
reuses one scratch workspace for the entire solve. Final samples were:

| sample | LP mean | MIP mean | nodes | objective | bound |
|---:|---:|---:|---:|---:|---:|
| 1 | 19,715 ns | 19,995 ns | 3 | 3 | 3 |
| 2 | 20,601 ns | 23,322 ns | 3 | 3 | 3 |
| 3 | 21,553 ns | 24,791 ns | 3 | 3 | 3 |

All samples used 5 LP iterations, 8 LP work units, 1 cold solve, and 2
reoptimizations. Compared with the immediately preceding implementation, the
node count fell from 5 to 3 and LP work from 12 to 8 while preserving status,
objective, and bound. As above, timings are local smoke measurements; the
structural counters are the reproducible evidence.

## 2026-09-10 pseudocost feedback result

Directional pseudocost state now records objective degradation per unit branch
distance when a child LP is solved. The default branching rule consumes a
solve-local read-only score buffer; sparse node metadata carries only the
branch origin needed for feedback. Three final samples were:

| sample | LP mean | MIP mean | nodes | LP work | observations |
|---:|---:|---:|---:|---:|---:|
| 1 | 12,582 ns | 16,461 ns | 3 | 8 | 2 |
| 2 | 16,308 ns | 21,146 ns | 3 | 8 | 2 |
| 3 | 14,577 ns | 18,310 ns | 3 | 8 | 2 |

Every sample remained optimal with objective and bound equal to 3, one cold LP
solve, and two reoptimizations. The two observations match the two evaluated
children, providing an end-to-end invariant for the learning path. This model
has only one useful branch candidate, so it validates feedback overhead and
correctness rather than candidate-quality improvements.

## 2026-09-10 reliability-branching comparison

`test/model/data/reliability_mip.mps` fixes a two-candidate case. The
`branching` profile disables propagation so reliability probing is observable;
`branching-no-probes` changes only the probing switch. Representative
ReleaseFast results over 1,000 solves were:

| profile | MIP mean | nodes | LP work | observations | probes |
|---|---:|---:|---:|---:|---:|
| branching-no-probes | 22,951--34,879 ns | 5 | 12 | 2 | 0 |
| branching | 38,481--46,918 ns | 3 | 19 | 4 | 4 |

Both profiles returned optimal objective and bound 2. Reliability probing
proved branch directions infeasible and reduced evaluated tree nodes from 5 to
3, but its seven probe iterations increased total work and wall time on this
tiny instance. This is a correctness and accounting baseline, not evidence
that probing should always run. Candidate caps and per-node iteration/work
budgets are public controls; corpus-scale tuning remains required.

Run the comparison with:

```text
mip-solver-bench test/model/data/reliability_mip.mps 1000 branching-no-probes
mip-solver-bench test/model/data/reliability_mip.mps 1000 branching
```

## 2026-09-11 root cover-cut comparison

The `default` profile enables root separation; `no-cuts` changes only
`cut_separation`. On `test/model/data/tiny_mip.mps`, representative
ReleaseFast runs over 1,000 solves produced:

| profile | LP mean | MIP mean | nodes | LP work | cold / reopt | cut rounds / added |
|---|---:|---:|---:|---:|---:|---:|
| no-cuts | 12,066 ns | 14,259 ns | 3 | 8 | 1 / 2 | 0 / 0 |
| default | 15,366 ns | 22,374 ns | 1 | 7 | 2 / 0 | 2 / 1 |

Both profiles returned optimal objective and bound 3. The generated binary
cover inequality closes the root integrality gap, reducing the evaluated tree
from three nodes to one. Adding the row forces a second cold root solve because
the LP structure changes. On this tiny instance that setup cost outweighs one
unit of LP work and makes the cut-enabled wall time slower. This is why the
default remains a structural baseline rather than a blanket performance claim;
adaptive activation requires a representative corpus. Status, objective,
bound, nodes, work, and cut counts are the reproducible invariants.

Run the comparison with:

```text
mip-solver-bench test/model/data/tiny_mip.mps 1000 no-cuts
mip-solver-bench test/model/data/tiny_mip.mps 1000 default
```
