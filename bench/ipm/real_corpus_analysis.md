# IPM 40-model corpus: end-to-end validation and comparison with HiGHS 1.12.0

Environment: Windows 10, zig 0.16.0, `ReleaseFast` `ipm-corpus` runner, HiGHS
1.12.0 (highspy) with `threads=1`, `parallel=off`, `run_crossover=off`,
presolve on for the simplex reference and off for the IPM reference. Corpus:
the 40 frozen SoPlex `check/instances` LPs (netlib-era models, 8x8 to
536x1775, hash-locked against `bench/simplex/end_to_end_corpus.lock.tsv`).
Fresh results: `real_corpus_results_fresh.json` (summary:
`python bench/ipm/analyze_results.py <json>`).

## Headline numbers

| Metric | Value |
|---|---|
| Models solved to `optimal` | 27/40 (HiGHS IPM: 38/40; two `Solve error`s are box1/ex72a) |
| Correct optimal status + full independent KKT validation | 18/40 |
| Reference simplex disagrees on any objective | 0/30 (all validated optima match HiGHS simplex objectives within 1e-9 relative) |
| zhighs/HiGHS-IPM time, 26 comparable optimal-optimal | geometric mean 0.63x (faster on 17/26), best 0.063x, worst 56.8x |
| zhighs solve time spent in the Cholesky factor | 93.4% |

## Failure taxonomy (13 no-answer + 9 accuracy-miss models)

1. **Infeasible/unbounded models (8): bgetam, forest6, galenet, gams10am,
   gas11, klein1, refinery, woodinfe.** HiGHS simplex reports 7 infeasible + 1
   unbounded; HiGHS IPM also detects every one of them. zhighs reports
   `numerical_failure` because the engine only has divergence heuristics, no
   infeasibility certificate and no self-dual embedding. This is the largest
   single capability gap.
2. **Iteration limit (2): box1, ex72a.** Infeasible models that neither
   converge nor trip the divergence bound within 100 iterations. HiGHS IPM
   errors on these too, so the bar is "detect", which zhighs currently
   cannot.
3. **Numerical failure on feasible models (3): brandy, capri, lotfi.** The
   reference is Optimal; zhighs exits via divergence detection after 13-23
   iterations. Regularization tops out at 1e-6 and there is no iterative
   refinement, so ill-conditioned late-phase normals are fatal.
4. **Optimal but outside strict validation thresholds (9): agg, beaconfd,
   bore3d, grow7, scagr25, scfxm1, share1b, shell, vtp-base.** All nine fail
   only `primal_rel` (1e-7 bar): raw primal violations land at 1e-7..1e-5
   scale. Objectives still match the reference within ~1e-9 relative. Without
   crossover this is the classic interior-point "slightly-interior" answer;
   HiGHS produces crossover-clean vertices by default.

## Where zhighs is faster

Every model with a friendly normal-equation envelope (m below ~250 or
near-banded structure) beats HiGHS IPM by 3-16x: afiro 0.07x, sc50a 0.06x,
scsd1 0.12x, share2b 0.09x, sc205 0.19x. These solves run in 0.1-3 ms — the
assembly and solve loops are tight, allocations are hoisted once per solve,
and there is no presolve or crossover overhead in the timed path.

## Where HiGHS is faster, and why

The profile (envelope) Cholesky stores the full row envelope of A·Aᵀ with no
ordering, so fill is decided by raw row order:

| model | rows | A nnz | stored normal entries | stored/m² |
|---|---|---|---|---|
| seba | 515 | 4,352 | 364,599 | 1.375 |
| grow7 | 140 | 2,612 | 63,035 | 3.216 |
| etamacro | 400 | 2,409 | 120,196 | 0.751 |
| bgetam | 400 | 2,409 | 120,196 | 0.751 |
| finnis | 497 | 2,310 | 113,064 | 0.458 |

seba stores 84x the nonzeros of A (and 1.375x a fully dense triangle: the
envelope exceeds m²/2 by a third). The factor phase then owns 93.4% of solve
time, and the same six slowest models are exactly the six worst fill ratios:
seba 56.8x, etamacro 21.3x, agg 22.3x, finnis 14.5x. Iteration counts on
these models are equal to HiGHS (seba 16 vs 17), so the entire deficit is
the linear-algebra kernel. HiGHS/IPX factors the same systems with an
approximate-minimum-degree ordering and a sparse up-looking LDLᵀ whose
factor is typically a small multiple of nnz.

## Conclusion vs open-source state of the art

The Mehrotra predictor-corrector iteration itself is competitive: comparable
iteration counts on every model that converges, and a 0.63x geometric-mean
solve time. The gaps, in priority order:

1. **P0 — ordering for the normal-equation factorization.** Replace the
   unordered envelope with an AMD/minimum-degree permutation + sparse
   LDLᵀ (supernodal or up-looking). This alone should move the worst six
   models from 14-57x slower to parity or better, since fill collapses from
   O(m²) to O(nnz·small factor).
2. **P0 — infeasible/unbounded detection.** A self-dual embedding (IPX's
   approach) or explicit infeasibility certificates. Eight models currently
   exit as `numerical_failure` that HiGHS resolves exactly.
3. **P1 — numerical robustness.** brandy/capri/lotfi need either iterative
   refinement of the normal-equation solve, a less conservative
   regularization schedule, or both.
4. **P1 — crossover (or tighter termination).** Nine models lose acceptance
   on 1e-7..1e-5 primal residuals only; a crossover to a basic solution
   would close the gap to the simplex-layer accuracy contract.