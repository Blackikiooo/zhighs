# Stage 7 end-to-end acceptance

This report tracks the long-running acceptance work from
`src/lp/simplex/todo.md` section 7. It is intentionally separate from the
40-model fast correctness gate: failures found here remain visible until they
are fixed or explicitly classified.

## Reproducible setup

- Date: 2026-07-19.
- Host: WSL2 Linux 6.6.87.2, Intel Core i9-10900KF (10 cores / 20 threads).
- zhighs: Zig 0.16.0, `ReleaseFast`, serial simplex, presolve unavailable.
- HiGHS: commit `de09bbad9fb7c5d39a1a464a7641bbb5531c6e9d`, `-O3
  -march=native -DNDEBUG -flto`, serial simplex, presolve disabled.
- C++ compiler: GCC 13.3.0.
- Initial process timeout: 10 seconds per solver/model. Exceptions were rerun
  with a 60-second process timeout.

`stage7_netlib.lock.tsv` pins SHA-256 for the 93 traditional compressed MPS
instances decoded with Netlib's official `lp/data/emps.c`. The official index
also lists `stocfor3` and `truss` as special generated/bundled entries and
QAP8/QAP12/QAP15 as “see NOTES”; those five acquisition exceptions are not yet
part of the 93-file lock and must not be silently counted as tested.

Run a selected acceptance set with:

```sh
TIMEOUT_SECONDS=60 OUTPUT_FILE=/tmp/stage7.tsv \
  bench/simplex/run_stage7_corpus.sh /path/to/netlib-mps afiro brandy
awk -f bench/simplex/summarize_stage7.awk /tmp/stage7.tsv
```

`MEMORY_LIMIT_KB` applies a virtual-memory cap, while the harness records GNU
time peak RSS for every solver. `CLP_RUNNER=/path/to/compatible-runner` enables
the third comparison once a pinned CLP runner is available. CLP was not
installed for this initial run, so no CLP result is claimed here.

## Parser gate found by the full corpus

The initial pass exposed two data-compatibility defects that the fast corpus
did not cover:

1. Traditional fixed MPS permits blanks inside an eight-character row or
   column name. Whitespace tokenization split FORPLAN names such as
   `DEDO3 1R`; section-aware fixed-field extraction now preserves them.
2. Fixed MPS may omit the BOUNDS vector name. Treating the column name as a set
   name silently discarded bounds in `gfrd-pnc` and `sierra`; the parser now
   recognizes value-requiring bound codes and the omitted-set layout.

After these fixes all 93 locked inputs parse. The corrected `gfrd-pnc` and
`sierra` objectives match HiGHS (`6902235.999548811` and
`15394362.183631929`). Unit tests cover both fixed-format cases.

## Initial Netlib result

With the 10-second cap:

| Solver | Optimal | Numerical failure | Timeout |
| --- | ---: | ---: | ---: |
| zhighs | 84 | 3 | 6 |
| HiGHS | 91 | 0 | 2 |

All 84 models completed as optimal by both solvers match objective within
`1e-7 * (1 + abs(reference))`; there are zero objective mismatches. Across
those common completions, the single-pass total-time distribution was:

| Metric | zhighs | HiGHS | zhighs / HiGHS |
| --- | ---: | ---: | ---: |
| median | 39.19 ms | 25.39 ms | 1.21x |
| p95 | 5.019 s | 0.831 s | 8.68x |

These are corpus-distribution percentiles from one run per model, not warmed
repeated-run benchmark percentiles. Stage 7's final performance table still
requires warmups and repeated samples.

### Timeout and failure isolation

- `d2q06c` completes in 21.35 s after 100,941 iterations; HiGHS completes
  inside 10 s.
- `d6cube` completes in 11.80 s after 101,442 iterations; HiGHS completes
  inside 10 s.
- `dfl001` and `fit2p` still exceed 60 s.
- `pilot87` reaches `numerical_failure/optimality_check` after 32,065
  iterations and 33.56 s. Its maximum dual violation is `2.47e-3`.
- `tuff` reaches the one-million-iteration limit entirely in Phase I. Of those
  pivots, 999,999 are classified degenerate; this is a cycling/anti-degeneracy
  failure rather than a sparse-kernel throughput result.
- `modszk1` stops at `pivot_factorization` after 18,683 iterations. It records
  355 growth-triggered and 136 solve-residual reinversions.
- `scsd8` stops at `optimality_check` with primal violation `2.905478` after
  5,028 iterations.
- `wood1p` stops at `pivot_factorization` after 2,281 iterations, with maximum
  FTRAN relative residual `4.12e-5`.

HiGHS also exceeds the initial 10-second cap on `dfl001` and `pilot87`; its
remaining 91 locked models complete as optimal.

### Degeneracy strategy A/B

The forced perturbation mode repairs two baseline failures without changing
their reference objective:

- `scsd8`: optimal in 4,648 iterations, objective `904.9999999254647`.
- `wood1p`: optimal in 598 iterations, objective `1.4429024115734086`.

It does not repair `modszk1`, `pilot87`, or the one-million-pivot `tuff` loop.
Forced taboo mode produces false infeasible conclusions on `tuff` and
`wood1p`, so this Stage 7 evidence explicitly rejects enabling taboo by
default. The next change must fix or guard its certificate path before any
automatic dispatch is considered.

## Validated bounded perturbation update

The initial A/B above is retained as the discovery baseline. A subsequent
implementation changed perturbation ranks from tie-order-only metadata into
bounded positive primal margins, while preserving three correctness guards:

1. taboo exhaustion is retried without exclusions before Phase I may stop;
2. perturbed Phase-I infeasibility restores the logical basis and reinstalls
   the exact same artificial-column construction used by the initial solve;
3. failed optimal/unbounded validation cold-restarts the baseline policy.

Automatic mode waits for 256 consecutive degenerate pivots before activating;
the explicit perturb mode retains the earlier trigger for controlled A/B.
This separates long faces from frequent short local ties.

Targeted ReleaseFast results are:

| Model | Previous default | Automatic bounded perturbation |
| --- | ---: | ---: |
| `tuff` | 1,000,000 iteration limit | optimal, 1,068 iterations, 29 ms |
| `modszk1` | factorization failure at 18,683 | optimal, 6,448 iterations, 242 ms |
| `brandy` | 3,384 iterations | optimal, 1,519 iterations, 23 ms |
| `scsd8` | optimality-check failure | optimal, 3,655 iterations, 133 ms |
| `wood1p` | factorization failure at 2,281 | optimal, 802 iterations, 90 ms |
| `d6cube` | 101,442 iterations, 11.9 s | 61,506 iterations, 6.39 s |
| `d2q06c` | 100,941 iterations, 21.3 s | 98,703 iterations, 17.9 s |

The 40-model status/objective/residual/ray gate passes in automatic mode. A
fresh 93-model run with a 10-second process cap reports:

| Solver/policy | Optimal | Numerical failure | Timeout |
| --- | ---: | ---: | ---: |
| previous zhighs baseline | 84 | 3 | 6 |
| zhighs automatic bounded perturbation | 88 | 0 | 5 |
| pinned HiGHS reference | 91 | 0 | 2 |

All 88 completed zhighs objectives match the pinned HiGHS results at `1e-7`
relative tolerance. Maximum published primal and dual residuals are
`9.93e-8` and `5.30e-8`. The remaining zhighs timeouts are `d2q06c`,
`dfl001`, `fit2p`, `pilot`, and `pilot87`; HiGHS also times out on `dfl001`
and `pilot87` under the same initial cap, leaving three zhighs-only long-tail
cases. Single-pass corpus median and p95 total times are 43.20 ms and 5.36 s.

An isolated attempt to update only the new nonbasic leaving-column Devex
weight increased `d6cube` to 136,676 iterations and `d2q06c` to 114,097. It
was rejected and removed: a production replacement must implement the full
Devex or projected steepest-edge recurrence, not a one-column approximation.

## Batched dual bound flips

The next long-tail pass replaces one FTRAN per accepted dual bound flip with
one allocation-free FTRAN of the accumulated equation-space displacement
`sum(A_j * delta_j)`. Variable states and basic values are committed only
after the aggregate solve succeeds. Dedicated tests verify that two flips
produce the same basic displacement while issuing one FTRAN, and the 40-model
status/objective/residual/ray gate passes in automatic degeneracy mode.

As a targeted exercised dual path, forced dual Phase I on `scsd1` remains
optimal at objective `8.666666674333364`, with primal residual `5.55e-17` and
dual residual `3.52e-8`. It records 36 flips in 33 batches, saving three
FTRAN calls. The modest reduction is expected because this model's batches
are predominantly single flips; the new counters permit corpus-wide
attribution before any further default-policy change.

## Incremental dual reduced costs

Phase-II `solveDual` now retains the old-basis pivotal tableau row and applies
the exact rank-1 recurrence `r' = r - (r_q / alpha_pq) alpha_p`. This removes
the previous full BTRAN and matrix reprice, followed by a complete dual
feasibility scan, after every dual pivot. The entering reduced cost is pinned
to zero, while logical and artificial entries are updated in the same scaled
coordinates as structural columns.

The incremental vector is never accepted as a final certificate. An exact
reprice runs after eight updates or immediately after a fresh factorization,
records normalized drift, and revalidates dual feasibility. Adaptive mode may
shorten the refresh period when drift exceeds its tolerance; terminal
`finishOptimal` independently recomputes basic values and reduced costs before
publishing the solution. Formula-level and imported-basis dual integration
tests pass, including deliberate drift injection and exact recovery. The
40-model default status/objective/residual/ray gate remains green. Counters
`dual_reduced_cost_updates` and `dual_exact_reprices` were added to the stats
runner. Forced dual Phase I on `scsd1` correctly reports zero for these
Phase-II-only counters because all 181 pivots occur before Phase II.

## Full primal Devex reference framework

An explicit `legacy/framework` policy now separates the frozen default from a
full primal Devex implementation. The framework stores one byte per internal
column for the frozen nonbasic reference set. Each basis replacement computes
the pivotal reference norm from the existing FTRAN direction, divides it by
the squared pivot, and applies the resulting lower bound to every nonbasic
weight in the complete old-basis tableau row. Four accumulated overestimated
entering weights trigger a deterministic framework rebuild after the pivot is
committed. This is the complete recurrence, not the previously rejected
leaving-column-only approximation.

Both policies pass the 40-model status/objective/residual/ray gate. Framework
mode gives the following 21-run ReleaseFast result on `brandy`:

| Metric | Legacy median / p95 | Framework median / p95 |
| --- | ---: | ---: |
| iterations | 1,519 | 498 |
| total | 22.87 / 23.50 ms | 9.70 / 10.08 ms |
| Phase I | 17.79 / 18.44 ms | 7.03 / 7.44 ms |
| Phase II | 7.39 / 7.68 ms | 3.44 / 3.64 ms |
| PRICE | 5.88 / 6.06 ms | 2.91 / 3.02 ms |
| FTRAN | 4.52 / 4.71 ms | 1.65 / 1.89 ms |
| BTRAN | 2.48 / 2.56 ms | 1.22 / 1.32 ms |
| requested bytes | 418,233 | 426,441 |
| peak RSS | 2,016 KiB | 2,012 KiB |

One instrumented framework solve uses seven frameworks, performs 372 full
weight updates, and records 23 bad entering weights. The 40-model run also
shows non-universal iteration changes: `bore3d` 252 -> 333, `scorpion`
512 -> 619, and `seba` 564 -> 704, while many larger paths improve sharply.
The mode therefore remains explicit until the missing `d2q06c/d6cube` inputs
and the full Stage 7 corpus are rerun; it is not yet the release default.

## Dual DSE to Devex fallback

The explicit `steepest-devex` dual edge-weight policy now starts each dual
phase with exact BTRAN-norm DSE weights. A Huangfu pivotal-weight rejection or
a deterministic recurrence-update budget switches the phase transactionally
to unit-weight dual Devex. Dual Devex applies the full HiGHS recurrence to all
row weights using the already hot FTRAN column; the fallback allocates no
memory and performs no additional basis solve. The inherited default pricing
policy is unchanged.

Unit coverage includes a real imported-basis dual pivot with a one-update
budget, direct invalid-framework recovery, and the full row-weight recurrence.
The default 40-model status/objective/residual/ray gate remains green. A forced
dual Phase-I `scsd1` correctness A/B produced:

| Policy | iterations | solve time | DSE updates | Devex updates | fallbacks |
| --- | ---: | ---: | ---: | ---: | ---: |
| inherited dual Devex | 181 | 6.59 ms | 0 | 0 | 0 |
| steepest-devex, budget 64 | 115 | 3.84 ms | 64 | 51 | 1 budget |

These are single ReleaseFast diagnostic runs, not a repeated timing claim.
The explicit result passed the locked status, objective and residual checks;
the policy remains non-default until a full Stage 7 A/B and repeated timing
run establish that its extra exact DSE initialization is worthwhile broadly.

## Open acceptance gates

- Acquire and lock the five official Netlib special/generated cases.
- Add a pinned CLP runner and repeat status/objective/certificate comparison.
- Reduce the three zhighs-only 10-second long tails: `d2q06c`, `fit2p`, and
  `pilot`. The former already improves from 100,941 to 98,703 iterations but
  remains outside the initial cap.
- Run the full Devex framework on the Stage 7 corpus, especially
  `d2q06c/d6cube`, before deciding whether it replaces the legacy default.
- Run the explicit dual `steepest-devex` policy on the full Stage 7 corpus and
  collect repeated median/p95 data before selecting a default update budget.
- Add reversible LP presolve (at minimum fixed columns, empty rows/columns and
  singleton rows) with primal/dual/ray/certificate postsolve validation.
- Finish 60-second classification for the large models, then choose the final
  Netlib timeout from evidence rather than converting timeout to failure.
- Acquire and lock the Mittelmann corpus, set both timeout and memory caps, and
  run the same raw-output harness.
- Produce repeated-run median/p95 timing and requested-bytes/peak-RSS tables
  only after the correctness gate is stable.
