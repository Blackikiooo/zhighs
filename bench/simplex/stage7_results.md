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

## Open acceptance gates

- Acquire and lock the five official Netlib special/generated cases.
- Add a pinned CLP runner and repeat status/objective/certificate comparison.
- Fix or explicitly classify the three baseline numerical failures and the
  `tuff` Phase-I cycle; rerun the entire locked corpus after every fix.
- Finish 60-second classification for the large models, then choose the final
  Netlib timeout from evidence rather than converting timeout to failure.
- Acquire and lock the Mittelmann corpus, set both timeout and memory caps, and
  run the same raw-output harness.
- Produce repeated-run median/p95 timing and requested-bytes/peak-RSS tables
  only after the correctness gate is stable.
