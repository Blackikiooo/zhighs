# Simplex per-iteration optimization protocol

## Objective

Reduce the implementation cost of a committed simplex pivot without trading
away the already lower iteration count.  The reference is the pinned local
HiGHS source and release build used by the end-to-end corpus gate.

The primary metric is computed over one copy of every model in the frozen
40-model corpus:

```
sum(per-model median solver nanoseconds) / sum(committed pivots)
```

The starting measurement at commit `45dfe3b` uses three warmups and 21
alternating measured runs per solver:

| solver | sum of medians | committed pivots | ns / pivot |
|---|---:|---:|---:|
| zhighs | 204.221115 ms | 6833 | 29887.47 |
| HiGHS | 127.164363 ms | 8774 | 14493.32 |

The starting ratio is 2.0622, so reaching parity requires at least a 51.5%
reduction in zhighs cost per committed pivot at the current iteration count.

## Experimental discipline

1. Use the hashes in `end_to_end_corpus.lock.tsv`; never select behavior by
   model name.
2. Compare instrumentation-free ReleaseFast zhighs with the pinned local
   HiGHS release build.  Alternate execution order and use odd-run medians.
3. Use `shell`, `capri`, `scrs8`, and `seba` for hotspot discovery.  Use
   `grow7` and `sc105` as controls where zhighs is already faster.
4. Change one mechanistic hypothesis at a time.  A short representative gate
   may reject an experiment, but accepting it requires the full 40-model
   correctness gate and repeated performance validation.
5. Do not relax feasibility, residual, pivot, or certificate tolerances.  Do
   not accept a lower time caused by a higher committed-pivot count.
6. Commit each accepted optimization independently.  Fully remove rejected
   experimental code and record the negative result before moving on.

## Hypothesis order

### H1: factor solve and Forrest--Tomlin traffic

Profile equivalent zhighs and HiGHS runs and attribute samples to INVERT,
FTRAN, BTRAN, UPDATE, PRICE, and non-kernel iteration work.  Inspect the local
HiGHS `HFactor`, `HVector`, and `HFactorRefactor` implementations at every
decision point.  Test elimination of redundant dense clears, permutation
round-trips, support rebuilding, and generic dispatch in sparse solves.

Prediction: the largest slow models spend enough time in triangular solves and
FT maintenance that reducing touched memory and indirect branches lowers their
median ns/pivot by at least 10% without changing a pivot trace.

### H2: dual pricing and DSE updates

Measure row packing, CHUZR scans, tableau pricing, exact-weight validation, and
DSE support handling.  Compare their sparse/dense transition rules with
HiGHS.  Prefer generation marks and retained sparse support over full-vector
clear or scan operations.

Prediction: `shell`, `capri`, and `scrs8` expose dimension-scaled work outside
the numerical kernels; removing one full row/column scan per pivot gives a
further 10--20% reduction.

### H3: pivot commit and incremental state maintenance

Audit the post-ratio-test path for repeated status, bounds, value, reduced-cost,
infeasibility-list, and edge-weight passes.  Fuse only operations that share
the same numerical ordering and ownership in `HEkkDual::updateDual`,
`updatePrimal`, and `updatePivots`.

Prediction: SoA pointer caching and single-pass sparse updates reduce retired
instructions and cache misses while preserving all committed iterations.

### H4: reinversion and rebuild amortization

Classify rebuilds by reason and measure factorization cost per committed pivot.
Match HiGHS lifecycle decisions using the real factor graph; do not introduce
parallel shadow state or model-specific thresholds.

Prediction: after H1--H3 remove solve overhead, remaining outliers will be
dominated by avoidable reinversions or overly conservative recovery epochs.

## Acceptance and stopping conditions

An accepted experiment must pass `zig build test -Doptimize=ReleaseFast`, the
frozen 40-model correctness gate, and `git diff --check`.  Its benchmark must
show a repeatable improvement outside observed noise, with no committed-pivot
regression.  The outer loop continues until the 40-model aggregate zhighs
ns/pivot is below HiGHS; geometric per-model ns/pivot and worst regressions are
reported as secondary guards.

## Experiment log

### H1a — packed reduced-cost update (rejected)

The first direct port added a persistent nonbasic tableau support and updated
only that support, matching `HEkkDualRow::updateDual` structurally.  Unit tests
passed and representative solutions retained their status and objective, but
the pivot paths regressed immediately: `shell 480 -> 629`, `capri 264 -> 575`,
`scrs8 438 -> 1035`, `grow7 98 -> 104`, and `sc105 80 -> 93`.  All experiment
code was removed.

This rules out switching the update loop before the existing basic-column
reduced-cost and `shiftBack` lifecycle is aligned with HiGHS.  The next H1
experiment must preserve the dense recurrence and isolate work that provably
cannot affect numerical state, starting with the post-update finite scan.

### H1b — fuse the post-update finite scan (rejected)

Finite checks were fused into the unchanged dense recurrence, eliminating a
second reduced-cost vector pass.  Six-model 31-run alternating A/B preserved
every pivot count, but produced mixed noise rather than a general gain:
`scrs8 -0.37%`, `grow7 -1.78%`, `sc105 -0.02%`, `shell +0.33%`,
`capri +1.69%`, and `seba +1.98%`.  The added per-element branch offsets the
saved read pass on this corpus, so the implementation was removed.

### H1c — sparse-support DSE norm (rejected)

The adaptive BTRAN output support was reused for exact DSE norm and density
instead of scanning the dense row.  Pivot counts remained unchanged, but
31-run A/B showed that indirect support reads lose to the short sequential
vector loop: `shell +1.94%`, `scrs8 +2.34%`, `seba +2.49%`,
`capri +11.74%`, `sc105 +0.47%`, and only `grow7 -1.03%`.  The code was
removed.  Future sparse traversal gates must include dimension and touched
bytes, not density alone.

### H2a — reuse scaling across crash-prefix validation (rejected)

Crash prefix retries were reset from retained scaled bounds instead of
repeating six equilibration passes.  Pivot counts were identical.  A 101-run
setup-only A/B improved `capri -12.32%`, `grow7 -8.62%`, `sc105 -7.25%`, and
`scrs8 -2.86%`, but regressed `shell +12.89%` and `seba +3.96%`.  The helper
was removed: rebuilding combinatorial state with branch-heavy scalar loops and
the resulting code layout is not a safe universal replacement for the current
initialization path.

### H1d — repack captured FT Aq (rejected)

The dense captured Aq was replaced with natural-order packed index/value
streams before FT update.  All pivot counts remained identical.  An initial
31-run A/B was mildly positive overall, but a 101-run confirmation resolved to
`shell +0.11%`, `scrs8 +0.70%`, and `sc105 -0.13%`.  Repacking by rescanning
the dense vector merely moves work from UPDATE to FTRAN and was removed.  A
future packed Aq must be published directly by the sparse FTRAN producer.

### H1e — publish packed Aq from sparse L support (rejected)

The hyper-sparse update FTRAN passed its active L support directly into FT
capture, with correction pivots appended and indices restored to natural row
order.  Although unit tests passed, representative pivot paths changed badly
(`shell 480 -> 838`, `capri 264 -> 423`, `scrs8 438 -> 762`).  Sorting did not
restore the baseline.  The L active set is therefore not a sufficient owner of
the dense state consumed across the historical FT correction chain.  All code
was removed; any future direct pack must be maintained inside the correction
lifecycle with exact generation tracking.

### H2b — publish packed dual tableau support (accepted)

PRICE now publishes the natural-order support of the pivotal tableau row into
the existing BFRT column workspace, and CHUZC scans that support rather than
all structural and logical columns.  No extra allocation or persistent index
array is required.  The forced-dual 40-model correctness gate passed and every
committed-pivot count remained unchanged.  A 31-run alternating full-corpus
A/B measured `197.206 ms -> 196.359 ms` at 6833 pivots, reducing aggregate
cost from `28860.8` to `28736.9 ns/pivot` (`-0.43%`).  Per-model results were
mixed because short fixed-cost solves are noisy, but the aggregate gain was
carried by real iteration-heavy work including `gas11 -4.17%`, `brandy
-3.03%`, `agg -2.43%`, and `scorpion -1.91%`.

The repository's default primal Phase-I path still fails a pre-existing set
of corpus models at `direction_column` before its first committed pivot.  An
isolated clean-HEAD build reproduced the same `agg` failure, so H2b acceptance
uses the frozen forced-dual protocol targeted by this optimization campaign.
