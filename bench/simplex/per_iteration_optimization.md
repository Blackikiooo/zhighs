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

### H2c — specialize PRICE after the coefficient audit (accepted)

HiGHS removes small matrix entries before simplex and its column PRICE inner
loop therefore has no coefficient-retention branch.  zhighs already audited
the complete matrix at solve entry, so the result is now retained and selects
one of two compile-time-specialized PRICE kernels: the common branch-free path
and the original filtered path for models that actually contain dropped
coefficients.  Numerical ordering is unchanged for both paths.

The forced-dual 40-model correctness gate passed with all 6833 pivot counts
unchanged.  A 31-run alternating full-corpus A/B measured `197.034 ms ->
195.402 ms`, or `28835.7 -> 28596.9 ns/pivot` (`-0.83%`).  The representative
gate showed the intended PRICE-heavy improvements, including `shell -3.12%`,
`scrs8 -2.52%`, `grow7 -2.17%`, and `capri -1.78%`.

### H1f — specialize entering-column assembly (rejected)

The same solve-entry audit was applied to compile-time-specialized dense and
sparse FTRAN column assembly, mirroring HiGHS `collectAj`.  Pivot counts stayed
unchanged, but the 10-model 31-run alternating gate regressed from `102.626`
to `103.033 ms` (`+0.40%`), including `gas11 +0.96%` and `scfxm1 +1.67%`.
Entering columns are short enough that the predictable retention branch is
cheaper than the extra helper dispatch and code footprint.  The experiment
was removed while the accepted PRICE specialization remains isolated.

### H2d — fuse DSE norm into logical tableau packing (rejected)

The exact DSE norm and row-density pass was fused into the later natural-order
logical-tableau packing pass, preserving arithmetic order and all thresholds.
Despite deleting one dense row read, the 10-model 31-run gate regressed from
`103.310` to `104.038 ms` (`+0.70%`), led by `scfxm1 +2.63%`, `scorpion
+1.70%`, and `capri +1.26%`.  Moving the substantial DSE validation block into
the PRICE/CHUZC tail harmed hot-code layout enough to dominate the saved read.
The original phase ordering was restored.

### H2e — fuse DSE norm with column-PRICE scaling (accepted)

Keeping the DSE validation block in its original location avoids H2d's code
layout regression.  On column-PRICE pivots, the mandatory exact-norm scan now
also materializes the row-scaled BTRAN vector that PRICE consumes next; Devex
and row-PRICE retain their prior loops.  The squared-norm accumulation order,
support threshold, and validation timing are unchanged.

The forced-dual 40-model gate passed with all 6833 pivot counts unchanged.  A
31-run alternating full-corpus A/B measured `193.680 ms -> 192.338 ms`, or
`28344.8 -> 28148.5 ns/pivot` (`-0.69%`).  The 10-model discovery gate was
also positive (`-0.42%`), including `scorpion -2.24%`, `sc105 -1.91%`, and
`shell -1.04%`.

### H3a — exact-nonzero sparse reduced-cost recurrence (rejected)

PRICE published every exactly nonzero tableau entry into a dedicated
natural-order support, allowing the reduced-cost recurrence to skip only
strict `alpha == 0` terms.  This corrected H1a's tolerance-truncated support:
all representative pivot counts and statuses remained identical.  Performance
still regressed sharply in the 10-model 31-run gate, from `101.666` to
`103.850 ms` (`+2.15%`), including `brandy +5.53%`, `capri +3.35%`, and
`scrs8 +3.34%`.  For these dimensions, contiguous dense recurrence beats
indirect sparse reads even when the support is exact.  The index allocation
and sparse update were removed; future work on this pass must preserve
contiguous traversal, for example by fusion or proven SIMD generation.
