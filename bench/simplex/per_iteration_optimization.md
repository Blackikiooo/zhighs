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

### H1g — validate FT diagonals at their write boundary (accepted)

Mutable-U diagonal values can change only in `reset` and `update`.  Updates
already reject zero/non-finite pivots; base diagonals are now checked while
`reset` copies their pivot IDs, allowing every subsequent `solveUpper` to
remove the redundant per-pivot diagonal validation.  The finite-result check
after division remains, and invalid reset input now fails earlier.

The forced-dual 40-model gate passed with all 6833 pivot counts unchanged.  A
31-run alternating full-corpus A/B measured `196.003 ms -> 194.164 ms`, or
`28684.7 -> 28415.7 ns/pivot` (`-0.94%`).  The independent 10-model discovery
gate was also positive at `-0.52%`, with improvements in `gas11`, `scorpion`,
`shell`, and `sc105`.

### H4a — short-circuit crash coefficient filtering (rejected)

The solve-entry coefficient audit was passed into the crash planner so normal
models could short-circuit four repeated coefficient-threshold checks without
duplicating the planner.  Matching decisions and pivot counts were unchanged,
but the 10-model 31-run gate regressed slightly from `102.256` to `102.431 ms`
(`+0.17%`): gains in `gas11`, `shell`, and `grow7` were outweighed by
`scorpion`, `scrs8`, and `seba`.  A runtime invariant still leaves a branch in
the deep loops and adds register pressure.  The change was removed; a future
crash optimization must use true outer specialization or reduce the greedy
planner's repeated full-column scans.

### H1h — sparse DSE tau FTRAN from retained Ep support (rejected)

The exact nonzero support of each pivotal BTRAN row was retained and passed to
the adaptive solve for `tau = B^-1 B^-T e_p`, matching HiGHS' indexed
`HVector` interface.  Representative trajectories stayed identical, but the
10-model 31-run gate regressed from `102.306` to `105.300 ms` (`+2.93%`),
including `scorpion +9.39%`, `capri +7.21%`, and `shell +4.16%`.  Despite
`shell`'s roughly 3% DSE Ep density, the current FT-aware hyper FTRAN pays more
for marking, scatter, and dense materialization than it saves in the L pass.
All support storage was removed.  Sparse DSE dispatch should be revisited only
after improving that lower-level FT hyper kernel.

### H3b — consume packed FTRAN support in pivot commit (rejected)

HiGHS `HEkkDualRHS::updatePrimal` updates the basic solution through the
indexed `col_aq` whenever its count is at most 40% of the row dimension, and
uses a dense pass only above that threshold.  zhighs already gathers the exact
nonzero FTRAN result into `changed_row_index`, but `performPivot` ignores that
support for the numerical update and then performs a second full-row pass to
scatter every basic value back into the global primal vector.

Change only the traversal ownership: use the existing natural-order support
for both operations when it contains at most 40% of the rows, otherwise retain
the current dense loops.  Every indexed update must use the same scalar
expression as the dense path, and the entering/leaving assignments remain in
their current order.  The prediction is unchanged pivot traces with a lower
`performPivot` cost on sparse-FTRAN models; the representative gate is
`shell`, `capri`, `scrs8`, `seba`, `grow7`, and `sc105`, followed by the frozen
40-model forced-dual gate if the aggregate signal is positive.

ReleaseFast tests passed and every representative pivot count, status, and
objective remained identical.  The 31-run alternating gate nevertheless
regressed by `+6.67%` in aggregate: `shell +17.88%`, `capri +15.19%`,
`sc105 +1.00%`, `grow7 -0.28%`, `scrs8 -1.40%`, and `seba -2.22%`.
Two contiguous row passes are cheaper than two indirect support traversals in
the dominant regimes, despite FTRAN sparsity.  The implementation was removed.
This rules out transplanting HiGHS' 40% consumer threshold onto a support that
zhighs reconstructs after the solve; the next alignment step must make sparse
support producer-owned across the factor solve and its consumers.

### H1i — publish exact entering support inside FT upper solve (protocol)

HiGHS `HFactor::ftranFT` extends `HVector.index` while applying the mutable
upper solve.  zhighs instead materializes the dense direction, scans all rows
to collect exact nonzeros, and scans the base pivot order again whenever the
result is at most 10% dense.  Extend only the indexed entering-column FT solve
so `solveUpper` publishes every exactly nonzero solved pivot while that value
is already hot.  For sparse outputs, restore the existing reverse base-pivot
order before exposing the support; for outputs above 10%, publish no support
and retain the current natural-order scan.  This preserves all consumer order
and numerical expressions.

The prediction is identical pivot traces with one or two dimension-sized
post-FTRAN scans removed on sparse entering directions.  Validate first on the
six representative models with 31-run alternating A/B, then run the frozen
40-model forced-dual gate if the signal is positive.
