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
