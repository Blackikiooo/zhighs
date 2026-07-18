# Simplex end-to-end acceptance

The runner includes MPS parsing, canonical CSC construction, and a complete
simplex solve. HiGHS 1.14.0 is run with serial simplex and presolve disabled.
Times below are one ReleaseFast smoke run and are diagnostic, not stable
performance claims.

| model | zhighs status / iterations | HiGHS status / iterations | objective | zhighs primal / dual residual |
|---|---:|---:|---:|---:|
| afiro | optimal / 16 | Optimal / 22 | -464.75314285714285 | 5.68e-14 / 9.99e-16 |
| adlittle | optimal / 129 | Optimal / 74 | 225494.9631623804 | 9.09e-13 / 7.42e-11 |
| sc50a | optimal / 48 | Optimal / 55 | -64.5750770585645 | 1.14e-13 / 1.67e-16 |
| sc105 | optimal / 118 | Optimal / 124 | -52.202061211707225 | 5.40e-13 / 1.39e-16 |
| brandy | optimal / 4561 | Optimal / 304 | 1518.5098964881295 | 2.41e-11 / 1.24e-14 |

`sc105` exposed an FT construction error: `captureEp` applied historical row
corrections after solving the current mutable `U^-T`. The correction was then
counted again by the next update. Separating the pure mutable-U solve from full
BTRAN fixes the default run; a pre-ratio-test `B*aq-a` check remains as a
backward-error safety gate.

The `brandy` pivot hook showed that the default FT and reinvert-after-every-
pivot oracle have an identical 1025-event prefix. Earlier divergences at
iterations 50 and 541 were unstable Bland leaving pivots and are prevented by
a direction-scaled pivot threshold. The remaining cycle came from applying a
standard-form Bland fallback to generalized bounded Phase I with artificial
columns. Phase I now retains Devex/Harris candidate selection and clears the
fallback state when transitioning objectives. This restores correctness;
closing the large iteration-count gap to HiGHS remains performance work.

Run the deterministic smoke corpus with:

```sh
bench/simplex/run_end_to_end_corpus.sh
```

## Expanded local Netlib pass

The 40-model SoPlex LP corpus exposed two fixed-MPS compatibility defects:
section keywords such as `RHS` are legal set names, and a traditional fixed
record may omit the RHS/RANGES set field. Both are now parsed without losing
the records. Models including `forest6`, `sc205`, `scagr7`, `scrs8`, and
`lotfi` consequently changed from false zero-RHS results to the same status
and objective as HiGHS.

After the parser and numerical fixes, all 40 models in the frozen corpus match
the HiGHS reference status, and all optimal objectives pass the acceptance
tolerance. The four former numerical failures now finish as follows:

| model | status | objective | iterations | primal / dual residual |
|---|---:|---:|---:|---:|
| blend | optimal | -30.81214984582824 | 1416 | 5.68e-14 / 8.60e-16 |
| grow7 | optimal | -47787811.8147115 | 394 | 9.37e-10 / 1.06e-13 |
| scsd1 | optimal | 8.666666674333367 | 559 | 5.55e-17 / 1.49e-8 |
| vtp-base | optimal | 129831.46246136136 | 718 | 1.23e-11 / 2.18e-11 |

`grow7` exposed an incorrectly scaled FTRAN residual: `||a_q-Bx||/||a_q||`
rejects a backward-stable solve when large `B*x` terms cancel. The fused
residual traversal now accumulates `|a_q| + |B||x|` without allocation.
`blend` and `scsd1` additionally require fresh-factorization validation for
pivotal components near the forward-accuracy boundary. Once encountered, the
solve retains fresh factorizations because returning to intermittent FT updates
can change later ratio-test decisions on these ill-conditioned bases.

`gas11` is now independently certified unbounded. Direct HiGHS 1.14.0 with
presolve disabled reports Unbounded in 699 iterations and about 0.05 seconds.
The zhighs run performs 866 iterations in about 0.022 seconds on the same host;
its published primal ray has maximum bound/row directional violation
`1.71e-13` and objective direction `-3.600008e7`.

The main speedup is structural rather than allocator-related. Phase-I reduced
costs previously cleared and dotted a rows-long dense vector for every
internal column, about 835 million scalar visits on gas11. A direct CSC
`c - A^T y` pass reduces this to O(nnz + rows + columns). Reduced costs are
then maintained incrementally with periodic exact sparse refreshes. Model
coefficients at or below `1e-9` are consistently dropped before scaling,
matching the HiGHS diagnostic policy for this model.

## Frozen pre-optimization gate

The 2026-07-18 baseline uses zhighs commit `ceeee07`, Zig 0.16.0 ReleaseFast,
and HiGHS commit `de09bbad9f` built with GCC 13.3.0 using `-O3 -march=native
-DNDEBUG -flto`. It was captured under WSL2 on an Intel i9-10900KF. Absolute
times are diagnostic and are not correctness thresholds.

- `end_to_end_corpus.lock.tsv` pins all 40 model SHA-256 digests, expected
  statuses, and optimal objectives.
- `end_to_end_trace.lock.tsv` pins pivot-event counts and trace digests for
  `gas11`, `brandy`, `sc105`, and `scsd1`.
- `end_to_end_baseline.tsv` records status, objective, iterations, residuals,
  factorization lifecycle, and complete solve time for zhighs and HiGHS.

Run the complete correctness gate with:

```sh
bench/simplex/run_end_to_end_corpus.sh
```

The command validates corpus hashes and the HiGHS commit, builds both runners,
checks every status and optimal objective, enforces `1e-7` primal/dual residual
limits for optimal zhighs solutions, validates the `gas11` unbounded ray, and
checks the four pinned pivot traces. Any unexpected result returns a non-zero
exit status. Set `VERIFY_TRACES=0` only for an explicitly requested status-only
diagnostic run; it is not valid for the commit gate.

## Phase-attributed baseline

Stage 6.1 adds opt-in phase and kernel statistics. Normal solver calls do not
read clocks or scan vectors for density. The benchmark deliberately runs each
model twice per sample: an instrumentation-free run supplies fair total time,
while a separate profiling run supplies Phase I/II, INVERT, FTRAN, BTRAN,
PRICE, UPDATE, rebuild, density, requested-byte, and peak-RSS measurements.
Profiling stage times must not be added to or compared directly with the
instrumentation-free total.

Run two warmups and seven measured repetitions with:

```sh
bench/simplex/run_end_to_end_benchmark.sh
```

The complete median/p95 table is stored in `end_to_end_phase_baseline.tsv`.
Median total solve results from the frozen host are:

| model | zhighs | HiGHS | relative result | main zhighs cost |
|---|---:|---:|---:|---|
| gas11 | 20.76 ms | 47.15 ms | zhighs 2.27x faster | Phase I / pricing |
| sc105 | 0.90 ms | 5.65 ms | zhighs 6.24x faster | Phase II |
| brandy | 57.25 ms | 21.62 ms | zhighs 2.65x slower | Phase I / pricing |
| scsd1 | 13.00 ms | 6.46 ms | zhighs 2.01x slower | Phase II / rebuild / pricing |

Retained simplex requested bytes are 629741 (`gas11`), 394588 (`brandy`),
145752 (`sc105`), and 183485 (`scsd1`). The runner also reports process
`ru_maxrss`; unlike requested bytes this includes runtime and parser memory and
is platform-specific.
