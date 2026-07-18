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

After the parser and numerical fixes, 39 of 40 models match the available
HiGHS reference. The four former numerical failures now finish as follows:

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

`gas11` remains classified as long-running/unresolved because the local HiGHS
reference also exceeded the 30-second diagnostic limit.
