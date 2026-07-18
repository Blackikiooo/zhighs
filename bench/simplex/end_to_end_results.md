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

After these fixes, 35 of 40 models match the available HiGHS reference.
`blend`, `grow7`, `scsd1`, and `vtp-base` remain numerical failures even with
reinversion after every pivot, which excludes FT update accumulation. `gas11`
is classified as long-running because the local HiGHS reference also exceeded
the 30-second diagnostic limit.
