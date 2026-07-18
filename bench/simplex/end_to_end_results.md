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
| sc105 | numerical_failure / 116 | Optimal / 124 | HiGHS -52.202061211707225 | failure |
| brandy | numerical_failure / 31 | Optimal / 304 | HiGHS 1518.5098964881281 | failure |

`sc105` solves to `-52.20206121170723` in 114 iterations when the diagnostic
runner forces reinversion after every pivot. Its default failure therefore
isolates cumulative Forrest--Tomlin update error rather than parsing, model
semantics, or base SparseLU factorization. `brandy` still fails after 5636
iterations under the same forced-reinversion diagnostic, so it has an
additional long-running simplex stability or pivot-path defect.

Run the deterministic smoke corpus with:

```sh
bench/simplex/run_end_to_end_corpus.sh
```
