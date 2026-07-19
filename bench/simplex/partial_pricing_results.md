# Partial/multiple pricing A/B

## Environment and method

- zhighs commit under test: `16cf3a88fa20a929246e938bb27e3b58da687c4a`
- Zig: `0.16.0`, `ReleaseFast`
- CPU: Intel Core i9-10900KF, 10 cores / 20 threads, 20 MiB L3
- Corpus: SHA-256-pinned Stage 7 Netlib MPS files
- Policies: inherited Devex versus explicit `PRIMAL_PRICING_STRATEGY=partial`
- Five independent, sequential, statistics-enabled processes per model and policy
- Reported p95 is the maximum of five samples; no warm process state is shared

`pilot` runs to a verified optimal solution. Historical Stage 7 runs show that
`dfl001` and `fit2p` exceed 60 seconds, so their A/B uses the same deterministic
20,000-iteration budget. Their `iteration_limit` status is a measurement
boundary, not a solver conclusion.

## Results

| Model | Policy | Status | Iterations | Solve median / p95 | PRICE median | PRICE share | Pool searches / refills | Reuse hit rate | Requested bytes | Peak RSS median |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| dfl001 | inherit | iteration limit | 20,000 | 13.820 / 13.863 s | 4.837 s | 35.0% | 0 / 0 | 0% | 8,658,701 | 9,232 KiB |
| dfl001 | partial | iteration limit | 20,000 | 13.564 / 13.663 s | 4.615 s | 34.0% | 18,953 / 9,477 | 50.0% | 8,658,701 | 9,236 KiB |
| fit2p | inherit | iteration limit | 20,000 | 5.297 / 5.351 s | 1.497 s | 28.3% | 0 / 0 | 0% | 11,329,413 | 8,964 KiB |
| fit2p | partial | iteration limit | 20,000 | 5.715 / 5.751 s | 1.448 s | 25.3% | 256 / 128 | 50.0% | 11,329,413 | 8,824 KiB |
| pilot | inherit | optimal | 28,241 | 10.416 / 10.452 s | 2.432 s | 23.3% | 0 / 0 | 0% | 5,944,430 | 6,968 KiB |
| pilot | partial | optimal | 19,895 | 7.700 / 7.859 s | 1.718 s | 22.3% | 19,897 / 9,951 | 50.0% | 5,944,430 | 6,812 KiB |

The reuse hit rate is `(pool_searches - full_refills) / pool_searches`. The
strict one-use policy intentionally caps it near 50%: one validated cached
search must be followed by one global refill. This policy uses the existing
`flip_columns` workspace and therefore adds no requested memory.

## Interpretation

- `pilot` is a real algorithm-path win: iterations fall 29.6%, PRICE time
  falls 29.4%, and total solve time falls 26.1% while status, objective and
  residual validation remain unchanged.
- `dfl001` gains only 1.9% in total time at fixed work. PRICE falls 4.6%, but
  pricing still occupies 34% of the run and Phase I remains highly degenerate.
- `fit2p` rejects default enablement: PRICE falls 3.3%, but total time rises
  7.9%. Only 256 searches use the pool before perturbation dominates, while
  the changed pivot path increases factorization and solve work.
- An eight-use cache was tested and rejected: stale candidates led `scrs8` to
  a factorization failure after 5,728 pivots. The retained one-use policy
  passes the 40-model status/objective/residual/ray gate and changes `scrs8`
  from 1,716 to 1,617 iterations.

Partial/multiple pricing therefore remains an explicit forcing mode. A future
default needs a deterministic activation rule based on width, degeneracy epoch
and measured refill benefit; model-name dispatch is prohibited.
