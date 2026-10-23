# IPM 40-model real-corpus validation

This report records the final run in `real_corpus_results_final.json`. The
corpus is the complete 40-model SoPlex `check/instances` continuous-LP set:
8x8 through 536x1775. Every input SHA-256 matches
`bench/simplex/end_to_end_corpus.lock.tsv`.

The candidate is the ReleaseFast sparse-Cholesky IPM. HiGHS 1.12.0 runs with
one thread and deterministic settings. HiGHS simplex with presolve supplies
the independent status/objective oracle; HiGHS IPM without presolve or
crossover supplies the timing reference. Timings are one sample per model and
are diagnostic rather than a stable benchmark claim.

Run it with:

```powershell
zig build build-ipm-corpus -Doptimize=ReleaseFast
python -I -S bench/ipm/validate_corpus.py `
  --corpus D:/CodesFiles/CppFiles/scipoptsuite-10.0.0/soplex/check/instances `
  --runner zig-out/bin/ipm-corpus.exe `
  --highspy-path .zig-cache/highs-1.12 `
  --highspy-path .zig-cache/bench-python `
  --output bench/ipm/real_corpus_results_final.json `
  --refinement-steps 0 --max-iterations 200 `
  --feasibility-tolerance 1e-7 --primal-correction-steps 5 `
  --max-regularization 1e-6
python -I -S bench/ipm/analyze_results.py bench/ipm/real_corpus_results_final.json
```

## Correctness result

| Metric | Current | Previous | Initial |
|---|---:|---:|---:|
| Independently KKT-validated models | 29/40 | 27/40 | 18/40 |
| Validated among reference-optimal models | 29/30 | 27/30 | 18/30 |
| Candidate status `optimal` | 29/40 | 27/40 | 27/40 |
| Hash-matched inputs | 40/40 | 40/40 | 40/40 |

All 29 published optima pass checks performed on the matrix independently
parsed by HiGHS: dimensions, finite vectors, original row and column bounds,
dual signs, stationarity, primal-dual gap, recomputed objective, and agreement
with the HiGHS simplex objective. Across accepted models, the worst relative
primal residual is 9.31e-8, stationarity is 4.24e-16, primal-dual gap is
9.96e-8, and reference-objective error is 7.41e-8.

The two new accepted models are `agg` and `vtp-base`. A failed unscaled solve
now retries once with max-norm equilibration within the original time budget;
models whose matrix coefficient range reaches 1e7 start with the scaled path.
The scaled convergence test evaluates primal residuals in unscaled row units;
the old denominator could report convergence while original-coordinate row
violations remained around 1e-6. An active-set projection on the scaled retry
removes correction variables that would cross zero.

`capri` is the only remaining reference-optimal failure. It reaches a late
normal-equation factorization failure after 17 iterations. Dense factorization,
equilibration, two refinement steps, 500 iterations, and regularization up to
1.0 were also tested and did not recover it.

The other ten failures match the non-optimal reference set: nine infeasible
models and the unbounded `gas11`. The current infeasible-start formulation has
no homogeneous self-dual embedding or certificate path, so it does not yet
publish those statuses. No false optimal result occurred.

## Performance result

The normal-equation backend now uses exact minimum-degree row ordering, with
RCM selected for row graphs above 60% density, and a row-oriented sparse
Cholesky pattern. Normal assembly positions and sparse
dot-product intersections are precomputed; intersection pairs use packed
32-bit indices. This replaced the skyline envelope that dominated the earlier
run. Examples:

| Model | Previous stored entries | Current stored entries |
|---|---:|---:|
| seba | 364,599 | 55,449 |
| etamacro | 120,196 | 12,208 |
| finnis | 113,064 | 5,973 |
| shell | 82,704 | 4,583 |
| grow7 | 63,035 | 5,293 |

Among the 28 models where both candidate and reference IPM report optimal,
the candidate/HiGHS-IPM solve-time geometric mean is 0.29x. The candidate is
faster on 25/28, compared with 19/26 previously. The best ratio is 0.060x.
Factorization now accounts for 59.7% of successful solve time, down from 94.3%.

The requested comprehensive performance win is not yet achieved. The three
measured counterexamples are `seba` (12.54x), `agg` (2.70x), and `israel`
(1.44x). `seba` performs 5,233,661 symbolic factor products per numeric
factorization; it needs a supernodal/block sparse kernel or a materially better
ordering. `agg` now starts scaled but its factorization remains slower than the
reference on this model.
These are the next concrete performance targets; the current evidence supports
a broad 25/28 win and a 3.45x geometric-mean speedup, not a universal win.
