# Presolve ablation benchmark

Measures presolve reductions and cold LP-relaxation solve time for five rule sets:

| mode | implied fixing | exact parallel rows | purpose |
|---|---|---|---|
| `off` | - | - | raw simplex baseline |
| `legacy` | no | no | pre-change rule set |
| `implied` | yes | no | isolate forcing/implied fixing |
| `parallel` | no | yes | isolate exact scaled rows |
| `full` | yes | yes | production defaults |

## Build and run

```powershell
zig build build-bench-presolve -Doptimize=ReleaseFast
python bench/presolve/generate_instances.py
python bench/presolve/compare.py --repetitions 10 --samples 3
```

`compare.py` writes `results/presolve_results.json` and
`results/presolve_results.md`. It also accepts arbitrary LP/MPS paths:

```powershell
python bench/presolve/compare.py --models path/to/model.mps path/to/model.lp
```

The binary itself can be used directly; it prints one TSV line per invocation:

```powershell
zig-out/bin/presolve-bench.exe bench/presolve/data/forcing_lp.lp full 10
```

## What is measured

- Parsing happens once and is outside the timed region.
- Each repetition builds a fresh engine and pipeline, so the solve is cold.
- Presolve, reduced solve, and postsolve are timed separately. `total_mean_ns`
  is their sum.
- The raw `off` objective is the reference when `off` reaches optimality. Every
  optimal mode must match it; recovered solutions are also checked against the
  original column bounds and row activities (`violation`).
- MIP models are sampled through their **LP relaxation**, with `nature` passed
  to the MIP-safe presolver so integer fixings are exercised. Full branch and
  bound performance belongs to `bench/mip`.

## Controlled instances

`generate_instances.py` writes deterministic LP files. Coefficients are small
integers and scale factors are `2`, `3`, `-1`, so exact parallel-row detection
is meaningful and reproducible.

- `forcing_lp.lp`: forced equalities plus a small core LP.
- `parallel_lp.lp`: base rows with exact scaled copies.
- `mip_forcing.lp`: binary forcing blocks plus knapsack rows.

The models are deliberately small enough that the raw simplex baseline reaches
optimality. Larger random sparse LPs can trigger raw-simplex numerical failures;
those are robustness observations, not presolve speedups.

## Interpretation

See `docs/presolve-strategy-review.md` for the literature review, correctness
argument, complete results and limitations. The short version: forcing and
parallel rows give large gains when the structure is present, but the added
propagation has a measurable fixed cost on small or already-clean models.
