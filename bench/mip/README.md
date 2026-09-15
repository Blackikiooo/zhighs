# MIP performance comparison

Build an optimized executable with `zig build build-bench-mip -Doptimize=ReleaseFast -j2`.
Before changing solver code, copy `zig-out/bin/mip-solver-bench.exe` to a separate
baseline path. Rebuild the candidate with the same Zig version, target and mode.

`compare_highs.py` requires Python and highspy. The recorded Windows run used
highspy 1.12.0; the 1.14.0 wheel crashed during extension import on this machine.
An isolated installation is sufficient:

```powershell
python -m pip install --target .zig-cache/highs-1.12 highspy==1.12.0
python bench/mip/compare_highs.py --baseline .zig-cache/perf-baseline/mip-solver-bench.exe --candidate zig-out/bin/mip-solver-bench.exe --highspy-path .zig-cache/highs-1.12 --output bench/mip/results/comparison.json --samples 5 --repetitions 30 test/model/data/warm_restart_mip.lp bench/mip/branching_compare.lp
```

Use `--repetitions 1000 --highs-repetitions 10` for microsecond-scale cases.
Run comparisons without other builds or benchmarks competing for the CPU.
The script rotates solver order and performs unmeasured warmup batches. Each
measured solve starts without retained solver state from the previous solve.
Model parsing and process startup are outside the timed region. Zhighs includes
solver initialization; HiGHS times `run()` after `clearSolver()`. Python call
overhead remains in the HiGHS samples.

HiGHS runs with one thread, parallel mode off, seed 0, relative gap 1e-4 and
absolute gap 1e-10. Both global presolve on and off are reported. Zhighs uses its
default MIP engine with root LP presolve, not the `Model.optimize` global MIP
presolve pipeline. These are comparisons of the supplied implementations, not
an assertion that algorithm settings or compiler flags are identical.

JSON records all samples, objective/bound checks, structural counters, model and
binary hashes, and reference version. A failed warmup is recorded and subsequent
models still run; any recorded failure makes the script exit nonzero. Runtime
status or objective mismatches are errors rather than performance samples.
The benchmark's LP counters currently exclude LP work inside primal heuristics;
use elapsed solve time when assessing diving improvements.

`rounding-only` and `no-heuristics` profiles in the Zig runner allow component
ablation. Keep `default` for the paired HiGHS comparison. Detailed results and
interpretation are recorded in `docs/mip-performance.md`.
