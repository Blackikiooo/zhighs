"""Paired simplex end-to-end comparison: zhighs vs HiGHS 1.14.0 (highspy).

Both sides solve from the same MPS file with presolve disabled and a single
thread. zhighs timing is the runner's internal solve_ns (parse excluded);
HiGHS timing wraps h.run() after readModel. Each side takes the median of
--repetitions samples. Status and objective agreement with the frozen corpus
lock is asserted before ratios are reported.

Usage:
  python bench/simplex/compare_e2e_highs.py [--repetitions 5] [--models a b c]
"""

import argparse
import json
import statistics
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CORPUS_DIR = Path(r"d:\CodesFiles\CppFiles\scipoptsuite-10.0.0\soplex\check\instances")
LOCK = Path(__file__).with_name("end_to_end_corpus.lock.tsv")
RUNNER = REPO / "zig-out" / "bin" / "simplex-end-to-end.exe"
HIGHSPY = REPO / ".zig-cache" / "bench-python"

RUNNER_ARGS = [
    "1000000", "100", "no-trace", "8", "2", "64", "no-stats", "32",
    "primal", "auto", "0", "auto", "inherit", "fixed", "auto", "framework",
    "inherit", "auto", "inherit", "auto",
]


def load_lock():
    rows = {}
    for line in LOCK.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        model, _sha, status, objective = line.split("\t")
        rows[model] = (status, float(objective))
    return rows


def zhighs_sample(model_path, repetitions):
    solve_ns, status, objective, iterations = [], None, None, None
    for _ in range(repetitions):
        out = subprocess.check_output(
            [str(RUNNER), str(model_path), *RUNNER_ARGS], text=True, timeout=600)
        fields = out.strip().split("\t")
        status, objective = fields[2], float(fields[4])
        iterations = int(fields[5])
        solve_ns.append(int(fields[17]))
    return statistics.median(solve_ns), status, objective, iterations


def highs_sample(highspy, model_path, repetitions):
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("threads", 1)
    h.setOptionValue("presolve", "off")
    h.setOptionValue("solver", "simplex")
    h.setOptionValue("random_seed", 0)
    if h.readModel(str(model_path)) != highspy.HighsStatus.kOk:
        raise RuntimeError(f"HiGHS cannot read {model_path}")
    elapsed = []
    status = objective = iterations = None
    for _ in range(repetitions):
        h.clearSolver()
        start = time.perf_counter_ns()
        h.run()
        elapsed.append(time.perf_counter_ns() - start)
        model_status = h.getModelStatus()
        status = str(model_status).replace("HighsModelStatus.k", "").lower()
        objective = h.getObjectiveValue() if model_status == highspy.HighsModelStatus.kOptimal else 0.0
        iterations = h.getInfo().simplex_iteration_count
    return statistics.median(elapsed), status, objective, iterations


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--models", nargs="*", default=None)
    args = parser.parse_args()

    import highspy

    lock = load_lock()
    models = args.models or sorted(lock)
    results = []
    for model in models:
        path = CORPUS_DIR / f"{model}.mps"
        z_ns, z_status, z_obj, z_iter = zhighs_sample(path, args.repetitions)
        h_ns, h_status, h_obj, h_iter = highs_sample(highspy, path, args.repetitions)
        expected_status, expected_obj = lock[model]
        ok = z_status == expected_status and (
            expected_status != "optimal" or abs(z_obj - expected_obj) <= 1e-8 * (1 + abs(expected_obj)))
        results.append({
            "model": model, "ok": ok,
            "z_ns": z_ns, "h_ns": h_ns, "ratio": z_ns / h_ns if h_ns else float("inf"),
            "z_status": z_status, "h_status": h_status,
            "z_iter": z_iter, "h_iter": h_iter,
            "iter_ratio": (z_iter / h_iter) if h_iter else None,
        })
        flag = "" if ok else "  <-- MISMATCH"
        print(f"{model:10s} z={z_ns/1e6:9.3f}ms h={h_ns/1e6:9.3f}ms "
              f"ratio={results[-1]['ratio']:6.3f} iter={z_iter:6d}/{h_iter:<6d}{flag}", flush=True)

    results.sort(key=lambda r: r["ratio"], reverse=True)
    print("\n=== sorted by z/h ratio ===")
    for r in results:
        print(f"{r['model']:10s} {r['ratio']:6.3f}")
    failed = [r["model"] for r in results if not r["ok"]]
    if failed:
        print("MISMATCH:", ", ".join(failed))
        sys.exit(1)
    out = Path(__file__).with_name("e2e_ratio_latest.json")
    out.write_text(json.dumps(results, indent=1))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()