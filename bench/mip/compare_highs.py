"""Paired cold-solve comparison; parsing and process startup are not timed.

Install highspy into an isolated directory and pass --highspy-path if needed.
Results retain every sample, model hash, solver counters and environment.
"""

import argparse
import hashlib
import json
import math
import os
import platform
from pathlib import Path
import statistics
import subprocess
import sys
import time


def zig_sample(executable, model, repetitions):
    output = subprocess.check_output(
        [str(executable), str(model), str(repetitions), "default"],
        text=True, timeout=180,
    )
    fields = dict(part.split("=", 1) for part in output.strip().split("\t")[2:])
    if fields["status"] != "optimal":
        raise RuntimeError(f"{executable}: {model}: {fields}")
    return {"ns": int(fields["mip_mean_ns"]), "objective": float(fields["objective"]),
            "bound": float(fields["bound"]), "counters": fields}


def highs_sample(highspy, model, repetitions, presolve):
    h = highspy.Highs()
    for name, value in {"output_flag": False, "threads": 1, "parallel": "off",
                        "random_seed": 0, "mip_rel_gap": 1e-4,
                        "mip_abs_gap": 1e-10, "presolve": presolve,
                        "time_limit": 120.0}.items():
        if h.setOptionValue(name, value) != highspy.HighsStatus.kOk:
            raise RuntimeError(f"HiGHS rejected option {name}")
    if h.readModel(str(model)) != highspy.HighsStatus.kOk:
        raise RuntimeError(f"HiGHS cannot read {model}")
    elapsed = 0
    for _ in range(repetitions):
        h.clearSolver()  # Discard the previous solution and basis outside timing.
        start = time.perf_counter_ns()
        h.run()
        elapsed += time.perf_counter_ns() - start
        if h.getModelStatus() != highspy.HighsModelStatus.kOptimal:
            raise RuntimeError(f"HiGHS: {model}: {h.getModelStatus()}")
    info = h.getInfo()
    objective = h.getObjectiveValue()
    return {"ns": elapsed // repetitions, "objective": objective,
            "bound": info.mip_dual_bound if h.getLp().integrality_ else objective,
            "nodes": info.mip_node_count, "iterations": info.simplex_iteration_count}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--highspy-path", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=5)
    parser.add_argument("--repetitions", type=int, default=20)
    parser.add_argument("--highs-repetitions", type=int, help="separate batch size for the reference solver")
    parser.add_argument("models", type=Path, nargs="+")
    args = parser.parse_args()
    if min(args.samples, args.repetitions) < 1:
        parser.error("sample and repetition counts must be positive")
    highs_repetitions = args.highs_repetitions if args.highs_repetitions is not None else args.repetitions
    if highs_repetitions < 1:
        parser.error("HiGHS repetition count must be positive")
    if args.highspy_path:
        sys.path.insert(0, str(args.highspy_path.resolve()))
    import highspy

    report = {"platform": platform.platform(), "cpu": platform.processor() or os.environ.get("PROCESSOR_IDENTIFIER", "unknown"),
              "python": sys.version, "highs": highspy.Highs().version(),
              "samples": args.samples, "repetitions": args.repetitions,
              "highs_repetitions": highs_repetitions,
              "zig_profile": "default", "highs_threads": 1,
              "highs_parallel": "off", "highs_random_seed": 0,
              "mip_relative_gap": 1e-4, "mip_absolute_gap": 1e-10,
              "binary_sha256": {name: hashlib.sha256(path.read_bytes()).hexdigest()
                                for name, path in (("baseline", args.baseline), ("candidate", args.candidate))},
              "models": []}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for model in args.models:
        runners = {
            "baseline": lambda: zig_sample(args.baseline.resolve(), model, args.repetitions),
            "candidate": lambda: zig_sample(args.candidate.resolve(), model, args.repetitions),
            "highs_on": lambda: highs_sample(highspy, model, highs_repetitions, "on"),
            "highs_off": lambda: highs_sample(highspy, model, highs_repetitions, "off"),
        }
        # Warm the executable pages/caches, but every measured solve stays cold.
        try:
            for run in runners.values():
                run()
        except (RuntimeError, subprocess.SubprocessError) as error:
            report["models"].append({"path": str(model),
                "sha256": hashlib.sha256(model.read_bytes()).hexdigest(), "error": str(error)})
            args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
            print(f"{model.name}: FAILED: {error}", flush=True)
            continue
        samples = {name: [] for name in runners}
        names = list(runners)
        for sample in range(args.samples):
            order = names[sample % len(names):] + names[:sample % len(names)]
            for name in order:
                samples[name].append(runners[name]())
        expected = samples["highs_on"][0]["objective"]
        for name, rows in samples.items():
            for row in rows:
                if not math.isclose(row["objective"], expected, rel_tol=1e-7, abs_tol=1e-6):
                    raise RuntimeError(f"objective mismatch: {model}: {name}: {row}")
                if not math.isclose(row["bound"], expected, rel_tol=1e-4, abs_tol=1e-6):
                    raise RuntimeError(f"bound mismatch: {model}: {name}: {row}")
        medians = {name: statistics.median(row["ns"] for row in rows) for name, rows in samples.items()}
        entry = {"path": str(model), "sha256": hashlib.sha256(model.read_bytes()).hexdigest(),
                 "samples": samples, "median_ns": medians,
                 "speedup": medians["baseline"] / medians["candidate"]}
        report["models"].append(entry)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(f"{model.name}: speedup={entry['speedup']:.3f} medians_ns={medians}", flush=True)
    if any("error" in entry for entry in report["models"]):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
