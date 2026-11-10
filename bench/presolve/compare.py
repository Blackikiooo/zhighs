#!/usr/bin/env python3
"""Run the presolve ablation binary and write a compact JSON/Markdown report."""

import argparse
import json
from pathlib import Path
import statistics
import subprocess
import sys

MODELS = [
    "bench/presolve/data/forcing_lp.lp",
    "bench/presolve/data/parallel_lp.lp",
    "bench/presolve/data/mip_forcing.lp",
    "bench/mip/branching_compare.lp",
    "test/model/data/presolve_mip.mps",
    "test/model/data/warm_restart_mip.lp",
    "test/model/data/pipeline_lp.lp",
]
MODES = ["off", "legacy", "implied", "parallel", "full"]


def parse_line(text: str) -> dict:
    fields: dict[str, str] = {}
    for part in text.strip().split("\t")[1:]:
        if "=" in part:
            key, value = part.split("=", 1)
            fields[key] = value
    return fields


def run(binary: Path, model: str, mode: str, repetitions: int) -> dict:
    output = subprocess.check_output(
        [str(binary), model, mode, str(repetitions)],
        text=True,
        cwd=Path(__file__).resolve().parents[2],
        timeout=180,
    )
    for line in output.splitlines():
        if line.startswith("presolve-bench\t"):
            fields = parse_line(line)
            if "total_mean_ns" in fields:
                return fields
    raise RuntimeError(f"no result line: {output}")


def median_sample(binary: Path, model: str, mode: str, repetitions: int, samples: int) -> dict:
    rows = [run(binary, model, mode, repetitions) for _ in range(samples)]
    result = dict(rows[0])
    for key in ("presolve_mean_ns", "solve_mean_ns", "postsolve_mean_ns", "total_mean_ns"):
        result[key] = str(int(statistics.median(int(row[key]) for row in rows)))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=Path("zig-out/bin/presolve-bench.exe"))
    parser.add_argument("--output", type=Path, default=Path("bench/presolve/results/presolve_results.json"))
    parser.add_argument("--markdown", type=Path, default=Path("bench/presolve/results/presolve_results.md"))
    parser.add_argument("--models", nargs="*", default=MODELS,
                        help="LP/MPS paths; defaults to the controlled cases and repository models")
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--samples", type=int, default=3)
    args = parser.parse_args()
    if not args.binary.exists():
        raise SystemExit(f"binary not found: {args.binary}")

    report = {"binary": str(args.binary), "repetitions": args.repetitions,
              "samples": args.samples, "models": []}
    for model in args.models:
        if not Path(model).exists():
            continue
        entry = {"model": model, "runs": {}}
        for mode in MODES:
            try:
                entry["runs"][mode] = median_sample(args.binary, model, mode,
                                                    args.repetitions, args.samples)
            except (subprocess.CalledProcessError, RuntimeError) as error:
                entry["runs"][mode] = {"error": str(error)}
        report["models"].append(entry)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# Presolve ablation results",
        "",
        f"Binary: `{args.binary}`; repetitions per sample: {args.repetitions}; samples: {args.samples}.",
        "Times are means over the repetitions, then a median over the samples. `off` disables presolve;",
        "`legacy` keeps the pre-existing rules; `implied` and `parallel` enable one new rule each;",
        "`full` enables both. Objective must match `off` whenever the raw solve is optimal.",
        "",
        "| model | mode | rows | cols | nnz | implied fix | parallel rows | presolve us | solve us | total us | iters | status |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for entry in report["models"]:
        name = Path(entry["model"]).name
        legacy = entry["runs"].get("legacy", {})
        for mode in MODES:
            row = entry["runs"][mode]
            if "error" in row:
                lines.append(f"| {name} | {mode} | - | - | - | - | - | - | - | - | - | error |")
                continue
            total_us = int(row["total_mean_ns"]) / 1000.0
            solve_us = int(row["solve_mean_ns"]) / 1000.0
            presolve_us = int(row["presolve_mean_ns"]) / 1000.0
            lines.append(
                f"| {name} | {mode} | {row['rows']} | {row['cols']} | {row['nnz']} | "
                f"{row['implied_fixed']} | {row['parallel_rows']} | {presolve_us:.1f} | "
                f"{solve_us:.1f} | {total_us:.1f} | {row['iterations']} | {row['status']} |"
            )
            if mode == "full" and "error" not in legacy:
                legacy_total = int(legacy["total_mean_ns"])
                if legacy_total > 0:
                    lines.append(f"| {name} | speedup vs legacy | | | | | | | | "
                                 f"{legacy_total / int(row['total_mean_ns']):.2f}x | | |")
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {args.output} and {args.markdown}", file=sys.stderr)


if __name__ == "__main__":
    main()
