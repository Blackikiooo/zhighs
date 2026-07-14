#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
HIGHS_SOURCE="${HIGHS_SOURCE:-}"
if [ -z "$HIGHS_SOURCE" ] || [ ! -d "$HIGHS_SOURCE/highs" ]; then
  echo "HIGHS_SOURCE must point to a pinned HiGHS source checkout" >&2
  exit 1
fi
for command in cmake g++ bc taskset git; do
  command -v "$command" >/dev/null || {
    echo "required command missing: $command" >&2
    exit 1
  }
done

export HIGHS_SOURCE
export RUNS="${MATRIX_HIGHS_RUNS:-1}"
export PERF_STAT=0
export FORCE_REBUILD="${MATRIX_HIGHS_FORCE_REBUILD:-1}"
export BUILD_ROOT="${MATRIX_HIGHS_BUILD_ROOT:-/tmp/zhighs-matrix-gate-highs}"
"$ROOT/bench/matrix/run_isolated_comparison.sh"

summary="$BUILD_ROOT/results/summary.csv"
test -s "$summary"
if awk -F, 'NR > 1 && ($9 != "yes" || $10 != "yes") { exit 1 }' "$summary"; then
  echo "HiGHS differential checksums and structural hashes match"
else
  echo "HiGHS differential mismatch in $summary" >&2
  exit 1
fi
