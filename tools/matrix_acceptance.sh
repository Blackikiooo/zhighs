#!/usr/bin/env bash
set -uo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PROFILE="${1:-quick}"
if [ "$PROFILE" != "quick" ] && [ "$PROFILE" != "full" ]; then
  echo "usage: $0 [quick|full]" >&2
  exit 2
fi

STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
REPORT_DIR="${MATRIX_ACCEPTANCE_REPORT_DIR:-/tmp/zhighs-matrix-acceptance/$STAMP}"
GLOBAL_CACHE="${MATRIX_ACCEPTANCE_ZIG_CACHE:-/tmp/zhighs-zig-cache}"
mkdir -p "$REPORT_DIR/logs"
SUMMARY="$REPORT_DIR/summary.tsv"
printf 'gate\tstatus\tlog\n' > "$SUMMARY"
failures=0

run_gate() {
  local gate="$1"
  shift
  local log="$REPORT_DIR/logs/$gate.log"
  echo "=== $gate ==="
  if "$@" >"$log" 2>&1; then
    printf '%s\tPASS\t%s\n' "$gate" "$log" >> "$SUMMARY"
    echo "PASS: $gate"
  else
    printf '%s\tFAIL\t%s\n' "$gate" "$log" >> "$SUMMARY"
    echo "FAIL: $gate (see $log)"
    failures=$((failures + 1))
  fi
}

gate_highs_differential() {
  "$ROOT/tools/matrix_gate/highs_differential.sh"
}

gate_structural_safety() {
  cd "$ROOT"
  zig build test-matrix-acceptance --global-cache-dir "$GLOBAL_CACHE"
  zig build test-matrix-acceptance -Dhighs-int-width=w64 --global-cache-dir "$GLOBAL_CACHE"
  zig build test-model --global-cache-dir "$GLOBAL_CACHE"
  zig build test-model -Dhighs-int-width=w64 --global-cache-dir "$GLOBAL_CACHE"
}

gate_large_datasets() {
  MATRIX_DATASET_REPORT="$REPORT_DIR/datasets.tsv" \
    "$ROOT/tools/matrix_gate/large_datasets.sh"
}

gate_configuration_regression() {
  cd "$ROOT"
  local modes=(Debug)
  if [ "$PROFILE" = "full" ]; then modes=(Debug ReleaseSafe ReleaseFast); fi
  local mode width
  for mode in "${modes[@]}"; do
    for width in w32 w64; do
      zig build test-matrix -Doptimize="$mode" -Dhighs-int-width="$width" \
        --global-cache-dir "$GLOBAL_CACHE"
    done
  done
}

if [ "$PROFILE" = "full" ]; then
  run_gate highs_differential gate_highs_differential
else
  printf 'highs_differential\tNOT_RUN_QUICK\t-\n' >> "$SUMMARY"
fi
run_gate structural_fuzz_oom gate_structural_safety
if [ "$PROFILE" = "full" ]; then
  run_gate large_real_datasets gate_large_datasets
else
  printf 'large_real_datasets\tNOT_RUN_QUICK\t-\n' >> "$SUMMARY"
fi
run_gate configuration_regression gate_configuration_regression

echo "Report: $SUMMARY"
if [ "$failures" -ne 0 ]; then
  echo "MATRIX ACCEPTANCE: FAIL ($failures gate(s))"
  exit 1
fi
if [ "$PROFILE" = "quick" ]; then
  echo "MATRIX ACCEPTANCE: QUICK PASS (not a production-ready verdict)"
else
  echo "MATRIX ACCEPTANCE: FULL PASS"
fi
