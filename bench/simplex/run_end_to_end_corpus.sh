#!/usr/bin/env bash
set -euo pipefail

readonly SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
readonly REPOSITORY_ROOT=$(cd -- "$SCRIPT_DIR/../.." && pwd)
readonly LOCK_FILE="$SCRIPT_DIR/end_to_end_corpus.lock.tsv"
readonly TRACE_LOCK_FILE="$SCRIPT_DIR/end_to_end_trace.lock.tsv"
readonly EXPECTED_HIGHS_COMMIT=de09bbad9fb7c5d39a1a464a7641bbb5531c6e9d

HIGHS_ROOT=${HIGHS_ROOT:-/home/godv/codefiles/cppfiles/HiGHS}
PHASE_ONE_STRATEGY=${PHASE_ONE_STRATEGY:-primal}
CRASH_STRATEGY=${CRASH_STRATEGY:-logical}
CRASH_MAX_COLUMNS=${CRASH_MAX_COLUMNS:-0}
DEGENERACY_STRATEGY=${DEGENERACY_STRATEGY:-baseline}
PHASE_ONE_PRICING=${PHASE_ONE_PRICING:-inherit}
ADAPTIVE_REPRICE=${ADAPTIVE_REPRICE:-fixed}
PRICING_KERNEL=${PRICING_KERNEL:-column}
DEVEX_STRATEGY=${DEVEX_STRATEGY:-legacy}
DUAL_EDGE_WEIGHT_STRATEGY=${DUAL_EDGE_WEIGHT_STRATEGY:-inherit}
DUAL_DSE_UPDATE_BUDGET=${DUAL_DSE_UPDATE_BUDGET:-64}
PRIMAL_PRICING_STRATEGY=${PRIMAL_PRICING_STRATEGY:-inherit}
CORPUS_DIR=${1:-/home/godv/codefiles/cppfiles/scipoptsuite-10.0.2/soplex/check/instances}
if (($#)); then shift; fi

if (($#)); then
  MODELS=("$@")
else
  mapfile -t MODELS < <(awk -F '\t' '!/^#/ && NF {print $1}' "$LOCK_FILE")
fi

fail() {
  echo "end-to-end corpus gate: $*" >&2
  exit 1
}

[[ -f "$LOCK_FILE" ]] || fail "missing $LOCK_FILE"
[[ -d "$CORPUS_DIR" ]] || fail "missing corpus directory $CORPUS_DIR"
[[ -d "$HIGHS_ROOT" ]] || fail "missing HiGHS checkout $HIGHS_ROOT"

if [[ ${ALLOW_HIGHS_VERSION_MISMATCH:-0} != 1 ]]; then
  actual_highs_commit=$(git -C "$HIGHS_ROOT" rev-parse HEAD)
  [[ "$actual_highs_commit" == "$EXPECTED_HIGHS_COMMIT" ]] ||
    fail "HiGHS commit $actual_highs_commit does not match pinned $EXPECTED_HIGHS_COMMIT"
fi

for model in "${MODELS[@]}"; do
  path="$CORPUS_DIR/$model.mps"
  [[ -f "$path" ]] || fail "missing model $path"
  expected_hash=$(awk -F '\t' -v model="$model" '$1 == model {print $2}' "$LOCK_FILE")
  [[ -n "$expected_hash" ]] || fail "model $model is not present in the corpus lock"
  actual_hash=$(sha256sum "$path")
  actual_hash=${actual_hash%% *}
  [[ "$actual_hash" == "$expected_hash" ]] ||
    fail "SHA-256 mismatch for $model: expected $expected_hash, got $actual_hash"
done

cd "$REPOSITORY_ROOT"
zig build build-simplex-end-to-end -Doptimize=ReleaseFast
c++ -O3 -march=native -DNDEBUG -flto \
  -I"$HIGHS_ROOT/highs" -I"$HIGHS_ROOT/build" \
  bench/simplex/highs_end_to_end_runner.cpp \
  -L"$HIGHS_ROOT/build/lib" -lhighs \
  -Wl,-rpath,"$HIGHS_ROOT/build/lib" -o zig-out/bin/highs-end-to-end

result_file=$(mktemp /tmp/zhighs-end-to-end-results.XXXXXX.tsv)
trap 'rm -f "$result_file"' EXIT
for model in "${MODELS[@]}"; do
  path="$CORPUS_DIR/$model.mps"
  zig-out/bin/simplex-end-to-end "$path" 1000000 100 no-trace 8 2 64 no-stats 32 \
    "$PHASE_ONE_STRATEGY" "$CRASH_STRATEGY" "$CRASH_MAX_COLUMNS" \
    "$DEGENERACY_STRATEGY" "$PHASE_ONE_PRICING" "$ADAPTIVE_REPRICE" "$PRICING_KERNEL" "$DEVEX_STRATEGY" \
    "$DUAL_EDGE_WEIGHT_STRATEGY" "$DUAL_DSE_UPDATE_BUDGET" "$PRIMAL_PRICING_STRATEGY" | tee -a "$result_file"
  zig-out/bin/highs-end-to-end "$path" | tee -a "$result_file"
done

# Timings and iteration counts are recorded for comparison but deliberately do
# not fail the correctness gate. Status, optimal objective, residuals and the
# unbounded certificate are deterministic acceptance conditions.
awk -F '\t' -v lock="$LOCK_FILE" '
  BEGIN {
    while ((getline line < lock) > 0) {
      if (line ~ /^#/ || line == "") continue
      split(line, field, "\t")
      expected_status[field[1]] = field[3]
      expected_objective[field[1]] = field[4]
      expected[field[1]] = 1
    }
  }
  function basename(path, result) {
    result = path
    sub(/^.*\//, "", result)
    sub(/\.mps$/, "", result)
    return result
  }
  function absolute(value) { return value < 0 ? -value : value }
  function reject(message) {
    print "end-to-end corpus gate: " message > "/dev/stderr"
    failed = 1
  }
  # Published residual upper bounds from the 93-model Stage 7 corpus:
  # max primal residual: 9.93e-8, max dual residual: 5.30e-8 (bounded
  # perturbation baseline). Assert at 2x the observed max as a regression
  # guard against residual drift.
  $1 == "zhighs" {
    model = basename($2)
    status = tolower($3)
    seen_zhighs[model]++
    if (!expected[model]) reject("unexpected zhighs model " model)
    if (status != expected_status[model])
      reject(model " zhighs status " status ", expected " expected_status[model])
    if ($4 != "none") reject(model " stopped at failure site " $4)
    if (status == "optimal") {
      tolerance = 1e-8 * (1 + absolute(expected_objective[model]))
      if (absolute($5 - expected_objective[model]) > tolerance)
        reject(model " objective " $5 ", expected " expected_objective[model])
      if ($7 > 2e-7) reject(model " primal residual " $7 " exceeds 2e-7 (baseline max 9.93e-8)")
      if ($8 > 1e-7) reject(model " dual residual " $8 " exceeds 1e-7 (baseline max 5.30e-8)")
    } else if (status == "unbounded") {
      if ($9 > 1e-7) reject(model " ray residual " $9 " exceeds 1e-7")
      if (model == "gas11" && $10 >= -1e-9)
        reject(model " ray objective " $10 " is not improving")
    }
  }
  $1 == "highs" {
    model = basename($2)
    status = tolower($3)
    seen_highs[model]++
    if (!expected[model]) reject("unexpected HiGHS model " model)
    if (status != expected_status[model])
      reject(model " HiGHS status " status ", expected " expected_status[model])
    if (status == "optimal") {
      tolerance = 1e-8 * (1 + absolute(expected_objective[model]))
      if (absolute($4 - expected_objective[model]) > tolerance)
        reject(model " HiGHS objective " $4 ", expected " expected_objective[model])
    }
  }
  END {
    for (model in expected) {
      selected = 0
      for (key in seen_zhighs) if (key == model) selected = 1
      for (key in seen_highs) if (key == model) selected = 1
      if (!selected) continue
      if (seen_zhighs[model] != 1) reject(model " has " seen_zhighs[model] " zhighs rows")
      if (seen_highs[model] != 1) reject(model " has " seen_highs[model] " HiGHS rows")
    }
    if (failed) exit 1
  }
' "$result_file"

if [[ ${VERIFY_TRACES:-1} == 1 ]]; then
  while IFS=$'\t' read -r model expected_count expected_digest; do
    [[ -n "$model" && ${model:0:1} != "#" ]] || continue
    selected=0
    for requested in "${MODELS[@]}"; do
      [[ "$requested" == "$model" ]] && selected=1
    done
    ((selected)) || continue

    trace_file=$(mktemp "/tmp/zhighs-$model-trace.XXXXXX.tsv")
    zig-out/bin/simplex-end-to-end "$CORPUS_DIR/$model.mps" 1000000 100 trace 8 2 64 no-stats 32 \
      "$PHASE_ONE_STRATEGY" "$CRASH_STRATEGY" "$CRASH_MAX_COLUMNS" \
      "$DEGENERACY_STRATEGY" "$PHASE_ONE_PRICING" "$ADAPTIVE_REPRICE" "$PRICING_KERNEL" "$DEVEX_STRATEGY" \
      "$DUAL_EDGE_WEIGHT_STRATEGY" "$DUAL_DSE_UPDATE_BUDGET" "$PRIMAL_PRICING_STRATEGY" \
      >/dev/null 2>"$trace_file"
    actual_count=$(awk -F '\t' '$1 == "pivot" {count++} END {print count + 0}' "$trace_file")
    # Hash the structural pivot path, not floating diagnostics that may vary
    # by target. Fields are phase, iteration, entering/leaving ids, row and
    # factor-update count.
    actual_digest=$(awk -F '\t' '$1 == "pivot" {print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $9}' "$trace_file" | sha256sum)
    actual_digest=${actual_digest%% *}
    rm -f "$trace_file"
    [[ "$actual_count" == "$expected_count" ]] ||
      fail "$model trace has $actual_count events, expected $expected_count"
    [[ "$actual_digest" == "$expected_digest" ]] ||
      fail "$model trace digest $actual_digest, expected $expected_digest"
  done < "$TRACE_LOCK_FILE"
fi

echo "end-to-end corpus gate: PASS (${#MODELS[@]} models)" >&2
