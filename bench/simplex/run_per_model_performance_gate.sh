#!/usr/bin/env bash
set -euo pipefail

readonly SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
readonly REPOSITORY_ROOT=$(cd -- "$SCRIPT_DIR/../.." && pwd)
readonly MANIFEST="$SCRIPT_DIR/per_model_performance.tsv"
readonly CORPUS_LOCK="$SCRIPT_DIR/end_to_end_corpus.lock.tsv"

CORPUS_DIR=${CORPUS_DIR:-/home/godv/documents/codefiles/cppfiles/scipoptsuite-10.0.1/soplex/check/instances}
HIGHS_ROOT=${HIGHS_ROOT:-/home/godv/documents/codefiles/cppfiles/HiGHS}
HIGHS_BUILD=${HIGHS_BUILD:-$HIGHS_ROOT/cmake-build-release}
RUNS=${RUNS:-11}
WARMUPS=${WARMUPS:-2}
DUAL_EDGE_WEIGHT_STRATEGY=${DUAL_EDGE_WEIGHT_STRATEGY:-steepest-devex}
PRICING_KERNEL=${PRICING_KERNEL:-column}
DUAL_DSE_UPDATE_BUDGET=${DUAL_DSE_UPDATE_BUDGET:-64}

if ! [[ "$RUNS" =~ ^[1-9][0-9]*$ ]] || ((RUNS % 2 == 0)); then
  echo "RUNS must be a positive odd integer" >&2
  exit 2
fi
[[ "$WARMUPS" =~ ^[0-9]+$ ]] || { echo "WARMUPS must be non-negative" >&2; exit 2; }
[[ -f "$MANIFEST" && -f "$CORPUS_LOCK" ]] || { echo "missing simplex performance manifest" >&2; exit 2; }

if (($#)); then
  MODELS=("$@")
else
  mapfile -t MODELS < <(awk -F '\t' '!/^#/ && NF {print $1}' "$MANIFEST")
fi

for model in "${MODELS[@]}"; do
  awk -F '\t' -v model="$model" '$1 == model {found=1} END {exit !found}' "$MANIFEST" || {
    echo "model $model is not in $MANIFEST" >&2
    exit 2
  }
done

cd "$REPOSITORY_ROOT"
zig build build-simplex-end-to-end -Doptimize=ReleaseFast
c++ -O3 -march=native -DNDEBUG -flto \
  -I"$HIGHS_ROOT/highs" -I"$HIGHS_BUILD" \
  bench/simplex/highs_end_to_end_runner.cpp \
  -L"$HIGHS_BUILD/lib" -lhighs \
  -Wl,-rpath,"$HIGHS_BUILD/lib" -o zig-out/bin/highs-end-to-end

readonly ZIG_RUNNER="$REPOSITORY_ROOT/zig-out/bin/simplex-end-to-end"
readonly HIGHS_RUNNER="$REPOSITORY_ROOT/zig-out/bin/highs-end-to-end"
readonly SOLVE_ARGS=(1000000 100 no-trace 8 2 64 no-stats 32 dual logical 0 baseline inherit fixed "$PRICING_KERNEL" framework "$DUAL_EDGE_WEIGHT_STRATEGY" "$DUAL_DSE_UPDATE_BUDGET" inherit highs)
samples=$(mktemp /tmp/zhighs-per-model-performance.XXXXXX.tsv)
trap 'rm -f "$samples"' EXIT

for model in "${MODELS[@]}"; do
  path="$CORPUS_DIR/$model.mps"
  [[ -f "$path" ]] || { echo "missing model $path" >&2; exit 2; }
  expected_hash=$(awk -F '\t' -v model="$model" '$1 == model {print $2}' "$CORPUS_LOCK")
  actual_hash=$(sha256sum "$path")
  actual_hash=${actual_hash%% *}
  [[ -n "$expected_hash" && "$actual_hash" == "$expected_hash" ]] || {
    echo "hash mismatch for $model" >&2
    exit 2
  }

  for ((run = 0; run < WARMUPS; run++)); do
    "$ZIG_RUNNER" "$path" "${SOLVE_ARGS[@]}" >/dev/null
    "$HIGHS_RUNNER" "$path" >/dev/null
  done
  for ((run = 1; run <= RUNS; run++)); do
    if ((run % 2)); then
      zig_output=$("$ZIG_RUNNER" "$path" "${SOLVE_ARGS[@]}")
      highs_output=$("$HIGHS_RUNNER" "$path")
    else
      highs_output=$("$HIGHS_RUNNER" "$path")
      zig_output=$("$ZIG_RUNNER" "$path" "${SOLVE_ARGS[@]}")
    fi
    awk -F '\t' -v model="$model" -v run="$run" '
      $1 == "zhighs" {OFS="\t"; print model,run,"zhighs",$18,$6,$3,$5}
    ' <<<"$zig_output" >>"$samples"
    awk -F '\t' -v model="$model" -v run="$run" '
      $1 == "highs" {OFS="\t"; print model,run,"highs",$9,$5,tolower($3),$4}
    ' <<<"$highs_output" >>"$samples"
  done
done

awk -F '\t' -v manifest="$MANIFEST" '
  function sort_values(key, count, i, j, value) {
    for (i = 1; i <= count; i++) ordered[i] = sample[key, i]
    for (i = 2; i <= count; i++) {
      value = ordered[i]; j = i - 1
      while (j >= 1 && ordered[j] > value) { ordered[j + 1] = ordered[j]; j-- }
      ordered[j + 1] = value
    }
  }
  BEGIN {
    while ((getline line < manifest) > 0) {
      if (line ~ /^#/ || line == "") continue
      split(line, field, "\t")
      class[field[1]] = field[2]
      target[field[1]] = field[5] + 0
    }
    OFS="\t"
    print "model","class","runs","zhighs_median_ms","highs_median_ms","z_over_h","zhighs_iterations","highs_iterations","result"
  }
  {
    key = $1 SUBSEP $3
    count[key]++
    sample[key, count[key]] = $4 + 0
    iterations[key] = $5 + 0
    status[key] = $6
    objective[key] = $7 + 0
    models[$1] = 1
  }
  END {
    failed = 0
    for (model in models) {
      zig_key = model SUBSEP "zhighs"
      highs_key = model SUBSEP "highs"
      sort_values(zig_key, count[zig_key]); zig_median = ordered[int((count[zig_key] + 1) / 2)]; delete ordered
      sort_values(highs_key, count[highs_key]); highs_median = ordered[int((count[highs_key] + 1) / 2)]; delete ordered
      ratio = zig_median / highs_median
      result = ratio < target[model] ? "PASS" : "FAIL"
      if (status[zig_key] != status[highs_key]) result = "STATUS_MISMATCH"
      if (result != "PASS") failed = 1
      print model,class[model],count[zig_key],zig_median/1e6,highs_median/1e6,ratio,iterations[zig_key],iterations[highs_key],result
    }
    exit failed
  }
' "$samples" | sort -t $'\t' -k1,1
