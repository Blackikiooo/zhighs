#!/usr/bin/env bash
set -euo pipefail

readonly SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
readonly REPOSITORY_ROOT=$(cd -- "$SCRIPT_DIR/../.." && pwd)
CORPUS_DIR=${CORPUS_DIR:-/home/godv/codefiles/cppfiles/scipoptsuite-10.0.2/soplex/check/instances}
RUNS=${RUNS:-7}
WARMUPS=${WARMUPS:-2}
PHASE_ONE_STRATEGY=${PHASE_ONE_STRATEGY:-primal}
CRASH_STRATEGY=${CRASH_STRATEGY:-logical}
CRASH_MAX_COLUMNS=${CRASH_MAX_COLUMNS:-0}
DEGENERACY_STRATEGY=${DEGENERACY_STRATEGY:-auto}
PHASE_ONE_PRICING=${PHASE_ONE_PRICING:-inherit}
ADAPTIVE_REPRICE=${ADAPTIVE_REPRICE:-fixed}
PRICING_KERNEL=${PRICING_KERNEL:-column}
DEVEX_STRATEGY=${DEVEX_STRATEGY:-legacy}
MODELS=("$@")
if (($# == 0)); then MODELS=(gas11 brandy sc105 scsd1); fi

[[ "$RUNS" =~ ^[1-9][0-9]*$ ]] || { echo "RUNS must be positive" >&2; exit 2; }
[[ "$WARMUPS" =~ ^[0-9]+$ ]] || { echo "WARMUPS must be non-negative" >&2; exit 2; }

cd "$REPOSITORY_ROOT"
VERIFY_TRACES=0 "$SCRIPT_DIR/run_end_to_end_corpus.sh" "$CORPUS_DIR" "${MODELS[@]}" >/dev/null
readonly COMMON_ARGS=(1000000 100 no-trace 8 2 64)
readonly POLICY_ARGS=(32 "$PHASE_ONE_STRATEGY" "$CRASH_STRATEGY" "$CRASH_MAX_COLUMNS" \
  "$DEGENERACY_STRATEGY" "$PHASE_ONE_PRICING" "$ADAPTIVE_REPRICE" "$PRICING_KERNEL" "$DEVEX_STRATEGY")

for model in "${MODELS[@]}"; do
  for ((run = 0; run < WARMUPS; run++)); do
    zig-out/bin/simplex-end-to-end "$CORPUS_DIR/$model.mps" "${COMMON_ARGS[@]}" no-stats "${POLICY_ARGS[@]}" >/dev/null
    zig-out/bin/simplex-end-to-end "$CORPUS_DIR/$model.mps" "${COMMON_ARGS[@]}" stats "${POLICY_ARGS[@]}" >/dev/null
    zig-out/bin/highs-end-to-end "$CORPUS_DIR/$model.mps" >/dev/null
  done
done

samples=$(mktemp /tmp/zhighs-end-to-end-benchmark.XXXXXX.tsv)
trap 'rm -f "$samples"' EXIT
for model in "${MODELS[@]}"; do
  for ((run = 1; run <= RUNS; run++)); do
    zig_output=$(zig-out/bin/simplex-end-to-end "$CORPUS_DIR/$model.mps" "${COMMON_ARGS[@]}" no-stats "${POLICY_ARGS[@]}")
    stats_output=$(zig-out/bin/simplex-end-to-end "$CORPUS_DIR/$model.mps" "${COMMON_ARGS[@]}" stats "${POLICY_ARGS[@]}")
    highs_output=$(zig-out/bin/highs-end-to-end "$CORPUS_DIR/$model.mps")
    awk -F '\t' -v model="$model" -v run="$run" '
      $1 == "zhighs" && total == "" { total = $18 }
      $1 == "stats" {
        for (i = 3; i <= NF; i++) {
          split($i, pair, "=")
          value[pair[1]] = pair[2]
        }
      }
      END {
        OFS = "\t"
        print model, run, "zhighs_total", total
        print model, run, "phase1", value["phase1_ns"]
        print model, run, "phase2", value["phase2_ns"]
        print model, run, "rebuild", value["rebuild_ns"]
        print model, run, "invert", value["invert_ns"]
        print model, run, "ftran", value["ftran_ns"]
        print model, run, "btran", value["btran_ns"]
        print model, run, "price", value["pricing_ns"]
        print model, run, "update", value["update_ns"]
        print model, run, "requested_bytes", value["requested_bytes"]
        print model, run, "peak_rss_kb", value["peak_rss_kb"]
      }
    ' <<<"$zig_output
$stats_output" >>"$samples"
    awk -F '\t' -v model="$model" -v run="$run" '$1 == "highs" {OFS="\t"; print model, run, "highs_total", $9}' \
      <<<"$highs_output" >>"$samples"
  done
done

awk -F '\t' '
  function sort_values(values, count, i, j, temporary) {
    for (i = 2; i <= count; i++) {
      temporary = values[i]
      j = i - 1
      while (j >= 1 && values[j] > temporary) {
        values[j + 1] = values[j]
        j--
      }
      values[j + 1] = temporary
    }
  }
  {
    key = $1 SUBSEP $3
    count[key]++
    sample[key, count[key]] = $4 + 0
    model[$1] = 1
    metric[$3] = 1
  }
  END {
    OFS = "\t"
    for (model_name in model) for (metric_name in metric) {
      key = model_name SUBSEP metric_name
      n = count[key]
      if (n == 0) continue
      delete values
      for (i = 1; i <= n; i++) values[i] = sample[key, i]
      sort_values(values, n)
      median_index = int((n + 1) / 2)
      p95_index = int((95 * n + 99) / 100)
      if (p95_index > n) p95_index = n
      print model_name, metric_name, n, values[median_index], values[p95_index], values[1], values[n]
    }
  }
' "$samples" | {
  printf 'model\tmetric\truns\tmedian\tp95\tminimum\tmaximum\n'
  sort -t $'\t' -k1,1 -k2,2
}
