#!/usr/bin/env bash
set -euo pipefail

readonly SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
readonly REPOSITORY_ROOT=$(cd -- "$SCRIPT_DIR/../.." && pwd)
readonly LOCK_FILE="$SCRIPT_DIR/end_to_end_corpus.lock.tsv"
CORPUS_DIR=${1:-/home/godv/documents/codefiles/cppfiles/scipoptsuite-10.0.1/soplex/check/instances}
if (($#)); then shift; fi

if (($#)); then
  MODELS=("$@")
else
  mapfile -t MODELS < <(awk -F '\t' '!/^#/ && NF {print $1}' "$LOCK_FILE")
fi

fail() {
  echo "solver pipeline corpus: $*" >&2
  exit 1
}

[[ -d "$CORPUS_DIR" ]] || fail "missing corpus directory $CORPUS_DIR"
for model in "${MODELS[@]}"; do
  path="$CORPUS_DIR/$model.mps"
  [[ -f "$path" ]] || fail "missing model $path"
  expected_hash=$(awk -F '\t' -v model="$model" '$1 == model {print $2}' "$LOCK_FILE")
  actual_hash=$(sha256sum "$path")
  actual_hash=${actual_hash%% *}
  [[ -n "$expected_hash" && "$actual_hash" == "$expected_hash" ]] || fail "SHA-256 mismatch for $model"
done

cd "$REPOSITORY_ROOT"
ZIG_GLOBAL_CACHE_DIR=${ZIG_GLOBAL_CACHE_DIR:-/tmp/zhighs-zig-cache} \
  zig build build-solver-pipeline -Doptimize=ReleaseFast

result_file=$(mktemp /tmp/zhighs-solver-pipeline.XXXXXX.tsv)
trap 'rm -f "$result_file"' EXIT
for model in "${MODELS[@]}"; do
  zig-out/bin/solver-pipeline "$CORPUS_DIR/$model.mps" 1000000 dual | tee -a "$result_file"
done

awk -F '\t' -v lock="$LOCK_FILE" '
  BEGIN {
    while ((getline line < lock) > 0) {
      if (line ~ /^#/ || line == "") continue
      split(line, field, "\t")
      expected_status[field[1]] = field[3]
      expected_objective[field[1]] = field[4]
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
    print "solver pipeline corpus: " message > "/dev/stderr"
    failed = 1
  }
  {
    model = basename($2)
    seen[model]++
    if (NF != 14) reject(model " returned an unavailable or malformed result")
    if ($3 != expected_status[model]) reject(model " status " $3 ", expected " expected_status[model])
    if ($3 == "optimal") {
      tolerance = 1e-8 * (1 + absolute(expected_objective[model]))
      if (absolute($4 - expected_objective[model]) > tolerance)
        reject(model " objective " $4 ", expected " expected_objective[model])
      if ($6 > 2e-7) reject(model " primal residual " $6)
      if ($7 > 1e-7) reject(model " dual residual " $7)
      if ($10 != "true" || $11 != "true" || $12 != "true")
        reject(model " optimal result validity flags are false")
    }
    if ($3 == "unbounded" && $8 > 1e-7) reject(model " ray residual " $8)
  }
  END {
    for (model in expected_status) if (seen[model] > 1) reject(model " appeared more than once")
    exit failed ? 1 : 0
  }
' "$result_file"

echo "solver pipeline corpus: PASS (${#MODELS[@]} models)" >&2
