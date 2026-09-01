#!/usr/bin/env bash
set -euo pipefail

readonly SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
readonly REPOSITORY_ROOT=$(cd -- "$SCRIPT_DIR/../.." && pwd)
readonly LOCK_FILE="$SCRIPT_DIR/end_to_end_corpus.lock.tsv"
CORPUS_DIR=${1:-/home/godv/documents/codefiles/cppfiles/scipoptsuite-10.0.1/soplex/check/instances}

mapfile -t MODELS < <(awk -F '\t' '!/^#/ && NF {print $1}' "$LOCK_FILE")
cd "$REPOSITORY_ROOT"
ZIG_GLOBAL_CACHE_DIR=${ZIG_GLOBAL_CACHE_DIR:-/tmp/zhighs-zig-cache} \
  zig build build-solver-session -Doptimize=ReleaseFast

result_file=$(mktemp /tmp/zhighs-solver-session.XXXXXX.tsv)
trap 'rm -f "$result_file"' EXIT
for model in "${MODELS[@]}"; do
  path="$CORPUS_DIR/$model.mps"
  [[ -f "$path" ]] || { echo "solver session corpus: missing $path" >&2; exit 1; }
  zig-out/bin/solver-session "$path" | tee -a "$result_file"
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
    result = path; sub(/^.*\//, "", result); sub(/\.mps$/, "", result); return result
  }
  function absolute(value) { return value < 0 ? -value : value }
  function reject(message) { print "solver session corpus: " message > "/dev/stderr"; failed = 1 }
  {
    model = basename($2)
    if ($3 != expected_status[model] || $4 != expected_status[model])
      reject(model " cold/warm status mismatch: " $3 "/" $4)
    if ($3 == "optimal") {
      tolerance = 1e-8 * (1 + absolute(expected_objective[model]))
      if (absolute($5 - expected_objective[model]) > tolerance || absolute($6 - expected_objective[model]) > tolerance)
        reject(model " cold/warm objective mismatch: " $5 "/" $6)
    }
    if ($9 != 1 || $10 != 1) reject(model " session counters are " $9 "/" $10)
    if ($8 > $7) reject(model " warm pivots " $8 " exceed cold pivots " $7)
  }
  END { exit failed ? 1 : 0 }
' "$result_file"

echo "solver session corpus: PASS (${#MODELS[@]} models)" >&2
