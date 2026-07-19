#!/usr/bin/env bash
set -euo pipefail

# Reproducible long-running acceptance harness for Stage 7. The script keeps
# every solver's native TSV payload and adds process outcome plus peak RSS, so
# incomplete runs remain distinguishable from solver-reported statuses.

readonly SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
readonly REPOSITORY_ROOT=$(cd -- "$SCRIPT_DIR/../.." && pwd)
readonly EXPECTED_HIGHS_COMMIT=de09bbad9fb7c5d39a1a464a7641bbb5531c6e9d

HIGHS_ROOT=${HIGHS_ROOT:-/home/godv/codefiles/cppfiles/HiGHS}
CORPUS_LOCK=${CORPUS_LOCK:-$SCRIPT_DIR/stage7_netlib.lock.tsv}
TIMEOUT_SECONDS=${TIMEOUT_SECONDS:-60}
MEMORY_LIMIT_KB=${MEMORY_LIMIT_KB:-0}
OUTPUT_FILE=${OUTPUT_FILE:-/tmp/zhighs-stage7-results.tsv}
BUILD_RUNNERS=${BUILD_RUNNERS:-1}
RUN_HIGHS=${RUN_HIGHS:-1}
CLP_RUNNER=${CLP_RUNNER:-}
DEGENERACY_STRATEGY=${DEGENERACY_STRATEGY:-auto}
ADAPTIVE_REPRICE=${ADAPTIVE_REPRICE:-fixed}

fail() {
  echo "stage7 corpus: $*" >&2
  exit 1
}

[[ $# -ge 1 ]] || fail "usage: $0 CORPUS_DIR [MODEL ...]"
CORPUS_DIR=$1
shift
[[ -d "$CORPUS_DIR" ]] || fail "missing corpus directory $CORPUS_DIR"
[[ "$TIMEOUT_SECONDS" =~ ^[1-9][0-9]*$ ]] || fail "TIMEOUT_SECONDS must be positive"
[[ "$MEMORY_LIMIT_KB" =~ ^[0-9]+$ ]] || fail "MEMORY_LIMIT_KB must be non-negative"
[[ -x /usr/bin/time ]] || fail "GNU /usr/bin/time is required for peak RSS"

if (($#)); then
  MODELS=("$@")
elif [[ -f "$CORPUS_LOCK" ]]; then
  mapfile -t MODELS < <(awk -F '\t' '!/^#/ && NF {print $1}' "$CORPUS_LOCK")
else
  mapfile -t MODELS < <(find "$CORPUS_DIR" -maxdepth 1 -type f -name '*.mps' -printf '%f\n' | sed 's/\.mps$//' | sort)
fi
((${#MODELS[@]} != 0)) || fail "no models selected"

if [[ -f "$CORPUS_LOCK" ]]; then
  for model in "${MODELS[@]}"; do
    path="$CORPUS_DIR/$model.mps"
    [[ -f "$path" ]] || fail "missing model $path"
    expected_hash=$(awk -F '\t' -v model="$model" '$1 == model {print $2}' "$CORPUS_LOCK")
    [[ -n "$expected_hash" ]] || fail "$model is not present in $CORPUS_LOCK"
    actual_hash=$(sha256sum "$path")
    actual_hash=${actual_hash%% *}
    [[ "$actual_hash" == "$expected_hash" ]] ||
      fail "SHA-256 mismatch for $model: expected $expected_hash, got $actual_hash"
  done
fi

cd "$REPOSITORY_ROOT"
if [[ "$BUILD_RUNNERS" == 1 ]]; then
  zig build build-simplex-end-to-end -Doptimize=ReleaseFast
fi
readonly ZHIGHS_RUNNER="$REPOSITORY_ROOT/zig-out/bin/simplex-end-to-end"
[[ -x "$ZHIGHS_RUNNER" ]] || fail "missing $ZHIGHS_RUNNER"

HIGHS_RUNNER="$REPOSITORY_ROOT/zig-out/bin/highs-end-to-end"
if [[ "$RUN_HIGHS" == 1 ]]; then
  [[ -d "$HIGHS_ROOT" ]] || fail "missing HiGHS checkout $HIGHS_ROOT"
  actual_highs_commit=$(git -C "$HIGHS_ROOT" rev-parse HEAD)
  if [[ ${ALLOW_HIGHS_VERSION_MISMATCH:-0} != 1 && "$actual_highs_commit" != "$EXPECTED_HIGHS_COMMIT" ]]; then
    fail "HiGHS commit $actual_highs_commit does not match pinned $EXPECTED_HIGHS_COMMIT"
  fi
  if [[ "$BUILD_RUNNERS" == 1 || ! -x "$HIGHS_RUNNER" ]]; then
    c++ -O3 -march=native -DNDEBUG -flto \
      -I"$HIGHS_ROOT/highs" -I"$HIGHS_ROOT/build" \
      bench/simplex/highs_end_to_end_runner.cpp \
      -L"$HIGHS_ROOT/build/lib" -lhighs \
      -Wl,-rpath,"$HIGHS_ROOT/build/lib" -o "$HIGHS_RUNNER"
  fi
fi

work_dir=$(mktemp -d /tmp/zhighs-stage7-run.XXXXXX)
trap 'rm -rf -- "$work_dir"' EXIT
: > "$OUTPUT_FILE"
{
  printf '# schema\trecord\tmodel\tsolver\tprocess_outcome\texit_code\tpeak_rss_kb\tnative_payload\n'
  printf '# zhighs_commit\t%s\n' "$(git rev-parse HEAD)"
  printf '# highs_commit\t%s\n' "${actual_highs_commit:-not-run}"
  printf '# timeout_seconds\t%s\n' "$TIMEOUT_SECONDS"
  printf '# memory_limit_kb\t%s\n' "$MEMORY_LIMIT_KB"
  printf '# degeneracy_strategy\t%s\n' "$DEGENERACY_STRATEGY"
  printf '# adaptive_reprice\t%s\n' "$ADAPTIVE_REPRICE"
  printf '# corpus_lock\t%s\n' "$CORPUS_LOCK"
} >> "$OUTPUT_FILE"

run_solver() {
  local model=$1
  local solver=$2
  shift 2
  local stdout_file="$work_dir/$model.$solver.stdout"
  local stderr_file="$work_dir/$model.$solver.stderr"
  local rss_file="$work_dir/$model.$solver.rss"
  local exit_code outcome rss

  set +e
  if ((MEMORY_LIMIT_KB == 0)); then
    /usr/bin/time -o "$rss_file" -f '%M' \
      timeout --signal=TERM --kill-after=5s "${TIMEOUT_SECONDS}s" "$@" \
      >"$stdout_file" 2>"$stderr_file"
  else
    (
      ulimit -v "$MEMORY_LIMIT_KB"
      /usr/bin/time -o "$rss_file" -f '%M' \
        timeout --signal=TERM --kill-after=5s "${TIMEOUT_SECONDS}s" "$@"
    ) >"$stdout_file" 2>"$stderr_file"
  fi
  exit_code=$?
  set -e

  rss=$(tail -1 "$rss_file" 2>/dev/null || true)
  [[ "$rss" =~ ^[0-9]+$ ]] || rss=0
  if ((exit_code == 0)); then
    outcome=completed
  elif ((exit_code == 124)); then
    outcome=timeout
  else
    outcome=error
  fi

  if [[ -s "$stdout_file" ]]; then
    while IFS= read -r line; do
      printf 'result\t%s\t%s\t%s\t%d\t%s\t%s\n' \
        "$model" "$solver" "$outcome" "$exit_code" "$rss" "$line" >> "$OUTPUT_FILE"
    done < "$stdout_file"
  else
    printf 'result\t%s\t%s\t%s\t%d\t%s\n' \
      "$model" "$solver" "$outcome" "$exit_code" "$rss" >> "$OUTPUT_FILE"
  fi
  if [[ -s "$stderr_file" ]]; then
    sed "s/^/stage7 corpus [$model $solver]: /" "$stderr_file" >&2
  fi
}

for model in "${MODELS[@]}"; do
  path="$CORPUS_DIR/$model.mps"
  echo "stage7 corpus: $model" >&2
  run_solver "$model" zhighs "$ZHIGHS_RUNNER" "$path" \
    1000000 100 no-trace 8 2 64 stats 32 primal logical 0 \
    "$DEGENERACY_STRATEGY" inherit "$ADAPTIVE_REPRICE" column
  if [[ "$RUN_HIGHS" == 1 ]]; then
    run_solver "$model" highs "$HIGHS_RUNNER" "$path"
  fi
  if [[ -n "$CLP_RUNNER" ]]; then
    [[ -x "$CLP_RUNNER" ]] || fail "CLP_RUNNER is not executable: $CLP_RUNNER"
    run_solver "$model" clp "$CLP_RUNNER" "$path"
  fi
done

echo "stage7 corpus: wrote $OUTPUT_FILE (${#MODELS[@]} models)" >&2
