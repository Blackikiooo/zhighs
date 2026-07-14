#!/usr/bin/env bash
set -euo pipefail

DATASET_DIR="${MATRIX_DATASET_DIR:-}"
RUNNER="${MATRIX_DATASET_RUNNER:-}"
REPORT="${MATRIX_DATASET_REPORT:-/tmp/zhighs-matrix-datasets.tsv}"
if [ -z "$DATASET_DIR" ] || [ ! -d "$DATASET_DIR" ]; then
  echo "MATRIX_DATASET_DIR must contain the pinned Matrix Market corpus" >&2
  exit 1
fi
if [ -z "$RUNNER" ] || [ ! -x "$RUNNER" ]; then
  echo "MATRIX_DATASET_RUNNER must be an executable dataset benchmark/validator" >&2
  exit 1
fi

dataset_count="$(find "$DATASET_DIR" -type f -name '*.mtx' -size +1M | wc -l)"
if [ "$dataset_count" -lt 3 ]; then
  echo "dataset corpus must contain at least three .mtx files larger than 1 MiB" >&2
  exit 1
fi

"$RUNNER" "$DATASET_DIR" "$REPORT"
test -s "$REPORT"
awk -F '\t' '
  NR == 1 {
    if ($1 != "dataset" || $2 != "rows" || $3 != "cols" || $4 != "nnz" ||
        $5 != "elapsed_ms" || $6 != "peak_rss_kb" || $7 != "status") exit 2
    next
  }
  $7 != "PASS" || $2 < 10000 || $3 < 10000 || $4 < 100000 || $5 <= 0 || $6 <= 0 { exit 3 }
  { passed += 1 }
  END { if (passed < 3) exit 4 }
' "$REPORT"
echo "large dataset report accepted: $REPORT"
