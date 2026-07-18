#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
HIGHS_ROOT=${HIGHS_ROOT:-/home/godv/codefiles/cppfiles/HiGHS}
HIGHS_BUILD=${HIGHS_BUILD:-$HIGHS_ROOT/build-release-zhighs}
REPEATS=${REPEATS:-101}

cd "$ROOT"
ZIG_GLOBAL_CACHE_DIR=${ZIG_GLOBAL_CACHE_DIR:-/tmp/zhighs-zig-cache} \
  zig build build-bench-sparse-lu -Doptimize=ReleaseFast -Dcpu=native
clang++ -O3 -march=native -DNDEBUG -flto -std=c++17 \
  -I"$HIGHS_ROOT/highs" -I"$HIGHS_BUILD" \
  bench/simplex/highs_sparse_lu_bench.cpp \
  -L"$HIGHS_BUILD/lib" -lhighs -Wl,-rpath,"$HIGHS_BUILD/lib" \
  -o zig-out/bin/highs-sparse-lu-bench

for dimension in 128 256 512 1024; do
  # Alternate order by dimension to reduce systematic thermal/order bias.
  if (( dimension == 256 || dimension == 1024 )); then
    zig-out/bin/highs-sparse-lu-bench "$dimension" "$REPEATS"
    zig-out/bin/sparse-lu-bench "$dimension" "$REPEATS"
  else
    zig-out/bin/sparse-lu-bench "$dimension" "$REPEATS"
    zig-out/bin/highs-sparse-lu-bench "$dimension" "$REPEATS"
  fi
done
