#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
HIGHS_ROOT=${HIGHS_ROOT:-/home/godv/codefiles/cppfiles/HiGHS}
HIGHS_BUILD=${HIGHS_BUILD:-$HIGHS_ROOT/build-release-zhighs}
SUITESPARSE_DIR=${SUITESPARSE_DIR:-/home/godv/codefiles/datasets/suitesparse/matrices}
NETLIB_DIR=${NETLIB_DIR:-/home/godv/codefiles/cppfiles/scipoptsuite-10.0.2/soplex/check/instances}
REPEATS=${REPEATS:-21}
LIMIT=${LIMIT:-128}
OUTPUT_DIR=${OUTPUT_DIR:-/tmp/zhighs-basis-corpus}

mkdir -p "$OUTPUT_DIR"
cd "$ROOT"
ZIG_GLOBAL_CACHE_DIR=${ZIG_GLOBAL_CACHE_DIR:-/tmp/zhighs-zig-cache} \
  zig build -Doptimize=ReleaseFast -Dcpu=native
clang++ -O3 -march=native -DNDEBUG -flto -std=c++17 \
  -I"$HIGHS_ROOT/highs" -I"$HIGHS_BUILD" \
  bench/simplex/highs_sparse_lu_bench.cpp \
  -L"$HIGHS_BUILD/lib" -lhighs -Wl,-rpath,"$HIGHS_BUILD/lib" \
  -o "$OUTPUT_DIR/highs-sparse-lu-bench"

for matrix in cage12 thermal1 webbase-1M; do
  zig build bench-suitesparse-basis -Doptimize=ReleaseFast -Dcpu=native -- \
    "$SUITESPARSE_DIR/$matrix.mtx" "$LIMIT"
done

for model in afiro sc50a sc105 adlittle brandy; do
  basis="$OUTPUT_DIR/$model.basis"
  "$OUTPUT_DIR/highs-sparse-lu-bench" --netlib "$NETLIB_DIR/$model.mps" "$REPEATS" "$basis"
  zig build bench-basis-file -Doptimize=ReleaseFast -Dcpu=native -- "$basis" "$REPEATS"
done
