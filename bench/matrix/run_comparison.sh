#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
HIGHS_SOURCE="${HIGHS_SOURCE:-/home/godv/documents/codefiles/cppfiles/HiGHS}"
BUILD_ROOT="${BUILD_ROOT:-/tmp/zhighs-matrix-comparison}"
CPU_CORE="${CPU_CORE:-2}"
RUNS="${RUNS:-7}"

cmake -S "$HIGHS_SOURCE" -B "$BUILD_ROOT/highs" \
  -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_CXX_FLAGS_RELEASE=-O3 -march=native -DNDEBUG' \
  -DBUILD_TESTING=OFF -DFAST_BUILD=ON
cmake --build "$BUILD_ROOT/highs" --config Release -j"$(nproc)"

zig build build-bench-matrix -Doptimize=ReleaseFast -Dcpu=native \
  --prefix "$BUILD_ROOT/zig-prefix"

g++ -std=c++17 -O3 -march=native -DNDEBUG -flto \
  -I"$BUILD_ROOT/highs" -I"$HIGHS_SOURCE/highs" \
  "$ROOT/bench/matrix/highs_matrix_bench.cpp" \
  -L"$BUILD_ROOT/highs/lib" -Wl,-rpath,"$BUILD_ROOT/highs/lib" -lhighs \
  -o "$BUILD_ROOT/highs-matrix-bench"

mkdir -p "$BUILD_ROOT/results"
for run in $(seq 1 "$RUNS"); do
  taskset -c "$CPU_CORE" "$BUILD_ROOT/zig-prefix/bin/matrix-bench" \
    2> "$BUILD_ROOT/results/zig-$run.csv"
  taskset -c "$CPU_CORE" "$BUILD_ROOT/highs-matrix-bench" \
    > "$BUILD_ROOT/results/cpp-$run.csv"
done

echo "Raw CSV files: $BUILD_ROOT/results"
