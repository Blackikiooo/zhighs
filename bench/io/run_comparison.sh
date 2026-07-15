#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 3 ]]; then
  echo "usage: $0 MODEL [iterations=7] [warmups=2]" >&2
  exit 2
fi

MODEL=$1
ITERATIONS=${2:-7}
WARMUPS=${3:-2}
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
HIGHS_ROOT=${HIGHS_ROOT:-/home/godv/documents/codefiles/cppfiles/HiGHS}
HIGHS_BUILD=${HIGHS_BUILD:-$HIGHS_ROOT/cmake-build-release}
HIGHS_LIB=$HIGHS_BUILD/lib

if [[ ! -e "$HIGHS_LIB/libhighs.so" && ! -e "$HIGHS_LIB/libhighs.a" ]]; then
  echo "release HiGHS library not found under $HIGHS_LIB" >&2
  echo "configure it with: cmake -S $HIGHS_ROOT -B $HIGHS_BUILD -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON" >&2
  exit 1
fi

zig build -Doptimize=ReleaseFast build-io-bench
c++ -O3 -DNDEBUG -std=c++17 \
  -I"$HIGHS_ROOT/highs" -I"$HIGHS_BUILD" "$ROOT/bench/io/highs_parser_bench.cpp" \
  -L"$HIGHS_LIB" -lhighs -Wl,-rpath,"$HIGHS_LIB" \
  -o "$ROOT/zig-out/bin/highs-parser-bench"

echo $'implementation\tfile\trows\tcolumns\tnonzeros\tbest_ns\tmedian_ms\tMiB_per_s\tpeak_rss_kb\tchecksum'
"$ROOT/zig-out/bin/io-parser-bench" "$MODEL" "$ITERATIONS" "$WARMUPS"
"$ROOT/zig-out/bin/highs-parser-bench" "$MODEL" "$ITERATIONS" "$WARMUPS"
