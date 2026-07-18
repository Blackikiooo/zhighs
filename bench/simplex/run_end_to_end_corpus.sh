#!/usr/bin/env bash
set -euo pipefail

HIGHS_ROOT=${HIGHS_ROOT:-/home/godv/codefiles/cppfiles/HiGHS}
CORPUS_DIR=${1:-/home/godv/codefiles/cppfiles/scipoptsuite-10.0.2/soplex/check/instances}
if (($#)); then shift; fi
if (($#)); then MODELS=("$@"); else MODELS=(afiro adlittle sc50a sc105 brandy); fi

zig build build-simplex-end-to-end -Doptimize=ReleaseFast
c++ -O3 -march=native -DNDEBUG -flto \
  -I"$HIGHS_ROOT/highs" -I"$HIGHS_ROOT/build" \
  bench/simplex/highs_end_to_end_runner.cpp \
  -L"$HIGHS_ROOT/build/lib" -lhighs \
  -Wl,-rpath,"$HIGHS_ROOT/build/lib" -o zig-out/bin/highs-end-to-end

for model in "${MODELS[@]}"; do
  path="$CORPUS_DIR/$model.mps"
  zig-out/bin/simplex-end-to-end "$path"
  zig-out/bin/highs-end-to-end "$path"
done
