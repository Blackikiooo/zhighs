#!/usr/bin/env bash
# Run an interleaved Zig-vs-C++ comparison for all matrix kernels.
# Outputs: $BUILD_ROOT/results/raw.csv  (all individual runs)
#          $BUILD_ROOT/results/summary.csv  (medians & derived stats)
#
# Usage:
#   CPU_CORE=2 RUNS=11 BUILD_ROOT=/tmp/zhighs-matrix-isolated \
#     HIGHS_SOURCE=/home/godv/documents/codefiles/cppfiles/HiGHS \
#     ./bench/matrix/run_isolated_comparison.sh
#
# Set FORCE_REBUILD=1 to rebuild both binaries even if they already exist.
# Set PERF_STAT=1 to collect 'perf stat' data for lagging kernels.

set -eo pipefail
# Intentionally NO set -u — some array expansions are empty.

# ── Configuration ─────────────────────────────────────────────────
CPU_CORE="${CPU_CORE:-2}"
RUNS="${RUNS:-11}"
BUILD_ROOT="${BUILD_ROOT:-/tmp/zhighs-matrix-isolated}"
HIGHS_SOURCE="${HIGHS_SOURCE:-/home/godv/documents/codefiles/cppfiles/HiGHS}"
FORCE_REBUILD="${FORCE_REBUILD:-1}"
PERF_STAT="${PERF_STAT:-0}"
PERF_BIN="${PERF_BIN:-perf}"
if [ "$PERF_BIN" = perf ]; then
  for candidate in /usr/lib/linux-tools-*/perf; do
    if [ -x "$candidate" ]; then PERF_BIN="$candidate"; break; fi
  done
fi
ZHIGHS_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
RESULTS_DIR="$BUILD_ROOT/results"

# ── Kernels (matching both perf_profile.zig and highs_perf_profile.cpp) ──
KERNELS=(
  clear_output
  csc_ax_dense
  csc_ax_sparse_skip
  csc_ax_sparse_view
  csc_sparse_add_no_clear
  csr_ax_dense
  csc_atx_dense
  csr_atx_dense
  alpha_ax_plus_y
  product_quad
  apply_scale
  csc_to_csr_into
  csc_to_csr_owning
  transpose_into
  transpose_owning
  builder_freeze_sorted
  builder_freeze_prepopulated
  builder_freeze_canonical
  builder_freeze_reusable
  builder_freeze_general
  sparse_accumulate
)

# Kernels where Zig is currently behind — collected with perf stat when PERF_STAT=1.
LAGGING_KERNELS=(
  clear_output
  csc_ax_dense
  csr_ax_dense
  csr_atx_dense
  transpose_into
  transpose_owning
  builder_freeze_sorted
)

# ── Helper: median from a list of numbers (newline-separated) ─────
median_of() {
  local arr=("$@")
  if [ ${#arr[@]} -eq 0 ]; then echo 0; return; fi
  local sorted
  sorted=($(printf '%s\n' "${arr[@]}" | sort -n))
  local n=${#sorted[@]}
  if (( n % 2 == 1 )); then
    echo "${sorted[$(( n / 2 ))]}"
  else
    awk -v a="${sorted[$n/2-1]}" -v b="${sorted[$n/2]}" \
      'BEGIN { printf "%.0f", (a + b) / 2 }'
  fi
}

# Median absolute deviation (MAD) as percentage of median
mad_pct() {
  local median="$1"
  shift
  local vals=("$@")
  if [ "${#vals[@]}" -eq 0 ] || [ "$(echo "$median == 0" | bc -l 2>/dev/null)" = "1" ]; then
    echo 0
    return
  fi
  local devs=()
  for v in "${vals[@]}"; do
    devs+=("$(echo "sqrt(($v - $median)^2)" | bc -l 2>/dev/null)")
  done
  local sorted_devs
  sorted_devs=($(printf '%s\n' "${devs[@]}" | sort -n))
  local n=${#sorted_devs[@]}
  local mad_median
  if (( n % 2 == 1 )); then
    mad_median="${sorted_devs[$(( n / 2 ))]}"
  else
    mad_median="$(awk -v a="${sorted_devs[$n/2-1]}" -v b="${sorted_devs[$n/2]}" \
      'BEGIN { printf "%.6f", (a + b) / 2 }')"
  fi
  echo "scale=2; 100 * $mad_median / $median" | bc -l 2>/dev/null
}

# ── Check environment ──────────────────────────────────────────────
echo "=== Environment check ==="
echo "CPU core: $CPU_CORE"
echo "Runs: $RUNS"
LOAD=$(cat /proc/loadavg | awk '{print $1}')
echo "Load average (1m): $LOAD"
if (( $(echo "$LOAD > 4.0" | bc -l 2>/dev/null) )); then
  echo "ERROR: Load average $LOAD is > 4.0. Refusing to produce noisy results."
  echo "Wait for system to quiet down and try again."
  exit 1
elif (( $(echo "$LOAD > 2.0" | bc -l 2>/dev/null) )); then
  echo "WARNING: Load average $LOAD is > 2.0. Results may be unreliable."
fi
lscpu | grep 'Model name\|CPU(s)\|Thread(s) per core' || true
echo "Zig: $(zig version)"
echo "G++: $(g++ --version | head -1)"
echo "HiGHS: $(git -C "$HIGHS_SOURCE" rev-parse HEAD 2>/dev/null || echo unknown)"
echo "Governor: $(cat /sys/devices/system/cpu/cpu$CPU_CORE/cpufreq/scaling_governor 2>/dev/null || echo unknown)"

# ── Build / locate binaries ───────────────────────────────────────
echo "=== Building ==="

# Always rebuild Zig perf-profile to pick up source changes
echo "Building Zig perf-profile..."
cd "$ZHIGHS_ROOT"
zig build build-perf-profile -Doptimize=ReleaseFast -Dcpu=native 2>&1

# Copy to BUILD_ROOT
mkdir -p "$BUILD_ROOT/bin"
cp zig-out/bin/perf-profile "$BUILD_ROOT/bin/perf-profile"

# Build or rebuild C++ perf-profile
if [ "$FORCE_REBUILD" = "1" ] || [ ! -x "$BUILD_ROOT/bin/highs-perf-profile" ]; then
  echo "Building HiGHS (if needed)..."
  if [ ! -f "$BUILD_ROOT/highs-build/lib/libhighs.so" ]; then
    cmake -S "$HIGHS_SOURCE" -B "$BUILD_ROOT/highs-build" \
      -DCMAKE_BUILD_TYPE=Release \
      '-DCMAKE_CXX_FLAGS_RELEASE=-O3 -march=native -DNDEBUG' \
      -DBUILD_TESTING=OFF -DFAST_BUILD=ON 2>&1
    cmake --build "$BUILD_ROOT/highs-build" --config Release -j"$(nproc)" 2>&1
  fi

  echo "Building C++ perf-profile..."
  rm -f "$BUILD_ROOT/bin/highs-perf-profile"
  g++ -std=c++17 -O3 -march=native -DNDEBUG -flto \
    -I"$BUILD_ROOT/highs-build" \
    -I"$BUILD_ROOT/highs-build/highs" \
    -I"$HIGHS_SOURCE/highs" \
    "$ZHIGHS_ROOT/bench/matrix/highs_perf_profile.cpp" \
    -L"$BUILD_ROOT/highs-build/lib" \
    -Wl,-rpath,"$BUILD_ROOT/highs-build/lib" -lhighs \
    -o "$BUILD_ROOT/bin/highs-perf-profile"
else
  echo "Using existing C++ perf-profile (set FORCE_REBUILD=1 to rebuild)."
fi

ZIG_BIN="$BUILD_ROOT/bin/perf-profile"
CPP_BIN="$BUILD_ROOT/bin/highs-perf-profile"
export LD_LIBRARY_PATH="$BUILD_ROOT/highs-build/lib:${LD_LIBRARY_PATH:-}"

# Verify both binaries work and record SHA-256
echo "Verifying binaries..."
ZHIGHS_PERF_KERNEL=clear_output ZHIGHS_PERF_ALLOCATOR=c "$ZIG_BIN" > /dev/null 2>&1 || {
  echo "ERROR: Zig binary failed to execute"; exit 1
}
ZHIGHS_PERF_KERNEL=clear_output "$CPP_BIN" > /dev/null 2>&1 || {
  echo "ERROR: C++ binary failed to execute"; exit 1
}
ZIG_SHA256=$(sha256sum "$ZIG_BIN" | awk '{print $1}')
CPP_SHA256=$(sha256sum "$CPP_BIN" | awk '{print $1}')
echo "Zig binary SHA-256:  $ZIG_SHA256"
echo "C++ binary SHA-256:  $CPP_SHA256"

# ── Run interleaved comparison ────────────────────────────────────
mkdir -p "$RESULTS_DIR"
RAW_CSV="$RESULTS_DIR/raw.csv"
# Write header with environment and binary hashes (SHA-256 now known)
cat > "$RAW_CSV" <<RAW_HEADER
# zhighs-matrix benchmark $(date -Iseconds)
# CPU core: $CPU_CORE  Runs: $RUNS  Load: $LOAD  Governor: $(cat /sys/devices/system/cpu/cpu$CPU_CORE/cpufreq/scaling_governor 2>/dev/null || echo unknown)
# Zig binary SHA-256:  $ZIG_SHA256
# C++ binary SHA-256:  $CPP_SHA256
# Owning allocator policy: Zig c_allocator, C++ malloc/std::vector
# Columns: run,kernel,impl,total_ns,ns_per_repeat,checksum,struct_hash
RAW_HEADER

for kernel in "${KERNELS[@]}"; do
  echo "=== Kernel: $kernel ==="
  declare -a zig_vals=()
  declare -a cpp_vals=()

  # Determine if we should run perf stat for this kernel
  RUN_PERF=0
  if [ "$PERF_STAT" = "1" ]; then
    for lag in "${LAGGING_KERNELS[@]}"; do
      if [ "$lag" = "$kernel" ]; then
        RUN_PERF=1
        break
      fi
    done
  fi

  for run in $(seq 1 "$RUNS"); do
    # Alternate starting order: odd runs Zig first, even runs C++ first
    if (( run % 2 == 1 )); then
      order=(zig cpp)
    else
      order=(cpp zig)
    fi

    for impl in "${order[@]}"; do
      if [ "$impl" = "zig" ]; then
        bin="$ZIG_BIN"
        impl_name="zig"
      else
        bin="$CPP_BIN"
        impl_name="cpp"
      fi

      # Run with taskset pinning — no || true, crashes must surface.
      raw_output="$(taskset -c "$CPU_CORE" env ZHIGHS_PERF_KERNEL="$kernel" ZHIGHS_PERF_ALLOCATOR=c "$bin" 2>&1)"
      _rc=$?
      if [ $_rc -ne 0 ]; then
        echo "ERROR: $impl_name/$kernel crashed with exit code $_rc at run $run"
        echo "  raw output was: $raw_output"
        exit 1
      fi
      output="$(echo "$raw_output" | grep '^zig,\|^cpp,' | head -1)"

      # Parse output: format is "impl,kernel,total_ns,ns_per_repeat,checksum,struct_hash"
      total_ns="$(echo "$output" | cut -d',' -f3)"
      ns_per_repeat="$(echo "$output" | cut -d',' -f4)"
      checksum_val="$(echo "$output" | cut -d',' -f5)"
      struct_hash_val="$(echo "$output" | cut -d',' -f6)"

      # Validate — abort on missing data
      if [ -z "$total_ns" ] || [ "$total_ns" = "" ]; then
        echo "ERROR: empty output for $impl_name/$kernel at run $run"
        echo "  raw output was: $raw_output"
        exit 1
      fi

      # Per-round checksum & struct-hash validation against first successful run
      cs_key="${kernel}__${impl_name}__checksum"
      sh_key="${kernel}__${impl_name}__struct_hash"
      if [ -z "${!cs_key}" ]; then
        eval "$cs_key=\"$checksum_val\""
        eval "$sh_key=\"$struct_hash_val\""
      else
        if [ "$checksum_val" != "${!cs_key}" ]; then
          echo "ERROR: $impl_name/$kernel checksum changed at run $run: expected ${!cs_key}, got $checksum_val"
          exit 1
        fi
        if [ "$struct_hash_val" != "${!sh_key}" ]; then
          echo "ERROR: $impl_name/$kernel struct hash changed at run $run: expected ${!sh_key}, got $struct_hash_val"
          exit 1
        fi
      fi

      echo "$run,$kernel,$impl_name,$total_ns,$ns_per_repeat,$checksum_val,$struct_hash_val" >> "$RAW_CSV"

      if [ "$impl_name" = "zig" ]; then
        zig_vals+=("$ns_per_repeat")
        _zig_cs_this_round="$checksum_val"
        _zig_sh_this_round="$struct_hash_val"
      else
        cpp_vals+=("$ns_per_repeat")
        _cpp_cs_this_round="$checksum_val"
        _cpp_sh_this_round="$struct_hash_val"
      fi
    done

    # Cross-implementation validation for this round: Zig vs C++ checksum & struct hash
    if [ -n "$_zig_cs_this_round" ] && [ -n "$_cpp_cs_this_round" ]; then
      # For non-structural kernels hash is 0 on both sides — that is expected and valid.
      # Compare checksum numerically via bc to handle formatting differences.
      _cs_diff=$(echo "$_zig_cs_this_round - $_cpp_cs_this_round" | bc -l 2>/dev/null || echo "1")
      if [ "$(echo "($_cs_diff < 0.0000000001) && ($_cs_diff > -0.0000000001)" | bc -l 2>/dev/null)" != "1" ]; then
        echo "ERROR: $kernel round $run cross-impl checksum mismatch: zig=$_zig_cs_this_round cpp=$_cpp_cs_this_round"
        exit 1
      fi
      if [ "$_zig_sh_this_round" != "$_cpp_sh_this_round" ]; then
        echo "ERROR: $kernel round $run cross-impl struct hash mismatch: zig=$_zig_sh_this_round cpp=$_cpp_sh_this_round"
        exit 1
      fi
    fi
  done

  # Compute stats — bail out if insufficient data
  if [ ${#zig_vals[@]} -lt "$RUNS" ] || [ ${#cpp_vals[@]} -lt "$RUNS" ]; then
    echo "ERROR: insufficient samples (zig: ${#zig_vals[@]}, cpp: ${#cpp_vals[@]}, expected: $RUNS)"
    exit 1
  fi

  zig_median=$(median_of "${zig_vals[@]}")
  zig_min=$(printf '%s\n' "${zig_vals[@]}" | sort -n | head -1)
  zig_max=$(printf '%s\n' "${zig_vals[@]}" | sort -n | tail -1)
  zig_mad=$(mad_pct "$zig_median" "${zig_vals[@]}")

  cpp_median=$(median_of "${cpp_vals[@]}")
  cpp_min=$(printf '%s\n' "${cpp_vals[@]}" | sort -n | head -1)
  cpp_max=$(printf '%s\n' "${cpp_vals[@]}" | sort -n | tail -1)
  cpp_mad=$(mad_pct "$cpp_median" "${cpp_vals[@]}")

  # Speedup: positive means Zig faster
  speedup=$(echo "scale=2; 100 * ($cpp_median - $zig_median) / $cpp_median" | bc -l 2>/dev/null)

  # Checksum and struct hash check
  last_zig_cs="$(grep ",$kernel,zig," "$RAW_CSV" 2>/dev/null | tail -1 | cut -d',' -f6)"
  last_cpp_cs="$(grep ",$kernel,cpp," "$RAW_CSV" 2>/dev/null | tail -1 | cut -d',' -f6)"
  last_zig_sh="$(grep ",$kernel,zig," "$RAW_CSV" 2>/dev/null | tail -1 | cut -d',' -f7)"
  last_cpp_sh="$(grep ",$kernel,cpp," "$RAW_CSV" 2>/dev/null | tail -1 | cut -d',' -f7)"
  cs_match="yes"
  sh_match="yes"
  if [ -n "$last_zig_cs" ] && [ -n "$last_cpp_cs" ]; then
    cs_diff="$(echo "$last_zig_cs - $last_cpp_cs" | bc -l 2>/dev/null || echo 1)"
    if [ "$(echo "($cs_diff < 0.0000000001) && ($cs_diff > -0.0000000001)" | bc -l 2>/dev/null)" != "1" ]; then
      cs_match="no (zig=$last_zig_cs cpp=$last_cpp_cs)"
    fi
  else
    cs_match="unknown"
  fi
  if [ -n "$last_zig_sh" ] && [ -n "$last_cpp_sh" ] && [ "$last_zig_sh" != "$last_cpp_sh" ]; then
    sh_match="no (zig=$last_zig_sh cpp=$last_cpp_sh)"
  fi

  echo "  Zig median: $zig_median ns  (min: $zig_min, max: $zig_max, MAD: $zig_mad%)"
  echo "  C++ median: $cpp_median ns  (min: $cpp_min, max: $cpp_max, MAD: $cpp_mad%)"
  echo "  Speedup: $speedup%"
  echo "  Checksum: $cs_match  StructHash: $sh_match"
done

# ── Optional: perf stat for lagging kernels ───────────────────────
if [ "$PERF_STAT" = "1" ]; then
  echo ""
  echo "=== Perf stat for lagging kernels ==="
  PERF_DIR="$RESULTS_DIR/perf"
  mkdir -p "$PERF_DIR"
  for kernel in "${LAGGING_KERNELS[@]}"; do
    echo "--- $kernel (zig) ---"
    # Single warmup run, then perf stat
    taskset -c "$CPU_CORE" env ZHIGHS_PERF_KERNEL="$kernel" ZHIGHS_PERF_ALLOCATOR=c "$ZIG_BIN" 2>&1 > /dev/null || true
    "$PERF_BIN" stat -e cycles,instructions,cache-references,cache-misses,L1-dcache-loads,L1-dcache-load-misses,LLC-loads,LLC-load-misses \
      taskset -c "$CPU_CORE" env ZHIGHS_PERF_KERNEL="$kernel" ZHIGHS_PERF_ALLOCATOR=c "$ZIG_BIN" \
      2>&1 | tee "$PERF_DIR/${kernel}_zig.perf"
    echo ""
    echo "--- $kernel (cpp) ---"
    taskset -c "$CPU_CORE" env ZHIGHS_PERF_KERNEL="$kernel" "$CPP_BIN" 2>&1 > /dev/null || true
    "$PERF_BIN" stat -e cycles,instructions,cache-references,cache-misses,L1-dcache-loads,L1-dcache-load-misses,LLC-loads,LLC-load-misses \
      taskset -c "$CPU_CORE" env ZHIGHS_PERF_KERNEL="$kernel" "$CPP_BIN" \
      2>&1 | tee "$PERF_DIR/${kernel}_cpp.perf"
    echo ""
  done
fi

# ── Produce summary.csv ───────────────────────────────────────────
SUMMARY_CSV="$RESULTS_DIR/summary.csv"
echo "kernel,zig_median_ns,cpp_median_ns,speedup_percent,zig_min_ns,cpp_min_ns,zig_mad_percent,cpp_mad_percent,checksum_match,struct_hash_match,zig_max_ns,cpp_max_ns" > "$SUMMARY_CSV"

for kernel in "${KERNELS[@]}"; do
  zig_ns=($(grep ",$kernel,zig," "$RAW_CSV" 2>/dev/null | cut -d',' -f5))
  cpp_ns=($(grep ",$kernel,cpp," "$RAW_CSV" 2>/dev/null | cut -d',' -f5))

  if [ ${#zig_ns[@]} -eq 0 ] || [ ${#cpp_ns[@]} -eq 0 ]; then
    echo "$kernel,,,,,,,,,," >> "$SUMMARY_CSV"
    continue
  fi

  zig_median=$(median_of "${zig_ns[@]}")
  cpp_median=$(median_of "${cpp_ns[@]}")
  zig_min=$(printf '%s\n' "${zig_ns[@]}" | sort -n | head -1)
  cpp_min=$(printf '%s\n' "${cpp_ns[@]}" | sort -n | head -1)
  zig_max=$(printf '%s\n' "${zig_ns[@]}" | sort -n | tail -1)
  cpp_max=$(printf '%s\n' "${cpp_ns[@]}" | sort -n | tail -1)
  zig_mad=$(mad_pct "$zig_median" "${zig_ns[@]}")
  cpp_mad=$(mad_pct "$cpp_median" "${cpp_ns[@]}")
  speedup=$(echo "scale=2; 100 * ($cpp_median - $zig_median) / $cpp_median" | bc -l 2>/dev/null)

  last_zig_cs="$(grep ",$kernel,zig," "$RAW_CSV" 2>/dev/null | tail -1 | cut -d',' -f6)"
  last_cpp_cs="$(grep ",$kernel,cpp," "$RAW_CSV" 2>/dev/null | tail -1 | cut -d',' -f6)"
  last_zig_sh="$(grep ",$kernel,zig," "$RAW_CSV" 2>/dev/null | tail -1 | cut -d',' -f7)"
  last_cpp_sh="$(grep ",$kernel,cpp," "$RAW_CSV" 2>/dev/null | tail -1 | cut -d',' -f7)"
  cs_match="yes"
  sh_match="yes"
  if [ -n "$last_zig_cs" ] && [ -n "$last_cpp_cs" ]; then
    cs_diff="$(echo "$last_zig_cs - $last_cpp_cs" | bc -l 2>/dev/null || echo 1)"
    if [ "$(echo "($cs_diff < 0.0000000001) && ($cs_diff > -0.0000000001)" | bc -l 2>/dev/null)" != "1" ]; then
      cs_match="no"
    fi
  else
    cs_match="unknown"
  fi
  if [ -n "$last_zig_sh" ] && [ -n "$last_cpp_sh" ] && [ "$last_zig_sh" != "$last_cpp_sh" ]; then
    sh_match="no"
  fi

  echo "$kernel,$zig_median,$cpp_median,$speedup,$zig_min,$cpp_min,$zig_mad,$cpp_mad,$cs_match,$sh_match,$zig_max,$cpp_max" >> "$SUMMARY_CSV"
done

echo ""
echo "=== Done ==="
echo "Raw CSV: $RAW_CSV"
echo "Summary CSV: $SUMMARY_CSV"
cat "$SUMMARY_CSV"
