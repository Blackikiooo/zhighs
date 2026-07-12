// C++ single-kernel profiling harness matching bench/matrix/perf_profile.zig.
// Build:
//   g++ -std=c++17 -O3 -march=native -DNDEBUG -flto \
//     -I<highs-build> -I<highs-build>/highs -I<highs-source>/highs \
//     bench/matrix/highs_perf_profile.cpp \
//     -L<highs-build>/lib -Wl,-rpath,<highs-build>/lib -lhighs \
//     -o highs-perf-profile
//
// Usage:
//   ZHIGHS_PERF_KERNEL=csc_to_csr_into ./highs-perf-profile

#define main highs_full_benchmark_main
#include "highs_matrix_bench.cpp"
#undef main

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <vector>

// kDimension, kNnz, checksum, makeCsc, cscAxSkippingZeros, cscAxSparseView,
// cscSparseAdd, Triplet, transposeCscInto, transposeCsc, buildFromSorted
// are all defined inside highs_matrix_bench.cpp (included above).

/// Structural hash covering starts + indices + values, so the comparison
/// harness can verify that structure-modifying kernels produce identical
/// layouts across implementations.  All integer fields are normalized to
/// uint64_t so the hash is independent of the native HighsInt width.
uint64_t structuralHash(const std::vector<HighsInt>& starts,
                        const std::vector<HighsInt>& indices,
                        const std::vector<double>& values) {
  // FNV-1a 64-bit
  uint64_t h = 0xcbf29ce484222325ULL;
  for (HighsInt s : starts) {
    uint64_t v = static_cast<uint64_t>(s);
    const auto* bytes = reinterpret_cast<const uint8_t*>(&v);
    for (size_t i = 0; i < sizeof(v); ++i) {
      h ^= bytes[i];
      h *= 0x100000001b3ULL;
    }
  }
  for (HighsInt idx : indices) {
    uint64_t v = static_cast<uint64_t>(idx);
    const auto* bytes = reinterpret_cast<const uint8_t*>(&v);
    for (size_t i = 0; i < sizeof(v); ++i) {
      h ^= bytes[i];
      h *= 0x100000001b3ULL;
    }
  }
  for (double v : values) {
    const auto* bytes = reinterpret_cast<const uint8_t*>(&v);
    for (size_t i = 0; i < sizeof(v); ++i) {
      h ^= bytes[i];
      h *= 0x100000001b3ULL;
    }
  }
  return h;
}

/// CSC → CSR via histogram + prefix sum + scatter (same algorithm as
/// zhighs fillCsrFromCscAssumeValid).  All output arrays must be pre-allocated.
void fillCsrFromCsc(const HighsSparseMatrix& csc,
                    std::vector<HighsInt>& row_starts,
                    std::vector<HighsInt>& col_indices,
                    std::vector<double>& values,
                    std::vector<HighsInt>& cursor) {
  // Clear row-starts (used as per-row counts during histogram)
  std::fill(row_starts.begin(), row_starts.begin() + csc.num_row_ + 1, 0);
  // Histogram
  for (HighsInt col = 0; col < csc.num_col_; ++col)
    for (HighsInt p = csc.start_[col]; p < csc.start_[col + 1]; ++p)
      ++row_starts[csc.index_[p] + 1];
  // Prefix sum → row starts
  for (HighsInt row = 0; row < csc.num_row_; ++row)
    row_starts[row + 1] += row_starts[row];
  // Copy cursor from row_starts
  std::copy(row_starts.begin(), row_starts.begin() + csc.num_row_, cursor.begin());
  // Scatter
  for (HighsInt col = 0; col < csc.num_col_; ++col) {
    for (HighsInt p = csc.start_[col]; p < csc.start_[col + 1]; ++p) {
      HighsInt dest = cursor[csc.index_[p]]++;
      col_indices[dest] = col;
      values[dest] = csc.value_[p];
    }
  }
}

int main() {
  const char* kernel = std::getenv("ZHIGHS_PERF_KERNEL");
  if (!kernel) {
    std::fprintf(stderr, "Usage: ZHIGHS_PERF_KERNEL=<name> %s\n", "highs-perf-profile");
    return 1;
  }

  // ── Build test matrix ──────────────────────────────────────────
  HighsSparseMatrix csc = makeCsc();
  HighsSparseMatrix csr = csc;
  csr.ensureRowwise();

  // Workspace buffers
  std::vector<double> dense_x(kDimension, 1.0);
  std::vector<double> sparse_x(kDimension, 0.0);
  std::vector<HighsInt> sparse_indices;
  std::vector<double> sparse_values;
  for (HighsInt i = 0; i < kDimension; i += 20) {
    sparse_x[i] = 1.0;
    sparse_indices.push_back(i);
    sparse_values.push_back(1.0);
  }
  std::vector<double> y(kDimension);

  // CSR fill buffers (pre-allocated)
  std::vector<HighsInt> csr_row_starts(kDimension + 1);
  std::vector<HighsInt> csr_col_indices(kNnz);
  std::vector<double> csr_values(kNnz);
  std::vector<HighsInt> csr_cursor(kDimension);

  // Transpose buffers (pre-allocated)
  HighsSparseMatrix transpose_target;
  transpose_target.format_ = MatrixFormat::kColwise;
  transpose_target.num_row_ = kDimension;
  transpose_target.num_col_ = kDimension;
  transpose_target.start_.resize(kDimension + 1);
  transpose_target.index_.resize(kNnz);
  transpose_target.value_.resize(kNnz);
  std::vector<HighsInt> transpose_cursor(kDimension);

  // Scale buffers
  HighsScale scale;
  scale.num_row = kDimension;
  scale.num_col = kDimension;
  scale.row.assign(kDimension, 1.0);
  scale.col.assign(kDimension, 1.0);

  // ── Select kernel ──────────────────────────────────────────────
  int repeats = 0;
  double result_checksum = 0.0;
  uint64_t result_struct_hash = 0;
  const auto start = Clock::now();

  if (!std::strcmp(kernel, "clear_output")) {
    repeats = 100000;
    for (int r = 0; r < repeats; ++r) {
      std::fill(y.begin(), y.end(), 0.0);
      asm volatile("" : : "g"(y.data()) : "memory");
    }
    result_checksum = checksum(y);
  } else if (!std::strcmp(kernel, "csc_ax_dense")) {
    repeats = 500;
    for (int r = 0; r < repeats; ++r) {
      csc.product(y, dense_x);
      asm volatile("" : : "g"(y.data()) : "memory");
    }
    result_checksum = checksum(y);
  } else if (!std::strcmp(kernel, "csc_ax_sparse_skip")) {
    repeats = 2000;
    for (int r = 0; r < repeats; ++r) {
      cscAxSkippingZeros(csc, sparse_x, y);
      asm volatile("" : : "g"(y.data()) : "memory");
    }
    result_checksum = checksum(y);
  } else if (!std::strcmp(kernel, "csc_ax_sparse_view")) {
    repeats = 4000;
    for (int r = 0; r < repeats; ++r) {
      cscAxSparseView(csc, sparse_indices, sparse_values, y);
      asm volatile("" : : "g"(y.data()) : "memory");
    }
    result_checksum = checksum(y);
  } else if (!std::strcmp(kernel, "csc_sparse_add_no_clear")) {
    repeats = 4000;
    std::fill(y.begin(), y.end(), 0.0);
    for (int r = 0; r < repeats; ++r) {
      cscSparseAdd(csc, sparse_indices, sparse_values, y);
      asm volatile("" : : "g"(y.data()) : "memory");
    }
    result_checksum = checksum(y);
  } else if (!std::strcmp(kernel, "csr_ax_dense")) {
    repeats = 500;
    for (int r = 0; r < repeats; ++r) {
      csr.product(y, dense_x);
      asm volatile("" : : "g"(y.data()) : "memory");
    }
    result_checksum = checksum(y);
  } else if (!std::strcmp(kernel, "csc_atx_dense")) {
    repeats = 500;
    for (int r = 0; r < repeats; ++r) {
      csc.productTranspose(y, dense_x);
      asm volatile("" : : "g"(y.data()) : "memory");
    }
    result_checksum = checksum(y);
  } else if (!std::strcmp(kernel, "csr_atx_dense")) {
    repeats = 500;
    for (int r = 0; r < repeats; ++r) {
      csr.productTranspose(y, dense_x);
      asm volatile("" : : "g"(y.data()) : "memory");
    }
    result_checksum = checksum(y);
  } else if (!std::strcmp(kernel, "alpha_ax_plus_y")) {
    repeats = 500;
    std::fill(y.begin(), y.end(), 0.0);
    for (int r = 0; r < repeats; ++r) {
      csc.alphaProductPlusY(1.0, dense_x, y, false);
      asm volatile("" : : "g"(y.data()) : "memory");
    }
    result_checksum = checksum(y);
  } else if (!std::strcmp(kernel, "product_quad")) {
    repeats = 500;
    for (int r = 0; r < repeats; ++r) {
      csc.productQuad(y, dense_x);
      asm volatile("" : : "g"(y.data()) : "memory");
    }
    result_checksum = checksum(y);
  } else if (!std::strcmp(kernel, "apply_scale")) {
    repeats = 500;
    for (int r = 0; r < repeats; ++r) {
      csc.applyScale(scale);
      asm volatile("" : : "g"(csc.value_.data()) : "memory");
    }
    result_checksum = checksum(csc.value_);
  } else if (!std::strcmp(kernel, "csc_to_csr_into")) {
    repeats = 100;
    for (int r = 0; r < repeats; ++r) {
      fillCsrFromCsc(csc, csr_row_starts, csr_col_indices, csr_values, csr_cursor);
      asm volatile("" : : "g"(csr_values.data()) : "memory");
    }
    result_checksum = checksum(csr_values);
    result_struct_hash = structuralHash(csr_row_starts, csr_col_indices, csr_values);
  } else if (!std::strcmp(kernel, "csc_to_csr_owning")) {
    repeats = 100;
    for (int r = 0; r < repeats; ++r) {
      HighsSparseMatrix fresh = csc;
      fresh.ensureRowwise();
      asm volatile("" : : "g"(fresh.value_.data()) : "memory");
    }
    result_checksum = checksum(csc.value_);
    // Build CSR via fillCsrFromCsc for struct hash verification (the owning path
    // produces CSR-ordered data; fillCsrFromCsc gives us the reference CSR layout).
    {
      std::vector<HighsInt> hash_starts(kDimension + 1);
      std::vector<HighsInt> hash_indices(kNnz);
      std::vector<double> hash_values(kNnz);
      std::vector<HighsInt> hash_cursor(kDimension);
      HighsSparseMatrix fresh = csc;
      fillCsrFromCsc(fresh, hash_starts, hash_indices, hash_values, hash_cursor);
      result_struct_hash = structuralHash(hash_starts, hash_indices, hash_values);
    }
  } else if (!std::strcmp(kernel, "transpose_into")) {
    repeats = 100;
    for (int r = 0; r < repeats; ++r) {
      transposeCscInto(csc, transpose_target, transpose_cursor);
      asm volatile("" : : "g"(transpose_target.value_.data()) : "memory");
    }
    result_checksum = checksum(transpose_target.value_);
    result_struct_hash = structuralHash(transpose_target.start_, transpose_target.index_, transpose_target.value_);
  } else if (!std::strcmp(kernel, "transpose_owning")) {
    repeats = 100;
    for (int r = 0; r < repeats; ++r) {
      HighsSparseMatrix t = transposeCsc(csc);
      asm volatile("" : : "g"(t.value_.data()) : "memory");
    }
    // Build one extra for checksum + structural hash
    {
      HighsSparseMatrix t = transposeCsc(csc);
      result_checksum = checksum(t.value_);
      result_struct_hash = structuralHash(t.start_, t.index_, t.value_);
    }
  } else if (!std::strcmp(kernel, "builder_freeze_sorted")) {
    repeats = 1000;
    for (int r = 0; r < repeats; ++r) {
      std::vector<Triplet> triplets;
      triplets.reserve(kNnz);
      size_t seq = 0;
      for (HighsInt col = 0; col < kDimension; ++col) {
        if (col != 0) triplets.push_back({col - 1, col, -1.0, seq++});
        triplets.push_back({col, col, 4.0, seq++});
        if (col + 1 < kDimension) triplets.push_back({col + 1, col, -1.0, seq++});
      }
      HighsSparseMatrix m = buildFromSorted(triplets);
      asm volatile("" : : "g"(m.value_.data()) : "memory");
    }
    // Build one extra for checksum
    {
      std::vector<Triplet> triplets;
      triplets.reserve(kNnz);
      size_t seq = 0;
      for (HighsInt col = 0; col < kDimension; ++col) {
        if (col != 0) triplets.push_back({col - 1, col, -1.0, seq++});
        triplets.push_back({col, col, 4.0, seq++});
        if (col + 1 < kDimension) triplets.push_back({col + 1, col, -1.0, seq++});
      }
      HighsSparseMatrix m = buildFromSorted(triplets);
      result_checksum = checksum(m.value_);
      result_struct_hash = structuralHash(m.start_, m.index_, m.value_);
    }
  } else if (!std::strcmp(kernel, "builder_freeze_prepopulated")) {
    repeats = 1000;
    // Pre-build sorted arrays once
    std::vector<Triplet> triplets;
    triplets.reserve(kNnz);
    size_t seq = 0;
    for (HighsInt col = 0; col < kDimension; ++col) {
      if (col != 0) triplets.push_back({col - 1, col, -1.0, seq++});
      triplets.push_back({col, col, 4.0, seq++});
      if (col + 1 < kDimension) triplets.push_back({col + 1, col, -1.0, seq++});
    }
    for (int r = 0; r < repeats; ++r) {
      HighsSparseMatrix m = buildFromSortedConst(triplets);
      asm volatile("" : : "g"(m.value_.data()) : "memory");
    }
    // Build one extra for checksum + struct hash
    {
      HighsSparseMatrix m = buildFromSortedConst(triplets);
      result_checksum = checksum(m.value_);
      result_struct_hash = structuralHash(m.start_, m.index_, m.value_);
    }
  } else if (!std::strcmp(kernel, "builder_freeze_canonical")) {
    repeats = 1000;
    // Pre-build canonical (no-duplicate) arrays
    std::vector<Triplet> ctriplets;
    ctriplets.reserve(kNnz);
    size_t cseq = 0;
    for (HighsInt col = 0; col < kDimension; ++col) {
      if (col != 0) ctriplets.push_back({col - 1, col, -1.0, cseq++});
      ctriplets.push_back({col, col, 4.0, cseq++});
      if (col + 1 < kDimension) ctriplets.push_back({col + 1, col, -1.0, cseq++});
    }
    for (int r = 0; r < repeats; ++r) {
      HighsSparseMatrix m = buildFromCanonical(ctriplets);
      asm volatile("" : : "g"(m.value_.data()) : "memory");
    }
    {
      HighsSparseMatrix m = buildFromCanonical(ctriplets);
      result_checksum = checksum(m.value_);
      result_struct_hash = structuralHash(m.start_, m.index_, m.value_);
    }
  } else if (!std::strcmp(kernel, "builder_freeze_reusable")) {
    repeats = 1000;
    // Pre-build canonical arrays once
    std::vector<Triplet> rtriplets;
    rtriplets.reserve(kNnz);
    size_t rseq = 0;
    for (HighsInt col = 0; col < kDimension; ++col) {
      if (col != 0) rtriplets.push_back({col - 1, col, -1.0, rseq++});
      rtriplets.push_back({col, col, 4.0, rseq++});
      if (col + 1 < kDimension) rtriplets.push_back({col + 1, col, -1.0, rseq++});
    }
    // Pre-allocate reusable buffers
    std::vector<HighsInt> r_starts(kDimension + 1);
    std::vector<HighsInt> r_indices(kNnz);
    std::vector<double> r_values(kNnz);
    // Pre-touch
    std::fill(r_starts.begin(), r_starts.end(), 0);
    std::fill(r_indices.begin(), r_indices.end(), 0);
    std::fill(r_values.begin(), r_values.end(), 0.0);
    for (int r = 0; r < repeats; ++r) {
      std::fill(r_starts.begin(), r_starts.end(), 0);
      for (const auto& t : rtriplets) ++r_starts[t.col + 1];
      for (HighsInt col = 0; col < kDimension; ++col)
        r_starts[col + 1] += r_starts[col];
      for (size_t i = 0; i < rtriplets.size(); ++i) {
        r_indices[i] = rtriplets[i].row;
        r_values[i] = rtriplets[i].value;
      }
      asm volatile("" : : "g"(r_values.data()) : "memory");
    }
    result_checksum = checksum(r_values);
    result_struct_hash = structuralHash(r_starts, r_indices, r_values);
  } else if (!std::strcmp(kernel, "builder_freeze_general")) {
    repeats = 100;
    for (int r = 0; r < repeats; ++r) {
      std::vector<Triplet> triplets;
      triplets.reserve(kNnz);
      size_t seq = 0;
      for (HighsInt col = kDimension; col-- > 0;) {
        if (col + 1 < kDimension) triplets.push_back({col + 1, col, -1.0, seq++});
        triplets.push_back({col, col, 4.0, seq++});
        if (col != 0) triplets.push_back({col - 1, col, -1.0, seq++});
      }
      std::sort(triplets.begin(), triplets.end(),
                [](const Triplet& a, const Triplet& b) {
                  if (a.col != b.col) return a.col < b.col;
                  if (a.row != b.row) return a.row < b.row;
                  return a.sequence < b.sequence;
                });
      HighsSparseMatrix m = buildFromSorted(triplets);
      asm volatile("" : : "g"(m.value_.data()) : "memory");
    }
    // Build one extra for checksum
    {
      std::vector<Triplet> triplets;
      triplets.reserve(kNnz);
      size_t seq = 0;
      for (HighsInt col = kDimension; col-- > 0;) {
        if (col + 1 < kDimension) triplets.push_back({col + 1, col, -1.0, seq++});
        triplets.push_back({col, col, 4.0, seq++});
        if (col != 0) triplets.push_back({col - 1, col, -1.0, seq++});
      }
      std::sort(triplets.begin(), triplets.end(),
                [](const Triplet& a, const Triplet& b) {
                  if (a.col != b.col) return a.col < b.col;
                  if (a.row != b.row) return a.row < b.row;
                  return a.sequence < b.sequence;
                });
      HighsSparseMatrix m = buildFromSorted(triplets);
      result_checksum = checksum(m.value_);
      result_struct_hash = structuralHash(m.start_, m.index_, m.value_);
    }
  } else if (!std::strcmp(kernel, "sparse_accumulate")) {
    repeats = 500;
    HighsSparseVectorSum accumulator(kDimension);
    for (int r = 0; r < repeats; ++r) {
      accumulator.clear();
      for (HighsInt i = 0; i < kDimension; ++i) {
        accumulator.add(i, 1.0);
        accumulator.add(i, -0.5);
      }
      asm volatile("" : : "g"(&accumulator) : "memory");
    }
    result_checksum = accumulator.getValue(kDimension / 2);
  } else {
    std::fprintf(stderr, "unknown ZHIGHS_PERF_KERNEL=%s\n", kernel);
    return 2;
  }

  const auto total =
      std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count();
  std::printf("cpp,%s,%lld,%.3f,%.17g,%llu\n", kernel, (long long)total,
              double(total) / repeats, result_checksum,
              (unsigned long long)result_struct_hash);
  return 0;
}
