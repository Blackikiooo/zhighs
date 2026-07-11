#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "lp_data/HStruct.h"
#include "util/HighsSparseMatrix.h"

namespace {
constexpr HighsInt kDimension = 50000;
constexpr HighsInt kNnz = 3 * kDimension - 2;
constexpr int kProductRepeats = 200;
constexpr int kQuadRepeats = 20;
constexpr int kTransformRepeats = 10;
constexpr int kAccumulatorRepeats = 20;

using Clock = std::chrono::steady_clock;

void report(const char* implementation, const char* kernel, int repeats,
            Clock::time_point start, double checksum) {
  auto total = std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count();
  std::printf("%s,%s,%d,%d,%d,%lld,%.3f,%.17g\n", implementation, kernel,
              int(kDimension), int(kNnz), repeats, (long long)total,
              double(total) / repeats, checksum);
}

inline void clobber(std::vector<double>& values) {
  asm volatile("" : : "g"(values.data()) : "memory");
}

double checksum(const std::vector<double>& values) {
  double sum = 0.0;
  for (double value : values) sum += value;
  return sum;
}

HighsSparseMatrix makeCsc() {
  HighsSparseMatrix matrix;
  matrix.format_ = MatrixFormat::kColwise;
  matrix.num_row_ = kDimension;
  matrix.num_col_ = kDimension;
  matrix.start_.resize(kDimension + 1);
  matrix.index_.reserve(kNnz);
  matrix.value_.reserve(kNnz);
  for (HighsInt col = 0; col < kDimension; ++col) {
    matrix.start_[col] = matrix.index_.size();
    if (col != 0) {
      matrix.index_.push_back(col - 1);
      matrix.value_.push_back(-1.0);
    }
    matrix.index_.push_back(col);
    matrix.value_.push_back(4.0);
    if (col + 1 < kDimension) {
      matrix.index_.push_back(col + 1);
      matrix.value_.push_back(-1.0);
    }
  }
  matrix.start_[kDimension] = matrix.index_.size();
  return matrix;
}

void cscAxSkippingZeros(const HighsSparseMatrix& matrix,
                        const std::vector<double>& x,
                        std::vector<double>& y) {
  std::fill(y.begin(), y.end(), 0.0);
  for (HighsInt col = 0; col < matrix.num_col_; ++col) {
    if (x[col] == 0.0) continue;
    for (HighsInt p = matrix.start_[col]; p < matrix.start_[col + 1]; ++p)
      y[matrix.index_[p]] += matrix.value_[p] * x[col];
  }
}

void cscAxSparseView(const HighsSparseMatrix& matrix,
                     const std::vector<HighsInt>& indices,
                     const std::vector<double>& values,
                     std::vector<double>& y) {
  std::fill(y.begin(), y.end(), 0.0);
  for (size_t k = 0; k < indices.size(); ++k) {
    HighsInt col = indices[k];
    for (HighsInt p = matrix.start_[col]; p < matrix.start_[col + 1]; ++p)
      y[matrix.index_[p]] += matrix.value_[p] * values[k];
  }
}

void cscSparseAdd(const HighsSparseMatrix& matrix,
                  const std::vector<HighsInt>& indices,
                  const std::vector<double>& values,
                  std::vector<double>& y) {
  for (size_t k = 0; k < indices.size(); ++k) {
    HighsInt col = indices[k];
    for (HighsInt p = matrix.start_[col]; p < matrix.start_[col + 1]; ++p)
      y[matrix.index_[p]] += matrix.value_[p] * values[k];
  }
}

struct Triplet {
  HighsInt row;
  HighsInt col;
  double value;
  size_t sequence;
};

void transposeCscInto(const HighsSparseMatrix& matrix, HighsSparseMatrix& result,
                      std::vector<HighsInt>& next) {
  std::fill(result.start_.begin(), result.start_.end(), 0);
  for (HighsInt row : matrix.index_) ++result.start_[row + 1];
  for (HighsInt col = 0; col < result.num_col_; ++col)
    result.start_[col + 1] += result.start_[col];
  std::copy(result.start_.begin(), result.start_.end() - 1, next.begin());
  for (HighsInt source_col = 0; source_col < matrix.num_col_; ++source_col) {
    for (HighsInt p = matrix.start_[source_col]; p < matrix.start_[source_col + 1]; ++p) {
      HighsInt destination = next[matrix.index_[p]]++;
      result.index_[destination] = source_col;
      result.value_[destination] = matrix.value_[p];
    }
  }
}

HighsSparseMatrix transposeCsc(const HighsSparseMatrix& matrix) {
  HighsSparseMatrix result;
  result.format_ = MatrixFormat::kColwise;
  result.num_row_ = matrix.num_col_;
  result.num_col_ = matrix.num_row_;
  result.start_.resize(result.num_col_ + 1);
  result.index_.resize(matrix.index_.size());
  result.value_.resize(matrix.value_.size());
  std::vector<HighsInt> next(result.num_col_);
  transposeCscInto(matrix, result, next);
  return result;
}

HighsSparseMatrix buildFromSorted(const std::vector<Triplet>& triplets) {
  HighsSparseMatrix result;
  result.format_ = MatrixFormat::kColwise;
  result.num_row_ = kDimension;
  result.num_col_ = kDimension;
  result.start_.assign(kDimension + 1, 0);
  result.index_.resize(triplets.size());
  result.value_.resize(triplets.size());
  for (const auto& entry : triplets) ++result.start_[entry.col + 1];
  for (HighsInt col = 0; col < kDimension; ++col)
    result.start_[col + 1] += result.start_[col];
  for (size_t i = 0; i < triplets.size(); ++i) {
    result.index_[i] = triplets[i].row;
    result.value_[i] = triplets[i].value;
  }
  return result;
}
}  // namespace

int main() {
  HighsSparseMatrix csc = makeCsc();
  HighsSparseMatrix csr = csc;
  csr.ensureRowwise();
  std::vector<double> dense_x(kDimension, 1.0);
  std::vector<double> sparse_x(kDimension, 0.0);
  std::vector<HighsInt> sparse_indices;
  std::vector<double> sparse_values;
  for (HighsInt i = 0; i < kDimension; i += 20) sparse_x[i] = 1.0;
  for (HighsInt i = 0; i < kDimension; i += 20) {
    sparse_indices.push_back(i);
    sparse_values.push_back(1.0);
  }
  std::vector<double> y(kDimension);

  std::puts("implementation,kernel,dimension,nnz,repeats,total_ns,ns_per_repeat,checksum");
  double sink = 0.0;
  auto start = Clock::now();
  for (int r = 0; r < kProductRepeats; ++r) clobber(y);
  sink = checksum(y);
  report("cpp_reference", "barrier_only", kProductRepeats, start, sink);

  start = Clock::now();
  for (int r = 0; r < kProductRepeats; ++r) {
    std::fill(y.begin(), y.end(), 0.0);
    clobber(y);
  }
  sink = checksum(y);
  report("cpp_reference", "clear_output", kProductRepeats, start, sink);

  start = Clock::now();
  for (int r = 0; r < kProductRepeats; ++r) {
    csc.product(y, dense_x);
    clobber(y);
  }
  sink = checksum(y);
  report("highs", "csc_ax_dense", kProductRepeats, start, sink);

  sink = 0.0;
  start = Clock::now();
  for (int r = 0; r < kProductRepeats; ++r) {
    cscAxSkippingZeros(csc, sparse_x, y);
    clobber(y);
  }
  sink = checksum(y);
  report("cpp_reference", "csc_ax_sparse_skip", kProductRepeats, start, sink);

  start = Clock::now();
  for (int r = 0; r < kProductRepeats; ++r) {
    cscAxSparseView(csc, sparse_indices, sparse_values, y);
    clobber(y);
  }
  sink = checksum(y);
  report("cpp_reference", "csc_ax_sparse_view", kProductRepeats, start, sink);

  std::fill(y.begin(), y.end(), 0.0);
  start = Clock::now();
  for (int r = 0; r < kProductRepeats; ++r) {
    cscSparseAdd(csc, sparse_indices, sparse_values, y);
    clobber(y);
  }
  sink = checksum(y);
  report("cpp_reference", "csc_sparse_add_no_clear", kProductRepeats, start, sink);

  sink = 0.0;
  start = Clock::now();
  for (int r = 0; r < kProductRepeats; ++r) {
    csr.product(y, dense_x);
    clobber(y);
  }
  sink = checksum(y);
  report("highs", "csr_ax_dense", kProductRepeats, start, sink);

  sink = 0.0;
  start = Clock::now();
  for (int r = 0; r < kProductRepeats; ++r) {
    csc.productTranspose(y, dense_x);
    clobber(y);
  }
  sink = checksum(y);
  report("highs", "csc_atx_dense", kProductRepeats, start, sink);

  start = Clock::now();
  for (int r = 0; r < kProductRepeats; ++r) {
    csr.productTranspose(y, dense_x);
    clobber(y);
  }
  sink = checksum(y);
  report("highs", "csr_atx_dense", kProductRepeats, start, sink);

  std::fill(y.begin(), y.end(), 0.0);
  start = Clock::now();
  for (int r = 0; r < kProductRepeats; ++r) {
    csc.alphaProductPlusY(1.0, dense_x, y, false);
    clobber(y);
  }
  sink = checksum(y);
  report("highs", "alpha_ax_plus_y", kProductRepeats, start, sink);

  sink = 0.0;
  start = Clock::now();
  for (int r = 0; r < kQuadRepeats; ++r) {
    csc.productQuad(y, dense_x);
    clobber(y);
  }
  sink = checksum(y);
  report("highs", "product_quad", kQuadRepeats, start, sink);

  HighsScale scale;
  scale.num_row = kDimension;
  scale.num_col = kDimension;
  scale.row.assign(kDimension, 1.0);
  scale.col.assign(kDimension, 1.0);
  start = Clock::now();
  for (int r = 0; r < kProductRepeats; ++r) csc.applyScale(scale);
  sink = checksum(csc.value_);
  report("highs", "apply_scale", kProductRepeats, start, sink);

  std::vector<HighsSparseMatrix> conversions(kTransformRepeats, csc);
  start = Clock::now();
  for (auto& matrix : conversions) matrix.ensureRowwise();
  sink = conversions[0].value_[0];
  report("highs", "csc_to_csr_scratch", kTransformRepeats, start, sink);

  HighsSparseMatrix reusable;
  reusable.format_ = MatrixFormat::kColwise;
  reusable.num_row_ = csc.num_col_;
  reusable.num_col_ = csc.num_row_;
  reusable.start_.resize(reusable.num_col_ + 1);
  reusable.index_.resize(csc.index_.size());
  reusable.value_.resize(csc.value_.size());
  std::vector<HighsInt> reusable_cursor(reusable.num_col_);
  start = Clock::now();
  for (int r = 0; r < kTransformRepeats; ++r)
    transposeCscInto(csc, reusable, reusable_cursor);
  sink = reusable.value_[0];
  report("cpp_reference", "csc_to_csr_into", kTransformRepeats, start, sink);

  std::vector<HighsSparseMatrix> transposes;
  transposes.reserve(kTransformRepeats);
  start = Clock::now();
  for (int r = 0; r < kTransformRepeats; ++r)
    transposes.push_back(transposeCsc(csc));
  sink = transposes[0].value_[0];
  report("cpp_reference", "transpose", kTransformRepeats, start, sink);

  start = Clock::now();
  for (int r = 0; r < kTransformRepeats; ++r)
    transposeCscInto(csc, reusable, reusable_cursor);
  sink = reusable.value_[0];
  report("cpp_reference", "transpose_into", kTransformRepeats, start, sink);

  std::vector<Triplet> sorted_triplets;
  sorted_triplets.reserve(kNnz);
  size_t sequence = 0;
  for (HighsInt col = 0; col < kDimension; ++col) {
    if (col != 0) sorted_triplets.push_back({col - 1, col, -1.0, sequence++});
    sorted_triplets.push_back({col, col, 4.0, sequence++});
    if (col + 1 < kDimension) sorted_triplets.push_back({col + 1, col, -1.0, sequence++});
  }
  start = Clock::now();
  HighsSparseMatrix sorted_matrix = buildFromSorted(sorted_triplets);
  sink = sorted_matrix.value_[sorted_matrix.value_.size() / 2];
  report("cpp_reference", "builder_freeze_sorted", 1, start, sink);

  std::vector<Triplet> triplets;
  triplets.reserve(kNnz);
  sequence = 0;
  for (HighsInt col = kDimension; col-- > 0;) {
    if (col + 1 < kDimension) triplets.push_back({col + 1, col, -1.0, sequence++});
    triplets.push_back({col, col, 4.0, sequence++});
    if (col != 0) triplets.push_back({col - 1, col, -1.0, sequence++});
  }
  start = Clock::now();
  std::sort(triplets.begin(), triplets.end(), [](const Triplet& a, const Triplet& b) {
    if (a.col != b.col) return a.col < b.col;
    if (a.row != b.row) return a.row < b.row;
    return a.sequence < b.sequence;
  });
  HighsSparseMatrix general_matrix = buildFromSorted(triplets);
  sink = general_matrix.value_[general_matrix.value_.size() / 2];
  report("cpp_std_sort", "builder_freeze_general", 1, start, sink);

  HighsSparseVectorSum accumulator(kDimension);
  start = Clock::now();
  for (int r = 0; r < kAccumulatorRepeats; ++r) {
    accumulator.clear();
    for (HighsInt index = 0; index < kDimension; ++index) {
      accumulator.add(index, 1.0);
      accumulator.add(index, -0.5);
    }
    asm volatile("" : : "g"(&accumulator) : "memory");
  }
  sink = accumulator.getValue(kDimension / 2);
  report("highs", "sparse_accumulate", kAccumulatorRepeats, start, sink);
  return 0;
}
