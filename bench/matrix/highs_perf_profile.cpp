#define main highs_full_benchmark_main
#include "highs_matrix_bench.cpp"
#undef main

#include <cstdlib>
#include <cstring>

int main() {
  const char* kernel = std::getenv("ZHIGHS_PERF_KERNEL");
  if (!kernel) kernel = "csc_ax_skip";

  HighsSparseMatrix csc = makeCsc();
  HighsSparseMatrix csr = csc;
  csr.ensureRowwise();
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

  if (!std::strcmp(kernel, "csc_ax_skip")) {
    for (int r = 0; r < 2000; ++r) {
      cscAxSkippingZeros(csc, sparse_x, y);
      clobber(y);
    }
  } else if (!std::strcmp(kernel, "csc_ax_sparse_view")) {
    for (int r = 0; r < 4000; ++r) {
      cscAxSparseView(csc, sparse_indices, sparse_values, y);
      clobber(y);
    }
  } else if (!std::strcmp(kernel, "csc_sparse_add")) {
    std::fill(y.begin(), y.end(), 0.0);
    for (int r = 0; r < 4000; ++r) {
      cscSparseAdd(csc, sparse_indices, sparse_values, y);
      clobber(y);
    }
  } else if (!std::strcmp(kernel, "csr_ax")) {
    for (int r = 0; r < 1000; ++r) {
      csr.product(y, dense_x);
      clobber(y);
    }
  } else if (!std::strcmp(kernel, "csr_atx")) {
    for (int r = 0; r < 1000; ++r) {
      csr.productTranspose(y, dense_x);
      clobber(y);
    }
  } else if (!std::strcmp(kernel, "product_quad")) {
    for (int r = 0; r < 500; ++r) {
      csc.productQuad(y, dense_x);
      clobber(y);
    }
  } else if (!std::strcmp(kernel, "csc_to_csr_into")) {
    HighsSparseMatrix result;
    result.format_ = MatrixFormat::kColwise;
    result.num_row_ = csc.num_col_;
    result.num_col_ = csc.num_row_;
    result.start_.resize(result.num_col_ + 1);
    result.index_.resize(csc.index_.size());
    result.value_.resize(csc.value_.size());
    std::vector<HighsInt> cursor(result.num_col_);
    for (int r = 0; r < 500; ++r) {
      transposeCscInto(csc, result, cursor);
      asm volatile("" : : "g"(result.value_.data()) : "memory");
    }
  } else if (!std::strcmp(kernel, "sparse_accumulate")) {
    HighsSparseVectorSum accumulator(kDimension);
    for (int r = 0; r < 500; ++r) {
      accumulator.clear();
      for (HighsInt i = 0; i < kDimension; ++i) {
        accumulator.add(i, 1.0);
        accumulator.add(i, -0.5);
      }
      asm volatile("" : : "g"(&accumulator) : "memory");
    }
  } else {
    std::fprintf(stderr, "unknown ZHIGHS_PERF_KERNEL=%s\n", kernel);
    return 2;
  }

  std::printf("%.17g\n", checksum(y));
  return 0;
}
