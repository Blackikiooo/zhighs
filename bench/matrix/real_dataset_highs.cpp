#include <sys/resource.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include <cstdlib>
#ifdef ZHIGHS_RSS_SPAN_PARSER
#include <span>
#include <spanstream>
#endif

#include "util/HighsSparseMatrix.h"

namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;

struct Triplet {
  HighsInt row;
  HighsInt col;
  double value;
};

static HighsSparseMatrix read_matrix_market(const fs::path& path) {
#ifdef ZHIGHS_RSS_SPAN_PARSER
  // Match the Zig benchmark parser's lifetime: retain the complete Matrix
  // Market text while triplets are built and sorted. This makes peak-RSS
  // comparisons about data structures rather than stream-vs-whole-file I/O.
  std::ifstream file(path, std::ios::binary);
  if (!file) throw std::runtime_error("cannot open " + path.string());
  file.seekg(0, std::ios::end);
  const std::streamsize content_size = file.tellg();
  file.seekg(0, std::ios::beg);
  std::string content(static_cast<std::size_t>(content_size), '\0');
  if (!file.read(content.data(), content_size))
    throw std::runtime_error("cannot read " + path.string());
  std::ispanstream input(std::span<const char>(content.data(), content.size()));
#else
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open " + path.string());
#endif
  std::string line, banner, object, format, field, symmetry;
  if (!std::getline(input, line)) throw std::runtime_error("missing banner");
  std::istringstream header(line);
  header >> banner >> object >> format >> field >> symmetry;
  if (banner != "%%MatrixMarket" || object != "matrix" || format != "coordinate")
    throw std::runtime_error("unsupported Matrix Market header");
  const bool pattern = field == "pattern";
  const bool symmetric = symmetry == "symmetric" || symmetry == "hermitian";
  if (!pattern && field != "real" && field != "integer")
    throw std::runtime_error("unsupported Matrix Market field");

  do {
    if (!std::getline(input, line)) throw std::runtime_error("missing dimensions");
  } while (line.empty() || line[0] == '%');
  std::istringstream dimensions(line);
  HighsInt rows, cols;
  std::size_t stored_nnz;
  dimensions >> rows >> cols >> stored_nnz;

  std::vector<Triplet> entries;
  entries.reserve(symmetric ? stored_nnz * 2 : stored_nnz);
  HighsInt row, col;
  double value;
  std::size_t seen = 0;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '%') continue;
    std::istringstream entry(line);
    entry >> row >> col;
    value = 1.0;
    if (!pattern) entry >> value;
    --row;
    --col;
    entries.push_back({row, col, value});
    if (symmetric && row != col) entries.push_back({col, row, value});
    ++seen;
  }
  if (seen != stored_nnz) throw std::runtime_error("stored nnz mismatch");

  std::sort(entries.begin(), entries.end(), [](const Triplet& lhs, const Triplet& rhs) {
    return std::tie(lhs.col, lhs.row) < std::tie(rhs.col, rhs.row);
  });
  std::vector<Triplet> canonical;
  canonical.reserve(entries.size());
  for (const auto& entry : entries) {
    if (!canonical.empty() && canonical.back().row == entry.row && canonical.back().col == entry.col)
      canonical.back().value += entry.value;
    else
      canonical.push_back(entry);
  }
  canonical.erase(std::remove_if(canonical.begin(), canonical.end(), [](const Triplet& entry) {
                    return entry.value == 0.0;
                  }),
                  canonical.end());

  HighsSparseMatrix matrix;
  matrix.format_ = MatrixFormat::kColwise;
  matrix.num_row_ = rows;
  matrix.num_col_ = cols;
  matrix.start_.assign(static_cast<std::size_t>(cols) + 1, 0);
  matrix.index_.resize(canonical.size());
  matrix.value_.resize(canonical.size());
  for (const auto& entry : canonical) ++matrix.start_[entry.col + 1];
  for (HighsInt index = 0; index < cols; ++index) matrix.start_[index + 1] += matrix.start_[index];
  for (std::size_t index = 0; index < canonical.size(); ++index) {
    matrix.index_[index] = canonical[index].row;
    matrix.value_[index] = canonical[index].value;
  }
  return matrix;
}

static double max_relative_difference(const std::vector<double>& lhs,
                                      const std::vector<double>& rhs) {
  double result = 0.0;
  for (std::size_t index = 0; index < lhs.size(); ++index) {
    const double scale = std::max({1.0, std::abs(lhs[index]), std::abs(rhs[index])});
    result = std::max(result, std::abs(lhs[index] - rhs[index]) / scale);
  }
  return result;
}

static double elapsed_ms(Clock::time_point started, std::size_t repeats = 1) {
  return std::chrono::duration<double, std::milli>(Clock::now() - started).count() /
         static_cast<double>(repeats);
}

static std::size_t current_rss_kb() {
  std::ifstream status("/proc/self/status");
  std::string line;
  while (std::getline(status, line)) {
    if (line.rfind("VmRSS:", 0) != 0) continue;
    std::istringstream fields(line.substr(6));
    std::size_t value;
    fields >> value;
    return value;
  }
  throw std::runtime_error("cannot read VmRSS");
}

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: real-dataset-highs DATASET_DIR REPORT.tsv\n";
    return 2;
  }
  std::vector<fs::directory_entry> datasets;
  for (const auto& entry : fs::directory_iterator(argv[1])) {
    if (entry.is_regular_file() && entry.path().extension() == ".mtx" &&
        entry.file_size() > 1024 * 1024)
      datasets.push_back(entry);
  }
  std::sort(datasets.begin(), datasets.end(), [](const auto& lhs, const auto& rhs) {
    return lhs.file_size() < rhs.file_size();
  });
  std::ofstream report(argv[2]);
  const bool rss_only = std::getenv("ZHIGHS_DATASET_RSS_ONLY") != nullptr;
  const char* filter = std::getenv("ZHIGHS_DATASET_FILTER");
  if (rss_only)
    report << "implementation\tdataset\trows\tcols\tnnz\trequested_bytes"
              "\tcurrent_rss_kb\tpeak_rss_kb\n";
  else
    report << "implementation\tdataset\trows\tcols\tnnz\tcsc_spmv_ms\tcsr_spmv_ms"
              "\tcsc_to_csr_ms\tpeak_rss_kb\n";
  report.setf(std::ios::fixed);
  report.precision(3);

  for (const auto& dataset : datasets) {
    if (filter && dataset.path().filename() != filter) continue;
    std::cerr << "benchmarking HiGHS " << dataset.path().filename() << "\n";
    HighsSparseMatrix matrix = read_matrix_market(dataset.path());
    if (rss_only) {
      rusage usage{};
      getrusage(RUSAGE_SELF, &usage);
      const std::size_t requested = matrix.start_.size() * sizeof(HighsInt) +
                                    matrix.index_.size() * sizeof(HighsInt) +
                                    matrix.value_.size() * sizeof(double);
      report << "HiGHS\t" << dataset.path().filename().string() << '\t' << matrix.num_row_
             << '\t' << matrix.num_col_ << '\t' << matrix.numNz() << '\t' << requested
             << '\t' << current_rss_kb() << '\t' << usage.ru_maxrss << '\n';
      continue;
    }
    std::vector<double> x(matrix.num_col_);
    for (std::size_t index = 0; index < x.size(); ++index)
      x[index] = 1.0 + static_cast<double>(index % 17) * 0.03125;
    std::vector<double> csc_result;
    const std::size_t repeats = std::clamp<std::size_t>(100000000 / std::max<HighsInt>(matrix.numNz(), 1), 5, 50);
    matrix.product(csc_result, x);
    auto started = Clock::now();
    for (std::size_t repeat = 0; repeat < repeats; ++repeat) matrix.product(csc_result, x);
    const double csc_ms = elapsed_ms(started, repeats);

    HighsSparseMatrix rowwise;
    constexpr std::size_t transform_repeats = 7;
    rowwise.createRowwise(matrix);
    started = Clock::now();
    for (std::size_t repeat = 0; repeat < transform_repeats; ++repeat)
      rowwise.createRowwise(matrix);
    const double conversion_ms = elapsed_ms(started, transform_repeats);
    std::vector<double> csr_result;
    rowwise.product(csr_result, x);
    if (max_relative_difference(csc_result, csr_result) > 1e-12)
      throw std::runtime_error("CSC/CSR semantic mismatch");
    started = Clock::now();
    for (std::size_t repeat = 0; repeat < repeats; ++repeat) rowwise.product(csr_result, x);
    const double csr_ms = elapsed_ms(started, repeats);

    rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    report << "HiGHS\t" << dataset.path().filename().string() << '\t' << matrix.num_row_ << '\t'
           << matrix.num_col_ << '\t' << matrix.numNz() << '\t' << csc_ms << '\t' << csr_ms
           << '\t' << conversion_ms << '\t' << usage.ru_maxrss << '\n';
  }
}
