// HiGHS HFactor INVERT benchmark with the same cyclic tridiagonal basis used
// by sparse_lu_bench.zig. Setup/allocation is outside the measured region.

#include "util/HFactor.h"
#include "Highs.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <vector>

static int netlib_basis(const char* model_path, int repeats,
                        const char* export_path) {
  Highs highs;
  highs.setOptionValue("output_flag", false);
  highs.setOptionValue("solver", "simplex");
  highs.setOptionValue("presolve", "off");
  if (highs.readModel(model_path) != HighsStatus::kOk ||
      highs.run() != HighsStatus::kOk)
    return 1;
  const HighsLp& lp = highs.getLp();
  const HighsBasis& basis = highs.getBasis();
  const HighsInt n = lp.num_row_;
  std::vector<HighsInt> basic;
  for (HighsInt col = 0; col < lp.num_col_; ++col)
    if (basis.col_status[col] == HighsBasisStatus::kBasic) basic.push_back(col);
  for (HighsInt row = 0; row < n; ++row)
    if (basis.row_status[row] == HighsBasisStatus::kBasic)
      basic.push_back(lp.num_col_ + row);
  if (basic.size() != static_cast<size_t>(n)) return 3;

  HFactor factor;
  factor.setup(lp.num_col_, n, lp.a_matrix_.start_.data(),
               lp.a_matrix_.index_.data(), lp.a_matrix_.value_.data(),
               basic.data());
  if (factor.build() != 0) return 4;
  std::vector<std::uint64_t> samples(repeats);
  for (int run = 0; run < repeats; ++run) {
    const auto begin = std::chrono::steady_clock::now();
    const HighsInt deficiency = factor.build();
    const auto end = std::chrono::steady_clock::now();
    if (deficiency != 0) return 4;
    samples[run] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
  }
  std::sort(samples.begin(), samples.end());

  // Export the exact final simplex basis in deterministic CSC text form.
  std::ofstream out(export_path);
  HighsInt basis_nnz = 0;
  for (HighsInt variable : basic)
    basis_nnz += variable < lp.num_col_
                     ? lp.a_matrix_.start_[variable + 1] - lp.a_matrix_.start_[variable]
                     : 1;
  out << "ZHIGHS_BASIS_V1 " << n << ' ' << basis_nnz << '\n';
  for (HighsInt variable : basic) {
    if (variable < lp.num_col_) {
      const HighsInt begin = lp.a_matrix_.start_[variable];
      const HighsInt end = lp.a_matrix_.start_[variable + 1];
      out << end - begin;
      std::vector<std::pair<HighsInt, double>> entries;
      for (HighsInt p = begin; p < end; ++p)
        entries.emplace_back(lp.a_matrix_.index_[p], lp.a_matrix_.value_[p]);
      std::sort(entries.begin(), entries.end());
      for (const auto& entry : entries)
        out << ' ' << entry.first << ' ' << std::setprecision(17) << entry.second;
      out << '\n';
    } else {
      out << "1 " << variable - lp.num_col_ << " 1\n";
    }
  }
  std::cout << "highs-netlib," << model_path << ',' << n << ',' << basis_nnz
            << ',' << samples[samples.size() / 2] << ',' << export_path << '\n';
  return out ? 0 : 5;
}

int main(int argc, char** argv) {
  if (argc >= 2 && std::string(argv[1]) == "--netlib") {
    if (argc != 5) return 2;
    return netlib_basis(argv[2], std::stoi(argv[3]), argv[4]);
  }
  const HighsInt n = argc > 1 ? std::stoll(argv[1]) : 512;
  const int repeats = argc > 2 ? std::stoi(argv[2]) : 21;
  if (n < 3 || repeats <= 0) return 2;
  std::vector<HighsInt> starts(n + 1), rows(n * 3), basic(n);
  std::vector<double> values(n * 3);
  HighsInt output = 0;
  for (HighsInt column = 0; column < n; ++column) {
    starts[column] = output;
    HighsInt local[3] = {column == 0 ? n - 1 : column - 1, column,
                         column + 1 == n ? 0 : column + 1};
    std::sort(local, local + 3);
    for (HighsInt row : local) {
      rows[output] = row;
      values[output] = row == column ? 4.0 + (column % 7) * 0.125 : -0.5;
      ++output;
    }
  }
  starts[n] = output;
  std::iota(basic.begin(), basic.end(), 0);
  HFactor factor;
  factor.setup(n, n, starts.data(), rows.data(), values.data(), basic.data());
  if (factor.build() != 0) return 1;
  std::vector<std::uint64_t> samples(repeats);
  for (int run = 0; run < repeats; ++run) {
    std::iota(basic.begin(), basic.end(), 0);
    const auto begin = std::chrono::steady_clock::now();
    const HighsInt deficiency = factor.build();
    const auto end = std::chrono::steady_clock::now();
    if (deficiency != 0) return 1;
    samples[run] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
  }
  std::sort(samples.begin(), samples.end());
  std::cout << "highs," << n << ',' << output << ',' << repeats << ','
            << samples[samples.size() / 2] << '\n';
}
