// HiGHS HFactor INVERT benchmark with the same cyclic tridiagonal basis used
// by sparse_lu_bench.zig. Setup/allocation is outside the measured region.

#include "util/HFactor.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

int main(int argc, char** argv) {
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
