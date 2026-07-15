// End-to-end HiGHS model reader benchmark with the same timing boundary and
// output columns as parser_bench.zig.

#include "Highs.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <sys/resource.h>
#include <vector>

int main(int argc, char** argv) {
  if (argc < 2 || argc > 4) {
    std::cerr << "usage: highs-parser-bench MODEL [iterations=7] [warmups=2]\n";
    return 2;
  }
  const std::string path = argv[1];
  const int iterations = argc >= 3 ? std::stoi(argv[2]) : 7;
  const int warmups = argc >= 4 ? std::stoi(argv[3]) : 2;
  if (iterations <= 0 || warmups < 0) return 2;
  std::vector<std::uint64_t> samples(iterations);
  HighsInt rows = 0, columns = 0, nonzeros = 0;
  double checksum = 0.0;
  for (int run = 0; run < warmups + iterations; ++run) {
    const auto started = std::chrono::steady_clock::now();
    Highs highs;
    highs.setOptionValue("output_flag", false);
    const HighsStatus status = highs.readModel(path);
    const auto stopped = std::chrono::steady_clock::now();
    if (status == HighsStatus::kError) {
      std::cerr << "HiGHS failed to read " << path << '\n';
      return 1;
    }
    const HighsLp& lp = highs.getLp();
    rows = lp.num_row_;
    columns = lp.num_col_;
    nonzeros = lp.a_matrix_.numNz();
    checksum = lp.offset_;
    for (double value : lp.col_cost_) checksum += value;
    for (double value : lp.a_matrix_.value_) checksum += value;
    if (run >= warmups) {
      samples[run - warmups] = std::chrono::duration_cast<std::chrono::nanoseconds>(stopped - started).count();
    }
  }
  std::sort(samples.begin(), samples.end());
  const auto best = samples.front();
  const auto median = samples[samples.size() / 2];
  const double mib = static_cast<double>(std::filesystem::file_size(path)) / (1024.0 * 1024.0);
  const double throughput = mib / (static_cast<double>(median) / 1e9);
  struct rusage usage {};
  getrusage(RUSAGE_SELF, &usage);
  std::cout << "highs\t" << path << '\t' << rows << '\t' << columns << '\t' << nonzeros << '\t'
            << best << '\t' << std::fixed << std::setprecision(3) << static_cast<double>(median) / 1e6
            << '\t' << throughput << '\t' << usage.ru_maxrss << '\t' << std::setprecision(17) << checksum << '\n';
}
