// End-to-end HiGHS MPS + serial simplex reference runner.
#include "Highs.h"
#include <chrono>
#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
  if (argc != 2) return 2;
  const auto total_started = std::chrono::steady_clock::now();
  Highs highs;
  // Keep the reference path unambiguously equivalent to zhighs' current
  // single-threaded, presolve-free dual simplex path. In particular, do not
  // leave the global thread count or simplex strategy at their automatic
  // values even though parallel=off currently selects the serial solver.
  if (highs.setOptionValue("output_flag", false) == HighsStatus::kError ||
      highs.setOptionValue("solver", "simplex") == HighsStatus::kError ||
      highs.setOptionValue("presolve", "off") == HighsStatus::kError ||
      highs.setOptionValue("parallel", "off") == HighsStatus::kError ||
      highs.setOptionValue("threads", 1) == HighsStatus::kError ||
      highs.setOptionValue("simplex_strategy", 1) == HighsStatus::kError) {
    return 5;
  }
  const auto read_started = std::chrono::steady_clock::now();
  // HiGHS may return kWarning for a usable model, for example when tiny
  // coefficients are ignored. Only a real read error invalidates the run.
  if (highs.readModel(argv[1]) == HighsStatus::kError) return 3;
  const auto solve_started = std::chrono::steady_clock::now();
  const auto read_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(solve_started - read_started).count();
  if (highs.run() == HighsStatus::kError) return 4;
  const auto finished = std::chrono::steady_clock::now();
  const auto& info = highs.getInfo();
  const auto solve_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(finished - solve_started).count();
  const auto total_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(finished - total_started).count();
  std::cout << std::setprecision(17) << "highs\t" << argv[1] << '\t'
            << highs.modelStatusToString(highs.getModelStatus()) << '\t' << info.objective_function_value << '\t'
            << info.simplex_iteration_count << '\t' << info.max_primal_infeasibility << '\t'
            << info.max_dual_infeasibility << '\t' << read_ns << '\t' << solve_ns << '\t' << total_ns << '\n';
}
