#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

#include "util/HighsCDouble.h"

static constexpr std::size_t n = 4096;
static constexpr std::size_t repeats = 20000;

static double make_hd(std::size_t i) {
  std::uint64_t x = static_cast<std::uint64_t>(i) * 6364136223846793005ull +
                    1442695040888963407ull;
  double mant = static_cast<double>((x >> 12) & ((1ull << 40) - 1)) /
                static_cast<double>(1ull << 40);
  double sign = (x & 1) == 0 ? 1.0 : -1.0;
  double scale = static_cast<double>((x >> 52) & 15) * 0.0625 + 0.5;
  return sign * (mant + 0.125) * scale;
}

static double make_small_hd(std::size_t i) { return make_hd(i) * 1e-12; }

static void fill_hd(std::vector<double>& values, bool small) {
  for (std::size_t i = 0; i < values.size(); ++i)
    values[i] = small ? make_small_hd(i + 1) : make_hd(i + 1);
}

static void fill_hcd(std::vector<HighsCDouble>& values, bool small) {
  for (std::size_t i = 0; i < values.size(); ++i) {
    double hi = small ? make_small_hd(i + 1) : make_hd(i + 1);
    values[i] = HighsCDouble(hi) + make_small_hd(i + 10001);
  }
}

template <typename F>
static void run(const char* name, F&& f) {
  auto start = std::chrono::steady_clock::now();
  double checksum = f();
  auto end = std::chrono::steady_clock::now();
  auto ns =
      std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
  double ns_per_op = static_cast<double>(ns) / static_cast<double>(n * repeats);
  std::cout << name << "," << ns << "," << std::fixed << std::setprecision(6)
            << ns_per_op << "," << std::scientific << std::setprecision(17)
            << checksum << "\n";
}

int main() {
  std::vector<double> hd_values(n);
  std::vector<double> small_hd_values(n);
  std::vector<HighsCDouble> hcd_values(n);

  fill_hd(hd_values, false);
  fill_hd(small_hd_values, true);
  fill_hcd(hcd_values, false);

  std::cout << "name,total_ns,ns_per_op,checksum\n";

  run("cpp.add_hd_assign", [&] {
    HighsCDouble x = HighsCDouble(1.0) + 1e-16;
    for (std::size_t r = 0; r < repeats; ++r)
      for (double v : hd_values) x += v;
    return static_cast<double>(x);
  });

  run("cpp.add_hd_ordered_assign", [&] {
    HighsCDouble x = HighsCDouble(1.0e8) + 1e-16;
    for (std::size_t r = 0; r < repeats; ++r)
      for (double v : small_hd_values) x += v;
    return static_cast<double>(x);
  });

  run("cpp.add_hcd_assign", [&] {
    HighsCDouble x = HighsCDouble(1.0) + 1e-16;
    for (std::size_t r = 0; r < repeats; ++r)
      for (const HighsCDouble& v : hcd_values) x += v;
    return static_cast<double>(x);
  });

  run("cpp.multiply_hd_assign", [&] {
    double checksum = 0.0;
    for (std::size_t r = 0; r < repeats; ++r) {
      HighsCDouble x =
          HighsCDouble(1.0000001 + static_cast<double>(r) * 1e-18) + 1e-16;
      for (double v : hd_values) x *= 1.0 + std::abs(v) * 1e-12;
      checksum += static_cast<double>(x);
    }
    return checksum;
  });

  run("cpp.multiply_hcd_assign", [&] {
    double checksum = 0.0;
    for (std::size_t r = 0; r < repeats; ++r) {
      HighsCDouble x =
          HighsCDouble(1.0000001 + static_cast<double>(r) * 1e-18) + 1e-16;
      for (const HighsCDouble& v : hcd_values)
        x *= HighsCDouble(1.0 + std::abs(static_cast<double>(v)) * 1e-12);
      checksum += static_cast<double>(x);
    }
    return checksum;
  });

  run("cpp.divide_hd_assign", [&] {
    double checksum = 0.0;
    for (std::size_t r = 0; r < repeats; ++r) {
      HighsCDouble x =
          HighsCDouble(1.0000001 + static_cast<double>(r) * 1e-18) + 1e-16;
      for (double v : hd_values) x /= 1.0 + std::abs(v) * 1e-12;
      checksum += static_cast<double>(x);
    }
    return checksum;
  });

  run("cpp.divide_hcd_assign", [&] {
    double checksum = 0.0;
    for (std::size_t r = 0; r < repeats; ++r) {
      HighsCDouble x =
          HighsCDouble(1.0000001 + static_cast<double>(r) * 1e-18) + 1e-16;
      for (const HighsCDouble& v : hcd_values)
        x /= HighsCDouble(1.0 + std::abs(static_cast<double>(v)) * 1e-12);
      checksum += static_cast<double>(x);
    }
    return checksum;
  });
}
