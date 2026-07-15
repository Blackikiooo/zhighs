# Model parser benchmark

The benchmark compares equal end-to-end work:

```text
open file -> read bytes -> parse syntax -> construct complete in-memory model
```

Solving, presolve, model destruction, and result printing are outside the timed
region. Both runners validate row, column, and nonzero counts and calculate the
same objective-plus-matrix checksum. Each process reports best time, median
time, input throughput, and `ru_maxrss`.

Build the Zig runner:

```sh
zig build -Doptimize=ReleaseFast build-io-bench
```

Run a fair Release HiGHS comparison:

```sh
bench/io/run_comparison.sh /path/to/model.mps 9 3
```

`HIGHS_ROOT` defaults to
`/home/godv/documents/codefiles/cppfiles/HiGHS`. `HIGHS_BUILD` defaults to its
`cmake-build-release` directory. The script refuses to compare against a
missing release library.

## Initial local baseline

HiGHS 1.14.0 was built with `CMAKE_BUILD_TYPE=Release`, shared libraries, and
IPO/LTO. zhighs used `ReleaseFast`. Results are medians from nine measured runs
after three warmups.

| dataset | format | zhighs | HiGHS | speedup | zhighs RSS | HiGHS RSS |
|---|---:|---:|---:|---:|---:|---:|
| 80bau3b | MPS | 5.902 ms | 26.625 ms | 4.51x | 5,076 KiB | 10,124 KiB |
| 2122 | LP | 0.934 ms | 3.409 ms | 3.65x | 1,696 KiB | 7,412 KiB |

These are machine-local baselines, not permanent acceptance thresholds. A
performance gate should use a pinned corpus, CPU affinity/governor policy, and
several large instances including ranges, integer markers, long names, and
adversarial duplicate coordinates.
