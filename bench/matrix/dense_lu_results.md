# DenseLU FTRAN/BTRAN benchmark

## Reproduction

Build and run the independent dimension sweep:

```sh
zig build bench-dense-lu -Doptimize=ReleaseFast -Dcpu=native
```

The executable prints CSV and accepts these environment variables:

- `ZHIGHS_LU_KERNEL=ftran|btran` selects exactly one solve path.
- `ZHIGHS_LU_DIMENSION=N` selects exactly one basis dimension.
- `ZHIGHS_LU_REPEATS=N` fixes the timed repetition count.

The factorization and workspace allocation happen before the timed solve loop.
Every repetition restores the same RHS, invokes one allocation-free solve, and
adds a result element to a checksum so the optimizer cannot discard the work.

## ReleaseFast timing snapshot

Recorded on 2026-07-13 with Zig 0.16.0, Linux 6.17, `-Dcpu=native`, and an AMD
Ryzen 5 3500X. These are single-run development-machine numbers, so they are a
regression snapshot rather than a stable cross-machine performance target.

| basis dimension | LU bytes | repeats | FTRAN ns/solve | BTRAN ns/solve |
| ---: | ---: | ---: | ---: | ---: |
| 32  | 8,192     | 20,000 | 2,521.435   | 608.819 |
| 64  | 32,768    | 16,384 | 10,529.049  | 2,417.523 |
| 128 | 131,072   | 4,096  | 44,168.225  | 8,996.061 |
| 256 | 524,288   | 1,024  | 183,832.660 | 31,984.861 |
| 512 | 2,097,152 | 256    | 735,482.082 | 122,643.098 |

Doubling the dimension increases FTRAN time by about 4.0--4.2x and BTRAN by
about 3.6--4.0x in this sweep, consistent with the quadratic triangular-solve
work. The row-contiguous BTRAN is substantially faster than the current FTRAN
on this machine; the benchmark intentionally reports the paths independently
so future loop-layout changes can be compared without mixing them.

## Hardware cache/TLB counters

Use a large repeat count so the one-time LU factorization is amortized in the
whole-process `perf stat` counters. Repeat the following for both kernels and
the desired dimensions:

```sh
ZHIGHS_LU_KERNEL=btran \
ZHIGHS_LU_DIMENSION=512 \
ZHIGHS_LU_REPEATS=4000 \
perf stat -x, -r 5 \
  -e cycles,instructions,cache-references,cache-misses,dTLB-loads,dTLB-load-misses \
  zig-out/bin/dense-lu-bench
```

Record both absolute counts and these ratios:

```text
cache miss rate = cache-misses / cache-references
dTLB miss rate  = dTLB-load-misses / dTLB-loads
```

Recorded on 2026-07-14 with `sudo perf stat -r 5`. Values are the five-run
means reported by perf. The six events were multiplexed at 82--83% running
time, and perf scaled the counts accordingly. The one-time factorization is
included in the process counters, but 20,000 solves at dimension 64 and 4,000
solves at dimensions 256/512 amortize it.

| kernel | dimension | repeats | cache references | cache misses | cache miss rate | dTLB loads | dTLB misses | dTLB miss rate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| FTRAN | 64  | 20,000 | 4,122,403   | 98,593     | 2.39%  | 10,167    | 1,428  | 14.05% |
| BTRAN | 64  | 20,000 | 3,360,185   | 91,675     | 2.73%  | 10,966    | 772    | 7.04%  |
| FTRAN | 256 | 4,000  | 70,594,794  | 4,190,986  | 5.94%  | 750,499   | 27,476 | 3.66%  |
| BTRAN | 256 | 4,000  | 68,246,330  | 5,122,700  | 7.51%  | 693,637   | 21,567 | 3.11%  |
| FTRAN | 512 | 4,000  | 339,958,110 | 24,531,788 | 7.22%  | 4,222,592 | 50,017 | 1.18%  |
| BTRAN | 512 | 4,000  | 334,027,954 | 40,506,701 | 12.13% | 4,090,994 | 43,911 | 1.07%  |

The generic cache events do not show a lower miss rate for row-contiguous
BTRAN: its miss rate is slightly higher at 64/256 and reaches 12.13% at 512.
This is partly a reminder that generic `cache-*` events on this AMD CPU are not
a direct L1-data-cache measurement. BTRAN nevertheless executes far fewer
instructions and sustains much higher IPC, which explains its elapsed-time
advantage:

| kernel | dimension | cycles | instructions | IPC | mean ns/solve |
| --- | ---: | ---: | ---: | ---: | ---: |
| FTRAN | 64  | 829,567,983    | 1,026,181,144  | 1.24 | 10,633.920 |
| BTRAN | 64  | 190,225,815    | 635,762,173    | 3.34 | 2,553.286  |
| FTRAN | 256 | 2,845,796,368  | 3,227,715,889  | 1.13 | 181,814.694 |
| BTRAN | 256 | 506,733,646    | 1,926,113,276  | 3.80 | 31,950.058  |
| FTRAN | 512 | 11,538,825,189 | 13,013,887,690 | 1.13 | 734,901.596 |
| BTRAN | 512 | 1,972,838,223  | 7,780,207,784  | 3.94 | 125,074.050 |

At 512 dimensions BTRAN is about 5.88x faster in the five-run measurement.
Its dTLB-load-miss count is also 12.2% lower than FTRAN, while the cache-miss
count is higher. More precise cache-level attribution would require AMD model-
specific L1D/L2/LLC events instead of only perf's portable generic events.
