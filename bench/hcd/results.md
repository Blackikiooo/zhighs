# HCD Benchmark Results

This benchmark compares `src/foundation/double.zig` against the local HiGHS
`HighsCDouble` implementation at:

`/home/godv/documents/codefiles/cppfiles/HiGHS/highs/util/HighsCDouble.h`

## Setup

- Date: 2026-07-10
- Workload: `4096 * 20000` operations per case
- Zig build: `ReleaseFast`
- C++ build: `g++ -O3 -march=native -DNDEBUG`
- Metric: nanoseconds per operation, lower is better

## Hardware And Toolchain

- OS: Ubuntu 24.04, Linux `6.17.0-35-generic`, x86_64
- CPU: AMD Ryzen 5 3500X 6-Core Processor
- Cores / threads: 6 cores, 6 threads
- Max / min frequency: 4120.3120 MHz / 2200.0000 MHz
- Frequency boost: enabled
- L1d / L1i cache: 192 KiB / 192 KiB, 6 instances each
- L2 cache: 3 MiB, 6 instances
- L3 cache: 32 MiB, 2 instances
- Relevant CPU features: FMA, AVX, AVX2, SSE4.1, SSE4.2, BMI1, BMI2
- Zig: 0.16.0
- C++ compiler: g++ 13.3.0

## Commands

```bash
g++ -O3 -march=native -DNDEBUG \
  -I/home/godv/documents/codefiles/cppfiles/HiGHS/highs \
  bench/hcd/highs_cdouble_bench.cpp \
  -o bench/hcd/highs_cdouble_bench

bench/hcd/highs_cdouble_bench
zig build bench-hcd -Doptimize=ReleaseFast
```

## Results

| Operation | HiGHS C++ ns/op | Zig HCD ns/op | Notes |
|---|---:|---:|---|
| add HD | 1.249518 | 1.240079 | Zig slightly faster |
| add HD ordered fast | 1.271427 | 1.007600 | Zig faster via `addHDOrderedFast` |
| add HCD | 1.565812 | 1.553363 | Zig slightly faster |
| multiply HD | 6.544251 | 3.894990 | Zig faster, likely helped by FMA `twoProduct` |
| multiply HCD | 6.180917 | 4.614730 | Zig faster |
| divide HD | 15.411998 | 10.961034 | Zig faster |
| divide HCD safe | 17.335735 | 19.578694 | Zig safe path is slower |
| divide HCD fast | N/A | 11.636426 | Zig fast one-correction path |

## Raw Output

### HiGHS C++

```text
name,total_ns,ns_per_op,checksum
cpp.add_hd_assign,102360515,1.249518,2.52968143199691694e+04
cpp.add_hd_ordered_assign,104155294,1.271427,1.00000000000000030e+08
cpp.add_hcd_assign,128271285,1.565812,2.52968143199332735e+04
cpp.multiply_hd_assign,536105040,6.544251,2.00000020494413620e+04
cpp.multiply_hcd_assign,506340756,6.180917,2.00000020494413620e+04
cpp.divide_hd_assign,1262550888,15.411998,2.00000019505624696e+04
cpp.divide_hcd_assign,1420143386,17.335735,2.00000019505624696e+04
```

### Zig

```text
name,total_ns,ns_per_op,checksum
zig.add_hd_fast,101587262,1.240079,25296.81431996917000000
zig.add_hd_ordered_fast,82542608,1.007600,100000000.00000003000000000
zig.add_hcd_fast,127251485,1.553363,25296.81431993327400000
zig.multiply_hd,319077572,3.894990,20000.00204944136200000
zig.multiply_hcd,378038699,4.614730,20000.00204944136200000
zig.divide_hd,897927889,10.961034,20000.00195056247000000
zig.divide_hcd,1603886625,19.578694,20000.00195056247000000
zig.divide_hcd_fast,953255998,11.636426,20000.00195056247000000
```

## Takeaways

- Zig HCD is competitive with HiGHS `HighsCDouble` on add operations.
- Zig multiplication is faster in this test, most likely because `twoProduct`
  can use hardware FMA while the current HiGHS header uses Dekker split.
- `divideHCD` remains the default stable path.
- `divideHCDFast` is the recommended hot-path option when one quotient
  correction is enough for the surrounding algorithm's error budget.
