# Matrix 性能审计与基准结果

## 1. 范围与结论

本轮检查覆盖 `src/matrix` 中可能进入求解器热路径的构建、格式转换、稀疏矩阵乘法、高精度乘法、缩放和稀疏累加。所有 Zig 数字来自 `ReleaseFast -Dcpu=native`，C++ 来自本地 HiGHS 1.14.0 的 Release 构建及同一程序中的等价 C++ 参考实现，均固定在 CPU 2 上运行。

核心结论：

- 已消除两个数量级级别的问题：Zig 0.16 `MultiArrayList.sort` 当前会落到插入排序；通用 builder 现改用 PDQ sort，并用 `sequence` 保持重复项合并顺序确定。逆序 149,998 个 triplet 从约 **53.9 s** 降到本轮中位数 **2.15 ms**。
- Zig 0.16 在本机生成的通用 byte clear 会调用逐字节 `compiler_rt.memset`。矩阵模块新增对齐后的向量清零，400 KiB 输出由 **102.90 us** 降到 **3.86 us**，与 C++ `3.87 us` 基本相同。
- `Ax` 应优先使用 CSR，`A^T x` 应优先使用 CSC；如果输入稀疏，应传 `SparseVectorView`，若调用者已有 generation/zero workspace，则使用不清空输出的 add API。
- 通用 builder 的 Zig PDQ 路径比 C++ `std::sort` 快约 60%；显式转置的 caller-owned 路径也有明显收益。
- 当前需继续优化的热点是 compensated product、scaling、稀疏累加和 CSC→CSR。部分内核出现明显双峰，说明 `schedutil` governor 与睿频会显著影响单次结果，因此本报告同时保留最小/最大值，不把一个百分比当作跨机器结论。

## 2. 测试环境

| 项目 | 内容 |
|---|---|
| 日期 | 2026-07-11 16:25（Asia/Shanghai） |
| OS / Kernel | Linux 6.17.0-35-generic x86_64 |
| CPU | AMD Ryzen 5 3500X，6 核 6 线程，boost 开启 |
| 频率策略 | `schedutil`，2200–4120 MHz；固定 CPU 2，不固定实际频率 |
| Cache | L1d 192 KiB，L1i 192 KiB，L2 3 MiB，L3 32 MiB |
| 内存 | 15 GiB RAM，4 GiB swap |
| Zig | 0.16.0，`ReleaseFast -Dcpu=native` |
| C++ | g++ 13.3.0，`-O3 -march=native -DNDEBUG -flto` |
| HiGHS | 1.14.0，git `dcc25308d8`，Release、FAST_BUILD、IPO/LTO，32-bit HighsInt |
| zhighs index | 默认 `HInt=w32`，与本次 HiGHS 构建一致 |

## 3. 公平性与数据集

- 矩阵为 `50000 x 50000` 的三对角 CSC，`nnz=149998`，对角值 4、相邻值 -1；Zig/C++ 使用相同的解析生成方式。
- 普通乘法重复 200 次，高精度乘法和 accumulator 重复 20 次，转换/转置重复 10 次，builder 每个样本构建一次。
- 共采集 11 轮，奇偶轮交替 Zig/C++ 的执行先后，以中位数为主，同时报告最小值和最大值。
- 每个 kernel 有 compiler barrier，并在计时循环后读取完整输出 checksum，避免编译器删除计算。checksum 在两端一致。
- HiGHS 有直接对应实现时调用 HiGHS；HiGHS 没有公开对应组件时使用同一 C++ harness 中的简单 reference、`std::sort` 或预分配数组实现。
- `barrier_only` 与 `clear_output_bytes` 是诊断项，不参与语言性能结论。

该合成矩阵能稳定复现访问模式，但不能替代大型、非规则工业矩阵；后者仍保留在 `todo.md`。

## 4. 基准结果

单位为 `ns / repeat`。区间是 11 轮的 `[min, max]`。最后一列为 `(C++ / Zig - 1)`：正数表示 Zig 更快，负数表示 Zig 更慢。双峰或宽区间的行只能用于定位热点，不能视为稳定排名。

| Kernel | Zig 中位数 [min, max] | C++/HiGHS 中位数 [min, max] | Zig 相对值 |
|---|---:|---:|---:|
| output clear | 3,860 [3,604, 3,964] | 3,872 [3,698, 4,128] | +0.3% |
| CSC `Ax` dense | 317,080 [150,707, 441,529] | 236,191 [185,548, 359,984] | -25.5% |
| CSC `Ax`, skip zero | 42,411 [41,581, 44,771] | 34,611 [34,214, 34,952] | -18.4% |
| CSC `Ax`, sparse view | 21,865 [21,772, 30,935] | 18,279 [18,118, 23,269] | -16.4% |
| CSC sparse add, no clear | 15,631 [15,498, 25,061] | 13,702 [13,518, 18,647] | -12.3% |
| CSR `Ax` dense | 169,246 [142,101, 170,368] | 138,417 [137,465, 148,703] | -18.2% |
| CSC `A^T x` | 102,346 [98,769, 103,528] | 139,528 [138,103, 145,428] | +36.3% |
| CSR `A^T x` | 531,468 [216,496, 579,595] | 188,462 [185,737, 444,120] | -64.5% |
| `alpha*A*x+y` | 141,106 [129,485, 312,182] | 232,705 [180,893, 353,524] | +64.9% |
| compensated product / HiGHS productQuad | 592,879 [494,266, 640,156] | 378,600 [375,795, 510,371] | -36.1% |
| row/column scaling | 373,420 [167,910, 399,499] | 163,614 [155,440, 168,020] | -56.2% |
| CSC→CSR, owning/scratch | 1,209,025 [1,077,009, 1,374,172] | 930,606 [816,696, 1,024,353] | -23.0% |
| CSC→CSR, caller-owned | 736,308 [612,213, 972,051] | 612,668 [438,398, 951,382] | -16.8% |
| transpose, owning | 1,320,980 [1,275,166, 1,430,061] | 1,243,090 [1,107,829, 1,422,145] | -5.9% |
| transpose, caller-owned | 707,826 [421,671, 718,116] | 926,156 [558,364, 1,015,602] | +30.8% |
| builder, sorted checked | 962,513 [896,878, 1,411,302] | 814,793 [777,382, 942,334] | -15.3% |
| builder, general unsorted | 2,152,064 [1,988,615, 2,696,545] | 3,451,314 [3,417,019, 3,599,004] | +60.4% |
| sparse accumulator | 200,261 [199,593, 211,689] | 138,819 [138,483, 140,964] | -30.7% |

附加 Zig 精度取舍结果：exact robust 为 806,569 ns，exact fast 为 594,230 ns，HiGHS 语义对应的 compensated product 为 592,879 ns。exact fast 仍计算 HCD 乘积，不能与只对 f64 乘积做 compensated sum 的 HiGHS `productQuad` 当作完全相同算法比较。

## 5. 热路径 API 取舍

| 场景 | 安全默认 API | 显式高性能 API | 调用者承担的条件 |
|---|---|---|---|
| CSC dense `Ax` | `multiply` | `multiplyAssumeValid` | 切片长度、矩阵规范已验证 |
| 大量零的 dense x | `multiplySkippingZeros` | `multiplySkippingZerosAssumeValid` | 同上；仍扫描全部 x |
| 真正稀疏的 x | `multiplySparse` | `multiplySparseAssumeValid` | sparse view 已规范化 |
| 已有清零/generation workspace | `addSparseProduct` | `addSparseProductAssumeValid` | y 不会自动清空，语义是累加 |
| 行方向 `Ax` | `CsrView.multiply` | `CsrView.multiplyAssumeValid` | cache revision 和尺寸有效 |
| 高精度乘法 | `multiplyHighPrecision` | `...AssumeValid` / `...FastAssumeValid` / `multiplyCompensatedAssumeValid` | 调用者明确选择 exact robust、exact fast 或 HiGHS-compatible 语义 |
| Triplet 冻结 | `freeze` | `freezeSorted` / `freezeSortedAssumeValid` | 后两者要求输入已按 `(col,row)` 排序；AssumeValid 不再扫描检查 |
| CSC→CSR | `CsrCache.build` | `buildWithScratchAssumeValid` / `fillFromCscAssumeValid` | scratch 和输出数组尺寸正确；`fill` 不分配 |
| 转置 | `transpose` | `transposeIntoAssumeValid` | 调用者复用 starts/rows/values/cursor |
| 稀疏累加 | checked add | `reserve` + `addAssumeValid` | 预留容量，Id 在范围内，避免循环内扩容/检查 |

`AssumeValid` 不代表“关闭所有安全性”：它只移除函数文档中列出的重复验证。分配失败、容量不足和不能由编译期证明的结构约束仍需由返回类型或调用者保证。

## 6. 后续优化优先级

1. 用 `performance` governor 或固定频率复测双峰项目，并增加 `perf stat` 的 cycles、instructions、branch-misses、L1/LLC misses；当前结果只能确认热点，不能精确归因。
2. 对 compensated product 检查 HCD 累加的内联和寄存器生命周期；保持 exact 与 HiGHS-compatible 两套语义，不用精度换取默认速度。
3. 对 scaling 拆分 row/column benchmark，确认慢点来自 scatter 写、检查还是本机降频。
4. 优化 `SparseAccumulator` 的 generation/index/value 三数组访问，并针对 cut/presolve 的真实稀疏度分布测试。
5. 为 CSC→CSR/transpose 增加可复用、同宽 offset storage；在证据充分前不引入 SIMD prefetch 或并行化。
6. 增加大于 L3 的随机和真实 LP 矩阵，记录吞吐量、峰值内存与 cache 行为。

## 7. 复现

脚本会重新构建 HiGHS、Zig benchmark 和 C++ harness，并把每轮原始 CSV 写到临时目录：

```bash
HIGHS_SOURCE=/home/godv/documents/codefiles/cppfiles/HiGHS \
CPU_CORE=2 RUNS=11 \
bash bench/matrix/run_comparison.sh
```

默认原始结果目录是 `/tmp/zhighs-matrix-comparison/results`。Zig benchmark 当前把 CSV 写到 stderr，脚本已按此行为重定向；直接运行时可使用：

```bash
zig build bench-matrix -Doptimize=ReleaseFast -Dcpu=native 2> zig-matrix.csv
```
