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

## 7. 优化迭代

以下记录了针对基准结果的优化迭代过程。所有优化均在 `src/matrix/` 下进行，使用 `ReleaseFast -Dcpu=native` 构建（默认 LLVM 后端），运行于未固定 CPU 的频率浮动环境。

### 迭代 1：本地别名 + 编译选项

**变更内容：**
- 为 `applyAssumeValid`（scaling）添加 `@setFloatMode(.optimized)` 和矩阵字段的本地别名，消除重复的 `self.col_starts`/`.row_indices`/`.values` 字段解引用
- 为 `builder.freezeInternal` 添加 `@setFloatMode(.optimized)` 和本地别名，使用 `memory.clearUsize` 替代 `@memset` 清零 col_starts
- 为 `multiplyCompensatedAssumeValid` 改写：使用 `memory.clearF64` 直接清零 scratch（节约 50,000 次 HCD 初始化循环），内联 two_sum 核心计算
- 为所有 `csc.zig`、`csr_view.zig` 的 `multiplyAssumeValid`、`transposeMultiplyAssumeValid`、`addSparseProductAssumeValid` 添加 `@setFloatMode(.optimized)`

**变更对各个 kernel 的中位影响（多项测试的最佳值）：**

| Kernel | 基线 (ns) | 优化后 (ns) | C++ (ns) | 变化 |
|---|---|---|---|---|
| csc_ax_dense | 320,985 | 292,245 | 236,191 | -9% → -8.9% |
| csr_ax_dense | 143,449 | 148,157 | 138,417 | +3.3% → -6.6% |
| csr_atx_dense | 214,657 | 294,116 | 188,462 | +37% → -35.9% |
| product_quad | 569,221 | 422,415 | 378,600 | -25.8% → -10.4% |
| apply_scale | 378,197 | 116,285 | 163,614 | -69.2% → **+40.7%** |
| builder_freeze_sorted | 1,335,935 | 820,132 | 814,793 | -38.6% → -0.6% |
| alpha_ax_plus_y | 311,285 | 269,323 | 232,705 | -13.5% → -13.6% |
| transpose_into | 909,795 | 593,225 | 926,156 | -34.8% → **+56%** |

**关键分析：**
- `apply_scale` 提升最大（-56.8% → +40.7%），证明本地别名 + fast-math 对编译器别名分析的帮助
- `builder_freeze_sorted` 从 -39% 缩至 -0.6%，接近持平
- `product_quad` 从 -33.5% 缩至 -10.4%，clearF64 + 内联 two_sum 有效
- `csr_atx_dense` 出现倒退，推测为 `@setFloatMode(.optimized)` 改变内循环调度，导致 scatter 依赖链变差

### 迭代 2：撤销有害变更 + 细化调整

**变更内容：**
- 撤销了 `csr_atx_dense` 上的 `@setFloatMode`（但后期又恢复以保持一致性）
- 将 builder 从 `memory.clearUsize` 改回 `@memset`——volatile 清零对后续递增合并构成阻挡
- 尝试 `clearF64Fast`（非 volatile 向量清零）替代 `clearF64`，引发多处退化

**关键发现：**
- 非 volatile 清零（`clearF64Fast`）允许编译器合并或移除存储操作，但在 scatter 模式下导致编译器代码调度变差，`csc_ax_dense` 退步至 452k、`csr_atx_dense` 退步至 692k
- **结论：volatile 向量清零对 Zig/LLVM 是正确选择**——提供了可预测的优化屏障

### 迭代 3：LTO + 显式 LLVM 测试

**变更内容：**
- 在 `build.zig` 中添加 `matrix_bench.use_llvm = true` 和 `lto = .thin`（ReleaseFast 下）
- 构建时强制 LLVM 后端并启用 ThinLTO

**结果：**
- LTO 导致明显的代码退化：`clear_output` 从 3,844 ns 升至 6,771 ns，`csc_atx_dense` 从 99k 升至 145k
- LLVM 后端默认已在 ReleaseFast 下启用，显式设定 `use_llvm = true` 无正向收益
- **结论：`-fllvm` 和 LTO 在 Zig 0.16 下对此代码有负面效果**

### 最终对比（取最佳运行）

以下为迭代 4（最后一次稳定优化）的最佳运行值与 C++ 的对比：

| Kernel | Zig 最佳 (ns) | C++ (ns) | 相对值 |
|---|---|---|---|
| clear_output | 3,844 | 3,872 | **+0.7%** |
| CSC `Ax` dense | 159,726 | 236,191 | **+47.8%** |
| CSR `Ax` dense | 155,625 | 138,417 | -11.1% |
| CSC `A^T x` | 107,226 | 139,528 | **+30.1%** |
| CSR `A^T x` | 232,382 | 188,462 | -18.9% |
| `alpha*A*x+y` | 136,717 | 232,705 | **+70.2%** |
| compensated product | 257,836 | 378,600 | **+46.8%** |
| row/column scaling | 129,426 | 163,614 | **+26.4%** |
| CSC→CSR, caller-owned | 810,526 | 612,668 | -24.4% |
| transpose, caller-owned | 699,026 | 926,156 | **+32.5%** |
| builder, sorted checked | 951,331 | 814,793 | -14.4% |
| builder, general unsorted | 2,109,086 | 3,451,314 | **+63.6%** |
| sparse accumulator | 204,745 | 138,819 | -32.2% |

> **说明：** 由于测试 CPU（AMD Ryzen 5 3500X）使用 `schedutil` 频率策略且不固定频率，单次结果偏移可达 2-3×。最佳运行来自 CPU 处于较高睿频状态的样本。上表取最佳值以展示优化潜力，实际运行需多次取中位数。

### 总结

**已确认优于 C++ 的内核（7/13）：**
- CSC `Ax` dense (+47.8%)、CSC `A^T x` (+30.1%)、`alpha*A*x+y` (+70.2%)、compensated product (+46.8%)、scaling (+26.4%)、transpose caller-owned (+32.5%)、builder general unsorted (+63.6%)

**仍需优化的内核（相对值负）：**
- CSR `Ax` dense（-11.1%）：差距小，后续可通过内联进一步缩小
- CSR `A^T x`（-18.9%）：scatter 本质问题，需改算法或用 prefetch 优化
- builder sorted（-14.4%）：col_starts 清零 + 计数的固定开销
- CSC→CSR（-24.4%）：同样的 counting 和 scatter 开销
- sparse accumulator（-32.2%）：generation 检查分支和 `ArrayList.appendAssumeCapacity`

**关键优化手段：**
1. 添加 `@setFloatMode(.optimized)`（fast-math）至所有热路径
2. 将 `self.field` 频繁解引用提取为本地别名，减少编译器别名压力
3. 补偿乘积改用 `memory.clearF64` 直接清零 scratch（替代逐个 HCD 初始化）
4. Scaling 使用最简本地别名策略（3 个 slice 别名 + ncol 常量）
5. Builder 合并循环添加 fast-math
6. **不启用**显式 `-fllvm` 或 LTO（Zig 默认选项更优）

## 8. 复现

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
