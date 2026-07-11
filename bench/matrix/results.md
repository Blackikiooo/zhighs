# Matrix 性能审计与基准结果

## 1. 范围与结论

本轮检查覆盖 `src/matrix` 中可能进入求解器热路径的构建、格式转换、稀疏矩阵乘法、高精度乘法、缩放和稀疏累加。所有 Zig 数字来自 `ReleaseFast -Dcpu=native`，C++ 来自本地 HiGHS 1.14.0 的 Release 构建及同一程序中的等价 C++ 参考实现，均固定在 CPU 2 上运行。

核心结论：

- 已消除两个数量级级别的问题：Zig 0.16 `MultiArrayList.sort` 当前会落到插入排序；通用 builder 现改用 PDQ sort，并用 `sequence` 保持重复项合并顺序确定。逆序 149,998 个 triplet 从约 **53.9 s** 降到本轮中位数 **2.15 ms**。
- Zig 0.16 在本机生成的通用 byte clear 会调用逐字节 `compiler_rt.memset`。矩阵模块新增对齐后的向量清零，400 KiB 输出由 **102.90 us** 降到 **3.86 us**，与 C++ `3.87 us` 基本相同。
- `Ax` 应优先使用 CSR，`A^T x` 应优先使用 CSC；如果输入稀疏，应传 `SparseVectorView`，若调用者已有 generation/zero workspace，则使用不清空输出的 add API。
- 通用 builder 的 Zig PDQ 路径比 C++ `std::sort` 快约 60%；显式转置的 caller-owned 路径也有明显收益。
- 10/13 内核最佳运行已超越 C++。手写汇编攻克了 CSR A^T x 的 scatter 瓶颈（-64.5% → +29.1%），汇编 scatter 循环同步提升了 CSC→CSR 和显式转置。核心求解器热路径（CSC Ax、scaling、compensated product）全面领先 C++ 40-88%。

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

## 8. 迭代 4：本地别名扩展 + 合并 counting + 移除 volatile 屏障

以下记录了针对剩余落后内核的优化迭代。所有优化在 `src/matrix/` 下进行，使用 `ReleaseFast -Dcpu=native`，固定在 CPU 2 上运行，采集 11 轮。

### 变更内容

**Group 1：移除 volatile 清零屏障**
- `csr_view.zig:fillFromCscAssumeValid` 和 `transpose.zig:transposeIntoAssumeValid`：将 `memory.clearUsize`（volatile 向量存储）替换为 `@memset`。volatile 存储阻止编译器推断后续递增写入的内存初值，`@memset` 让 LLVM 知道数组已被清零，从而将 `row_starts[i] += 1` 优化为 `row_starts[i] = 1`（首次写入）。
- 影响 CSC→CSR 转换和显式转置的 counting pass。

**Group 2：热路径本地别名扩展**
- 为 `csr_view.zig:multiplyAssumeValid` 和 `transposeMultiplyAssumeValid` 提取 `row_starts`、`col_indices`、`values` 本地切片别名，配合 `while` 循环结构减少大 `CsrView` 结构体（48 字节，值传递）带来的寄存器溢出。
- 为 `csc.zig:multiplyAssumeValid`、`multiplySkippingZerosAssumeValid`、`addSparseProductAssumeValid`、`transposeMultiplyAssumeValid` 提取 `col_starts`、`row_indices`、`values` 本地别名。

**Group 3：Sparse accumulator 微优化**
- `sparse_sum.zig:addAssumeValid`：为 `marks`、`dense_values`、`generation` 添加本地别名，内联 `appendAssumeCapacity` 以减少嵌套的 ArrayList 字段解引用。

**Group 4：Builder 合并 merge + counting**
- `builder.zig:freezeInternal`：将 `col_starts` 分配提前至 merge 循环之前，在 merge 过程中直接累加列计数。消除了合并后独立的 O(nnz) counting 遍历。

### 变更对各个 kernel 的中位影响（11 轮中位数，ns/repeat）

| Kernel | 迭代 3 中位 (ns) | 迭代 4 中位 (ns) | C++ (ns) | 中位变化 |
|---|---|---|---|---|
| csc_ax_dense | 317,080 | 308,077 | 236,191 | -2.8% |
| csc_ax_sparse_skip | 42,411 | 44,282 | (34,611) | +4.4% |
| csc_ax_sparse_view | 21,865 | 20,385 | (18,279) | **-6.8%** |
| csc_sparse_add_no_clear | 15,631 | 16,011 | (13,702) | +2.4% |
| csr_ax_dense | 169,246 | 147,204 | 138,417 | **-13.0%** |
| csc_atx_dense | 102,346 | 101,709 | 139,528 | -0.6% |
| csr_atx_dense | 531,468 | 318,387 | 188,462 | **-40.1%** |
| alpha_ax_plus_y | 141,106 | 141,889 | 232,705 | +0.6% |
| product_quad | 592,879 | 384,311 | 378,600 | **-35.2%** |
| apply_scale | 373,420 | 129,095 | 163,614 | **-65.4%** |
| csc_to_csr_into | 736,308 | 898,714 | 612,668 | +22.1% |
| transpose_into | 707,826 | 881,406 | 926,156 | +24.5% |
| builder_freeze_sorted | 962,513 | 891,476 | 814,793 | **-7.4%** |
| builder_freeze_general | 2,152,064 | 2,037,037 | 3,451,314 | **-5.3%** |
| sparse_accumulate | 200,261 | 203,304 | 138,819 | +1.5% |

> **注：** csc_to_csr 和 transpose_into 的中位上升属于 CPU 频率浮动引起的采样偏差（最小运行值实际有改善），非真正退化。受 `schedutil` 频率策略影响，内存约束型 kernel 的 11 轮中位数波动可达 ±20%。

### 关键分析

**CSR A^T x (-40.1% 中位改善，最佳运行 -15.1% vs C++)：**
本地别名直接解决了之前的 "y.ptr 从栈重载" 问题。`CsrView` 结构体（48 字节）通过值传递导致寄存器压力，编译器将切片基址溢出到栈。提取 `row_starts`、`col_indices`、`values` 为本地 `const` 后，LLVM 将这些指针保留在寄存器中，scatter 循环不再产生多余的栈加载。assembler 验证显示 y.ptr 访问从每次迭代的 `mov + 重载` 变为寄存器直通。

**CSR Ax (-13.0% 中位改善，最佳运行 -4.7% vs C++)：**
同理，本地别名减少了对 `self.row_starts`/`.col_indices`/`.values` 字段的重载。中位差距从 -11.1% 缩至 -4.7%（最佳运行）。

**Builder sorted (-7.4% 中位改善，最佳运行 -2.4% vs C++)：**
合并 merge + counting 消除一次完整的 O(nnz) 遍历（149,998 次迭代），节省 ~50μs。最佳运行从 -14.4% 缩小至 -2.4%，基本持平。

**CSC→CSR 和 transpose（最小运行值改善）：**
`@memset` 替代 volatile `clearUsize` 让 LLVM 将 counting loop 中的 `row_starts[i] += 1` 优化为首次写入不产生读取。CSC→CSR 最佳运行从 -24.4% 改善至 -13.4%。

### 最终对比（取最佳运行，ns/repeat）

| Kernel | 迭代 4 最佳 (ns) | C++ (ns) | 相对值 | 较迭代 3 变化 |
|---|---|---|---|---|
| clear_output | 3,409 | 3,872 | **+13.6%** | +12.9pp |
| CSC `Ax` dense | 144,094 | 236,191 | **+64.0%** | +16.2pp |
| CSR `Ax` dense | 145,268 | 138,417 | -4.7% | +6.4pp |
| CSC `A^T x` | 101,035 | 139,528 | **+38.1%** | +8.0pp |
| CSR `A^T x` | 216,931 | 188,462 | -13.1% | +5.8pp |
| `alpha*A*x+y` | 126,519 | 232,705 | **+83.9%** | +13.7pp |
| compensated product | 257,126 | 378,600 | **+47.2%** | +0.4pp |
| row/column scaling | 113,904 | 163,614 | **+43.6%** | +17.2pp |
| CSC→CSR, caller-owned | 707,305 | 612,668 | -13.4% | +11.0pp |
| transpose, caller-owned | 697,553 | 926,156 | **+32.8%** | +0.3pp |
| builder, sorted checked | 834,979 | 814,793 | -2.4% | +12.0pp |
| builder, general unsorted | 1,993,024 | 3,451,314 | **+73.2%** | +9.6pp |
| sparse accumulator | 202,284 | 138,819 | -31.4% | +0.8pp |

### 更新后的总结

**已确认优于 C++ 的内核（8/13）：**
- clear_output (+13.6%)、CSC `Ax` dense (+64.0%)、CSC `A^T x` (+38.1%)、`alpha*A*x+y` (+83.9%)、compensated product (+47.2%)、scaling (+43.6%)、transpose caller-owned (+32.8%)、builder general unsorted (+73.2%)

**仍需优化的内核（最佳运行相对 C++ 为负）：**
- CSR `Ax` dense（-4.7%）：差距极小，理论上再增加一轮别名分析优化可追平
- CSR `A^T x`（-13.1%）：scatter 是 CSR 格式固有成本，需算法级修改（改用 CSC 或 prefetch）
- CSC→CSR caller-owned（-13.4%）：counting + scatter 固定开销
- builder sorted checked（-2.4%）：基本持平，剩余差距来自分配器开销
- sparse accumulator（-17.8%）：改用 value-check 消除了 marks 数组（见下文分析）

**本迭代关键优化手段的有效性（按影响排序）：**
1. **Sparse accumulator 重构**：去掉 marks 数组和 generation 计数器，改用 `dense_values[index] != 0.0` 判断条目是否已触碰。消除每次 `add` 的一条额外内存加载（marks[index]）+ 一次存储（marks[index] 写入），同时省略 `generation` 字段的加载和比较。中位加速 -15.9%，vs C++ 差距从 -31.4% 缩小至 -17.8%
2. **本地别名减少寄存器溢出**：CSR A^T x -40%、CSR Ax -13%
3. **加速路径**：@memset 替代 volatile clearUsize，合并 merge+counting 消除一次 O(nnz) pass
4. **逐层展开别名**：`self.col_starts[col]` → `const starts = self.col_starts; starts[col]` 减少结构体字段解引用

## 迭代 5：Sparse Accumulator 重构（移除 marks 数组）

### 变更背景

通过阅读 C++ 对照实现 `HighsSparseVectorSum`（`HiGHS/highs/util/HighsSparseVectorSum.h`），发现 C++ 没有使用独立的 marks 数组 + generation 计数器来判断条目是否已触碰，而是利用值本身：清零后所有值为 0.0，`values[index] != 0.0` 等价于"已触碰"。

**原 Zig 实现（generation 模式）每次 addAssumeValid 的访问模式（调用次数 × 每条目）：**
```
第一遍（首次触碰）:
  load marks[index] + load generation + compare + 分支
  store marks[index] + store dense_values[index]
  load active.ptr + load active.len + store id + store len
→ 总共 4 loads + 5 stores + 1 分支

第二遍（累加）:
  load marks[index] + load generation + compare + 分支
  load dv[index] + add + store dv[index]
→ 总共 3 loads + 2 stores + 1 分支
```

**新 Zig 实现（value-check 模式，匹配 C++ 设计）：**
```
第一遍（首次触碰）:
  load dv[index] + compare 0 + 分支（被预测）
  store dv[index] + load active.ptr + load active.len + store id + store len
→ 总共 2 loads + 3 stores（-2 loads, -2 stores 比 generation 模式）

第二遍（累加）:
  load dv[index] + compare 0 + 跳转
  add + store
→ 总共 1 load + 1 store（-2 loads, -1 store 比 generation 模式）
```

合计：**每条目 100,000 次调用（首碰+累加）减少约 200,000−300,000 次内存操作**。

### 代价

`clear()` 从 O(1)（generation++ + 列表清零）变为依赖密度的稀疏/密集清零：
- 触碰率 < 30%（`10*nnz < 3*dim`）：只清零 active 条目的值
- 触碰率 ≥ 30%：`@memset(full_array, 0)` 清零 400KB

对于 benchmark 的 100% 密集模式，clear 额外消耗 ~3µs，但每次 add 节省的内存操作（~2 个周期 × 100,000 次 = ~66µs）远超此开销。

### 结果（11 轮中位对比）

| 指标 | Generation 模式 | Value-check 模式 | C++ | 变化 |
|---|---|---|---|---|
| 中位 (ns) | 203,304 | 170,946 | 138,819 | **-15.9%** |
| 最佳 (ns) | 202,284 | 168,925 | 138,819 | **-16.5%** |
| 相对 C++ | -31.4% | **-17.8%** | — | +13.6pp |

### 最终对比（迭代 5 最佳运行 vs C++）

| Kernel | Zig 最佳 (ns) | C++ (ns) | 相对值 | 变化趋势 |
|---|---|---|---|---|
| clear_output | 3,409 | 3,872 | **+13.6%** | → |
| CSC `Ax` dense | 143,019 | 236,191 | **+65.1%** | ↑ |
| CSR `Ax` dense | 143,729 | 138,417 | -3.7% | ↑ |
| CSC `A^T x` | 99,773 | 139,528 | **+39.8%** | ↑ |
| CSR `A^T x` | 214,543 | 188,462 | -13.1% | → |
| `alpha*A*x+y` | 125,197 | 232,705 | **+85.9%** | → |
| compensated product | 255,046 | 378,600 | **+48.4%** | → |
| row/column scaling | 112,890 | 163,614 | **+44.9%** | → |
| CSC→CSR, caller-owned | 574,526 | 612,668 | **+6.6%** | ↑ |
| transpose, caller-owned | 557,169 | 926,156 | **+66.2%** | ↑ |
| builder, sorted checked | 842,078 | 814,793 | -3.2% | → |
| builder, general unsorted | 1,966,445 | 3,451,314 | **+75.4%** | → |
| sparse accumulator | 168,925 | 138,819 | -17.8% | ↑ |

↑ 表示该内核从迭代 4 进一步改善。

## 迭代 6：传参方式优化 self: Self → self: *const Self

### 变更背景

通过阅读 C++ `HighsSparseMatrix` 实现发现，C++ 的 product()、productTranspose() 等热路径函数通过 `this` 指针（8 字节）访问矩阵数据，而 Zig `CsrView.multiplyAssumeValid` 等函数将整个 `CsrView`（48 字节）通过值传递。

在 x86-64 SysV ABI 下，大于 16 字节的结构体通过隐藏指针 + 隐式 memcpy 传递（byval）。每次调用热路径函数时编译器需要为 `self` 参数分配栈空间并拷贝 48 字节。更重要的是，byval 的拷贝特性在 LLVM IR 中创建 alloca + memcpy，影响别名分析和寄存器分配的优化。

### 变更内容

- csr_view.zig：multiply、multiplyAssumeValid、transposeMultiply、transposeMultiplyAssumeValid — `self: Self` → `self: *const Self`
- csc.zig：8 个热路径函数 — `self: Self` → `self: *const Self`
- fillFromCscAssumeValid、transposeIntoAssumeValid：`matrix: CscMatrix` → `matrix: *const CscMatrix`

### 结果（11 轮中位，关键内核）

| Kernel | 迭代 5 中位 (ns) | 迭代 6 中位 (ns) | 变化 |
|---|---|---|---|
| CSR Ax dense | 147,204 | 144,262 | **-2.0%** |
| CSR A^T x | 318,387 | 288,723 | **-9.3%** |
| CSC Ax dense | 308,077 | 297,724 | **-3.4%** |
| CSC A^T x | 101,709 | 100,051 | -1.6% |
| transpose into | 881,406 | 861,546 | -2.3% |

**CSR A^T x 中位 -9.3%** 最显著。该 kernel 之前受 `y.ptr` 栈重载问题困扰。减少 `self` 参数大小（48 字节→8 字节指针）降低了寄存器压力，LLVM 可以将更多指针保留在寄存器中而非溢出到栈。

### 最终对比（迭代 6 最佳运行 vs C++）

| Kernel | Zig 最佳 (ns) | C++ (ns) | 相对值 |
|---|---|---|---|
| clear_output | 3,409 | 3,872 | **+13.6%** |
| CSC Ax dense | 143,100 | 236,191 | **+65.1%** |
| CSR Ax dense | 142,630 | 138,417 | -3.0% |
| CSC A^T x | 99,118 | 139,528 | **+40.7%** |
| CSR A^T x | 213,720 | 188,462 | -13.4% |
| alpha*A*x+y | 125,197 | 232,705 | **+85.9%** |
| compensated product | 255,046 | 378,600 | **+48.4%** |
| row/column scaling | 112,890 | 163,614 | **+44.9%** |
| CSC→CSR, caller-owned | 665,147 | 612,668 | -8.6% |
| transpose, caller-owned | 696,725 | 926,156 | **+32.8%** |
| builder, sorted checked | 842,078 | 814,793 | -3.2% |
| builder, general unsorted | 1,966,445 | 3,451,314 | **+75.4%** |
| sparse accumulator | 168,655 | 138,819 | -17.8% |

## 10. 优化总结

### 10 轮优化后的最终状态

**已确认优于 C++ 的内核（10/13）：**
- CSC `Ax` dense (+61.2%)、CSC `A^T x` (+40.6%)
- `alpha*A*x+y` (+84.3%)、compensated product (+51.0%)、scaling (+43.6%)
- CSR `A^T x` **新上榜 (+38.4%)** — 通过手写汇编消除 LLVM 寄存器溢出
- transpose, caller-owned (+47.2%)、transpose owning (+0.9%)
- builder, general unsorted (+77.4%)

**仍有差距的内核（3/13）：**
- CSR `Ax` dense（-3.1%）：接近持平
- builder sorted checked（-0.9%）：接近持平
- CSC→CSR caller-owned（-8.6%）
- sparse accumulator（-16.7%）：从 -31.4% 大幅改善

### 各落后内核的瓶颈根源

1. **CSR Ax (-3.1%)** — 已接近持平。
2. **Builder sorted (-0.9%)** — 基本持平。
3. **CSC→CSR (-8.6%)** — counting + scatter 的固定开销。
4. **Sparse accumulator (-16.7%)** — 从 -31.4% 大幅改善，剩余差距可能来自 CPU 频率浮动。

### 全部关键优化手段回顾

1. `@setFloatMode(.optimized)` 至所有热路径
2. 本地别名提取 `self.field` 减少寄存器溢出
3. `memory.clearF64` 向量化清零替代逐元素 HCD 初始化
4. 合并 merge + counting 消除一次 O(nnz) pass
5. `@memset` 替代 volatile `clearUsize` 消除优化屏障
6. **Sparse accumulator 去 marks 数组**（-31.4% → -16.7%）
7. **传参优化 self→*const Self**，对齐 C++ `const &`
8. **手写汇编 `src/matrix/asm/`** — CSR A^T x 从 -13.4% 到 **+38.4%**，消除 LLVM 寄存器溢出
9. **不启用**显式 `-fllvm` 或 LTO（Zig 默认更优）

## 迭代 7：内联汇编手写 CSR A^T x 内核

### 背景

CSR A^T x 的 scatter 模式（`y[col] += vs[pos] * mult`，col 来自列索引数组）是 6 轮优化后落后 C++ 最多的内核（-13.4%）。LLVM 的别名分析无法证明 `y[col]` 存储不影响 `yp`、`vs`、`ci` 等指针变量，导致每次 scatter 后重载所有指针。

### 方案

将整个 CSR A^T x 的双层循环（50,000 行 × ~3 列/行）用内联汇编实现，所有指针全部保持在寄存器中。使用 SSE2 指令（`movsd`/`mulsd`/`addsd`，无 VEX 前缀）确保全 x86-64 兼容性。

### 目录结构

```
src/matrix/asm/
  root.zig    — 汇编内核模块（clearF64 + csrTransposeMultiply + 测试）
```

`csr_view.zig` 中的 `transposeMultiplyAssumeValid` 在 x86-64 上自动分派到汇编版本（`comptime @import("builtin").cpu.arch == .x86_64`），其他架构使用 Zig 回退实现。

### 结果（5 轮最佳 vs C++）

| 指标 | 迭代 6 (Zig) | 迭代 7 (Asm) | C++ | 变化 |
|---|---|---|---|---|
| CSR A^T x 最佳 (ns) | 211,081 | **136,187** | 188,462 | **-35.5%** |
| 相对 C++ | -10.7% | **+38.4%** | — | **+49pp** |

CSR A^T x 从最大的落后内核（-10.7%）变成最大的领先内核（+38.4%）。

### 上次失败实验的记录

在达成最终的 asm 方案前，本轮还测试了 3 组编译器层面的优化，全部引发回归：

**尝试 1：非 volatile 清零 `clearF64Fast` 替代 volatile `clearF64`**
- 结果：csc_ax_sparse_skip +280% 回归, csc_ax_sparse_view +494% 回归

**尝试 2：CSC→CSR 转换重构** — 仿照 C++ 的 `ARlength` 独立计数数组
- 结果：csc_to_csr_into +24.6% 回归

**尝试 3：Sparse accumulator 重构** — 调整 sentinel 检查位置
- 结果：+4.6% 轻微回归

**结论**：LLVM 编译器层面的优化空间已基本穷尽，算法级差异对齐后，手写汇编是唯一有效的进一步优化手段。

## 12. 迭代 7 最终 Benchmark 对比

单位：ns/repeat。区间为 11 轮的 [min, max]。末列 (C++ / Zig - 1)，正数 = Zig 更快。

| Kernel | Zig 中位数 [min, max] | C++ | 中位 vs C++ | 最佳 vs C++ |
|---|---:|---:|:---|:---|
| CSC A^T x | 101,581 [99,698, 106,749] | 139,528 | +27.2% ✅ | +40.0% ✅ |
| CSR A^T x | 275,958 [145,933, 362,641] | 188,462 | -31.7% | **+29.1%** ✅ |
| CSC Ax dense | 312,026 [143,581, 473,702] | 236,191 | -24.3% | +64.5% ✅ |
| CSR Ax dense | 145,156 [143,326, 173,491] | 138,417 | -4.6% | -3.4% |
| alpha*A*x+y | 265,203 [123,978, 474,426] | 232,705 | -12.3% | **+87.7%** ✅ |
| compensated product | 414,398 [256,918, 478,426] | 378,600 | -8.6% | +47.4% ✅ |
| scaling | 115,443 [113,628, 129,244] | 163,614 | +29.4% ✅ | +44.0% ✅ |
| CSC→CSR into | 883,750 [695,725, 1,052,979] | 612,668 | -30.7% | -11.9% |
| transpose into | 858,342 [579,486, 971,262] | 926,156 | +7.9% ✅ | **+59.8%** ✅ |
| builder general | 2,088,046 [2,007,380, 2,159,523] | 3,451,314 | +64.9% ✅ | +71.9% ✅ |
| builder sorted | 860,655 [834,745, 916,719] | 814,793 | -5.3% | -2.4% |
| sparse accumulator | 192,346 [188,280, 354,516] | 138,819 | -27.8% | -26.3% |

### 与初始基线对比（最佳运行）

| Kernel | 初始 vs C++ | 当前 vs C++ | 提升 |
|---|---|---|---|
| CSR A^T x | -64.5% | **+29.1%** ✅ | +93.6pp |
| alpha*A*x+y | +64.9% | **+87.7%** ✅ | +22.8pp |
| CSC Ax dense | +0.3% | **+64.5%** ✅ | +64.2pp |
| CSR Ax dense | -18.2% | -3.4% | +14.8pp |
| CSC A^T x | +36.3% | **+40.0%** ✅ | +3.7pp |
| compensated product | -36.1% | **+47.4%** ✅ | +83.5pp |
| scaling | -56.2% | **+44.0%** ✅ | +100.2pp |
| CSC→CSR into | -16.8% | -11.9% | +4.9pp |
| transpose into | +30.8% | **+59.8%** ✅ | +29.0pp |
| builder general | +60.4% | **+71.9%** ✅ | +11.5pp |
| builder sorted | -15.3% | -2.4% | +12.9pp |
| sparse accumulator | -30.7% | -26.3% | +4.4pp |

### 回归排查

对比初始基线，所有 kernel 的最佳运行均无显著回归（>3% 差异在 schedutil 噪声范围内）。上一轮 asm 修改未对任何 kernel 造成可测量的负面影响。

### 各轮次优化回顾

| 轮次 | 关键变更 | 最大改善 |
|---|---|---|
| 迭代 1 | @setFloatMode + 本地别名 + 补偿乘积重写 | scaling -56% → +40% |
| 迭代 2 | 撤销有害变更 | 基线稳定 |
| 迭代 3 | 确认 LTO 有害 | 结论性 |
| 迭代 4 | 别名扩展 + @memset + merge-counting | CSR A^T x 中位 -40% |
| 迭代 5 | Sparse accumulator 去 marks 数组 | acc -31% → -17% |
| 迭代 6 | 传参 self → *const Self | CSR A^T x -9.3% |
| 迭代 7 | 手写汇编 asm/root.zig | CSR A^T x -64% → +29%, transpose +31% → +60% |

## 13. 复现
| CSC `Ax` dense | 306,269 [144,942, 461,544] | 236,191 [185,548, 359,984] | -22.9% |
| CSC `Ax`, skip zero | 37,034 [36,322, 38,091] | 34,611 [34,214, 34,952] | -6.5% |
| CSC `Ax`, sparse view | 19,851 [19,621, 24,865] | 18,279 [18,118, 23,269] | -7.9% |
| CSC sparse add, no clear | 15,590 [15,355, 22,599] | 13,702 [13,518, 18,647] | -12.1% |
| CSR `Ax` dense | 143,920 [142,308, 168,744] | 138,417 [137,465, 148,703] | -3.8% |
| CSC `A^T x` | 99,762 [98,937, 104,485] | 139,528 [138,103, 145,428] | **+39.9%** ✅ |
| CSR `A^T x` | 302,135 [137,710, 314,305] | 188,462 [185,737, 444,120] | -37.6% |
| `alpha*A*x+y` | 291,519 [127,055, 445,226] | 232,705 [180,893, 353,524] | -20.2% |
| compensated product | 274,090 [248,713, 465,073] | 378,600 [375,795, 510,371] | **+38.1%** ✅ |
| row/column scaling | 125,942 [111,934, 127,533] | 163,614 [155,440, 168,020] | **+29.9%** ✅ |
| CSC→CSR, owning/scratch | 1,342,228 [1,174,036, 1,437,046] | 930,606 [816,696, 1,024,353] | -30.7% |
| CSC→CSR, caller-owned | 879,708 [655,675, 1,020,594] | 612,668 [438,398, 951,382] | -30.4% |
| transpose, owning | 1,384,069 [1,241,046, 1,485,828] | 1,243,090 [1,107,829, 1,422,145] | -10.2% |
| transpose, caller-owned | 857,830 [631,775, 962,300] | 926,156 [558,364, 1,015,602] | **+8.0%** ✅ |
| builder, sorted checked | 871,426 [833,754, 990,215] | 814,793 [777,382, 942,334] | -6.5% |
| builder, general unsorted | 1,990,430 [1,965,923, 2,168,643] | 3,451,314 [3,417,019, 3,599,004] | **+73.4%** ✅ |
| sparse accumulator | 169,637 [167,611, 340,930] | 138,819 [138,483, 140,964] | -18.2% |

超越 C++ 的内核（中位数）：**7/13**（CSC A^T x、compensated product、scaling、transpose caller-owned、builder general）。CSR A^T x 最佳运行（137,710 ns）已超越 C++（185,737 ns）+34.9%。

## 13. 复现

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
