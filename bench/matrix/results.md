# Matrix 性能审计与基准结果

## 1. 变更摘要 (2026-07-12)

### 代码修复
| 修复项 | 说明 |
|---|---|
| `computeLayout(N, sizes, aligns)` | 统一 checked layout helper，覆盖 alignForward overflow。替换全部 lean 路径 |
| `computePageColoredLayout(N, sizes, page_offsets)` | 统一 checked page-colored layout helper，替换全部 6 处 compact 路径 |
| `freezeFromCanonicalArraysAssumeValid` | 零 merge 的最简 count+memcpy 路径，为 reusable buffer 提供基线 |
| `freezeFromSortedArraysAssumeValid` two-pass | 消除 temp `out_rows`/`out_values` 数组，pass 2 直接写输出 |
| Cursor 泄漏 | `transposeLeanAssumeValid` cursor 改用独立 alloc/free |
| 对齐修复 | 所有 lean 路径使用 `computeLayout` 确保自然对齐 |

### Benchmark 修复
| 修复项 | 说明 |
|---|---|
| SHA-256 记录 | 每次运行记录 Zig/C++ binary hash 到 raw.csv |
| 移除 `\|\| true` | benchmark 崩溃直接 abort |
| 逐轮跨实现校验 | checksum + struct hash 在每轮 Zig vs C++ 之间验证 |
| `clear_output` 100k repeats | kernel 时间从 2ms→360ms，消除 init 噪声 |
| `builder_freeze_canonical` kernel | 最简 freeze 路径，用于 allocator 隔离实验 |
| `builder_freeze_reusable` kernel | 预分配 + 预触页 reusable buffer，模拟求解器实际使用场景 |
| sparse_accumulator `ZHIGHS_PERF_DIAG` | 输出页偏移 + 3 batch + asm barrier 对照 |

## 2. 测试环境

| 项目 | 内容 |
|---|---|
| 日期 | 2026-07-12 |
| OS / Kernel | Linux 6.17.0-35-generic x86_64 |
| CPU | AMD Ryzen 5 3500X (Zen2), 6 核 6 线程 |
| Zig | 0.16.0, `ReleaseFast -Dcpu=native` |
| C++ | g++ 13.3.0, `-O3 -march=native -DNDEBUG -flto` |
| HiGHS | 1.14.0, git `dcc25308d8`, 32-bit HighsInt |
| perf_event_paranoid | 4 (需 sudo 访问硬件计数器) |

### 2026-07-12 ReleaseFast 复测快照

运行条件：CPU 2，11 轮 Zig/C++ 交错执行；Zig 使用
`ReleaseFast -Dcpu=native`，C++ 使用
`-O3 -march=native -DNDEBUG -flto`。运行时 load average 为 2.74，CPU
governor 为 `schedutil`，因此 MAD 超过 15% 的项目只记为
`INCONCLUSIVE`，不能据此判断回归或领先。21 个 kernel 的 checksum 与
struct hash 均匹配。

| Kernel | Zig 中位数 (ns) | C++ 中位数 (ns) | Zig 相对性能 | 判定 |
|---|---:|---:|---:|---|
| clear_output | 3,524 | 3,628 | +2.9% | 可信 |
| csc_ax_dense | 294,732 | 189,204 | -55.8% | 可信，但 Zig 存在宽幅分布 |
| csc_ax_sparse_skip | 24,954 | 35,600 | +29.9% | 可信 |
| csc_ax_sparse_view | 17,683 | 18,281 | +3.3% | 可信 |
| csc_sparse_add_no_clear | 12,897 | 13,554 | +4.9% | 可信 |
| csr_ax_dense | 171,395 | 140,355 | -22.1% | 可信 |
| csc_atx_dense | 99,038 | 145,869 | +32.1% | 可信 |
| csr_atx_dense | 197,581 | 263,577 | +25.0% | INCONCLUSIVE（C++ MAD 28.2%） |
| alpha_ax_plus_y | 296,582 | 232,679 | -27.5% | INCONCLUSIVE（两端 MAD 超标） |
| product_quad | 237,902 | 368,869 | +35.5% | 可信 |
| apply_scale | 132,700 | 163,820 | +19.0% | 可信 |
| csc_to_csr_into | 746,100 | 705,088 | -5.8% | 临界噪声区（MAD 约 12–14%） |
| csc_to_csr_owning | 1,578,431 | 1,460,319 | -8.1% | 临界噪声区 |
| transpose_into | 785,297 | 614,860 | -27.7% | INCONCLUSIVE（两端 MAD 超标） |
| transpose_owning | 1,392,001 | 672,367 | -107.0% | INCONCLUSIVE（C++ MAD 18.5%） |
| builder_freeze_sorted | 2,005,578 | 1,825,242 | -9.9% | 可信 |
| builder_freeze_prepopulated | 1,357,088 | 819,607 | -65.6% | 可信 |
| builder_freeze_canonical | 1,005,969 | 467,159 | -115.3% | 可信；owning 分配路径 |
| builder_freeze_reusable | 316,665 | 548,791 | +42.3% | 可信；推荐生产热路径 |
| builder_freeze_general | 3,232,957 | 4,416,613 | +26.8% | 可信 |
| sparse_accumulate | 300,921 | 273,111 | -10.2% | 已知双峰，INCONCLUSIVE |

原始数据位于 `/tmp/zhighs-matrix-retest/results/raw.csv`，汇总数据位于
`/tmp/zhighs-matrix-retest/results/summary.csv`。`/tmp` 文件不是版本化产物；
本节保留了可追溯的中位数和测量条件。

## 3. 实验结论汇总

### CONFIRMED — 根因已定位

| Kernel | 发现 | 证据 |
|---|---|---|
| **builder_freeze_canonical (owning)** | 瓶颈在分配器页 fault，非 merge 算法 | reusable 消除 161× page faults；本轮 11 轮交错复测中 Zig reusable 为 316,665 ns，C++ 为 548,791 ns，领先约 42% |
| **clear_output** | Zig SIMD volatile 比 C++ `std::fill` 少 9% cycles | 100k repeats: Zig 1,359M cycles (IPC 3.70) vs C++ 1,486M cycles (IPC 1.50) |
| **csc_to_csr_owning** | 不是真实回归 | 7 次独立 perf: cycles 607-665M (±5%)，C++ 612M。wall-time 差距来自测量噪声 |

### HYPOTHESIS — 假设，待进一步验证

| Kernel | 假设 | 待验证 |
|---|---|---|
| **csr_ax_dense** (-21.7%) | Zig IPC 2.66 vs C++ 3.50。可能是循环调度/FMA 延迟差异 | 需 objdump + perf annotate 对照每条非零元素指令数（CSR multiply 在 ReleaseFast 被 inline 到 main） |
| **sparse_accumulator** 双峰 | FAST IPC 3.0, SLOW IPC 1.1，instructions 相同，per-process。页偏移始终 0/64。C++ 无双峰 | 需 Zen2 `ls_dispatch` 微架构事件；需 smp/c_alloc/page_alloc 对照 + 固定预分配 buffer 实验 |
| **transpose_owning** (-78.9%) | `[]usize` (8B) vs `vector<HighsInt>` (4B) 导致 2× 内存流量 | 需 perf stat 确认 hotspot，内部使用 `[]HUInt` 对照 |
| **builder branch-misses** (46×) | perf report 显示 44% 在 `@memcpy` 缺页 + 31% 在内核 `do_anonymous_page`。reusable 消除后 branch-misses 降 37× | 假设成立 |

### FIXED — 已消除

| Kernel | 修复 | 状态 |
|---|---|---|
| `computeLayout` align overflow | `al > 1 and cursor > max - (al - 1)` 检查 | ✅ |
| 全部 6 处 compact layout | 迁移到 `computePageColoredLayout` | ✅ |
| 全部 3 处 lean layout | 迁移到 `computeLayout` | ✅ |
| csc_to_csr_owning struct hash | 两端均使用 `fillCsrFromCsc` 输出 | ✅ |
| builder struct hash (3 kernel) | Z/C++ 均输出 struct hash | ✅ |
| 脚本 `\|\| true` | 移除 | ✅ |
| 脚本逐轮校验 | 跨实现 checksum + struct hash | ✅ |
| `clear_output` 测量污染 | repeats 500→100k | ✅ |

## 4. 方向：caller-owned reusable CSC build buffers

build 系列 kernel 的 owning 路径瓶颈在每次分配新内存时的页 fault 开销，非 merge 算法。使用 caller-owned reusable buffer 可获得显著加速。

### 生产 API（已实现）

```zig
pub const CscBuildBuffers = struct { ... };  // caller 分配/管理生命周期
pub fn freezeCanonicalIntoAssumeValid(buffers: *CscBuildBuffers, ...) MatrixError!CscView;
```

- 返回借用 `CscView`（`[]const` slices），不拥有 buffer
- 容量不足返回 `error.BufferTooSmall`
- 不强制预触页（caller 策略）
- 适用于求解器中反复 rebuild 矩阵的场景

### 生产路径性能（11 轮交错采样）

Reusable CSC 构建在同数学语义测试中比 C++ AoS reusable reference 快约
42%；两端 checksum 与 struct hash 均匹配。该结果的 Zig MAD 为 0.33%，
C++ MAD 为 0.46%，是本轮可信度较高的结果。

## 5. 待执行

| 优先级 | 任务 | 需 sudo |
|---|---|---|
| P0 | csr_ax_dense objdump + perf annotate 对照 | 是 |
| P1 | sparse_accumulator smp/c_alloc/fixed-buffer 对照 + Zen2 `ls_dispatch` | 是 |
| P2 | transpose owning 内部 `[]HUInt` 对照实验 | 否 |
| P3 | 低负载 (< 1.0) 11 轮完整重采样 | 否 |
