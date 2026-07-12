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

## 3. 实验结论汇总

### CONFIRMED — 根因已定位

| Kernel | 发现 | 证据 |
|---|---|---|
| **builder_freeze_canonical (owning)** | 瓶颈在分配器页 fault，非 merge 算法 | reusable 消除 161× page faults → owning→reusable 3.4× faster。reusable 单次对比 Z/C++: Zig 458k ns vs C++ 545k ns（约 19% 优势）。正式生产 API 稳定性与 11 轮验证进行中 |
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

### 原型性能（单次采样）

Reusable CSC 构建原型在同数学语义测试中比 C++ AoS reusable reference 快约 19%；正式生产 API 与稳定性尚待 11 轮交错验证。

## 5. 待执行

| 优先级 | 任务 | 需 sudo |
|---|---|---|
| P0 | csr_ax_dense objdump + perf annotate 对照 | 是 |
| P1 | sparse_accumulator smp/c_alloc/fixed-buffer 对照 + Zen2 `ls_dispatch` | 是 |
| P2 | transpose owning 内部 `[]HUInt` 对照实验 | 否 |
| P3 | 低负载 (< 1.0) 11 轮完整重采样 | 否 |
