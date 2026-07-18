# Matrix 模块状态与路线

状态说明：

- `[x]` 已实现并通过 acceptance gate。
- `[ ]` 尚未实现。
- `[*]` 阻塞中 — 等待依赖项就绪。

## 文件总览

24 个源文件，全部已实现。acceptance gate 通过（4/4）。

| 文件 | 行数 | 职责 |
| --- | --- | --- |
| `root.zig` | 171 | 模块入口，公共 API re-export，API stability 标注 |
| `target_policy.zig` | 76 | 编译期架构特化：cache line / vector lanes / prefetch / unroll |
| `csc.zig` | 585 | Canonical CSC：`CscMatrix`（owning）+ `CscView`（borrowed）+ compact starts |
| `csr_view.zig` | 422 | CSR cache：`CsrView` + `CsrCache` + reusable `CsrBuffers` |
| `sparse_vector.zig` | 303 | Canonical `SparseVectorView` / `SparseVector`（SoA packed） |
| `sparse_vector_builder.zig` | 214 | SoA `SparseVectorBuilder`：MultiArrayList append → freeze |
| `builder.zig` | 810 | `MatrixBuilder`：SoA triplet → sort/merge → canonical CSC freeze |
| `store.zig` | 209 | `MatrixStore`：owning CSC + revision + lazy CSR cache 一致性 |
| `ops.zig` | 404 | SpMV / norms / statistics / `@mulAdd` FMA / 显式 `@Vector` |
| `memory.zig` | 244 | SIMD clear（volatile + fast）、computeLayout、pageColoredLayout |
| `dense_lu.zig` | 388 | Dense LU：partial pivot + row-major cache-friendly BTRAN |
| `transpose.zig` | 389 | CSC→CSCᵀ：owning / lean(compact) / reusable `Into` |
| `slice.zig` | 281 | 行列范围提取：owning + reusable `Into` |
| `scaling.zig` | 286 | Row/col scaling：transactional preflight + `@Vector` 显式向量化 |
| `permutation.zig` | 211 | Row/col 置换：owning + reusable `Into` |
| `edit.zig` | 631 | 结构编辑：append/delete rows/cols，owning + `Into` + Csr-based append |
| `dynamic_rows.zig` | 203 | Append-only dynamic rows：CSR-like SoA + O(1) checkpoint/rollback |
| `transform_buffers.zig` | 81 | `CscTransformBuffers`：edit/slice/permute 共享 scratch，64B 对齐 |
| `sparse_sum.zig` | 268 | `SparseAccumulator`：dense workspace + active set + generation mark |
| `sparse_basis.zig` | 223 | Simplex basis CSC 组装：SoA streams + 64B-aligned capacity reuse |
| `sparse_kernel.zig` | 1003 | Mutable sparse elimination kernel：SoA 池 + 侵入式行列链表 |
| `sparse_symbolic.zig` | 489 | 符号 pivot 分析：singleton 队列 + adaptive Markowitz 搜索 |
| `sparse_lu.zig` | 1194 | Packed sparse LU：Markowitz 排序 + FT 更新 + pivot trace replay/repair |
| `sparse_ft.zig` | 402 | Sparse Forrest–Tomlin：双链 U factor + row correction workspace |

## P0：已知性能差距（与 HiGHS Release 对照）

以下数据来自 2026-07-14 fair-allocator synthetic benchmark（同一 CPU、交替进程）。精确数据会随代码变更过期，应以 `bench/matrix/results.md` 和 `bench/matrix/real_dataset_results.md` 为准。

- [ ] **`alpha_ax_plus_y`** — 慢约 23%，当前最大单项差距。原因：CSR SpMV 后独立 `y += alpha*x` 循环，两次遍历而非融合。修复方向：单次 CSR traversal 中写入 `y[i] = alpha * sum + beta * y[i]`，但需测量 store-forwarding 影响。
- [ ] **Owning transpose** — 慢约 6%。拆为 allocation + histogram + scatter + compact restore 逐段对照 HiGHS，确认 bottleneck 在哪个阶段。
- [ ] **Owning CSR** — 慢约 5%。可能来自 `CsrBuffers` cursor 与 persistent output 共享 allocation 的 layout 代价。
- [ ] **双 offsets 存储** — `CscMatrix` 同时持有 `col_starts: []usize` 和 `compact_col_starts: []HUInt`。热 kernel 已消费 compact starts，但双份 offset 在 webbase-1M 多占约 7.63 MiB。需要统一为单份 offset 表示（tagged union 或默认 HUInt + 超大模型 fallback），且不增加每个 nnz 的分支。

## P1：矩阵编辑管线收尾

当前 `ModelEditPlan` 已实现 coefficient 合并 + bounds/obj/RHS/sense/type + added rows/columns + deleted ID。剩余：

- [ ] **Coefficient add 语义** — 公开 API 只有 set，需 add（`chgCoeff(row, col, delta)` 叠加而非替换）。
- [ ] **Name 更新** — pending change 尚无 name-edit variant。
- [ ] **Add-then-delete 消除** — 同一 segment 内加后即删的无效操作应被折叠。
- [ ] **Stable ID 集中解析** — 避免执行阶段为每个 edit 重复 lookup handle table。
- [ ] **细分 plan type 缓存** — `bounds_only` / `coefficients_only` / `append_rows` / `append_columns` / `general_structural_merge`，当前只有 `scalar_only` / `coefficients_only` / `mixed_nonstructural` / `structural` 粗粒度分类。

## P1：将复用缓冲区接入调用链

已有 `CsrBuffers`、`TransposeBuffers`、`CscTransformBuffers`、`CscBuildBuffers` 均提供 caller-owned 复用路径（CSC→CSR 约 3 倍加速，transpose 约 3 倍加速）。但当前仅 `MatrixStore` 持有 `CsrBuffers` 跨 revision 复用。

- [ ] compile/presolve session 接入 `TransposeBuffers` + `CscBuildBuffers`（阻塞：presolve 仍为空模块，LP compile 走 zero-copy view 不需要转换）。
- [ ] buffer capacity 增长策略：保留高水位，提供显式 shrink/release。
- [ ] `SparseVectorBuilder.freezeInto` — cut generation / presolve 复用输出缓冲。
- [ ] `DynamicRowMatrix` 批量追加 + rollback benchmark（MIP cut 场景）。

## P1：Compile/Presolve pass fusion

- [ ] Pass 1：统计 surviving dims/nnz，生成 row/col remap + scaling，不物化中间矩阵。
- [ ] Pass 2：一次写出 final packed CSC，合并 remap + scaling + 数值过滤 + compact offset。
- [ ] compile plan 按结构 revision 缓存；仅 bounds/objective revision 变化时复用结构分析。

## P2：按需的 fused kernel

仅在 profiler 确认单次 traversal 能同时消除多次 dense temporary 且不因寄存器压力退化时才实现：

- [ ] CSR `y = alpha*A*x + beta*y`（`alpha_ax_plus_y` 修复）。
- [ ] CSC `reduced_cost = c - A^T*pi` + row activity（simplex pricing 双输出遍历合并）。
- [ ] scaling + remap + CSC copy/freeze 最终写出 traversal 合并。

每个 fused kernel 分别测量 bytes read/write、temporary bytes、cache miss、总耗时。遍历次数减少但寄存器溢出增加时不得合入。热循环内不引入动态 dispatch 或 allocator。

## P2：冷路径清理

- [ ] `scaleRow` 已标注为冷路径（单次调用扫描 CSC 两次）。在 simplex/presolve 接入后检查是否有意外的逐行调用路径，评估是否需要批量接口或 assertion guard。

## 阻塞项（依赖其他模块）

- [*] Presolve 接入 `TransposeBuffers` / `CscBuildBuffers` — 等 presolve 模块有实际消费端。
- [*] `DynamicRowMatrix` batch merge benchmark — 等 MIP cut generation 链路可用。
- [*] `freezeInto` for `SparseVectorBuilder` — 等 presolve/cut 有具体调用方，否则设计的 contract 对不齐真实需求。

## 已验证且拒绝的建议

以下建议已通过 benchmark 验证后明确拒绝。除非有新 perf 证据，不再重新评估：

- `@memset` 替代 volatile SIMD clear — 慢约 12 倍。
- CSC 改 CSR authoritative — CSR 可通过 cache 获得，无需改 authoritative。
- `page_allocator` 默认 — 小矩阵慢 14–23 倍。
- 4 KiB stack fallback 替代 alloc/free — 无稳定收益。
- dense SpMV 手工展开 — `@mulAdd` + noinline leaf 已优于 HiGHS。
- CSC scatter / transpose scatter / row scaling gather 强制 SIMD — 间接寻址和写依赖无收益。
- `Arena` retain/reset 硬编码 — 收益显著（sorted build 降至 smp 的 0.17–0.28 倍），但会将返回对象生命周期绑定到 session reset。保留为调用方可选策略，不在 matrix 内部默认。

## 开发纪律

- 每项只改一个变量。固定 CPU、交替顺序，报告 median/MAD/min/max。
- 热路径不新增 allocator 调用。
- 性能结论使用 ReleaseFast + 固定数据集 + 足够重复。
- 新增 API 必须明确 ownership、借用失效条件、revision 关系。
- 结构变更保持 canonical CSC：偏移合法、行索引严格递增、无显式零、无重复。
- 重大优化同时记录正确性差分、耗时、吞吐、峰值内存。

## 验收

`tools/matrix_acceptance.sh full`（fail-closed，4 gate）— 当前 **PASS**。所有矩阵 kernel 的 correctness、structural fuzz/OOM、三组真实大数据集 (SuiteSparse)、w32/w64 × Debug/ReleaseSafe/ReleaseFast 全部通过。
