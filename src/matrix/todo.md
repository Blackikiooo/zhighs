# Matrix performance roadmap

本文件记录 `src/matrix` 及其上层集成路径中，与 DOD、CPU cache 利用率、
zero-copy 和可复用工作区相关的后续工作。任务按性能影响和调用热度排序；结构
变换中不可避免的输出复制不以“零复制”为目标，而以单次遍历、单次提交、连续
布局和可复用内存为目标。

状态说明：

- `[x]` 已实现并通过当前测试或 benchmark。
- `[ ]` 尚未实现或仍需补充验证。

## 已具备的性能基线

- [x] CSC 使用 `row_indices` / `values` SoA 连续布局。
- [x] `CscView`、`ProblemView` 等借用视图不会复制模型矩阵。
- [x] `MatrixBuilder` 使用 `MultiArrayList`，并提供预留容量的可信构建路径。
- [x] `SparseAccumulator` 使用 dense workspace + active IDs，并支持容量复用。
- [x] `DynamicRowMatrix` 使用 append-only CSR-like SoA 和 O(1) checkpoint/rollback。
- [x] `MatrixStore` 按 revision 延迟构建并复用 CSR cache。
- [x] CSC/CSR/transpose/builder 已提供部分 caller-owned `...Into`/buffer API。
- [x] 已建立 ReleaseFast `bench-matrix` 基准。

## P0：求解热路径

### Cache-friendly BTRAN

- [ ] 重写 row-major `DenseLU.solveTranspose`，消除
  `lu[j * n + i]` 形式的按列跨步读取。
- [ ] 使用 row-contiguous scatter/update 形式完成转置三角求解，不额外长期保存
  一份转置 LU，避免用双份存储换取 cache 命中率。
- [ ] 保持 RHS/workspace caller-visible 生命周期和现有 allocation-free solve 语义。
- [ ] 增加不同 basis 维度下的 ReleaseFast FTRAN/BTRAN 独立基准，并记录
  cache miss、TLB miss 和耗时变化。
- [ ] 对照现有实现做随机矩阵数值差分，验证 pivot permutation 和浮点容差。

说明：当前 dual simplex 的 tableau row、reduced cost 重算和 dual steepest-edge
初始化都会调用 BTRAN，因此这是当前最优先的 matrix 热路径问题。Dense LU 仍只应
作为小 basis fallback 和正确性 oracle；商业规模最终依赖 sparse LU、稀疏
FTRAN/BTRAN 和 Forrest--Tomlin 类更新。

## P0：矩阵变更与 revision 提交

### 批量系数修改

- [ ] 禁止每个 `chg_coeff` 单独复制完整 `col_starts`、`row_indices` 和 `values`。
- [ ] 对只修改既有非零数值的 batch，最多复制/写入一次 values storage。
- [ ] 对插入或删除非零的 batch，收集按 `(column, row)` 排序的 delta，并在一次
  canonical merge 中生成新 CSC。
- [ ] 一个 batch 只递增一次 matrix revision、只使派生 cache 失效一次。
- [ ] 为重复位置、删除为零、插入新非零和混合修改建立差分测试与吞吐基准。

### 统一 packed CSC 输出

- [ ] 为 `edit`、`slice`、`permutation` 和 model rebuild 建立共同的 packed CSC
  输出分配器，保持 64-byte alignment、page coloring 和一致的 ownership。
- [ ] 提供 owning 返回和 caller-owned `...Into` 两套边界；热路径优先复用容量。
- [ ] 保持 strong failure safety：构建成功后再原子替换 authoritative matrix。
- [ ] 避免结构变换后退化为三个彼此独立、未对齐且没有 compact metadata 的分配。

## P1：将复用缓冲区接入生产路径

- [ ] 让 `MatrixStore` 或上层 solve/compile session 持有可复用 CSR row cursor，
  避免每次 revision rebuild 都分配临时 `next`。
- [ ] 将 `CsrBuffers`、`TransposeBuffers`、`CscBuildBuffers` 接入 presolve、模型编译
  和重复矩阵转换路径，而不只用于测试和 benchmark。
- [ ] 为 buffer capacity 建立增长策略，保留高水位容量并提供显式 shrink/release。
- [ ] 所有返回借用内存的类型继续使用 `View` 后缀，并在 doc comment 中说明失效
  条件和 owning buffer 生命周期。

2026-07-13 当前机器上的 ReleaseFast 基准结果：

| kernel | 每次耗时 |
| --- | ---: |
| CSC -> CSR，重复构建输出 | 约 741 us |
| CSC -> CSR，复用 `Into` buffers | 约 246 us |
| transpose，重复构建输出 | 约 958 us |
| transpose，复用 `Into` buffers | 约 301 us |

完整复用在这两个转换中约有 3 倍收益，因此应优先把现有 API 接入真实调用链，
再考虑更复杂的 SIMD 或 prefetch。

## P1：动态行与批量结构编辑

- [ ] 为 `DynamicRowMatrix` 增加直接消费其 `row_starts`、`col_indices`、`values`
  的 `appendRowsFromCsrView`/等价接口。
- [ ] 删除 `appendToCsc` 中按动态行数量分配和填充 `SparseVectorView` 数组的中间层。
- [ ] 直接预计算各列新增 nnz，并一次写入最终 packed CSC。
- [ ] 增加 MIP cut/presolve 场景的批量追加、rollback 和同步 benchmark。

## P2：构建期 DOD 与单行操作

### SparseVectorBuilder

- [ ] 将 `Entry { id, value }` AoS 评估并迁移为 `MultiArrayList`/SoA，降低排序时
  无关字段搬运和 freeze 时的 field split 成本。
- [ ] 增加 `freezeInto`，允许 cut generation/presolve 复用 indices/values 输出缓冲区。
- [ ] 保留重复项浮点合并的确定性，并对 stable sort 与带 sequence 的 unstable
  sort 做基准后再选择。

### Row scaling

- [ ] 将 `scaleRow` 明确标记为冷路径，因为单次调用会扫描整个 CSC 两次。
- [ ] 批量缩放统一使用 `applyRowFactors`；已有 CSR cache 时评估 row-contiguous
  的单行缩放接口。
- [ ] 增加逐行调用与批量 factors 的复杂度和性能回归测试，防止出现
  `O(num_rows * nnz)` 的调用方式。

### Model rebuild 中间存储

- [ ] 删除新增变量/约束路径中只累计但未参与构建的 `col_nz_counts`。
- [ ] 删除变量/约束时先计算 surviving nnz，再直接写最终 packed storage，避免
  `ArrayList` 构建后再次 `dupe` 全部 row/value 数据。

## P2：CSC offset 宽度和视图设计

- [ ] 测量 simplex 主路径中 `usize` column offsets 与 `HUInt` compact offsets 的
  cache、带宽和整体迭代差异。
- [ ] 避免默认长期维护两套同步 offsets，而主 `CscView` 只能暴露 `usize` 的状态。
- [ ] 根据实测选择一种统一设计：常规模型使用 `HUInt`、超大模型切换宽 offset，
  或使用泛型/带标签的 view 暴露真实 offset 类型。
- [ ] 在没有端到端收益证据前，不为了减少 offset 字节数增加每个 nnz 的分支。

## 已验证且暂不修改

- [x] 保留当前 `clearF64` volatile SIMD 默认策略。当前 ReleaseFast 基准中，
  `clearF64` 清零 50,000 个 `f64` 约 7.66 us，而
  `@memset(sliceAsBytes(...))` 约 93.4 us；本机上前者约快 12 倍。
- [x] 保留 CSC 作为 authoritative storage、CSR 作为按 revision 缓存的设计。
  CSR 的 row-contiguous 访问收益通常高于仅保存 position map 后随机读取 CSC values。
- [x] 不要求结构变换实现不可能的严格 zero-copy；目标是避免重复中间复制、复用
  scratch，并让最终 owning 输出保持最适合后续热循环的布局。

## 验收要求

- 所有优化必须通过 `zig build test` 和 w32/w64 配置测试。
- 性能结论必须使用 ReleaseFast、固定数据集和足够重复次数验证。
- 热路径不得在 simplex iteration 内新增 allocator 调用。
- 新增 owning/view API 必须明确 ownership、借用失效条件和 revision 关系。
- 结构变更必须保持 canonical CSC：列偏移合法、行索引严格递增、无显式零、无重复项。
- 每个重大优化同时记录正确性差分、耗时、吞吐量、峰值内存及适用矩阵规模。
