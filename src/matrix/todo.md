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
- [x] `CscMatrix`、`CscView`、`CsrView` 和 `SparseVectorView` 使用命名构造入口；
  owning/borrowed/packed/trusted 语义不再依赖业务代码中的 struct literal。
- [x] `CscView`、`ProblemView` 等借用视图不会复制模型矩阵。
- [x] `MatrixBuilder` 使用 `MultiArrayList`，并提供预留容量的可信构建路径。
- [x] `SparseAccumulator` 使用 dense workspace + active IDs，并支持容量复用。
- [x] `DynamicRowMatrix` 使用 append-only CSR-like SoA；column/value 由一个
  `MultiArrayList` allocation 管理，并支持 O(1) checkpoint/rollback。
- [x] builder owning freeze 直接在最终 packed/page-colored storage 中生成
  `col_starts`；有序输入不再分配 starts/cursor scratch。
- [x] `MatrixStore` 按 revision 延迟构建并复用 CSR cache。
- [x] CSC/CSR/transpose/builder 已提供部分 caller-owned `...Into`/buffer API。
- [x] 已建立 ReleaseFast `bench-matrix` 基准。

## 2026-07-14 后续执行顺序

1. 补齐 packed CSC 结构变换的 caller-owned `...Into` 边界，以共同的 capacity
   contract 覆盖 edit、slice 和 permutation；这是把复用内存接入上层的前置条件。
2. 将已有 `CsrBuffers`、`TransposeBuffers`、`CscBuildBuffers` 接入 compile/presolve
   session，先兑现已测得的约 3 倍格式转换收益，并补 high-water/shrink 生命周期。
3. 分阶段实现 `ModelEditPlan`：先做 `existing_coefficients_only` 和 scalar-only plan，
   再做 append，最后才做 general structural merge；每阶段都保留小 batch 直接路径。
4. 在 compile plan 和 revision cache 稳定后实现两遍式 presolve/compile fusion，避免
   operation stream 与 compile 中间 IR 同时大改而难以定位语义回归。
5. specialized fused kernels、compact offsets、SIMD/prefetch 继续坚持 benchmark gate；
   只有端到端 bytes moved、cache 和总耗时同时给出收益证据才进入生产路径。

当前执行批次（2026-07-14）：`ModelEditPlan` 第二阶段已完成并接入 `applyPending`。
coefficient 与 bounds/objective/RHS/sense/type streams 使用 `MultiArrayList` DOD，
按 target/sequence 规范化并执行 last-write-wins；不超过 8 条的 scalar-only segment
走 allocation-free direct path。deleted variable/constraint IDs 也已连续收集、排序和
去重。added column/row metadata、offsets、indices、values 已迁入连续 plan streams，
属性初始化及单次 packed CSC rebuild 不再反复扫描原 `PendingChange`；新增对象紧随的
bounds/objective/type/RHS/sense scalar edits 会直接折叠进最终 metadata。w32/w64
`test-model` 及 model-level metadata/name/matrix 集成回归已通过；下一批次为
add-then-delete 消除、stable-ID 集中解析和更细的 append/general structural plan 分类。

本轮 ReleaseFast 基线（dimension=4096）：一次 flush 4096 rows、4096 columns 和
4096 个可折叠 objective edits 共 12,288 operations，耗时约 0.834 ms、约 14.74 M
operations/s。`bench-coefficient-edits` 已包含
`append_rows_columns_folded_scalars`，供 general structural merge 前后直接对比；无新增
对象的 plan 会跳过 folding passes，避免 scalar-only 路径多做线性扫描。

Structural plan 前置审计发现 `addConstr(num_nz > 0)` 的 CSR row payload 未进入 CSC；
已修复为与 added columns 在同一次 packed rebuild 中按列计数、scatter 和提交，并补
“existing variables + appended nonempty row”回归测试。同一 batch 的 added row/column
描述相同新坐标时采用 row payload last-set 覆盖，不生成重复项；显式回归入口位于
`test/model/root.zig`，w32/w64 `test-model` 均已通过。

## P0：求解热路径

### Cache-friendly BTRAN

- [x] 重写 row-major `DenseLU.solveTranspose`，消除
  `lu[j * n + i]` 形式的按列跨步读取。
- [x] 使用 row-contiguous scatter/update 形式完成转置三角求解，不额外长期保存
  一份转置 LU，避免用双份存储换取 cache 命中率。
- [x] 保持 RHS/workspace caller-visible 生命周期和现有 allocation-free solve 语义。
- [x] 增加不同 basis 维度下的 ReleaseFast FTRAN/BTRAN 独立基准，并记录耗时变化；
  基准入口和当前结果见 `bench/matrix/dense_lu_bench.zig` 与
  `bench/matrix/dense_lu_results.md`。
- [x] 使用 `perf stat -r 5` 补录 64/256/512 维下的 cache miss、TLB miss、
  cycles、instructions 和 IPC，并记录通用硬件事件在当前 AMD CPU 上的解释限制。
- [x] 对照原 column-gather BTRAN 实现做确定性随机矩阵数值差分，验证实际发生的
  pivot permutation、FTRAN/BTRAN 已知解和可配置浮点 pivot 容差。

说明：当前 dual simplex 的 tableau row、reduced cost 重算和 dual steepest-edge
初始化都会调用 BTRAN，因此这是当前最优先的 matrix 热路径问题。Dense LU 仍只应
作为小 basis fallback 和正确性 oracle；商业规模最终依赖 sparse LU、稀疏
FTRAN/BTRAN 和 Forrest--Tomlin 类更新。

## P0：矩阵变更与 revision 提交

### 批量系数修改

- [x] 禁止每个 `chg_coeff` 单独复制完整 `col_starts`、`row_indices` 和 `values`。
- [x] 对只修改既有非零数值的 batch，最多复制/写入一次 values storage。
- [x] 对插入或删除非零的 batch，收集按 `(column, row)` 排序的 delta，并在一次
  canonical merge 中生成新 CSC。
- [x] 一个 batch 只递增一次 matrix revision、只使派生 cache 失效一次。
- [x] 为重复位置、删除为零、插入新非零和混合修改建立差分测试与吞吐基准；
  确定性 dense-oracle 测试见 `src/model/model.zig`，ReleaseFast 吞吐入口为
  `zig build bench-coefficient-edits -Doptimize=ReleaseFast`。

### 统一 packed CSC 输出

- [x] 为 `edit`、`slice`、`permutation` 建立共同的 packed CSC 输出分配器，保持
  64-byte alignment 和一致的 ownership。
- [x] 将 model linear deletion rebuild 迁移到共同分配器，避免临时增长列表和最终
  数组的双份峰值；natural packing 与 page coloring 的适用阈值仍需 benchmark。
- [x] 提供 owning 返回和 caller-owned `...Into` 两套边界；`CscTransformBuffers`
  统一 edit、slice、permutation 的 64-byte aligned output/scratch capacity contract，
  `...Into` 返回的 `CscView` 在下一次 buffer 写入或 deinit 时失效。
- [x] 保持 strong failure safety：构建成功后再原子替换 authoritative matrix。
- [x] 避免结构变换后退化为三个彼此独立、未对齐的分配。

## P0：借鉴 Burn 的有界操作流合并

借鉴 Burn Tensor Operation Streams 的目标是减少中间存储、内存读写、分配和重复
调度，而不是在求解器中引入通用 tensor JIT。zhighs 应采用带有明确同步屏障的
延迟执行：只在 model update、compile/presolve 和少数专用 matrix kernel 边界捕获
操作，生成可缓存的计划类型，然后一次执行和提交。

推荐的数据流：

```text
用户建模操作
    -> ModelEditStream
    -> 规范化、合并、消除无效操作
    -> ModelEditPlan
    -> 一次矩阵物化与 revision 提交
```

### ModelEditStream 和 ModelEditPlan

- [x] coefficient、scalar 和 deleted-ID pending changes 已视为有限操作流，在 flush
  前按对象/字段或 dense target 分组；每类只执行规范化后的 stream。
- [x] added rows/columns（含新增对象 name metadata）已进入有限操作流；独立 name-edit
  API 尚未定义，后续随 name 更新语义补齐。
- [x] 定义第一阶段 DOD `ModelEditPlan`，分别连续保存 coefficient、bounds、
  objective、RHS、sense 和 type；执行阶段只遍历实际存在的 stream。
- [x] deleted variable/constraint IDs 已迁入连续 plan storage，并在执行前排序去重。
- [x] 将 added rows/columns 的 metadata 与 CSR/CSC payload 迁入 plan；执行阶段只从
  连续 offsets/indices/values streams 构建 packed CSC，原 pending 仅保留资源 ownership。
- [x] 合并同一 coefficient 的连续 set，采用按 sequence 的 last-write-wins。
- [ ] 为 coefficient add 语义定义并实现保持规定浮点顺序的合并；当前公开 API 只有 set。
- [x] 合并同一变量的 lower/upper/objective/type，以及同一约束的 RHS/sense 更新。
- [x] `addVar`/`addConstr` 后紧随的 bounds/objective/type/RHS/sense 修改直接折叠进
  新对象最终初始化数据。
- [ ] 定义并合并独立 name 更新；当前公开 pending change 尚无 name-edit variant。
- [ ] 消除同一未发布对象在一个 segment 内的 add-then-delete 等无可观察效果操作。
- [ ] 将 stable ID 到 dense position 的解析集中在 plan 阶段，避免执行阶段重复查找。
- [ ] 只缓存不含具体句柄的计划类型，例如 `bounds_only`、`objective_only`、
  `existing_coefficients_only`、`append_rows`、`append_columns` 和
  `general_structural_merge`。
- [x] 已定义不含具体句柄的第一阶段 `PlanKind`：`scalar_only`、
  `coefficients_only`、`mixed_nonstructural` 和 `structural`。
- [ ] 在 structural plan 完成后缓存细分的 append/general plan type；当前只生成
  轻量 discriminant，不缓存包含具体 edit 数据的 plan instance。
- [x] 为不超过 8 条的 scalar-only batch 保留 allocation-free 直接路径；ReleaseFast
  分界基准为 1/4/8/16 条约 167/280/329/113 M edits/s，超过阈值进入 DOD plan。

### Compile/presolve pass fusion

- [ ] 将 validate、problem classification、presolve analysis 和 scaling analysis
  尽可能组合在分析阶段，但不在每一步物化一份完整中间矩阵。
- [ ] 采用两遍式 fusion：Pass 1 统计 surviving dimensions/nnz、生成 row/column
  remap 和 scaling；Pass 2 一次写出最终 packed CSC。
- [ ] 在最终写出中合并 row/column remap、coefficient scaling、删除固定或冗余项、
  数值过滤和 compact offset 生成。
- [ ] compile plan 按结构 revision 和相关配置缓存；仅 bounds/objective revision
  变化时复用结构分析、映射和矩阵 storage。
- [ ] 不以不安全的严格单遍为目标；优先保证最终尺寸已知、单次物化和 strong
  failure safety。

### 专用 fused matrix kernels

- [ ] 评估并补充 `y += alpha * A * x`、`y = alpha * A * x + beta * y` 等直接
  kernel，避免一次性 dense temporary。
- [ ] 评估在一次 CSC traversal 中完成 `reduced_cost = c - A^T*pi`、row activity、
  residual 或统计量，但只合并同一阶段必然同时需要的输出。
- [ ] 评估将 scaling、remap 和 CSC copy/freeze 合并到最终写出 traversal。
- [ ] fused kernel 使用静态 Zig 函数、comptime specialization 和 caller-owned
  workspace；simplex iteration 内不引入通用动态 scheduler、哈希或 allocator。
- [ ] 分别测量减少的 bytes read/write、temporary bytes、cache miss 和总耗时；遍历
  次数减少但寄存器压力、随机写或无用输出增加时不得合入。

### ExpressionGraph execution plan

- [ ] 在 NLP 路径启用前，将递归 `ExpressionGraph.evaluate` 编译为拓扑排序的 SoA
  execution plan，确保共享 DAG 节点每次 evaluation 只计算一次。
- [ ] 编译阶段实现 constant folding、common-subexpression elimination、identity
  elimination 和安全的 strength reduction。
- [ ] 根据节点引用次数分配并复用连续 value slots，避免每次 evaluation 分配。
- [ ] 为 value、gradient、Jacobian sparsity 和后续 Hessian evaluation 设计共享分析，
  避免分别递归遍历同一 DAG。
- [ ] 只执行不会改变规定浮点语义的代数化简；禁止未经证明的浮点重结合。

### 操作流同步屏障和语义

- [ ] `optimize`、导出/写文件、读取 pending 影响的数据、获取 matrix/basis view、
  callback 暴露状态和删除已发布稳定句柄之前必须 flush。
- [ ] segment 内发生 write-after-read 或 read-after-write 且读取结果对用户可观察时，
  在读取点结束当前 segment。
- [ ] deferred validation 必须保持清晰的错误归属；不能让非法操作只在很晚的
  `optimize` 中以无法定位的方式失败。
- [ ] operation fusion 必须保持 duplicate merge 顺序、NaN/Inf/overflow 行为、
  canonical CSC、stable ID 生命周期和 revision/basis/cache 失效语义。
- [ ] 每个计划执行必须先构建完整 replacement，成功后再一次提交，失败时原模型
  和所有当前 view/cache 仍保持有效。

## P1：将复用缓冲区接入生产路径

- [x] 让 `MatrixStore` 或上层 solve/compile session 持有可复用 CSR row cursor，
  避免每次 revision rebuild 都分配临时 `next`。
- [x] `MatrixStore` 已持有并跨 revision 复用完整 `CsrBuffers`，matrix-values 或结构
  revision 失效只清除 cache stamp；容量足够时不释放 row/column/value/cursor storage。
- [ ] 在 compile/presolve 出现真实 transpose/build 重复转换调用点后，将
  `TransposeBuffers`、`CscBuildBuffers` 接入对应 session；当前 presolve 仍为空模块，
  continuous LP compile 已为 zero-copy view，不引入没有消费者的短命 buffer。
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

- [x] 为 `DynamicRowMatrix` 增加直接消费其 `row_starts` 和 MultiArrayList
  column/value field streams 的 `appendRowsFromCsrView`/等价接口。
- [x] 删除 `appendToCsc` 中按动态行数量分配和填充 `SparseVectorView` 数组的中间层。
- [x] 直接预计算各列新增 nnz，并一次写入最终 packed CSC。
- [ ] 增加 MIP cut/presolve 场景的批量追加、rollback 和同步 benchmark。

## P2：构建期 DOD 与单行操作

### SparseVectorBuilder

- [x] 将 `Entry { id, value }` AoS 迁移为 `MultiArrayList`/SoA，稳定排序直接操作
  field streams，并提供预留容量的可信 append 路径。
- [x] owning freeze 将 indices/values 放入一个 64-byte aligned packed allocation，
  同时保留旧的独立切片 ownership 兼容路径。
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

- [x] 删除新增变量/约束路径中只累计但未参与构建的 `col_nz_counts`。
- [x] 删除变量/约束时先计算 surviving nnz，再直接写最终 packed storage，避免
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

生产级终止条件已固化为 `tools/matrix_acceptance.sh full`。该工具 fail-closed：
HiGHS 自动差分、structural property/fuzz 与 failing-allocator、真实大型稀疏数据集、
w32/w64 多配置回归四个 gate 必须全部通过；缺少外部依赖或数据集按失败处理，不记为
skip-success。日常使用 `tools/matrix_acceptance.sh quick`，其成功不代表生产准入。

2026-07-14 当前 gate 状态：

- [ ] HiGHS differential：provider 和 checksum/struct-hash 校验已接入；当前缺少
  `HIGHS_SOURCE` pinned checkout，因此 full gate 正确失败。
- [x] Structural property/fuzz/OOM：100 组确定性随机结构变换和 owning 路径 exhaustive
  failing-allocator 已通过 w32/w64。
- [ ] Large real datasets：已定义严格 runner/report contract；当前缺少 pinned Matrix
  Market corpus 与 runner，因此 full gate 正确失败。
- [x] Configuration regression：Debug/ReleaseSafe/ReleaseFast x w32/w64 已全部通过。

当前 full verdict：`FAIL (2/4 gates passed)`；在四项全部为 `[x]` 前，不得将 matrix
标记为 production-ready。

- 所有优化必须通过 `zig build test` 和 w32/w64 配置测试。
- 性能结论必须使用 ReleaseFast、固定数据集和足够重复次数验证。
- 热路径不得在 simplex iteration 内新增 allocator 调用。
- 新增 owning/view API 必须明确 ownership、借用失效条件和 revision 关系。
- 结构变更必须保持 canonical CSC：列偏移合法、行索引严格递增、无显式零、无重复项。
- 每个重大优化同时记录正确性差分、耗时、吞吐量、峰值内存及适用矩阵规模。
