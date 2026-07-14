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

- [x] HiGHS differential：使用本地 pinned checkout
  `/home/godv/codefiles/cppfiles/HiGHS`，commit
  `de09bbad9fb7c5d39a1a464a7641bbb5531c6e9d`；21 个 matrix kernels 的 checksum 与
  structural hash 全部匹配。
- [x] Structural property/fuzz/OOM：100 组确定性随机结构变换和 owning 路径 exhaustive
  failing-allocator 已通过 w32/w64。
- [x] Large real datasets：SuiteSparse `thermal1`、`cage12`、`webbase-1M` 已固定来源和
  SHA-256；in-tree ReleaseFast runner 验证 canonical CSC、CSC/CSR/transpose 语义、
  scaling、slice/permutation、运行时间和 peak RSS，三组均 PASS。完整数据与 HiGHS
  同机对比见 `bench/matrix/real_dataset_results.md`。
- [x] Configuration regression：Debug/ReleaseSafe/ReleaseFast x w32/w64 已全部通过。

当前 full verdict：`PASS (4/4 gates passed)`。2026-07-14 的 fail-closed full run 已
实际验证 HiGHS differential、structural fuzz/OOM、large real datasets 和
configuration regression 全部 PASS。按约定硬门槛，matrix 当前已达到其已实现功能
范围内的 production-candidate 终止点；这不等价于完整 solver 已满足商业化要求。

- 所有优化必须通过 `zig build test` 和 w32/w64 配置测试。
- 性能结论必须使用 ReleaseFast、固定数据集和足够重复次数验证。
- 热路径不得在 simplex iteration 内新增 allocator 调用。
- 新增 owning/view API 必须明确 ownership、借用失效条件和 revision 关系。
- 结构变更必须保持 canonical CSC：列偏移合法、行索引严格递增、无显式零、无重复项。
- 每个重大优化同时记录正确性差分、耗时、吞吐量、峰值内存及适用矩阵规模。

## 2026-07-14 矩阵内核专项：布局、公平性、perf 与反汇编

本专项暂时高于 model/solver 集成路线；在下列证据闭环完成前，不根据单次 wall-clock
结果修改热路径。每项优化必须记录优化前后 `perf stat`、热点符号/指令、正确性回归和
真实数据集结果。

### A. Zig 内存布局审计

- [x] 建立可重复运行的 layout audit，记录 `@sizeOf`、`@alignOf`、字段 offset、切片/可选
  切片开销，以及 w32/w64 差异；覆盖 `CscMatrix`、`CsrView`、`CsrBuffers`、builder
  triplet、transpose buffers 和稀疏向量类型。
- [ ] 按实际分配而非仅按 struct header 统计 CSC/CSR authoritative bytes、scratch bytes、
  page-color padding 和 allocator high-water；禁止把完整验收流水线 RSS 当成单矩阵内存。
- [ ] 复核所有热数据结构的 AoS/SoA、字段宽度和生命周期；结构体仅用于控制块/借用视图，
  大规模元素优先连续 field streams，并以 cache-line/bytes-per-nnz 证据决定布局。
- [ ] 消除默认 compact CSC 同时长期保存 `usize` 与 `HUInt` 两套 offsets、但 dense CSC
  SpMV 只读取 `usize` 的无收益状态；先完成 perf/汇编 A/B，再选择统一 offset 表示。
- [ ] 将构建期内存与运行期内存分开：审计 Matrix Market 整文件缓冲、24-byte/nnz
  builder triplet（含 `sequence: usize`）及 freeze 同时存活峰值。

### B. zhighs / HiGHS 公平性审计

- [ ] 固定同一 CPU core、SMT 状态、数据、矩阵顺序、warm-up、重复次数及交错执行顺序；
  报告 median、MAD 和进程级重复，不以单进程内部均值代替统计稳定性。
- [x] 核对 Zig `ReleaseFast -Dcpu=native` 与 C++ `-O3 -march=native -DNDEBUG -flto` 的
  ISA、LTO、链接方式和断言/边界检查；记录编译器版本与二进制哈希。
- [ ] 确保双方 kernel contract 一致：是否包含 output clear、是否复用输出容量、offset/index
  宽度、输入稀疏模式、checksum barrier 和 owning/borrowed 转换必须逐项相同。
- [ ] 将真实数据集 runner 拆成 kernel-only 与 lifecycle/memory 两类；HiGHS 与 zhighs 的
  CSC-only、CSR-only、CSC+CSR、转换 scratch 和构建峰值必须同口径。
- [x] 初步代码审计确认：双方 CSC/CSR dense SpMV 算法同阶且循环结构接近；当前验收 RSS
  不公平，zhighs 同时持有 transpose、values 副本、scaling、slice/permutation 等对象，
  HiGHS runner 主要持有 CSC+CSR，现有 RSS 不得用于宣称单矩阵内存差距。

### C. perf 与汇编定位

- [x] 对 `csc_ax_dense`、`csr_ax_dense`、`clear_output`、`csc_to_csr_into` 先执行固定核
  `perf stat -r`，至少记录 cycles、instructions、IPC、branches、branch-misses、
  cache-references/misses、page-faults 和 context-switches；WSL/PMU 不支持的事件明确标记。
- [x] 使用 `perf record`/`perf report` 确认采样落在目标 kernel 而非 harness、allocator、
  checksum 或动态链接层；必要时添加稳定的 noinline profiling boundary。
- [x] 对存在稳定差距的 kernel，用 `objdump`/`llvm-objdump` 按符号提取双方汇编，比较
  清零实现、offset 扩展、循环分支、地址生成、load/store 数量和寄存器 spill。
- [ ] 为 `usize offsets`、`HUInt offsets`、volatile clear、fast clear 建立独立 A/B；每次
  只改变一个变量，汇编和 perf 同时证明原因后才修改生产默认路径。

### D. 已识别候选修复（等待 perf 证据）

- [x] CSC dense SpMV 当前忽略已存在的 `compact_col_starts`，导致双份 offset 内存没有转化
  为热循环收益；评估 compact 专用 kernel 或统一 32-bit/64-bit tagged view。
- [ ] `CsrBuffers` 将 row cursor 与常驻 CSR 输出放在同一 allocation，重复转换有利但
  常驻多 `4 * num_rows`（w32）；拆分 persistent view 与 reusable workspace 的生命周期。
- [x] 重新验证 `clearF64` 的 volatile SIMD 默认策略。既有“比 memset 快 12 倍”结论需在
  同编译器、同输出语义和汇编下复测，排除旧 benchmark 的 memset 类型/优化消除问题。
- [ ] builder 的 `sequence: usize` 令 w32 triplet field streams 达 24 bytes/nnz；比较稳定
  sort、32-bit sequence、分块/流式构建，保证重复项浮点合并顺序后再选择。

专项完成标准：公平基准中主要矩阵 kernel 的差距有硬件事件和汇编解释；修复后 quick/full
gate 通过，真实数据与 synthetic profile 方向一致，且 `todo.md` 保留失败实验而不只记录
成功结果。

### 2026-07-14 专项实验记录

- Layout audit 已加入 `zig build audit-matrix-layout`。w32 下 `CscMatrix`/`CscView`/
  `CsrView` 控制块分别为 104/64/64 bytes；slice 与 optional slice 均为 16 bytes。
  普通 Zig struct 为 auto layout：builder triplet 的实际字段 offset 为 value=0、
  sequence=8、row=16、col=20，总计 24 bytes；w64 为 32 bytes。控制块大小不是百万级
  矩阵 RSS 主因，外部 SoA arrays、双 offsets 和 scratch 生命周期才是主因。
- 已安装并验证 perf 6.8.12；WSL2 用户态 PMU 可用。初始 synthetic 结果：CSC Zig/HiGHS
  约 363M/318M cycles，CSR 约 355M/284M cycles；branch miss 均极低，cache miss 没有
  解释差距。CSC->CSR Zig 约 136M cycles，HiGHS 约 175M cycles，证明不是 Zig/ID/SoA
  的普遍开销。
- 反汇编确认 HiGHS/GCC 使用 `vfmadd*sd`，原 Zig 使用分离 `vmulsd+vaddsd`；现已在
  CSC/CSR dense kernels 显式使用 `@mulAdd`，w32/w64 matrix tests 通过，生成 FMA 已核实。
- 更大的差距来自 benchmark/代码生成边界：Zig kernel 被内联进包含 21 个分支的巨大
  profiling `main`，导致 values/indices/y bases 在每个 nnz 从栈重载；HiGHS shared-library
  `product` 是独立 leaf，bases 常驻寄存器。给 CSC/CSR 建立稳定 noinline leaf 后，栈重载
  消失；synthetic CSR 降到约 222M cycles，优于 HiGHS 约 283M cycles。
- CSC dense kernel 已按一次性 tag dispatch 使用 `compact_col_starts`，每个 nnz 不增加
  分支；synthetic CSC 约 282M cycles，HiGHS 约 318M cycles。5 个独立真实数据进程的
  zhighs median 为 thermal/cage/web CSC 0.561/1.774/4.712 ms，CSR
  0.526/1.648/4.314 ms；5 个 HiGHS 进程 median 分别为 CSC
  0.602/1.842/5.162 ms，CSR 0.540/1.845/5.094 ms。当前为分组重复而非逐次交错，最终
  公平报告仍需完成真实数据 interleaved runner 后固化。
- `clear_output` perf 中 Zig 约 3.28B cycles/100k，HiGHS 约 3.46B cycles/100k；HiGHS
  `std::fill` 主要降为 `rep stos`，因此 retired instructions 不可直接比较。清零不是当前
  CSC/CSR 差距来源，保留 volatile clear 的结论暂时成立。
- 修复后的 `zig build test`（w32/w64）、quick gate 和 fail-closed full gate 均通过；
  full report 为 `/tmp/zhighs-matrix-acceptance/20260714T065632Z/summary.tsv`，四项全部 PASS。
  首次 full 尝试因未传 `MATRIX_DATASET_RUNNER` 而非代码失败，补齐 runner 路径后通过；
  该环境配置失败保留记录，不计作内核回归。

## 2026-07-14 Matrix 布局、显式 SIMD 与 allocator 专项

本批次坚持单变量实验：auto-layout/SoA 审计、显式 `@Vector` 和 allocator 策略分别测量，
不得用一个组合改动宣称其中任一项有效。返回 owning 对象的 API 必须继续由调用方 allocator
分配和释放；特殊 allocator 只用于生命周期边界清晰的 session 或短命 scratch。

### E. 布局与生命周期

- [x] 逐项计算 CSC、CSR、transpose、builder、sparse vector 的 authoritative bytes、scratch
  bytes、padding 和 header，覆盖空/小/中/百万 nnz；结果写入 `bench/matrix/layout_results.md`。
  authoritative owning data 若 page coloring 无可重复 cycles 收益且增加超过一页或 5%，则拒绝。
- [x] 将默认 `transposeAssumeValid` 从“最终结果永久携带 cursor scratch”的布局切换为已经
  存在的 lean compact 路径，前提是时间、结构 hash 和峰值字节实验通过。
- [x] 评估 `CsrBuffers`/`TransposeBuffers` 是否应拆为 output storage 与 workspace；复用型
  buffers 可以组合持有，但发布的 owning matrix/cache 不得保留无消费者的 scratch。
- [x] 审计 Zig auto-layout 控制块中 slice 数量和按值传递边界；只在测得调用/寄存器开销时
  缩减 header，不为了几十字节破坏清晰 ownership。

布局结论：Zig 普通 struct 是 auto-layout，字段顺序不能作为 ABI/padding 优化手段；当前
`CscMatrix` 104 B、view 64 B，真正的成本在 SoA 数据流而非控制头。w32 triplet builder 的
`MultiArrayList` 为 24 B/entry（w64 为 32 B），SoA 方向正确，但 sequence 流仍为 8 B/entry。
默认 compact builder 已由 page-colored 改为自然 packed + 64 B stream alignment：16 列/32 nnz
由 12,736 B 降为 704 B，synthetic CSC sampled cycles 288.5M -> 284.5M。默认 owning transpose
改为复用 compact starts 作为临时 cursor，scatter 后恢复最终 offsets；50k 列案例减少约
200 KiB（约旧结果 8%）常驻空间且不再发生第二次分配，五进程 synthetic median
约 982 us -> 906 us（快 7.7%），checksum/structural hash 相同。复用型
`CsrBuffers`/`TransposeBuffers` 继续允许 output+workspace 组合持有；发布的
`CsrCache` 和 owning transpose 均不再把 cursor 当 authoritative bytes。双 wide/compact starts
在 webbase-1M 仍增加约 7.63 MiB，因 compact starts 已被热内核消费，留作后续公共 API 迁移项。

### F. 显式 `@Vector` 候选

- [x] 对连续列 value scaling 建立 scalar/auto-vector/manual `@Vector` A/B，并检查 LLVM 汇编；
  checked preflight 与 mutation 分开测量，保持 failure atomicity。
- [x] 对 `absoluteRange`/`maxAbs`/`assessValues` 评估显式 abs/min/max/compare reduction；验证
  NaN、Inf、signed zero 语义和短 slice 尾部。
- [x] 明确拒绝无收益候选：CSC scatter SpMV、transpose/CSR scatter、row scaling 的 gather/
  scatter 或循环依赖若 perf 不支持，不引入手工 lane 构造。
- [x] 仅当 perf cycles、instructions 与真实调用基准稳定改善且代码尺寸可接受时替换生产实现；
  生产实现已替换，本批次 w32/w64、完整工程测试与 full acceptance 已收口。

SIMD 结论：ReleaseFast scalar 列循环反汇编为 `vmulsd`，显式路径为 `vmulpd`。8/16/32/128
nnz/列的 manual/scalar 时间比为 0.659/0.636/0.624/0.773，3 nnz/列为 0.993（中性）；因此
短列自然走 scalar tail，较长连续列使用 native vector。`absoluteRange` 在 4K/131K/2M values
的比值为 0.249/0.252/0.358，并由 `maxAbs` 复用。canonical matrix 禁止 NaN/Inf；新增测试覆盖
非对齐 borrowed slice、signed zero、vector bulk 和 tail。`assessValues` 的计数循环、CSC scatter
SpMV、transpose/CSR scatter、row scaling 均因早退、间接寻址或写依赖不强制向量化。

### G. allocator 场景实验

- [x] 比较 `smp_allocator`、`page_allocator`、Arena retain/reset 对短命批量 matrix build 的
  时间、逻辑分配次数和 Arena retained capacity；Arena 只允许绑定 compile/presolve session 生命周期。
- [x] 对 CSR/transpose 短命 cursor 比较普通 alloc/free、caller-owned reusable scratch 与小型
  stack/fixed-buffer fallback；必须覆盖小矩阵和大矩阵回退。
- [x] 保持 SpMV、norm、scaling 等迭代热循环 allocation-free；不得用 allocator 替代工作区复用。
- [x] 特殊 allocator 只有在端到端场景净收益且 ownership/deinit 明确时进入生产路径；否则记录
  “调用方可选策略”，不在 matrix 内核硬编码。

allocator 结论：4 KiB stack fallback 在 64/256/1024/4096 个 HUInt scratch 上分别为普通
alloc/free 的约 1.69/1.20/1.04/1.03 倍，拒绝替换；page allocator 在 64 维 build 慢 13.7--22.6 倍，
512/4096 维无稳定收益，也拒绝。Arena retain/reset 对重复 512/4096 维 sorted build 降至
`smp_allocator` 的 0.17--0.28/约 0.28，但会将返回对象生命周期绑定到 session reset，因此不在 matrix
内部硬编码。`MatrixBuilder.freeze` 已明确 construction allocator 与 owning output allocator 可
不同；compile/presolve 可由调用方用 Arena 管 triplet，再用长生命周期 allocator 输出矩阵。
现有 reusable build/CSR/transpose-into buffers 是 ownership 清晰的生产方案。
实验中的预留 sorted build 每轮有两次逻辑分配（triplet SoA 与 owning CSC）；Arena warm-up
之后 reset-retain 不再请求 backing allocation，64/512/4096 维 retained capacity 分别约
10.6/81.5/648.5 KiB。这里不把 Arena retention 当成 independently owned matrix 的峰值内存优势。

本批次验证：`zig build test-matrix -Doptimize=ReleaseFast` 的 w32/w64 均通过，完整
`zig build test -Doptimize=ReleaseFast` 通过，quick gate 通过；复用 compact-starts cursor 的
最终版本 fail-closed full report 为
`/tmp/zhighs-matrix-acceptance/20260714T085054Z/summary.tsv`，HiGHS differential、structural/
OOM、三份真实大数据、Debug/ReleaseSafe/ReleaseFast × w32/w64 四项全部 PASS。此前
`20260714T083447Z` 首轮仅因默认 Zig global cache 只读导致 differential 构建失败，使用
`/tmp/zhighs-zig-cache` 重跑后通过，该环境失败不计作内核回归。

### H. 修复后再次对比 HiGHS（2026-07-14）

- [x] synthetic 使用 CPU 2、11 进程、Zig/C++ 交替先后顺序；真实 SuiteSparse 使用 7 进程
  交替顺序。Zig `ReleaseFast -Dcpu=native`，C++ `-O3 -march=native -DNDEBUG -flto`，HiGHS
  commit `de09bbad9fb7c5d39a1a464a7641bbb5531c6e9d`；全部 checksum/structural hash 一致。
- [x] 真实数据中，thermal/cage/web 的 CSC SpMV 分别比 HiGHS 快 7.3%/1.0%/13.7%，CSR
  分别快 0.8%/5.9%/16.7%；CSC->CSR reusable 分别快 16.4%、慢 0.6%、快 4.2%。低于 1%
  只判定为持平，不宣称领先。完整表写入 `bench/matrix/real_dataset_results.md`。
- [x] synthetic 稳定项：CSR dense 快 20.4%，full scale 快 9.3%，CSR reusable conversion
  快 18.5%，owning CSR 快 2.2%，transpose-into 慢 1.3%（持平），reusable builder 快 10.7%，
  general builder 快 12.4%，sparse accumulate 快 66.6%。
- [ ] owning transpose 仍为 920 us vs HiGHS 337 us（Zig 慢 2.73x）；本次 cursor 复用已令
  Zig 自身快约 7.7%，但 reusable transpose 已基本持平，差距集中于 owning allocation 与
  双 wide/compact starts 输出。下一批用 perf/汇编单独定位，不再修改 scatter 内核。
- [ ] canonical owning builder 仍为 745 us vs 256 us（慢 2.91x），但 reusable builder 已快
  10.7%；下一批拆分 output allocation、wide starts 填充、输入扫描/memcpy 做单变量 perf。
- [ ] 建立等生命周期的 RSS runner。目前 acceptance Zig 进程额外保留 CSR/transpose/scaling/
  permutation 验证对象，不能与窄 C++ runner RSS 直接比较。按 authoritative w32 CSC 数组计算，
  Zig 相对 HiGHS 在 thermal/cage/web 多 9.2%/4.2%/19.4%；webbase 的双 offsets 多 7.63 MiB。

本轮结论：不能说所有 matrix 路径已经完全超过 HiGHS。SpMV、reusable CSR conversion、通用
builder 和 sparse accumulator 已达到持平或领先；主要剩余差距已收敛到 owning transpose、
canonical/sorted owning build 以及双 offsets 的内存成本。高 MAD 的 synthetic CSC/product
排名不作为结论，真实数据交错 medians 才作为 SpMV 判据。原始结果位于
`/tmp/zhighs-matrix-after-layout/results` 与 `/tmp/zhighs-matrix-after-layout/real`。
