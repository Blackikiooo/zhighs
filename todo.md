# zhighs 实现路线

任务按依赖顺序推进。每个阶段必须先满足验收条件，再开始依赖它的模块。

状态说明：

- `[x]` 已实现并通过当前测试。
- `[~]` 有骨架/部分实现，尚不可用于生产。
- `[ ]` 尚未实现或仍需补充验证。
- `[*]` 暂不实施 — 验收条件未满足或无紧迫需求（标注原因）。

## 0. 工程基线

- [x] 建立 `build.zig` 和 Zig 包配置。
- [x] 支持 `HInt` 32/64 位构建选项。
- [x] 实现 `HD`、`HCD` 及单元测试。
- [x] 建立 HCD 与 HiGHS C++ 的对照基准。
- [x] 将基础类型迁移到 `src/foundation/`。
- [x] 建立 `src/matrix/` 模块边界。
- [x] 为当前目录和规划模块补充 README。

验收命令：

```bash
zig build test
zig build test -Dhighs-int-width=w64
zig build bench-hcd -Doptimize=ReleaseFast
```

## 1. 稀疏数据结构

- [x] 定义强类型 `RowId` 和 `ColId`，集中实现 checked 边界转换。
- [x] 实现规范 `SparseVectorView/SparseVector` 和 `SparseVectorBuilder`。
- [x] 实现接受 triplet 的 SoA `MatrixBuilder`。
- [x] 实现构建后冻结的规范 CSC：排序、合并重复项、删除显式零。
- [x] 实现带 revision 检查的按需 CSR cache/view，并提供可信热路径。
- [x] 实现 CSC 的 `Ax`、`A^T y`，并区分 checked/AssumeValid 热路径。
- [x] 添加 CSC 随机矩阵乘法 property test 与结构边界测试。
- [x] 实现矩阵范数、绝对值统计、显式转置和规范行列切片。
- [x] 实现可逆 scaling、行列 permutation 及失败不产生半修改。
- [x] 实现函数式追加/删除行列和动态行 checkpoint/rollback。
- [x] 实现 generation-mark `SparseAccumulator`，服务于聚合、cuts 和 presolve。
- [x] 建立 `bench-matrix` ReleaseFast 专项基准。
- [x] 深度审计 matrix 热路径，为验证、排序、清零、稀疏乘法、格式转换和高精度计算提供显式性能取舍 API。
- [x] 建立 HiGHS/C++ Release 对照 benchmark，并在 `bench/matrix/results.md` 记录环境、方法、结果和复现命令。
- [x] 对照本地 HiGHS 实现 `range/assessment/alphaProductPlusY/productQuad` 等通用接口。
- [x] 将 partitioned row-wise、pricing、collectAj 和 basis update 明确归属 `lp/simplex`。
- [x] 实现 `DenseLU` 参考分解、`SparseLU`（Markowitz + FT 更新）、`SparseBasisBuffers`、`SparseForrestTomlin`。
- [x] 实现 `target_policy.zig` — 编译期架构特化（cache line、vector lanes、prefetch distance、unroll factor）。
- [ ] 建立自动化 HiGHS C++ 矩阵差分测试，覆盖构建、乘法、scaling、切片和格式转换。
- [*] 增加大型稀疏矩阵数据集 benchmark，记录时间、吞吐量、峰值内存和 cache 行为。（当前矩阵运算已在对照基准和 simplex e2e 中间接覆盖，专门的大数据集 benchmark 可在 presolve/MIP 阶段按需补。）
- [*] 根据 benchmark 结果评估 SIMD、prefetch 和并行矩阵乘法；没有实测收益时不引入复杂度。（SIMD 清零内核已就位，并行乘法留待 `parallel` 模块成熟后统一评估。）

## 2. 模型层

- [x] 建立 `src/model/` 模块边界。
- [x] 定义 ProblemClass 分类体系（8 类：LP/MILP/QP/MIQP/QCP/MIQCP/NLP/MINLP）。
- [x] 实现 LinearModel（原 LpModel 重构） — SoA 布局，owning MatrixStore。
- [x] 实现 LinearModelBuilder — 收集坐标、validate、freeze 为 LinearModel。
- [x] 实现 Integrality 枚举（continuous/integer/semi-continuous/semi-integer）。
- [x] 实现 Hessian 存储（triangular/diagonal，1/2 xᵀQx 语义）。
- [x] 实现 QuadraticModel（LinearModel + Hessian + QuadraticConstraint 数组）。
- [x] 实现 ExpressionGraph（DAG，SoA 布局，reference evaluator）。
- [x] 实现 NonlinearModel（QuadraticModel + ExpressionGraph 组合）。
- [x] 实现 CompiledModel（union(enum) 分派：linear/quadratic/nonlinear）。
- [x] 实现 SolverCapability 查询类型。
- [x] 实现 validate 分层框架（Linear → Hessian → Quadratic → Expression → Nonlinear → Compiled）。
- [x] 模型层 w32/w64 独立测试（`zig build test-model`）。
- [x] 实现 `Solution`、`Basis` 求解状态和信息结构（`model/solution.zig`、`model/types.zig`）。
- [x] 实现 primal/dual residual 与 KKT 检查（`model/residual.zig`）。
- [x] 实现 `Model -> CompiledModel` 编译路径（`model/compile_model.zig`）。
- [x] 实现 Gurobi 风格 `Model` API：lazy-update 语义、addVar/addConstr/chgCoeff、属性读写（`"LB"`, `"Obj"`, `"X"`, `"Status"` 等）、generation-checked entity handles（`model/model.zig` + 15 个子模块）。
- [x] 实现 `CompiledModelView` 零拷贝 LP 视图缓存，避免热路径重复编译。
- [x] 实现 `model_solve.zig` — optimize 入口、LP 快速路径、问题类 dispatch。
- [x] 实现 `model_update.zig` — 批量应用 pending changes、矩阵重分配、revision 推进。
- [ ] 建立 HiGHS C API 小模型差分测试（可随 LP 正确性闭环一起完成）。

## 3. LP 正确性闭环

- [x] 实现 dense LU 和 sparse LU 参考分解（`matrix/dense_lu.zig`、`matrix/sparse_lu.zig`）。
- [x] 实现 FTRAN/BTRAN、Forrest-Tomlin 更新、迭代 refinement（`lp/simplex/factorization.zig`）。
- [x] 实现 primal revised simplex 主循环（primal + dual pricing、Harris ratio test、bound flipping）。
- [x] 实现 dual revised simplex 主循环（dual 选行、primal 选列、bound shifting）。
- [x] 实现 Phase I（人工变量 + 成本扰动）、infeasibility/unbounded 检测。
- [x] 实现 basis snapshots、warm start import/export、reoptimization。
- [x] 实现 logical crash basis 和 column/row scaling equilibration。
- [x] 实现 Devex/steepest-edge 加权 pricing（`lp/simplex/pricing.zig`）。
- [x] 实现 anti-cycling（Bland's rule fallback）、退化检测、pivot tolerance。
- [x] 实现奇异 basis 检测与 rank repair。
- [x] 实现 `SimplexEngine` 求解入口、`LpSolveSession`（cold solve / reoptimize）、迭代回调。
- [x] 实现 `bench/simplex/end_to_end_runner.zig` — 与 HiGHS 对照的端到端测试。
- [x] 实现 simplex 专项 benchmark：sparse LU、sparse FT、sparse basis、Suitesparse 对照、session bench。
- [ ] 对 Netlib LP 回归集执行 HiGHS 差分测试，state/objective/value 在统一容差内一致。
- [ ] 补齐 dual Phase I（当前仅在 primal infeasible + dual feasible 路径可用 warm-basis repair；纯 dual infeasible 的逻辑 crash 回退到 primal Phase I）。

## 4. Presolve/Postsolve

- [~] 建立 `src/presolve/` 模块骨架。（`root.zig` 和 `rules/root.zig` 已占位，无实际 reduction）
- [ ] 先定义 reduction 记录和 `PostsolveStack`。
- [ ] 实现 empty row/column、fixed column、singleton row。
- [ ] 实现基础 bound tightening 与 redundant row。
- [ ] 为每条规则实现解和 basis 恢复。

验收：presolve 开关前后的状态和目标一致，恢复后的原模型解满足容差。

> 优先级说明：presolve 是超越开源求解器性能的关键分水岭。当前 simplex 引擎已就位，presolve 补齐后即可在 LP 性能上与 HiGHS 正面竞争。

## 5. SCIP 式组件框架

- [~] 建立 `Stage` 枚举（`framework/stage.zig`）和 `Kind` 枚举（`plugin/kind.zig` — 15 种组件类型）。
- [~] 建立 `src/framework/`、`src/plugin/`、`src/plugins_builtin/` 模块骨架。（均为空壳）
- [ ] 实现 `Registry`、`Scheduler`、`Event` 和 `Services`。
- [ ] 定义统一插件生命周期和显式执行结果类型。
- [ ] 实现 presolver、branching、heuristic 三类最小接口。
- [ ] 添加 priority、frequency、max-depth 和 stage 合法性测试。
- [ ] 保证插件不能直接持有完整 `Solver`。

验收：组件可预测调度，非法阶段调用返回明确错误，热循环不使用动态分发。

> 优先级说明：框架是 MIP 和后续扩展的基础设施。需在 presolve 验证通过后启动。

## 6. 最小 MILP

- [~] 建立 `src/mip/` 模块骨架。（`root.zig` 为空壳）
- [ ] 实现整数可行性、变量域修改和可回滚的 domain stack。
- [ ] 实现 node、node queue、LP relaxation 和 incumbent。
- [ ] 实现 most-infeasible 与 pseudocost branching。
- [ ] 实现 rounding heuristic、cut pool 和一种基础割。
- [ ] 实现 node/time/gap limit。
- [ ] 建立小型 MIPLIB 回归集。

验收：小型 MILP 的状态和目标值与 HiGHS 一致，节点回溯无状态污染。

## 7. I/O 与文件格式

详见 [`src/io/todo.md`](src/io/todo.md)。

### 7a. 基础设施 ✅

- [x] `types.zig` — IoError 错误集、Format/Compression/FileKind、ReadOptions（资源限制 + 原子取消）、WriteOptions、ModelView。
- [x] `format.zig` — 文件后缀检测（.lp/.rlp/.mps/.rew + .gz/.bz2/.zip/.7z/.xz）。
- [x] `input.zig` — FileInput（mmap/buffered/empty 三态，自动选择策略）。
- [x] `output.zig` — Names（确定性名称 + 重复检测）、print/write。
- [x] `string_arena.zig` — 非 owning 字符串池，单次 name interning。
- [x] `builder.zig` — 共享 Builder：LP 列链快速路径、MPS 有序路径、通用排序合并、资源限制。
- [x] `model_data.zig` — 64 字节对齐紧凑存储（属性 + 名称池 + owning CSC），keep_names=false 模式。

### 7b. LP 解析器

- [x] 零复制 streaming lexer + Location 追踪（`lp/lexer.zig`、`lp/token.zig`）。
- [x] Section header 识别（Min/Max/ST/Bounds/Binaries/Generals/Semi-Continuous/Semi-Integer/End）。
- [x] 目标函数解析（线性项 + 系数 + 跨行常数 + 常量项）。
- [x] 约束解析（label/无 label、<=/>=/==、双边范围 a <= expr <= b）。
- [x] Bounds 解析（简单/free/fixed/双边/±inf）。
- [x] 变量类型 section（Binary/Integer/Semi-Continuous/Semi-Integer）。
- [x] 列链直通 CSC 快速路径 + 资源限制 + 协作取消。
- [x] LP writer — 紧凑输出、双边约束、类型 section、确定性名称。
- [ ] **【优先】** 约束续行 — 当前约束必须在同一物理行内，大型生成模型的多行表达式无法导入。
- [ ] 二次表达式 — lexer 已预留 `*` `^` `[` `]`，但 parser 未接入 Hessian 构建。
- [ ] 完整名称规则验证 — 特殊字符、关键字冲突、quoted/escaped names。
- [ ] 结构化诊断 — line/column/错误 token/期望 token/源码片段。
- [ ] 健壮性验证 — fuzz testing、malformed LP corpus、跨求解器差分测试。

### 7c. MPS 解析器

- [x] Free + fixed 字段格式兼容。
- [x] ROWS/COLUMNS/RHS/RANGES/BOUNDS section 完整支持。
- [x] 整数 marker（'MARKER' 'INTORG'/'INTEND'）。
- [x] 多 RHS/RANGES/BOUNDS set（取第一个）。
- [x] 列有序快速路径 + 回退排序路径。
- [x] MPS writer — 完整 OBJSENSE/ROWS/COLUMNS/RHS/RANGES/BOUNDS 输出。
- [ ] OBJNAME 多目标 — header 已识别但未消费。
- [ ] 二次 MPS 扩展 — QUADOBJ / QSECTION / QCMATRIX / QMATRIX。
- [ ] SOS section（SETS type 1/2）。

### 7d. 格式覆盖缺口

- [ ] `.rlp` / `.rew` — 可复用现有 parser/grammar。
- [ ] `.dua` / `.dlp` / `.ilp` / `.opb` — 格式枚举已有，无实现。
- [ ] 压缩输入 — 后缀检测已支持，但 readFile 遇到压缩文件直接返回 UnsupportedCompression。

### 7e. 跨切面

- [ ] 流式/chunked 输入 — 当前解析器需要完整内存缓冲。
- [ ] 直接从 ModelData 构造 public Model — 去除当前 pending-change copy。
- [ ] Hessian/Quadratic 数据通道 — ModelData 当前无 Q 项存储。

## 8. 后续扩展

### 8a. LP 引擎深化

- [ ] 实现 hyper-sparse pricing（大规模稀疏 LP）。
- [ ] 实现 PDLP 首阶方法（`src/pdlp/` 已占位 — 骨架）。
- [ ] 实现 IPM 内点法（`src/ipm/` 已占位 — 骨架）。
- [ ] 实现 IPM→simplex crossover。
- [ ] 自适应 scaling 策略与数值稳定性改进。

### 8b. MIP 深化

- [ ] reliability/strong branching。
- [ ] clique、implication 和 conflict analysis。
- [ ] RINS、diving、feasibility pump。
- [ ] 更多 separator（Gomory、MIR、flow cover）与 constraint handler。
- [ ] 并行 branch-and-bound（依赖 `src/parallel/`）。

### 8c. QP / QCP / NLP

- [~] QP 骨架已占位（`src/qp/root.zig`）。
- [ ] 实现 active-set QP 求解器。
- [ ] 实现 barrier QP 求解器。
- [ ] NLP 表达式求值与自动微分。

### 8d. 并行与 GPU

- [~] 并行框架骨架已占位（`src/parallel/root.zig`）。
- [ ] 实现 task scheduler、线程池、work stealing。
- [ ] 实现并行矩阵-向量乘法与并行 pricing。
- [ ] 评估 GPU 加速方案（CUDA/Vulkan/ROCm），优先在稀疏 LU 和最耗时的 simplex 阶段试验。
- [ ] GPU 加速模块独立于核心求解器，通过 feature flag 控制编译。

### 8e. API 与生态

- [x] C ABI 骨架已占位（`src/bindings/root.zig`、`src/api/root.zig`）。
- [ ] 实现 HiGHS 兼容 C API，使现有 HiGHS 用户可直接替换。
- [ ] Python 绑定（通过 Zig C ABI + ctypes/cffi）。
- [ ] 求解日志、性能 profiling 和可视化工具。

## 当前实际状态总览

```text
模块             状态     代码量      说明
─────────────────────────────────────────────────────────────
foundation/      ✅       ~1K        强类型索引、HD/HCD、double 工具
matrix/          ✅       ~9.5K      CSC/CSR、Builder、LU、FT、ops、memory
model/           ✅      ~11K        Model API、compile、validate、solution、residual
lp/simplex/      ✅       ~4K        SimplexEngine（primal+dual+Phase I+warm）
io/              ✅       ~2.9K      MPS/LP 读写、lexer、builder、ModelData
io: LP parser     ✅        620      线性目标/约束/bounds/types、列链 CSC
io: MPS parser    ✅        492      ROWS/COLUMNS/RHS/RANGES/BOUNDS
io: LP 续行        ⬜         —       约束跨行、二次表达式
io: 压缩/其他      ⬜         —       .gz/.bz2、.dua/.dlp/.ilp/.opb
solver/          ✅       ~140       求解分发、LpSolveSession
presolve/        🟡         15       骨架占位（rules 为空）
framework/       🟡         20       Stage 枚举、空壳
plugin/          🟡         25       Kind 枚举（15 种组件类型）
plugins_builtin/ 🟡          6       空壳
mip/             🟡          9       空壳
pdlp/            🟡          8       空壳
ipm/             🟡          8       空壳
qp/              🟡          8       空壳
parallel/        🟡          8       空壳
nla/             🟡          9       空壳（因子分解已在 matrix/ 和 simplex/ 实现）
analysis/        🟡          1       空壳
diagnostics/     🟡          1       空壳
bindings/        🟡          1       空壳（C ABI 预留）
api/             🟡          1       空壳
─────────────────────────────────────────────────────────────
✅  生产就绪    🟡  骨架/占位    ⬜  未开始
```

## 当前建议执行顺序

1. ✅ 完成稀疏数据结构、模型 IR、Gurobi 风格 Model API、MPS/LP 读写。
2. ✅ 完成 revised simplex 引擎（primal/dual/Phase I/warm start/reoptimize）。
3. [ ] **【当前优先级】** Netlib 回归测试 + dual Phase I 补齐 → 确认 LP 正确性闭环。
4. [ ] 实现 presolve/postsolve（empty row/col、fixed col、singleton row、bound tightening）。
5. [ ] 实现 SCIP 式组件框架（Registry、Scheduler、Event、最小插件接口）。
6. [ ] 实现最小 branch-and-bound（node、domain、branching、heuristic、cut pool）。
7. [ ] MIP 深化（strong branching、conflict analysis、更多 cuts、propagation）。
8. [ ] LP 引擎深化（hyper-sparse pricing、IPM、crossover、自适应 scaling）。
9. [ ] QP 支持 → 并行框架 → GPU 加速。

> 隐式目标：每一步都必须在对标测试中达到或超越开源求解器的性能水平，并为后续对标商业求解器留出架构扩展空间。GPU 加速作为远期差异化能力，在 CPU 路径稳定后启动。
