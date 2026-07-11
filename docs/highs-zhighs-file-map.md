# HiGHS 与 zhighs 功能文件对应表

本文用于在移植、测试和代码审查时快速定位“两个项目中承担相同功能的文件”。

## 基线与状态说明

- HiGHS 源码：`/home/godv/documents/codefiles/cppfiles/HiGHS`
- 本次核对版本：`v1.14.0-4-gdcc25308d8-dirty`
- zhighs 源码：本仓库 `src/`
- 对应关系按功能划分，不表示逐行翻译，也不要求一对一文件布局。

| 标记 | 含义 |
|---|---|
| ✅ | Zig 文件中已有实际实现和测试 |
| 🧱 | Zig 模块入口已经创建，具体算法尚未实现 |
| 📝 | 建议的目标文件，后续实现时创建 |
| ➖ | 不计划直接移植，功能由其他层吸收或不属于首期范围 |

## 1. 公共 API 与求解编排

| 功能 | HiGHS C++ 文件 | zhighs Zig 文件 | 状态与说明 |
|---|---|---|---|
| 顶层求解器对象 | `highs/Highs.h` | `src/api/root.zig`, `src/solver/root.zig` | 🧱 拆分公共 facade 与内部编排，避免形成单个上帝类。 |
| 顶层对象核心实现 | `highs/lp_data/Highs.cpp` | `src/solver/solver.zig` | 📝 负责模型生命周期、presolve、求解、postsolve 和结果发布。 |
| 模型修改和查询 API | `highs/lp_data/HighsInterface.cpp` | `src/api/solver.zig`, `src/model/lp_model_builder.zig` | 📝 待真正的 LP 数据结构完成后实现；API 做参数检查，builder 负责模型变更。 |
| LP 求解器路由 | `highs/lp_data/HighsSolve.h/.cpp` | `src/solver/solve_lp.zig` | 📝 在 Simplex、IPM、PDLP 和无约束快速路径之间选择。 |
| LP 求解调用上下文 | `highs/lp_data/HighsLpSolverObject.h` | `src/lp/context.zig` | 📝 用显式 context 传递模型、basis、solution、info、options 和 timer。 |
| API 返回状态 | `highs/lp_data/HighsStatus.h/.cpp` | `src/api/status.zig` | 📝 区分调用状态、模型状态、presolve/postsolve 状态。 |
| 选项定义与校验 | `highs/lp_data/HighsOptions.h/.cpp` | `src/api/options.zig` | 📝 使用 typed options；解析文件由 `io/` 负责。 |
| 求解信息 | `highs/lp_data/HighsInfo.h/.cpp` | `src/api/info.zig`, `src/diagnostics/statistics.zig` | 📝 稳定公开字段与内部统计分开。 |
| 用户回调 | `highs/lp_data/HighsCallback.h/.cpp`, `HighsCallbackStruct.h` | `src/api/callback.zig`, `src/framework/event.zig` | 📝 API 回调数据与内部事件分开。 |
| 多目标调度 | `Highs::multiobjectiveSolve`（`Highs.cpp`） | `src/solver/multiobjective.zig` | 📝 多次求解、优先级、容差和 objective fixing。 |
| 旧 API 兼容 | `highs/lp_data/HighsDeprecated.cpp` | `src/api/compat.zig` | ➖ 只在 zhighs 发布稳定 API 后按需提供。 |

## 2. 基础类型、容器与数值工具

| 功能 | HiGHS C++ 文件 | zhighs Zig 文件 | 状态与说明 |
|---|---|---|---|
| 可配置整数类型 | `highs/util/HighsInt.h` | `src/foundation/int.zig` | ✅ 已支持 32/64 位构建配置。 |
| 强类型行列标识 | HiGHS 中通常直接使用 `HighsInt` 和 `-1` 哨兵 | `src/foundation/index.zig` | ✅ 使用 `comptime` 工厂生成 `RowId/ColId`，并以 `maxInt(HUInt)` 实现同宽度紧凑 OptionalId。 |
| 双双精度累加 | `highs/util/HighsCDouble.h` | `src/foundation/double.zig` | ✅ 已实现并有单元测试及 C++ 对照基准。 |
| 常量和公共枚举 | `highs/lp_data/HConst.h`, `simplex/SimplexConst.h` | `src/foundation/constants.zig`, `src/lp/simplex/constants.zig` | 📝 通用常量与算法常量分层。 |
| 排序工具 | `highs/util/HighsSort.h/.cpp` | `src/foundation/sort.zig` | 📝 仅保留标准库不足的索引排序。 |
| 随机数 | `highs/util/HighsRandom.h` | `src/foundation/random.zig` | 📝 必须支持确定性 seed 和可复现测试。 |
| 哈希和集合 | `HighsHash*`, `HSet*`, `HighsRbTree.h`, `HighsSplay.h` | `src/foundation/collections/` | 📝 先评估 Zig 标准库，再移植专用结构。 |
| 并查集 | `highs/util/HighsDisjointSets.h` | `src/foundation/disjoint_sets.zig` | 📝 供 clique、冲突和图算法使用。 |
| 通用数据栈 | `highs/util/HighsDataStack.h` | `src/foundation/data_stack.zig` | 📝 与 presolve 的语义化 postsolve stack 区分。 |
| 计时 | `HighsTimer.h`, `FactorTimer.h`, `MipTimer.h` | `src/foundation/timer.zig`, `src/diagnostics/timing.zig` | 📝 时钟原语在 foundation，指标命名在 diagnostics。 |
| 字符串工具 | `highs/util/stringutil.h/.cpp` | `src/io/text.zig` | 📝 归入输入输出层。 |
| 通用数学工具 | `highs/util/HighsUtils.h/.cpp` | `src/foundation/math.zig` | 📝 禁止重新形成无边界的 util 大杂烩。 |

## 3. 模型、矩阵、解与 Basis

| 功能 | HiGHS C++ 文件 | zhighs Zig 文件 | 状态与说明 |
|---|---|---|---|
| LP 数据 | `highs/lp_data/HighsLp.h/.cpp` | `src/model/lp_model.zig` | 📝 必须同时包含目标、行列界、integrality 和矩阵；尚未创建。 |
| 通用模型容器 | `highs/model/HighsModel.h/.cpp` | `src/model/model.zig` | 📝 未来组合 `LpModel`、Hessian 和多目标；当前尚未创建，避免与 LP 数据结构混名。 |
| Hessian | `HighsHessian.h/.cpp`, `HighsHessianUtils.*` | `src/model/hessian.zig` | 📝 只接受凸 QP 所需的规范存储。 |
| 多线性目标 | `HighsLinearObjective`（`HStruct.h`） | `src/model/objective.zig` | 📝 保存 weight、priority 和容差。 |
| 解结构 | `HighsSolution`（`HStruct.h`）, `HighsSolution.*` | `src/model/solution.zig` | 📝 primal/dual validity 和行列值分开表达。 |
| Basis | `HighsBasis`（`HStruct.h`） | `src/model/basis.zig` | 📝 保存行列 basis status，不持有 factorization。 |
| Scaling 数据 | `HighsScale`（`HStruct.h`） | `src/matrix/scaling.zig` | ✅ 正因子校验、事务式应用/移除和 max equilibration。 |
| 模型合法性检查 | `HighsLpUtils.*`, `HighsModelUtils.*` | `src/model/validate.zig` | 📝 builder freeze 时执行。 |
| 半连续/半整数 reformulation | `HighsLpUtils.*`, `HighsLpMods` | `src/model/reformulation.zig` | 📝 形成可逆 transformation，不藏在顶层 run 内。 |
| 稀疏矩阵 | `highs/util/HighsSparseMatrix.h/.cpp` | `src/matrix/csc.zig`, `src/matrix/csr_view.zig` | ✅ 权威 CSC、列 view、`Ax/A^T y`，以及带 revision 失效检查的按需 CSR cache/view。 |
| 矩阵构建 | `HighsSparseMatrix::addRows/addCols` 等 | `src/matrix/builder.zig` | ✅ SoA triplet 输入，按 `(col,row,ordinal)` 排序，稳定合并、去零并冻结为 CSC。 |
| 规范稀疏行/列及构建 | `HighsMatrixSlice.h`、矩阵行列访问逻辑 | `src/matrix/sparse_vector.zig`, `sparse_vector_builder.zig` | ✅ sorted unique `(Id, f64)`、借用 view、合并重复和去零。 |
| 矩阵算法 | `HighsMatrixUtils.*` | `src/matrix/ops.zig`, `transpose.zig`, `slice.zig`, `permutation.zig`, `edit.zig` | ✅ `Ax/A^T y`、稳定范数、转置、切片、置换及函数式结构编辑。 |
| 矩阵乘加与高精度乘法 | `HighsSparseMatrix::alphaProductPlusY/productQuad/productTransposeQuad` | `src/matrix/ops.zig` | ✅ `y += alpha*A*x`、转置版本及基于 `HCD` 的高精度乘法。 |
| 矩阵数值评估 | `HighsSparseMatrix::range/assessSmallValues/hasLargeValue` | `src/matrix/ops.zig` | ✅ 绝对值范围、small/large 计数和阈值查询；日志由上层 diagnostics 负责。 |
| Simplex pricing 矩阵视图 | `createRowwisePartitioned`, `priceByRow/Column`, `collectAj`, `update` | `src/lp/simplex/pricing_matrix.zig` | 📝 依赖 basis partition、`HVector` 和 pricing 策略，不进入通用 matrix。 |
| Simplex 稀疏工作向量 | `HVector.h`, `HVectorBase.h/.cpp` | `src/nla/sparse_work_vector.zig` | 📝 后续支持 dense values + active index 集合，不与矩阵切片混用。 |
| 稀疏向量和 | `HighsSparseVectorSum.h` | `src/matrix/sparse_sum.zig` | ✅ dense value + generation mark + active ID，用于割、聚合和 O(1) 逻辑清空。 |
| 动态 MIP 行矩阵 | `mip/HighsDynamicRowMatrix.*` | `src/matrix/dynamic_rows.zig` | ✅ append-only 连续行存储、checkpoint/rollback 和批量合并；cuts 不直接修改基础 CSC。 |

## 4. 数值线性代数与 Revised Simplex

| 功能 | HiGHS C++ 文件 | zhighs Zig 文件 | 状态与说明 |
|---|---|---|---|
| Basis factorization 主体 | `util/HFactor.h/.cpp` | `src/nla/factorization.zig` | 🧱 `nla` 入口已建，具体实现待完成。 |
| Refactor / update | `HFactorRefactor.cpp`, `HFactorExtend.cpp` | `src/nla/refactor.zig`, `src/nla/update.zig` | 📝 Forrest–Tomlin/产品形式按验证顺序实现。 |
| Factor 辅助与调试 | `HFactorUtils.cpp`, `HFactorDebug.*`, `HFactorConst.h` | `src/nla/factor_utils.zig`, `src/diagnostics/factor_check.zig` | 📝 算法与调试检查分开。 |
| Simplex 总控状态 | `simplex/HEkk.h/.cpp`, `HEkkControl.cpp` | `src/lp/simplex/solver.zig`, `state.zig`, `control.zig` | 🧱 Simplex 模块入口已建。 |
| Simplex 公共接口 | `HSimplex.h/.cpp`, `HEkkInterface.cpp` | `src/lp/simplex/interface.zig` | 📝 供顶层 LP router 和 MIP relaxation 使用。 |
| Primal simplex | `HEkkPrimal.h/.cpp` | `src/lp/simplex/primal.zig` | 📝 Phase I/II、entering/leaving 和更新。 |
| Dual simplex | `HEkkDual.h/.cpp`, `HEkkDualMulti.cpp` | `src/lp/simplex/dual.zig` | 📝 MIP reoptimization 主线。 |
| Dual RHS | `HEkkDualRHS.h/.cpp` | `src/lp/simplex/dual_rhs.zig` | 📝 管理 primal infeasibility 与 leaving row。 |
| Dual row / pricing | `HEkkDualRow.h/.cpp` | `src/lp/simplex/pricing.zig` | 📝 Devex、DSE 和候选选择。 |
| NLA 桥接 | `HSimplexNla.h/.cpp`, `HSimplexNla*.cpp` | `src/lp/simplex/nla.zig` | 📝 管理 factor、FTRAN/BTRAN、freeze 和 product-form。 |
| Simplex 数据结构 | `SimplexStruct.h` | `src/lp/simplex/state.zig` | 📝 工作 cost/bounds/value、basis map、edge weights 和 rays。 |
| Simplex 计时与统计 | `SimplexTimer.h`, `HighsSimplexAnalysis.*` | `src/diagnostics/simplex_stats.zig` | 📝 不放入核心 pivot 逻辑。 |
| Simplex 日志 | `HSimplexReport.h/.cpp` | `src/diagnostics/simplex_report.zig` | 📝 只读 solver snapshot。 |
| 应用入口 | `simplex/HApp.h` | `src/lp/simplex/solve.zig` | 📝 `solveLpSimplex` 对应入口。 |

## 5. Presolve 与 Postsolve

| 功能 | HiGHS C++ 文件 | zhighs Zig 文件 | 状态与说明 |
|---|---|---|---|
| Presolve 主循环 | `presolve/HPresolve.h/.cpp` | `src/presolve/engine.zig` | 🧱 模块入口已建。 |
| Presolve 分析 | `HPresolveAnalysis.h/.cpp` | `src/diagnostics/presolve_stats.zig` | 📝 规则效果与耗时统计。 |
| Postsolve 恢复栈 | `HighsPostsolveStack.h/.cpp` | `src/presolve/postsolve_stack.zig` | 📝 reduction 必须以 LIFO 方式恢复解、basis 和证书。 |
| 组件包装 | `PresolveComponent.h/.cpp` | `src/plugin/presolver.zig`, `src/framework/registry.zig` | 📝 用统一插件契约替代单独组件基类。 |
| 具体 reduction | `HPresolve.cpp` 中各规则 | `src/presolve/rules/*.zig` | 🧱 rules 入口已建，规则按功能拆文件。 |
| Symmetry | `HighsSymmetry.h/.cpp` | `src/plugins_builtin/presolve/symmetry.zig` | 📝 可选高级组件。 |
| Crash | `ICrash.h/.cpp`, `ICrashUtil.*`, `ICrashX.*` | `src/lp/simplex/crash.zig` | 📝 从 presolve 目录移到 simplex 初始化职责。 |

## 6. MIP Branch-and-Cut

| 功能 | HiGHS C++ 文件 | zhighs Zig 文件 | 状态与说明 |
|---|---|---|---|
| MIP 顶层 | `mip/HighsMipSolver.h/.cpp` | `src/mip/solver.zig` | 🧱 `mip` 入口已建。 |
| MIP 共享状态 | `HighsMipSolverData.h/.cpp` | `src/mip/state.zig` | 📝 按所有权继续拆分，避免共享大对象。 |
| 搜索循环 | `HighsSearch.h/.cpp` | `src/mip/search.zig` | 📝 节点处理与组件调度。 |
| 节点队列 | `HighsNodeQueue.h/.cpp` | `src/mip/node_queue.zig` | 📝 node selection 可通过插件策略替换。 |
| Domain 与回滚 | `HighsDomain.h/.cpp`, `HighsDomainChange.h` | `src/mip/domain.zig`, `domain_stack.zig` | 📝 所有边界修改必须可回溯。 |
| LP relaxation | `HighsLpRelaxation.h/.cpp` | `src/mip/relaxation.zig` | 📝 持有节点 LP 状态和 simplex warm start。 |
| 变换后的 LP | `HighsTransformedLp.h/.cpp` | `src/mip/transformed_lp.zig` | 📝 MIP 内部变量/行映射。 |
| Cut pool | `HighsCutPool.h/.cpp` | `src/mip/cut_pool.zig` | 📝 管理有效域、age、重复和删除。 |
| Cut generation | `HighsCutGeneration.h/.cpp` | `src/plugins_builtin/separator/cut_generation.zig` | 📝 具体 cut 算法作为 separator。 |
| Separation 调度 | `HighsSeparation.*`, `HighsSeparator.*` | `src/plugin/separator.zig`, `src/framework/scheduler.zig` | 📝 接口与调度分开。 |
| Tableau cuts | `HighsTableauSeparator.*` | `src/plugins_builtin/separator/tableau.zig` | 📝 Gomory/MIR 等 tableau 派生割。 |
| Path / mod-k cuts | `HighsPathSeparator.*`, `HighsModkSeparator.*` | `src/plugins_builtin/separator/path.zig`, `mod_k.zig` | 📝 延后实现。 |
| Primal heuristics | `HighsPrimalHeuristics.*` | `src/plugins_builtin/heuristic/*.zig` | 📝 rounding、diving、RINS/RENS 等拆分。 |
| Feasibility jump | `HighsFeasibilityJump.cpp`, `feasibilityjump.hh` | `src/plugins_builtin/heuristic/feasibility_jump.zig` | 📝 独立启发式组件。 |
| Pseudocost | `HighsPseudocost.h/.cpp` | `src/mip/pseudocost.zig` | 📝 数据归 MIP，branching policy 归插件。 |
| Implications | `HighsImplications.h/.cpp` | `src/mip/implications.zig` | 📝 服务 propagation、cuts 和 branching。 |
| Clique table | `HighsCliqueTable.h/.cpp` | `src/mip/clique_table.zig` | 📝 0/1 变量逻辑结构。 |
| Conflict pool | `HighsConflictPool.h/.cpp` | `src/mip/conflict_pool.zig` | 📝 冲突约束、age 和传播。 |
| LP aggregation | `HighsLpAggregator.h/.cpp` | `src/mip/aggregation.zig` | 📝 生成聚合行和 cut 输入。 |
| Reduced-cost fixing | `HighsRedcostFixing.h/.cpp` | `src/plugins_builtin/propagator/reduced_cost.zig` | 📝 作为 propagation 组件。 |
| 目标管理 | `HighsObjectiveFunction.h/.cpp` | `src/mip/objective.zig` | 📝 incumbent、cutoff 与 objective integral scale。 |
| MIP 调试与统计 | `HighsMipAnalysis.*`, `HighsDebugSol.*` | `src/diagnostics/mip_stats.zig`, `mip_check.zig` | 📝 与搜索状态分离。 |

## 7. IPM、PDLP 与 QP

| 功能 | HiGHS C++ 文件 | zhighs Zig 文件 | 状态与说明 |
|---|---|---|---|
| IPX 包装 | `ipm/IpxWrapper.h/.cpp`, `IpxSolution.h` | `src/ipm/ipx_adapter.zig`, `solution.zig` | 🧱 IPM 模块入口已建，延后实现。 |
| IPX 算法 | `ipm/ipx/*` | `src/ipm/ipx/` | 📝 KKT、迭代、crossover 和线性系统。 |
| HiPO | `ipm/hipo/*` | `src/ipm/hipo/` | 📝 作为可选 backend。 |
| IPM LU | `ipm/basiclu/*` | `src/nla/basic_lu/` 或外部依赖 | 📝 先评估复用统一 NLA。 |
| cuPDLP 适配 | `pdlp/CupdlpWrapper.*` | `src/pdlp/cupdlp_adapter.zig` | 📝 可选外部/GPU backend。 |
| HiPDLP 适配 | `pdlp/HiPdlpWrapper.*` | `src/pdlp/solve.zig` | 🧱 PDLP 模块入口已建。 |
| HiPDLP 算法 | `pdlp/hipdlp/*` | `src/pdlp/pdhg.zig`, `restart.zig`, `scaling.zig` | 📝 延后于 revised simplex。 |
| QP 顶层 | `qpsolver/quass.h/.cpp`, `a_quass.*` | `src/qp/solver.zig` | 🧱 QP 模块入口已建。 |
| Active-set | `qpsolver/a_asm.*` | `src/qp/active_set.zig` | 📝 凸 QP 工作集主循环。 |
| QP basis/factor | `basis.*`, `factor.hpp` | `src/qp/basis.zig`, `src/nla/qp_factor.zig` | 📝 复用统一 NLA 契约。 |
| 梯度与 reduced gradient | `gradient.hpp`, `reducedgradient.hpp`, `reducedcosts.hpp` | `src/qp/gradient.zig` | 📝 QP 一阶最优性量。 |
| QP pricing/ratio | `*pricing.hpp`, `ratiotest.*` | `src/qp/pricing.zig`, `ratio_test.zig` | 📝 策略可静态组合。 |
| QP scaling/perturbation | `scaling.*`, `perturbation.*` | `src/qp/scaling.zig`, `perturbation.zig` | 📝 保持与 LP scaling 语义分离。 |

## 8. 分析、I/O、接口与并行

| 功能 | HiGHS C++ 文件 | zhighs Zig 文件 | 状态与说明 |
|---|---|---|---|
| KKT 检查 | `HighsSolution.*`, `test_kkt/*` | `src/analysis/kkt.zig`, `test/differential/kkt.zig` | 🧱 analysis 入口已建。 |
| IIS | `HighsIis.h/.cpp`, `lp_data/Iis.md` | `src/analysis/iis.zig` | 📝 统一 LP/MIP IIS 结果。 |
| Ranging | `HighsRanging.h/.cpp` | `src/analysis/ranging.zig` | 📝 依赖有效 simplex basis。 |
| Ray/certificate | `Highs.cpp`, `HighsSolution.*`, `SimplexStruct.h` | `src/analysis/certificate.zig` | 📝 primal ray、dual ray 和 infeasibility certificate。 |
| MPS | `io/FilereaderMps.*`, `HMPSIO.*`, `HMpsFF.*` | `src/io/mps.zig` | 🧱 I/O 入口已建。 |
| LP 格式 | `io/FilereaderLp.*`, `io/filereaderlp/*` | `src/io/lp.zig` | 📝 parser 写入 `LpModelBuilder`。 |
| 文件读取分派 | `io/Filereader.*` | `src/io/reader.zig` | 📝 可通过 reader plugin 扩展。 |
| 选项文件 | `io/LoadOptions.*` | `src/io/options.zig` | 📝 解析结果交给 typed options 校验。 |
| 日志 | `io/HighsIO.*` | `src/diagnostics/logger.zig` | 🧱 diagnostics 入口已建。 |
| C API | `interfaces/highs_c_api.h/.cpp` | `src/bindings/c_api.zig` | 🧱 bindings 入口已建。 |
| C#/Fortran | `interfaces/highs_csharp_api.cs`, `highs_fortran_api.f90` | 独立 bindings 包 | 📝 基于稳定 C ABI，不进入核心库。 |
| Python | `highspy/*`, `highs_bindings.cpp` | 独立 Python binding | 📝 不直接暴露 Zig 内部布局。 |
| 任务执行器 | `parallel/HighsTaskExecutor.*`, `HighsTask.h` | `src/parallel/executor.zig`, `task.zig` | 🧱 parallel 入口已建。 |
| 同步原语 | `HighsMutex.h`, `HighsSpinMutex.h`, `HighsBinarySemaphore.h` | `src/parallel/sync.zig` | 📝 优先采用 Zig 标准库。 |
| 调度数据结构 | `HighsSplitDeque.h`, `HighsCombinable.h` | `src/parallel/work_queue.zig` | 📝 保证串行 fallback 和可取消性。 |

## 9. SCIP 式组件化新增文件

这些文件在 HiGHS 中通常以内嵌策略或具体类存在，zhighs 将其抽成稳定边界：

| zhighs 文件 | 吸收的 HiGHS 功能 | 状态 |
|---|---|---|
| `src/framework/stage.zig` | `Highs::run()` 隐含生命周期和合法调用时机 | ✅ 已定义初始阶段枚举。 |
| `src/plugin/kind.zig` | MIP 中 separator/heuristic/branching 等散布组件 | ✅ 已定义项目范围内组件种类。 |
| `src/framework/registry.zig` | 各组件的集中持有和生命周期 | 📝 |
| `src/framework/scheduler.zig` | priority、frequency、depth、timing 调度 | 📝 |
| `src/framework/services.zig` | 对组件暴露受限能力，阻止反向依赖 Solver | 📝 |
| `src/framework/event.zig` | callback、变量域、节点、incumbent 等事件 | 📝 |
| `src/plugin/*.zig` | 每种组件的 context、metadata、callback 和 result | 📝 |
| `src/plugins_builtin/**` | HiGHS 内置 cuts、heuristics、branching、propagation | 🧱 根模块已建，具体实现待创建。 |

## 维护规则

1. 移植一个 HiGHS 文件前，先在本表确认它应落入哪个职责层。
2. 一个 HiGHS 文件拆到多个 Zig 文件时，应在同一行或相邻行记录全部目标。
3. Zig 文件实现后，把状态从 📝/🧱 更新为 ✅，并附对应测试位置。
4. 如果上游 HiGHS 重命名或拆分文件，记录核对的 commit，不静默改表。
5. 本表描述功能等价，不承诺 ABI、内部数据布局或逐行代码等价。
