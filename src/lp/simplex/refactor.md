# Simplex refactor tracking

不阻塞功能开发。以下各项在 dual simplex Phase I 修复完成、全路径 A/B 通过后
再执行。每项完成后勾选。

## Engine split

engine.zig 拆分前 4652 行，承载了全部 simplex 变体的主循环、pricing 分派、
退化处理、Devex/DSE 权重、warm-start 修复、trace 诊断。目标：engine.zig 只
保留 solveProblem、solvePrimal、solveDual 三个入口函数 + 公共基础设施；
各变体独立为子模块。

- [x] **dual_pivot.zig**（实际落地为 `engine_dual.zig`）：从 engine.zig 提取
  `solveDual` 和 `solveDualPhaseOne` 主循环。包含 `chooseDualLeavingRow`、
  `computeDualTableauRow`、`updateReducedCostsAfterDualPivot`、
  `repairWarmBasisWithDual`，以及 bound flip 组（`applyBoundFlips`、
  `boundFlipResidualAcceptable`、`accumulateBoundFlipRhs`）和
  `buildDualPhaseOneCosts`、`deterministicUnit`
- [x] **dual_edge_weight.zig**（实际落地为 `engine_dual_edge_weight.zig`）：
  从 engine.zig 提取 dual Devex/DSE 权重逻辑。包含
  `beginDualEdgeWeightPhase`、`switchDualDseToDevex`、
  `ensureExactDualEdgeWeights`、`updateDualSteepestEdgeWeights`、
  `updateDualDevexWeights` 和 `DualDseFallbackReason`。primal Devex 权重
  随 primal 主循环提取到了 `engine_primal.zig`
  （`initializePrimalDevexFramework`、`updatePrimalDevexFramework`、
  `updateLegacyDevexWeights`）
- [x] **dual_candidate.zig**（实际落地为 `engine_dual_candidate.zig`）：
  从 engine.zig 提取 dual 候选集维护逻辑。包含 `dualCandidate`、
  `dualCandidateScore`、`bestDualCandidate`、`rebuildDualCandidateList`
- [x] **dual_diagnostic.zig**：`recordDualPhaseOneNoEntering` 随 dual Phase I
  一并进入 `engine_dual.zig`；trace 事件类型作为公共 API 类型保留在
  engine.zig

拆分采用 `src/model/` 的先例：engine.zig 保留全部 pub 类型、
`SimplexEngine` 字段、`init`/`deinit`/`requestedBytes`/`solve` 占位、
`solveProblem`/`reoptimizeProblem` 入口，方法以
`pub const foo = @import("engine_foo.zig").foo;` 重导出。其余落地文件：
`engine_setup.zig`（问题存储/缩放/basis 安装）、`engine_basis.zig`
（basis 导入/导出/修复）、`engine_primal.zig`（primal 主循环 + primal
Phase I + primal Devex）、`engine_pivot.zig`（pivot 机制、reduced cost
更新、factorization 生命周期）、`engine_degeneracy.zig`（退化处理驱动）、
`engine_progress.zig`（求解控制、统计计时、密度观测）、
`engine_finish.zig`（收尾与解校验）。

## Pricing consolidation

dual leaving 选择当前分散在三处：

- engine_dual.zig `chooseDualLeavingRow()`（分派：Dantzig / weighted / Bland）
- pricing.zig `chooseDualLeavingWeighted()`（实际计算）
- engine_dual.zig `chooseDualLeavingBland()`（Bland fallback）

- [x] 将 `chooseDualLeavingBland` 移入 pricing.zig，统一接口为
  `pricing.chooseDualLeaving(engine)` 一次分派（engine 通过 anytype 传入，
  避免循环依赖；Bland 和 weighted 分支统一定义在 pricing.zig 中，hyper-sparse
  候选通过 engine 方法回调）

## engine.zig size target

- 拆分前：4652 行，~30 dual 函数内联在 engine 中
- 拆分后（实际结果）：engine.zig 855 行；engine_primal.zig 988 行、
  engine_dual.zig 847 行、engine_pivot.zig 559 行、
  engine_basis.zig 368 行、engine_progress.zig 356 行、
  engine_setup.zig 335 行、engine_degeneracy.zig 229 行、
  engine_dual_edge_weight.zig 203 行、engine_finish.zig 195 行、
  engine_dual_candidate.zig 103 行。
  primal 循环、dual 循环、Devex/DSE 权重、退化策略各自独立文件，engine
  保留类型定义、字段、入口与重导出
