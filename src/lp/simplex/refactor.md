# Simplex refactor tracking

不阻塞功能开发。以下各项在 dual simplex Phase I 修复完成、全路径 A/B 通过后
再执行。每项完成后勾选。

## Engine split

engine.zig 当前 240k 行，承载了全部 simplex 变体的主循环、pricing 分派、
退化处理、Devex/DSE 权重、warm-start 修复、trace 诊断。目标：engine.zig 只
保留 solveProblem、solvePrimal、solveDual 三个入口函数 + 公共基础设施；
各变体独立为子模块。

- [ ] **dual_pivot.zig**：从 engine.zig 提取 `solveDual` 和 `solveDualPhaseOne`
  主循环。包含 `chooseDualLeavingRow`、`computeDualTableauRow`、
  `updateReducedCostsAfterDualPivot`、`repairWarmBasisWithDual`
- [ ] **dual_edge_weight.zig**：从 engine.zig 提取 dual Devex/DSE 权重逻辑。
  包含 `beginDualEdgeWeightPhase`、`switchDualDseToDevex`、
  `ensureExactDualEdgeWeights`、`updateDualSteepestEdgeWeights`、
  `updateDualDevexWeights`。注意：primal Devex 权重也在 engine.zig 中，
  可一并提取到同一文件或 `devex_weight.zig`
- [ ] **dual_candidate.zig**：从 engine.zig 提取 dual 候选集维护逻辑。
  包含 `dualCandidateScore`、`bestDualCandidate`、`rebuildDualCandidateList`
- [ ] **dual_diagnostic.zig**：从 engine.zig 提取 `recordDualPhaseOneNoEntering`
  和相关 trace 事件类型

## Pricing consolidation

dual leaving 选择当前分散在三处：

- engine.zig `chooseDualLeavingRow()`（分派：Dantzig / weighted / Bland）
- pricing.zig `chooseDualLeavingWeighted()`（实际计算）
- engine.zig `chooseDualLeavingBland()`（Bland fallback）

- [ ] 将 `chooseDualLeavingBland` 移入 pricing.zig，统一接口为
  `pricing.chooseDualLeaving(engine, problem)` 一次分派

## engine.zig size target

- 当前：~240k 行，~30 dual 函数内联在 engine 中
- 目标：engine.zig < 80k 行。primal 循环、dual 循环、Devex 权重、退化策略
  各自独立文件，engine 保留 solve/dispatch/finish/infeasibility 公共入口
