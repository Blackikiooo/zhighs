# Dual simplex 对齐路线

基准源码：本地 HiGHS `de09bbad9fb7c5d39a1a464a7641bbb5531c6e9d`。

目标：先让 Zig 的 serial/plain dual simplex 在状态机、工作数组、阶段切换、选主元、更新和
cleanup 上与 `HEkkDual` 源码同构；验证口径为相同模型、presolve off、相同 simplex 选项下的
committed simplex iterations。稀疏矩阵、primal simplex 和 tasks/multi dual 不在本路线首轮范围。

## 1. HiGHS serial dual simplex 算法路径

1. `HEkkDual::solve`
   - `initialiseSolve`；无扰动 cost → `computeDual` → 原始 dual infeasibility。
   - 计算 `force_phase2/near_optimal/perturb_costs`；重新初始化 work cost。
   - 初始化 DSE/Devex；扰动后重新 `computeDual`，fixed/boxed flip 后决定 Phase I/II。
   - 在同一个 major phase loop 内运行 Phase I、Phase II、cleanup 和 phase re-entry。
2. `solvePhase1/solvePhase2`
   - 安装对应 work bounds/value/move。
   - outer loop: `rebuild`；inner loop: `iterate`，由 `rebuild_reason` 统一中断和恢复。
3. `rebuild`
   - 按 reason 决定 reinvert；随后固定顺序执行 `computeDual → correctDualInfeasibilities →
     computePrimal → rebuild primal infeasibility list → dual objective`。
4. `iterate`
   - `chooseRow` → `chooseColumn` → bad-basis guard。
   - `FTRAN(BFRT)` → `FTRAN(pivot)` → `FTRAN(DSE)` → row/column pivot verify。
   - `updateDual` → `updatePrimal + weights` → `updatePivots/factor/matrix/infeas-list`。
5. `chooseRow`
   - weighted CHUZR；`BTRAN(e_p)`；exact DSE 只拒绝低估超过 4 倍的权重并重选。
6. `chooseColumn`
   - CHUZC0 free move；CHUZC1 pack row；CHUZC2 candidate；CHUZC3 large-step；
     CHUZC4 breakpoint groups/large-alpha/strict prior-group flips。
   - 小 pivot 时 iterative-refinement BTRAN + quad PRICE；仍小时删除 pivot 并重跑 CHUZC。
7. `updateDual/updatePrimal/updatePivots`
   - 增量维护 `workDual/workShift/workValue/baseValue/infeasibility list/DSE`；basis 更新最后提交。
8. Phase-I exit 和 cleanup
   - fresh rebuild 后按 Phase-I dual objective、原始 LP dual infeasibility判断。
   - remove/re-perturb cost、free-variable shift、Phase I/II re-entry 都在同一 work state 上完成。

## 2. Zig 源码差异矩阵

| 路径 | Zig 现状 | 判定 | 必须对齐的内容 |
|---|---|---|---|
| solve controller | shifted dual、独立 Phase-I、snapshot retry、primal fallback 分支拼接 | 未对齐 | 单一 `solve_phase + rebuild_reason` major loop |
| work state | model bounds 与 `DualPhaseOneWorkspace` 临时数组混用 | 未对齐 | 持久 `workCost/workShift/workLower/workUpper/workRange/workValue/nonbasicMove` |
| cost init | 扰动公式、启用时机、scaling 与 HiGHS 不同 | 未对齐 | 两次 cost init、near-optimal policy、精确 RNG、cost scale |
| phase selection | 预先调用多个 solver 尝试并 fallback | 未对齐 | fixed/boxed flip 后在同一 controller 决定 Phase I/II |
| Phase-I bounds | bounds 形状基本存在，但生命周期独立 | 部分 | 由统一 work state 安装/恢复，禁止重建另一套求解路径 |
| rebuild | 分散在周期 reprice/refactor/fallback | 未对齐 | reason-driven 固定 rebuild 顺序和 fresh 标志 |
| CHUZR/DSE | exact row 校验、DSE/Devex recurrence 已有 | 部分 | framework 生命周期、scaled norm、拒绝/重选和 reinvert 保留语义 |
| PRICE | BTRAN/tableau row 已有 | 部分 | pack 顺序、scale、small-pivot refinement/quad retry |
| CHUZC0-4 | CHUZC2/3 和 gated CHUZC4 已有 | 部分 | 统一 work state 输入；free move、完整 retry、启用 CHUZC4 |
| BFRT/FTRAN | Zig 先直接 apply flips，再计算 entering FTRAN | 未对齐 | 先形成 BFRT RHS/FTRAN，再 pivot FTRAN，再统一更新 primal |
| dual update | Zig pivot 后 recurrence + 周期 exact reprice | 未对齐 | `theta_dual`、`shiftCost/shiftBack`、workDual 增量顺序 |
| primal update | basic value 更新与 BFRT 分散 | 部分 | BFRT correction、pivot column update、infeasibility list 同序 |
| pivot verify | 有 residual guard，但不是 row/column alpha 对照 | 未对齐 | `alpha_row/alpha_col` 数值 trouble → rebuild reason |
| factor/basis update | eta update 和 basis commit 已有 | 部分 | 必须作为 iteration 最后一步，失败不提交半状态 |
| Phase-I assessment | 恢复 original 后重新分类并可能切 primal | 未对齐 | objective/cleanup/unperturbed assessment/exit reset duals |
| Phase-II cleanup | recursion、fallback 与证书分支 | 未对齐 | remove perturbation → dual count → optimal/primal cleanup/phase re-entry |
| status/certificate | correctness gate 已有 | 部分 | 只允许 fresh rebuild 后形成最终结论，匹配 HiGHS phase exit |

## 3. 实现顺序和验收门槛

每项只能在前置项完成后启用；实验失败必须回退行为开关，但源码差异仍保留为未完成。

### A. 统一 controller 和 work state

- [x] A1a 定义 HiGHS 同构的 `DualSolvePhase`、`DualRebuildReason`、fresh rebuild 和
  perturb/shift/cleanup flags；接入 root reset、Phase-I/II entry 与 refactor fresh lifecycle。
  ReleaseFast PASS；四模型路径保持 `114/293/1738/1079`，correctness 4/4 PASS。
- [x] A1b 删除 original-cost cleanup 中递归 `solveDual` 的第二次 solver 入口；在同一
  controller/iteration 生命周期切换 Phase II，并用 `dual_control.costs_shifted` 取代独立的
  `shifted_dual_accounting_active`。ReleaseFast 与 highs-strategy 40-model correctness
  **40/40 PASS**，四模型迭代数不变。
- [x] A2a 将 engine 字段重命名为持久 `dual_work`，补齐 `workShift/workRange`，并用
  `DualWorkView` 将 BasisState 的 live lower/upper/value/dual 作为唯一数组暴露；exact work
  reprice 已按 `workCost + workShift` 计算；新增数组计入 `requestedBytes`。ReleaseFast PASS，
  未启用半迁移行为时四模型 correctness 4/4 PASS、迭代数不变。
- [x] A2b 枚举并迁移 basis status/bound/value mutation，使 pivot 和 bound flip 在 Phase I/II
  都持续同步 `nonbasicMove`；move 初始化改为由 work bounds 判定 fixed/free，不再信任可能陈旧的
  status 标签。实现源码同构的 `correctDualInfeasibilities` 后，Phase-II CHUZC 已正式切换 explicit
  move。此前 `bore3d` 24 iterations 错判 infeasible 的根因是 `lower==upper` 变量从 at-lower/
  at-upper 标签错误得到非零 move；修复后四模型恢复 `114/293/1738/1079`，highs-strategy
  forced-dual correctness **40/40 PASS**，无模型特例。
- [ ] A2c 全部 iteration kernel 只读取 `DualWorkView`；原始 bounds/costs 只作初始化源。
  CHUZC 的 dual/move/lower/upper/value 以及 leaving target 已迁移；其余 update/rebuild 调用点待迁移。
- [ ] A3 用一个 major phase loop 替换 shifted-dual/Phase-I/snapshot retry 的求解器拼接；
  暂保留旧 controller 作为 correctness fallback，但不计入“对齐路径”。
- 验收：logical basis 上 Phase 选择、work bounds/value/move、cost/dual snapshot 与 HiGHS
  iteration 0 逐字段一致；40-model correctness 40/40。

### B. 初始化和 rebuild

- [x] B1a 显式实现第一次无扰动 cost initialization：Phase-I 入口现在先把原始 scaled structural
  cost、零 logical cost 和零 workShift 写入持久 work arrays，再执行 `computeDual` 和
  unperturbed dual/primal infeasibility、`force_phase2/near_optimal` assessment，不再依赖上层碰巧
  留下的 reduced costs。ReleaseFast PASS；九模型迭代数和状态完全不变。
- [ ] B1 对齐两次 cost initialization、near-optimal、perturbation 和 cost scale。
  精确移植单次 HiGHS 扰动公式的实验已完成但未启用：仅替换为
  `(1+random)*(abs(cost)+1)*base`、logical cost=0 和全变量 boxed-rate 后，`etamacro`
  `1079 → 2073`。第一次 unperturbed cost/computeDual 和 infeasibility assessment 已由 B1a 补齐；
  当前剩余工作是把第二次初始化按 `near_optimal` 原子切换为原始 cost 或精确 perturbation，并使其
  phase selection 与 B2 生命周期消费同一份 assessment。已保留 `force_phase2`、`near_optimal` 状态，
  并对齐 logical row perturbation `(0.5-random)*1e-12`；禁止按模型筛选。
  B1a 完成后再次隔离精确第二初始化：`blend 114→117`、`etamacro 1079→2073`，其余七个目标
  模型不变，故再次撤回。固定 HiGHS 日志确认 `etamacro` 参考路径为 `DuPh1 152 + DuPh2 372 +
  PrPh2 cleanup 8 = 532`；当前安全路径则含 `dual Phase-I 384 + primal Phase-I 398 + primal
  Phase-II 297`，差距主要来自 dual Phase-I fallback，而非单次迭代内核。
  2026-07-23 用 gdb 截取固定 HiGHS 首次 Phase-I rebuild 后发现：cold solve 会连续调用两次
  `initialiseSimplexLpRandomVectors`，且中间不重置 `random_`；原 `highs_random.zig` 只模拟一轮。
  现已修正双轮 RNG 生命周期、第二轮 permutation/random vector 和 correction RNG 起点，并加入固定
  HiGHS 数值测试。精确 perturbation 实验下 `workCost` 与 permutation 已二进制完全一致，workDual
  最大差仅 `6.35e-14`，dual objective 为 HiGHS `-29.562269393289178`、Zig
  `-29.562269393289181`。由于后续 CHUZR sparse list 尚未对齐，精确公式行为开关仍暂时撤回。
- [ ] B2 对齐 `initialiseBound + initialiseNonbasicValueAndMove` 的 Phase I/II 生命周期。
  通用 nonbasic value/move 初始化已按 work bounds 对齐 fixed/free/one-sided/boxed-invalid。完整 phase
  transition 实验确认两个剩余源码差异：HiGHS 进入 Phase I 前不翻 boxed，退出后 move=0 boxed
  固定到 lower；孤立启用使 `bore3d 293→361`，已撤销行为并等待统一 controller/rebuild。
  初始 `force_phase2` 分支实验通过 40/40 correctness，产生 `brandy 1738→894`、`grow7
  993→413`，但也产生 `klein1 105→265`、`scorpion 358→434`、`scsd1 104→115`、
  `shell 670→686`。拆分确认全部来自 force-phase2，说明 Phase II 更新链仍未对齐；分支已撤销，
  禁止按模型筛选。
  gdb 全数组快照进一步确认 `etamacro` 首次 Phase-I rebuild 的 `basicIndex/workLower/workUpper/
  nonbasicMove/baseValue/baseLower/baseUpper/workShift` 已与 HiGHS 二进制一致。46 个 workValue 差异
  全部是 basic 变量在 Zig full-primal 镜像中的值；HiGHS basic 值只存 baseValue，且两边 baseValue
  完全一致，不构成 CHUZR/CHUZC 输入差异。B2 的首次 rebuild 验收已完成，剩余是 phase exit 生命周期。
- [ ] B3 实现 reason-driven rebuild：reinvert、computeDual、correct dual、computePrimal、
  infeasibility list、objective、fresh flag 固定顺序。已完成初始 Phase-II 的
  `computeDual → correctDualInfeasibilities → computePrimal → fresh`，并接入与 HiGHS 连续的 RNG
  状态；已增加 rebuild primal-infeasibility dense state、汇总和 work dual objective。periodic reinvert
  已统一 exact reprice，但 full computePrimal 暂未启用：完整实验使 `blend 114→111`、`bore3d
  293→280`、`brandy 1738→1710`，同时使 `etamacro 1079→1777`；拆分确认奇异值来自
  computePrimal，而非 correctDual。须先完成 B2 的 nonbasic value 生命周期和 rebuild reason 对齐，
  再统一启用，不能保留三优一退的行为。
  追加 no-entering reason-driven fresh rebuild 实验：`etamacro` 首次失败在 iteration 155，仍有 10 个
  符号合格候选，但 CHUZC 返回 non-positive expanded step/no entering。按 HiGHS 请求 fresh rebuild
  后没有结束 Phase-I，反而使 `etamacro 1079→2367`（dual Phase-I 累计 1674），其余八模型不变；
  已撤回。证明 B3 不能在陈旧的 B2 work value/move 状态上孤立启用。
- 验收：`blend/brandy/bore3d/etamacro` 首次 rebuild 全量 snapshot 一致。

### C. CHUZR 和 edge weights

- [ ] C1 对齐 DSE 初始化、scaled-space norm、1/4 低估拒绝和候选重选。数值 recurrence 和低估重选
  已完成。`HEkkDualRHS` 的 serial-corpus sparse infeasibility `workIndex/workMark/workCount`、平方
  infeasibility、20% dense 切换、动态追加、连续 RNG 随机起点已实现；FTRAN 输出索引按 HiGHS
  `ftranU` 的 reverse-U-pivot 顺序发布，而不是对 dense 结果做自然行号扫描。尚缺 `workCutoff`
  的 >500-candidate hyper-sparse nth-element 路径，故 C1 暂不整项勾选。
- [ ] C2 对齐 DSE recurrence、Devex framework reset、phase/reinvert 生命周期。
- 验收：四模型前 20 次 leaving row、weight、computed weight 一致。

  `etamacro` 精确初始状态 + CHUZC4 实验的首个 pivot 已与 HiGHS 一致（entering 44、leaving row
  284/variable 972）。首次分叉在 iteration 2：row 190 与 row 270 都是 `value=±4、weight=1、
  squared merit=16`；HiGHS 的 64 项 sparse list 从随机位置 9 扫描，动态顺序使 row 270 先被访问，
  Zig dense row-order 扫描选择 row 190。禁止用任意 tie-break 补丁替代，下一步完整移植 RHS list。

  2026-07-23 C1 serial 路径验收：精确 B1 cost + C1 list + D3 CHUZC4 下，Zig 与固定 HiGHS 的
  `etamacro` 前 **60/60** 次 `(entering, leaving variable, step)` 逐项一致。iteration 13
  曾是第二个分叉点：两边 `workCount=105/randomStart=56` 且集合相同，但 Zig 正向 pivot-column
  收集使 row 223 先于 row 222；改为与 `HFactor::ftranU` 一致的 reverse-U 顺序后，第 13--24 次也
  全部对齐；随后又补齐 `HVector::reIndex` 超过 10% density 时恢复自然行序的规则。该组合完整运行时
  `blend=109`，已等于 HiGHS；但 `grow7` 仍被后续未对齐的 B/E/F 路径
  错判 infeasible，其他模型也在 25 次之后分叉。因此 exact B1/C1/D3 的联合行为开关仍关闭，安全
  基线保持不变；下一首分叉已移动到单次迭代更新/phase controller。
  默认安全开关恢复后，ReleaseFast 全测通过；按既定公平口径
  `phase_one=dual + dual_initialization=highs + dual_edge_weight=steepest-devex` 的端到端正确性门禁
  **40/40 PASS**。若把 edge-weight 留为 `inherit`，`agg/grow7` 会走另一条已知未对齐路径并失败，
  因而该配置不能混入本轮与 HiGHS DSE 的迭代数比较。

### D. PRICE 和 CHUZC

- [x] D1a 对齐 CHUZC0 free move：显式 `nonbasicMove` 存在时，`move=0` 只允许真正
  free/superbasic 且双侧无限的列按 `alpha*move_out` 获得本次 ratio test 的临时方向；fixed 和
  bound-active 零 move 不再由 status 隐式复活，且不修改持久 move。ReleaseFast PASS；九模型
  `114/293/1738/1079/413/105/358/104/670` 与安全基线完全一致。诊断 trace 证实 `bore3d`
  没有 bound-active missing move。CHUZC1 pack 顺序与 row scale 仍待逐项 trace 对齐。
- [ ] D2 对齐 CHUZC2 dynamic Ta/Td 和 CHUZC3 large-step 的输入/输出。
- [ ] D3 启用已实现的 CHUZC4 grouping/large-alpha/strict flip set/permutation tie。
- [ ] D4 实现 small-pivot iterative refinement、quad PRICE、删除 pivot后重跑。
- 验收：四模型前 20 次 entering column、alpha、theta、flip set 一致。

  CHUZC4 隔离记录（2026-07-22）：完整开启 grouping 后仍出现 `bore3d 293→521`、
  `brandy 1738→1078`、`scorpion 358→1418`、`scsd1 104→1431`、`shell 670→1313`，且
  `grow7` 被错误判 infeasible，因此行为开关保持关闭。源码核对确认 CHUZC4 选中列满足
  `dual*move<=0` 时 HiGHS 设置 `theta_dual=0`，随后走 `shiftCost(variable_in,-workDual)`；Zig
  当前统一使用 `r_q/alpha` 更新。曾把所有路径直接改成 ratio-test theta，导致多个模型在几十次
  迭代内误判 infeasible；撤回后安全基线完全恢复。结论：D3 不能孤立启用，必须先完成 E3 的
  `theta_dual==0` shiftCost/shiftBack 和 workShift 生命周期，再重新验收 CHUZC4。

  后续已将 `zero_dual_step` 显式接入现有 `shiftCost/shiftBack` 存储，并保证普通路径仍使用已验证的
  tableau recurrence；重新开启 CHUZC4 后上述九模型结果逐项完全不变，说明这些模型在首次分叉前
  没有走零 theta 分支。进一步源码核对得到当前更早的结构差异：HiGHS CHUZC1 按动态 row-wise
  matrix 的 `row_ap.index` 顺序、再按 `row_ep.index + num_col` 顺序形成 pack；Zig 当前按全局列号
  dense 扫描。CHUZC3/4 的原位 swap/reduce 对 pack 顺序敏感，因此 D3 的下一前置项是对齐
  CHUZC1 pack 顺序，而不是继续修改 grouping 公式。CHUZC4 行为开关再次关闭。

  追加隔离（2026-07-22）：曾把 CHUZC4 theta 改成 `dual/(tableau*move_out*move)`，会使
  `blend/bore3d/brandy/grow7/scorpion/shell` 很快误判 infeasible。复核源码确认 CHUZC2 写入
  `workData.second` 时已经乘过一次 `move_out*move`，CHUZC4 形成 `workAlpha` 时再乘一次，符号平方
  后恢复原始 tableau alpha；因此 `dual/tableau` 原实现正确，effective-alpha theta 方案已撤回。
  随后的静态 CSR pack-order 实验在 CHUZC4 开启时仅改变 `blend 114→132`、`etamacro
  1079→1229`，其余退化和 `grow7` 错判完全不变；关闭 CHUZC4 后八模型全部恢复原基线。HiGHS
  的真实 pack 还依赖 HVector 非零顺序及动态 partition swap，静态 CSR 不是源码等价实现，相关
  实验代码已撤回。当前结论：在 B1/B2 初始 work state 尚未对齐前，不能用 CHUZC4 轨迹判定
  kernel 差异，执行顺序回到 B1→B2→B3，再重新开启 D3。

  2026-07-23 在双轮 RNG、精确 cost 和首次 rebuild 全数组对齐后重开 D3：`etamacro` 首个 CHUZC4
  选择已正确从 236 变为 HiGHS 的 44，前 25 个 pivot 的变量集合高度一致；但由于 C1 sparse
  CHUZR 顺序尚缺，完整路径仍产生 `grow7/shell` 错判，故 CHUZC4 开关再次撤回，安全路径不变。

### E. 单次迭代更新顺序

- [ ] E1 BFRT collect/FTRAN 与 pivot-column FTRAN 同构。
- [ ] E2 row/column pivot verify 产生 rebuild reason，不在半更新状态直接退出。
- [ ] E3 对齐 incremental dual shift/update、primal/BFRT update、DSE/Devex update。当前首要项：
  `theta_dual==0 → shiftCost` 和 `shiftBack(out)` 已接入持久 `workShift`；剩余源码差异是非零 theta
  packed-row update、`workDual[in]=0/workDual[out]=-theta` 的严格更新顺序。禁止把 CHUZC theta
  直接塞入现有 tableau recurrence。
  2026-07-23 已删除 dual Phase-I 的“每个 pivot 后 full reprice”旁路：Phase-I 现在与 Phase-II
  共用增量 dual/shiftBack 更新，只在真实 reinversion 后 exact reprice，符合
  `HEkkDual::updateDual` 生命周期。ReleaseFast 全测及固定配置 40-model correctness **40/40 PASS**；
  八模型 committed iterations 均不变，因此该项只记为结构和单次迭代开销修复，不记作迭代收益。
  当前 40 模型横向结果为 Zig 胜 7、平 3、负 30，总迭代 `13509 vs 8774`（1.540x）；
  最大差距 `brandy 1735 vs 304`。E3 尚未完成，下一步是消除 dense tableau recurrence 与
  HiGHS packed-row `updateDual(theta)` 的数据流差异。
- [ ] E4 basis/factor/matrix/infeasibility-list 最后原子提交。
- 验收：四模型前 20 次 pivot 后 workDual/baseValue/weights/basis head 一致。

### F. Phase exit、cleanup 和结论

- [ ] F1 对齐 Phase-I objective assessment、remove perturbation 和 Phase-I re-entry。
- [ ] F2 对齐 `exitPhase1ResetDuals`、free dual shift、Phase-II transition。
- [ ] F3 对齐 Phase-II optimal/possibly-unbounded fresh rebuild、cleanup levels和最终状态。
- 验收：禁止 snapshot retry 和非 HiGHS primal fallback 后，40-model correctness 40/40。

### G. 迭代数验收

- [ ] G1 四模型完整 pivot trace 与 HiGHS 对齐。
- [ ] G2 完整 40-model 逐模型记录 Zig/HiGHS committed iterations 和奇异值。
- [ ] G3 所有未对齐模型回到首次分叉所属模块，不允许跨模块补丁。

当前执行点：**E1 BFRT/FTRAN → E3 incremental update → F1/F2 phase exit**。只接受两类成果：
40-model correctness 不退化；相同配置下完整 committed iterations 降低并最终不高于 HiGHS。
局部 tick、前缀轨迹和源码相似度只用于定位，不计入验收结果。

#### 已拒绝并删除的实验（2026-07-23）

曾为复现 HiGHS synthetic-clock reinversion 增加独立 U/UR shadow graph。它能让 `etamacro` 与
HiGHS 同在 update 54 触发 rebuild，但累计 tick 仍差 30，未改善任何完整模型迭代数，并使
ReleaseFast 热路径约退化 10%。该实验层、字段、开关和额外接口参数已全部删除；仅保留结论：
reinversion 时点是后续差异，不是当前最大收益点。若未来重做，必须直接复用真实 factor graph，
不得再引入平行 shadow 状态。

清理后 ReleaseFast 全测、`git diff --check` 均通过。固定公平配置的八模型安全基线为：
`blend 114 / bore3d 293 / brandy 1735 / etamacro 1095 / grow7 413 / scorpion 358 /
scsd1 104 / shell 670`，状态全部 optimal。后续改动必须按 E1/E3/F 的模块边界完成，禁止
模型特例、固定 iteration 阈值或跨模块补偿。

## 4. 可维护性注释补全进度

本节只追踪文档工作，不改变上述算法对齐状态。注释必须说明数据所有权、缩放坐标、
状态生命周期、数值失败条件和热路径约束，不能用重复函数名的空泛描述代替。

- [x] simplex 公共模型与基础状态：`problem.zig`、`solution.zig`、`basis.zig`、
  `basis_snapshot.zig`、`numerical.zig`、`root.zig`。
- [x] simplex 定价与比率测试：`pricing.zig`、`pricing_workspace.zig`、
  `ratio_test.zig`。
- [x] dual 工作状态与算法内核：`dual_state.zig`、`dual_phase_one.zig`、
  `highs_random.zig`、`engine_dual.zig`、`engine_dual_candidate.zig`、
  `engine_dual_edge_weight.zig`。
- [x] primal、pivot、basis 与 solve controller：`engine_primal.zig`、
  `engine_pivot.zig`、`engine_basis.zig`、`engine_setup.zig`、
  `engine_finish.zig`、`engine_progress.zig`、`engine_degeneracy.zig`。
- [x] factorization、crash 与退化工作区：`factorization.zig`、`crash.zig`、
  `degeneracy.zig`。
- [x] 2026-07-23 simplex 注释批次验证：`zig fmt`、ReleaseFast 全量测试和
  `git diff --check` 全部通过；未修改求解行为。
- [x] matrix 注释补全。
  - [x] 核心表示与乘法/缓存首批：`csc.zig`、`csr_view.zig`、`dense_lu.zig`；
    已补齐所有权、packed/page-colored 布局、字段语义、checked/assume-valid 边界和
    核心乘法/求解路径说明。
  - [x] 首批验证：`zig fmt`、ReleaseFast 全量测试和 `git diff --check` 全部通过。
  - [x] builder/vector/transform/edit/store。
  - [x] sparse symbolic/kernel/LU/Forrest--Tomlin：已说明所有 SoA 平行数组、
    intrusive list/free-list、有效前缀、置换坐标、symbolic/numerical 分界及更新链生命周期。
- [x] io 注释补全：通用输入输出与资源限制、语义 builder、packed `ModelData`、
  string arena、LP lexer/parser/writer、MPS parser/writer 均已覆盖。
- [x] 最终机械审计：`simplex`、`matrix`、`io` 中每个函数声明前均存在作用说明；
  三个模块的主要结构体字段均说明所有权、有效范围或状态含义。
- [x] 最终验证：全模块 `zig fmt`、ReleaseFast 全量测试及 `git diff --check` 通过。
