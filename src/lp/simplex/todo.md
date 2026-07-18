# Simplex roadmap

完成的工作使用 `[x]` 标记；尚未完成或仍需端到端验收的工作使用 `[ ]`。

## 1. 完整 simplex iteration

- [x] Pricing，包括 bound-aware primal pricing 和 dual leaving-row pricing。
- [x] Entering/leaving variable 选择。
- [x] Harris primal ratio test 和 dual ratio test。
- [x] Boxed-variable bound flip。
- [x] Basis replacement 及 factorization update。
- [x] Primal/dual infeasibility 增量更新与重新分类。

## 2. Factorization 生命周期

- [x] FTRAN `aq` 与 update capture。
- [x] BTRAN `ep`。
- [x] Eta/Forrest--Tomlin update。
- [x] Growth、update limit 和 residual 触发 reinversion。
- [x] Reinversion 后恢复 pricing、Devex/steepest-edge 状态。

## 3. 数值鲁棒性

- [x] Rank-deficient imported basis repair：确定性地以非基 logical columns
  替换 structural basics，并保留能够继续使用的 warm-basis 前缀。
- [x] Allocation-free iterative refinement。
- [x] Threshold/small-pivot rejection。
- [x] 最优解发布前的 primal/dual residual 与人工变量校验。
- [x] Stalled/degenerate pivot 监控；连续零步长触发 Bland-style
  entering/leaving 与字典序扰动 tie-break，有效移动后退出 fallback。
- [x] 二次幂 row/objective scaling，并仅在 row-scaled matrix 动态范围超过
  `1e6` 时启用 column scaling；solution view 还原 primal/dual/reduced cost。
- [x] 原始坐标下统一 `1e-9` model-coefficient dropping policy，basis assembly、
  pricing、FTRAN residual 和 certificate validation 使用同一结构。
- [x] Unbounded primal ray 构造、原坐标发布及 bounds/row/objective 独立验收。

## 4. Phase I / Phase II 与状态输出

- [x] Artificial variables 和 Phase-I objective。
- [x] Infeasible/unbounded 判定。
- [x] Phase-I cleanup 与 Phase-II objective transition。
- [x] Phase-I reduced cost 使用一次 CSC 扫描，消除逐列 dense clear/dot；
  primal pivot 后增量维护并每 8 次精确稀疏刷新，病态模式每次刷新。
- [ ] 实现带 cost perturbation 的专用 dual Phase I，再允许 neither-feasible
  cold logical crash 走 dual -> primal；现有 warm-basis repair 不满足该前提。
  具体实现和验收顺序见第 6.3 节，完成对应阶段前不得提前启用自动算法切换。
- [x] Basis status、primal、dual、reduced cost 和模型 solution 回写。

## 5. 阶段零：冻结最小端到端验收基线

本节是进入第 6 节性能优化前必须完成的前置关卡。这里不要求先跑完 Netlib 和
Mittelmann，而是冻结可重复的快速回归集合、当前已知失败和基础结果，确保后续
每项优化都能判断收益来源并及时发现正确性回归。

- [x] 建立 MPS -> canonical CSC -> simplex 的 zhighs/HiGHS serial-simplex
  差分 runner，统一输出 objective、status、iterations、residual 和耗时。
- [x] 修复 `sc105` 暴露的 FT `ep` 历史 correction 重复应用；增加
  allocation-free FTRAN backward-error gate 和 caller-owned pivot trace。
- [x] 修复 `brandy` 在相同 FT/reinvert trace 下复现的 Phase-I 零步长循环：
  artificial/bounded Phase I 保留 Devex/Harris，目标切换时清除 Bland fallback。
- [x] 修复 fixed-MPS 中 `RHS`/`RANGES` 同名 set 与省略 set field 导致的
  约束静默丢失；本地 40 模型扩展 corpus 已有 35 个匹配 HiGHS。
- [x] 修复 `grow7`、`blend`、`scsd1` 和 `vtp-base` 的非 FT 数值失败：
  FTRAN 使用 `|a_q| + |B||x|` 后向误差尺度；小 pivot 在换基前以 fresh
  factorization 重算，并在强病态 basis 上切换到本次 solve 的安全模式。
- [x] 冻结现有 40 模型快速 corpus、输入文件校验值、构建模式、随机种子和
  HiGHS 对照版本；逐例标记当前通过、失败或超时状态，不要求先修复全部已知失败。
- [x] 为快速 corpus 保存当前 objective、status、iteration count、residual 和
  完整 solve 时间，作为第 6 节所有阶段共同使用的 pre-optimization baseline。
- [x] 固定 `gas11` 的 unbounded ray 验收，以及 `brandy`、`sc105`、`scsd1` 的
  basis/pivot trace；连续运行结果必须可重复。
- [x] 确认单元测试和快速 corpus runner 能以单一命令执行，并在任一非预期
  status、objective、residual 或 certificate 变化时返回失败。

## 6. 下一阶段性能优化（严格执行顺序）

只有第 5 节全部勾选后才能进入本节。以下 6.1 -> 6.2 -> 6.3 -> 6.4 -> 6.5
按顺序实施，每一阶段必须完成正确性测试、固定 corpus 差分和指标记录，再勾选并
进入下一阶段；
第 6 节完成后进入第 7 节完整验收。借鉴 HiGHS 的状态机和数值策略，但保持 zhighs
的 SoA、可复用 workspace、allocation-free hot path 和零拷贝稀疏接口。启发式仅
用于选择搜索路径；Phase 切换、optimal、infeasible 和 unbounded 结论必须由
无扰动状态下的 fresh rebuild/reprice 验证。

### 6.1 阶段一：建立可归因的性能基线

- [x] 增加 Phase I、Phase II 和 cleanup 的独立耗时与 iteration 计数。
- [x] 分别统计 INVERT、FTRAN、BTRAN、PRICE、UPDATE 和 rebuild 的耗时与调用次数。
- [x] 记录退化/零步长 pivot、bound flip、reinversion 原因、FT chain 长度、
  factor growth 和 solve residual。
- [x] 记录 RHS、`aq`、`ep` 和 pricing 向量密度，以及 dense/hyper、row/column
  kernel 的选择次数；kernel 仅在一次粗粒度操作开始时分发，禁止在非零元素
  内循环中动态判断。
- [x] 为 `gas11`、`brandy`、`sc105`、`scsd1` 建立固定输入、固定随机种子和
  pivot/basis trace；输出 warmup 后多轮 median、p95、requested bytes 和 peak RSS。
- [x] 同一 corpus 上记录 zhighs 与 HiGHS 的 status、objective、iterations、
  residual、阶段耗时和总耗时，形成后续阶段共同使用的基线报告。

基线确定的优化优先级：

1. **P0 `scsd1` rebuild 归因与消除**：当前 347 次 rebuild 中只有 3 次能够由
   update-limit/growth、solve-residual 或 small-pivot policy 解释，必须先闭合原因计数。
2. **P1 `brandy` Phase I**：2987 次 Phase-I iteration、3101 次退化 pivot，
   Phase I 和 pricing 是其相对 HiGHS 慢约 2.65 倍的主要来源。
3. **P2 共享 pricing/kernel 优化**：`brandy` 和 `scsd1` 的 profiling PRICE
   分别约 16.0 ms 和 4.7 ms，但须在 pivot/rebuild 路径稳定后再改变 representation。

### 6.2 阶段二：`scsd1` rebuild 原因闭环与更新路径修复（P0）

- [ ] 扩展 rebuild reason，至少区分 policy reinversion、FT update rejection、
  `direction_requires_reinversion`、`fresh_factorization_mode`、cleanup 和 edge-weight
  reset；任何 rebuild 都必须且只能归入一个原因。
- [ ] 对 dense Eta 和 sparse FT update 的失败错误分类计数，保留首个失败位置和
  固定 pivot-trace suffix；诊断热路径不得格式化字符串或分配内存。
- [ ] 定位 `scsd1` 约 344 次未归因 rebuild 的主导来源，修复 update capture、
  pivot acceptance 或 safe-mode 生命周期；禁止仅通过放宽 residual/pivot 容差掩盖错误。
- [ ] 避免已知不可更新的 basis 在每轮重复尝试同一 FT 路径；fallback 必须具有
  basis-epoch 生命周期，并在 fresh rebuild 验证后才能恢复更新。
- [ ] A/B 验收 `scsd1` 的 rebuild 次数、INVERT/rebuild 时间和完整 solve 时间；
  同时保持 objective、residual 和固定 trace 可解释，并保证 40 模型 corpus 不回退。

### 6.3 阶段三：专用 dual Phase I 与算法选择（P1）

- [ ] 增加持久复用的 `DualPhaseOneWorkspace`，以 SoA 保存工作 bounds、cost、
  dual infeasibility 和 perturbation；迭代热路径不得分配或复制模型矩阵。
- [ ] 根据变量 bound 类型构造 dual Phase-I 工作边界和 infeasibility objective，
  不再用零目标 primal warm-basis repair 代替 dual Phase I。
- [ ] 实现与 basis epoch 绑定的确定性 cost perturbation，统一 perturbation 的
  启用、移除和恢复流程。
- [ ] dual Phase I 退出时移除 perturbation，执行 fresh rebuild、BTRAN/reprice
  和 primal/dual infeasibility 重算；只有重验通过才允许进入 Phase II。
- [ ] 根据 logical/crash basis 的 primal/dual infeasibility 数量与幅度、free/boxed
  比例以及 coefficient/cost range，确定性选择 primal 或 dual 路径。
- [ ] 先以显式选项启用并完成 `gas11`、`brandy`、`sc105`、`scsd1` 回归；确认
  status/objective/ray 或 certificate 后，再启用自动 dual -> primal 选择。
- [ ] 以 `brandy` 为主要性能验收：分别记录 primal/dual Phase-I iterations、
  退化 pivot、PRICE 时间和总耗时；只有显著优于冻结的 primal Phase-I baseline
  且其他 corpus 无回退时，dual Phase I 才能参与默认算法选择。

### 6.4 阶段四：稀疏 crash basis（P1）

- [ ] 实现 LTSSF 风格的 singleton/near-singleton peeling，优先形成结构三角、
  低 fill 的初始 basis。
- [ ] 使用 column nnz、可接受 pivot 幅度、bound/cost compatibility 的确定性评分，
  并复用 crash workspace，禁止候选扫描期间逐列分配。
- [ ] 对 crash basis 执行 factorization、rank、growth 和 residual 验证；失败时
  确定性回退到 logical basis，不得发布部分损坏的 basis。
- [ ] 在 LTSSF 基线稳定后，再评估 Bixby 风格评分；必须通过 A/B 数据证明其
  Phase-I iterations 或总耗时收益后才能成为默认策略。
- [ ] 比较 crash 前后初始 fill、INVERT 时间、Phase-I iterations、退化 pivot
  数量和完整 solve 时间，重点验收 `brandy`。
- [ ] 若 crash 没有降低 `brandy` 的 Phase-I iterations 或总耗时，则保留为显式
  策略而非默认路径，并记录失败原因，禁止仅凭初始 fill 下降判定成功。

### 6.5 阶段五：退化 pivot 与 pricing 路径优化（P2）

- [ ] 将现有全局 Bland fallback 升级为 basis-epoch 级确定性 bound/cost
  perturbation，并以连续零步长和目标改善停滞作为触发条件。
- [ ] 记录坏的 `(entering, leaving, bound_direction)` basis change，增加短期
  taboo、suffix repair 或 backtracking；有效移动或 reinversion 后按策略失效。
- [ ] perturbation cleanup 后强制 fresh rebuild/reprice，并在原始无扰动坐标下
  重验 primal、dual、reduced cost 和最终状态。
- [ ] 根据实测 `ep`/pricing 密度，在一次 pricing 操作边界选择 row pricing 或
  column pricing；保留可强制指定 kernel 的选项用于固定 trace A/B 对照。
- [ ] 为 column pricing 增加可复用的 row-wise/partitioned matrix view，只有当
  `brandy`/`scsd1` 的 PRICE median 和总耗时同时受益时才启用；view 构建成本、
  requested bytes 和 cache 行为必须计入完整 solve。
- [ ] 以 `brandy` 为主验收退化 pivot 和总 iterations 的下降，以 `scsd1` 验收
  rebuild 修复后的 Phase-II pricing 收益，同时保证其他
  corpus 的 objective、status、residual、ray/certificate 和性能不回退。

### 6.6 每阶段提交与验收规则

- [x] 阶段一提交：`bench: add phase-level simplex counters`。
- [ ] 阶段二建议提交：`simplex: classify and reduce rebuild fallbacks`。
- [ ] 阶段三建议提交：`simplex: implement dedicated dual phase one`。
- [ ] 阶段四建议提交：`simplex: add sparse triangular crash basis`。
- [ ] 阶段五建议提交：`simplex: add epoch perturbation and adaptive pricing`。
- [ ] 每次提交前运行单元测试、固定 corpus 差分和对应 benchmark；将结果写入
  可由 git 追踪的报告，报告同时保留测试环境、编译模式和 HiGHS 版本。
- [ ] 未达到本阶段正确性门槛时不得以性能提升为理由进入下一阶段；未经 fresh
  rebuild 验证的启发式结果不得用于最终状态或 certificate。

## 7. 优化完成后的完整端到端验收

本节只在第 6 节全部完成后执行，避免 dual Phase I、crash basis 和退化策略改变
pivot 路径后重复生成大规模报告。快速 corpus 仍须在第 6 节每次提交前运行。

- [ ] Netlib 完整求解结果与 HiGHS/CLP 对比。
- [ ] Mittelmann 完整求解结果与 HiGHS/CLP 对比，并单独记录超时和内存上限。
- [ ] 汇总 objective、status、iteration count、primal/dual residual、ray 和
  infeasibility certificate。
- [ ] 汇总 Phase I/II iterations、退化 pivot、reinversion 次数与原因、FT chain
  长度、factor growth 和 solve residual。
- [ ] 汇总 INVERT、FTRAN、BTRAN、PRICE、UPDATE、rebuild 和完整 solve 时间，
  同时报告 warmup 后多轮 median、p95、requested bytes 和 peak RSS。
- [ ] 对第 5 节 baseline、各阶段提交和最终结果生成可由 git 追踪的对比报告，
  明确性能改善、回退、已知失败和下一轮优化候选。
