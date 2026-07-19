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
- [x] 实现带 cost perturbation 的专用 dual Phase I，并以显式策略验收
  neither-feasible cold logical/crash basis 的 dual -> primal 路径；现有实现已通过
  40 模型正确性 gate，但 `brandy` 未达到性能门槛，因此不参与当前默认算法选择。
  具体实现、拒绝依据和后续 scale-aware 重设计见第 6.3 节。
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

- [x] 扩展 rebuild reason，至少区分 policy reinversion、FT update rejection、
  `direction_requires_reinversion`、fresh-factorization epoch、cleanup 和 edge-weight
  reset；任何 rebuild 都必须且只能归入一个原因。
- [x] 对 dense Eta 和 sparse FT update 的失败错误分类计数，保留首个失败位置和
  固定 pivot-trace suffix；诊断热路径不得格式化字符串或分配内存。
- [x] 定位 `scsd1` 约 344 次未归因 rebuild 的主导来源，修复 update capture、
  pivot acceptance 或 safe-mode 生命周期；禁止仅通过放宽 residual/pivot 容差掩盖错误。
- [x] 避免已知不可更新的 basis 在每轮重复尝试同一 FT 路径；fallback 必须具有
  basis-epoch 生命周期，并在 fresh rebuild 验证后才能恢复更新。
- [x] A/B 验收 `scsd1` 的 rebuild 次数、INVERT/rebuild 时间和完整 solve 时间；
  同时保持 objective、residual 和固定 trace 可解释，并保证 40 模型 corpus 不回退。

### 6.3 阶段三：专用 dual Phase I 与算法选择（P1）

- [x] 增加持久复用的 `DualPhaseOneWorkspace`，以 SoA 保存工作 bounds、cost、
  dual infeasibility 和 perturbation；迭代热路径不得分配或复制模型矩阵。
- [x] 根据变量 bound 类型构造 dual Phase-I 工作边界和 infeasibility objective，
  不再用零目标 primal warm-basis repair 代替 dual Phase I。
- [x] 实现与 basis epoch 绑定的确定性 cost perturbation，统一 perturbation 的
  启用、移除和恢复流程。
- [x] dual Phase I 退出时移除 perturbation，执行 fresh rebuild、BTRAN/reprice
  和 primal/dual infeasibility 重算；只有重验通过才允许进入 Phase II。
- [x] 根据 logical/crash basis 的 primal/dual infeasibility 数量与幅度、free/boxed
  比例以及 coefficient/cost range，确定性选择 primal 或 dual 路径。
- [x] 先以显式选项启用并完成 `gas11`、`brandy`、`sc105`、`scsd1` 回归；确认
  status/objective/ray 或 certificate 后，再启用自动 dual -> primal 选择。
- [x] 以 `brandy` 为主要性能验收并记录 primal/dual Phase-I iterations、退化 pivot、
  PRICE 时间和总耗时；验收结论为“未通过默认启用门槛”，因此 dual Phase I 不参与
  当前默认算法选择。此勾选表示 A/B 和决策闭环完成，不表示性能门槛已经通过。

6.3 当前验收记录（ReleaseFast，7 轮 median）：

- `scsd1` 自动选择 dual Phase I：完整 solve `9.55 ms -> 6.21 ms`，总 iterations
  `587 -> 179`，PRICE `3.82 ms -> 1.73 ms`；dual Phase-I 179 次，Phase II 0 次，
  fresh cleanup 后 objective/primal/dual residual 通过验收。
- `gas11` 自动保留 primal Phase I，无界 ray 验收通过；`sc105` 初始 primal
  feasible，不进入 Phase I；自动策略的 40 模型 corpus status/objective/residual/
  ray gate 全部通过。
- `brandy` 的显式 dual Phase-I 在首轮 ratio test 无可用 entering column，事务性
  恢复后可完整复现冻结 primal 路径，但不构成 dual Phase-I 性能收益。自动策略
  因此仍选择 primal，`SolveControl` 默认也继续保持 `.primal`。第 6.4 节 LTSSF/Bixby
  crash basis A/B 同样未解除该限制，当前阶段以“实现完成、默认启用被拒绝”闭环。

#### 6.3 后续独立优化：scale-aware dual Phase I（不阻塞第 7 阶段）

以下工作属于下一轮算法研究，不参与当前默认选择，也不阻塞第 7 阶段对已验证默认路径
的完整验收。禁止使用 `brandy` 名称特判或通过放宽 feasibility/pivot tolerance 绕过问题。

- [ ] 重新设计 scale-aware working bounds，使基本变量 violation、boxed bound-flip
  capacity 和 row/column scaling 处于一致坐标；替代当前固定单位 `[0,1]` 工作边界。
- [ ] 为 boxed/free 变量建立明确的 nonbasic move 与剩余 violation 表示；bound flips
  无法消除 leaving-row violation 时，必须继续产生数值有效的 pivot 候选或给出可验证拒绝。
- [ ] 在原坐标下验证 Phase-I objective 单调性、cleanup 后 primal/dual feasibility 和
  certificate；失败必须恢复最近 verified logical/crash basis epoch。
- [ ] 以 `brandy` 和第 7 阶段暴露的同类模型做显式 A/B；只有完整 corpus 无正确性和
  性能回退，且显著优于当前 primal Phase-I/taboo 路径时，才重新评估自动启用。

### 6.4 阶段四：稀疏 crash basis（P1）

- [x] 实现 LTSSF 风格的 singleton/near-singleton peeling，优先形成结构三角、
  低 fill 的初始 basis。
- [x] 使用 column nnz、可接受 pivot 幅度、bound/cost compatibility 的确定性评分，
  并复用 crash workspace，禁止候选扫描期间逐列分配。
- [x] 对 crash basis 执行 factorization、rank、growth 和 residual 验证；失败时
  确定性回退到 logical basis，不得发布部分损坏的 basis。
- [x] 在 LTSSF 基线稳定后，再评估 Bixby 风格评分；必须通过 A/B 数据证明其
  Phase-I iterations 或总耗时收益后才能成为默认策略。
- [x] 比较 crash 前后初始 fill、INVERT 时间、Phase-I iterations、退化 pivot
  数量和完整 solve 时间，重点验收 `brandy`。
- [x] 若 crash 没有降低 `brandy` 的 Phase-I iterations 或总耗时，则保留为显式
  策略而非默认路径，并记录失败原因，禁止仅凭初始 fill 下降判定成功。

6.4 当前 LTSSF 验收记录（ReleaseFast，7 轮 median）：

- 新增复用型 `CrashWorkspace`，构建 CSR companion view，并按 active column nnz、
  bound/cost compatibility、scaled pivot 幅度和 row degree 做确定性匹配；候选扫描和
  peeling 热路径不分配。near-singleton 整批 basis 若不可逆，则二分缩短到最长可验证
  前缀，并保留至少 1/4 logical columns 作为数值锚点。
- 每个候选 basis 必须通过 fresh LU、rank、pivot growth/condition 和 basic residual
  验证；后续 dual Phase-I 若失败，也会恢复完整 logical basis epoch。显式 LTSSF 的
  40 模型 status/objective/residual/ray gate 全部通过。
- `brandy` 安装 41 个结构列后仍在 dual Phase-I 首轮 ratio test 回落，最终完整复现
  logical 路径：iterations `3384 -> 3384`、Phase-I `2987 -> 2987`、退化 pivot
  `3101 -> 3101`，完整 solve `62.39 ms -> 63.65 ms`，PRICE `15.94 ms -> 16.07 ms`。
  因此 `.automatic` 和默认策略都继续固定为 logical，LTSSF 仅保留显式 A/B 选项。
- `scsd1` 的 LTSSF 显式路径可将 iterations `587 -> 259`、退化 pivot `542 -> 7`、
  PRICE `3.88 ms -> 2.53 ms`，但 crash/INVERT 开销使完整 solve `9.74 ms -> 9.76 ms`，
  仍不足以支持默认启用。
- Bixby 风格评分复用相同 workspace、peeling 和数值验收，只将候选优先级改为
  scaled pivot quality 相对于 row/column fill、bound penalty 和 cost penalty 的确定性
  merit。ReleaseFast 21 轮 median 中，`scsd1` 的显式 Bixby 路径将 structural columns
  限制为 7 个，iterations `587 -> 173`、退化 pivot `542 -> 7`、PRICE
  `3.77 ms -> 1.66 ms`，完整 solve `9.51 ms -> 6.48 ms`；这是可复现的局部收益。
- 同一验收中 `brandy` 安装 82 个结构列后仍回落到 primal Phase I，iterations 和
  退化 pivot 不变，完整 solve `62.14 ms -> 63.01 ms`。显式 Bixby 的 40 模型
  status/objective/residual/ray gate 全部通过，但跨模型收益不稳定，因此保留为显式
  A/B 策略，`.automatic` 和默认路径继续使用 logical basis。

### 6.5 阶段五：退化 pivot 与 pricing 路径优化（P2）

- [x] 先冻结 `brandy` 的退化归因 trace：按 bound tie、ratio tie、零 primal step、
  零 Phase-I objective improvement、repeated basis signature、small pivot rejection 和
  bound flip 分类；计数必须互斥，热路径只写 caller-owned POD buffer，不格式化或分配。
- [x] 单独追踪 `brandy` 显式 dual Phase-I 首轮“无 entering column”：记录 leaving row、
  `ep`、候选列的 bound/status、signed pivot、reduced cost 和 rejection reason，确认是
  数学上不适用还是候选过滤/符号约定错误；在闭环前不得把该回退固化为模型特判。
- [x] 将现有全局 Bland fallback 升级为 basis-epoch 级确定性 bound/cost
  perturbation，并以连续零步长和 Phase-I objective 改善停滞作为双触发条件；使用
  复用型 SoA workspace 保存原始值、shift 和 generation，pivot 热路径不得分配。
- [x] perturbation 幅度必须由 row/column/objective scale 和可行性容差共同约束，使用
  basis epoch 与变量 id 生成跨平台确定顺序；增加启用迟滞、最大存活 epoch 和幅度
  升级上限，禁止扰动随迭代无界累积或依赖遍历哈希表的非确定顺序。
- [x] 记录坏的 `(entering, leaving, bound_direction)` basis change，增加短期
  taboo、suffix repair 或 backtracking；使用固定容量 generation-marked ring/cache，
  有效目标改善、fresh reinversion 或 basis epoch 切换后按策略失效。
- [x] 增加低成本 basis fingerprint/repetition window，只用于触发诊断和退化策略，
  不作为正确性判据；对命中记录 pivot trace suffix，以区分真正 cycling、长退化链和
  浮点 tie，不允许在每轮复制完整 basis。
- [x] 在 perturbation 下重新验收 Harris two-pass ratio candidate set、boxed bound flip
  和稳定 tie-break；显式提供 baseline/perturbed/taboo 三种强制模式，保证固定 trace
  A/B 能分别归因，禁止多个启发式同时启用后只比较最终总耗时。
- [x] perturbation cleanup 后强制 fresh rebuild/reprice，并在原始无扰动坐标下
  重验 primal、dual、Phase-I artificial objective、reduced cost 和最终状态；cleanup
  失败必须恢复最近一个已验证 basis epoch，不能携带 perturbed reduced cost 发布结果。
- [x] 为 Phase I 独立比较 Dantzig、Devex 和近似 steepest-edge pricing，维护
  generation-marked edge weights 和定期精确刷新；仅当 `brandy` 的 pivot 数及完整 solve
  同时下降才允许替换默认 Devex，禁止用更昂贵 pricing 换取无意义的 iteration 降低。
- [x] 为增量 reduced cost 增加抽样 drift 指标和自适应 exact reprice：按 basis update、
  退化链长度和误差阈值触发，而不是固定每 8 次刷新；刷新频率、PRICE 时间和错误上界
  必须纳入 A/B 报告。
- [x] 根据实测 `ep`/pricing 密度，在一次 pricing 操作边界选择 row pricing 或
  column pricing；保留可强制指定 kernel 的选项用于固定 trace A/B 对照。
- [x] 为 column pricing 增加可复用的 row-wise/partitioned matrix view，只有当
  `brandy`/`scsd1` 的 PRICE median 和总耗时同时受益时才启用；view 构建成本、
  requested bytes 和 cache 行为必须计入完整 solve。
- [x] 以 `brandy` 为主验收退化 pivot 和总 iterations 的下降，以 `scsd1` 验收
  rebuild 修复后的 Phase-II pricing 收益，同时保证其他
  corpus 的 objective、status、residual、ray/certificate 和性能不回退。

6.5 严格按“归因 trace -> epoch perturbation -> taboo/repetition -> cleanup ->
Phase-I pricing -> row/column kernel”的顺序实施。每个策略先以显式控制单独通过 40 模型
corpus，再进行组合；任何策略成为默认前，`brandy` ReleaseFast 21 轮 median 至少满足
退化 pivot 或总 iterations 下降 30%、完整 solve 下降 10%，且 p95、requested bytes、
peak RSS 和其他 corpus 不出现不可解释回退。最终目标不是仅降低 PRICE 单轮耗时，而是
将 `brandy` 的 2987 次 Phase-I iteration 和 3101 次退化 pivot 显著压低，并逐步缩小
与冻结 HiGHS 约 304 次总迭代的算法路径差距。

6.5 当前退化归因记录（ReleaseFast，logical basis，primal Phase I）：

- 新增 caller-owned `DegeneracyTraceEvent` 和互斥 scalar counters。普通求解只执行 O(1)
  分类；ratio tie 只在 statistics/trace 模式下读取已有 ratio workspace，basis fingerprint
  只在显式 trace 模式扫描，生产热路径不分配、不复制 basis，也不改变 ratio-test 决策。
- `brandy` 仍为 objective `1518.509896488128`、3384 iterations、2987 次 Phase-I
  iteration 和 3101 次退化 pivot；分类总数严格闭合为 3101：bound tie `871`、ratio tie
  `1780`、普通零 primal step `232`、Phase-I objective stall `218`，small-pivot retry 和
  bound flip 均为 `0`。
- 32-entry basis fingerprint window 在完整 3101-event trace 中没有 repeated-basis 命中，
  表明当前主导问题是连续 ratio/bound tie 形成的长退化链，而不是短周期 cycling。下一步
  优先诊断 dual Phase-I 首轮候选拒绝，再以 ratio/bound tie 作为 epoch perturbation 的
  主要触发信号。
- dual Phase-I 首轮失败已由 caller-owned candidate/`ep` trace 闭环：leaving row `4`
  在 `.at_upper` 方向违反 `131.5`，`ep` 只有 `ep[4] = 1`。469 列中 461 列因 tableau
  小元素被排除、1 列为 basic/fixed、没有 wrong-sign rejection；其余 7 列全部是合法的
  boxed bound-flip 候选，因此不是候选符号、状态或 tolerance 实现错误。
- 这 7 列的 working bounds 均为 `[0,1]`、tableau/signed pivot 均为 `1`，每列最大
  correction 为 `1`，总容量仅 `7`；全部 flip 后仍有 `124.5` 的基本上界违反且不存在
  可 pivot 的非 boxed 候选。根因是固定单位 dual Phase-I working bounds 与 logical
  basis 的基本值尺度不匹配，当前 ratio test 的空 entering 结论对该工作问题是正确的。
  禁止强制选择最后一个 boxed 列或放宽 pivot tolerance；若要使该路径可用，必须重新
  设计 scale-aware Phase-I bounds/nonbasic move 表示，或先构造能降低基本违反的 crash
  basis。默认继续使用 primal Phase I，不增加 `brandy` 模型特判。
- basis-epoch virtual perturbation 只改变 Harris/weighted-pricing 容差内的稳定排序，不改变
  原始 step 或 reduced cost；row/column rank 由 tolerance、scale、epoch 和 id 确定，幅度
  上限为对应容差的 2 倍，单 epoch 最多存活 256 个退化 pivot。64-entry generation-marked
  taboo ring 记录 `(entering, leaving, direction)`，有效移动和 epoch 失效会清空 exclusion。
- cleanup 先移除 artificial basis，再 fresh rebuild、精确 Phase-I reprice 并重验 artificial
  objective；失败时事务性恢复 logical basis epoch，并关闭 perturbation 重跑 Phase I。
  `perturb`、`taboo` 两种显式模式分别通过 40 模型 status/objective/residual/ray gate。
- `brandy` 21 轮 ReleaseFast：baseline 为 3384 iterations、3101 degenerate pivots、
  64.65 ms median / 66.03 ms p95；taboo-column 为 1237、947、23.29 ms / 23.81 ms，
  分别下降 63.4%、69.5% 和 64.0%。达到默认候选的 brandy 门槛，但其他模型存在不稳定
  iteration 回退，因此 `.automatic` 暂时仍保持 baseline，禁止按模型名称特判。
- Phase-I A/B 中 Devex 保持默认：Dantzig 在 brandy 增至 9195 次 Phase-I iteration；
  approximate steepest-edge 不降低 900 次 Phase-I iteration，却将 BTRAN 从 899 增至
  4199。三种强制模式均通过 40 模型 gate，保留用于固定 trace A/B。
- adaptive reprice 记录 exact snapshot drift；brandy taboo 单轮最大归一化 drift
  `5.21e-10`。它只允许缩短已验证的 8-update refresh period，且 Phase-I 周期不会泄漏到
  Phase II；显式 adaptive 模式通过 40 模型 gate。
- reusable CSR row view 在 solve 边界构建，pricing 操作按 touched CSR entries 一次分发。
  brandy taboo-row 21 轮为 22.20 ms median / 22.45 ms p95；scsd1 row pricing median 从
  4.03 ms 降至 3.24 ms、总耗时 10.05 ms 降至 9.83 ms，但 iteration 587 -> 652 且 p95
  10.19 ms -> 11.24 ms。因此 row/auto/column 均保留强制模式，默认继续 column；三者和
  taboo+adaptive+row 组合均通过 40 模型正确性 gate。

### 6.6 每阶段提交与验收规则

- [x] 阶段一提交：`bench: add phase-level simplex counters`。
- [x] 阶段二提交：`simplex: classify and reduce rebuild fallbacks`。
- [x] 阶段三提交：`simplex: add dedicated dual phase one`（`1fb004a`）。
- [x] 阶段四提交：`simplex: add sparse triangular crash basis`（`d27edcc`），并以
  `simplex: evaluate bixby crash scoring`（`29f4cbe`）完成评分策略 A/B。
- [x] 阶段五建议提交：`simplex: add epoch perturbation and adaptive pricing`。
- [x] 每次提交前运行单元测试、固定 corpus 差分和对应 benchmark；将结果写入
  可由 git 追踪的报告，报告同时保留测试环境、编译模式和 HiGHS 版本。
- [x] 未达到本阶段正确性门槛时不得以性能提升为理由进入下一阶段；未经 fresh
  rebuild 验证的启发式结果不得用于最终状态或 certificate。

## 7. 优化完成后的完整端到端验收

本节只在第 6 节全部完成后执行，避免 dual Phase I、crash basis 和退化策略改变
pivot 路径后重复生成大规模报告。快速 corpus 仍须在第 6 节每次提交前运行。
第 6.3 节 scale-aware dual Phase-I 重设计已明确移入后续独立优化，不属于本阶段前置
条件；第 7 阶段只验收当前已通过快速 corpus 的默认 primal Phase-I 路径。

- [x] 冻结官方 Netlib 传统压缩 MPS 中可直接解码的 93 个模型及其 SHA-256，新增
  带逐模型 timeout、memory limit、peak RSS、原始 TSV 保留和 HiGHS/可选 CLP 对照的
  可复现 runner；首轮结果记录于 `bench/simplex/stage7_results.md`。
- [x] 完成 93 模型解析前置关卡：修复 fixed MPS 名称内空格及省略 BOUNDS set name
  两类兼容问题，并以 `gfrd-pnc`、`sierra` objective 对照确认未再静默丢失 bounds。
- [ ] 获取并锁定 `stocfor3`、`truss` 与 QAP8/QAP12/QAP15 五个特殊生成模型；在此之前
  “93/93 可解析”不得表述为完整 Netlib 已通过。
- [ ] 处理首轮暴露的 `modszk1`、`scsd8`、`wood1p` 数值失败、`tuff` Phase-I 循环及
  大模型超时；强制 perturbation 已修复 `scsd8`/`wood1p`，但 taboo 在 `tuff`/`wood1p`
  产生错误 infeasible，未通过默认启用门槛。
- [ ] 增加固定版本 CLP runner；当前环境没有 CLP，因此只允许报告 zhighs/HiGHS 的
  阶段性数据，不能勾选三方完整对照。

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
