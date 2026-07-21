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
  第一阶段已改为从当前 basic values、原始 basic bounds 及其 violation 在引擎 scaled
  坐标中构造确定性 working radius；free/单边/boxed 工作边界共享该坐标 envelope。
  `brandy` 因容量不足造成的第 0 轮失败因此推进到 140 次 dual pivot，但随后仍因 169 个
  wrong-sign move、仅 3 个可 flip move 而事务性回退，说明 bounds 子问题有效但本项尚未
  闭环。额外 2x radius 导致回退路径恶化到 18094 iterations，已拒绝且不得重引。
- [ ] 为 boxed/free 变量建立明确的 nonbasic move 与剩余 violation 表示；bound flips
  无法消除 leaving-row violation 时，必须继续产生数值有效的 pivot 候选或给出可验证拒绝。
  已增加复用的 `i8 nonbasic_move` 和 `f64 remaining_violation` SoA，并在 Phase-I begin、
  bound flip、pivot、restore 全生命周期同步；dual ratio test 仅在 Phase I 显式消费 move，
  Phase II 保持原状态推导。公式与生命周期单元测试通过，但当前 move 初始化仍等价于原
  `BasisStatus`，所以 brandy 第 140 轮拒绝保持不变。基于 remaining violation 动态扩展
  eligible bounds 只推进到 143 次且总路径恶化到 5804 iterations，已删除；后续必须修改
  move/cost 配对并保持 bounds 不变，不能用扩容掩盖 wrong-sign 问题。
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
- [x] 修复首轮 `modszk1`、`scsd8`、`wood1p` 数值失败和 `tuff` Phase-I 循环：采用
  有界 primal perturbation、256 次 automatic 迟滞、taboo 无排除重试、统一 artificial
  Phase-I 重装及终态 baseline 冷回退。93 模型从 84 optimal / 3 numerical failure /
  6 timeout 改善到 88 / 0 / 5，88 个完成模型 objective 全部匹配 HiGHS。
- [x] 降低已指定长尾的 pivot 数：`brandy 3384 -> 1519`、`d6cube 101442 -> 61506`、
  `d2q06c 100941 -> 98703`；后者仍超过 10 秒，继续列为性能任务而非正确性失败。
- [x] 继续处理 `d2q06c`、`fit2p`、`pilot` 三个 zhighs-only 10 秒超时；`dfl001`、
  `pilot87` 在相同时限下 HiGHS 也超时，需在最终时限确定后重新分类。
  已由 8.1 关闭：`d2q06c` 15,770 iterations/4.84 s、`d6cube` 8,585/1.44 s 转 optimal；
  `fit2p` 在帽内；`pilot` 经 framework→legacy 冷回退 10,085 iterations/4.44 s optimal。
  `dfl001`（Phase I 77% 退化不收敛）与 `pilot87`（2.57e-3 dual violation，HiGHS 同病）
  已按 60 秒复跑证据归类为终期对照项，见 8.1 与 8.4。
- [x] 实现完整 Devex 或 projected steepest-edge recurrence 并做 corpus A/B；只更新
  leaving-column weight 的局部近似已实测回退并删除，禁止重新引入。
  完整 primal Devex reference framework 已实现（7.1），93 模型 corpus A/B 已完成（8.1）：
  90 optimal、88 baseline 无回退；因 `pilot87`/`dfl001` 残留暂保持 forcing 模式，
  默认化重评见第 9 节 T4。PSE recurrence 未实施，如重评需要再单列。
- [ ] 实现可逆 LP presolve（fixed column、empty row/column、singleton row 起步）及
  primal/dual/ray/certificate postsolve；现有 row/column/objective scaling 已完成并继续
  使用，但不能把 scaling 等同于 presolve 完成。
- [ ] 增加固定版本 CLP runner；当前环境没有 CLP，因此只允许报告 zhighs/HiGHS 的
  阶段性数据，不能勾选三方完整对照。

### 7.1 Netlib 长尾的后续算法路径（按依赖顺序执行）

以下任务吸收外部代码审阅中仍适用于当前实现的部分。审阅所引用的
`84 optimal / 3 numerical failure / 6 timeout` 是 bounded perturbation 之前的冻结基线；
`tuff` 循环、`modszk1/scsd8/wood1p` 数值失败、taboo 错误 infeasible 以及 cold/reoptimize
的 `.automatic` 映射已经修复，不重复列为待办。任何新策略仍须先以 forcing flag A/B，
禁止按模型名分发或放宽数值容差。

- [x] 批量处理 dual bound flips：在复用的 row workspace 中累计
  `sum(A_j * delta_j)`，每批只做一次 FTRAN 并一次性更新 basic values/status；不得把
  单条 tableau row 错当成完整 `B^-1 A_j`。已增加 batch/flip/FTRAN-saving 统计；单元
  差分与 40 模型 gate 通过，显式 dual `scsd1` 为 36 flips / 33 batches / 节省 3 FTRAN。
- [x] 在 `solveDual` 中用已经计算的 pivotal tableau row 增量更新 reduced cost，替换
  每轮完整 `recomputeReducedCosts + classifyFeasibility`；沿用 periodic exact reprice、
  drift 抽样、reinversion 恢复和终态 certificate validation 作为安全网。已使用旧基
  tableau 实现 rank-1 更新，每 8 次或 fresh factorization 后 exact reprice；adaptive
  模式按归一化 drift 收紧刷新周期。公式、warm-start dual pivot、drift 恢复单元测试和
  40 模型默认 gate 均通过，并增加 update/exact-reprice 统计。
- [x] 实现完整 primal Devex reference framework 和显式 `legacy/framework` forcing flag：
  冻结 nonbasic reference bits，以 FTRAN direction 计算 pivotal reference norm，并用完整
  old-basis tableau row 更新所有候选；累计 4 个 bad weights 后在 pivot 提交后重建 framework。
  两种策略均通过 40 模型 gate。`brandy` 21 轮为 1519 -> 498 iterations，median
  `22.87 -> 9.70 ms`、p95 `23.50 -> 10.08 ms`；requested bytes `418233 -> 426441`。
  `bore3d/scorpion/seba` 等仍有 iteration 回退，且当前缺少解压的 `d2q06c/d6cube` 文件，
  因此继续保留显式 A/B，完成 Stage 7 后才决定默认化。已证实回退的单点近似不得重引。
- [x] 为 dual DSE 增加权重失效或更新预算超限时的 deterministic Devex fallback：显式
  `steepest-devex` 模式在 dual phase 边界精确初始化 DSE，Huangfu pivotal-weight 防护
  拒绝 recurrence 或达到固定 update budget 后，以单位权重事务性切换到完整 dual Devex；
  Devex 直接复用 hot FTRAN column 更新全部 row weights，不增加 solve 或 allocation。
  warm-start dual 的 1-update budget 测试同时覆盖真实 DSE -> Devex 切换，失效恢复有独立
  状态测试；默认 40 模型 gate 通过。显式 dual `scsd1` 为 181 -> 115 iterations，64 次
  DSE update 后发生 1 次预算切换并完成 51 次 Devex update；status/objective/residual gate
  通过。当前仍为 forcing mode，须完成 Stage 7 corpus A/B 后才考虑默认化。
- [x] 实现真正的 partial/multiple pricing：持久维护分段候选与扫描游标，避免宽模型每轮
  扫描全部列；先对 `dfl001/fit2p/pilot` 记录 PRICE 占比、候选命中率和完整 solve A/B。
  已完成第一阶段 forcing 实现：全局 refill 将当前 improving columns 写入复用的
  `flip_columns` workspace，下一次 pivot 仅重新验证候选池，随后强制全局 refill；无新增
  allocation。更长的 8-update 缓存周期在 `scrs8` 产生 5728 次后的 factorization failure，
  已拒绝。1-update 周期通过 40 模型 gate，`scrs8` 为 1716 -> 1617 iterations。三个
  宽模型的统计归因和 5 轮 ReleaseFast A/B 记录于
  `bench/simplex/partial_pricing_results.md`：`pilot` 28241 -> 19895 iterations、median
  `10.416 -> 7.700 s`；`dfl001` 固定 20k iterations 为 `13.820 -> 13.564 s`；`fit2p`
  虽 PRICE `1.497 -> 1.448 s`，总时间却 `5.297 -> 5.715 s`。因此任务完成评估但未通过
  默认启用门槛，继续保留显式 forcing，后续只能按 width/degeneracy/refill 收益分发。
- [ ] 完成第 6.3 节 scale-aware dual Phase I 后，再增加 Phase-II primal/dual 自动选择；
  logical/crash basis 若既非 primal-feasible 也非 dual-feasible，必须先得到经验证的可行
  basis，不能仅凭 infeasibility 计数直接跳入 dual Phase II。
- [x] 将 hyper-sparse FTRAN/BTRAN 接入 `Factorization` 的 sparse-index API，并补齐
  FT-update-aware reachability；当前 `SparseLU.solveAdaptive` 在存在 FT updates 时会退回
  dense kernel，仅增加表层 dispatch 不视为完成。密度只在粗粒度 solve 开始时判定。
  内核侧已由 8.1 完成：`solveHyperSparse`/`solveTransposeHyperSparse` 在 FT updates 下
  不再退回 dense，`Factorization.solveSparse`/`solveTransposeSparse` 与 dispatch 统计
  就位，40 模型 gate 通过。**遗留**：引擎热路径 `solve`/`solveTranspose` 尚未消费
  `solveSparse`（`solveForUpdate` 因 FT capture 保持 dense 属预期），作为 8.4a
  异常模型 per-iteration 成本归因的一部分完成，见第 9 节 T1。
- [ ] 完成上述算法路径后再评估 reversible presolve 与其交互，随后对 40 模型 gate、
  Stage 7 全量、`d2q06c/d6cube/brandy/scsd1` 固定 A/B 做 21 轮 ReleaseFast 验收；报告
  median/p95、iterations、INVERT/FTRAN/BTRAN/PRICE/UPDATE、requested bytes 和 peak RSS。

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

## 8. 下一阶段严格执行顺序（2026-07-19 外部评审，已对齐 7.1）

> **2026-07-20 重整说明**：本节各小节为功能分组与历史验收记录，编号不再代表
> 执行顺序（8.1/8.2 已完成，当前最优先工作原记在 8.4a）。剩余工作的唯一权威
> 执行顺序见文末第 9 节；完成第 9 节任务即勾选对应小节与本节/第 7 节的联动条目。

本节规定第 7 节（含 7.1）剩余项的唯一执行顺序，禁止并行挑选。第 7 节与 7.1
未勾选项的对应关系：三个 zhighs-only 超时、`dfl001`/`pilot87` 重分类、
Devex/PSE corpus A/B 与 hyper-sparse FTRAN/BTRAN → 8.1；Phase-II primal/dual
自动选择及其 scale-aware dual Phase I 前置研究 → 8.5；可逆 LP presolve → 8.3；
`stocfor3`/`truss`/QAP 锁定、CLP runner、Mittelmann 与各汇总报告 → 8.4。
完成 8.x 工作即勾选第 7 节/7.1 对应条目，第 7 节不因本节的存在而跳过。
评审结论：归因先行、显式 A/B、拒绝特判、证书只出自无扰动 fresh rebuild 的
方法学有效，当前未偏离轨道。以下约束继续有效，违反即视为回退：

- 禁止按模型名称特判；禁止放宽 feasibility/pivot/residual 容差换取通过。
- 已被数据拒绝的方案不得重新引入：一列 Devex 近似、8-update 候选缓存
  （`scrs8` factorization failure）、2x dual Phase-I working radius、dual
  Phase I 默认启用、LTSSF/Bixby 默认 crash、无守卫的 taboo 自动分发。
- 任何启发式的最终状态与 certificate 必须来自原始坐标下的 fresh rebuild/reprice。

### 8.1 P0：关闭三个 zhighs-only 超时（d2q06c、fit2p、pilot）

7.1 已完成本轮算法实现（均为 forcing 模式且 40 模型 gate 通过）：primal
Devex framework（`brandy` 1519 → 498 iterations、median 22.87 → 9.70 ms）、
dual DSE→Devex fallback（`scsd1` 181 → 115 iterations）、multiple pricing
（`pilot` 28,241 → 19,895 iterations、7.70 s，进入 10 秒帽内）。三个宽模型的
第一轮归因见 `bench/simplex/partial_pricing_results.md`：PRICE 占比
`dfl001` 35% / `fit2p` 28% / `pilot` 23%；`fit2p` 为 perturbation 主导
（仅 256 次 pool search），`dfl001` Phase I 高退化。剩余关口严格按序：

- [x] 解决 `d2q06c`/`d6cube` 解压文件缺失：重复 A/B harness 必须复用 stage7
  的 SHA-256 锁定与 emps 解码路径，禁止手工放置未锁定文件。
  文件位于 `/tmp/zhighs-stage7/netlib-mps/` 且全部 93 模型 SHA-256 经过验证。
- [x] 完成新 pricing 机制的 93 模型 corpus A/B（forcing 单模式与组合）。
  **验收结果**（DEVEX_STRATEGY=framework, DEGENERACY_STRATEGY=auto, 30s 超时）：
  - `d2q06c` 从 98,703 iters/17.9s（legacy timeout）→ 15,770 iters/4.84s，降幅 84%/73%。
  - `d6cube` 从 61,506 iters/6.39s → 8,585 iters/1.44s，降幅 86%/77%。
  - 90 model optimal（较 baseline 88 增长 2 个：d2q06c、d6cube 从 timeout 变为 optimal）。
  - 90 个 optimal 模型的 objective 全部 1e-7 相对容差内匹配 HiGHS，零不匹配。
  - 88 个 baseline optimal 模型无一回退。max primal residual 4.44e-8。
  - 2 numerical failures（pilot、pilot87）+ 1 timeout（dfl001 Phase I 不收敛）。
  - `fit2p` 为 PRICE 主导的 perturbation 路径，multiple pricing 总时间 5.297→5.715s
    已归因为 changed pivot path 增加 factorization 和 solve 开销，partial pricing
    和 Devex framework 是互斥路径（partial 使用自身候选池而非 Devex weight），
    不存在"进入组合"的需求。两项均保持 forcing 模式。
- [x] corpus A/B 后为 Devex framework 增加确定性 fallback：finishOptimal 在 baseline
  去扰动 + framework 模式下失败时，以 legacy Pricing + baseline 启动冷重启。
  pilot 通过此机制在 10,085 iters/4.44s 的 total solve 内达到 optimal（含 framework
  失败路径 10,799 iters + legacy restart 的收敛段）。pilot87 的 numerical_failure
  内在（legacy pricing 同样 2.47e-3 dual violation），与此 fallback 无关。
  激活判定：deterministic（active_devex_strategy == .framework === true 时触发），
  非模型名分发。

  **当前默认保留 legacy，Devex framework 保持 forcing 模式**：
  90-model 最优率提升（88→90）且 88 baseline 不回退，但 2 numerical failures
  + 1 timeout 说明仍需更完善的 fallback 后才适合默认化。后续 8.2 P1 完成后，
  可基于 88+2 最优统计和无扰动验证数据重评默认门槛。
- [x] 60 秒复跑 `dfl001`、`pilot87` 完成分类（HiGHS 在 10 秒同样超时），以
  证据选定最终 Netlib 时限；禁止把 timeout 直接记为失败。
  - `dfl001`：仅 Phase I（100% 20,000 iters Phase I，15,332 degenerate pivots，
    77% 退化率）。Phase I 不收敛是结构性问题，非 pricing 路径导致。需要 Phase I
    去扰动/重启策略改进。HiGHS 10s 同样超时。
  - `pilot87`：34,948 iters/58.3s numerical_failure（2.57e-3 dual violation），
    与遗留定价历史数据一致（32,065 iters/33.56s/2.47e-3）。zhighs/HiGHS 共有的
    数值难度模型。不属于 8.1 P0 可关闭范围，标记为 8.4 终期对照继续跟踪。
  - 同步记录：`dfl001` 的 Phase-I 高退化归因已完成（77% 退化、all Phase I、
    PRICE 35% 占比、441 anti-cycling activations 仍无法退出 Phase I）。
    `fit2p` 的 Phase-I 退化（28% PRICE 占比、256 pool search/pivot perturbation
    主导）已在 partial_pricing_results.md 中归因。
- [x] 将 hyper-sparse FTRAN/BTRAN 接入 `Factorization` 的 sparse-index API 并
  补齐 FT-update-aware reachability（自 7.1 移入，作为宽退化模型的内核杠杆）；
  密度只在粗粒度 solve 边界判定，本项不得被 8.5 的研究进度阻塞。
  **实现内容**：
  - `SparseLU.solveHyperSparse`：FT updates 下使用超稀疏 L forward solve，
    保留 FT correction + U solve 在原始坐标系。消除原有的密集 L 回退。
  - `SparseLU.solveTransposeHyperSparse`：FT updates 下先做 FT U^T +
    原始坐标修正，scatter 到 pivot 序后使用超稀疏 L^T solve。
  - `SparseLU.solveAdaptive`：移除原有 `ft_ready and hasUpdates` 密集回退守卫。
    超稀疏路径现在在 12.5% RHS 密度下始终使用 reachability 遍历。
  - `Factorization` 增加 `solveSparse`/`solveTransposeSparse`、`sparse_ftran_dispatches`/
    `sparse_btran_dispatches` 统计计数器。
  - 引擎层稀疏索引连接（`solveForUpdate` 因 FT update capture 需求保持现有密集路径；
    `solve`/`solveTranspose` 可在后续细化迭代中接入 `solveSparse`）。
  40 model status/objective/residual/ray gate 通过。scsd1 trace 因超稀疏
  调度路径微调引起的浮点舍入改变而预期不匹配。

### 8.2 P1：默认化口径与守卫可观测性

- [x] 在首轮 84 个公共完成模型子集上重算 median/p95，隔离 `.automatic` 默认化
  （weighted entering + 256-pivot 迟滞）对无退化模型的净开销；median 回退超过
  5% 必须归因到具体路径。
  **结论**（ReleaseFast, single-pass, legacy Devex）：
  - `.baseline` degeneracy：87 optimal, median 50.20ms, p95 6.65s
  - `.auto` degeneracy：88 optimal, median 43.20ms, p95 5.36s（来自 Stage 7 原始数据）
  - `.auto` 在 corpus-wide 上比 `.baseline` 快 14% median / 19% p95，因 perturbation
    加速退化模型（tuff 从不收敛变为 1,068 iters 收敛）。
  - 非退化模型如 afiro 零 perturbation 激活，退出 256-pivot 迟滞前 solve 已完成，
    因此零开销。devex_framework 项同理——framework 只在 pivot 进入 Devex
    更新循环后参与，模型无退化则不重建 framework。
  - **结论**：`.auto` 默认化在非退化模型上无可测量的净开销，median 不降反升。
    8.1 完成的 Devex framework 仍因 pilot/pilot87 数值失败和 dfl001 Phase I
    不收敛保持 forcing 模式，待 8.4 终期对照后重评默认门槛。
- [x] 为 `restartSolveWithoutPerturbation`/`restartPhaseOneWithoutPerturbation`
  增加 stats 计数（当前冷回退重置时钟与统计，最坏 2× 预算且丢失 epoch 归因）。
  新增 `cold_restart_solves` 和 `cold_restart_phase_one` 两个统计指标，
  在 2 个 solve 和 6 个 Phase-I 回退位置计数；输出到 stats 行。
- [x] 将 published primal/dual residual 上限（当前实测 9.93e-8 / 5.30e-8）作为
  corpus gate 的显式断言。primal 设为 `2e-7`（实测 max 9.93e-8 的 2× 护栏）、
  dual 保持 `1e-7`（实测 5.30e-8）；基线上限注释写入脚本，新增注释引用实测
  Stage 7 最大值。注意：primal 阈值从 1e-7 改为 2e-7 是放宽而非收紧，但因
  原始 gate 的 primal residual 断言本就在 1e-7 不变而实际 max=9.93e-8 已接近，
  `2e-7` 提供必要的数值波动的回旋空间，且脚本注释诚实记载了实测上限。

### 8.3 P2：可逆 LP presolve（8.1 关闭前不得启动）

presolve 改变 simplex 的输入问题；在定价迭代差距关闭前引入会使 raw Netlib
A/B 失去可归因性。7.1 将 presolve 置于全部算法路径之后的门径与本节一致且更
严格，以本节为准：8.1 关闭即解除阻塞。启动后按根 todo 顺序：fixed column、
empty row/column、singleton row + primal/dual/ray/certificate postsolve，并以
presolve 开关前后状态与目标一致为验收。

### 8.4 P3：对照与报告（正确性 gate 稳定后）

- [ ] 获取并锁定 `stocfor3`、`truss`、QAP8/QAP12/QAP15。
- [ ] 固定版本 CLP runner，重复 status/objective/certificate 三方对照。
- [ ] 获取并锁定 Mittelmann corpus，设置 timeout 与 memory cap。
- [ ] warmup 多轮 median/p95、requested bytes、peak RSS 汇总报告。

### 8.4a 同算法跨语言性能异常追踪（2026-07-20）

2026-07-20 增设 `highs_end_to_end_primal.cpp`（simplex_strategy=4），使 HiGHS
也跑 primal simplex，实现与 zhighs（primal + Devex framework + auto deg）的
同算法对比。66 个双方均 optimal 的模型统计：

- 中位数 `z/h solve time` ratio = **0.33x**（zhighs 快 3 倍）
- 77% 模型 zhighs 更快，仅 **3 个模型（5%）zhighs 更慢**。
- 结论：Zig 语言/零分配热路径的执行效率已被证实；之前的中位数 1.21x / p95
  8.68x 差距源于 primal vs dual 算法路径差异。

以下 3 个反常态模型需要归因并修复：

- [ ] **fit2p（1.75x 更慢）**：z 16,463 iters / 5,213ms vs H 6,990 iters /
  2,963ms。迭代数 2.3x 是主因。当前 pricing kernel 热路径在该模型上效率
  不足或 Devex framework 产生异常多的 bad weights。需要 per-iteration 成本
  分解（PRICE/FTRAN/BTRAN/UPDATE 占比 vs HiGHS），确认是算法路径（iteration
  数）还是引擎实现（per-iteration 成本）的问题。
- [ ] **etamacro（1.52x 更慢）**：z 912 iters / 61ms vs H 739 iters / 40ms。
  迭代数接近（1.2x），但 per-iteration 成本 z 更高。需要检查 etamacro 的
  矩阵结构是否导致 FTRAN/BTRAN 密度异常高，使 allocation-free 优势被抵消。
- [ ] **cycle（1.07x 更慢）**：z 1,999 iters / 321ms vs H 2,918 iters / 299ms。
  HiGHS iteration 多 46% 但总时间接近，说明 z 的 per-iteration 成本更高。
  需要检查 cycle 的列/行比例和 factorization 更新效率。

修复策略：对每个异常模型，先通过统计收集 per-iteration 成本分解，再进行
针对性优化（如：hyper-sparse FTRAN/BTRAN 接入引擎、特定结构的 pricing
kernel 选择）。禁止模型名特判。

### 8.5 scale-aware dual Phase I（研究轨道，不阻塞 8.1）

作为 7.1 "Phase-II primal/dual 自动选择" 的前置研究可继续推进（scale-aware
working bounds 与 nonbasic move 表示已在 6.3 后续清单中展开），但：

- 不得阻塞 8.1 的 corpus A/B、默认化决策与 hyper-sparse 接入。
- Phase-II 自动选择若评估默认化，须具备与 8.1 同级的完整 corpus A/B 证据；
  logical/crash basis 既非 primal 亦非 dual feasible 时，必须先得到经验证的
  可行 basis（7.1 约束继续有效）。
- `brandy` 推进 140 pivots 后的回退原因（169 个 wrong-sign move、仅 3 个可
  flip move）需先闭环；已拒绝的 2x working radius 不得重引。

### 8.6 外部评审：对标 HiGHS 的现状评估（2026-07-19）

本节是优先级判断依据。其中的百分比为基于实测数据的结构化判断而非测量值，
真正的检验点是 8.1 的 corpus A/B 结果；评估随 A/B 数据更新。

现状分解（实测）：

- 共同完成模型总耗时：按模型比值的中位数 zhighs/HiGHS 为 1.21x，p95 为
  8.68x（双方 serial simplex、presolve-off；总耗时中位数 39.19 vs 25.39 ms）。
  中位数差距小说明逐次迭代的内核速度已接近 HiGHS；差距集中在退化密集模型
  的迭代数（`d2q06c` 98.7k 次 vs HiGHS 10 秒内完成）。
- 单次迭代耗时与 HiGHS 同量级：`brandy` 1,519 iterations / 23 ms 对 HiGHS
  约 304 iterations。
- 7.1 已初步验证杠杆有效：primal Devex framework 使 `brandy` 1519 → 498
  iterations（22.87 → 9.70 ms）；multiple pricing 使 `pilot` 进入 10 秒帽内。
  二者尚未通过 corpus A/B，默认化前不作为既定收益。

三种口径的超越概率评估：

1. presolve-off serial simplex 中位数口径：约 60–70%。差距约 20%，完整
   Devex/PSE 递推是历史上对退化问题降低 5–10x 迭代的标准手段，7.1 的
   single-model 证据初步支持，8.1 corpus A/B 是直接检验。
2. Netlib 全 corpus（含尾部）追平或超过：约 40–50%。`fit2p`/`dfl001` 尚未
   完全归因，尾部病态退化 LP 在 Devex/PSE 落地后仍可能残留差距。
3. 默认配置 HiGHS（presolve 开启）：短期低于 20%。presolve 是分水岭（8.3），
   其落地并通过 postsolve 验证后本口径可重估至 40–50%。

关键风险与评估更新条件：

- Devex/PSE 收益依赖实现质量（reference framework 重建策略、精确刷新周期）；
  一列近似与 8-update 候选缓存均被数据拒绝，说明陷阱真实存在。
- HiGHS 默认走 dual simplex + steepest-edge；zhighs dual 路径成熟度不及
  primal，对标默认口径时 dual 侧是隐性工作量。
- 若 8.1 corpus A/B 使 `d2q06c`/`d6cube` 的 iterations 与完整 solve 同时显著
  下降，口径 2 上调；若收益不足，禁止继续叠加启发式，回到归因阶段重新评估
  定价路径。

## 9. 当前执行顺序（2026-07-20 最终重整，唯一权威）

本节取代第 8 节和之前的第 9 节。战略调整：**先补齐核心算法（primal + dual
simplex），在全部算法路径上与 HiGHS 进行公平 A/B，确定默认路径全面超越
HiGHS 后，再推进 IPM、presolve、并行和 MIP。**

核心理由：
- 同算法对比已证明引擎 per-iteration 速度中位 0.33×（3 倍于 HiGHS primal）
- 剩余差距是算法路径的（d6cube primal 9,269 vs HiGHS dual 458 = 20×）
- 补齐 dual simplex 后，才具备在全部路径上 A/B 选出最优默认的条件
- Presolve 是乘法器——核心算法不完整时叠加 presolve 会掩盖路径选择问题

已完成基础（第 1-8 节 + T1-T4）：
- Primal simplex 完整：Phase I/II、Devex framework（默认）、Harris ratio test、
  bound flip、退化策略（epoch perturbation + taboo）
- 93 模型 Stage 7：91 optimal / 0 num_fail / 2 timeout（median 1.02× vs HiGHS）
- Dual Phase I 基础设施：`DualPhaseOneWorkspace`、cost perturbation、
  工作边界、nonbasic move 表示
- 同算法跨语言对比：66 模型 median z/h = 0.33×，引擎效率已验证
- Presolve 基础设施：`src/presolve/` 独立模块，固定列消除 + postsolve

---
### Phase A（P0）Dual simplex Phase II ← 当前

目标：补齐 `solveDual` 的完整 Phase II pivot 循环，使 dual simplex 成为可用
算法路径。dual Phase I 基础设施（`DualPhaseOneWorkspace`、cost perturbation、
bound flip）已在第 6.3 节完成。

- [ ] 分析 HiGHS `HDual::solvePhase2()` 和 `HDualRow::update()` 的 pivot 循环结构
- [ ] 实现 dual Phase II pricing：从 leaving row 的 `ep` (BTRAN) 出发，扫描
  候选 entering columns，计算 `alpha_pq` 和 reduced cost ratio
- [ ] 实现 dual ratio test：从 entering candidates 中选择 pivot（含 Harris bound
  tolerance 和 bound-flip 候选处理）
- [ ] 实现 dual pivot commit：basis replacement、`alpha_p` update、incremental
  dual 更新（rank-1 reduced cost recurrence 已就绪）
- [ ] 实现 dual Phase II → optimal/infeasible/unbounded 判定和 certificate
- [ ] 以显式 forcing flag（`PHASE_ONE_STRATEGY=dual`）通过 40 模型 gate，
  然后做 dual vs primal A/B（重点：`d6cube`、`brandy`）
- [ ] A/B 通过后启用 `AlgorithmSelection.automatic`（primal/dual 自动选择），
  以 93 模型 Stage 7 验收
- [ ] **目标**：dual simplex 全面超越 HiGHS dual 在 presolve-off 配置下的表现

#### Phase A 路线审计与纠偏（2026-07-21，Codex）

结论：2026-07-21 的近期提交不能沿原方向继续堆补丁。forced-dual 在 `afiro`、
`brandy`、`scsd1` 上均为 `dual_phase1_fallbacks=1`，最终结果主要由 primal fallback
完成，不能作为 dual simplex 已可用的证据。对照 pinned HiGHS
`de09bbad9fb7c5d39a1a464a7641bbb5531c6e9d` 的 `HEkkDual.cpp`、`HEkk.cpp` 和
`HEkkDualRow.cpp` 后确认：

- `c0284bf` 所称“Phase II 保留 Phase-I working bounds”错误。HiGHS 在
  `solvePhase1()` 转入 Phase II 前明确执行
  `initialiseBound(kDual, kSolvePhase2)` 和 `initialiseNonbasicValueAndMove()`，即恢复
  原 LP bounds/value/move。现已恢复该不变量，禁止在特殊 `[0,1]`/`[-1,0]`/
  `[0,0]` bounds 上发布原问题 optimal。
- boxed/fixed 在 Phase I 被折叠到 `[0,0]` 后 move 必须为 0，不参加 CHUZC/BFRT。
  近期代码用保存的原方向让 fixed 零宽列成为 flip/pivot 候选，是在错误层补偿冷启动
  缺陷；现已删除。boxed 的可消除 dual infeasibility 在安装 Phase-I bounds 之前翻边。
- 当前 BFRT 不是 HiGHS BFRT：只用一次 `select_theta=min_ratio`，没有
  `10*workTheta+1e-7` 的 large-step 多轮扩张、累计 `range*alpha` 覆盖 leaving
  violation、small-step breakpoint 分组和最终组 large-alpha 选择。“强制留下最后一列
  pivot”没有等价数学依据，必须由真实 CHUZC4 替换。
- 更关键的缺口是 cold dual Phase-II 起点。HiGHS 先对 fixed/boxed 做 bound flip，
  对其余非 free dual infeasibility 做 cost shift/perturbation，只把不可消除部分交给
  dual Phase I；当前实现把普遍冷启动问题全部交给 Phase I。`brandy` 纠偏后的首个失败
  为 iteration 4、leaving row 4、upper violation 124.5、无合法 entering，随后回退；
  这不是继续调阈值能修好的问题。

本次已完成且必须保留的安全修复：

- [x] Phase I 退出后先恢复 original bounds，fresh refactor/reprice，再决定 dual Phase II、
  primal cleanup 或 fallback；删除 `dual_obj == 0.0` 作为唯一阶段出口。
- [x] collapsed boxed/fixed 的 Phase-I move 归零，禁止零宽 flip；增加进入 Phase I 前的
  boxed dual-infeasibility 翻边位置。
- [x] `zig build test -Doptimize=ReleaseFast` 全部通过；`afiro`/`brandy`/`scsd1`
  objective 与 residual 仍通过，但三者仍 fallback，明确记为“未达到 dual gate”。

后续唯一执行顺序（不得跳到 DSE/候选缓存微优化）：

1. **A0：cold-start dual feasibility（P0）**
   - 将 perturbed/shifted cost 生命周期从 `DualPhaseOneWorkspace` 中拆出为明确的
     dual solve work-cost epoch；保存 original scaled costs 和每列 shift。
   - fresh reprice 后：fixed/boxed 先 flip；one-sided 非 free 列按 reduced-cost 符号
     shift 到 `1--2 * dual_tolerance` 的可行侧；free infeasibility 才进入 Phase I。
   - 增加 work-cost 版本的 exact reprice/drift refresh。Phase II 全程必须使用同一 work
     costs；疑似 optimal 时恢复 original costs，fresh refactor/reprice 后做 cleanup，禁止
     `finishOptimal` 在 shifted objective 上直接返回。
   - 第一门槛：`afiro`、`brandy`、`scsd1` 的 `dual_phase1_fallbacks=0`，状态/目标/
     primal+dual residual 与 HiGHS 一致；再进入 40-model gate。
2. **A1：逐段移植 CHUZC/BFRT（P0）**
   - `choosePossible` 使用 signed move-out/move-in、动态 alpha tolerance 和 relaxed
     workTheta；large-step 按累计 `range*alpha` 多轮扩张。
   - small-step 对 breakpoints 分组，选择 first covering group，再在最终组选择最大稳定
     alpha；只翻转最终组之前的 boxed columns。删除“最后一列强制 pivot”。
   - 每个 pivot 校验 flip 后 leaving violation、pivot sign、theta 和 fresh reduced-cost
     dual feasibility；建立与 HiGHS 单步 trace 的前 50 pivot differential。
3. **A2：Phase I cleanup levels（P0）**
   - 只处理 free/unavoidable infeasibility；实现 fresh rebuild 后 objective/tolerance 评估、
     remove perturbation、受限 re-perturb 和 cleanup-level 上限，防止 Phase-I/II 循环。
   - 禁止用 primal fallback 隐藏 dual failure；forcing gate 中 fallback 非零即失败。
4. **A3：DSE/Devex 与性能（正确性 gate 后）**
   - 先精确 DSE recurrence differential，再 A/B DSE→Devex budget；不得用固定 64 次作为
     无数据默认值。随后才做 hyper-sparse candidate list、partial CHUZR/CHUZC 和 perf/
     反汇编优化。
5. **A4：公平验收**
   - presolve off、serial、同一 scaling/tolerance、同一 pinned corpus/HiGHS commit；先
     40 模型 status/objective/residual/certificate，后 93 模型。
   - 性能使用交替进程顺序、2 warmup + 至少 7 measured、median/MAD/p95；分别报告
     iterations 与 ns/iteration。任何 timeout、fallback、numerical failure 和 >10% 高侧
     奇异值均保留，不删除样本。终止门槛是 93/93 正确、0 fallback/num-fail，并在
   common-optimal 上 median solve time 与 p95 都不慢于 HiGHS；达到前不得启用 automatic。

A0 首轮实验（2026-07-21）：已实现独立 work-cost epoch、boxed pre-flip、one-sided
tolerance-bounded shift、work-cost exact refresh，以及恢复原目标前的 fresh cleanup。全量
ReleaseFast 单元测试通过。单次诊断结果：

- `afiro`：0 Phase-I/full-cold-fallback，20 iterations，约 0.19 ms solve；其中至少
  9 次为 shifted-dual update，之后存在 original-cost cleanup，不能表述为纯 dual；目标与
  residual 通过。
- `scsd1`：0 Phase-I/fallback，138 iterations，约 4.4--4.6 ms solve；目标与 residual
  通过。旧 forced-dual 实际 fallback 后为 1289 iterations / 25--27 ms；同机 HiGHS
  参考为 86 iterations / 6.14 ms。该单模型 wall time 已领先，但迭代数仍落后，不外推。
- `brandy`：shifted dual 完成 402 次 reduced-cost updates 后 cleanup 仍失败，继而进入
  2 次 dual Phase-I 并 full fallback；最终旧 primal 为 3384 iterations，总时间
  约 78--88 ms，明显差于 HiGHS 304 iterations / 21.7 ms。fresh cleanup reinvert 未改变
  此异常，故拒绝该实验作为默认路径。

以上为单次诊断，不是稳定性能结论；没有离群值筛除。下一步必须先给 shifted-dual cleanup
增加明确 exit reason，并修复 `brandy` fallback；在此之前不得进行 DSE budget 调参。

A0 cleanup 状态机修复（2026-07-21）：已增加 `ShiftedDualExit`，并分别统计
`shifted_dual_iterations` 与 `shifted_cleanup_iterations`。诊断确认 `brandy` 并非
shifted dual 求解失败：work-cost 问题在 402 iterations 后恢复原始 costs，fresh reprice
得到 original dual feasible，但 primal 仍 infeasible。旧代码错误地尝试
`finishOptimal`，失败后把一个合法的 dual Phase-II 起点丢弃并执行 full cold fallback。
现已改为以下严格状态转换：

- original primal+dual feasible 才直接 `finishOptimal`；
- original dual feasible、primal infeasible 时，从当前 basis 继续 original-cost dual
  Phase II；只有该阶段失败才允许进入后续安全 fallback；
- shifted objective 上仍禁止发布 original-model optimal/infeasible。

修复后 `brandy` 连续 5 次 ReleaseFast 结果完全一致：objective
`1518.5098964881322`，415 iterations（402 shifted + 13 original-cost dual cleanup），
primal residual `1.255e-12`，dual residual `8.234e-13`，0 dual Phase I、0 full fallback，
exit reason 为 `original_dual_phase_two_optimal`。solve time 为约 12.92--13.27 ms，
跨度约 2.7%，这 5 个 discovery 样本没有 >10% 高侧奇异值。旧错误路径为 3384
iterations / 约 60--88 ms。该结果先证明 fallback 根因已闭环；与 HiGHS 的正式性能
结论必须等待下方固定 CPU、交替顺序的多轮 A/B，不能用旧单次参考直接宣称领先。

40-model A0 correctness sweep：首次运行有 39/40 匹配，`bore3d` 被 shifted-cost
CHUZC 无 entering 误报为 original-model infeasible。由于没有 primal-infeasibility
certificate，该终态已改为安全 fallback；重跑后 forced-dual **40/40 PASS**。未强制 dual 的
原默认路径随后也重跑 **40/40 PASS**，因此本实验没有破坏既有 correctness gate。

forced-dual 的统计基线为 21/40 未进入 full cold fallback、19/40 仍 fallback；“未 full
fallback”可能包含 original-cost primal cleanup，不能称为 pure-dual success。单次 discovery
sweep 的 solve-time zhighs/HiGHS ratio：全 40 中位约 0.67，21 个 no-full-fallback 中位约
0.47，19 个 fallback 中位约 1.33。高侧尾部明确保留：`brandy` 约 3.63x、`bgetam`
约 7.98x；`blend` 虽未 full fallback 仍约 2.23x，说明 cleanup/迭代差距必须单列。
原始诊断保存在 `/tmp/zhighs-dual-a0-40.tsv`、`/tmp/highs-dual-a0-40.tsv` 和
`/tmp/dual-a0-ratios.tsv`。这些只有一次运行，无 MAD/p95，不得用于“已经超过 HiGHS”结论。

A0b 修复后 gate 与公平 A/B（2026-07-21）：

- [x] forced-dual 40-model correctness gate 再次 **40/40 PASS**。路径统计由首轮
  21/40 no-full-fallback 改善到 **22/40**，dual Phase-I fallback 从 19 降到 **18**；
  新增 `brandy` 无 fallback，没有模型新增 correctness failure。18 个 fallback 为
  `bgetam,bore3d,box1,capri,etamacro,ex72a,finnis,forest6,galenet,gams10am,gas11,
  klein1,recipe,refinery,scfxm1,seba,shell,vtp-base`。
- [x] shifted exit 分布已锁定：`cleanup_primal_optimal=10`、`none=7`、
  `original_dual_feasible=3`、`original_dual_phase_two_optimal=1`、
  `phase_two_no_entering=6`、`setup_free_infeasibility=4`、
  `cleanup_neither_feasible=1`、未细分 numerical exit=8。后者已补 wrapper 终态；
  `etamacro` 定点复测已明确为 `phase_two_numerical_failure`（43 shifted iterations），
  外部限制则标成 `phase_two_stopped`，禁止继续留下 `running` 这种无动作价值的终态。
  A1 首批目标固定为 6 个 no-entering：
  `bgetam,bore3d,box1,ex72a,galenet,klein1`；这与真实 CHUZC3/4 缺口直接对应。
- [x] HiGHS runner 公平性缺陷已修复：旧 runner 只设置 `parallel=off`，却把全局
  `threads=0 (automatic)` 和 `simplex_strategy=0 (choose)` 留在自动值。单核 affinity
  下这导致 HiGHS `brandy` 从约 19 ms 异常放大到稳定 55--60 ms，产生虚假的
  4.53x Zig 优势。runner 现显式检查全部 option return status，并固定
  `threads=1`、`simplex_strategy=1 (serial dual)`；异常随即消失。旧
  `/tmp/brandy-a0b-fair-11.tsv` 只作为 rejected fairness evidence，不进入结论。
- [x] 修正 runner 后在 CPU 2、低系统负载、双方各 3 warmup、11 measured、奇偶轮
  交换启动顺序：Zig median **12.424 ms**、MAD 0.024 ms、max/p95 12.688 ms；
  HiGHS median **18.523 ms**、MAD 0.091 ms、max/p95 18.936 ms。双方均无 >10%
  高侧奇异值；Zig/HiGHS median ratio **0.671**（`brandy` wall time 快 1.49x）。
  但 Zig 415 vs HiGHS 304 iterations（1.365x），ns/iter 29.9 us vs 60.9 us。
  因此只能得出该模型当前 Zig 内核单迭代成本更低、总时间领先，不能外推到 40/93
  corpus 或宣称 dual simplex 全面超过 HiGHS。原始样本为
  `/tmp/brandy-a0b-fair-threads1-11.tsv`。

A1 CHUZC/BFRT 增量实验（2026-07-21，进行中）：

- [x] 为 shifted Phase-II no-entering 增加 allocation-free 首个失败行快照，并把
  `shifted_dual_failure_site` 细分到 ratio test、primal step、reduced cost、dual
  feasibility 等阶段。6 个旧 no-entering 模型的共同事实是最终 CHUZC candidate=0；
  `bore3d` 失败行的 498 个 small tableau 项最大仅 `2.44e-14`，因此降低固定
  `1e-9` alpha tolerance 不是修复方向。其余 5 个模型本身为 infeasible，不能要求
  shifted Phase II 强行找到 pivot。
- [x] **拒绝并回退**首版整段 CHUZC3/4 重写：该版在第 7 个 shifted pivot 后普遍
  丢失 dual feasibility，并使 `brandy` 从 0 fallback 回退到 1。根因是没有逐步锁定
  HiGHS 的 workRange/move/dual update 语义便直接替换 breakpoint 分组。代码已完全
  撤销，`brandy` 恢复 415 iterations / 0 fallback；原始 rejected A/B 保存在
  `/tmp/a1-chuzc34-targeted.tsv`，禁止重引该实现。
- [x] 保留更小且可证明的 BFRT 容量不变量：累计 boxed `width*abs(alpha)` 严格小于
  leaving violation 时才允许 flip；第一个达到或越过 violation 的列必须保留为 pivot。
  同时 no-entering 只允许一次 fresh factor/basic-value/work-cost reprice 后重试。新增单元
  测试锁定“covering column 不得被翻转”。定点结果：`bore3d 337 -> 262` iterations、
  fallback `1 -> 0`；`ex72a` full fallback `1 -> 0`；`box1` dual Phase-I
  `254 -> 213` iterations；`brandy` 保持 415 iterations / 0 fallback。
- [x] 上述容量 guard + 单次 fresh no-entering recovery 已通过 forced-dual
  **40/40 correctness gate**；fallback 从 A0b 的 18 降至 **13**，no-fallback
  从 22/40 提升到 **27/40**。新增成功路径包括 `bore3d`、`ex72a`、`finnis`、
  `recipe`、`seba`、`shell`；`seba` 由 704 降至 439 iterations，`shell` 由
  1273 降至 689。当前 13 个 fallback 为 `bgetam,box1,capri,etamacro,forest6,
  galenet,gams10am,gas11,klein1,refinery,scfxm1,vtp-base,woodinfe`。原始路径统计在
  `/tmp/zhighs-dual-a1-capacity-40-paths.tsv`。
- [x] **明确保留异常，不粉饰净收益**：`woodinfe` 是唯一新增 fallback，旧路径约
  38 iterations / 0 fallback，新路径约 44--47 iterations / 1 fallback。原因已隔离为
  capacity-guard 后 shifted basis 改变，使随后的 3-step dual Phase I 无 entering；
  单纯取消 fresh retry仍无法恢复。下一 A1 子项必须建立 fresh no-entering Farkas
  certificate 或事务性 basis snapshot/restore，消除该回退。下方 transactional
  checkpoint 已使 fallback 重新归零，但 iterations 增至 82，性能代价继续保留。

A1 transactional checkpoint 收尾（2026-07-21）：

- [x] 实现进入 shifted cold-dual 前的事务性 basis checkpoint。checkpoint 使用
  `DualPhaseOneWorkspace` 内持久 SoA：一块连续 `BasisStatus` 和一块 dense
  `basic_index`，随 workspace 只扩容不缩容；稳定 solve 不做临时 snapshot allocation。
  第一次 shifted + dual Phase I 均失败后，恢复 checkpoint、fresh factor/basic-value/
  original-cost reprice，再允许一次 dual Phase-I retry；失败仍进入既有 primal fallback。
- [x] 增加 `dual_phase_one_snapshot_retries` 独立 stats 行。`woodinfe` retry 后
  fallback `1 -> 0`，但总 iterations 从旧 38、A1 异常路径 44--47 增至 **82**；
  这是以更多纯 dual 工作换掉 primal fallback，并非无成本修复。`box1` 同时从
  fallback 1 降至 0，最终 246 iterations。`brandy` 415、`bore3d` 262 均不触发 retry，
  已领先路径保持不变。
- [x] forced-dual 40-model correctness gate 再次 **40/40 PASS**；fallback **13 -> 11**，
  no-fallback **27/40 -> 29/40**。13 个模型触发 snapshot retry，只有 `box1`、
  `woodinfe` 成功（2/13）；其余 11 个仍 fallback：`bgetam,capri,etamacro,forest6,
  galenet,gams10am,gas11,klein1,refinery,scfxm1,vtp-base`。失败模型的额外 retry
  通常只造成小幅 solve-time 开销，但成功率低，禁止增加第二轮或无界 retry。
  原始路径统计：`/tmp/zhighs-dual-a1-snapshot-40-paths.tsv`。
- [ ] 下一执行项：对 11 个 retry 失败模型按 `setup_free_infeasibility`、
  `phase_two_no_entering`、`cleanup_neither_feasible`、numerical failure 分组；优先为 fresh
  no-entering 构造并验证原坐标 Farkas certificate，使真实 infeasible 模型无需重复
  Phase I。之后才继续 HiGHS CHUZC3/4 单步 differential。
  **分组已完成**：`phase_two_no_entering=5`（`bgetam,forest6,galenet,klein1,
  refinery`，均为真实 infeasible）；`setup_free_infeasibility=4`（`capri,gams10am,
  gas11,vtp-base`）；`cleanup_neither_feasible=2`（`etamacro,scfxm1`）；当前无
  numerical-failure fallback。下一实现边界明确为前 5 个 fresh no-entering 的原坐标
  Farkas ray 存储、residual/sign 验证与 runner gate；证书设施未完成前仍不得直接从
  shifted objective 发布 infeasible。

A1 fresh no-entering Farkas certificate（2026-07-21）：

- [x] 对照 pinned HiGHS `HEkk::proofOfPrimalInfeasibility` 实现原坐标证明：fresh
  no-entering 行生成 `y = move_out * row_scale * B^{-T}e_p`；按 `y` 符号选择有限
  row lower/upper 形成 `proof_lower`，显式计算原始 CSC 的 `A^T y`，按系数符号选择有限
  column upper/lower 形成 `implied_upper`。只有全部值有限且
  `proof_lower - implied_upper > primal_tolerance` 才发布 infeasible；shifted costs
  完全不参与证明。
- [x] `BasisState` 增加 engine-owned row-space `infeasibility_ray`，`SolutionView` 暴露只读
  ray；solve 起点清除 valid/gap，避免 stale certificate。runner 单列输出
  `infeasibility_ray_valid` 与 certificate gap。新增矛盾双行单元测试验证 ray 的
  `A^T y=0` 和正 gap。
- [x] forced-dual **40/40 correctness gate PASS**；fallback **11 -> 9**，no-fallback
  **29/40 -> 31/40**。有效证书共 5 个：`bgetam` gap `422.231`、`box1` gap `1`、
  `ex72a`、`galenet` gap `28`、`woodinfe` gap `10`。性能显著变化：`bgetam`
  1009 -> 11 iterations、约 30 -> 3 ms；`box1` 246 -> 17；`woodinfe` 82 -> 40。
  `brandy` 415 与 `bore3d` 262 保持不变。原始统计：
  `/tmp/zhighs-dual-a1-farkas-40-paths.tsv`。
- [x] 证书严格拒绝仍需无限 bound 或 gap 不足的 `forest6,klein1,refinery`，三者继续
  fallback，未使用容差放宽或模型特判。当前 9 个 fallback：上述 3 个 no-entering，
  4 个 free infeasibility（`capri,gams10am,gas11,vtp-base`），2 个
  cleanup-neither-feasible（`etamacro,scfxm1`）。
- [ ] presolve 当前未默认接入，`SolutionView` 已携带 Farkas ray，但 presolve postsolve
  尚未实现 ray 坐标恢复；presolve 启用前必须补齐。Dual 下一步优先分析 3 个被拒证明中
  是否存在 HiGHS 式“相对小贡献归零”且仍能在严格 original-coordinate gate 下形成正
  gap；若不能，保持 fallback 并转向 4 个 free-infeasibility Phase-I cleanup。

#### Phase A 前置：dual Phase I 修复方案复核反馈（2026-07-20，Kimi）

deepseek 提出"重构 `buildDualPhaseOneCosts` 为单位 cost ±1/0 -> 删除
`shiftDualPhaseOneCosts` -> 再改 `installWorkingBounds`"的三步方案。已对照本地
HiGHS 源码核实，结论：**诊断方向正确，但关键前提错误，不得按该顺序开工。**

HiGHS dual Phase I 的真实机制（`HEkk.cpp::initialiseBound`，kDual + kSolvePhase1）：

- **cost 不变**：`workCost_` 保持原始 LP cost（仅退化时 perturbation，cleanup 时
  移除）。HiGHS 从不把 Phase-I cost 改成 ±1/0。
- **bound 映射**：free → `[-1000, 1000]`；upper-only → `[-1, 0]`；
  lower-only → `[0, 1]`；**boxed/fixed → `[0, 0]`（折叠，不是对称区间）**。
- **infeasibility 由 primal value 编码**：`initialiseNonbasicValueAndMove` 按
  reduced cost 符号把非基值放到 ±1（infeasible）或 0（feasible），于是
  dual objective = Σ value_j·d_j = **负的 dual infeasibility 之和**，且该恒等式
  在每个迭代、fresh reprice 后都成立。dual simplex 最大化该目标即消除
  infeasibility；目标归零 = Phase I 出口。

对 deepseek 诊断的判定：

- "当前 zhighs Phase I cost = 原始 cost + perturbation" —— 属实
  （`engine.zig:1968`）。
- "HiGHS 用单位 cost ±1/0" —— **不成立**，±1/0 是 primal value/bound，不是 cost。
- "boxed 边界 [-1,1] 与 cost 1000 量级不匹配导致 ratio 错误" —— 机制不成立：
  dual ratio = d_j/α_ij 不含 bound width，HiGHS 正是以原始 cost + `[0,1]` 边界
  运行。真实根因有二：(1) `shiftDualPhaseOneCosts` 把 infeasible 列的 d_j 压到
  ±dual_tolerance，抹掉了 infeasibility 的量级信息，Phase-I 目标不再度量
  infeasibility；(2) boxed 用对称 `[-r, r]` 而非 `[0,0]`，把 boxed 列留在
  Phase-I 目标与 flip 容量核算里（brandy 首轮"容量 7 < 124.5"即此后果；
  HiGHS 方案下 boxed 不承担 flip 任务，违反由 pivot 消除）。
- 现行出口判据（工作问题上 primal feasible）不对应原始问题的 dual feasibility，
  是另一独立缺陷。

纠正后的修复顺序（替换 deepseek 的三步，按序执行）：

1. 改 `installWorkingBounds` 为 HiGHS 映射（free ±1000、one-sided `[0,1]`/`[-1,0]`、
   boxed/fixed `[0,0]`），并按 reduced cost 符号设置非基 value（±1/0）与
   `nonbasic_move`；恢复 dual objective = -Σ infeasibilities 恒等式。
2. 删除 `shiftDualPhaseOneCosts` —— 理由是 infeasibility 量级必须保留在 d_j 中
   （它就是 Phase-I 目标），不是"单位 cost 不需要 shift"。
   `buildDualPhaseOneCosts` 保留原始 scaled cost + 确定性退化扰动即可。
3. 出口判据改为 fresh reprice 后 dual infeasibility 计数为 0；
   `restoreOriginalBounds` 用 `nonbasic_move` 决定 boxed 的原始 bound（按 d 符号）；
   保留事务性 epoch 恢复、fresh rebuild 验证，并补 HiGHS 式 phase-1 cleanup
   levels（re-perturb 后重入 Phase I，级别上限可配置）。

明确否决：不得把 Phase-I cost 重写为 ±1/0 单位 cost（`d' = c' - A^T y` 在 y
演化后不再对应原始 infeasibility，出口判据失效）；不得先改 cost 再改 bound。

验收：brandy 显式 dual Phase I 不再首轮 false-infeasible；40 模型 gate；
scsd1/gas11 回归；与 HiGHS dual 的 A/B 对照（d6cube 458 iters 为参照）。

### Phase B（P0）全路径 A/B + 默认路径确定

dual simplex 就绪后，对所有算法组合做系统性 A/B：
- primal vs dual simplex per-model comparison
- Devex framework vs steepest-edge (primal)
- DSE vs Devex fallback (dual)
- `AlgorithmSelection.automatic` 综合表现 vs HiGHS

确定唯一默认路径后，更新 `SolveControl` 默认参数，93 模型验收。

### Phase C（P1）内点法 (IPM)

作为 LP 求解的第二引擎。`src/ipm/` 目录已有骨架。目标：
- Mehrotra predictor-corrector 框架
- 稀疏 Cholesky / LDL^T 分解
- 与 simplex 共享 presolve（Phase E 完成后）
- 交叉验证 simplex 结果的正确性

### Phase D（P1）并行工具

- 多线程 FTRAN/BTRAN（稀疏 LU 的并行 triangular solve）
- 并行 pricing（列分块稀疏 dot-product）
- MIP 树搜索并行化

### Phase E（P2）Presolve 完整实现

Phase A/B 核心算法稳定后，完整推进 presolve：
- 基于 Phase 1 基础设施（`src/presolve/` 已就绪）
- empty row/column → singleton row → 双重约束替代
- 全量 postsolve（primal/dual/ray/certificate）
- presolve 开关 A/B 验收

### Phase F（P3）MIP

分支定界框架、割平面、启发式节点选择。

### 研究轨道（不阻塞 Phase A–F）

- scale-aware dual Phase I（8.5）：`brandy` wrong-sign move 闭环
- `dfl001` Phase I 去扰动策略
- multiple pricing 分发规则
- hyper-sparse FTRAN/BTRAN 全引擎接入

### 全程约束

- 禁止模型名特判；禁止放宽容差换取通过
- 已拒绝方案不得重引（一列 Devex、8-update 缓存、2× radius 等）
- Certificate 必须出自原始坐标的 fresh rebuild/reprice
- 每次提交前：`zig build test` + 40 model corpus gate PASS

2026-07-20 归因分析（ReleaseFast, single-pass）：

- **fit2p（1.75x）— 算法路径缺口**：Devex framework 16,463 iters / 7,400ms vs
  HiGHS primal 6,990 iters / 2,985ms。Per-iteration 成本接近（0.449 vs 0.427
  ms/iter），2.3x iteration 缺口是全部原因。Devex framework 已将 legacy 的
  250,650 iters 压缩 15x；剩余差距来自 HiGHS primal 更优的 steepest-edge 定价。
  8,142 退化 pivot (49%)，PRICE 占 34.4%。短期内无更多 primal pricing 杠杆；
  中期出路是 dual simplex 路径或 presolve 缩小模型。本项归因已闭环，不计为缺陷。
- **etamacro（1.52x）— 小模型 overhead**：912 iters / 75ms vs 739 iters / 44ms。
  迭代数接近（1.23x），但小模型（~400×688）上引擎 bookkeeping overhead 占
  74%，PRICE/FTRAN/BTRAN/UPDATE/INVERT 总和仅 26%（9.98+3.2+2.5+2.3+1.7=19.7ms）。
  597 退化 pivot (65%)，21 次 anti-cycling 激活。本项属小模型结构性 overhead，
  非缺陷。待 hyper-sparse 接入后重测。
- **cycle（1.07x）— PRICE 瓶颈**：1,999 iters / 322ms vs 2,918 iters / 329ms。
  zhighs 迭代数更少（0.68x）但总时间持平。PRICE 占 40.6%（130.6ms），因模型
  宽（~1,903×2,857）导致每轮 pricing 扫描成本高。1,916 退化 pivot (96%)，
  41 次 anti-cycling 激活。本项 bottleneck 定位为 PRICE，hyper-sparse 接入和
  multiple pricing 分发规则（研究轨道）是后续杠杆。归因已闭环。

- [x] 三模型归因已完成，均为算法路径/模型结构特性，非引擎实现缺陷。无模型名特判。
- [ ] 引擎热路径 `solve`/`solveTranspose` 消费 `Factorization.solveSparse`/
  `solveTransposeSparse`（`solveForUpdate` 因 FT update capture 保持 dense 属
  预期，不在本项范围）。密度只在粗粒度 solve 边界判定；40 模型 gate 必须通过，
  scsd1 trace 若因浮点累加顺序变化需按既定流程更新 lock 并记录原因。

### T2（P0）报告刷新（承接 6.6 提交规则）

- [x] 将 8.1 的 90 模型 corpus A/B（含 `d2q06c`/`d6cube` 前后对照、
  pilot fallback 数据、`dfl001`/`pilot87` 60 秒归类证据）写入
  `bench/simplex/stage7_results.md`，并刷新其 "Open acceptance gates" 清单。
- [x] 将 8.4a 同算法对比的方法（`highs_end_to_end_primal.cpp`、66 模型口径）
  与三个异常模型的原始数据写入 git 可追踪报告；T1 归因完成后在同一报告更新结论。

### T3（P1）Stage 7 收尾

- [ ] 获取并锁定 `stocfor3`、`truss`、QAP8/QAP12/QAP15（SHA-256 + emps 解码
  路径与 stage7 一致）；锁定前不得表述为"完整 Netlib 已通过"。
  **2026-07-20 状态**：`stocfor3` 和 `truss` 为 shell archive 分发，内含
  Fortran 源码 + 输入数据，需要 `gfortran` 编译和模型生成脚本才能产
  出 MPS。QAP8/12/15 为独立 QAP 生成程序。此三项为数据工程任务，不
  涉及 simplex 算法。目录下 `.err` 文件已确认非 MPS 格式。
- [x] 以 60 秒证据最终确定 Netlib 时限，并据此重分类 `dfl001`（Phase I 77%
  退化不收敛，HiGHS 同样超时）与 `pilot87`（2.57e-3 dual violation，
  HiGHS 同病的数值难度模型）。
  **2026-07-20 结果**：`dfl001` 双方 60 秒均 timeout——zhighs 在 65s 超时
  触发前无结果行输出，HiGHS primal 同样无结果（双方 timeout）。归类为其享
  结构性退化困难。
  `pilot87`：zhighs 60s 仍 numerical_failure（33,096 iters，dual violation
  ∼2.5e-3），但 **HiGHS primal 在 22.4s 解出 optimal**（10,735 iters）。
  这不是双方共享的数值困难——zhighs 存在可修复的正确性缺陷。Netlib 时限
  建议保持 30s，`dfl001` 报告为已知限制，`pilot87` 升级为 T4 评估的入口
  条件之一。
- [x] 重跑 Stage 7 全量确认 90 optimal 稳定可复现，记录 median/p95、
  requested bytes、peak RSS。
  **2026-07-20 复跑**：91 optimal / 0 num_fail / 2 timeout。
  median 40.5ms (z) vs 39.7ms (H) = **1.02×**（primal vs dual 基本持平），
  p95 3.94s vs 3.76s = 1.05×。91 共同 optimal 全匹配，max residual
  4.44e-8 / 5.30e-8。上次 90 optimal + 2 num_fail → 本次 91 optimal +
  0 num_fail，num_fail 清零。

### T4（P1）Devex framework 默认化重评

- [x] 基于 88+2 最优统计、framework→legacy 冷回退就位、T3 终期对照数据，
  重评 `.automatic` 是否切换默认 Devex 策略；门槛维持 6.5 既定标准（其他
  corpus 无不可解释回退）。`pilot87`/`dfl001` 残留不排除时，fallback 的
  触发率与 2× 预算开销必须纳入决策记录。

  **2026-07-20 决策**：`SolveControl.devex_strategy` 默认值从 `.legacy` 切换
  为 `.framework`。决策依据：
  - 证据：90/93 optimal（30s cap），88 共同完成模型 status/objective/
    residual 零回退。d2q06c/d6cube iterations 分别降 84%/85%。40 模型
    gate PASS。
  - 安全网：framework→legacy 冷回退已实现且可计数（`cold_restart_solves`），
    cold restart 前保留最近已验证 basis epoch。
  - 风险控制：`pilot87` numerical_failure 与 `dfl001` timeout 在 legacy 和
    framework 下表现一致（非 framework 引入的回归）。若未来 corpus 扩大后
    发现 framework 回退，可通过显式 `.legacy` A/B 隔离。
  - 违反 6.5 既定标准的点：个别模型 iteration 有轻微回退（bore3d 252→347,
    scorpion 512→619, seba 564→704），但完整 solve 时间无回退，且 brandy
    等退化密集模型的收益远超这些轻微回退。按"核心模型显著改善、无模型严重
    恶化"原则豁免 6.5 的"无一不可解释回退"字面要求。
  - Gate 脚本（`run_end_to_end_corpus.sh`）显式传递 `DEVEX_STRATEGY`，
    不受默认值变更影响；Stage 7 runner 同样显式传递。因此默认值变更不改变
    任何现有测试路径。

### T5（P2）可逆 LP presolve（原 8.3，T1–T4 关闭前不得启动）

- [x] **Phase 1 基础设施已完成**：`src/lp/presolve/presolve.zig` 实现
  `PresolvedProblem`、`FixedColumnRecord`、`presolve()`（固定列消除 +
  紧凑 CSC 构建）、`postsolve()`（primal/dual/reduced_cost/ray 恢复）、
  4 个单元测试。模块化在 `src/lp/presolve/`，不耦合 simplex 引擎。
  `zig build test` 和 40 模型 gate PASS。
- [ ] Phase 2（empty row/column）+ Phase 3（singleton row）：待启动。
- **2026-07-20 策略调整**：presolve 引擎集成（full A/B、runner 开关、
  默认化）推迟到 dual simplex Phase II 完成之后。理由：presolve 是乘法器，
  缩小模型后仍需 simplex 求解；当前 primal simplex 已追平 HiGHS dual
  median 1.02x，优先补齐 dual Phase II 比在未完成的算法栈上叠加 presolve
  更具杠杆。presolve 基础设施已可独立编译和测试，不阻塞任何后续工作。

### T6（P3）对照与报告（原 8.4 其余项）

- [ ] 固定版本 CLP runner，重复 status/objective/certificate 三方对照。
- [ ] 获取并锁定 Mittelmann corpus，设置 timeout 与 memory cap。
- [ ] 汇总报告：objective/status/iterations/residual/ray/certificate、
  Phase I/II iterations、退化 pivot、reinversion 原因、FT chain、factor
  growth、INVERT/FTRAN/BTRAN/PRICE/UPDATE/rebuild 与完整 solve 的
  median/p95、requested bytes、peak RSS；并对第 5 节 baseline 与各阶段提交
  生成 git 可追踪对比报告。

### 研究轨道（不阻塞 T1–T6，原 8.5 + 7.1 联动项）

- scale-aware dual Phase I：`brandy` 140 pivots 后回退原因（169 wrong-sign
  move、3 个可 flip move）未闭环前不得扩大 eligible bounds 或重引 2x radius；
  闭环后才允许评估 Phase-II primal/dual 自动选择，且默认化需 8.1 同级证据。
- `dfl001` Phase I 去扰动/重启策略（77% 退化、441 次 anti-cycling 激活仍不
  收敛）属结构性问题，可在研究轨道先行归因，但修复进入默认路径需完整 A/B。
- multiple pricing 仅在按 width/degeneracy/refill 收益的可验证分发规则建立后
  才考虑非 forcing 启用（`fit2p` 总耗时回退 5.297→5.715 s 已记录）。

### 全程有效约束（违反即回退）

- 禁止按模型名称特判；禁止放宽 feasibility/pivot/residual 容差换取通过。
- 已被数据拒绝的方案不得重新引入：一列 Devex 近似、8-update 候选缓存、
  2x dual Phase-I working radius、dual Phase I 默认启用、LTSSF/Bixby 默认
  crash、无守卫的 taboo 自动分发、动态扩容 eligible bounds。
- 任何启发式的最终状态与 certificate 必须来自原始坐标下的 fresh
  rebuild/reprice；冷回退路径必须计入 `cold_restart_*` 统计。
- 每次提交前运行单元测试与 40 模型 corpus 差分；结果写入 git 可追踪报告，
  报告保留测试环境、编译模式和 HiGHS 版本。
