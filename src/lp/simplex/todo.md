# Dual simplex 对齐路线

基准源码：本地 HiGHS `de09bbad9fb7c5d39a1a464a7641bbb5531c6e9d`。

目标：先让 Zig 的 serial/plain dual simplex 在状态机、工作数组、阶段切换、选主元、更新和
cleanup 上与 `HEkkDual` 源码同构；验证口径为相同模型、presolve off、相同 simplex 选项下的
committed simplex iterations。稀疏矩阵、primal simplex 和 tasks/multi dual 不在本路线首轮范围。

## 0. 框架收敛原则（2026-08-11 重置）

当前目标从局部性能追赶切回框架完整性：算法、状态转移和数值顺序以本地 HiGHS 为准，Zig 只在
表达所有权、泛型专门化、错误联合、内存布局和编译期消除分支上发挥语言优势。完成结构等价前，
不再接受通过模型尺寸、密度区间或固定迭代阈值选择算法行为的新分支。

实施边界：

1. 建立泛型 `IndexedVector(Real)`，覆盖 `HVectorBase<T>` 的 indexed/dense 状态、稀疏清零、
   `tight/reIndex/pack/norm2/saxpy` 和 synthetic tick；先作为独立基础设施，不接入热路径。
2. 建立 `SimplexNla` 边界，让 collect/FTRAN/BTRAN/UPDATE 原生接收同一 indexed vector，索引由
   producer 维护，不允许 solve 后补扫 dense 数组伪造 support。
3. 建立 serial-dual iteration workspace，对应 `row_ep/row_ap/col_aq/col_BFRT`、`dualRow` 和
   `dualRHS`；CHUZR 到 UPDATE 的整个迭代共享它，不再把临时向量散落在 `BasisState`。
4. 用单一 major phase/rebuild controller 接管 B--F 各阶段；新路径通过逐字段/逐 pivot 门禁后，
   删除旧的模型形状开关和 fallback 拼接，而不是继续为它们增加补偿条件。
5. 框架验收顺序固定为：单元语义 → 首次 rebuild 全状态 → pivot 前缀 → 完整模型 → 40 模型。
   在完整 40 模型与 HiGHS 算法路径对齐以前，运行时间只作回归保护，不作设计选择依据。

### F0a — Zig 泛型 indexed vector（完成）

直接参照本地 `HVectorBase.{h,cpp}` 实现拥有内存的 `IndexedVector(Real)`。用显式 enum 表示
HiGHS 的 `count < 0` dense/unknown sentinel，以类型化状态取代带符号魔数；用 comptime `Real`
生成 f64 和扩展精度版本；串行框架不携带 PAMI `next` 链。其余清零阈值、reindex 阈值、索引顺序、
可能非零语义和算术顺序必须与 HiGHS 一致。验收为 ReleaseFast 全测、独立语义测试和零现有求解
行为变化；接入 iteration workspace/NLA 属于后续 F0b/F0c，不在本提交混合完成。

实现已落在 `indexed_vector.zig`：显式 `.indexed/.dense` 状态替代负 count，f64/f128 由同一
comptime 类型生成；保留 HiGHS 的 30% clear、10% reindex、possible-nonzero、`1e-14/1e-50`、
pack snapshot、扩展精度运算顺序、synthetic tick 及 UPDATE scratch 所有权。独立测试覆盖 sparse/
dense clear、自然序 reindex、tight、pack、norm、跨精度 copy 和 saxpy cancellation；ReleaseFast
全量测试及 `git diff --check` 通过，现有 solver 尚未接入，因此求解行为为零变化。

### F0b — serial dual iteration workspace（完成）

参照 `HEkkDual` 的成员生命周期，用 `DualIterationWorkspace(Real)` 唯一拥有 `row_ep`、
`row_ap`、`col_aq`、`col_BFRT` 和 CHUZC 的组合 pack。行向量维度固定为 `num_row`，结构 tableau
维度为 `num_col`，组合 pack 容量为 `num_col + num_row`；`beginIteration` 只按各自产生者 support
清零，不触碰跨迭代 factor 状态。扩展精度由同一 comptime 定义生成，不复制 workspace 代码。

本阶段只建立所有权、维度和 pack 顺序语义并完成独立测试，不接入现有 dual kernel；F0c 再让
SimplexNla 原生消费这些 vector。验收仍是 ReleaseFast 全测和零现有求解行为变化。

实现已落在 `dual_iteration_workspace.zig`：四个 vector 的维度和生命周期互不混用，CHUZC pack
严格按 `row_ap.index` 后接 `row_ep.index + num_col` 的 producer 顺序形成，f64/f128 共用同一
定义。独立测试覆盖非自然 support 顺序、logical offset 和逐迭代 sparse clear；ReleaseFast 全测
通过，尚未接入现有 engine，求解行为零变化。下一步 F0c 建立 NLA 的 vector-native contract。

### F0c — vector-native SimplexNla contract（完成）

建立 concrete `SimplexNla`（无 vtable）并让其拥有 `Factorization`。FTRAN/BTRAN/UPDATE 的公共
参数改为 `IndexedVector(f64)`：空 indexed RHS 保持空；SparseLU hyper kernel 已发布的 support
原样保留；dense backend 或尚未支持 producer index 的 FT 路径必须显式把结果标成 `.dense`，禁止
在 NLA 内全量扫描补造索引。`ftranForUpdate` 保持 captured-Aq 生命周期，UPDATE 必须消费同一
vector。先用 identity 和非平凡 dense basis 验证数值、状态降级和零 RHS；现有 engine 接入留到
F0d，因此本阶段仍要求求解行为零变化。

实现已落在 `simplex_nla.zig`：concrete NLA 拥有 factorization，空 indexed RHS 零操作，SparseLU
hyper FTRAN/BTRAN 的 producer support 保持 indexed，dense/FT 未发布 support 时显式降级 dense。
同时补齐 indexed RHS 的 stale-complement 边界：只有真正 hyper-sparse dispatch 可以直接消费
indices，其余 dense consumer 先用 vector 自有 packed scratch allocation-free 物化零补集。identity、
非平凡 dense basis、stale complement、sparse support 和 captured-Aq UPDATE 测试全部通过；现有
engine 尚未接入，求解行为零变化。

### F0d — engine NLA ownership migration（完成）

将 `SimplexEngine.factorization` 替换为 `SimplexEngine.nla: SimplexNla`，所有 factor build/solve/
update、reinversion policy、统计和测试观测显式经 `nla.factorization` 访问。本阶段是机械所有权迁移：
现有 engine 数值向量仍调用 compatibility slice API，不同时切换 dual iteration 数据流。验收为
ReleaseFast 全测、40 模型 forced-dual correctness 和 committed pivot 数逐项不变。通过后，F0e 才
把 `DualIterationWorkspace` 装入 engine 并迁移一条完整的 BTRAN→PRICE→CHUZC 只读链。

`SimplexEngine` 现只拥有 concrete `SimplexNla`，factor build/solve/update、reinversion policy、
统计、benchmark runner 和 session 内存复用测试均从同一 NLA 所有权边界访问 factorization；旧的
公开 `engine.factorization` 路径已清零。ReleaseFast 全测通过，forced-dual 40 模型 correctness
**40/40 PASS**，committed pivots 合计仍为 **6833**，本阶段未改变任何迭代 kernel 的数值输入或顺序。

### F0e — engine-owned serial dual workspace（完成）

把 `DualIterationWorkspace(f64)` 安装为 `SimplexEngine` 的模型形状资源，生命周期与当前 LP 实例
一致而不再依附 `BasisState`。同尺寸重复 solve 必须复用分配；尺寸变化时以先成功构造、后替换的
transactional 顺序更新，OOM 不得破坏旧 workspace；`deinit/requestedBytes` 必须覆盖完整所有权。

本阶段不读取或写入 iteration vector，求解数值行为应为零变化。验收为生命周期/内存复用单测、
ReleaseFast 全测、forced-dual 40 模型 correctness **40/40** 且 committed pivots 合计仍为 **6833**。
通过后 F0f 才原子迁移 BTRAN→PRICE→CHUZC 的整条只读链，禁止同时保留两套具有不同 support
顺序的 active tableau row。

`SimplexEngine.dual_iteration` 现拥有唯一的模型形状 workspace；同尺寸请求直接保留 vector 和 pack
分配，尺寸变化使用 allocate-before-release，solve/reoptimize、deinit 和 requested memory accounting
均覆盖该生命周期。复用/换形单测及 ReleaseFast 全测通过；forced-dual 40 模型 correctness
**40/40 PASS**，committed pivots 合计仍为 **6833**，确认尚未改变 iteration 数值路径。

### F0f — workspace-native BTRAN→PRICE→CHUZC chain（完成）

一次性迁移完整的 pivotal-row 只读生命周期：`nla.btran(row_ep)` 产生 BTRAN 结果，PRICE 直接写
`row_ap`，workspace 按 structural 后 logical 的顺序建立唯一 CHUZC pack。ratio test 改收
`PivotalRowView`（两段 slice、inline column lookup），不把两段值复制回 `BasisState.tableau`；证书、
Phase-I 诊断、row/column pivot verify 和 dual reduced-cost update 全部读取同一 workspace。

双单纯形不得再把 `BasisState.dual_row` 或 `BasisState.tableau` 当作该迭代的权威值；前者在所有
辅助 DSE 路径迁移后删除，后者仅保留 primal pricing。CHUZC candidate/flip prefix 由 workspace pack
存储，BFRT 不再借用 primal 的 `flip_columns`。验收顺序为 view/pack 单测、ReleaseFast 全测、
`blend/bore3d/etamacro` pivot trace 前缀、forced-dual 40 模型 correctness 及 **6833** aggregate pivots。

`computeDualTableauRow` 现由 `nla.btran(row_ep)` 和直接写入 `row_ap` 的 PRICE 构成，分段
`PivotalRowView` 为 scalar consumer 提供 inline global-column lookup，CHUZC/BFRT 只使用 workspace
pack。证书、Phase-I 诊断、pivot verify、DSE 和 dual reduced-cost update 已迁移至同一权威行；
`BasisState.dual_row/tableau_count` 已删除，`BasisState.tableau/flip_columns` 仅供 primal 使用。

ReleaseFast 全测与完整构建通过，forced-dual correctness **40/40 PASS**、aggregate pivots
**6833**。与 F0e 提交逐模型对比的结构 trace 完全一致：`blend 105` digest `15d11cf0…`、
`bore3d 173` digest `c9af8318…`、`etamacro 484` digest `db4e80c0…`。这确认迁移同时保持了
数值结果、迭代数和 pivot 顺序。

### F0g — single-owner entering column and UPDATE（完成）

将 dual entering direction 从 `BasisState.pivot_direction` 迁移到 workspace `col_aq`：collectAj 直接
发布 RHS support，`nla.ftranForUpdate(col_aq)` 保留 SparseLU 的 partial-Aq capture，row/column
pivot verify、BFRT 后 primal update、DSE/Devex recurrence 和 changed-row gather 均消费该 vector，
最后由 `nla.update(col_aq, ...)` 提交同一个 captured column。负 entering move 仍原地缩放 vector，
保持既有算术顺序和 factor update 的 column-scale 约定。

`BasisState.pivot_direction` 暂只保留给 primal simplex 和 residual scratch；公共 `performPivot` 以
comptime-free algorithm dispatch 绑定一次 direction slice，之后循环不得重复分支。验收为 NLA
captured-Aq 单测、ReleaseFast 全测、F0f 三模型 trace digest 不变、forced-dual correctness 40/40
和 aggregate pivots **6833**。

dual iteration 现由 `fillInternalColumnSparse` 直接发布 `col_aq` RHS support，随后
`nla.ftranForUpdate`、pivot verify、DSE/Devex、basic-value update 和 `nla.update` 始终消费同一
vector；primal direction 继续由 `BasisState` 独立拥有。初版曾在共享 pivot 内依据可变
`engine.algorithm` 选择方向，但 `scsd1` 的 primal cleanup 证明 phase kernel 与顶层算法标签并非
同义（`77→88` pivots、反复 reinvert），因此已改为显式 `performDualPivot/performPivot` 入口，再
汇合到一个 direction-bound transaction；修复后该负结果未保留。

ReleaseFast 全测、runner 完整构建和 forced-dual correctness **40/40 PASS**，aggregate pivots
恢复 **6833**；F0f 三模型的 counts/digests 全部不变（`105/15d11cf0…`、`173/c9af8318…`、
`484/db4e80c0…`）。

### F0h — workspace-native BFRT vector（完成）

将 boxed flips 的 equation-space RHS 和 FTRAN 结果统一放入 `col_bfrt`：CHUZC flip prefix 按顺序
collect 到 vector，producer 用 byte mark 保证行 support 唯一并保留 cancellation，随后
`nla.ftran(col_bfrt)`；basic value correction、changed-row gather 和 infeasibility-list update 直接
消费结果。删除 `residual_work → rhs_work` 的全向量复制，不改变 flip 提交的 transactional 边界。

验收为重复行/抵消 support 单测、ReleaseFast 全测、F0g 三模型 trace digest 不变、forced-dual
correctness 40/40 和 aggregate pivots **6833**。

BFRT 现直接在 `col_bfrt` 中按 CHUZC flip 顺序聚合；byte mark 让共享 matrix row 只发布一次
support，同时允许数值抵消继续保留 possible-nonzero 语义。`nla.ftran`、basic correction、changed-row
gather 和 infeasibility update 消费同一 vector，原先两次 row-sized materialize/copy 已移除。
重复行 support 单测、ReleaseFast 全测和 runner build 通过；forced-dual correctness **40/40 PASS**、
aggregate pivots **6833**，三模型 trace counts/digests 继续保持 `105/15d11cf0…`、
`173/c9af8318…`、`484/db4e80c0…`。

### F0i — complete serial workspace with col_DSE（完成）

修正 F0b 的初始范围遗漏：本地 `HEkkDual::initialiseInstance` 除 `row_ep/row_ap/col_aq/col_BFRT`
外还显式 setup `col_DSE`。在 `DualIterationWorkspace(Real)` 增加同维 vector；DSE recurrence 由
`col_dse.copyFrom(row_ep) → nla.ftran(col_dse)` 产生 tau，禁止继续借用 `BasisState.residual_work`。
copy 必须保留 row_ep producer support 和 synthetic metadata，dense backend 由 NLA 显式降级。

验收为 workspace clear/bytes 与 DSE recurrence 单测、ReleaseFast 全测、F0h 三模型 trace digest
不变、forced-dual correctness 40/40 和 aggregate pivots **6833**。

workspace 现完整拥有本地 `HEkkDual` serial path 的五个 vector；DSE tau 由 indexed `row_ep`
直接 `copyFrom` 到 `col_dse` 后经 `nla.ftran` 产生，recurrence 不再占用 basis residual scratch。
f64/f128 lifecycle/bytes、DSE exact recurrence、ReleaseFast 全测和 runner build 通过；forced-dual
correctness **40/40 PASS**、aggregate pivots **6833**，三模型 trace digests 保持不变。

### F0j — complete HEkkDualRow ownership（完成）

将 CHUZC 的 per-column ratio/group 和 selected move direction 从 `BasisState` 移入
`DualIterationWorkspace`，与 pivotal-row pack/candidate prefix 组成同一 `HEkkDualRow` 对应物。
数组容量固定为 `num_col + num_row`，ratio test 获得 contiguous mutable slices；iteration clear
只重置有效计数，不全量初始化下一轮必然覆写的 work entries。删除 basis 中两个 dual-only allocation。

验收为 workspace shape/bytes 单测、ReleaseFast 全测、三模型 trace digest 不变、forced-dual
correctness 40/40 和 aggregate pivots **6833**。

CHUZC ratio/group 与 move direction 现和 pack/candidate prefix 一起由 iteration workspace 拥有，
`BasisState.dual_ratio/dual_direction` 及其分配、初始化、释放和 memory accounting 已删除。数组保持
SoA 连续布局且只覆写 active candidates。ReleaseFast 全测、runner build、forced-dual **40/40**、
**6833** pivots 及三模型 trace digest 全部保持。

### F0k — complete dual RHS candidate ownership（完成）

将旧 hyper-sparse CHUZR 的 retained candidate row/score 从 `BasisState` 移入持久
`DualWorkState`，使所有 row infeasibility、index/mark、changed-row support 和候选筛选状态均由
同一个 dual RHS 对应物拥有。`ensureCapacity` 必须事务式扩展候选 SoA，`deinit/requestedBytes`
覆盖完整内存；basis 删除 dual-only 字段、分配和初始化。

本阶段只修正所有权，不改变候选容量、评分、遍历顺序或启用条件。验收为 candidate lifecycle/
memory accounting 单测、ReleaseFast 全测、F0j 三模型 trace digest 不变、forced-dual correctness
40/40 和 aggregate pivots **6833**。通过后，major controller 才能在不依赖 basis scratch 的前提下
统一 CHUZR 生命周期，并最终删除旧 selector 与模型形状开关。

候选 row/score 现由 `DualWorkState.ensureCapacity` 与其余 RHS SoA 一次性事务扩展，重复小尺寸请求
保持地址稳定，释放和 requested-memory accounting 均覆盖；`BasisState` 已删除对应 dual-only
字段。ReleaseFast 全测和 runner build 通过，forced-dual correctness **40/40 PASS**、aggregate
pivots **6833**；三模型 trace 保持 `105/15d11cf0…`、`173/c9af8318…`、`484/db4e80c0…`。

### F0l — model-independent HEkkDualRHS lifecycle（完成）

删除 `useHighsRhsInfeasibilityList/useHighsBreakpointGroups` 对 CHUZR 状态生命周期的模型行列范围
门控。只要选择 `.highs` 初始化，cold rebuild、Phase-I bounds 安装、reinversion 和 Phase-II bounds
恢复都必须无条件重建同一 `DualWorkState` infeasibility list，`chooseRow` 也必须始终消费它；模型
密度只能决定 indexed/dense 表示，不能决定是否采用 HiGHS 算法。

本阶段不同时启用 CHUZC4 grouping、改变 update limit 或删除 fallback。验收首先要求 ReleaseFast
全测和 forced-dual correctness 40/40；随后逐模型记录 committed pivot 变化并定位首次分叉。若发生
正确性退化，不恢复模型尺寸门控，而是补齐 `HEkkDualRHS` 缺失的通用 cutoff/rebuild 语义。

`.highs` 路径现在在 cold/phase/reinversion/restore 四类边界无条件创建或重建 RHS list；旧 selector
只作为非 HiGHS 初始化的兼容路径存在。ReleaseFast 全测和 forced-dual correctness **40/40 PASS**。
aggregate pivots `6833→8261`，仍低于 HiGHS 的 `8774`，但个别模型暴露出后续未对齐：`refinery
341→1117`、`blend 105→173`（dual trace）。本地 debug HiGHS 的逐 pivot 检查证明 refinery 新路径
前 7 个 `(entering,row_out)` 与 HiGHS 完全一致，旧路径从第 0 个就分叉；新路径首次分叉在第 8 个：
HiGHS `(9,60)`，Zig `(100,118)`。因此保留通用生命周期，下一步从 incremental primal/DSE update
定位该分叉，不允许恢复尺寸门控。

### F0m — model-independent perturbation and CHUZC grouping（完成）

F0l 后用本地 debug HiGHS 对 `refinery` 逐字段比较：前 5 次 pivot 的 entering/row/alpha/theta 完全
一致；第 6 次仍选相同 `(92,167)`，但 HiGHS `theta=0.8639937870`、BFRT count 0，Zig
`theta=0.1438501116`、flip count 1，随后第 8 次 CHUZR 才分叉。根因边界因此确定为 CHUZC/BFRT，
不是 RHS list 或 FTRAN。

删除 `.highs` 路径 cost perturbation 公式和 CHUZC3/4 breakpoint grouping 的模型形状门控：初始化
策略决定算法，密度只决定 kernel 表示。保留 update-limit、zero-pivot phase shortcut 和 fallback 的
独立实验边界。本阶段验收为 refinery 前 8 次 alpha/theta/flip 与 HiGHS、ReleaseFast 全测、forced-
dual correctness 40/40；若失败，必须修正通用 ratio/update 语义，不得恢复尺寸条件。

`.highs` 现在无条件使用其随机向量 perturbation 公式和 CHUZC breakpoint grouping；baseline 策略
仍保留原比较路径。`refinery` 的前 12 个 `(entering,row_out)` 已与本地 HiGHS 完全一致，第 6 次
alpha/theta 从错误的 `1.60953125/0.14385011 + 1 flip` 修正为 HiGHS 的
`1.60953125/0.86399379 + 0 flip`；首次变量分叉推迟到第 13 次（HiGHS `(107,98)`，Zig
`(50,68)`）。ReleaseFast 全测和 forced-dual correctness **40/40 PASS**。aggregate pivots
`8261→8753`，接近 HiGHS `8774`；该总数不是完成标准，`refinery 1609 vs 396` 和 `ex72a 594 vs
139` 表明后续 incremental update/phase exit 仍须逐字段对齐。

### F0n — HVector-tight FTRAN support before RHS update（否定）

F0m 的第 13 次 refinery CHUZR 中，HiGHS rows 68/98 的 infeasibility 和 DSE weight 各自完全相等，
分叉只来自 retained list/random start。此时 HiGHS `workCount=103`，Zig `infeasibility_count=209`。
源码边界是 HiGHS FTRAN 的 HVector support 已执行 tiny-value 清理，而 Zig 因 FT backend 尚未发布
index，随后以 `value != 0` 把数值填充全部当作 changed rows。

在 `col_aq/col_BFRT` 离开 NLA、进入 primal infeasibility update 前应用泛型 `IndexedVector.tight()`，
再由现有 factor-order gather 产生 changed support；禁止使用模型阈值或 feasibility tolerance 替代
HVector 的固定 `1e-14` 语义。本阶段验收 refinery 首次分叉、RHS list count、40/40 correctness 和
aggregate pivots；下一阶段再让 sparse FT 原生发布 index，删除 reconstruction scan。

隔离结果：在 `col_aq/col_BFRT` FTRAN 后执行 `tight()`，refinery 首次分叉仍为第 13 次，完整模型
`1609→1629`；debug Zig 更证明首次 CHUZR 前 list 已是 `211/dense`，而非迭代更新后逐渐膨胀，
HiGHS 同点为 `30/sparse`。因此 tiny fill 假设错误，实验代码已撤回。差异发生在第一次 FTRAN
之前的 Phase-I bounds/nonbasic-value/rebuild 状态，下一步转入 B2 初始化顺序对齐。

### F0o — HiGHS-default cold basis policy（完成）

继续逐字段检查推翻了“同一初始 basis”的隐含前提：本地 HiGHS `simplex_crash_strategy` 默认值为 0，
refinery 首次 CHUZR 的 basis head 为 logical（row 1 → column 465）；Zig `.automatic` 通过模型尺寸/
密度规则安装 structural crash（同一 row → column 32）。因此 F0m 的 pivot 前缀只是局部巧合，不能
作为完整算法等价证据。

将 Zig `.automatic` crash 解析为 HiGHS 默认 logical basis；仅 `.ltssf/.bixby` 显式请求可安装
structural crash。删除 automatic 的模型形状选择、覆盖率/condition 验收、Bixby 二次尝试、DSE
预算补偿和 recursive logical fallback；显式 crash 的失败仍原子恢复 logical basis。验收为首次
basis head/RHS list 与 HiGHS、ReleaseFast 全测、forced-dual 40/40，以及逐模型 pivot 变化。此项
优先级高于保留既有低迭代数，因为目标是先建立相同算法起点。

`.automatic` 现在严格解析为 logical basis，模型尺寸/密度 crash 选择、coverage/condition 验收、
Bixby 二次尝试、crash 相关 DSE 补偿及递归 logical fallback 共删除 242 行；`.ltssf/.bixby` 仍是
独立显式策略。ReleaseFast 全测和 forced-dual correctness **40/40 PASS**。统一逻辑基后 aggregate
pivots 为 `10702 vs HiGHS 8774`，真实暴露出 `refinery 1377 vs 396`、`brandy 813 vs 304`、
`ex72a 516 vs 139` 等框架差异，不再由数据集 crash 掩盖。

`refinery` 的前 12 个 pivot 仍与 HiGHS 一致；第 13 次 CHUZR 两侧 RNG state、random start、list
count、rows 68/98 的 basic value、squared infeasibility 和 DSE weight 均逐位一致。首次分叉只来自
FTRAN HVector support 的插入顺序：HiGHS 在相同 103-row list 中先遇 row 98，Zig 当前 sparse-FT
reconstruction 先遇 row 68。下一阶段应让 NLA/FT 原生发布 producer-order support，不能再用基维度
或模型名称规则修补 tie-break。

### F0p — FT-native FTRAN producer support（完成）

F0o 首次分叉的全部数值状态和 RNG 已排除，只剩 `col_aq.index` 的 producer order。HiGHS
`HFactor::ftranU` 在 mutable U 上按 reverse logical pivot order 求解并直接写 HVector index；Zig
`SparseForrestTomlin::solveUpper` 只写 dense values，engine 随后错误地用 base LU `pivot_columns`
重建索引。FT 更新追加的 pivotal rows 因而丢失其真实顺序，CHUZR 在精确 merit tie 时分叉。

让 `SparseForrestTomlin` 在 upper solve 内按实际遍历顺序发布 support，`SparseLU/Factorization/NLA`
逐层借用该 slice，并让 `IndexedVector` 保持 `.indexed`；engine 只消费 producer-owned support，不再
推测 factor order。使用固定 `1e-14` HVector tiny 语义在 producer 边界清零。验收为 sparse-FT 单测、
refinery 第 13 次 pivot、ReleaseFast 全测、forced-dual 40/40 和 aggregate pivots；禁止在 CHUZR
中加入人为 row tie-break。

`SparseForrestTomlin::solveUpper` 现在按 mutable-U reverse logical traversal 直接发布经过 `1e-14`
清理的 support；`SparseLU → Factorization → SimplexNla` 保持借用链，`col_aq` 离开 NLA 时仍是
indexed vector，engine 删除 base-LU `pivot_columns` 重建。新增单测覆盖 FT pivot 追加后的顺序和
tiny 清零；ReleaseFast 全测、forced-dual correctness **40/40 PASS**。aggregate pivots
`10702→10571`（HiGHS `8774`），其中 `brandy 813→767`、`vtp-base 175→167`、`lotfi 337→263`；
`refinery 1377→1391`，第 13 次 tie 仍分叉，因为 Zig 与 HiGHS 的 base LU pivot permutation 不同，
而非 support 所有权仍错误。后续不应追求 factor-order 逐项相同，应转向影响 Phase-I exit 的状态机
和 update 语义；该 producer-owned support 作为通用架构修正保留。

### F0q — model-independent HiGHS initialization（完成）

`.automatic` 当前仍由 `shouldUseBaselineDualInitialization` 按 row/column 数、equality 比例、cost 和
matrix range 为 4 个模型选择旧 baseline。显式 `highs` A/B 已通过 correctness **40/40**：
`brandy 767→312`（HiGHS 304），同时 `gas11 595→1075`、`klein1 119→436`、`agg 198→225`，
aggregate `10571→10940`。这证明 baseline 门控既掩盖了成功，也掩盖了待修的统一生命周期缺陷。

删除 `shouldUseBaselineDualInitialization`，令 `.automatic` 无条件选择 `.highs`；`.baseline` 仅保留为
显式诊断策略。验收为 ReleaseFast 全测、forced-dual 40/40，并以显式 `highs` A/B 数值作为路径锁。
本阶段允许 aggregate 暂时回退，后续必须在统一 HiGHS 路径修复 gas11/klein1，不得重新增加模型判定。

`shouldUseBaselineDualInitialization` 及其整套尺寸/scaling 扫描已删除，`.automatic` 现在恒等于
`.highs`，`.baseline` 只能显式请求。ReleaseFast 全测和 forced-dual correctness **40/40 PASS**；
自动路径逐模型与预先运行的显式 `highs` A/B 完全一致，aggregate `10940 vs 8774`。因此后续所有
修复都在同一初始化状态上进行，brandy 的接近对齐和 gas11/klein1 的回退都不会再被策略分流隐藏。

### F0r — source-equivalent dual phase entry（完成）

HiGHS 在 perturb 后若无原始 dual infeasibility 会直接进入 Phase II；否则进入 Phase I，由统一 major
loop 和 rebuild reason 决定退出。Zig 仍仅对两个模型形状开放 zero-pivot Phase-II，按第三种形状缩短
update chain，并用 `shouldPreflightDualPhaseOneFallback` 在 row/flip/mass 阈值上提前放弃 dual Phase I；
零 pivot 后还对一个形状强制 Devex。这些都不是 HiGHS 算法状态。

删除 `useHighsRhsInfeasibilityList/useHighsUltrawideLifecycle/useHighsBreakpointGroups/`
`useShortDualUpdateChain` 和 preflight helper：zero-infeasibility transition 对 `.highs` 路径普遍生效，
其余模型一律进入同一个 dual Phase-I loop；update limit 和 DSE rule 只由通用控制/观测量决定。
验收为 ReleaseFast 全测、forced-dual 40/40、逐模型 phase iteration 与 aggregate；失败只能修统一
phase/rebuild 语义，不得恢复形状阈值。

上述 5 个 helper、preflight fallback、shape-specific update-limit 和 zero-pivot Devex override 已删除；
无原始 dual infeasibility 现在统一进入 perturbed Phase II。ReleaseFast 全测和 forced-dual correctness
**40/40 PASS**。aggregate `10940→11936`，但三组强证据说明统一入口正确：`ex72a 516→136`
（HiGHS 139）、`klein1 436→134`（HiGHS 134）、`etamacro 557→532`（HiGHS 532）。回退集中为
`shell 568→1367`、`scorpion 335→833`、`brandy 312→695`，共同指向 zero-pivot Phase-II 入口缺少
HiGHS `computeDual → correctDualInfeasibilities → computePrimal → rebuild RHS` 顺序，而不是 phase
选择错误。下一阶段在该统一入口补齐 rebuild，不回滚模型分流。

### F0s — zero-pivot Phase-II rebuild（完成）

源码核对确认 F0r 的 regression 边界：HiGHS 决定 Phase II 后，`solvePhase2` 首先执行完整 `rebuild`，
顺序为 perturbed `computeDual → correctDualInfeasibilities(force_phase2) → computePrimal → RHS list`。
Zig zero-pivot 分支仅 `createPrimalInfeasibilityList` 后直接 iterate，漏掉 fixed/boxed flip 和 one-sided
cost shift，导致 Phase II 虽选择正确却不是 dual-feasible work state。

zero-pivot 分支改用现有 `rebuildDualState(use_work_costs=true, correct_dual=true)`，禁止复制另一套 correction；
只有 correction 后 free infeasibility 为零且 basis dual feasible 才进入 Phase II，否则沿统一失败边界返回。
验收重点为 shell/scorpion/brandy，另跑 ReleaseFast 全测与 forced-dual 40/40。

zero-pivot 分支现复用 `rebuildDualState`，并在 free infeasibility/dual feasibility 检查后进入 Phase II。
ReleaseFast 全测和 forced-dual correctness **40/40 PASS**，aggregate 及 shell/scorpion/brandy pivot
数均不变。这一否定结果说明这些逻辑基在 perturbed cost 下无需 correction，缺失 rebuild 不是长路径
根因；调用仍保留，因为它补齐了 source-equivalent phase boundary 和未来非逻辑 warm basis 的正确性。
下一步比较 Phase-II DSE recurrence/切换，而不是继续修改初始化。

### F0t — model-independent DSE lifecycle（完成）

当前 DSE 有两层数据集补丁：recurrence 仅在 65–95 row/超宽形状使用 HiGHS 的 expanded evaluation
order；automatic budget 又按 10 组 row/column/flip/crash 阈值设为 0/1/4/8/16/20/32。HiGHS 对所有
模型使用同一 recurrence，并仅由 costly-DSE density/frequency 或 weight error 切换 Devex。
显式 budget=0 A/B 已 40/40，aggregate `11936→12044`，其中 `refinery 1391→1091`、
`etamacro 532→514`，但 `grow7 387→851`，表明 grow7 依赖人为 20-pivot cap，必须由通用 costly-DSE
观测修复。

删除全部 initial-flip/shape budget 改写；automatic 只运行已有 HiGHS density/error criterion，显式数字
budget 保留为诊断控制。DSE recurrence 对所有模型固定使用 HiGHS expanded expression，禁止按尺寸选择
浮点表达式。验收为 ReleaseFast、40/40、DSE switch 统计和 aggregate；grow7 回退作为下一步观测量
对齐目标，不恢复模型 cap。

所有 shape/flip budget 分支已删除；automatic 不再读取 hard cap，显式数字 budget 的既有诊断语义
保留。DSE recurrence 对所有模型固定为 HiGHS expanded evaluation order。ReleaseFast 全测和
forced-dual correctness **40/40 PASS**，aggregate `12021 vs 8774`；相对 F0r，`refinery
1391→1091`、`etamacro 532→514`、`agg 225→205`，主要回退仅 `grow7 387→851`。下一步校正
costly-DSE density 的 HVector count（当前错误地按 simplex `1e-12` 扫 dense values），使 grow7 由
通用观测触发 Devex。

### F0u — HiGHS HVector density streams（完成）

源码核对进一步定位到密度的维度错误，而不仅是 tiny threshold：HiGHS 独立维护
`row_ep.count / num_row`、`col_aq.count / num_row`、`row_DSE.count / num_row` 和
`row_ap.count / num_col`；Zig 却把 `row_ep + row_ap` 合并后除以 `num_row + num_col`，再作为
`row_ap_density` 输入 costly-DSE。该值在长宽不平衡模型上没有相同物理意义，会错误延迟或提前
DSE→Devex。

恢复四条独立的 producer-owned density stream：优先读取 `IndexedVector.activeIndices().len`，仅在 dense
producer 撤销 support 时用 HiGHS `1e-14` 计数；`row_ap` 只除以结构列数，hyper-sparse PRICE 的算法
判据仍保留其现有 tableau support 口径，避免把观测修复混入定价策略。costly-DSE frequency 和 5%/10%
切换条件保持源码常数。验收为 ReleaseFast、40/40、grow7 DSE 统计及全 corpus aggregate；不得恢复任何
模型预算。

四条 density stream 现已按各自 HVector 维度独立更新，并统一使用 `kHighsTiny=1e-14` 的 support
语义；`row_ap` 不再混入逻辑列，cost switch 也只在当前 iteration costly 时触发，匹配
`HEkk::switchToDevex`。ReleaseFast 全测和 forced-dual correctness **40/40 PASS**，aggregate 保持
`12021 vs 8774`，grow7 仍为 851 且其 `row_DSE / max(Ep,Aq,Ap)` 不满足 costly 条件。这一否定结果排除
了 DSE 切换作为 grow7 长路径根因，同时保留了正确的通用控制框架；下一步转向 basis update/reinvert
lifecycle，其 18 次 update-growth reinversion 明显偏离 HiGHS。

### F0v — remove foreign growth-tier reinversion（完成）

HiGHS 的正常 update lifecycle 只有显式 `simplex_update_limit`、solve/pivot numerical trouble 与 factor
synthetic clock；它不按 `max(abs(FTRAN column))/abs(pivot)` 的 `1e2/1e4/1e7/1e10` 四档缩短 update
chain。Zig 这套档位是早期稳定性补丁：grow7 的 `max_update_growth=1.10e7` 将 100-update chain 强制缩到
16，产生 18 次非源码 reinversion，并改变 rebuild 后的定价路径。

删除 growth-tier 对 reinversion decision 的影响，但保留 maximum growth 作为只读诊断；稀疏 FT 依
显式 update limit，dense Eta 依 retained capacity，residual/pivot warning 仍可提前 rebuild。后续单独实现
producer synthetic tick，不能把 growth ratio 冒充 tick。先验收 ReleaseFast、grow7 correctness/pivots/
reinversion，再跑 40-model corpus；若失败只能修 factor update/residual，不能恢复尺寸或模型阈值。

四档 growth policy 已从 decision path 删除，maximum ratio 仅保留为统计；update limit、Eta capacity 和
numerical warning 继续生效。ReleaseFast 全测和 forced-dual correctness **40/40 PASS**。grow7
`851→482`（HiGHS 339），factor reinversion `20→5`，其中 2 次 update-limit、3 次 solve-residual；全 corpus
aggregate `12021→11725`（HiGHS 8774）。另有 `scorpion 833→824`、`scfxm1 482→480`，主要回退
`israel 139→155`。结果确认外加 growth 档位同时增加 rebuild 和 pivot，不应保留；下一步对齐 HiGHS
synthetic tick 与 numerical-trouble 的 producer 边界，并单独追踪 israel 的 residual-trigger divergence。

### F0w — abort pivot on dual numerical trouble（完成）

HiGHS `updateVerify` 比较 `alpha_col`/`alpha_row` 后，一旦 numerical trouble 便设置 rebuild reason；同次
iteration 的 `updateDual/updatePrimal/updatePivots` 全部跳过，outer loop 完整 rebuild 后重新 CHUZR/CHUZC。
Zig 当前却立刻 refactor、重新 FTRAN，然后沿用 rebuild 前选出的 row/tableau/ratio candidate 提交 pivot。
这违反 reason-driven rebuild 边界，也解释了 F0v 后 grow7 的 3 次 solve-residual 会改变长路径。

引入内部 typed direction result `ready/rebuild_required/numerical_failure`。pivot mismatch 只发布 reason，不在
direction producer 内继续旧 transaction；各 dual loop 在自己的 work-cost/bounds 语境完成 rebuild 并
重新选择，warm zero-cost repair 保持其零 reduced-cost 状态。验收 ReleaseFast、grow7/israel 的
solve-residual 路径、40-model corpus；禁止在 mismatch 后复用旧 pivotal row。

dual direction producer 现返回 typed internal status；mismatch 不再自行 refactor/FTRAN 并复用旧 row，而由
Phase I、Phase II 或 zero-cost repair controller 在各自 work state 中 rebuild 后重新选择。ReleaseFast
全测和 forced-dual correctness **40/40 PASS**，aggregate 与逐模型 pivot 保持 `11725 vs 8774`；grow7
的 3 次 solve-residual 均未经过该 dual mismatch 分支，确认它们来自 primal `B*x-a_q` 周期检查。该修复
保留为 source-equivalent transaction boundary；下一步处理 primal 路径中 HiGHS 不存在的每 8 update
全列残差 gate。

### F0x — HiGHS primal pivot verification（完成）

HiGHS primal 不每 8 updates 构造 `B*x-a_q`；它在 CHUZR 确定 row 后执行一次 unit BTRAN + PRICE，并用
`alpha_col` 对比 `alpha_row`，阈值同 dual 为 `1e-7`。mismatch 只请求 rebuild，禁止提交旧 pivot。
Zig 当前的 whole-column `1e-10` gate 既早于 CHUZR，又立即 refactor 后复用旧 entering candidate；grow7
F0v 后的 3 次 solve-residual 正来自这里。

删除 periodic dense residual/refinement state，在现有 primal reduced-cost pivotal-row producer 中加入 O(nnz
entering column) 的 alpha-row 独立校验，并返回 typed `ready/rebuild_required/failure`；Phase I/II controller
各自 exact reprice 后重新选择。终局 residual validation 不变。验收 ReleaseFast、grow7 路径和 40-model
correctness/aggregate；数值 trouble 不得通过调宽模型阈值规避。

periodic whole-column residual/refinement 及其跨 pivot 状态已删除；primal pivotal-row producer 现在先完整
生成 BTRAN/PRICE row，按 HiGHS `1e-7` 比较 `alpha_col/alpha_row`，通过后才更新 reduced costs，失败则
由 Phase I/II controller rebuild、exact reprice 并重新选择。ReleaseFast 全测和 forced-dual correctness
**40/40 PASS**。grow7 `482→486`（HiGHS 339），solve-residual `3→1`、总 reinversion `5→4`；全 corpus
aggregate `11725→11729`（HiGHS 8774）。4-pivot 增量来自唯一真实 mismatch 后重新选择，是正确 transaction
语义；下一步清除 pivot commit 中剩余的模型形状分支，并继续对齐 basis move/status lifecycle。

### F0y — bound-derived leaving status（完成）

HiGHS `updatePivots` 对所有模型先检查 outgoing variable 的 `workLower == workUpper`：相等时 value 取唯一
bound、nonbasic move 为 0；否则才按 `move_out` 选择 lower/upper。Zig 仅对 65–95 row 且超宽模型把相等
bounds 提交为 `.fixed`，其他模型错误保留 `.at_lower/.at_upper`，使后续 pricing/flip 把 fixed variable 当作
可移动 nonbasic。

删除尺寸判断，对任意 finite/infinite 表示下的相等 bound 统一提交 `.fixed`；非 fixed 继续使用 ratio test
给出的 leaving side。复用 `nonbasicStatusForBounds` 的状态语义，不增加 hot-loop 动态分派。验收
ReleaseFast、40-model correctness/aggregate 和 fixed-bound unit coverage。

模型尺寸分支已删除，所有 outgoing fixed variable 统一提交 `.fixed` 和唯一 bound；非 fixed 仍保留 ratio
test side，并增加纯状态单测。ReleaseFast 全测和 forced-dual correctness **40/40 PASS**。aggregate
`11729→9930`（HiGHS 8774，差距 `33.7%→13.2%`）；`shell 1367→558`（HiGHS 623）、
`scorpion 824→322`（精确等于 HiGHS）、`brandy 695→304`（精确等于 HiGHS），另有
`lotfi 254→229`、`etamacro 514→495`。这确认 fixed move 的模型补丁是三条最大长路径的共同根因；下一步
删除 setup scaling 中同形状的 unit-scale override，恢复统一 HiGHS scaling policy。

### F0z — range-derived scaling decision（完成）

HiGHS `scaleLp` 的默认入口先对所有模型计算 matrix range；若全部绝对值位于 `[0.2, 5]`，直接跳过
matrix scaling，否则执行 forced equilibration。Zig 已复现六轮 equilibration，却只对 65–95 row 且
超宽模型执行该 range gate，仍是数据形状补丁而非算法策略。

把原始 matrix min/max audit 提升为所有模型共享的 setup decision，完全删除 ultrawide 条件；range 内统一
保留 unit row/column scale，range 外运行现有 HiGHS 六轮 power-of-two equilibration。该分支只在 solve
setup 执行，不进入 hot loop。验收 ReleaseFast、40-model correctness/aggregate，并逐模型记录路径变化。

matrix range audit 现对所有模型执行，`[0.2,5]` 内跳过 scaling，范围外才运行六轮 equilibration；原
ultrawide 条件及先计算后覆盖的无效 setup 工作均已删除。ReleaseFast 全测和 forced-dual correctness
**40/40 PASS**，aggregate `9930→9926`（HiGHS 8774），`forest6 80→77`、`sc50b 49→48` 均与
HiGHS 完全相同，其余模型稳定。下一步审计剩余按列/行比例选择行为的分支，区分合法 work-estimate 与
模型 patch，并继续删除后者。

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
  2026-08-08 已对齐 `HEkkDual::updateFtranBFRT` 的 residual ownership：BFRT 仅聚合 RHS 并执行
  一次 FTRAN，紧随其后的 entering-column FTRAN 统一负责本迭代的 residual/reinversion guard；删除
  每个 flip batch 额外的完整基矩阵乘积。forced-dual 40-model correctness **40/40 PASS**，默认策略
  40-model + 四条锁定 trace **40/40 PASS**，八个目标模型 committed iterations 完全不变。
  31 次交错 ReleaseFast A/B 中，`etamacro`（229 个 BFRT batch）中位数
  `39.93 → 39.34 ms`（-1.5%），`shell`（29 batch）`43.95 → 43.74 ms`（-0.5%）；
  无 BFRT 的 `brandy/scorpion` 为噪声对照。E1 仍待把 flip value/status 的提交边界与后续
  pivot-column FTRAN/rebuild reason 合并成完整原子 iteration，因此暂不整项勾选。
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
  2026-08-12 源码协议核对（HiGHS `HEkkDual::solvePhase1`、`rebuild`、
  `assessPhase1Optimality`）：CHUZR 返回 no-row 只产生 `possibly_optimal` rebuild reason，不能直接
  退出 Phase I。控制器必须先以当前 Phase-I bounds/cost 完成一次完整 fresh rebuild（必要时
  reinvert、重算 dual/base value/infeasibility list 和 dual objective），然后再次 CHUZR；仅 fresh
  状态仍无 leaving row 才允许移除 cost perturbation 并判断 Phase-I objective。Zig 当前在首次
  no-row 后立即恢复原始 bounds，随后把双/原可行性组合映射成 cold primal fallback，是
  `refinery 1087 vs 396` 和 `gas11 1127 vs 699` 的主要控制器差异。F1a 将先引入独立
  `possibly_optimal` 数值 rebuild reason 和“首次 no-row 请求 rebuild、fresh no-row 才提交”的
  原子边界；定向验收 refinery/gas11/grow7，随后跑 forced-dual 40-model correctness 和 aggregate。
  F1a 实现结果：新增的 reason 已贯通数值 rebuild 分类和 benchmark 诊断；Phase-I no-row 在
  `has_fresh_rebuild=false` 时请求完整 rebuild，只有 fresh no-row 才提交。ReleaseFast 全测通过，
  定向 3/3 和 forced-dual 40-model correctness **40/40 PASS**。refinery 确实新增 1 次
  `possibly_optimal` reinversion，但 committed iterations 保持 `1087`；gas11/grow7 仍为
  `1127/486`，全语料 aggregate 仍为 `9926 vs 8774`。因此旧状态提前提交已被排除，F1 的下一根因
  是 fresh no-row 后缺少 `assessPhase1Optimality` 的 objective/去扰动/重入协议，而不是 rebuild
  本身。
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
