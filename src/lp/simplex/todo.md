# Simplex 实施看板

更新时间：2026-08-21
范围：以 serial revised simplex 为当前交付主线；目标不是逐行翻译 C++，而是在相同模型、
presolve off、串行及相同策略下，得到可验证的状态、目标、证书、残差和稳定性能。

## 1. 当前结论

### 1.1 现在进行到哪一步

```text
模型/矩阵输入
    ↓ 已有
scaling + logical/crash basis
    ↓ 已有
SparseLU INVERT + FTRAN/BTRAN + Forrest–Tomlin UPDATE
    ↓ 已有
primal/dual Phase I → CHUZR/PRICE/CHUZC/BFRT → pivot commit
    ↓ 已有
Phase II → cleanup → optimal/infeasible/unbounded certificate
    ↓ 当前工作区
大语料稳定性 + HiGHS 路径差分 + 稀疏热路径性能
    ↓ 尚未形成产品闭环
presolve/router/public API 接线，以及 IPM/MIP/QP 等 HiGHS 其余求解器族
```

结论：**serial simplex 已经是能真实求解 MPS 模型的实现，不是骨架；当前约处于“算法闭环完成，
正确性和性能收口”阶段。** 不能据此宣称“HiGHS 已重写完成”：multi/task dual、完整 presolve、
顶层算法路由、IPM、PDLP、QP、MIP 和稳定 C ABI 仍未完成。

### 1.2 已有的真实实现

| 能力 | 主要 Zig 文件 | 当前状态 |
|---|---|---|
| basis、缩放、crash、warm start | `basis.zig`, `basis_snapshot.zig`, `engine_setup.zig`, `crash.zig` | 已实现并接入 solve |
| 稀疏数值线性代数 | `factorization.zig`, `simplex_nla.zig`, `matrix/sparse_lu.zig`, `matrix/sparse_ft.zig` | SparseLU + FT 更新已进入主路径 |
| primal revised simplex | `engine_primal.zig`, `pricing.zig`, `ratio_test.zig` | Phase I/II 已实现；大语料仍有历史缺陷待重验 |
| serial dual revised simplex | `engine_dual.zig`, `dual_state.zig`, `dual_phase_one.zig` | Phase I/II、rebuild、cleanup 已形成闭环 |
| dual iteration 数据流 | `dual_iteration_workspace.zig`, `indexed_vector.zig` | `row_ep/row_ap/col_aq/col_bfrt/col_dse` 已单一所有者化 |
| dual pricing 与更新 | `engine_dual_candidate.zig`, `engine_dual_edge_weight.zig` | CHUZR、CHUZC、BFRT、DSE/Devex 已接入 |
| pivot 与最终结果 | `engine_pivot.zig`, `engine_finish.zig`, `solution.zig` | 原子 pivot、残差检查、ray/certificate 已实现 |
| 可观测性 | `engine_progress.zig`, `bench/simplex/end_to_end_runner.zig` | callback、stats、pivot trace 和 corpus runner 可用 |

静态规模仅用于判断实现不是占位：`src/lp/simplex/` 当前约 13.8k 行，含 97 个具名测试；
其中 `engine_dual.zig` 约 2.2k 行，`engine_primal.zig` 约 1.2k 行，factorization 约 0.9k 行。

### 1.3 2026-08-20 实测基线

| 门禁 | 本次结果 | 解释 |
|---|---:|---|
| `zig build test -Doptimize=ReleaseFast` | 166/166 PASS | 当前工作树通过全量单元/模块测试 |
| forced-dual 40-model 状态/目标/残差/ray | 40/40 PASS | 包含 DSE sparse support、两级 SIMD clear、automatic PRICE 与 warm Phase-I/certificate rescue |
| committed pivots | zhighs 8479 / HiGHS 8774 | 约 0.966x；胜/平/负 13/26/1，唯一负例为 `scsd1` 87/86 |
| patch A/B 完整 pivot trace | 40/40 byte-exact | 保存的改动前 runner 与当前 runner；正式策略专属 lock 仍由 P0.2 补齐 |

同机 36-model 性能门禁（31 轮、3 次预热、2026-08-20）在 SIMD clear 前为 zhighs
`227.539 ms`、HiGHS `118.909 ms`，即 **1.914x**；SIMD clear 后、column PRICE 为
`219.369 / 119.727 ms = 1.832x`；当前 automatic PRICE 工作树为
`201.946 / 118.351 ms = 1.706x`；接通 warm Phase-I fallback 后的独立门禁为
`190.500 / 126.986 ms = 1.500x`，6/36 模型达到 `<0.98x`。接着 `klein1` certificate rescue 与
IndexedVector SIMD clear 后，独立门禁为 `165.327 / 118.718 ms = 1.393x`；UPDATE growth support
优化并修复 CPU pinning 后，当前严格门禁为 **`152.334 / 114.807 ms = 1.327x`**，8/36 达标。
该门禁跨轮次受系统频率变化影响，
改动收益以下述同会话、固定 CPU 的 A/B 为准。
当前仍明显没有达到“性能超过 HiGHS”的项目目标。两次 HiGHS 数据都受下述 dirty reference 限制，
仅用于本机开发排序，不是可发布 benchmark。

本次 forced-dual 命令口径：

```sh
ALLOW_HIGHS_VERSION_MISMATCH=1 VERIFY_TRACES=0 \
PHASE_ONE_STRATEGY=dual \
DUAL_EDGE_WEIGHT_STRATEGY=steepest-devex \
DUAL_INITIALIZATION_STRATEGY=highs \
bench/simplex/run_end_to_end_corpus.sh
```

注意：以上对比只可作为当前 Zig 回归结果。本地 HiGHS checkout 是
`dcc25308d890237531422fc4d0e2e8a4740c8608` 且为 dirty，corpus 脚本锁定的却是
`de09bbad9fb7c5d39a1a464a7641bbb5531c6e9d`。在 P0 完成前，不发布新的正式性能结论。

历史 Stage 7（2026-07-20）记录为 zhighs 90 optimal / 2 numerical failure / 1 timeout，
HiGHS 为 93 optimal；其中 `pilot87` 是已确认的 zhighs 特有正确性缺陷。该结果时间较早，必须在
固定参考版本后重跑，不能直接当作当前结果。

## 2. 执行顺序

```text
P0 固定参考和门禁
 ├── P1 收口当前 DSE sparse-support 改动
 └── P2 重跑 93-model，建立当前缺陷清单
          ├── P3 dual 路径差分：gas11 / klein1
          └── P4 primal 数值缺陷：pilot87
                    ↓
          P5 权威 U/UR + FT 内核（由证据决定是否启动）
                    ↓
          P6 simplex 产品化接线与发布门禁
```

规则：任何优化先过策略专属 trace 和 40-model 正确性；任何数值修复再过 93-model；禁止模型名分支、
固定迭代号补丁、为了追平 pivot 数而放宽残差或证书阈值。

## 3. 可执行 TODO

### P0 — 固定唯一 HiGHS 参考与可重复门禁（最高优先级）

- [ ] **P0.1 选择并冻结参考提交。** 优先恢复一个干净的
  `de09bbad9fb7c5d39a1a464a7641bbb5531c6e9d` worktree；若决定升级到 `dcc25308…`，必须在独立
  clean worktree 重建 CMake Release runner，并一次性更新脚本、文档和 lock 元数据。
  - 修改：`bench/simplex/run_end_to_end_corpus.sh`, `end_to_end_*.tsv`, `README.md`。
  - 产物：输出 Zig commit、HiGHS commit、dirty 标记、Zig/C++ 编译参数和模型 hash。
  - 验收：不设置 `ALLOW_HIGHS_VERSION_MISMATCH` 也能运行；dirty reference 必须立即失败。

- [ ] **P0.2 把 trace lock 按策略分开。** 至少建立 `default` 与
  `forced-dual+highs-init+steepest-devex` 两套 lock，lock 头记录全部策略参数。
  - 修改：`run_end_to_end_corpus.sh`, `end_to_end_trace.*.tsv`。
  - 验收：两种策略各自 40/40 PASS；错误使用另一套 lock 时在运行模型前失败，而不是末尾假阴性。

- [ ] **P0.3 生成机器可读摘要。** corpus gate 结束时输出
  `status_pass/objective_pass/max_residual/pivots/win-tie-loss`，并保存 TSV；不再人工从日志相加。
  - 验收：本看板当前的 8479/8774 可由单条命令重现，并与逐模型表求和一致。
  - 已完成子项：`run_per_model_performance_gate.sh` 现真正执行 `CPU=<id>` pinning，并输出 aggregate
    time/ratio、总 pivots 和 pass count；当前 36-model 摘要为
    `152.334 / 114.807 ms = 1.32687x`、`8272 / 8566` pivots、`8/36`。

### P1 — 收口当前 `engine_dual_edge_weight.zig` 未提交改动

当前改动把 DSE recurrence 从整行扫描改成：`col_aq` 为 indexed 时只遍历 producer support，
`col_dse` 的非零统计也优先读 support。它已通过 ReleaseFast 和 40-model 数值门禁，但尚不能直接落地。

- [x] **P1.1 增加 dense/indexed 等价单测。** 同一 `col_aq/col_dse` 数值分别以 `.dense`、自然序
  support、乱序 support 表示，断言 edge weights、`aq/tau` density 统计、fallback 决策完全一致。
  - 文件：`engine_dual_edge_weight.zig`。
  - 必测边界：空 support、仅 leaving row、低于 `highs_tiny`、精确零、support 超过 clear/reIndex
    阈值和非有限更新拒绝。
  - 结果：新增 representation/support-order invariant test；作为独立 test root 执行为 1/1 PASS，
    并在 `simplex/root.zig` 显式纳入 engine fragment 测试发现。

- [x] **P1.2 锁定 producer contract。** 证明 `ftranForUpdate(col_aq)` 在 `.indexed` 状态下的
  `activeIndices()` 覆盖 recurrence 所需的所有可能非零；在 `simplex_nla.zig` 测试 FT correction
  产生的新 support，防止漏更新 edge weight。
  - 结果：新增 64 阶 sparse identity + FT correction 测试，输入 support 仅 row 0，更新后新增的
    row 5 数值也必须出现在 published support；直接执行该文件时该测试 PASS。

- [ ] **P1.3 做 HEAD vs patch A/B。** 同一 clean 基线、同一语料、交替执行；先比逐模型 pivot
  digest，再比 8 个代表模型各 31 轮 `cycles/instructions`，最后比完整 40 模型 solve 中位数之和。
  - 保留条件：40/40、策略专属 trace lock-exact、逐模型 pivots 不意外漂移，并有稳定硬件工作量收益。
  - correctness：40/40 的结果字段与结构 pivot digest 均 0 mismatch，HEAD/patch aggregate 都是
    9360；旧文档的 9365 已确认是陈旧记录。
  - 40% possible-support 阈值：完整 40 模型、每项 31 轮、固定 CPU、交替顺序，中位数之和
    `215.434 → 214.565 ms`（**-0.40%**），29 胜 11 负，最坏 `seba +0.91%`。
  - 负实验：无阈值 sparse loop 为 +1.47%；20% 阈值四模型收益从 -0.61% 降为 -0.29%；10%
    阈值为 +0.004%，均拒绝。保留与 HiGHS `sparseLoopStyle` 一致的 40% 版本。
  - 未完成门禁：宿主 `perf_event_paranoid=4`，普通和受控非 sudo `perf` 均无 CAP_PERFMON；
    `sudo -n` 无缓存凭据。不得把密码写入命令/stdin/脚本，待用户在外部预热 sudo 或调整 sysctl 后
    补测 31 轮 `cycles,instructions,branches,branch-misses,cache-misses`。

- [x] **P1.4 保留清零语义，替换 SparseLU 的标量 `memset` 热点。** `matrix/memory.zig` 已有显式、
  目标相关的 volatile SIMD clear；把 `solveForUpdateSparse`、hyper FTRAN scatter 和
  `solveAdaptive` materialization 的 4 个 dense `f64` clear 改为该 kernel，不改变 support、遍历次序或
  FT 数值路径。
  - 因果验证：`shell` Callgrind 从 `167,236,646` 降到 `146,608,500` instructions
    （**-12.36%**）；此前仅给 runner/core module 设置 `link_libc` 仍调用
    `compiler_rt.memset`，instructions 基本不变，因此链接方案已撤销。
  - correctness：40/40 结果签名一致，40/40 完整 stderr pivot trace 与保存的改动前二进制逐字节一致；
    ReleaseFast 全量测试 PASS。
  - 性能：固定 CPU、3 次预热、交替顺序、40 模型各 31 轮，中位数之和
    `215.568 → 203.111 ms`（**-5.78%**），31 胜 9 负；最差 `klein1 +2.81%`。
    4 个热点模型各 101 轮合计为 `42.0427 → 38.8289 ms`（**-7.51%**）。
  - 后续门禁：取得 perf 权限后补硬件 counters；若不同 CPU/target 上 aggregate 不再为正，允许按
    `target_policy` 选择 SIMD clear 或 `@memset`，不得按模型选择。

- [x] **P1.5 将 PRICE 默认策略从强制 column 改为 automatic。** automatic 按 `row_ep` support 和 CSR
  实际触达量，在一次 PRICE 开始时选择行式或列式 kernel；两种 coefficient inner loop 内均无表示分支。
  HiGHS 同样根据密度在 `priceByColumn` 与 `priceByRowWithSwitch` 间选择，显式 `column` 仍保留作诊断。
  - correctness：165/165 tests PASS；forced-dual 40-model 状态、目标、残差、ray 全 PASS。38 个模型与
    column 的结果/计数相同；`israel 158 → 143`、`lotfi 240 → 238` pivots，目标不变且残差过门禁。
  - 性能：固定 CPU、3 次预热、交替顺序、40 模型各 31 轮，column/auto 中位数之和
    `202.705 → 197.115 ms`（**-2.76%**），21 胜 19 负、最差 `sc105 +1.66%`；5 个热点模型各
    101 轮合计 **-3.15%**，其中 `seba -17.17%`、`shell -3.74%`、`gas11 -1.72%`。
  - Callgrind：`shell` 从 SIMD-column 的 `146,608,500` 降至 `140,210,390` instructions
    （再降 **4.36%**）；row/column dispatch 为 421/149。

- [x] **P1.6 将 IndexedVector 的 dense `f64` clear 接入同一个 SIMD kernel。** `init`、packed values、
  dense `clear` 和 `materializeZeroComplement` 统一经 `clearValues` 分派；非 `f64` 泛型仍使用
  `@memset`，不改变 indexed support 的所有权和发布顺序。
  - correctness：165/165 tests PASS；40/40 结果字段及 40/40 完整 stderr pivot trace 与保存的改动前
    runner 逐字节一致。
  - Callgrind：`shell` 从 `140,209,274` 降至 `128,042,115` instructions（**-8.68%**）。
  - 性能：固定 CPU、3 次预热、交替顺序、40 模型各 31 轮，中位数之和
    `166.976 → 157.672 ms`（**-5.57%**），39 胜 1 负；4 个热点模型各 101 轮合计 **-6.11%**。

- [x] **P1.7 用 `col_aq` producer support 计算 UPDATE growth。** `PivotUpdateView` 可携带 exact
  `direction_indices`；仅在 dimension ≥128 且 support ≤10% 时做 indexed gather，小维度或较稠密列仍走
  可向量化的连续 dense scan。
  - contract：新增 stale NaN/Inf complement 与越界 support 单测；当前 166/166 tests PASS，40/40 结果与
    完整 pivot trace byte-exact。
  - 因果：`shell` Callgrind `128,042,115 → 125,718,085`（**-1.82%**）。
  - 同代码布局严格 A/B：40 模型各 31 轮 `157.233 → 154.808 ms`（**-1.54%**），32 胜 8 负，
    最坏 `sc50b +0.84%`；`shell -2.65%`、`seba -4.09%`。4 模型各 101 轮合计 **-2.61%**。
  - 阈值实验：无 dimension gate 时 `klein1 +4.08%`、`adlittle +2.23%`；加入 dimension ≥128 后分别
    `-0.40%/+0.07%`，因此拒绝无条件 indexed gather。

- [x] **P1.8 让 dual row PRICE 消费 BTRAN producer support。** `computeDualTableauRow` 已在 BTRAN 后
  `reIndex()` 得到自然序 `row_ep` support；新增 support-aware row/column selector，并让 dual row scatter
  只遍历 support。阈值和自然列 support 顺序不变，因此不改变 pivot 决策。
  - correctness：ReleaseFast 全测 PASS；forced-dual 40-model corpus PASS，逐模型 committed pivots 与
    改动前 adaptive runner 一致。
  - A/B：改动前 adaptive runner vs 当前、CPU 2、自动 PRICE、40 模型各 31 轮、3 次预热、交替顺序，
    中位数总和 `148.512 → 143.372 ms`（**-3.46%**），34 胜 6 负；`shell -5.27%`、`seba -8.97%`、
    `gas11 -4.91%`，最坏 `galenet +3.32%`（5 pivots tiny model）。
  - 当前 zhighs/HiGHS 36-model gate（auto PRICE）：aggregate `141.106 / 114.907 ms = 1.228x`，
    **11/36 PASS**，committed pivots `8272 / 8566`；主要剩余热点仍是 `shell/seba/scagr25` 单次迭代。
  - 负实验：row PRICE scatter 期间标记 touched structural columns，并在 touched ≤ 50% 时排序后只 finalize
    touched 列，保持 40-model correctness 和 pivots 不变，但完整 40 模型各 31 轮为
    `143.247 → 159.712 ms`（**+11.49%**），1 胜 39 负；`shell +29.49%`、`scrs8 +25.17%`。
    mark/sort 和额外写流成本明显超过省下的自然顺序 dense scan，因此实现已撤销。
  - 负实验：在 `RatioTest.chooseDualEntering` 中对 ≤8 个 CHUZC candidates 使用手写 insertion sort，并跳过
    0/1 个 BFRT flip 的通用 sort，保持 40-model correctness 和 pivots 不变。4 热点模型各 101 轮为
    `41.466 → 41.391 ms`（**-0.18%**），但完整 40 模型各 31 轮仅
    `140.001 → 139.923 ms`（**-0.056%**），最坏 `scsd1 +1.79%`；收益小于噪声和尾部风险，因此实现已撤销。
  - 负实验：把 `sparse_ft_adaptive.zig` 的 tagged-union thin wrappers 强制 `inline`，保持 40-model
    correctness 和 pivots 不变，但 4 热点模型各 101 轮为 `41.697 → 42.631 ms`（**+2.24%**），
    `seba +4.96%`、`shell +2.74%`。代码膨胀/指令局部性损失超过分发开销，因此实现已撤销。

### P2 — 重建当前 93-model correctness 矩阵

- [ ] **P2.1 用冻结参考重跑 Netlib 93-model。** 分别跑默认策略、forced dual、forced primal，记录
  status、objective、primal/dual residual、certificate、failure site、iterations 和 timeout。
  - 文件：`run_stage7_corpus.sh`, `stage7_netlib.lock.tsv`。
  - 验收：结果写到新的日期化 TSV；`stage7_results.md` 只引用机器产物。
  - 2026-08-21 当前机器只找到 30/93 个锁定 MPS（`scipoptsuite/soplex/check/instances`），hash 全部匹配；
    其余 63 个（含 `pilot87/dfl001/tuff`）本地不可见，完整 93 门禁暂被数据集缺失阻塞。
  - 30-model forced-dual 子集：zhighs **30/30 optimal**，HiGHS **30/30 optimal**；目标 30/30 在
    `1e-7` 相对容差内匹配；committed pivots `7109 / 7241`，最大残差 `1.54e-9 / 3.73e-9`。
  - 同一 30-model 子集用 Stage7 脚本默认 `primal logical legacy` 路径时，zhighs 为
    **5 optimal / 25 numerical_failure**，HiGHS 为 **30 optimal**；这属于 P4 primal 数值收口，不归因于
    adaptive FT forced-dual 路径。

- [ ] **P2.2 将失败分为三类。** `wrong status/certificate`、`numerical_failure`、`timeout`；每个失败
  保存最后 fresh rebuild snapshot 和末 64 个 pivot，不允许只记模型名。

- [ ] **P2.3 建立发布底线。** 先达到 93/93 与 HiGHS 状态一致、共同 optimal objective 在
  `1e-7` 相对容差内、published residual/ray 不超过现有 gate，再做全局性能优化。

### P3 — dual 路径剩余差异：先解决额外 fallback，不追逐任意 pivot 相等

当前 40-model 总量已比 HiGHS 少约 295 pivots；这不是“更正确”或“更快”的充分条件，仍需按模型锁定
状态、证书和生命周期。当前重点回归模型为：

| 模型 | zhighs / HiGHS | 已知现象 | 下一观测点 |
|---|---:|---|---|
| `gas11` | 593 / 699 | 已消除 cold fallback；Zig 现少 106 pivots | 锁定 dual→warm-primal 生命周期与 unbounded ray |
| `klein1` | 134 / 134 | raw ray 仅因绝对尺度触发 infinite-column tail；严格失败后归一化 rescue 已消除 299-pivot cold Phase-I | 锁定 raw-first 证书语义，禁止对已有效 ray 做无条件缩放 |
| `scsd1` | 87 / 86 | 仅差 1 pivot | 只做回归样本，不为它写特例 |
| `etamacro` | 526 / 532 | Zig 已少 6 pivots | 保持为 phase transition/cleanup 正向门禁 |
| `israel` | 143 / 146 | automatic row PRICE 比 column 少 15 pivots，且已少于 HiGHS | 锁作 automatic PRICE 正向门禁 |

- [ ] **P3.1 为 HiGHS runner 加同构 trace。** 每 pivot 输出 phase、enter/leave、row、alpha、theta、
  update count、DSE updated/computed weight、rebuild reason、dual objective；浮点同时输出 hex bits。
- [ ] **P3.2 自动找首次分叉。** 新增脚本对两个 trace 做 longest-common-prefix，并把分叉归属到
  CHUZR、PRICE、CHUZC、FTRAN、UPDATE 或 phase transition；输出可直接贴进 issue 的最小报告。
- [ ] **P3.3 只修首次分叉所属模块。** `gas11` 和 `klein1` 的已知 fallback 已收口；后续新差异仍按
  trace 首次分叉选择模块，不在 cleanup 中增加补偿分支或模型特例。

- [x] **P3.4 接通 dual Phase-I 退出后的 warm primal Phase-I。** 本地 HiGHS 开发日志证明 `gas11`
  在 dual Phase-I iteration 484 清除 cost perturbation 后仍有 2 个 Phase-I dual infeasibilities，随后
  保留 basis 进入 primal cleanup；zhighs 此前在相同诊断（2 个、max `1.069e-5`）处返回
  `not_implemented`，外层丢弃 483 个 pivots 并从 logical basis 冷启动。
  - 修改：`engine_dual.zig` 在恢复原始 bounds、确认 basis 同时 primal/dual infeasible 后调用已有但未
    接线的 `solveWarmPhaseOne`；成功后直接进入 primal Phase II，失败仍保留原 cold fallback。
  - `gas11`：`1157 → 593` pivots（HiGHS 699），warm primal Phase-I 80 pivots；状态仍 unbounded，
    improving ray residual `9.85e-14`。forced-dual 40-model certificate gate PASS，165/165 tests PASS。
  - 严格 A/B：保存的改动前二进制 vs 当前、固定 CPU、3 次预热、交替顺序、40 模型各 31 轮，
    中位数总和 `212.926 → 168.702 ms`（**-20.77%**），39 胜 1 负；`gas11 -69.54%`。
  - 参考一致性：HiGHS 日志为 `DuPh1 485; PrPh1 97; PrPh2 109; PrSwap 88; Total 691`，随后验证
    cleanup 到 699；Zig 不要求逐阶段 pivot 相等，但必须保留证书和 warm-basis 生命周期。

- [x] **P3.5 以 raw-first、normalized-rescue 修正 `klein1` 的 full-dual Farkas 验证。** 临时
  `Highs::getDualRay` 诊断证明 Zig pivotal-row ray 与 HiGHS 几乎逐分量一致（例如 row 0
  `217.052134376 / 217.052134686`，row 50 `669251.273741 / 669251.274701`）；失败来自最大值约
  `6.69e5` 时累计的绝对舍入尾，而不是 ray 来源、符号或约束语义错误。
  - 实现：先用原始 ray 运行严格 validator；仅在失败原因是 `infinite_column_bound` 时重建 pivotal ray，
    且只向下归一化到 max 1 后重试。不会放大小 ray，也不会掩盖其他 certificate failure。
  - 负实验：无条件归一化虽修好 `klein1`，却令 `refinery 361 → 1369` pivots，因此撤销；当前
    raw-first 设计保持 `refinery 361`。
  - 回归：新增 `1e12` scale-invariance validator 用例；165/165 tests 与 40-model certificate gate PASS。
    `klein1 434 → 134`，primal Phase-I `299 → 0`，与 HiGHS 同为 134 pivots。
  - 严格 A/B：40 模型各 31 轮 `169.613 → 166.568 ms`（**-1.79%**），`klein1 -62.91%`。

### P4 — primal simplex 正确性收口

- [ ] **P4.1 首先复现 `pilot87`。** 固定 forced-primal、presolve off、60 s/迭代双限制，锁定
  `direction_column` 或 `optimality_check` 的实际 failure site；若当前版本已通过，则把模型加入回归。
- [ ] **P4.2 保存失败前数值证据。** pivotal column residual、row/column alpha、condition estimate、
  FT update count 和 fresh-factorization 重试结果必须进入 trace。
- [ ] **P4.3 修复后跑三层门禁。** 最小回归模型 → 40-model 三策略 → 93-model forced primal；不能只让
  `pilot87` 返回 optimal 而不验证 KKT。

### P5 — 直接替换 mutable U/UR 权威表示（下一性能主线）

现在已有足够证据启动本阶段：`shell` 的 FT UPDATE self cost 为 12.14M instructions，BTRAN 中
`solveUpperTranspose` 单独约 9.03M inclusive；多个“在 intrusive 表示外围再加 support/heap/marker”的
实验都降低不了真实时间。本地 HiGHS `HFactor::updateFT` 的生产表示不是逐 entry 双向链，而是连续 U
列加带余量的 UR 行 slab。

- [x] **P5.1 写清 HiGHS 对照 contract。** 本地源码位置：`highs/util/HFactor.h` 和
  `HFactor.cpp:updateFT`。
  - U：`u_start/u_last_p/u_index/u_value` 为 append-only 连续列；旧 logical pivot 标 `-1`，新 pivotal
    column 追加到尾部，`u_pivot_lookup/index/value` 保持 stable row ↔ logical position。
  - UR：`ur_start/ur_lastp/ur_space/ur_index/ur_value` 为每行连续 slab；空间用尽时按
    `row_count * 1.1 + 5` 搬到尾部。删除 reciprocal entry 用 swap-with-last，不维持 row 内排序。
  - UPDATE：先从所有 U 列删除 pivotal row，再从所有 UR 行删除 pivotal column；随后按 Aq pack 顺序
    追加新 U 列并写回 UR；FT correction 仍按 Ep pack 顺序追加，pivot 更新为 `old_pivot * alpha`。
  - FTRAN/BTRAN：分别顺序遍历 U column segment / UR row segment；logical holes 跳过，浮点遍历顺序由
    segment 内存顺序决定。

- [x] **P5.2 先实现 test-only packed FT oracle。** `src/matrix/sparse_ft_packed.zig` 已实现
  `reset/prepareFtran/solveUpper/solveUpperTranspose/captureEp/update`，U 为连续列、UR 为带余量行 slab。
  固定 4 次 replacement 和确定性随机 32 次 replacement 均与 intrusive oracle 逐位一致；`validate`
  同时检查 U→UR、UR→U reciprocal entry 唯一性。首次差分曾发现 reset 借用 `work` 计数后未清零，
  已改为复用未发布的 `row_space` 元数据计数并由测试锁定。
  原验收要求：只接 `reset/solveUpper/solveUpperTranspose/update`
  四个接口；用固定和随机 replacement 序列逐次对比当前 intrusive 实现的向量、pivot lookup、U/UR
  reciprocal invariants。测试阶段允许两份表示，production solve 禁止双写。

- [x] **P5.3 切换 production 权威表示。** 纯 packed 版本已实际接入并验证，不是 shadow/double-write；
  FTRAN、BTRAN、`captureEp`、UPDATE 全部只读写 packed streams。严格数据证明小维度 dense FT 不适合该
  表示，因此 production 已改为 tagged-union adaptive backend：`n < 128` 用 linked，`n >= 128` 用 packed，
  任一时刻只存活并更新一个后端，仍禁止每 pivot shadow synchronization。
  - 新文件：`src/matrix/sparse_ft_adaptive.zig`；`SparseLU` 当前源码已导入它，并把测试改为
    `updateCount()` accessor。
  - 验证：`sparse_ft_adaptive.zig` 单测、linked/packed differential、ReleaseFast 全测、runner build
    和 forced-dual 40-model correctness 均通过；40 模型逐模型 committed pivots 与 linked 基线一致。
  - 修正：linked FT 的 producer-order 单测不再使用会被新 pivotal column 精确消掉的全 1 spike，避免
    测试期待不存在的 support。

- [ ] **P5.4 门禁与删除旧路径。** 先要求 randomized differential 及 40/40 byte-exact；再测
  `shell/seba/scagr25` 101 轮和完整 40 模型 31 轮。保留条件：aggregate 至少 -3%、热点三模型无
  >1% 稳定回退、
  `requested_bytes` 不超过当前 1.5x；随后跑 93-model 并删除旧 production 路径。
  - 已过门禁：linked 基线 vs adaptive、CPU 2、交替顺序。4 热点模型各 101 轮为
    `46.362 → 45.366 ms`（**-2.15%**），`shell -2.31%`、`scagr25 -5.86%`、`gas11 -0.39%`、
    `seba +1.34%`。完整 40 模型各 31 轮为 `157.789 → 151.481 ms`（**-4.00%**），32 胜 8 负，
    最坏 `klein1 +3.01%`；迭代数逐模型一致。
  - 内存：40 模型 stats 的 `requested_bytes` 合计 `20.479M → 20.003M`（0.977x），最坏
    `etamacro 1.006x`，低于 1.5x；peak RSS 基本不变。
  - 剩余：补齐缺失 Netlib MPS 后跑完整 93-model correctness；linked 后端作为 small-dimension backend
    仍是 production 的一部分，后续只能删除旧的公开绕行入口，不能删除 adaptive 需要的小维度实现。

#### P5 暂停检查点（2026-08-21）

已完成且有可复现证据：

- 纯 packed production 在每个实验版本上均通过 40/40 结果字段和完整 pivot trace byte-exact；adaptive
  wrapper 加入后已重跑 `sparse_ft_adaptive.zig` 单测、ReleaseFast 全测、runner build 和 40-model
  correctness，均通过。
- 初版 UPDATE 对 `captured_aq` 做三次全维扫描：4 热点模型各 101 轮为 **+3.72%**，`shell` Callgrind
  `125.72M → 130.97M` instructions，否决该实现。
- 把“筛选新列、计算 UR relocation、预写未发布列”融合为一次 dense scan 后，`shell` 降到
  `125.49M` instructions，但 4 模型仍约 **+3.84%**。cache simulation 给出根因：D1 misses
  `2,094,190 → 2,188,081`（+4.48%），I1 misses `50,743 → 59,678`（+17.60%），LLd 基本不变。
- 当 abandoned U entries 超过 live entries 的 50% 时原地压缩连续列，完整 40 模型、31 轮、3 次预热、
  CPU 2、交替顺序为 **`154.447 → 149.904 ms`（-2.94%）**，29 胜 11 负；`agg -9.42%`、
  `brandy -9.51%`、`israel -8.37%`，但最坏 `klein1 +5.39%`，尚未满足 P5.4 尾部门禁。
- 阈值负实验：25% 在 4 模型 101 轮仅 -0.05%，100% 为 +0.05%；50% 为 -0.72%，因此保留 50%。
  给压缩加 dimension ≥128 gate 和 `noinline` 后，4 模型仍为 -0.77%，`klein1 +4.72%`；组件计时显示
  packed 使该 55 阶模型的 FTRAN 约慢 13%，故下一步采用维度自适应而非继续调压缩阈值。
- 尚未跑：完整 93-model correctness（当前缺 63 个锁定 MPS），以及当前 zhighs/HiGHS 36-model 比率。
  adaptive 40-model A/B、30-model forced-dual Stage7 子集和 `requested_bytes` 1.5x 门禁已通过，纯 packed
  的 -2.94% 不再作为保留结论。

恢复工作的第一组命令：

```sh
zig fmt --check src/matrix/sparse_ft_adaptive.zig src/matrix/sparse_ft_packed.zig src/matrix/sparse_lu.zig
zig test src/matrix/sparse_ft_adaptive.zig -OReleaseSafe
zig build test -Doptimize=ReleaseFast
zig build build-simplex-end-to-end -Doptimize=ReleaseFast
bash /tmp/zhighs-ft-correctness.sh
CPU=2 RUNS=101 WARMUPS=5 bash /tmp/zhighs-linklibc-ab.sh
CPU=2 RUNS=31 WARMUPS=3 bash /tmp/zhighs-full-ab.sh
```

其中 `/tmp/zhighs-btran-support-baseline` 是 retained linked-FT 基线；临时脚本和二进制可能在重启后消失，
若不存在，应先从当前记录重建等价 baseline，不能拿 stale `zig-out` 代替。

2026-08-20 SIMD clear 前的 instrumentation 画像（31 轮中位数，stats 运行本身含计时开销）：

| model | total ms | FTRAN | BTRAN | PRICE | UPDATE | rebuild |
|---|---:|---:|---:|---:|---:|---:|
| `shell` | 15.475 | 4.016 | 2.114 | 2.275 | 1.588 | 0.395 |
| `seba` | 8.577 | 1.776 | 1.145 | 1.715 | 1.012 | 0.382 |
| `scagr25` | 16.734 | 4.960 | 3.349 | 0.854 | 1.778 | 1.307 |
| `gas11` | 55.394 | 8.434 | 8.415 | 13.659 | 2.805 | 17.008 |

SIMD clear + automatic PRICE + IndexedVector clear 后，`shell` 整体 Callgrind 为 128.04M instructions；
当前最大 self hotspot 是 FT UPDATE 路径的 `applyForrestTomlinUpdate`，约 12.14M（9.48%）。行级审计
显示 `captureEp` 的逻辑后缀三角遍历是必要工作，但它把精确零也记录进 `ep_touched`，导致后续 clear
和 correction support 二次遍历。

2026-08-20 已完成两个 sparse-`ep_touched` 负实验，均不保留：分支版把 `shell` Callgrind
`128.04M → 125.07M`（-2.32%），但 40 模型 31 轮实际时间 **+0.94%**（8 胜/32 负）；无分支计数版
Callgrind 为 `125.70M`（-1.83%），实际时间 **+1.71%**（3 胜/37 负）。这说明当前 CPU 上新增的
逐逻辑 pivot 数据依赖/比较比省掉后续 support 访问更贵。后续 FT 优化应优先减少三角遍历本身，或让
FTRAN producer 直接发布 `captured_aq` support；不能仅凭 Callgrind 指令下降保留改动。

同日也否决了 sparse `captured_aq` 镜像：直接使用 factor-order support 会令 `agg` trace 分叉；升序排序
虽恢复 40/40 byte-exact，却使 `shell` Callgrind **+0.13%**、4 模型 101 轮 **+6.25%**。排序/marker
成本超过省下的 memcpy 与 dense scan，因此实现已完整撤销。

另外三项布局/遍历负实验也已撤销：BTRAN output support + pivot-order sort 虽 byte-exact、Callgrind
`-0.54%`，但 4 模型 **+2.31%**（`shell +3.94%`）；max-heap sparse upper byte-exact，却令 Callgrind
`+1.55%`、4 模型 **+10.47%**（`seba +25.42%`）；新列优先 high-water 连续分配仅 **+0.11%** 且增加
链内存占用。这些结果共同否定“给 intrusive pool 继续叠索引”的方向。

### P6 — simplex 产品闭环

- [ ] **P6.1 接通 solver → presolve → simplex → postsolve。** 当前 `solver.solveLp` 直接调用 simplex；
  `presolve.zig` 已有 fixed-column reduction，但尚未进入统一 solve pipeline。
- [ ] **P6.2 公开稳定结果语义。** 对 optimal/infeasible/unbounded/limit/numerical failure，统一 basis、
  primal/dual、reduced cost、ray、iteration 与 callback 生命周期。
- [ ] **P6.3 补齐 basis/solution I/O 和 warm reoptimization corpus。** 除一次 cold solve 外，验证导出 basis、
  导入 basis、修改 bounds/cost 后 reoptimize 与 HiGHS 状态/目标一致。
- [ ] **P6.4 明确非目标。** serial simplex 发布不等待 multi/task dual、IPM、PDLP、QP、MIP；这些进入
  项目级 roadmap，不能继续塞进 simplex engine。

## 4. 每次提交的 Definition of Done

```sh
zig fmt --check src build.zig
zig build test -Doptimize=ReleaseFast
git diff --check
```

算法或数值路径改动还必须满足：

1. 最小语义测试覆盖本次 contract，而不是只靠 corpus 碰巧经过；
2. default 与 forced-dual 各使用自己的 trace lock；
3. 40-model 状态、目标、残差和 ray 全通过；
4. 可能影响鲁棒性的改动再过 93-model；
5. 性能提交记录同会话交替 A/B 的 commit、命令、样本数、median、cycles、instructions；
6. 负实验删除实现，只在简短实验日志中保留假设、数据和拒绝原因。

## 5. 下一次动手顺序

- [ ] 先做 P0.1/P0.2，消除参考版本和 trace 策略混用。
- [ ] 取得 perf 权限后补完 P1.3 硬件 counters；DSE 40% 阈值、SparseLU SIMD clear 与 automatic
  PRICE、IndexedVector SIMD clear 暂保留；不使用或记录用户密码。
- [ ] 跑 P2.1 得到当前 93-model 真相；若仍有错误状态/数值失败，暂停全局 perf。
- [ ] correctness 达标后推进 P4 的 primal `pilot87`；P3 的 `gas11`/`klein1` fallback 已收口。
- [ ] 从“P5 暂停检查点”恢复：先让 `sparse_ft_adaptive.zig` 编译并通过边界切换测试，再跑全测和
  40/40 byte-exact；失败时先修 wrapper，不要用磁盘上 stale 的纯 packed runner 做结论。
- [ ] adaptive 正确后先跑 4 模型 101 轮；只有 `klein1/adlittle` 回归显著收窄且 aggregate 为正，才跑
  完整 40 模型和 HiGHS gate。若 tagged-union dispatch 抵消收益，则撤销 adaptive wrapper并重新设计
  backend selection，不保留未过门禁的 production 切换。

完成上述五步后的可验收结果不是“又写了一批代码”，而是：参考可复现、两套策略门禁无假阴性、
当前 DSE 优化有明确结论、93-model 缺陷清单可归零，并且每个剩余迭代差距都能定位到首次分叉模块。
