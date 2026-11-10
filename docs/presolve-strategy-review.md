# Presolve 策略核对、增强实现与基准报告

日期：2026-11-08。范围限定在 `src/presolve` 及 LP/MIP 求解入口，不改变 simplex、
branch-and-bound 或 postsolve 的坐标恢复约定。

## 1. 结论摘要

现有 presolve 已覆盖经典 LP presolve 的核心行/列删除规则，但缺少 HiGHS/SCIP
规则集中最有价值的**通用行活动传播（forcing / implied bound）**与**按比例平行行**
两类规则。本次实现：

1. **行活动强制固定**（`implied_fixing`）：对每行计算其余项的活动区间，向外取整后
   求列的可隐含上下界；只有当区间收缩为单点（MIP 含整数取整后收缩）时才固定变量，
   固定后的变量走既有的精确列替换记录，因此 postsolve 的对偶/约化成本恢复逻辑无需
   修改。
2. **精确比例平行行消除**（`parallel_rows`）：按列模式分桶，要求行向量满足
   `A_k = scale * A_i` 且两个边界在比例变换后逐位相等；负比例自动交换上下界。
   被删行的对偶乘子取 0，代表行保留。
3. 在 `hasPresolveReductions` 门控中加入 forcing 潜力与平行行检测，保证新增规则
   在 LP 路径上真正被调用，而不是被门控跳过。
4. 新增可复现的消融基准：`bench/presolve/presolve_bench.zig`、生成器、
   `compare.py`，可比较 `off / legacy / implied / parallel / full` 五种模式。

在受控算例上，新增规则带来的收益是真实且可测量的：

| 算例 | 关键结构 | legacy 总时间 | full 总时间 | 加速 |
|---|---|---:|---:|---:|
| `forcing_lp.lp` | 60 个变量被强制固定 | 717.2 us | 36.9 us | **19.46x** |
| `parallel_lp.lp` | 30 行精确比例副本 | 197.7 us | 116.9 us | **1.69x** |
| `mip_forcing.lp` | 80 个二元变量被强制固定 | 560.3 us | 106.2 us | **5.28x** |

反过来，在没有任何可约简结构的小模型上，`full` 会多出数微秒到数十微秒的传播
开销；`presolve_mip.mps` 上 `full/legacy = 0.75x`。因此不应无条件宣称“presolve
一定更快”，而应按模型结构/规模决定是否启用更强传播。受控算例规模较小，重复运行
的绝对时间波动可达 20%–50%，上表用于验证方向和数量级，不应外推为通用性能结论。

## 2. 文献核对

### 2.1 核对对象

| 文献/实现 | 与本项目相关的规则 |
|---|---|
| Andersen & Andersen, *Presolving in linear programming*, Mathematical Programming 71:221–245, 1995, DOI [10.1007/BF01586000](https://doi.org/10.1007/BF01586000) | 经典 LP presolve：空行/空列、固定变量、singleton 行/列、forcing 行、冗余行、平行行/列、doubleton 等式、dominated columns、bound tightening |
| Achterberg, Bixby, Gu, Rothberg, Weninger, *Presolve Reductions in Mixed Integer Programming*, INFORMS Journal on Computing 32(2):473–506, 2020, DOI [10.1287/ijoc.2018.0857](https://doi.org/10.1287/ijoc.2018.0857) | Gurobi presolve 分类法；MIP 侧的 probing、coefficient strengthening、clique、implied integer 等 |
| Gleixner, Gottwald, Hoen, *PaPILO: A Parallel Presolving Library for Integer and Linear Optimization with Multiprecision Support*, INFORMS Journal on Computing, 2023, DOI [10.1287/ijoc.2022.0171](https://doi.org/10.1287/ijoc.2022.0171) | 并行 presolve 架构、规则可组合性、多精度/精确算术 |
| Savelsbergh, *Preprocessing and probing techniques for mixed integer programming problems*, ORSA Journal on Computing 6(4):445–454, 1994, DOI [10.1287/ijoc.6.4.445](https://doi.org/10.1287/ijoc.6.4.445) | MIP probing 与预处理 |
| HiGHS 源码 `highs/presolve/HPresolve.h`（master，2026-11-08 访问，<https://github.com/ERGO-Code/HiGHS>） | 可执行规则清单：`removeRowSingletons`、`presolveColSingletons`、`removeDoubletonEquations`、`detectDominatedCol`、`detectParallelRowsAndCols`、`strengthenInequalities`、`detectImpliedIntegers`、`sparsify`、`aggregator` 等 |

> 说明：论文给出的“规则类别”与 HiGHS 的“函数名”并非一一对应；上表按语义对齐。
> 本次没有外网下载 MIPLIB 实例，算例全部来自仓库内模型和可复现生成器，避免把
> 网络/机器差异混入结论。

### 2.2 规则覆盖矩阵

| 规则 | 文献/实现 | 本项目状态 |
|---|---|---|
| 空行、空列、固定变量 | Andersen 1995; HiGHS `kEmptyRow`/fixed col | ✅ 已有 |
| singleton 行 | Andersen 1995; HiGHS `removeRowSingletons` | ✅ 已有 |
| 零成本自由 singleton 列 | HiGHS `kFreeSingCol` | ✅ 已有 |
| 行活动冗余行 | HiGHS `kRedundantRow` | ✅ 已有（向外取整） |
| 精确重复行 | HiGHS dependent rows | ✅ 已有 |
| **forcing / 行活动强制固定** | Andersen 1995; HiGHS `kForcingRow` | ✅ **本次新增** |
| **精确比例平行行** | Andersen 1995; HiGHS `detectParallelRowsAndCols` | ✅ **本次新增（仅精确倍数）** |
| 通用 implied bound tightening（发布收紧后的列界） | Andersen 1995; HiGHS implied bounds | ❌ 未实现；需要同时实现对偶 postsolve |
| singleton 列 stuffing | HiGHS `singletonColStuffing` | ❌ 未实现 |
| dominated columns / dual fixing | HiGHS `detectDominatedCol`/`dualFixing` | ❌ 未实现 |
| doubleton equation / aggregator | HiGHS `removeDoubletonEquations`/`aggregator` | ❌ 未实现 |
| 平行列、dependent free columns | HiGHS `detectParallelRowsAndCols`/`removeDependentFreeCols` | ❌ 未实现 |
| coefficient strengthening（MIP） | Achterberg et al. 2020; HiGHS `strengthenInequalities` | ❌ 未实现 |
| probing / clique / implied integer | Achterberg et al. 2020; Savelsbergh 1994; HiGHS `kProbing` | ❌ 未实现（MIP 引擎已有部分节点传播，但不是 presolve 规则） |
| sparsify、scaling、reduction limit、规则日志 | HiGHS `sparsify`、`setRelaxedImpliedBounds`、presolve log | ❌ 未实现 |

因此，**当前策略仍不全面**，但已经从“无 forcing、无平行行”推进到与经典 LP
presolve 主体规则对齐；剩余缺口主要是 MIP 强化、通用 bound tightening 的对偶
postsolve、以及列级消除（doubleton/aggregator/dominated column）。


## 3. 新规则的数学与 postsolve 安全性

### 3.1 行活动强制固定

对行

```text
l <= a_j x_j + S <= u,    S = sum_{k != j} a_k x_k
```

用当前列界得到 `S in [S_l, S_u]`。向外取整后有：

```text
a_j > 0:  (l - S_u) / a_j <= x_j <= (u - S_l) / a_j
a_j < 0:  (u - S_l) / a_j <= x_j <= (l - S_u) / a_j
```

实现细节（`src/presolve/rules/implied_fixing.zig`）：

- 乘、加、除全部向外取整；除法统一用 `nextAfter(q, direction)`，对正负分母都安全。
- 从行总活动区间中排除当前项时，用该项**最小值的上界**和**最大值的下界**，
  即 `total_lower - minUpper`、`total_upper - maxLower`，保证剩余项界仍然严格。
- 连续变量只有区间宽度不超过 `1e-12 * max(1, |l|, |u|)` 才固定；整数/二元变量用
  `ceil(l)`/`floor(u)` 收缩，唯一整数才固定。二元变量先夹到 `[0,1]`。
- 传播最多 8 轮，链式固定（一个变量固定后使相邻行继续收缩）在同一调用内完成。
- **只返回固定值，不返回收紧后的列界**。固定变量在 `presolveImpl` 中复用既有的
  `ColumnSubstitutionRecord`（`.fixed`），目标偏移与行界调整、postsolve 的对偶和
  约化成本计算都不需要新逻辑。这样避免了“收紧后的列界在 postsolve 中被当成原界”
  的对偶互补性问题。
- MIP 模式（`nature`）下，整数舍入只在临时传播界上使用；输出的固定值保证是整数。

正确性边界：如果向外取整后仍有间隙（`new_lower > new_upper`），规则选择**跳过**，
不宣称不可行；不可行仍由 simplex 的 Farkas 证书负责。

### 3.2 精确比例平行行

对同列模式的两行 `A_i, A_k`，令 `scale = A_k[0] / A_i[0]`，要求：

```text
A_k[m] == scale * A_i[m]              （逐位相等）
scale > 0:  scale*l_i == l_k, scale*u_i == u_k
scale < 0:  scale*u_i == l_k, scale*l_i == u_k
```

满足时两行定义同一可行域，删除其中一行。实现位于
`src/presolve/rules/duplicate_row.zig` 的 `patternHash`/`parallelScale`/
`boundsMatchScale`，集成在 `presolveImpl` 的重复行检测之后。

- 只接受**逐位精确**的倍数；近似平行（例如 `2` 与 `2.000000000000001`）不会被删，
  避免用容差误删约束。这是对“更激进”的 HiGHS 平行行检测的保守替代。
- 被删行的对偶乘子为 0；代表行保留原对偶。Farkas 证书中同理，等价约束的组合
  系数取 0 不改变不可行性证明。
- 负比例交换上下界；无界侧通过 `±inf` 的精确乘法比较处理。

## 4. 代码改动

| 文件 | 改动 |
|---|---|
| `src/presolve/rules/implied_fixing.zig` | 新增：项区间、隐含界、强制固定值、向外除法；含单元测试 |
| `src/presolve/rules/redundant_row.zig` | 将 `outwardAdd`/`outwardProduct` 暴露给新规则复用 |
| `src/presolve/rules/duplicate_row.zig` | 新增列模式哈希、精确比例、比例边界匹配 |
| `src/presolve/rules/root.zig` | 导出 `implied_fixing` |
| `src/presolve/presolve.zig` | `PresolveOptions`；`presolveWithOptions`/`presolveMipWithOptions`；`computeImpliedFixings`；固定值接入列替换；平行行删除；新计数器 |
| `src/presolve/root.zig` | 导出新选项/入口 |
| `src/solver/root.zig` | `hasPresolveReductions` 增加 forcing 潜力与平行行检测，保证 LP 路径调用新规则 |
| `build.zig` | `bench-presolve` / `build-bench-presolve` |
| `bench/presolve/presolve_bench.zig` | 五种模式的冷启动消融基准，输出 TSV |
| `bench/presolve/generate_instances.py` | 生成受控 forcing/parallel/MIP 算例 |
| `bench/presolve/compare.py` | 采样、聚合、JSON/Markdown 报告 |
| `bench/presolve/data/*.lp` | 受控算例 |

## 5. 基准方法与结果

### 5.1 方法

```powershell
zig build build-bench-presolve -Doptimize=ReleaseFast
python bench/presolve/generate_instances.py
python bench/presolve/compare.py --repetitions 10 --samples 3
```

- `off`：直接对原问题调用 simplex，不做 presolve。
- `legacy`：保留本次之前的规则，关闭 `implied_fixing` 与 `parallel_rows`。
- `implied`：只开启强制固定；`parallel`：只开启精确比例平行行。
- `full`：两者都开。
- 每个模式重复 10 次取均值，重复 3 个样本取中位数；解析在计时外，求解为冷启动。
- 目标值必须与 `off` 一致（`off` 成功时）；postsolve 后按原模型检查列界与行活动
  最大相对违反，报告 `violation`。
- LP 与 MIP 都只求解 **LP relaxation**，用于隔离 presolve 本身；MIP 的完整
  branch-and-bound 另有 `bench-mip` 入口，本次没有把 B&B 变量混入 presolve 结论。

### 5.2 结果

完整表见 `bench/presolve/results/presolve_results.md`，原始数据见
`bench/presolve/results/presolve_results.json`。关键行（Windows，ReleaseFast）：

| 算例 | 模式 | rows | cols | nnz | implied fixed | parallel rows | presolve us | solve us | total us | iters |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| forcing_lp | legacy | 40 | 60 | 80 | 0 | 0 | 51.6 | 660.4 | 717.2 | 61 |
| forcing_lp | implied | 0 | 0 | 0 | 60 | 0 | 27.4 | 10.5 | **40.6** | 0 |
| forcing_lp | full | 0 | 0 | 0 | 60 | 0 | 25.1 | 10.2 | **36.9** | 0 |
| parallel_lp | legacy | 40 | 13 | 80 | 0 | 0 | 25.4 | 170.9 | 197.7 | 15 |
| parallel_lp | parallel | 10 | 13 | 20 | 0 | 30 | 22.0 | 33.4 | **56.8** | 15 |
| parallel_lp | full | 10 | 13 | 20 | 0 | 30 | 53.5 | 62.4 | **116.9** | 15 |
| mip_forcing | legacy | 40 | 120 | 280 | 0 | 0 | 28.1 | 531.6 | 560.3 | 94 |
| mip_forcing | implied | 20 | 40 | 200 | 80 | 0 | 48.5 | 28.2 | **78.0** | 0 |
| mip_forcing | full | 20 | 40 | 200 | 80 | 0 | 68.6 | 35.9 | **106.2** | 0 |

解释：

- `forcing_lp` 的 60 个变量被两两/链式约束强制固定；`implied` 与 `full` 把模型降到
  0×0，总时间约为 legacy 的 1/19，且 `violation=0`。
- `parallel_lp` 中 30 行被识别为精确比例副本；删除后 nnz 从 80 降到 20，solve 时间
  约降为 1/5，总时间约 1/1.7。
- `mip_forcing` 的 80 个二元变量在整数舍入后唯一确定；`implied` 将行/列降到
  20/40，solve 时间从 531.6 us 降到 28.2 us，总时间约 1/5.3。
- 小模型/无可约简模型上，传播有固定开销：`presolve_mip.mps` 的
  `full/legacy = 0.75x`；`branching_compare.lp` 和 `warm_restart_mip.lp` 的比值在
  1.35x/0.87x 附近波动，主要是小样本噪声。因此后续需要 reduction limit/规模阈值。
- `off` 在较大随机稀疏 LP 上出现 `numerical_failure`（例如曾尝试的 80 行/160 列
  版本），说明当前直接 simplex 对某些结构仍不稳健；受控算例已缩小到 `off` 可解，
  以免把“presolve 让问题可解”和“presolve 提速”混为一谈。`parallel_lp` 的 40 行
  版本中 `off/legacy` 仍失败，而 `parallel/full` 成功，这是鲁棒性收益，但不计入
  上表的加速倍数。

## 6. 局限与后续工作

按优先级：

1. **通用 implied bound tightening 的对偶 postsolve**：当前只固定变量，不发布收紧
   列界。若要发布，需要在 postsolve 中记录每个界的来源行，并在最终解位于收紧界时
   把约化成本转移到对应行乘子；多行来源需要小型线性系统。实现前不应贸然开启。
2. **MIP 强化**：coefficient strengthening、probing、clique table、implied integer
   检测（Achterberg et al. 2020；Savelsbergh 1994）。这些通常比 LP 规则更能改善
   分支定界，但需要与现有 `src/mip` 的 domain propagation/clique 代码对齐，避免
   重复实现。
3. **列级消除**：singleton column stuffing、dominated column、doubleton equation、
   aggregator、parallel columns/dependent free columns（HiGHS 对应函数见第 2 节）。
4. **平行行容差策略**：当前要求逐位精确。可以引入相对容差 + 事后验证（例如在
   postsolve 检查代表行对被删行的最大违反），但这需要明确的数值风险政策。
5. **presolve 控制面**：reduction limit、规则开关/日志、按模型规模自动关闭强传播，
   以及把 postsolve 统计（对偶/约化成本/证书）纳入回归。
6. **更大规模基准**：在允许下载数据的环境中接入 MIPLIB 2017 / Netlib，并按
   “有无新规则”做配对统计；当前受控算例用于规则正确性和方向性验证，不应外推为
   通用性能结论。

## 7. 复现与验证

```powershell
# 单元测试（包含新增规则的数学、postsolve、pipeline 测试）
zig build test --summary all

# 构建基准
zig build build-bench-presolve -Doptimize=ReleaseFast

# 生成受控算例并跑消融
python bench/presolve/generate_instances.py
python bench/presolve/compare.py --repetitions 10 --samples 3

# 在任意 LP/MPS 集合上跑（例如自备 MIPLIB 路径）
python bench/presolve/compare.py --models path/to/a.mps path/to/b.lp
```

本次验证结果：`zig build test --summary all` 为 `332/332 tests passed`（Debug）；
ReleaseFast 下 `zhighs` CLI 对 `forcing_lp.lp`、`parallel_lp.lp`、`mip_forcing.lp`
均返回 `OPTIMAL`，目标值与基准一致（MIP 为 393，LP relaxation 为 403，属于不同
问题的目标值）。

