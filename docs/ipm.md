# 连续 LP 内点法

`Method=2` 将连续线性模型交给 Mehrotra 预测—校正内点法，贯通 LP/MPS 读取、标准形变换、迭代求解、原坐标解恢复、模型属性和解文件导出。

```powershell
zig build run -Doptimize=ReleaseFast -- --param Method=2 --solution solution.sol test/model/data/pipeline_lp.lp
zig build run -Doptimize=ReleaseFast -- --param Method=2 test/model/data/ipm_lp.mps
zig build test-ipm
zig build test-end-to-end
zig build bench-ipm -Doptimize=ReleaseFast
```

库调用使用 `zhighs.ipm.solve(allocator, problem, options)`，结果由调用者 `deinit()`。计时或时间限额使用 `solveWithStats(allocator, problem, options, io, stats)`；没有时钟却指定时间限额会返回 `ClockRequired`。非法数值选项返回 `InvalidOptions`，非法模型返回 `invalid_model` 状态。

模型 API 设置 `env.setIntParam("Method", 2)` 后调用 `model.optimize()`。发布 `X`、`Pi`、`RC`、`Slack`、目标、解数量和 `BarIterCount`；结果不包含 simplex basis，不伪造 crossover。重新求解失败或达到限额时清除旧解。支持 `IterationLimit`、`BarIterLimit`、`TimeLimit`、`WorkLimit`、中断标志及 `.barrier` 回调；默认内点迭代上限为 100，显式 `BarIterLimit` 可覆盖。此路径的 WorkLimit 每次预测—校正迭代记一单位。时间与中断在迭代边界检查，不能抢占单次分解，时间限额不含标准形构造和内存分配。

## 性能实现

- 从标准形 CSC 列结构分析 `A Θ Aᵀ` 下三角的行轮廓，分析一次，重复使用存储和工作向量。迭代不分配内存。
- 只装配下三角。Cholesky 在轮廓内计算，前代和回代均按连续行访问；预测步和校正步复用同一次数值分解。
- 保留 `Options.dense_normal=true` 的原稠密 Cholesky 内核，便于对照。两后端共享标准形、残差与步长逻辑；这是后端比较，不是未经修改的整仓历史版本比较。
- 正则化重试有次数上限；校验增长因子，避免无穷重试。NaN/Inf 残差不能满足收敛条件。
- 修复原有计时入口引用缺失字段的问题，统计装配、分解、方向、残差、迭代总耗时及正规矩阵存储项数。

轮廓存储依赖无需选主元时消元不会产生轮廓外填充的性质，参见 [Netlib Skyline Storage](https://www.netlib.org/utk/people/JackDongarra/etemplates/node378.html)。每行从最左结构元素存到对角线，包含内部零元素；它不等于通用压缩稀疏 Cholesky。

## 可复现性能证据

环境：本次 Windows x64 工作区、Zig 0.16.0、ReleaseFast、单线程；CPU 型号未获取。命令：`zig build bench-ipm -Doptimize=ReleaseFast --global-cache-dir .zig-cache/global`。

基准为偶数 n 的 `max Σxᵢ`，约束 `xᵢ+xᵢ₊₁≤2, x≥0`，解析最优值 n。每个后端预热一次，随后取三次完整求解耗时中位数，包含标准形与分配，排除编译和进程启动。分阶段列来自最后一次测量，并非各阶段中位数。

| n | 稠密中位数 ms | 轮廓中位数 ms | 速度比 | 正规矩阵项数：稠密 / 轮廓 |
| --- | ---: | ---: | ---: | ---: |
| 128 | 2.153 | 0.145 | 14.8× | 16,129 / 253 |
| 512 | 166.717 | 0.764 | 218.2× | 261,121 / 1,021 |

以上为本次第一组测量，后续重复运行会随调度和频率变化。四组均 6 次迭代；目标分别约 127.99999999 和 511.99999875，最大原始可行性违反为 0。基准对目标和可行性执行断言。存储项数只包括一个正规矩阵缓冲区，因子还需要同等大小的浮点缓冲区，另有索引与迭代向量。

这组模型具有窄轮廓，能明显展示跳过结构零的收益，不能外推为任意 LP 的速度提升。

## 验证与边界

内点法测试覆盖上下界、范围行、自由/固定变量、最大化与目标偏移、对偶符号、稠密/轮廓目标对照、具有内部零及填充的线性系统残差、多右端求解、常量约束、非法选项、限额、中断和所有分配失败位置的清理。模型测试检查对偶/松弛/无基结果、旧解清理和回调中断；端到端检查 LP 与 MPS 的状态、目标及导出解。

最终执行 `zig build test -j2 --global-cache-dir .zig-cache/global --summary all`：50/50 构建步骤、317/317 测试通过，其中内点法 17 项，模型/端到端测试入口 12 项，另含命令行输出及解文件检查。修改文件的 `zig fmt --check` 与 `git diff --check` 均通过。

另外使用完整的 40 模型 SoPlex 真实 LP 集进行独立验证，规模从 8x8 到
536x1775。验证器通过 HiGHS 单独解析原始 MPS，再检查原始坐标下的边界、
stationarity、对偶符号、互补性、primal-dual gap 及参考目标。最终单一配置
严格通过 27/40，即参考最优模型中的 27/30；初始实现为 18/40。全部输入哈希
与冻结清单匹配，没有 false-optimal。完整命令、逐项结果和剩余失败见
[`bench/ipm/real_corpus_analysis.md`](../bench/ipm/real_corpus_analysis.md)。

当前仍有限制：

- 面向连续 LP；MIP 保持原有分支定界路径，QP 不支持。
- 没有 AMD/嵌套剖分排序或超节点分解；宽轮廓时存储仍可达到 O(m²)、分解 O(m³)，列上界转换出的辅助行也可能扩大轮廓。
- 没有自动缩放、齐次自对偶嵌入、通用不可行/无界证书或 crossover。一般发散返回数值失败或迭代限额，不将启发式现象当作证明。可直接证明的常量矛盾及无约束负成本方向单独处理。
- `Method=2` 当前直接求解原模型，不接入 simplex presolve、warm start 或日志回调；`Presolve`、`Crossover`、`BarOrder`、`BarHomogeneous` 在此路径尚未实现。

因此本次实现提供了可验证的性能改进和可用端到端路径，还不能据此宣称达到 HiGHS/IPX 的通用大规模鲁棒性。
