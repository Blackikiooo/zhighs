# zhighs 实现路线

任务按依赖顺序推进。每个阶段必须先满足验收条件，再开始依赖它的模块。

状态说明：

- `[x]` 已实现并通过当前测试。
- `[ ]` 尚未实现或仍需补充验证。

当前主线：完成真正的 `LpModelBuilder -> LpModel`，然后建立小规模 LP 正确性闭环。

## 0. 工程基线

- [x] 建立 `build.zig` 和 Zig 包配置。
- [x] 支持 `HInt` 32/64 位构建选项。
- [x] 实现 `HD`、`HCD` 及单元测试。
- [x] 建立 HCD 与 HiGHS C++ 的对照基准。
- [x] 将基础类型迁移到 `src/foundation/`。
- [x] 建立 `src/matrix/` 模块边界。
- [x] 为当前目录和规划模块补充 README。

验收命令：

```bash
zig build test
zig build test -Dhighs-int-width=w64
zig build bench-hcd -Doptimize=ReleaseFast
```

## 1. 稀疏数据结构

- [x] 定义强类型 `RowId` 和 `ColId`，集中实现 checked 边界转换。
- [x] 实现规范 `SparseVectorView/SparseVector` 和 `SparseVectorBuilder`。
- [x] 实现接受 triplet 的 SoA `MatrixBuilder`。
- [x] 实现构建后冻结的规范 CSC：排序、合并重复项、删除显式零。
- [x] 实现带 revision 检查的按需 CSR cache/view，并提供可信热路径。
- [x] 实现 CSC 的 `Ax`、`A^T y`，并区分 checked/AssumeValid 热路径。
- [x] 添加 CSC 随机矩阵乘法 property test 与结构边界测试。
- [x] 实现矩阵范数、绝对值统计、显式转置和规范行列切片。
- [x] 实现可逆 scaling、行列 permutation 及失败不产生半修改。
- [x] 实现函数式追加/删除行列和动态行 checkpoint/rollback。
- [x] 实现 generation-mark `SparseAccumulator`，服务于聚合、cuts 和 presolve。
- [x] 建立 `bench-matrix` ReleaseFast 专项基准。
- [x] 对照本地 HiGHS 实现 `range/assessment/alphaProductPlusY/productQuad` 等通用接口。
- [x] 将 partitioned row-wise、pricing、collectAj 和 basis update 明确归属 `lp/simplex`。
- [ ] 建立自动化 HiGHS C++ 矩阵差分测试，覆盖构建、乘法、scaling、切片和格式转换。
- [ ] 增加大型稀疏矩阵数据集 benchmark，记录时间、吞吐量、峰值内存和 cache 行为。
- [ ] 根据 benchmark 结果评估 SIMD、prefetch 和并行矩阵乘法；没有实测收益时不引入复杂度。

当前验收状态：随机矩阵运算与朴素实现一致，CSC/CSR 转换不改变矩阵语义；外部 HiGHS 自动差分和大型数据集性能验证待补。

## 2. 模型与解

- [x] 建立 `src/model/` 模块边界。
- [ ] 定义目标方向、变量上下界、行上下界和变量类型。
- [ ] 实现包含目标函数、行列上下界和变量类型的 `LpModelBuilder -> LpModel` 冻结流程。
- [ ] 实现 `Solution`、`Basis`、求解状态和信息结构。
- [ ] 实现 primal/dual residual 与 KKT 检查。
- [ ] 建立 HiGHS C API 小模型差分测试。

验收：内存模型能够表达 `L <= Ax <= U, l <= x <= u`，非法输入有明确错误。

## 3. LP 正确性闭环

- [ ] 实现只服务于小模型测试的 dense/reference simplex。
- [ ] 实现 Phase I、最优、不可行和无界判定。
- [ ] 对手工模型和随机小模型执行 HiGHS 差分测试。

验收：状态、目标值和可行性在统一容差内与 HiGHS 一致。

## 4. Revised Simplex 主线

- [ ] 建立 `src/nla/` 和 `src/lp/simplex/`。
- [ ] 定义 factorization 接口和 dense LU 参考实现。
- [ ] 实现 sparse LU、FTRAN/BTRAN 和 refactor。
- [ ] 实现 crash basis 与 primal revised simplex。
- [ ] 实现 dual revised simplex 和 basis warm start。
- [ ] 实现 Devex/steepest-edge pricing 与 Harris ratio test。
- [ ] 添加奇异 basis 检测、数值回退和性能基准。

验收：通过选定的 Netlib LP 回归集，warm start 对重复 LP 有稳定收益。

> PDLP 和 IPM 暂不进入主线。MIP 节点需要 basis 和快速 reoptimization，
> 因此先完成 revised simplex。

## 5. Presolve/Postsolve

- [ ] 先定义 reduction 记录和 `PostsolveStack`。
- [ ] 实现 empty row/column、fixed column、singleton row。
- [ ] 实现基础 bound tightening 与 redundant row。
- [ ] 为每条规则实现解和 basis 恢复。

验收：presolve 开关前后的状态和目标一致，恢复后的原模型解满足容差。

## 6. SCIP 式组件框架

- [ ] 建立 `Stage`、`Registry`、`Scheduler`、`Event` 和 `Services`。
- [ ] 定义统一插件生命周期和显式执行结果类型。
- [ ] 实现 presolver、branching、heuristic 三类最小接口。
- [ ] 添加 priority、frequency、max-depth 和 stage 合法性测试。
- [ ] 保证插件不能直接持有完整 `Solver`。

验收：组件可预测调度，非法阶段调用返回明确错误，热循环不使用动态分发。

## 7. 最小 MILP

- [ ] 实现整数可行性、变量域修改和可回滚的 domain stack。
- [ ] 实现 node、node queue、LP relaxation 和 incumbent。
- [ ] 实现 most-infeasible 与 pseudocost branching。
- [ ] 实现 rounding heuristic、cut pool 和一种基础割。
- [ ] 实现 node/time/gap limit。
- [ ] 建立小型 MIPLIB 回归集。

验收：小型 MILP 的状态和目标值与 HiGHS 一致，节点回溯无状态污染。

## 8. 后续扩展

- [ ] reliability/strong branching。
- [ ] clique、implication 和 conflict analysis。
- [ ] RINS、diving、feasibility pump。
- [ ] 更多 separator 与 constraint handler。
- [ ] MPS/LP 文件读写。
- [ ] PDLP、IPM、QP 和并行搜索。

## 当前建议执行顺序

1. 定义 LP 目标、行列上下界和变量类型。
2. 实现 `LpModelBuilder -> LpModel` 的验证、冻结和所有权测试。
3. 实现 `Solution/Basis` 与 primal/dual residual、KKT 检查。
4. 建立 dense/reference simplex 和 HiGHS C API 小模型差分测试。
5. 实现 NLA factorization、FTRAN/BTRAN，再进入 primal/dual revised simplex。
6. 在 LP 状态流稳定后实现 presolve/postsolve 和 SCIP 式组件框架。
7. 最后建立最小 branch-and-bound，并逐步加入 cuts、传播和启发式。
