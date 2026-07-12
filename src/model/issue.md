# 模型层已知问题与待办

## 已完成

- [x] objects 里面承载了太多的功能，需要 vars、models 等全部拆分出来，单独实现，避免成为单文件实现。
- [x] attrs 在 .zig 里面的实现要变得更加 zig-like。
- [x] 建立核心模型 IR 分类体系：
  - LinearModel（LP/MILP）
  - QuadraticModel + Hessian（QP/MIQP/QCP/MIQCP）
  - ExpressionGraph + NonlinearModel（NLP/MINLP）
  - CompiledModel（union 分派）
  - ProblemClass 自动分类
  - SolverCapability 能力查询
  - 分层 Validation（6+ 个 validator）
  - w32/w64 独立测试（`zig build test-model`）

## 待办（下一阶段）

- [ ] `Model -> CompiledModel` 编译路径（`compile_model.zig`）。
- [ ] Solution/Basis/Info 数据结构。
- [ ] primal/dual residual 与 KKT 检查。
- [ ] 求解 dispatcher 接入（目前 optimize() 返回 FeatureNotAvailable）。
- [ ] ExpressionGraph 当前能力边界：节点操作集（16 种 opcode）覆盖基本数学运算，但缺少循环、条件、自定义函数等表达。自动微分、区间传播、求值优化尚未实现。

## API 稳定性

- 所有核心模型类型标记为 **Experimental API**。
- 数学语义（bounds 语义、目标方向编码、1/2 xᵀQx 约定）已固定。
- 字段布局、分配策略和 convenience 函数可能在 presolve/simplex 集成时修正。
- `NonlinearModel` 即使没有 nonlinear root 也能构造（problemClass 返回 NLP）；调用者应通过 `hasNonlinearRoot()` 检查。

## 所有权风险

- `CompiledModel` 是 move-only；复制 union 后调用 `deinit()` 导致 double-free。
- `clone()` 未实现（记录在案，待下一阶段）。
- `LinearModelBuilder.freeze()` 后将数组移出，`deinit()` 检查 `frozen` 标志避免 double-free。

## 曲率分析

- `Curvature` 枚举已定义，默认值为 `.unknown`。
- 当前没有任何代码执行凸性检查（特征值分析未实现）。
- `validateQuadraticModel` 不会因为 curvature 为 `.unknown` 而拒绝模型。
