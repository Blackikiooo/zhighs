好，这个想法很有意义。让我先深入了解下 SCIP 的架构设计，再做对比分析。

好的，让我深入对比 HiGHS 和 SCIP 的架构差异，然后为你的 Zig 项目设计一个新架构。

好，现在我有足够的信息了。让我做一个完整的分析。

---

## 一、HiGHS 架构的问题

HiGHS 的 C++ 代码有几个明显的工程缺陷：

### 1. 循环依赖严重

最典型的是 `HighsSparseMatrix`——它放在 `util/` 却包含了 `simplex/SimplexStruct.h`：

```cpp
// highs/util/HighsSparseMatrix.h:19
#include "simplex/SimplexStruct.h"  // util → simplex 反向依赖！
```

这说明 HiGHS 的分层是**假的**——名义上 util → lp_data → solver，实际有突破依赖方向的代码。

### 2. 没有插件体系

所有约束类型全是硬编码的 enum：

```cpp
enum PresolveRuleType : int {
  kPresolveRuleSingletonRow,  // 写死
  kPresolveRuleDominatedCol,  // 加新规则得改 enum  + switch
  // ...
};
```

想加一种新约束类型（比如 SOS2、Indicator），你得改至少 5 个文件。**SCIP 加一个约束处理器只需要写一个插件并注册进去。**

### 3. "上帝类" 问题

- `HEkk` 包揽了 simplex 的所有事情——iteration control、density tracking、DSE/Devex 切换
- `Highs` 类既是公共 API 又是内部调度器，职责太多
- 没有清晰的 stage 生命周期管理

### 4. 文件组织是平的，不是分层的

`highs/util/` 里有 40+ 个文件，从 `HFactor.cpp`（几千行的 LU 分解）到 `stringutil.cpp`（字符串工具）全混在一起。一个目录包含了好几个"层级"的代码。

### 5. 算法逻辑和数据结构耦合

`HighsSparseMatrix.cpp` 里混合了大量算法（addVec、addRows、scaling），不只是一个数据结构。

---

## 二、SCIP 架构为什么好

SCIP 的核心设计哲学：

```
┌───────────────────────────────────────────┐
│              SCIP Core Framework           │
│   (Tree management, stage control, LP     │
│    interface, event dispatch)              │
├───────────────────────────────────────────┤
│  Plugin Type          │  Interface         │
├───────────────────────────────────────────┤
│ Constraint Handler    │  CONSCHECK/ENFO... │  ← 领域逻辑
│ Separator             │  SEPEXECLP         │  ← 切平面
│ Branching Rule        │  BRANCHEXECLP      │  ← 分支策略
│ Primal Heuristic      │  HEUREXEC          │  ← 启发式
│ Presolver             │  PRESOLEXEC        │  ← 预处理
│ Event Handler         │  EVENTEXEC         │  ← 事件订阅
│ Propagator            │  PROPAEXEC         │  ← 域传播
└───────────────────────────────────────────┘
```

**关键设计决策：**

1. **插件模式**：所有算法组件都通过接口注册到核心框架
2. **回调驱动**：求解过程中的每个决策点都通过回调委托给插件
3. **分阶段管理**：SCIP 有明确的 stage transition（Init → Problem → Transform → Presolve → Solve → Free）
4. **职责分离**：核心框架只负责编排；所有领域逻辑在插件里
5. **Constraint Handler 是核心**：不是"求解器支持 N 种约束"，而是"每个约束类型自己定义 check/propagate/separate/enforce"

---

## 三、建议的 zhighs 架构

这是适合 Zig 的混血设计——**SCIP 的分层 + 插件接口 + HiGHS 的高效算法内核**：

```
src/
├── lib.zig                  # 库入口
├── main.zig                 # CLI
│
├── core/                    ★ 核心框架层 ★
│   ├── types.zig            # 基础类型 (HInt/HUInt — ✅ 已有)
│   ├── cdouble.zig          # 双倍精度 (HCD — ✅ 已有)
│   ├── constants.zig        # 常量、枚举 (HConst)
│   ├── problem.zig          # 问题表示 (约束+变量+目标)
│   ├── solution.zig         # 解 + basis
│   ├── setting.zig          # 参数系统 (HighsOptions)
│   ├── message.zig          # 日志 (HighsLogOptions)
│   ├── clock.zig            # 计时/时间限制
│   ├── statistics.zig       # 求解统计
│   ├── event.zig            ★ 事件系统 ★
│   └── stage.zig            ★ 阶段机 (SCIP stage) ★
│
├── interface/               ★ 插件接口定义层 ★
│   ├── constraint.zig       # 约束处理器接口
│   ├── branching.zig        # 分支规则接口
│   ├── separator.zig        # 切平面接口
│   ├── heuristic.zig        # 启发式接口
│   ├── presolver.zig        # 预处理器接口
│   ├── propagator.zig       # 传播器接口
│   ├── event_handler.zig    # 事件处理器接口
│   └── lp_solver.zig        # LP 求解器接口 (抽象层)
│
│
├── lp/                      ★ LP 求解引擎 ★
│   ├── factor.zig           # LU 分解 (原 HFactor)
│   ├── simplex.zig          # Simplex (原 HEkk)
│   ├── ipm.zig              # 内点法
│   └── pdlp.zig             # 一阶方法
│
├── mip/                     ★ MIP 求解框架 ★
│   ├── relax.zig            # LP relaxation 管理
│   └── branch_and_cut.zig   # 主 B&C 循环
│
├── plugins/                 ★ 内置插件实现 ★
│   ├── constraint/
│   │   ├── linear.zig       # 线性约束
│   │   ├── knapsack.zig
│   │   ├── set_packing.zig
│   │   ├── set_partition.zig
│   │   ├── set_cover.zig
│   │   ├── sos1.zig
│   │   ├── indicator.zig
│   │   ├── cardinality.zig
│   │   └── quadratic.zig
│   ├── branching/
│   │   ├── most_infeasible.zig
│   │   ├── pseudo_cost.zig
│   │   ├── strong.zig
│   │   ├── reliability.zig
│   │   ├── inference.zig
│   │   └── hybrid.zig       ★ 根据深度自动切换 ★
│   ├── separator/
│   │   ├── gomory.zig
│   │   ├── clique.zig
│   │   ├── mir.zig
│   │   ├── flow_cover.zig
│   │   ├── implied_bound.zig
│   │   ├── zero_half.zig
│   │   ├── rlt.zig
│   │   └── proximity.zig
│   ├── heuristic/
│   │   ├── trivial.zig
│   │   ├── rounding.zig
│   │   ├── diving.zig
│   │   ├── feaspump.zig
│   │   ├── rins.zig
│   │   ├── ren.zig
│   │   ├── local_branching.zig
│   │   ├── mutation.zig
│   │   ├── repair.zig
│   │   └── shift_propagate.zig
│   ├── presolver/
│   │   ├── empty_row_col.zig
│   │   ├── singleton.zig
│   │   ├── dominated_col.zig
│   │   ├── free_col_sub.zig
│   │   ├── doubleton_eq.zig
│   │   ├── forcing_chain.zig
│   │   ├── dependent_eq.zig
│   │   ├── sparsify.zig
│   │   └── probing.zig
│   └── event_handler/
│       ├── log_incumbent.zig
│       └── empathy_gap.zig
│
├── data/                    ★ 数据结构层 ★
│   ├── vector.zig           # HVectorBase<T>
│   ├── sparse_matrix.zig    # HighsSparseMatrix
│   ├── sparse_vector.zig
│   ├── hash.zig             # HighsHash
│   ├── sort.zig             # HighsSort
│   ├── rbtree.zig           # HighsRbTree
│   ├── splay.zig            # HighsSplay
│   ├── disjoint_sets.zig
│   ├── matrix_utils.zig
│   └── random.zig
│
└── io/                      ★ 文件读写 ★
    ├── mps_reader.zig
    ├── lp_reader.zig
    └── file_reader.zig      # 基础文件解析
```

---

## 四、Zig 接口设计模式

关键是怎么用 Zig 的 features 实现 SCIP 式的插件系统。核心思路是用 **vtable struct** + **type-erased context pointer**：

### constraint 接口示例

```zig
// interface/constraint.zig

pub const CheckResult = enum {
    feasible,
    infeasible,
    unknown,
};

pub const EnforceResult = enum {
    feasible,
    branched,
    cut_added,
    domain_reduced,
    infeasible,
};

pub const ConstraintHandler = struct {
    name: []const u8,
    description: []const u8,
    priority: i32 = 0,

    // vtable — 所有回调
    // 核心: 检查可行性 (必须实现)
    check: *const fn (ctx: *anyopaque, sol: *const Solution) CheckResult,

    // 可选回调: 每个都有默认 null
    presolve: ?*const fn (ctx: *anyopaque) ?PresolveResult = null,
    propagate: ?*const fn (ctx: *anyopaque) ?PropagateResult = null,
    separate_lp: ?*const fn (ctx: *anyopaque, sol: *const Solution, pool: *CutPool) void = null,
    enforce_lp: ?*const fn (ctx: *anyopaque, sol: *const Solution) EnforceResult = null,
    enforce_ps: ?*const fn (ctx: *anyopaque) EnforceResult = null,
    active: ?*const fn (ctx: *anyopaque, node: *const Node) void = null,
    deactive: ?*const fn (ctx: *anyopaque, node: *const Node) void = null,
};

/// 将具体类型包装为 ConstraintHandler 的辅助函数
pub fn wrap(comptime T: type, ctx: *T, comptime impl: T) ConstraintHandler {
    return .{
        .name = impl.name,
        .description = impl.description,
        .priority = impl.priority,
        .check = struct {
            fn f(c: *anyopaque, sol: *const Solution) CheckResult {
                return @as(*T, @ptrCast(@alignCast(c))).check(sol);
            }
        }.f,
        .presolve = if (impl.presolve) |_| struct {
            fn f(c: *anyopaque) ?PresolveResult {
                return @as(*T, @ptrCast(@alignCast(c))).presolve();
            }
        }.f else null,
        // ... 其他回调同样包装
    };
}
```

### 使用方式

```zig
// plugins/constraint/knapsack.zig

const KnapsackHandler = struct {
    name: []const u8 = "knapsack",
    description: []const u8 = "Knapsack constraint handler",
    priority: i32 = 100,

    weights: []const f64,
    capacity: f64,

    fn check(self: *@This(), sol: *const Solution) CheckResult {
        var sum: f64 = 0;
        for (self.weights, sol.col_value) |w, x| {
            if (x > 0.5) sum += w;
        }
        return if (sum <= self.capacity + 1e-9) .feasible else .infeasible;
    }

    fn separate_lp(self: *@This(), sol: *const Solution, pool: *CutPool) void {
        // 生成 lifted cover cut
    }
};

pub fn create(weights: []const f64, capacity: f64) ConstraintHandler {
    return ConstraintHandler.wrap(KnapsackHandler, .{
        .weights = weights,
        .capacity = capacity,
    });
}
```

### branching 接口

```zig
// interface/branching.zig

pub const BranchResult = union(enum) {
    branched: struct { children: []ChildNode },
    cutoff,
    did_not_run,
};

pub const BranchingRule = struct {
    name: []const u8,
    priority: i32,
    max_depth: i32 = -1, // -1 = unlimited

    execute_lp: *const fn (ctx: *anyopaque, sol: *const Solution) BranchResult,
    execute_ps: ?*const fn (ctx: *anyopaque) BranchResult = null,
};
```

### 事件系统

```zig
// core/event.zig

pub const EventType = enum {
    node_processed,
    lp_solved,
    incumbent_found,
    branch_performed,
    separator_found_cut,
    presolve_finished,
    /// ...
};

pub const EventHandler = struct {
    name: []const u8,
    interests: []const EventType,

    execute: *const fn (ctx: *anyopaque, event: EventType, data: *anyopaque) void,
};

// 事件分发
pub const EventSystem = struct {
    handlers: std.ArrayList(RegisteredHandler),

    pub fn fire(self: *Self, event: EventType, data: anytype) void {
        for (self.handlers.items) |handler| {
            handler.fire(event, &data);
        }
    }
};
```

---

## 五、对比分析

| 维度 | HiGHS (原版) | SCIP | zhighs (新架构) |
|---|---|---|---|
| **依赖方向** | 循环依赖 | 严格单向 | 严格单向 `core → interface → {data, lp} → plugins` |
| **约束系统** | 硬编码 enum + switch | 约束处理器插件 | 约束处理器接口 + 内置实现 |
| **分支规则** | 内嵌在 HEkk 里 | 分支规则插件 | 分支规则接口 + 多个策略 |
| **可扩展性** | 改核心代码 | 写新插件注册 | 写新插件，编译期注册 |
| **状态管理** | 隐式 | 显式 stage 机 | 显式 stage 机 |
| **文件组织** | 扁平 | 按功能分层 | 按职责分层 |
| **Zig 特色利用** | 无 | N/A | comptime 泛型、vtable、tagged union |

---

## 六、你现在已有的代码怎么调整

你现在有 `src/types/root.zig`、`types/double.zig`、`types/int.zig`。迁移路径：

**第一步**（立即）：移动现有代码到新结构
```
// 当前                    →  新位置
types/int.zig              →  core/types.zig
types/double.zig           →  core/cdouble.zig
types/root.zig             →  合并进 lib.zig
builds/types_build.zig     →  保留，调整路径
src/sparse/root.zig        →  data/sparse_matrix.zig
```

**第二步**：写 `core/constants.zig`（HConst 迁移）

**第三步**：写 `interface/constraint.zig` 和 `interface/branching.zig`

**第四步**：实现 `data/vector.zig`（HVectorBase）

**第五步**：实现 `data/sparse_matrix.zig`

**第六步**：实现 `core/problem.zig`

**第七步**：实现 `lp/pdlp.zig` 作为第一个完整求解器

---

## 七、要不要我帮你落地？

如果你对我的方案认可，我可以直接开始帮你：

1. **重构现有代码**——把你已有的 `types/int.zig` 和 `types/double.zig` 迁移到新目录结构
2. **写核心接口**——`interface/constraint.zig`、`interface/branching.zig` 等
3. **写第一个求解器**——PDLP 或 Simplex

你觉得这个方向怎么样？有什么想调整的地方吗？