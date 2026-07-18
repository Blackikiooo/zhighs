# Model 层 — 商业求解器 API 对标进度

本文档以主流商业求解器 C API 为参照，逐项记录 `src/model/` 的实现进度。

状态说明：

- `✅` 已实现，通过测试，可用于生产
- `🧪` 已实现，但仅对 continuous LP 路径生效
- `🪣` 已有 stub（`FeatureNotAvailable`），数据存储就绪但求解路径未接入
- `❌` 尚未实现或完全缺失

## 1. 环境管理 (Environment)

| C API | 状态 | zhighs 对应 | 说明 |
| --- | --- | --- | --- |
| `loadenv` | ✅ | `Env.initSimple` | 支持基础参数、日志回调 |
| `emptyenv` | ❌ | — | 空环境构造未实现 |
| `startenv` | ❌ | — | 延迟启动未实现 |
| `freeenv` | ✅ | `Env.deinit` | |
| `getconcurrentenv` | ❌ | — | 并发优化环境未实现 |
| `getmultiobjenv` | ❌ | — | 多目标环境未实现 |
| `discardconcurrentenvs` | ❌ | — | |
| `discardmultiobjenvs` | ❌ | — | |

## 2. 模型创建与销毁

| C API | 状态 | zhighs 对应 | 说明 |
| --- | --- | --- | --- |
| `newmodel` | ✅ | `Model.init` | 创建空模型 |
| `loadmodel` | ✅ | `Model.read` / `Model.readModel` | 通过 `io` 模块加载 MPS/LP |
| `copymodel` | ✅ | `Model.copy` | 独立深拷贝 |
| `copymodeltoenv` | 🪣 | `model_advanced.copyModelToEnv` | Stub，返回 `FeatureNotAvailable` |
| `freemodel` | ✅ | `Model.deinit` | |

## 3. 变量与线性约束 (Linear Modelling)

| C API | 状态 | zhighs 对应 | 说明 |
| --- | --- | --- | --- |
| `addvar` | ✅ | `Model.addVar` | 单变量追加，lazy-update |
| `addvars` | ✅ | `Model.addVars` | 批量追加 |
| `addconstr` | ✅ | `Model.addConstr` | |
| `addconstrs` | ✅ | `Model.addConstrs` | |
| `addrangeconstr` | ✅ | `Model.addRangeConstr` | 双边范围约束 |
| `addrangeconstrs` | ✅ | `Model.addRangeConstrs` | |
| `chgcoeffs` | ✅ | `Model.chgCoeff` / `Model.chgCoeffs` | 单元素 + 批量，last-write-wins |
| `delvars` | ✅ | `Model.delVars` | 支持 SOS/QC/GenConstr/PWLObj remap |
| `delconstrs` | ✅ | `Model.delConstrs` | |
| `Xaddvars` | ❌ | — | 高级变量追加（多列整块）未实现 |
| `Xaddconstrs` | ❌ | — | |
| `Xaddrangeconstrs` | ❌ | — | |
| `Xchgcoeffs` | ❌ | — | |
| `Xloadmodel` | ❌ | — | 高级模型加载未实现 |
| `updatemodel` | ✅ | `Model.updateModel` / `Model.applyPending` | lazy-update flush |

## 4. 变量/约束属性 (Scalar Attribute Access)

当前 `Model` 支持 string-keyed 属性系统（通过 `attrs.zig` 的 `Attr` 枚举），`get/set*Attr` 系列完整实现。

| 属性族 | get | set | getElement | setElement | getArray | setArray | getList | setList |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Int | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| Dbl | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| Str | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| Char | — | — | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |

覆盖的核心属性：`LB`, `UB`, `Obj`, `RHS`, `Sense`, `VType`, `VarName`, `ConstrName`, `X`, `RC`, `Pi`, `Slack`, `VBasis`, `CBasis`, `Status`, `ObjVal`, `ObjBound`, `IterCount`, `NodeCount`, `BarIterCount`, `ModelName`, `ModelSense`, `NumVars`, `NumConstrs`, `NumNZs`, `Start`, `PStart`, `DStart`, `QVal`, `PWLObjCvx`, `Fingerprint`。

## 5. 二次目标与约束 (Quadratic)

| C API | 状态 | zhighs 对应 | 说明 |
| --- | --- | --- | --- |
| `addqpterms` | ✅ | `Model.addQPterms` / `addQPtermsExpr` | 追加 Q 项到目标，支持 Var handle |
| `delq` | ✅ | `Model.delQ` | 清空所有 Q 项 |
| `getq` | ✅ | `Model.getQ` | 读取 Q 项到 caller buffer |
| `addqconstr` | ✅ | `Model.addQConstr` | SoA 存储（Q term + linear term） |
| `delqconstrs` | ✅ | `Model.delQConstrs` | 压缩删除 |
| `getqconstr` | ✅ | `Model.getQConstr` | |
| `getqconstrbyname` | ✅ | `Model.getQConstrByName` | |

> 注意：QP 求解路径未实现。`Model.optimize` 遇到 QP/MIQP 返回 `FeatureNotAvailable`。二次数据的存储、查询和删除已完整就绪。

## 6. 一般约束 (General Constraints)

共 20 种一般约束类型，zhighs 全部实现了 **add + get**，存储在 `genconstr_*` SoA 中。数据层完整。

| C API | add | get | 数据 |
| --- | --- | --- | --- |
| `addgenconstrMax` | ✅ | ✅ | ✅ |
| `addgenconstrMin` | ✅ | ✅ | ✅ |
| `addgenconstrAbs` | ✅ | ✅ | ✅ |
| `addgenconstrAnd` | ✅ | ✅ | ✅ |
| `addgenconstrOr` | ✅ | ✅ | ✅ |
| `addgenconstrIndicator` | ✅ | ✅ | ✅ (packed bitcast 存储) |
| `addgenconstrPWL` | ✅ | ✅ | ✅ (packed bitcast 存储) |
| `addgenconstrPoly` | ✅ | ✅ | ✅ |
| `addgenconstrExp` | ✅ | ✅ | ✅ |
| `addgenconstrExpA` | ✅ | ✅ | ✅ (base-a) |
| `addgenconstrLog` | ✅ | ✅ | ✅ |
| `addgenconstrLogA` | ✅ | ✅ | ✅ (base-a) |
| `addgenconstrLogistic` | ✅ | ✅ | ✅ |
| `addgenconstrPow` | ✅ | ✅ | ✅ |
| `addgenconstrSin` | ✅ | ✅ | ✅ |
| `addgenconstrCos` | ✅ | ✅ | ✅ |
| `addgenconstrTan` | ✅ | ✅ | ✅ |
| `addgenconstrNorm` | ✅ | ✅ | ✅ |
| `addgenconstrNL` | ✅ | ✅ | ✅ (generic container) |
| `delgenconstrs` | ✅ | — | 对应 `Model.delGenConstrs`（清空全部） |

> 注意：全部 20 种一般约束的 solver 处理路径均未实现。`Model.optimize` 遇到有 genconstr 的模型返回 `FeatureNotAvailable`。

## 7. SOS 约束 (Special Ordered Sets)

| C API | 状态 | zhighs 对应 | 说明 |
| --- | --- | --- | --- |
| `addsos` | ✅ | `Model.addSOS` | SOS1/SOS2 |
| `delsos` | ✅ | `Model.delSOS` | 压缩删除 |
| `getsos` | ✅ | `Model.getSOS` | |

> 求解路径未接入。

## 8. 多目标与 PWL 目标

| C API | 状态 | zhighs 对应 | 说明 |
| --- | --- | --- | --- |
| `setobjectiven` | 🪣 | `model_advanced.setObjectiveN` | Stub |
| `setpwlobj` | ✅ | `Model.setPWLObj` | PWL 目标设置（SoA 存储） |
| `getpwlobj` | ✅ | `Model.getPWLObj` | |

## 9. 模型求解 (Optimization)

| C API | 状态 | zhighs 对应 | 说明 |
| --- | --- | --- | --- |
| `optimize` | 🧪 | `Model.optimize` | **仅 continuous LP** 路径完整。MILP/QP/MIQP/NLP 均返回 `FeatureNotAvailable` |
| `optimizeasync` | ❌ | — | 异步求解未实现 |
| `presolvemodel` | 🪣 | `model_advanced.presolveModel` | Stub（仅 set status） |
| `computeIIS` | 🪣 | `model_advanced.computeIIS` | Stub（仅 set status） |
| `converttofixed` | ✅ | `model_advanced.convertToFixed` | 修复整数变量到当前解 |
| `feasrelax` | 🪣 | `model_advanced.feasRelax` | Stub |
| `fixmodel` | ✅ | `model_advanced.fixModel` | 别名，委托到 `convertToFixed` |
| `reset` | ✅ | `Model.reset` | 清除解状态 |
| `sync` | ❌ | — | 多线程同步未实现 |

## 10. 回调系统 (Callbacks)

| C API | 状态 | 说明 |
| --- | --- | --- |
| `setcallbackfunc` | 🧪 | 已实现，simplex 迭代中触发。仅 `.simplex` where |
| `getcallbackfunc` | ✅ | |
| `setlogcallbackfunc` | ✅ | `model_params.setLogCallbackFunc` |
| `cbget` | 🪣 | Stub — MIP 回调信息查询未实现 |
| `cbcut` | 🪣 | Stub — 延迟 cut pool 添加未实现 |
| `cblazy` | 🪣 | Stub — 延迟约束添加未实现 |
| `cbsolution` | 🪣 | Stub — MIP 启发式解提交未实现 |
| `cbproceed` | 🪣 | Stub |
| `cbstoponemultiobj` | 🪣 | Stub |
| `cbsetdblparam` | 🪣 | Stub |
| `cbsetintparam` | 🪣 | Stub |
| `cbsetstrparam` | 🪣 | Stub |
| `cbsetparam` | 🪣 | Stub |
| `terminate` | ✅ | `Model.terminate` | 设置 interrupt 原子标志 |

## 11. 参数管理 (Parameters)

| C API | 状态 | zhighs 对应 | 说明 |
| --- | --- | --- | --- |
| `getintparam` | ✅ | `Model.getIntParam` | |
| `setintparam` | ✅ | `Model.setIntParam` | |
| `getdblparam` | ✅ | `Model.getDblParam` | |
| `setdblparam` | ✅ | `Model.setDblParam` | |
| `getstrparam` | ✅ | `Model.getStrParam` | |
| `setstrparam` | ✅ | `Model.setStrParam` | |
| `setparam` | ✅ | `Model.setParam` | 统一入口（ParamValue union） |
| `getdblparaminfo` | ✅ | `Model.getDblParamInfo` | |
| `getintparaminfo` | ✅ | `Model.getIntParamInfo` | |
| `getstrparaminfo` | ✅ | `Model.getStrParamInfo` | |
| `readparams` | ✅ | `Model.readParams` | |
| `writeparams` | ✅ | `Model.writeParams` | |
| `resetparams` | ✅ | `Model.resetParams` | |

已支持的参数：`TimeLimit`, `WorkLimit`, `IterationLimit`, `OutputFlag`, `SimplexPricing`, `SimplexLogInterval`, `Threads`。其余参数注册为 unknown。

## 12. 模型查询 (Model Queries)

| C API | 状态 | zhighs 对应 | 说明 |
| --- | --- | --- | --- |
| `getcoeff` | ✅ | `Model.getCoeff` | |
| `getconstrbyname` | ✅ | `Model.getConstrByName` | |
| `getconstrs` | ✅ | `Model.getConstrs` | |
| `getvarbyname` | ✅ | `Model.getVarByName` | |
| `getvars` | ✅ | `Model.getVars` | |
| `getenv` | ✅ | `Model.getEnv` | |
| `getjsonsolution` | ✅ | `model_advanced.getJSonSolution` | 基础 JSON 序列化 |
| `singlescenariomodel` | 🪣 | `model_advanced.singleScenarioModel` | Stub |
| `Xgetconstrs` | ❌ | — | |
| `Xgetvars` | ❌ | — | |

## 13. 文件读写 (I/O)

| C API | 状态 | 说明 |
| --- | --- | --- |
| `readmodel` | ✅ | `Model.read` / `Model.readModel`，支持 LP/MPS |
| `read` | ✅ | 别名，委托到 `io` 模块 |
| `write` | ✅ | `Model.write` / `Model.writeModel`，支持 LP/MPS |

## 14. 高级 Simplex 例程

| C API | 状态 | 说明 |
| --- | --- | --- |
| `FSolve` | ❌ | 单次 FTRAN，未暴露 |
| `BSolve` | ❌ | 单次 BTRAN，未暴露 |
| `BinvColj` | ❌ | 基逆的列 j 提取，未暴露 |
| `BinvRowi` | ❌ | 基逆的行 i 提取，未暴露 |
| `getBasisHead` | ✅ | `Model.getBasisHead` |

## 15. 批处理 (Batch API)

| C API | 状态 | 说明 |
| --- | --- | --- |
| 全部 12 个 Batch API | ❌ | 未实现。`abortbatch`, `discardbatch`, `freebatch`, `getbatch`, `getbatchenv`, `getbatchintattr`, `getbatchjsonsolution`, `getbatchstrattr`, `optimizebatch`, `retrybatch`, `updatebatch`, `writebatchjsonsolution` |

## 16. 参数调优 (Tuning)

| C API | 状态 | 说明 |
| --- | --- | --- |
| `tunemodel` | 🪣 | `Model.tuneModel` — stub |
| `gettuneresult` | 🪣 | `Model.getTuneResult` — stub |

## 17. 错误处理

| C API | 状态 | 说明 |
| --- | --- | --- |
| `geterrormsg` | ✅ | `Model.getErrormsg` |

## 18. 版本

| C API | 状态 | 说明 |
| --- | --- | --- |
| `version` | ✅ | `Model.version` |

---

## 进度总结

### 数据层（存储 + add/get/del）

```
Environment            ████████░░ 80%
线性建模                ██████████ 100%
变量/约束属性 (Attr)     ██████████ 100%
二次目标/约束 (QP/QC)    ██████████ 100%
一般约束 (GenConstr)×20  ██████████ 100%  全部 add+get 已实现
SOS 约束                ██████████ 100%
PWL 目标                ██████████ 100%
参数管理 (Param)         ██████████ 100%
文件 I/O                ██████████ 100%
模型查询 (Query)         ██████████ 100%
```

### 求解层（optimize 实际处理）

```
Continuous LP          ██████████ 100%  simplex 引擎完整接入
MILP (branch&bound)    ░░░░░░░░░░   0%  返回 FeatureNotAvailable
QP / MIQP              ░░░░░░░░░░   0%  返回 FeatureNotAvailable
QCP / MIQCP            ░░░░░░░░░░   0%  返回 FeatureNotAvailable
NLP / MINLP            ░░░░░░░░░░   0%  返回 FeatureNotAvailable
```

### 高级特性

```
回调系统                ████░░░░░░  40%  iteration 回调 ✓, MIP 回调 ✗
Presolve               █░░░░░░░░░  10%  model 入口 stub 存在, 无实际 reduction
IIS                    █░░░░░░░░░  10%  stub
FeasRelax              █░░░░░░░░░  10%  stub
Multi-objective        ░░░░░░░░░░   0%  stub
Concurrent env         ░░░░░░░░░░   0%  未实现
Batch API              ░░░░░░░░░░   0%  未实现
Tuning                 ░░░░░░░░░░   0%  stub
Advanced Simplex       ██░░░░░░░░  20%  getBasisHead ✓, FTRAN/BTRAN 未暴露
```

### 关键技术债

1. **19 个一般约束的 solver 处理路径全部缺失** — 数据存储完整但 optimizer 无法消费。优先级：低（需要 MIP/NLP solver 先就位）。

2. **MIP 回调全部 stub** — `cbCut`, `cbLazy`, `cbSolution`, `cbGet`, `cbProceed` 均返回 `FeatureNotAvailable`。优先级：中（MIP solver 就位前无法实现）。

3. **`setObjectiveN` stub** — 多目标优化的模型层入口缺失。优先级：低。

4. **`FSolve`/`BSolve`/`BinvColj`/`BinvRowi` 未暴露** — Simplex 引擎内部已有这些能力（`factorization.solve`/`solveTranspose`/`solveForUpdate`），但未通过公开 API 暴露。优先级：低（高级用户功能）。

5. **Batch API 完全缺失** — 12 个 batch API 函数均未实现。优先级：低（需先有稳定求解路径）。

6. **`optimizeasync` / `sync` 缺失** — 异步求解 + 线程同步未实现。优先级：低（需并行框架先就位）。

### 当前建议执行顺序

1. **MILP 求解路径**（branch-and-bound + cuts + heuristics）— 解锁 `Model.optimize` 对 MIP 模型的完整支持。这是从 "仅 LP" 到 "可对标 HiGHS" 的关键步骤。

2. **Presolve 实现** — 对 LP 和 MILP 同时生效的性能提升。依赖 `src/presolve/` 模块的 reduction rule 实现。

3. **MIP 回调接入** — `cbCut`（lazy constraints）、`cbSolution`（heuristic solutions）、`cbGet`（callback-time query）。MILP solver 就位后可以接入。

4. **QP 求解路径** — Active-set 或 barrier QP solver。数据层已完整，仅缺 solver 实现。

5. **一般约束求解** — 将 GenConstr 翻译为 MIP 等价形式（big-M、SOS、PWL 近似）或接入 NLP solver。
