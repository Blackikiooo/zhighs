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

## 4. Phase I / Phase II 与状态输出

- [x] Artificial variables 和 Phase-I objective。
- [x] Infeasible/unbounded 判定。
- [x] Phase-I cleanup 与 Phase-II objective transition。
- [x] Basis status、primal、dual、reduced cost 和模型 solution 回写。

## 5. 端到端验收

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
- [ ] Netlib 和 Mittelmann 完整求解结果与 HiGHS/CLP 对比。
- [ ] 汇总 objective、status、iteration count 和 residual。
- [ ] 汇总 reinversion 次数、原因和 FT chain 长度。
- [ ] 比较完整 solve 时间，而不只是 INVERT/FTRAN/BTRAN。
