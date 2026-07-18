# Matrix 模块性能审计

本文件记录 `src/matrix/` 中所有已识别的性能优化机会。每条发现标注位置、热度评估、预期收益和实验建议。标记 `[建议实验]` 的条目需要 benchmark 验证后决定是否合入。本审计由deepseek-v4-pro完成，需要结合审计内容认真做测试案例和benchmark评测，好的采纳，不好的建议要给出合理的建议并写在audit_answer.md中.

## 热度定义

| 等级   | 含义                                      |
| ------ | ----------------------------------------- |
| 🔴 top | 每次 simplex pivot / FTRAN / BTRAN 都执行 |
| 🟠 mid | 每次 factorization / reinversion 执行一次 |
| 🟡 low | 诊断、范数、统计等非迭代路径              |

---

## 🔴 热路径：pivot 应用

### P1. `applyPivot` 中 fill-in 查找的 generation-mark 条件过于保守

**文件**：[`sparse_kernel.zig:616`](sparse_kernel.zig#L616)

**当前代码**：

```zig
const use_lookup = u_count >= 2;
```

**问题**：generation-mark dense lookup 仅在 U-factor 产生 ≥2 列时才启用。但实际上大多数 pivot 的 `u_count == 1`（U 因子只有一列非零元），此时对每个 L 行执行一次 `find()` → 侵入式链表扫描。对于 `cache_column_maxima` 已启用（即 `nnz/n > 4`，kernel 密度较高）的场景，单次 generation-lookup 的 overhead（一次 `nextLookupGeneration` + 遍历目标行的所有列索引填充 lookup 表）可能比反复 `find()` 更低。关键洞察：`find()` 是 O(row_count) 的线性链表扫描，而 generation-lookup 建表是 O(row_count) 但分摊了所有 u_count 次查询。

**预期收益**：对高密度 kernel（nnz/n > 4）的 pivot throughput 可能有 5-15% 提升。

**实验方案**：改为 `const use_lookup = u_count >= 2 or self.cache_column_maxima;`，用 `sparse_lu_bench` 对比 `factorizeAssumeValid` 时间。风险低——generation-mark 路径已有测试覆盖。

**[建议实验]**

### P2. `find()` 没有利用"最近访问"的局部性

**文件**：[`sparse_kernel.zig:711-722`](sparse_kernel.zig#L711-L722)

**问题**：`find(row, column)` 沿较短链表扫描。在 pivot 应用的内层循环中，`scratch_columns` 通常按 pivot row 的 entry 顺序访问，这意味着同一 row 可能被多次 `find()`。但每次 `find()` 都从 `row_head[row]` 开始全新扫描。一个极低成本（~10 行）的改进：在 `MutableSparseKernel` 中为每行缓存"上次访问的 entry 和 column"，如果下一个查询的 column 恰好更大（因为 scratch_columns 保持 source row 的顺序），则从缓存位置继续而非每次都从 head 开始。

**预期收益**：对 fill-in 密集的 pivot（large l_count × u_count）可能减少链表遍历步数。

**实验方案**：添加 `last_find_entry: []u32` 和 `last_find_column: []u32` 两个数组（均 dimension 大小），在 `find()` 中先检测缓存是否命中或落在后面。注意：pivot 应用期间 entry 会被删除和插入，缓存失效逻辑要配套。中等实现复杂度。

**[建议实验]**

### P3. `ensureFactorCapacity` 在 pivot 循环内逐次检查

**文件**：[`sparse_lu.zig:235`](sparse_lu.zig#L235)

```zig
try self.ensureFactorCapacity(@max(self.l_nonzeros + pivot.l_rows.len, ...));
```

**问题**：每个 pivot 都分支检查容量 + 可能触发 `realloc`。因子的非零元增长完全可预测——不会超过当前 kernel nnz。在 pivot 循环前做一次保守扩容（`factor_capacity >= kernel_nnz * 1.5`），循环内替换为 `std.debug.assert`。

**预期收益**：微小但确定的 CPU 分支预测改善，每次 pivot 省去一个函数调用和条件分支。更容易量化的是：消除了 pivot 循环内任何可能触发 allocator 的代码路径，简化了 failure safety 分析。

**实验方案**：在 `factorizeImpl` 进入 pivot 循环前加入 `try self.ensureFactorCapacity(kernel_nnz + kernel_nnz / 2);`，替换循环内的 `try` 为 `std.debug.assert`。风险极低。

**[建议实验]**

### P4. `removePivotColumnEntry` 在 retire pivot column 时对每个 entry 调 `releaseEntry`

**文件**：[`sparse_kernel.zig:823-831`](sparse_kernel.zig#L823-L831)

**问题**：pivot column 的所有 entry 被逐个 `releaseEntry`——这个 inline 函数写入 freelist。但 retire 完整个 column 后，这些 entry 马上又可能被 fill-in 的 `insert()` 重新分配出来。对于 fill-in 密集的 pivot（inserted_fill ≈ removed_entries），这导致 freelist 的 push/pop ping-pong。更好的策略：批量释放——记录 column 的 entry 链表头尾，一次性拼接到 freelist 头部，避免逐 entry 的 `free_next` 写入。

**预期收益**：对 fill-heavy pivot 减少约 (removed_entries - 1) 次 `free_next` 写操作。

**实验方案**：在 `removePivotColumnEntry` 中仅解除 column/row 链接并减少 count，将 entry 加入一个 thread-local "to-free" 链表，在 pivot 应用完成后批量拼接到 freelist。中等复杂度。

**[建议实验]**

---

## 🔴 热路径：FTRAN / BTRAN solve

### P5. Hyper-sparse solve 在 FT 更新后完全退化为稠密求解

**文件**：[`sparse_lu.zig:402-416`](sparse_lu.zig#L402-L416)

**问题**：已在之前讨论中详述。核心：`solveHyperSparse` 在 `ft.hasUpdates()` 为 true 时执行 `@memset(work, 0, n)` + `@memcpy(rhs, work, n)` + `solve(rhs)` + 全扫描收集非零元——四个 O(n) 操作。而 `solveAdaptive` 的 12.5% density 阈值会将稀疏 RHS 主动导向这个伪装路径。

**短期修复**：在 FT 更新存在时，`solveAdaptive` 直接走 dense path。一行改动：

```zig
if (self.ft_ready and self.ft.hasUpdates()) return if (transpose) self.solveTranspose(rhs) else self.solve(rhs);
```

这避免了 hyper-sparse 的 active list 维护开销 + `@memset`/`@memcpy` + 全扫描收集——而实际上本来就做的是完全稠密求解。net effect：同样的稠密求解但没有额外的内存搬运和扫描。

**长期修复**：FT correction 的依赖关系是显式的——`correction_starts/correction_pivots/correction_indices` 形成了一个修正图。Hyper-sparse FTRAN 可以从 input_indices 出发，沿 FT correction → L factor 做 reachability 传播，只访问被填充的位置。`SparseForrestTomlin` 已有 `correction_indices/correction_values` 存储这个依赖图。

**预期收益**：短期修复对已更新 basis 的每个 pricing pass 节省 O(n) 的内存搬运（约 `2*n*sizeof(f64)` + n 次 branch）。长期修复使 hyper-sparse 在 FT 更新后真正发挥作用。

**实验方案（短期）**：改一行，用 simplex e2e benchmark 在多次 FT 更新后测试 pricing 时间。

**[建议实验]**

### P6. `solveTranspose` 和 `solveTransposeFt` 有重复的 L^T 回代循环

**文件**：[`sparse_lu.zig:342-363`](sparse_lu.zig#L342-L363)、[`sparse_lu.zig:365-379`](sparse_lu.zig#L365-L379)

**问题**：`solveTranspose` 中 U^T forward + L^T backward 共两个循环；`solveTransposeFt` 中 FT upper transpose solve + U^T forward + L^T backward 共三个循环。两个函数的 L^T 回代循环完全相同（逐字重复 ~15 行）。这不是性能 bug，但是 copy-paste 导致 BTRAN 路径的优化需要同步两处。

**建议**：提取 `fn backwardSubstitutionL(self, rhs)` 消除重复。不直接提升性能但降低 BTRAN 路径同步遗漏的风险。

### P7. `solve` 的 scatter 写入忽略 `rhs` 已为零的可能性

**文件**：[`sparse_lu.zig:316`](sparse_lu.zig#L316)

```zig
for (0..n) |position| rhs[self.pivot_columns[position]] = self.work[position];
```

**问题**：标准 FTRAN 在 simplex 中通常用于求解 `B * d = a_q`（进入方向），其中 `a_q` 只有一列非零。此时 work 中大部分位置为零，但上面的 scatter 无条件写入所有 n 个位置。对于有 5,000 行的 basis，这意味着 5,000 次 store——其中可能只有 ~50 次是非零的。

**建议**：利用 `self.active` 列表（已在 hyper-sparse 中维护了 active 位置），仅 scatter 被标记为 active 的位置。这需要 `solve` 在求解过程中维护 active 列表——当前只有 `solveHyperSparse` 做了这个。

**[建议实验]**

---

## 🔴 热路径：pivot 选择

### P8. `choosePivot` (DOD Markowitz) 缺少 early termination by favorable side

**文件**：[`sparse_kernel.zig:340-371`](sparse_kernel.zig#L340-L371) vs [`sparse_kernel.zig:394-474`](sparse_kernel.zig#L394-L474)

**问题**：`choosePivotHighs` 在列桶扫描中找到 `row_count[row] < column_count[column]` 的候选时立即返回（该候选在行侧肯定不如列侧好，不需要再扫行桶）。DOD 后端缺少这条提前终止规则，会扫满 `markowitz_candidate_limit` 个候选。

**预期收益**：对非对称稀疏 pattern 的 kernel，可能减少 ~30% 的候选检查。

**建议**：从 `choosePivotHighs` 移植规则到 DOD 后端。改动约 5 行。

**[建议实验]**

### P9. Symbolic planner 的 `chooseMarkowitz` 每次候选检查重算 column maximum

**文件**：[`sparse_symbolic.zig:261-269`](sparse_symbolic.zig#L261-L269)

```zig
for (begin..end) |entry| {
    const row = basis.rows[entry].toUsize();
    if (self.row_active[row]) maximum = @max(maximum, @abs(basis.values[entry]));
}
```

**问题**：symbolic planner 的 Markowitz 阶段（仅在 `planSingletonPrefix` 返回后调用一次，选中一个 kernel pivot 即返回），然而这个单次扫描中的 column maximum 计算是从 CSC values 逐 entry 读取——OK，因为 planner 不持有数值的预计算缓存。但对于有大量候选列的 dense kernel，这可能扫描很多列。由于 planner 只选一个 kernel pivot 就返回，实际开销很小。**标记为低优先级，仅记录。**

---

## 🟠 中频：factorization / reinversion

### P10. `loadImpl` 在 `validate=true` 时逐 entry 验证

**文件**：[`sparse_kernel.zig:270`](sparse_kernel.zig#L270)

```zig
if (validate and (row >= n or !std.math.isFinite(value) or value == 0.0)) return error.InvalidBasis;
```

**问题**：`validate=true` 时在内层循环逐 entry 检查数值有效性。但 `SparseBasisView` 的数据源（`SparseBasisBuffers`）在构建时已经过验证。checked 入口 `load()` 可以通过在循环前做一次集中验证（如验证 starts non-decreasing + rows sorted + values finite）来消除这个逐 entry 分支。这不仅提升了 `load()` 本身的性能（虽然仅在 factorization 开始时调用一次），更重要的是消除了 `loadImpl` 中 `comptime validate` 参数存在的理由——之后 checked 和 AssumeValid 入口可以共享完全相同的 hot loop。

**实验方案**：提取 `validateBasisView(basis) !void` 独立函数，在 `load()` 开头调用。之后 `loadImpl` 去掉 `comptime validate`，`load` 和 `loadAssumeValid` 共享同一个实现。

**[建议实验]**

### P11. `buildHyperViews` 中的 reverse column/view 构建可以用单次遍历完成

**文件**：[`sparse_lu.zig:627-656`](sparse_lu.zig#L627-L656)

**问题**：`buildHyperViews` 为 U column view 和 L row view 各执行一次 histogram + prefix + scatter（共 6 次数组遍历）。这两个视图的构建可以合并为：一次遍历 L+U factors，同时填充 U column 和 L row 的计数；一次 prefix sum；一次遍历同时 scatter 两者。

**预期收益**：减少 ~2 次 factor 数组遍历。但由于 `buildHyperViews` 只在首次 hyper-sparse solve 或 FT 初始化时调用（每个 factorization 一次），收益不会体现在迭代级 benchmark 上。**低优先级。**

### P12. `factorizeImpl` 中的 pivot trace replay 全扫描

**文件**：[`sparse_lu.zig:204-230`](sparse_lu.zig#L204-L230)

**问题**：在 trace replay 模式中（`factorizeWithTraceRepairAssumeValid`），每个 pivot 调用 `chooseRecordedPivotThreshold` 检查 recorded pivot 是否仍然满足 column 阈值——这需要 `columnMaximum()` 扫描当前列的全部 entry。如果 trace 的前 N 个 pivot 都有效（常见于邻近 basis），前 N 次 column 扫描是浪费的——因为阈值检查可以延迟到第一个失败的 pivot 之后。

**建议**：replay 阶段去掉阈值检查（trace 已来自同一求解过程的已验证 pivot），仅在第一个 `chooseRecordedPivot` 失败进入 repair 后，对 repair 阶段执行阈值检查。这消除了 trace replay 中的列扫描开销。

**[建议实验]**

---

## 🟡 低频：SpMV / 范数 / 统计

### P13. CSC dense SpMV 未使用显式 `@Vector`

**文件**：[`csc.zig:383-401`](csc.zig#L383-L401)

**问题**：`multiplyDenseKernel` 是标量 `@mulAdd` 循环。`absoluteRange`、`scaleColumn`、`columnOneNorms` 均已使用显式 `@Vector(VW, f64)` 并获得 0.25-0.77x 加速。SpMV 的列内循环与 scaling 高度相似（连续 values × scalar multiplier 累加到 y[row]）。对于平均 nnz/col > 4 的矩阵（大多数 SuiteSparse 矩阵满足），显式向量化可能减少 load/store 指令。

与 scaling 的关键区别：SpMV 的 y 写入是 scatter（y[row] += ...），而 scaling 是 contiguous store。Scatter 的 SIMD 需要 gather/scatter 指令（AVX2 `vgatherdpd`/`vscatterdpd`），这些指令的吞吐量通常低于连续 load/store。**因此 SIMD SpMV 的收益空间比 SIMD scaling 小，需要实测确认。**

**实验方案**：在列内 nnz >= lanes 时用 `@Vector` 读取 values[row:row+lanes] 和对应的 row indices，做 SIMD 乘法但 scatter 回 y 时回退到标量（因为 scatter 无收益）。这至少减少了 values 的 load 指令数。

**[建议实验]**

### P14. CSR SpMV 可以用 `@Vector` 水平求和替代标量 `@mulAdd` 链

**文件**：[`csr_view.zig:93-106`](csr_view.zig#L93-L106)

**问题**：CSR SpMV 的内循环是 `sum = @mulAdd(f64, vs[pos], x[ci[pos]], sum)`——FMA 链。对于 row 内 nnz >= 4 的行，用 `@Vector` 做 2/4 路并行累加后水平求和，可以减少 FMA 延迟链的长度。这类似于循环展开但交给 SIMD 单元。

**实验方案**：在 `row_nnz >= 2*lanes` 的行使用 `@Vector` 累加 + 水平求和。需要实测确认延迟链缩短的收益是否超过额外寄存器的 spill 成本。

**[建议实验]**

### P15. `addTransposeProductAssumeValid` 的乘积累加用标量而非常量

**文件**：[`ops.zig:200-208`](ops.zig#L200-L208)

**问题**：与 CSR SpMV 相同的问题——CSC^T 相当于 CSR 模式（axis=row 变为 inner loop），内层 `sum += vs[pos] * x[ri[pos].toUsize()]` 是标量累加，可以用 `@mulAdd` 展开。

**当前状态**：未使用 `@mulAdd`。直接用 `+=` 交给 LLVM 自行决定是否 FMA。对于 LLVM 18+，通常能识别这个模式并生成 FMA，但不保证。建议显式改为 `@mulAdd`。

---

## 🟡 低频：SparseAccumulator

### P16. `SparseAccumulator.clear` 的阈值选择可以自适应

**文件**：[`sparse_sum.zig:110-119`](sparse_sum.zig#L110-L119)

```zig
if (10 * self.active_len < 3 * self.dimension) {
    // 逐个清零 touched entry
} else {
    memory.clearF64(self.dense_values);  // 全量 SIMD 清零
}
```

**问题**：硬编码阈值 `active_len * 10 < dimension * 3`。对于实际使用场景（如 simplex pricing 中 aggregate 20 列到 5000 行），`active_len ≈ 100` 且 `dimension = 5000`，条件为 true，走逐 entry 清零。但逐 entry 清零在 `active_len` 较大的情况下可能比 SIMD 全量清零更慢，因为逐 entry 零值是随机写（cache miss），而全量清零是连续 streaming write。

**实验方案**：对不同的 `active_len/dimension` 比例做 microbenchmark，找到 SIMD 全量清零的实际盈亏平衡点。可能把阈值从 30% 调低到 10% 或更低。

**[建议实验]**

### P17. `SparseAccumulator.init` 的热启动问题：第一次分配后 cold cache

**文件**：[`sparse_sum.zig:56`](sparse_sum.zig#L56)

**问题**：`initWithCapacity` 分配合并存储后用 `memory.clearF64` 清零。但这个清零在第一次使用前发生，如果调用方紧接着在同一函数中调用 `addAssumeValid`，这些 cache line 会被立即写入。对于大 dimension（5000+），清零的 streaming write 可能浪费带宽。

**建议**：提供一个 `initWithCapacityUninitialized` 变体，由调用方在首次 `clear()` 时显式清零。或者至少用 non-volatile clear 而非 volatile，允许编译器消除死写。

---

## 🟡 低频：其他文件

### P18. `transposeIntoAssumeValidCompact` 中的 cursor reuse 可以消除一次 @memset

**文件**：[`transpose.zig`](transpose.zig) 的 `TransposeBuffers` cursor 用法

**问题**：每次 transpose 调用都 `@memset(cursor, 0, num_cols)` 清零 cursor。但 cursor 在 scatter 完成后等于最终的 `compact_starts`——下次 transpose 如果矩阵结构相同，cursor 的初始值可以是上一次的 `compact_starts`。这省去了一次 O(num_cols) 的 memset。不过这需要调用方保证矩阵结构不变，适用场景有限。

### P19. `DenseLU.factorize` 每次重新分配而非复用

**文件**：[`dense_lu.zig:33-79`](dense_lu.zig#L33-L79)

**问题**：`factorize` 每次分配新的 `next_lu`、`next_pivots`、`next_work`，用完后释放旧数组。`factorizeOwned` 和 `refactorizeInPlace` 已经提供了零分配路径。但 `factorize` 作为最通用的入口，在重复调用时的 alloc/free 对开销被 `refactorizeInPlace` 和 `factorizeOwned` 修复。**不是问题——已有优化路径。**

### P20. `scaleColumn` 的 preflight 检查可以合并到应用循环

**文件**：[`scaling.zig`](scaling.zig) 的 `scaleColumn`

**问题**：`scaleColumn` 先 preflight 检查所有乘积 finite（一次遍历），再 apply（第二次遍历）。可以合并为一次遍历：一旦发现无效乘积，回滚已修改的值。这需要保存原始值或使用 generation mark 跟踪已修改位置。**低优先级——scaleColumn 是冷路径。**

---

## 总结表

| ID  | 文件            | 热度 | 描述                               | 预期收益            | 实验难度       |
| --- | --------------- | ---- | ---------------------------------- | ------------------- | -------------- |
| P1  | sparse_kernel   | 🔴   | generation-lookup 阈值降低         | 中 (5-15%)          | 低 (1 行改动)  |
| P2  | sparse_kernel   | 🔴   | find() per-row caching             | 中                  | 中 (~30 行)    |
| P3  | sparse_lu       | 🔴   | ensureFactorCapacity 提到循环外    | 低                  | 低 (~5 行)     |
| P4  | sparse_kernel   | 🔴   | 批量释放 pivot entry               | 低                  | 中             |
| P5  | sparse_lu       | 🔴   | hyper-sparse FT 回退修复           | 高 (消除 O(n) 搬运) | 低 (短期 1 行) |
| P6  | sparse_lu       | 🔴   | L^T 回代重复消除                   | 无 (可维护性)       | 低             |
| P7  | sparse_lu       | 🔴   | FTRAN scatter 仅写 active          | 中                  | 中             |
| P8  | sparse_kernel   | 🔴   | DOD Markowitz early term           | 中 (~30% 候选)      | 低 (~5 行)     |
| P9  | sparse_symbolic | 🟠   | symbolic column max 扫描           | 极低                | —              |
| P10 | sparse_kernel   | 🟠   | validate 提取到循环外              | 低                  | 低             |
| P11 | sparse_lu       | 🟠   | buildHyperViews 合并遍历           | 极低                | 低             |
| P12 | sparse_lu       | 🟠   | trace replay 跳过阈值检查          | 中 (trace heavy)    | 低             |
| P13 | csc             | 🟡   | CSC SpMV @Vector                   | 低-中 (需实测)      | 中             |
| P14 | csr_view        | 🟡   | CSR SpMV @Vector 水平求和          | 低-中 (需实测)      | 中             |
| P15 | ops             | 🟡   | addTransposeProduct @mulAdd        | 低                  | 低             |
| P16 | sparse_sum      | 🟡   | SparseAccumulator clear 自适应阈值 | 低                  | 低             |
| P17 | sparse_sum      | 🟡   | initWithCapacity 跳过预清零        | 低                  | 低             |
| P18 | transpose       | 🟡   | cursor @memset 消除                | 低                  | 低             |
| P19 | dense_lu        | 🟡   | factorize 分配（已有优化路径）     | 无                  | —              |
| P20 | scaling         | 🟡   | scaleColumn preflight 合并         | 低                  | 中             |

20 条发现中，**P1、P5、P8、P12** 四条改动最小（1-5 行）、风险最低、且落在最热的路径上。建议优先实验验证这四条。
