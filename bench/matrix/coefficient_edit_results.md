# Batched coefficient edit benchmark

测量日期：2026-07-14。入口：

```bash
zig build bench-coefficient-edits -Doptimize=ReleaseFast
```

固定数据集为 4096 x 4096 三对角 CSC。`existing_values` 每个 batch 修改 4096
个已有对角非零；`structural_mixed` 每个 batch 交替插入/删除 1024 个远离带宽的
非零，并删除/恢复 1024 个已有对角非零。计时包含 enqueue、排序/合并、矩阵提交
和 pending queue 清理，不包含模型初始构建。

| kernel | batch size | repeats | total | throughput |
| --- | ---: | ---: | ---: | ---: |
| existing values | 4096 | 100 | 11.36 ms | 36.06 M edits/s |
| structural mixed | 2048 | 40 | 6.80 ms | 12.05 M edits/s |

接入第一阶段 `ModelEditPlan` 后的 2026-07-14 复测：

| kernel | batch size | throughput |
| --- | ---: | ---: |
| existing values | 4096 | 35.08 M edits/s |
| structural mixed | 2048 | 12.91 M edits/s |
| scalar last-write-wins | 4096 | 50.02 M edits/s |
| scalar direct path | 1 | 166.94 M edits/s |
| scalar direct path | 4 | 280.11 M edits/s |
| scalar direct path | 8 | 329.35 M edits/s |
| scalar DOD plan | 16 | 112.80 M edits/s |

因此 scalar-only direct-path 阈值固定为 8；更大的 segment 通过 target/sequence
排序消除重复写入。单次测量会受 CPU 频率影响，阈值调整必须重新运行同一入口。

结果用于防止批处理退化回逐项 CSC 复制。绝对耗时依赖 CPU、allocator 和 Zig
版本；比较优化前后时必须使用相同构建和机器，并重复运行以检查波动。
