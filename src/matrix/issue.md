## 人工审核时，记录目前实现的问题:
1. [已解决] `matrix/builder.zig` 已改用 Zig 官方 `std.MultiArrayList(Triplet)`。注意 Zig 0.16 的 MultiArrayList stable context sort 暂时是 insertion sort；当前使用 PDQ unstable sort，并以 sequence 字段作为最终排序键，在保持重复项确定性求和的同时避免 O(n²) 退化。
2. [已解决] MatrixBuilder.freezeInternal里面为什么需要把const rows = fields.items(.row) 解构出来，我看代码都是一起用的，是不是不结构也行？
3. [已解决] memory.zig中，有很多重复逻辑的代码，需要利用zig的comptime特性重构，减少重复代码。
