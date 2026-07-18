# I/O 模块实现状态

状态说明：

- `[x]` 已实现并通过当前测试。
- `[ ]` 尚未实现或仍需补充验证。
- `[*]` 暂不实施 — 优先级低或无紧迫需求（标注原因）。

## 1. 基础设施

- [x] `types.zig` — IoError 错误集、Format/Compression/FileKind 枚举、ObjectiveSense/RowSense/VariableType、Diagnostic、InputMode、ReadOptions（文件/行/token/行列/name-byte 等资源限制 + 原子取消标志）、WriteOptions、ModelView、ParseControl。
- [x] `format.zig` — 文件后缀检测（.lp/.rlp/.mps/.rew/.dua/.dlp/.ilp/.opb + .gz/.bz2/.zip/.7z/.xz）。
- [x] `input.zig` — FileInput（mmap/buffered/empty 三态，自动模式大文件 mmap 小文件 buffered，取消检查）。
- [x] `output.zig` — Names（确定性名称 + 重复检测）、print/write 包装。
- [x] `string_arena.zig` — 非 owning 连续字符串池，ModelData 的单次 name interning。
- [x] `builder.zig` — 共享语义 Builder，LP 列链模式（addColumnTerm 合并追加 + finishColumnTerms 直接 CSC），通用模式（addTerm 排序合并 + finishTerms），有序 MPS 快速路径（finishColumnOrdered），资源限制强制。
- [x] `model_data.zig` — 64 字节对齐的紧凑 ModelData（属性 + 名称池 + owning CSC），keep_names=false 紧凑模式。
- [x] `root.zig` — readFile（format detect → FileInput → parse → ModelData）、writeFile（ModelView → format writer）、LP/MPS round-trip 测试。

## 2. LP 解析器

- [x] 零复制 streaming lexer（`lexer.zig`）— 廉价 checkpoint、Location 追踪（byte_offset/line/column）、newline token 作为语法断点、注释跳过。
- [x] Token 定义（`token.zig`）— 21 种 tag（identifier/number/operator/relation/punctuation/stream boundary）。
- [x] Section header 识别 — Minimize/Maximize、Subject To/Such That/ST/S.T.、Bounds/Binaries/Generals/Semi-Continuous/Semi-Integer/End，以及 Semi 的二义消解（probe 后缀）。
- [x] 目标函数解析 — 线性项 + 系数 + 符号翻转，跨行续行（pending_objective_coefficient），常量项。
- [x] 约束解析 — 带/不带 label、双边范围（a <= expr <= b）、等式/不等式、常量项归一化到 RHS。
- [x] Bounds 解析 — 简单（x <= 4 / x >= 0）、free、fixed、双边范围（0 <= x <= 1）、±inf。
- [x] 变量类型 section — Binary/Integer/Semi-Continuous/Semi-Integer，直接修改列类型。
- [x] 列链直通 CSC 快速路径 — 行有序合并到列链 node pool，无全局排序，无逐列分配，最后精确 CSC 分配。
- [x] 资源限制逐层强制 — FileTooLarge、行长限制、token 长度限制、行列上限、matrix term 上限、name byte 上限。
- [x] 协作取消 — 基于 atomic flag 的 ParseControl。
- [x] LP writer — 紧凑输出、双边范围约束、变量类型 section、deterministic name 生成。

- [ ] 约束续行 — 当前约束必须位于一个物理行内。大型生成模型经常跨多行写一个线性表达式，这是线性 LP 兼容性的首要缺口。
- [ ] 二次表达式 — Lexer 已预留 `*` `^` `[` `]` `(` `)`，但 parser 尚未处理 `[ x^2 + x*y ] / 2` 和 `x^2 + y^2 <= 10`。
- [ ] 完整名称规则验证 — 特殊字符、关键字冲突、超长名称、数字开头、quoted/escaped names、不同求解器字符差异（. $ _ 等）、大小写语义。
- [ ] 高级 LP section — SOS、indicator constraints、general constraints（MIN/MAX/ABS/EXP/LOG/SIN/COS 等）、piecewise-linear、multiple objectives、lazy constraints/user cuts、hierarchy/cone sections。
- [ ] 结构化诊断 — Token 已有 line/column/byte_offset，但错误目前主要是 `error.InvalidSyntax`/`error.InvalidNumber`。需补充：文件名、错误 token、所在 section、期望 token、附近源码片段。
- [ ] 健壮性验证 — fuzz testing、malformed LP corpus、超长行/token、极值数值、大量重复坐标、与 Gurobi/CPLEX/HiGHS 差分测试、LP 写出后由其他求解器重新读取的 round-trip 测试。
- [ ] LP writer 峰值内存优化 — 当前按行输出 CSC 构造 O(nnz) row_columns/row_values 临时数组，写超大模型时可优化。

## 3. MPS 解析器

- [x] Free + fixed 字段格式兼容 — 空白符分割字段解析。
- [x] Section keyword 识别 — NAME/OBJSENSE/OBJNAME/ROWS/COLUMNS/RHS/RANGES/BOUNDS/ENDATA。
- [x] ROWS section — N/E/L/G 行类型，支持 N 行作为 objective row。
- [x] COLUMNS section — 列名解析、row-value 对、整数 marker（'MARKER' 'INTORG'/'INTEND'）。
- [x] RHS section — 多 set 支持（取第一个）、objective constant（通过 N 行名识别）、等式行 RHS 正确设置双边。
- [x] RANGES section — 多 set 支持（取第一个）、正确将 range 加到对应行 bounds。
- [x] BOUNDS section — LO/UP/FX/FR/MI/PL/BV/LI/UI/SC/SI 全部 MPS bound type。
- [x] 列有序 CSC 快速路径 — finishColumnOrdered：先检查 terms 是否有序，若是则直接建立 CSC 而无需排序；否则回退到排序路径。
- [x] 语义 name 表释放 — ROWS/COLUMNS 完成后释放 hash table，不重叠最终 CSC/属性分配。
- [x] MPS writer — OBJSENSE/ROWS/COLUMNS（含 INTORG/INTEND marker）/RHS（含 objective offset）/RANGES/BOUNDS 完整输出。
- [x] Section name 在数据 section 中仍为有效 set name（如 `RHS ROW VALUE` 使用 "RHS" 作为 set 名）。

- [ ] OBJNAME 支持 — header 已识别，但未解析多目标名称。标准 MIP 多目标场景需要。
- [ ] 二次 MPS 扩展 — QUADOBJ / QSECTION / QCMATRIX / QMATRIX section 未实现。
- [ ] SOS section — `SETS`/`SOS` section（SOS type 1/2）。
- [ ] 高级 bound type — `SC`/`SI`（semi-continuous/semi-integer upper bound）当前 writer 支持，parser 也支持。

## 4. 格式覆盖缺口

当前仅支持 `.lp` 和 `.mps` 的实际解析和写入。以下后缀在 `format.zig` 中识别但无实现：

- [ ] `.rlp` — reduced/reformatted LP（同 LP 语法，可复用 LP parser）。
- [ ] `.rew` — reweighted MPS。
- [ ] `.dua` — DUA 格式。
- [ ] `.dlp` — dense LP 格式。
- [ ] `.ilp` — indexed LP 格式。
- [ ] `.opb` — OPB（pseudo-Boolean）格式。

## 5. 跨切面能力

- [ ] 压缩输入 — suffix 检测已支持 `.gz/.bz2/.zip/.7z/.xz`，但 `readFile` 遇到压缩文件直接返回 `UnsupportedCompression`。需接入流式解压。
- [ ] 流式/chunked 输入 — 当前解析器需要完整的内存缓冲区。对于超大文件或流式源（如管道），需要窗口化输入。
- [ ] 直接从 ModelData 构造 public Model — 当前 adapter 需经过 pending-change copy，去除这层拷贝可减少导入时的分配峰值。
- [ ] 与 public Model 的 Hessian/Quadratic 数据通道 — ModelData 当前无 Q 项存储，LP/MPS 解析器无法返回二次模型。
- [ ] 解析性能 benchmark — 当前 parser_bench 仅覆盖 LP 路径，缺少 MPS 和大型 corpus benchmark。

## 6. 当前状态总览

```
组件             状态     代码量      说明
─────────────────────────────────────────────────────────────
types.zig         ✅        186      错误集、枚举、ReadOptions、ModelView
format.zig        ✅         51      后缀检测（含压缩后缀枚举）
input.zig         ✅        176      mmap/buffered FileInput
output.zig        ✅         66      Names + 确定性名称生成
string_arena.zig  ✅         70      非 owning 字符串池
builder.zig       ✅        593      列链+通用模式、合并排序、CSC冻结
model_data.zig    ✅        152      对齐紧凑存储、name pool
root.zig          ✅        166      readFile/writeFile + round-trip
lp/lexer.zig      ✅        257      零复制 lexer + Location
lp/token.zig      ✅         43      Token/Tag/Location
lp/root.zig       ✅        620      完整 LP parser + writer
mps/root.zig      ✅        492      完整 MPS parser + writer
─────────────────────────────────────────────────────────────
✅  生产就绪  ⬜  未实现
```

## 下一步执行优先级

1. **【当前最优先】** 约束续行 — 线性 LP 兼容性的首要缺口，大型生成模型的多行约束无法导入。
2. 结构化诊断 — 文件名、行号、列号、错误 token、期望 token、源码片段。
3. 健壮性验证 — fuzz testing、malformed LP corpus、与 Gurobi/CPLEX/HiGHS 差分测试。
4. 压缩输入支持 — gz/bz2 流式解压，大量公共数据集以 .mps.gz 形式分发。
5. 二次表达式 — LP 文件中的 `[ x^2 + x*y ] / 2` 和 `x^2 + y^2 <= 10`，需 Hessian builder 接入。
6. 高级 LP 结构 — SOS、indicator/general constraints、piecewise-linear（长期）。
