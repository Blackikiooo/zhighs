# LP 文件语法与解析约定

本文描述 `src/io/lp/` 当前实际接受的 LP 文本语法，以及 lexer 与 parser 之间的边界。
它面向维护者，不试图复刻某个商业求解器的完整格式手册。新增语法时，应同时更新本文、
源码测试和与 HiGHS 的差分测试。

## 1. 解析流程

```text
文件字节
  -> 单一零拷贝 Lexer（文件首部到 EOF）
  -> 行首 section/header 识别
  -> 流式线性表达式 parser
  -> semantic Builder
  -> 规范化 CSC ModelData
```

lexer 不分配 token 字符串：`Token.lexeme` 直接借用输入文件的字节切片。parser 也不构建
AST，而是在读取项时直接写入 builder。输入缓冲区必须至少存活到 `parse` 返回；最终的
`ModelData` 独立拥有发布后的模型数据。

`Lexer` 是一个只包含输入切片和游标的小值类型。复制 lexer 就能建立 checkpoint，适合
探测约束中 `lower <= expression <= upper` 之类的歧义前缀，不需要 token 数组或回滚缓存。

解析资源由 `ReadOptions` 限制，包括文件字节、物理行字节、单 token 字节、行数、列数、
临时矩阵项数和累计对象名称字节。检查发生在对应容器扩容之前。调用方还可以提供一个
线程安全的原子中断标志；lexer/parser 和 CSC 构建阶段按配置间隔协作式轮询，取消时返回
`Cancelled` 并清理所有临时存储。默认单行限制为 16 MiB，单 token 限制为 1 MiB。

## 2. 词法规则

下列规则使用近似 EBNF 表示，`?` 表示可选，`*` 表示重复。

```ebnf
digit       = "0" ... "9" ;
digits      = digit, { digit } ;
fraction    = ".", { digit } ;
exponent    = ("e" | "E"), ("+" | "-")?, digits ;
number      = (digits, fraction? | ".", digits), exponent? ;

relation    = "<=" | ">=" | "=" ;
punctuation = "+" | "-" | "*" | "^" | ":" | ","
            | "(" | ")" | "[" | "]" ;
comment     = "\\", { any-byte-except-newline } ;
```

具体约定：

- 空格、制表符和 `\r` 被跳过；`\n` 会产生 `newline` token。
- 反斜杠 `\` 从当前位置注释到行尾，换行本身不属于注释。
- 正负号是独立 token，不属于数字 token，例如 `-1.2e-3` 被拆成 `minus` 和
  `number("1.2e-3")`。
- 裸 `<` 和 `>` 非法，必须写成 `<=` 或 `>=`。
- 关键字首先统一识别为 identifier，由 parser 根据 section 和上下文判断含义。
- identifier 持续到空白、注释符或运算符分隔符。名称还会经过 parser 的额外校验。
- 每个 token 保存绝对 `byte_offset`，以及从 1 开始的 `line` 和 `column`。

lexer 会识别 `*`、`^`、括号和方括号，是为了给二次项、分段线性表达式等后续语法保留
明确边界；识别 token 不代表当前 parser 已支持相应模型能力。

## 3. 文件结构

推荐按照以下顺序组织文件：

```ebnf
lp-file          = objective-section,
                   constraint-section?,
                   bounds-section?,
                   type-section*,
                   end-header ;

objective-section  = objective-header, objective-line+ ;
constraint-section = constraint-header, constraint-line+ ;
bounds-section     = bounds-header, bound-line+ ;
type-section       = type-header, identifier+ ;
```

section header 不区分大小写，当前别名如下：

| Section | 可接受 header |
|---|---|
| 最小化目标 | `Minimize`, `Minimum`, `Min` |
| 最大化目标 | `Maximize`, `Maximum`, `Max` |
| 约束 | `Subject To`, `Such That`, `ST`, `S.T.` |
| 边界 | `Bounds` |
| 0-1 变量 | `Binary`, `Binaries`, `Bin` |
| 整数变量 | `General`, `Generals`, `Gen` |
| 半连续变量 | `Semi-Continuous`, `Semis`, `Semi` |
| 半整数变量 | `Semi-Integer` |
| 文件结束 | `End` |

目标方向 header 和 `End` 是必需的。空行和注释行可以出现在 section 之间。parser 使用
贯穿整个输入的单一 lexer 流，并只在行首通过 token checkpoint 识别 header；header
仍必须单独占一行。

## 4. 目标函数

```ebnf
objective-line = (name, ":")?, linear-expression ;
linear-expression = signed-term, { ("+" | "-"), term } ;
signed-term     = ("+" | "-")?, term ;
term            = number?, identifier | number ;
```

系数和变量之间可以有空白，也可以直接相邻：`2 x` 与 `2x` 等价。目标函数允许跨多行，
并允许常数偏移量。

```lp
Minimize
 obj: 2x - 3 y
      + z + 5
```

上例产生列费用 `(2, -3, 1)`，目标常数为 `5`。为了兼容历史的多行系数写法，某行末尾
单独出现的数字会暂存：若下一目标行以变量名开始，则数字是该变量的系数；若随后进入
其他 section，则数字作为目标常数提交。生成器会把常数与变量写在同一行，以避免歧义。

## 5. 线性约束

```ebnf
constraint-line = (name, ":")?,
                  (number, "<=")?,
                  linear-expression,
                  ("<=" | ">=" | "="),
                  number ;
```

支持单边、不等式、等式和双边约束：

```lp
Subject To
 capacity: 2x + y <= 10
 demand:   x - y >= -3
 balance:  x + y = 4
 ranged:   1 <= x + 2y + 3 <= 9
```

表达式中的常数会在构建模型时移到边界一侧。例如最后一行规范化为
`-2 <= x + 2y <= 6`。内部统一使用 `row_lower` 和 `row_upper`，不保留 sense/RHS 二元表示。

当前每条约束必须位于同一个物理行。全文件 lexer 已保留换行 token，后续可以在语法层
增加约束续行状态，而不修改 token API 或重新切分输入缓冲区。

## 6. 变量边界

```ebnf
bound-line = identifier, "free"
           | identifier, "=", scalar
           | identifier, "<=", scalar
           | identifier, ">=", scalar
           | scalar, "<=", identifier, "<=", scalar ;

scalar     = ("+" | "-")?, (number | "inf" | "infinity") ;
```

`inf` 和 `infinity` 不区分大小写，只能在 Bounds section 使用。

```lp
Bounds
 x free
 y = 2
 z >= -1
 w <= 10
 -inf <= q <= 4
```

如果没有 Bounds 记录，builder 的默认边界语义保持不变。互相矛盾的双边界会返回
`InvalidBounds`。

## 7. 变量类型

类型 section 中只包含由空白分隔的变量名，可以一行写一个或多个：

```lp
Binaries
 x y
Generals
 count
Semi-Continuous
 flow
Semi-Integer
 batch
```

变量第一次出现在目标、约束、边界或类型 section 时创建；后续引用通过名称哈希表复用
同一列。

## 8. 名称规则

当前名称校验采用以下约束：

- 长度为 1 到 255 字节；
- 首字节不能是数字；
- 不能含空白或 `+ - * ^ : \\ < > = ( ) [ ] ,`；
- 不能与 `min`、`maximize`、`bounds`、`binary`、`general`、`end`、`free` 等保留字
  发生不区分大小写的冲突。

lexer 按字节计算列位置，语法面向 ASCII。非 ASCII 名称当前不会在 lexer 层转码，但也
尚未承诺完整 Unicode 标识符语义。

## 9. 当前不支持的语法

以下内容目前会返回 `UnsupportedFeature` 或相应格式错误，而不会静默降级：

- 二次目标和二次约束，包括方括号二次块、`*` 和 `^` 表达式；
- SOS、indicator、general constraints 和分段线性约束；
- 多目标；
- 约束跨物理行；
- 压缩 LP 文件；
- 名称中嵌入空格。

扩展这些能力时，应优先增加新的 parser 状态或语义事件，不要让 lexer 直接依赖模型
builder。lexer 只负责稳定、快速地划分字节流。

## 10. 性能与测试约束

- lexer 热路径不得为 token 分配内存；
- `peek` 和歧义探测使用值复制 checkpoint；
- parser 直接产生 builder 事件，不建立通用 AST；
- 新语法必须增加最小单元测试、完整模型解析测试和错误输入测试；
- 性能改动使用 `bench/io/` 与同功能 HiGHS reader 比较，必须同时校验行数、列数、
  非零元和数值 checksum，不能只比较耗时。

当前实现入口：

- `src/io/lp/token.zig`：Token、Tag 和源码位置；
- `src/io/lp/lexer.zig`：零拷贝 lexer；
- `src/io/lp/root.zig`：section、语法和语义构建；
- `bench/io/parser_bench.zig`：端到端读取与解析基准。
