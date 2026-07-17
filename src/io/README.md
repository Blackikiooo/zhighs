# zhighs model I/O

The I/O module is independent of the public `Model` type. Format parsers return
an owning `ModelData` with canonical CSC storage; writers consume a borrowed
`ModelView`. `src/model/model_io.zig` is only an API and ownership adapter.

## Dependency and data flow

```text
file + std.Io
    -> suffix dispatch
    -> LP or MPS parser
    -> semantic Builder
    -> canonical ModelData (CSC + row bounds)
    -> model adapter / presolve / compiler
```

Row constraints use lower and upper bounds rather than a sense/RHS pair. This
represents ordinary one-sided rows, equalities, and MPS ranges without loss.

## Performance contracts

- Parsers borrow tokens and names from the input buffer during parsing.
- LP uses an allocation-free, zero-copy streaming lexer. Tokens retain absolute
  byte, line, and column locations, and cheap lexer copies provide checkpoints
  for ambiguous grammar prefixes without materializing token arrays.
- No syntax tree is built for linear expressions.
- LP terms enter per-column chains backed by one compact shared node pool;
  duplicate coordinates merge at append time and freeze directly into one
  exact-size CSC allocation without a global sort or per-column allocations.
- General and unordered format terms use compact strong row/column IDs, are
  stably sorted and duplicate-merged in place, then frozen directly into CSC.
  No secondary matrix-builder triplet copy is created.
- Published column/row attributes share one 64-byte-aligned packed allocation;
  the model and all retained row/column names share one contiguous string pool.
- LP/MPS name-index hash tables are released as soon as syntax resolution
  finishes, before CSC and published attribute allocations. With
  `keep_names=false`, no row or column name bytes enter the final string pool,
  and the per-row/per-column nullable name-reference tables are omitted.
- MPS attempts a direct ordered-column CSC freeze and falls back to the general
  deterministic duplicate-merging path for unordered input.
- Writers use a 128 KiB buffered `std.Io.Writer`; numeric formatting does not
  allocate per token.
- `ReadOptions` enforces file, physical-line, token, row, column, temporary
  matrix-term, and cumulative object-name-byte limits before parser storage is
  grown. Limit failures are reported as `FileTooLarge` or
  `ResourceLimitExceeded`, never as a misleading allocation failure.
- An optional caller-owned `std.atomic.Value(bool)` interrupt flag provides
  cooperative cancellation. Parsers and Builder finalization poll it at a
  configurable work interval and return `Cancelled` with normal cleanup.
- `ReadOptions.input_mode` selects automatic, buffered, or read-only
  memory-mapped input. Automatic mode maps large files without prefaulting and
  keeps small files on the lower-overhead buffered path.
- The injected `std.Io` backend keeps filesystem policy outside the parser and
  leaves room for async I/O without duplicating grammar code.

## Current format coverage

LP supports linear objectives, linear and ranged rows, bounds, binary, integer,
semi-continuous and semi-integer variables, multiline objectives, comments,
and generated or explicit names. The accepted grammar and lexer contracts are
documented in [`../../docs/lp-format.md`](../../docs/lp-format.md).

MPS supports the common free/fixed whitespace-compatible representation,
`OBJSENSE`, `ROWS`, `COLUMNS`, integer markers, `RHS`, `RANGES`, and the common
`BOUNDS` records. Multiple RHS, range, or bound vectors use the first vector,
matching the usual single-model import behavior.

Compressed files, quadratic expressions, SOS, general constraints, multiple
objectives, and fixed-field names containing embedded spaces are intentionally
reported as unsupported rather than silently changed.

## Next performance steps

1. Add chunked/windowed input for streams where whole-file mmap is unavailable.
2. Construct public `Model` storage directly from `ModelData`, removing the
   current pending-change copy in the compatibility adapter.
3. Pipeline decompression, scanning, and CSC finalization where profiling shows
   useful overlap; keep deterministic reduction and diagnostics.
