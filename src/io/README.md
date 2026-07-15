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
- No syntax tree is built for linear expressions.
- MPS attempts a direct ordered-column CSC freeze and falls back to the general
  deterministic duplicate-merging path for unordered input.
- Writers use a 128 KiB buffered `std.Io.Writer`; numeric formatting does not
  allocate per token.
- `ReadOptions.max_file_bytes` enforces an explicit input resource limit.
- The injected `std.Io` backend keeps filesystem policy outside the parser and
  leaves room for async I/O without duplicating grammar code.

## Current format coverage

LP supports linear objectives, linear and ranged rows, bounds, binary, integer,
semi-continuous and semi-integer variables, multiline objectives, comments,
and generated or explicit names.

MPS supports the common free/fixed whitespace-compatible representation,
`OBJSENSE`, `ROWS`, `COLUMNS`, integer markers, `RHS`, `RANGES`, and the common
`BOUNDS` records. Multiple RHS, range, or bound vectors use the first vector,
matching the usual single-model import behavior.

Compressed files, quadratic expressions, SOS, general constraints, multiple
objectives, and fixed-field names containing embedded spaces are intentionally
reported as unsupported rather than silently changed.

## Next performance steps

1. Replace whole-file buffering with windowed/memory-mapped input for very
   large uncompressed files.
2. Add a pooled name store so source text can be released before model publish.
3. Construct public `Model` storage directly from `ModelData`, removing the
   current pending-change copy in the compatibility adapter.
4. Pipeline decompression, scanning, and CSC finalization where profiling shows
   useful overlap; keep deterministic reduction and diagnostics.
