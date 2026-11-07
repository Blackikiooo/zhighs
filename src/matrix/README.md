# Sparse matrix module

Every file has exactly one job. This page is the map: what each module does,
what shape it produces, and which module to reach for.

## Ownership rule

`CscMatrix` is the only authoritative representation. Everything else is either
a builder (produces a `CscMatrix`), a borrowed view (`CscView`, `CsrView`,
`SparseVectorView`), a disposable cache (`CsrCache`), or a caller-owned buffer
(`CscTransformBuffers`, `CsrBuffers`, `TransposeBuffers`).

## Shape cheat sheet

| Module | Operation | Shape |
|---|---|---|
| `transpose.zig` | transpose | **m×n → n×m** |
| `slice.zig` | extract columns | m×n → **m×k** |
| `slice.zig` | extract rows | m×n → **k×n** |
| `slice.zig` | extract column range | m×n → **m×(to−from)** |
| `edit.zig` | append columns | m×n → **m×(n+k)** |
| `edit.zig` | append rows | m×n → **(m+k)×n** |
| `edit.zig` | delete columns | m×n → **m×(n−k)** |
| `edit.zig` | delete rows | m×n → **(m−k)×n** |
| `permutation.zig` | reorder rows/columns | m×n → **m×n** (indices relabelled) |
| `scaling.zig` | scale rows/columns | m×n → **m×n** (values only, in place) |
| `dynamic_rows.zig` | buffer then merge rows | 0×n → k×n → **(m+k)×n** |
| `csr_view.zig` | derive row-major view | m×n CSC → **m×n CSR** |
| `ops.zig` | statistics / products | m×n → **scalar or vectors** |
| `builder.zig` | canonicalise coordinates | triplets → **m×n** |
| `store.zig` | own + version the matrix | **m×n** (no shape change) |
| `sparse_vector.zig` | one sparse row/column | **1×n** or **m×1** payload |
| `sparse_sum.zig` | accumulate a sparse vector | **1×n** or **m×1** scratch |

## Who owns what

| Module | Returns | Mutates caller data |
|---|---|---|
| `builder.zig` | owning `CscMatrix` | no |
| `csc.zig` | owning `CscMatrix` / borrowed `CscView` | no |
| `ops.zig` | numbers | no |
| `transpose.zig` | owning `CscMatrix` / borrowed `CscView` | no |
| `slice.zig` | owning `CscMatrix` / borrowed `CscView` | no |
| `edit.zig` | owning `CscMatrix` / borrowed `CscView` | no |
| `permutation.zig` | owning `CscMatrix` / borrowed `CscView` | no |
| `scaling.zig` | factors or `void` | **yes, values in place** |
| `dynamic_rows.zig` | owning `CscMatrix` via `appendToCsc` | no |
| `csr_view.zig` | `CsrCache` / `CsrView` | no |
| `store.zig` | `*const CscMatrix` / `CsrView` | **yes, via its mutators** |

## Data flow

```text
triplets (i, j, v)
    └─ builder.zig ──freeze──► csc.zig (CscMatrix)
                                   │
              ┌────────────────────┼─────────────────────┐
              ▼                    ▼                     ▼
        csr_view.zig          ops.zig              store.zig
        (row access)          (read-only)          (authority + revision)
              │                                          │
              ▼                                          ▼
         CsrCache                       edit / slice / permutation /
                                        transpose / scaling  ──► new CSC
                                        dynamic_rows ──appendToCsc──►
```

## Choosing a module

| I want to... | Use |
|---|---|
| build from arbitrary (row, col, value) | `builder.zig` |
| build from already sorted arrays | `builder.freezeFromSortedArrays*` |
| read rows repeatedly | `csr_view.zig` (`CsrCache`) |
| compute norms or products | `ops.zig` |
| transpose | `transpose.zig` |
| take a submatrix | `slice.zig` |
| add or remove rows/columns | `edit.zig` |
| relabel indices | `permutation.zig` |
| rescale coefficients | `scaling.zig` |
| accumulate MIP cuts or presolve rows | `dynamic_rows.zig` |
| own the matrix a model solves | `store.zig` |

## Design conclusions

These are the decisions this module settled on. Treat them as the review
checklist for new code in this directory.

### 1. CSC is the only authoritative representation

- `CscMatrix` owns the canonical streams; `CsrCache` is a disposable projection
  keyed by the owner's revision.
- A second owning matrix type (say `CsrMatrix`) would create two sources of
  truth and force a double write on every edit.
- Consequence: row-oriented consumers take a `CsrView`, never a mutable CSR.

### 2. Ownership is expressed by type, not by flags

- `CscMatrix` (owning) vs `CscView` (borrowed); `CsrCache` (owning cache) vs
  `CsrView` (borrowed); `SparseVector` vs `SparseVectorView`.
- `CscMatrix.storage` is a tagged union (`separate` / `single_block` /
  `borrowed`), so "single block that also owns three slices" cannot be
  represented at all.
- Do not add `owns_*` booleans, `_`-prefixed fields, or a stored allocator:
  Zig has no private fields, and none of those change what `deinit` must free.

### 3. Builders accept disorder; transforms are pure

- `builder.zig` is the only module that accepts arbitrary `(row, col, value)`
  order; it sorts, merges duplicates, and drops numerical zeros at `freeze`.
- `transpose` / `slice` / `edit` / `permutation` never mutate their input; they
  return a new matrix, so a failed allocation leaves the old revision intact.
- `slice.zig` can only shrink; `edit.zig` can grow or shrink.
- `scaling.zig` and `MatrixStore`'s mutators are the only in-place writers.

### 4. Row-at-a-time growth needs its own buffer

- `dynamic_rows.zig` exists because MIP cuts and presolve rows arrive one at a
  time; `edit.appendRows` would rebuild the whole CSC per row.
- It is append-only with O(1) checkpoint/rollback; the merge happens once in
  `appendToCsc`, which delegates to `edit.appendRowsFromCsrAssumeValid`.

### 5. Reads take pointers, writes stay explicit

- `ops.zig` takes `*const CscMatrix`: the struct is a 104-byte header, so a
  by-value call spills the whole header at every call site. The pointer also
  matches `MatrixStore.csc()`, which already returns `*const CscMatrix`.
- `*const` does not make coefficients immutable (`values` is still `[]f64`).
  Type-level read-only requires `CscView` (`[]const f64`).

### 6. Layout is part of the contract

- Owning blocks are one 64-byte-aligned allocation on x86 and 128-byte-aligned
  on Apple Silicon / POWER (`memory.owning_block_alignment`).
- Multi-stream blocks use page coloring with cache-line-sized in-page offsets
  (`memory.lineOffsets`), so two streams never share a cache line.
- `target_policy.zig` derives vector widths from
  `std.simd.suggestVectorLength`, so `-Dcpu=...` and cross-compilation get the
  right width without source changes.

### 7. Known gaps

- `permutation.zig` has no production caller yet; only the matrix dataset
  runner round-trips it. Keep it, but do not assume it is on a hot path.
- `ops.zig` and the other `matrix/*.zig` unit tests are not collected by
  `zig build test` (only `matrix/root.zig`'s single test runs), so matrix-level
  coverage currently comes from the higher-level suites.
- `permuteAssumeValid` always re-sorts every destination column, even for
  column-only permutations; optimize only if it becomes hot.
- `transposeIntoAssumeValidCompact` and `permuteInto*` share
  `CscTransformBuffers`; a view returned by one is invalidated by the next call
  that writes the same buffers.

## Notes

- `permutation.zig` only relabels indices; `transpose.zig` changes the shape.
- `dynamic_rows.zig` batches row additions; `edit.appendRows` rebuilds once per call.
- `csr_view.zig` never owns the matrix; CSR is derived from CSC and keyed by the
  owner's revision.

