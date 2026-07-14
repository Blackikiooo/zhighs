# Matrix layout audit (2026-07-14)

Reproduce with:

```bash
zig build audit-matrix-layout -Doptimize=ReleaseFast
zig build audit-matrix-layout -Doptimize=ReleaseFast -Dhighs-int-width=w64
```

Zig ordinary structs use auto layout; source field order is not an ABI or padding
guarantee. Matrix element storage is SoA outside the control structs, so both the
control-block ABI and the allocated bytes per row/column/nonzero must be audited.

| item | w32 | w64 |
| --- | ---: | ---: |
| `usize` | 8 B | 8 B |
| `HUInt`, `RowId`, `ColId` | 4 B | 8 B |
| slice or optional slice header | 16 B | 16 B |
| `CscMatrix` control block | 104 B | 104 B |
| `CscView` / `CsrView` | 64 B | 64 B |
| `CsrBuffers` control block | 80 B | 80 B |
| `TransposeBuffers` control block | 96 B | 96 B |
| builder triplet field streams per entry | 24 B | 32 B |

For a w32 canonical matrix, authoritative entry streams cost 12 bytes/nnz
(`RowId` or `ColId` plus `f64`). Offsets cost 4 or 8 bytes per compressed
dimension entry. The compact CSC present at the time of this audit stores both
4-byte and 8-byte column offsets, so it pays 12 bytes per column offset. This is
material for matrices with many sparse columns and is tracked separately from
the 104-byte `CscMatrix` header.

`CsrBuffers` additionally owns a 4-byte-per-row cursor in w32. That scratch is
useful for repeated conversion but is not part of a read-only CSR view and must
not be counted as authoritative CSR storage in a fair HiGHS comparison.
