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

## Allocated bytes and lifetime findings

For w32, the current compact CSC owning representation is approximately
`12 * (num_cols + 1) + 12 * nnz` bytes before at most 63 bytes of alignment
padding per stream: both `usize` and `HUInt` starts are retained, followed by
4-byte row ids and 8-byte values. A conventional single-`HUInt` CSC is
`4 * (num_cols + 1) + 12 * nnz`. On `webbase-1M` (1,000,005 columns and
3,105,536 entries), these are about 46.98 MiB and 39.35 MiB respectively; the
dual offsets therefore cost about 7.63 MiB. The compact offsets are currently
consumed by measured hot kernels, so removing the wide public offsets requires
an API migration rather than struct field reordering.

Representative compact-CSC allocation sizes (excluding the 104-byte control
block) are:

| shape / nnz | naturally packed w32 bytes | observation |
| --- | ---: | --- |
| 0 columns / 0 nnz | 128 | fixed stream alignment dominates; empty-matrix specialization remains possible |
| 16 columns / 32 nnz | 704 | old page-colored builder used 12,736 bytes |
| 50,000 columns / 149,998 nnz | 2,400,112 | padding is only 124 bytes |
| webbase-1M / 3,105,536 nnz | 49,266,560 | only 56 padding bytes; dual starts dominate avoidable bytes |

The layout policy now rejects page coloring for authoritative owning data when
it adds more than one page or 5% without a repeatable end-to-end cycle win.

The corresponding authoritative w32 CSR is
`4 * (num_rows + 1) + 12 * nnz`. `CsrBuffers` adds a reusable
`4 * num_rows` cursor, but `CsrCache` does not publish or retain that cursor.
Likewise, `TransposeBuffers` deliberately combines output and scratch for
caller-owned reuse, while the default owning `transposeAssumeValid` now reuses
its final compact-starts stream as the cursor and restores it after scatter.
For the 50,000 by 50,000 synthetic transpose case this removes about 200 KiB of
persistent cursor plus page padding (roughly 8% of the old returned allocation);
at one million output columns it removes about 3.81 MiB plus padding. Avoiding
the separate cursor allocation also reduced the five-process synthetic median
from about 982 us to 906 us (7.7%), with the same checksum and structural hash.

Page coloring was not beneficial for the default compact builder allocation.
For a 16-column, 32-entry matrix, the old four-stream page layout occupied
12,736 bytes, versus 704 bytes with naturally packed 64-byte stream alignment.
On the synthetic CSC benchmark, natural packing used about 1.4% fewer cycles
(284.5M versus 288.5M across the sampled run), with identical instructions and
no stable cache-miss regression. The default builder therefore uses natural
packing; page-colored reusable workspaces remain separate experimental policy.

## Explicit SIMD and allocator decisions

ReleaseFast assembly showed that the ordinary dynamic inner column loop used
scalar `vmulsd`, while the explicit implementation used packed `vmulpd`.
Manual SIMD reduced column-scaling time to 62-66% of scalar for 8-32 entries per
column and to 77% at 128 entries; three-entry columns were neutral. Explicit
abs/min/max reduction used about 25% of scalar time on 4K-131K values and 36%
at about two million values. The production paths retain unaligned vector
loads/stores and scalar tails. Indirect scatter/gather paths (CSC SpMV,
transpose/CSR scatter and row scaling) were rejected as manual-vector
candidates.

Allocator experiments found no matrix-internal replacement that preserves the
existing independent ownership contract. A 4 KiB stack fallback was 2-70%
slower than `smp_allocator` for 64-4096-element cursor scratch; page allocation
was 13.7-22.6x slower for tiny builds and had no stable larger-build gain. A
retained Arena reduced repeated 512/4096-dimension sorted-build time to about
17-28%/28% of ordinary allocation after warm-up, but its lifetime couples every
returned object to a session reset. It is therefore documented as a caller-selected
compile/presolve-session policy, not hardcoded into matrix APIs. Reusable
`CscBuildBuffers`, `CsrBuffers`, and transpose-into workspaces remain the safe
allocation-free production routes.

The reserved sorted-build experiment makes two logical allocations per
iteration (triplet SoA and owning CSC). After Arena warm-up, retain/reset serves
both without another backing allocation; retained capacities for dimensions
64/512/4096 were about 10.6/81.5/648.5 KiB. Those bytes are intentionally
session-retained and are not reported as an independently owned matrix memory
reduction.

A later libc-linked A/B corrected the C++ fairness boundary for repeated large
owning objects. The identical Zig owning-transpose code measured about 905 us
with `smp_allocator`, 357 us with `c_allocator`, and 340 us in the C++ malloc
runner. Perf counted about 61,867 page faults for smp versus 3,853 for c and
4,038 for C++; retired instructions were unchanged between the Zig allocators.
Therefore libc-linked solver sessions that repeatedly create and destroy large
owning matrices should select `c_allocator` (or an explicitly retained session
allocator). Matrix APIs remain allocator-parametric and do not hardcode libc;
reusable buffers remain preferred when ownership permits.
