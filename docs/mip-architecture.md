# MIP architecture and performance contract

This document defines the module and ownership boundaries for the zhighs MIP
solver. HiGHS is used as an architectural reference, not as a requirement for
one-to-one source translation.

## Module map

| Responsibility | HiGHS reference area | zhighs module | Contract |
|---|---|---|---|
| Public solve facade | `HighsMipSolver` | `src/mip/root.zig`, `types.zig` | Stable public types and one solve entry point; no search implementation. |
| Search orchestration | `HighsSearch` | `src/mip/search.zig` | Limits, LP dispatch, plugin scheduling, and result publication only. |
| Open nodes | `HighsNodeQueue` | `src/mip/node_queue.zig` | Own queued nodes; default DFS is O(1); plugin selectors receive metadata views. |
| Node state | `HighsNodeQueue`, `HighsDomainChange` | `src/mip/node.zig` | Own sparse bound-change paths, never complete model copies. |
| Domain reconstruction | `HighsDomain`, `HighsDomainChange` | `src/mip/domain.zig` | Normalize the root domain once and replay sparse changes into reusable arrays. |
| Incumbents | `HighsMipSolverData` incumbent handling | `src/mip/incumbent.zig` | Validate integrality, bounds, rows, and objective before ownership transfer. |
| Global bound and gap | MIP primal/dual bound handling | `src/mip/bounds.zig` | Centralize objective-sense comparisons and stopping tolerances. |
| LP relaxation | `HighsLpRelaxation` | `src/solver/root.zig` | Reuse an LP session and publish cold/reoptimization counters. |
| Algorithms | branching and primal heuristics | `src/plugins_builtin/*.zig` | Algorithms implement restricted plugin contracts and do not own solver state. |

## Ownership invariants

1. `ProblemView` and plugin contexts borrow model arrays for the duration of a
   solve or callback. They never free or retain those arrays.
2. `Domain` owns the normalized root bounds, semi-variable active bounds, and
   one mutable reconstruction buffer pair.
3. A solve-local persistent arena owns all `BoundChange` entries. A `Node`
   stores only an optional, strongly typed `ChangeId` into that arena.
4. `NodeQueue` owns every queued node value. A node removed by `take` transfers
   to the caller; change entries remain stable until the solve ends.
5. `SolveResult` owns only its primal array. All search workspaces are released
   before `solve` returns.

## Hot-path rules

- Default DFS must not allocate or materialize plugin views when selecting a
  node. Custom selectors may pay O(open nodes) for the restricted view.
- Branching appends to a parent-linked change arena instead of cloning complete
  bound arrays or ancestor paths. Child creation is amortized O(1).
- LP bound arrays use stable addresses throughout a solve so the LP session can
  distinguish value revisions from structural revisions.
- Diagnostics and statistics are opt-in. Performance builds use `ReleaseFast`
  and benchmark parsing separately from repeated solves.
- No optimization is accepted from elapsed time alone. Status, objective,
  global bound, node count, and LP work must be checked in the same run.

## Performance workflow

Build and run the repeatable benchmark:

```text
zig build bench-mip -Doptimize=ReleaseFast -- model.mps 25
zig build build-bench-mip -Doptimize=ReleaseFast
```

The installed `mip-solver-bench` binary is the sampling-profiler target. Its TSV
output includes separate LP and MIP mean times, objective, bound, nodes, LP
iterations, work, cold solves, and reoptimizations. Comparisons require a pinned
model corpus, identical compiler mode, warmups, multiple samples, and
objective/status parity. Recorded smoke measurements live in
`docs/mip-performance.md`.

## Deliberately incomplete HiGHS-aligned modules

The current solver is branch-and-bound, not a complete HiGHS-class
branch-and-cut implementation. Add the following as separate modules rather
than expanding `search.zig`:

1. `propagation.zig` and reversible domain stack;
2. `pseudocost.zig` plus reliability/strong branching;
3. `cut_pool.zig`, `separation.zig`, and cut aging;
4. `conflict.zig` and implication/clique tables;
5. `symmetry.zig` and restart policy;
6. deterministic parallel node processing;
7. MIP event/callback bridging and structured timing statistics.

Each new algorithm needs a unit test, an end-to-end correctness corpus, and a
before/after benchmark or profiler trace before becoming a default component.
