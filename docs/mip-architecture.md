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
| Domain propagation | `HighsDomain::propagate` | `src/mip/propagation.zig` | Tighten linear-row bounds before LP solves and persist deductions for descendants. |
| Binary conflict graph | `HighsCliqueTable` | `src/mip/clique_table.zig` | Encode `x` and `~x` as typed vertices, extract sparse conflicts from both row sides, and propagate implications transitively. |
| Reduced-cost fixing | `HighsRedcostFixing` | `src/mip/reduced_cost_fixing.zig` | Convert incumbent cutoff gaps and LP reduced costs into persistent node-local integer bound tightenings. |
| Root reduced-cost cache | `HighsRedcostFixing` lurking bounds | `src/mip/root_reduced_costs.zig` | Retain the root dual certificate and strengthen the shared root domain when later incumbents improve the cutoff. |
| Learned conflicts | `HighsConflictPool` | `src/mip/conflict_pool.zig` | Canonicalize infeasible integer paths, remove root-implied conditions, deduplicate, age, and propagate through column occurrence worklists. |
| Pseudocost state | `HighsPseudocost` | `src/mip/pseudocost.zig` | Learn directional unit objective gains from child LPs and prepare read-only plugin scores. |
| Reliability branching | `HighsSearch::selectBranchingCandidate` | `src/mip/strong_branching.zig` | Probe unreliable directions in an isolated LP session, enforce budgets, and cache child feasibility/bounds. |
| Cut storage and aging | `HighsCutPool` | `src/mip/cut_pool.zig` | Own canonical sparse global cuts, deduplicate them, and age inactive rows. |
| Root separation | `HighsSeparator` family | `src/mip/separation.zig`, `cut_generation.zig` | Generate valid cover cuts and coordinate rounds without putting algorithms in the search loop. |
| Augmented relaxation | LP cut-row synchronization | `src/mip/cut_relaxation.zig` | Materialize active cuts only when the pool revision changes while retaining the immutable base matrix. |
| Incumbents | `HighsMipSolverData` incumbent handling | `src/mip/incumbent.zig` | Validate integrality, bounds, rows, and objective before ownership transfer. |
| Global bound and gap | MIP primal/dual bound handling | `src/mip/bounds.zig` | Centralize objective-sense comparisons and stopping tolerances. |
| LP relaxation | `HighsLpRelaxation` | `src/mip/relaxation.zig`, `src/solver/root.zig` | Keep MIP node policy separate from the reusable LP session and publish cold/reoptimization counters. |
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
6. `Pseudocost` owns directional statistics and its score buffer. Nodes retain
   only the compact branch origin required to submit a later LP observation;
   plugins borrow the prepared score slice for one selection call.
7. `StrongBranching` owns a separate LP session, reusable bound buffers, and
   generation-stamped probe results. It never mutates the main relaxation
   session; its iterations and work still count against global limits.
8. `CutPool` owns cut coefficients, metadata, and membership revisions.
   `CutRelaxation` owns its augmented CSC matrix and row bounds; the LP session
   borrows those arrays only until the next synchronization.
9. `CoverSeparator` scans the immutable base CSR view and owns reusable
   candidate/canonicalization buffers. Generated cuts are globally valid and
   never retain references to a node primal vector.
10. `CliqueTable` owns a compact adjacency structure, generation-stamped work
    queue, and literal marks derived from immutable global rows. Node
    propagation borrows reconstructed bounds and appends deductions to the same
    persistent domain-change arena.
11. Reduced-cost fixing owns no buffers: it borrows an optimal LP result and
    reconstructed domain, then stores every accepted tightening in the domain
    arena before the result view can be invalidated.
12. `ConflictPool` owns flat condition storage, SoA metadata, and one retained
    canonicalization buffer. It learns a path only when every condition has an
    exact integer negation; unsupported paths are rejected conservatively.
    Conflict identifiers never escape the pool, so compaction may rebuild all
    ranges and occurrence indices at a node boundary.

## Hot-path rules

- Default DFS must not allocate or materialize plugin views when selecting a
  node. Custom selectors may pay O(open nodes) for the restricted view.
- Branching appends to a parent-linked change arena instead of cloning complete
  bound arrays or ancestor paths. Child creation is amortized O(1).
- LP bound arrays use stable addresses throughout a solve so the LP session can
  distinguish value revisions from structural revisions.
- The root relaxation may use presolve. Later bound-only nodes solve in original
  coordinates and reuse the same basis/session through LP reoptimization.
- Propagation owns one solve-local workspace for row activities and round
  snapshots; processing another node performs no scratch-array allocation.
- Cover separation reuses solve-local candidate and cut buffers. The augmented
  CSC matrix is rebuilt only after a cut-pool revision, not for every LP solve.
- Clique extraction is bounded by `max_clique_row_entries` and
  `max_clique_edges`; node propagation is allocation-free apart from persistent
  deductions that descendants inherit.
- Reduced-cost fixing scans only ordinary binary/integer columns and allocates
  only when a useful persistent bound change is found; semi-variable
  disjunctions are deliberately excluded.
- Conflict propagation scans all active conflicts once, then revisits only
  ranges indexed by sparse node-path changes using two retained worklists.
  A shared-root revision change explicitly invalidates this shortcut and
  schedules every active conflict once.
  Learning reuses canonicalization capacity; entries, length, rounds, and age
  are all bounded by controls. Physical compaction runs only after aging and
  only when at least half of the conflict slots are inactive.
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
iterations, work, cold solves, reoptimizations, pseudocost observations, and
strong-branch probe/work counters, cut rounds/additions/expirations, and clique
edges/tightenings.
Comparisons require a pinned model corpus, identical compiler mode, warmups,
multiple samples, and objective/status parity. Recorded smoke measurements live in
`docs/mip-performance.md`.

## Deliberately incomplete HiGHS-aligned modules

`HeuristicContext.lp_basis` is an optional borrowed, read-only basis in original
model coordinates. `search.zig` supplies it only for a valid basis with no extra
cut rows. The diving plugin imports it into its own `LpSolveSession` for the first
fix, then reoptimizes that private session for subsequent bound changes. The
session is destroyed at the end of the invocation; neither the parent's basis
nor its solution buffers can be mutated by the heuristic. No solver ownership
or mutable engine reference crosses the plugin contract.

The current solver is branch-and-bound, not a complete HiGHS-class
branch-and-cut implementation. Add the following as separate modules rather
than expanding `search.zig`:

1. incremental/reversible propagation queues;
2. inference/cutoff scores and adaptive strong-branch budgets;
3. MIR/Gomory/clique separators, adaptive cut selection, and node-local cuts;
4. compressed threshold maps for very large root integer domains;
5. reason-frontier/recursive conflict minimization, non-binary implications, and
   two-watched conditions;
6. `symmetry.zig` and restart policy;
7. deterministic parallel node processing;
8. MIP event/callback bridging and structured timing statistics.

Each new algorithm needs a unit test, an end-to-end correctness corpus, and a
before/after benchmark or profiler trace before becoming a default component.
