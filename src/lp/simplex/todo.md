# Simplex roadmap

## Implemented baseline

- General column bounds: lower, upper, boxed, free, and fixed variables.
- General row bounds: less-than, greater-than, equality, ranged, and free rows.
- Primal revised simplex with bound-aware entering directions and bound flips.
- Phase I artificial objective, infeasibility detection, and stable artificial-basic cleanup.
- Dense LU base factorization with allocation-free Eta FTRAN/BTRAN updates and periodic reinversion.
- Allocation-free final primal/dual feasibility validation.
- Borrowed `ProblemView` input and borrowed `SolutionView` output.
- Monotonic wall-clock limits and caller-owned atomic interruption.
- `Model.optimize` dispatch with stable public solution/status publication.
- Validated zero-copy `BasisView` and owning `BasisSnapshot` import/export.
- Warm basis refactorization with primal/dual feasibility classification.
- Dual revised simplex pricing, ratio test, and boxed-variable BFRT.
- Automatic RHS/bounds dual reoptimization and objective primal reoptimization.
- Fine-grained structure/matrix/bounds/objective revisions and a persistent
  owning solve session that reuses its basis, LU factors, Eta chain, and SoA
  workspaces across compatible solves.
- Allocation-free basis residual evaluation and iterative refinement directly
  from CSC plus basis membership.
- Contiguous Devex reference-weight storage and weighted primal/dual scans.
- Pivot-spread stability monitoring and factorization FTRAN/BTRAN/update stats.
- ReleaseFast cold/warm, FTRAN, BTRAN, and BFRT microbenchmarks.

Artificial basics are retained only for rank-redundant rows that have no stable
non-artificial replacement. Presolve row removal will eliminate those rows.

## Completed solve controls

- Deterministic work limits shared across Phase I, primal, and dual iterations.
- Borrowed structured `ProgressEventView` reporting with independently
  controlled callback and logging intervals.
- Iteration callbacks with no allocation or synchronization and only a nullable
  branch in disabled-callback hot paths.
- Model-level `WorkLimit`, `OutputFlag`, and `SimplexLogInterval` integration,
  including adapters for existing simplex and message callbacks.

## Priority 2: advanced dual and MIP reoptimization

- Piecewise-linear dual Phase I maximizing the dual-infeasibility sum. Imported
  bases that are neither primal nor dual feasible already retain their basis
  and factorization through a zero-auxiliary-objective dual feasibility repair,
  followed by primal or dual Phase II.
- Incremental Forrest--Goldfarb dual steepest-edge updates are implemented with one
  allocation-free FTRAN-DSE per pivot, exact BTRAN initialization after
  reinversion, selected-row accuracy correction, and full reset after severe
  weight underestimation. `SimplexPricing=2` enables the strategy.
- Hyper-sparse dual leaving-row pricing uses a fixed-capacity top-attractiveness
  candidate list, cutoff-based stale-list rebuild, and automatic activation
  only after tableau density falls below 10%. Candidate rows and scores are
  engine-owned SoA arrays with no iteration-time allocation.
- MIP-node warm starts.

## Priority 3: sparse factorization

- [x] DOD compact sparse-basis assembly with retaining SoA buffers,
  compile-time target prefetch policy, w32/w64 tests, and a standalone
  perf/disassembly benchmark. Design and fair HiGHS comparison rules are in
  `docs/sparse-basis-design.md`.
- [x] Reusable row-entry symbolic workspace, count buckets, fill-free singleton
  elimination, and deterministic threshold Markowitz selection for the first
  kernel pivot. Later pivots correctly wait for numerical fill updates.
- [x] Mutable SoA kernel matrix, recycled fill slots, row/column count-bucket
  updates, repeated numerical threshold Markowitz, packed L/U, and
  allocation-free FTRAN/BTRAN MVP.
- [x] Non-allocator INVERT tuning: dead-pivot bucket suppression, direct CSC
  pool load, redundant entry-state removal, theoretical Markowitz lower-bound
  stop, explicit trusted zero-copy entry, and warm no-allocation test.
- [x] Dual sparse-ordering backends with reduced-kernel dispatch, HiGHS-style
  bounded row/column Markowitz search, dynamic pre-kernel singleton peeling,
  and fixed pivot-trace separation of ordering from numerical update cost.
- [x] Numerically validated pivot-trace prefix replay with automatic suffix
  repair, cost-gated adaptive bounded-Markowitz search windows, and a
  shape-gated single-candidate frontier for compact high-fill peeled kernels.
- [x] Generation-marked dense lookup for sparse Schur row accumulation, with
  O(n) clearing required only on u32 generation wraparound.
- Integrate sparse LU behind the factorization backend with dense fallback,
  rank-deficiency repair, iterative refinement, and production statistics.
- [x] SparseLU is integrated as the large-basis factorization backend with a
  small dense fallback. Mutable intrusive U/UR, captured partial FTRAN `aq`,
  partial BTRAN `ep`, pivotal row/column deletion, spike insertion, FT row
  corrections, growth-based reinversion, and repeated-replacement dense-oracle
  differential tests are implemented. Updated hyper-sparse solves currently
  fall back to the sequential FT kernel until update-graph reachability lands.
- [x] Reachability-based hyper-sparse FTRAN/BTRAN with dense-output adaptive
  dispatch and an explicit sparse-output API.
- Retain dense LU only as a small-basis fallback and correctness oracle.

## Priority 4: pricing and numerical robustness

- Exact Devex recurrence and primal/dual steepest-edge weights. The current
  reference-weight implementation provides the owning layout and weighted
  scan, but intentionally does not claim the exact recurrence.
- Partial and hyper-sparse pricing.
- Row/column scaling. Residual recomputation and iterative refinement are done.
- Formal condition estimation, perturbation, anti-cycling, and stronger
  small-pivot rejection. A pivot-spread warning signal is already integrated.

## Validation and performance gates

- Degenerate, cycling, rank-deficient, and ill-conditioned unit tests.
- Netlib and Mittelmann result comparison against HiGHS/CLP.
- Pricing-only benchmark and larger representative model corpus. Cold/warm,
  FTRAN, BTRAN, and BFRT microbenchmarks are available via `bench-simplex`.
