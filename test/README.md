# Test suites

Source-local `test` blocks cover small unit tests. This directory contains:

- `matrix/`: matrix-only tests and acceptance gates (`zig build test-matrix`, `test-matrix-acceptance`);
- `model/`: model core tests (`zig build test-model`);
- `zig build test-end-to-end`: model lifecycle regressions plus CLI LP/MPS
  parsing, LP/MILP optimization, presolve/postsolve, objective and status
  reporting, exported primal values, and stopping before the first incumbent.
  Included in `zig build test`; also supports `-Dhighs-int-width=w64`.
- `mip/`: MIP known-optimum regression suite (`zig build test-mip-regression`) —
  solves 8 small/medium instances under 4 solver profiles and asserts
  status=optimal with objective matching the reference value (1e-6).

Planned additions:

- `differential/`: comparisons against the pinned local HiGHS C API;
- `fuzz/`: model, matrix, parser, and presolve fuzz targets;
- `instances/`: small redistributable LP/MIP/QP fixtures and provenance.

Performance measurements belong in `bench/`, not here.

The diving tests cover both cold initialization and a borrowed parent basis,
including preservation of the parent primal vector. HiGHS performance and
objective comparisons are a separate, explicitly invoked workflow documented
in `bench/mip/README.md`; they require the reference package and optional public
instances and are not dependencies of the hermetic Zig tests.
