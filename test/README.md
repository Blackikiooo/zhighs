# Test suites

Source-local `test` blocks cover small unit tests. As implementation grows,
this directory will contain:

- `differential/`: comparisons against the pinned local HiGHS C API;
- `fuzz/`: model, matrix, parser, and presolve fuzz targets;
- `regression/`: minimal instances for previously fixed defects;
- `instances/`: small redistributable LP/MIP/QP fixtures and provenance.

Performance measurements belong in `bench/`, not here.
