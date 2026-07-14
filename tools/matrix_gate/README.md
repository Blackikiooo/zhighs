# Matrix production acceptance gate

Daily development check:

```bash
tools/matrix_acceptance.sh quick
```

Production verdict:

```bash
HIGHS_SOURCE=/path/to/pinned/HiGHS \
MATRIX_DATASET_DIR=/path/to/pinned/matrix-market-corpus \
MATRIX_DATASET_RUNNER=/path/to/dataset-runner \
tools/matrix_acceptance.sh full
```

`full` is fail-closed. Missing HiGHS, commands, datasets, runner, malformed
reports, or a failed test all produce a non-zero exit status. `quick` only runs
the internal structural/OOM and Debug w32/w64 regression gates; its success is
never a production-ready verdict.

The dataset runner receives two arguments: dataset directory and output report.
It must validate matrix construction and canonical CSC, CSC/CSR/transpose
semantic equality, scaling round trips, slice/permutation transforms, runtime,
and peak RSS. It must write tab-separated columns:

```text
dataset  rows  cols  nnz  elapsed_ms  peak_rss_kb  status
```

At least three `PASS` datasets with 10,000 rows, 10,000 columns, 100,000 nnz,
positive runtime, and positive peak RSS are required. Corpus provenance and
checksums must be pinned alongside the external dataset directory.

The pinned local corpus used on 2026-07-14 is recorded in
`tools/matrix_gate/suitesparse-corpus.lock.tsv`. Build the in-tree runner with
`zig build build-matrix-dataset-runner -Doptimize=ReleaseFast -Dcpu=native` and
point `MATRIX_DATASET_RUNNER` at `zig-out/bin/matrix-dataset-runner`. The large
dataset gate verifies every matrix SHA-256 against this lock before executing;
`MATRIX_DATASET_LOCK` may select a different explicitly pinned corpus.
