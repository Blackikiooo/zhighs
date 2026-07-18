const std = @import("std");
const buildFoundation = @import("build/foundation.zig").buildFoundation;

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const foundation = buildFoundation(b, .{
        .target = target,
        .optimize = optimize,
    });

    const matrix_module = b.createModule(.{
        .root_source_file = b.path("src/matrix/root.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{.{ .name = "foundation", .module = foundation.module }},
    });
    const lp_module = b.createModule(.{
        .root_source_file = b.path("src/lp/root.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "foundation", .module = foundation.module },
            .{ .name = "matrix", .module = matrix_module },
        },
    });
    const solver_module = b.createModule(.{
        .root_source_file = b.path("src/solver/root.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "foundation", .module = foundation.module },
            .{ .name = "matrix", .module = matrix_module },
            .{ .name = "lp", .module = lp_module },
        },
    });
    const io_module = b.createModule(.{
        .root_source_file = b.path("src/io/root.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "foundation", .module = foundation.module },
            .{ .name = "matrix", .module = matrix_module },
        },
    });
    const model_module = b.createModule(.{
        .root_source_file = b.path("src/model/root.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "foundation", .module = foundation.module },
            .{ .name = "matrix", .module = matrix_module },
            .{ .name = "solver", .module = solver_module },
            .{ .name = "io", .module = io_module },
        },
    });

    const mod = b.addModule("zhighs", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });
    mod.addImport("foundation", foundation.module);
    mod.addImport("matrix", matrix_module);
    mod.addImport("model", model_module);
    mod.addImport("lp", lp_module);
    mod.addImport("solver", solver_module);
    mod.addImport("io", io_module);

    const exe = b.addExecutable(.{
        .name = "zhighs",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "zhighs", .module = mod },
            },
        }),
    });

    b.installArtifact(exe);
    const run_step = b.step("run", "Run the app");

    const run_cmd = b.addRunArtifact(exe);
    run_step.dependOn(&run_cmd.step);

    run_cmd.step.dependOn(b.getInstallStep());

    if (b.args) |args| {
        run_cmd.addArgs(args);
    }

    const mod_tests = b.addTest(.{
        .root_module = mod,
    });

    const run_mod_tests = b.addRunArtifact(mod_tests);

    const exe_tests = b.addTest(.{
        .root_module = exe.root_module,
    });

    // A run step that will run the second test executable.
    const run_exe_tests = b.addRunArtifact(exe_tests);

    // A top level step for running all tests. dependOn can be called multiple
    // times and since the two run steps do not depend on one another, this will
    // make the two of them run in parallel.
    const test_step = b.step("test", "Run tests");
    test_step.dependOn(foundation.test_step);
    test_step.dependOn(&run_mod_tests.step);
    test_step.dependOn(&run_exe_tests.step);

    const hcd_bench = b.addExecutable(.{
        .name = "hcd-bench",
        .root_module = b.createModule(.{
            .root_source_file = b.path("bench/hcd/hcd_bench.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "foundation", .module = foundation.module },
            },
        }),
    });
    const run_hcd_bench = b.addRunArtifact(hcd_bench);
    const hcd_bench_step = b.step("bench-hcd", "Run the HCD microbenchmark");
    hcd_bench_step.dependOn(&run_hcd_bench.step);

    const matrix_bench_module = b.createModule(.{
        .root_source_file = b.path("bench/matrix/root.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "foundation", .module = foundation.module },
            .{ .name = "matrix", .module = matrix_module },
        },
    });

    const matrix_bench = b.addExecutable(.{
        .name = "matrix-bench",
        .root_module = b.createModule(.{
            .root_source_file = b.path("bench/matrix/matrix_bench.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "zhighs", .module = matrix_bench_module },
            },
        }),
    });
    const run_matrix_bench = b.addRunArtifact(matrix_bench);
    const matrix_bench_step = b.step("bench-matrix", "Run sparse matrix microbenchmarks");
    matrix_bench_step.dependOn(&run_matrix_bench.step);

    const simplex_session_bench = b.addExecutable(.{
        .name = "simplex-session-bench",
        .root_module = b.createModule(.{
            .root_source_file = b.path("bench/simplex/session_bench.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "zhighs", .module = mod }},
        }),
    });
    const run_simplex_session_bench = b.addRunArtifact(simplex_session_bench);
    const simplex_session_bench_step = b.step("bench-simplex", "Run cold and warm simplex-session benchmarks");
    simplex_session_bench_step.dependOn(&run_simplex_session_bench.step);

    const install_matrix_bench = b.addInstallArtifact(matrix_bench, .{});
    const build_matrix_bench_step = b.step("build-bench-matrix", "Build the matrix benchmark without running it");
    build_matrix_bench_step.dependOn(&install_matrix_bench.step);

    const dense_lu_bench = b.addExecutable(.{
        .name = "dense-lu-bench",
        .root_module = b.createModule(.{
            .root_source_file = b.path("bench/matrix/dense_lu_bench.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "zhighs", .module = matrix_bench_module }},
        }),
    });
    const run_dense_lu_bench = b.addRunArtifact(dense_lu_bench);
    const dense_lu_bench_step = b.step("bench-dense-lu", "Run multi-dimension DenseLU FTRAN/BTRAN benchmarks");
    dense_lu_bench_step.dependOn(&run_dense_lu_bench.step);
    const install_dense_lu_bench = b.addInstallArtifact(dense_lu_bench, .{});
    const build_dense_lu_bench_step = b.step("build-bench-dense-lu", "Build the DenseLU benchmark for perf stat");
    build_dense_lu_bench_step.dependOn(&install_dense_lu_bench.step);

    const coefficient_edit_bench = b.addExecutable(.{
        .name = "coefficient-edit-bench",
        .root_module = b.createModule(.{
            .root_source_file = b.path("bench/matrix/coefficient_edit_bench.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "zhighs", .module = mod }},
        }),
    });
    const run_coefficient_edit_bench = b.addRunArtifact(coefficient_edit_bench);
    const coefficient_edit_bench_step = b.step("bench-coefficient-edits", "Run batched coefficient-edit throughput benchmarks");
    coefficient_edit_bench_step.dependOn(&run_coefficient_edit_bench.step);

    const io_parser_bench = b.addExecutable(.{
        .name = "io-parser-bench",
        .root_module = b.createModule(.{
            .root_source_file = b.path("bench/io/parser_bench.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "io", .module = io_module }},
        }),
    });
    const install_io_parser_bench = b.addInstallArtifact(io_parser_bench, .{});
    const build_io_bench_step = b.step("build-io-bench", "Build the model I/O parser benchmark");
    build_io_bench_step.dependOn(&install_io_parser_bench.step);
    const run_io_parser_bench = b.addRunArtifact(io_parser_bench);
    if (b.args) |args| run_io_parser_bench.addArgs(args);
    const io_bench_step = b.step("bench-io", "Benchmark model parsing and canonical construction");
    io_bench_step.dependOn(&run_io_parser_bench.step);

    const perf_profile = b.addExecutable(.{
        .name = "perf-profile",
        .root_module = b.createModule(.{
            .root_source_file = b.path("bench/matrix/perf_profile.zig"),
            .target = target,
            .optimize = optimize,
            // Enables an explicit c_allocator A/B against C++ malloc. The
            // profiling default remains smp_allocator for continuity.
            .link_libc = true,
            .imports = &.{
                .{ .name = "zhighs", .module = matrix_bench_module },
            },
        }),
    });
    const install_perf_profile = b.addInstallArtifact(perf_profile, .{});
    const build_perf_profile_step = b.step("build-perf-profile", "Build the perf profiling binary");
    build_perf_profile_step.dependOn(&install_perf_profile.step);

    const matrix_layout_audit = b.addExecutable(.{
        .name = "matrix-layout-audit",
        .root_module = b.createModule(.{
            .root_source_file = b.path("bench/matrix/layout_audit.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "zhighs", .module = matrix_bench_module }},
        }),
    });
    const run_matrix_layout_audit = b.addRunArtifact(matrix_layout_audit);
    const matrix_layout_audit_step = b.step("audit-matrix-layout", "Print matrix type sizes, alignments, and field offsets");
    matrix_layout_audit_step.dependOn(&run_matrix_layout_audit.step);

    const matrix_vector_experiment = b.addExecutable(.{
        .name = "matrix-vector-experiment",
        .root_module = b.createModule(.{
            .root_source_file = b.path("bench/matrix/vector_experiment.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });
    const run_matrix_vector_experiment = b.addRunArtifact(matrix_vector_experiment);
    const matrix_vector_experiment_step = b.step("experiment-matrix-vector", "Compare scalar and explicit @Vector matrix loops");
    matrix_vector_experiment_step.dependOn(&run_matrix_vector_experiment.step);

    const matrix_allocator_experiment = b.addExecutable(.{
        .name = "matrix-allocator-experiment",
        .root_module = b.createModule(.{
            .root_source_file = b.path("bench/matrix/allocator_experiment.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "zhighs", .module = matrix_bench_module }},
        }),
    });
    const run_matrix_allocator_experiment = b.addRunArtifact(matrix_allocator_experiment);
    const matrix_allocator_experiment_step = b.step("experiment-matrix-allocator", "Compare matrix scratch and session allocators");
    matrix_allocator_experiment_step.dependOn(&run_matrix_allocator_experiment.step);

    const matrix_dataset_runner = b.addExecutable(.{
        .name = "matrix-dataset-runner",
        .root_module = b.createModule(.{
            .root_source_file = b.path("bench/matrix/dataset_runner.zig"),
            .target = target,
            .optimize = optimize,
            .link_libc = true,
            .imports = &.{.{ .name = "zhighs", .module = matrix_bench_module }},
        }),
    });
    const install_matrix_dataset_runner = b.addInstallArtifact(matrix_dataset_runner, .{});
    const build_matrix_dataset_runner_step = b.step("build-matrix-dataset-runner", "Build the real Matrix Market validator and benchmark");
    build_matrix_dataset_runner_step.dependOn(&install_matrix_dataset_runner.step);

    // ── Matrix-only test (no model/API dependency) ──────────────
    // Used during matrix performance work so tests stay green while
    // higher layers (model, API) are being edited independently.

    const matrix_test_root = b.createModule(.{
        .root_source_file = b.path("test/matrix/root.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "matrix", .module = matrix_module },
        },
    });

    const matrix_tests = b.addTest(.{
        .root_module = matrix_test_root,
    });
    const run_matrix_tests = b.addRunArtifact(matrix_tests);
    const matrix_test_step = b.step("test-matrix", "Run matrix-only tests (no model/API)");
    matrix_test_step.dependOn(&run_matrix_tests.step);

    const matrix_acceptance_root = b.createModule(.{
        .root_source_file = b.path("test/matrix/acceptance.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "matrix", .module = matrix_module },
            .{ .name = "foundation", .module = foundation.module },
        },
    });
    const matrix_acceptance_tests = b.addTest(.{ .root_module = matrix_acceptance_root });
    const run_matrix_acceptance_tests = b.addRunArtifact(matrix_acceptance_tests);
    const matrix_acceptance_step = b.step("test-matrix-acceptance", "Run matrix structural property and failing-allocator gates");
    matrix_acceptance_step.dependOn(&run_matrix_acceptance_tests.step);

    // ── Model core test (no API/solver/presolve dependency) ─────
    // Tests the solver-internal model IR layer independently of the
    // Gurobi-style user Model, API bindings, presolve, and simplex.

    const model_test_root = b.createModule(.{
        .root_source_file = b.path("test/model/root.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "model", .module = model_module },
        },
    });

    const model_tests = b.addTest(.{
        .root_module = model_test_root,
    });
    const run_model_tests = b.addRunArtifact(model_tests);
    const model_test_step = b.step("test-model", "Run model core tests (no API/solver)");
    model_test_step.dependOn(&run_model_tests.step);

    const io_tests = b.addTest(.{ .root_module = io_module });
    const run_io_tests = b.addRunArtifact(io_tests);
    const io_test_step = b.step("test-io", "Run model file parser and writer tests");
    io_test_step.dependOn(&run_io_tests.step);
    test_step.dependOn(&run_io_tests.step);
}
