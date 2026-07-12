const std = @import("std");
const buildFoundation = @import("build/foundation.zig").buildFoundation;

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const foundation = buildFoundation(b, .{
        .target = target,
        .optimize = optimize,
    });

    const mod = b.addModule("zhighs", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });
    mod.addImport("foundation", foundation.module);

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

    const matrix_module = b.createModule(.{
        .root_source_file = b.path("src/matrix/root.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "foundation", .module = foundation.module },
        },
    });
    // Make matrix available to model files (used via @import("matrix")).
    mod.addImport("matrix", matrix_module);

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

    const install_matrix_bench = b.addInstallArtifact(matrix_bench, .{});
    const build_matrix_bench_step = b.step("build-bench-matrix", "Build the matrix benchmark without running it");
    build_matrix_bench_step.dependOn(&install_matrix_bench.step);

    const perf_profile = b.addExecutable(.{
        .name = "perf-profile",
        .root_module = b.createModule(.{
            .root_source_file = b.path("bench/matrix/perf_profile.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "zhighs", .module = matrix_bench_module },
            },
        }),
    });
    const install_perf_profile = b.addInstallArtifact(perf_profile, .{});
    const build_perf_profile_step = b.step("build-perf-profile", "Build the perf profiling binary");
    build_perf_profile_step.dependOn(&install_perf_profile.step);

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

    // ── Model core test (no API/solver/presolve dependency) ─────
    // Tests the solver-internal model IR layer independently of the
    // Gurobi-style user Model, API bindings, presolve, and simplex.

    const model_module = b.createModule(.{
        .root_source_file = b.path("src/model/root.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "foundation", .module = foundation.module },
            .{ .name = "matrix", .module = matrix_module },
        },
    });

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
}
