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
}
