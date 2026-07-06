const std = @import("std");

const HighsIntWidth = enum {
    w32,
    w64,
};

pub const TypesBuild = struct {
    module: *std.Build.Module,
    test_step: *std.Build.Step,
};

pub fn buildTypes(
    b: *std.Build,
    target: std.Build.ResolvedTarget,
) TypesBuild {
    const mod = b.createModule(.{
        .root_source_file = b.path("src/types/root.zig"),
        .target = target,
    });

    const highs_int_width = b.option(HighsIntWidth, "highs-int-width", "Integer width for HighsInt: i32 or i64") orelse .w32;
    const config = b.addOptions();
    config.addOption(HighsIntWidth, "highs_int_width", highs_int_width);

    mod.addOptions("config", config);

    const tests = b.addTest(.{
        .root_module = mod,
    });
    const run_tests = b.addRunArtifact(tests);

    return .{
        .module = mod,
        .test_step = &run_tests.step,
    };
}
