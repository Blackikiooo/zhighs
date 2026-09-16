//! Command-line entry point: read a model file (MPS/LP) and solve it.
//!
//! Usage:
//!   zhighs [options] <model.mps|model.lp>
//!
//! Options:
//!   --param Key=Value   Set a solve parameter (repeatable), e.g.
//!                       --param TimeLimit=10 --param NodeLimit=100
//!   --solution <file>   Write the primal solution to <file>, one
//!                       "variable value" pair per line.
//!   --help              Show this usage text.
//!
//! Prints the solve status, objective value, MIP bound/gap/node count when
//! the model contains integer variables, wall-clock solve time, and the
//! non-zero entries of the primal solution.

const std = @import("std");
const zhighs = @import("zhighs");

const api = zhighs.api;

const usage_text =
    \\usage: zhighs [options] <model.mps|model.lp>
    \\
    \\options:
    \\  --param Key=Value   set a solve parameter (repeatable), e.g.
    \\                      --param TimeLimit=10 --param NodeLimit=100 --param MIPGap=0.01
    \\                      --param Presolve=0
    \\                      --param Method=2 (continuous LP interior-point)
    \\  --solution <file>   write the primal solution to <file>
    \\                      (one "variable value" pair per line)
    \\  --help, -h          show this usage text
    \\
;

fn nowNs() i128 {
    const timestamp = std.Io.Clock.now(.awake, std.Io.Threaded.global_single_threaded.io());
    return timestamp.nanoseconds;
}

fn fail(io_context: std.Io, comptime fmt: []const u8, args: anytype) noreturn {
    const stderr = std.Io.File.stderr();
    var buffer: [1024]u8 = undefined;
    const message = std.fmt.bufPrint(&buffer, fmt, args) catch "error";
    stderr.writeStreamingAll(io_context, "zhighs: ") catch {};
    stderr.writeStreamingAll(io_context, message) catch {};
    stderr.writeStreamingAll(io_context, "\n") catch {};
    std.process.exit(1);
}

/// Known parameter names, mirrored from `api.params` so the CLI rejects
/// typos instead of silently storing a no-op entry.
const known_param_names = blk: {
    const decls = @typeInfo(api.params).@"struct".decls;
    var names: [decls.len][]const u8 = undefined;
    for (decls, 0..) |decl, i| names[i] = @field(api.params, decl.name);
    break :blk names;
};

fn isKnownParam(name: []const u8) bool {
    for (known_param_names) |known| if (std.mem.eql(u8, name, known)) return true;
    return false;
}

/// Apply one `Key=Value` pair to the environment, dispatching on the value's
/// syntax: integer, floating-point, or plain string.
fn applyParam(env: *api.Env, pair: []const u8) !void {
    const eq = std.mem.indexOfScalar(u8, pair, '=') orelse return error.InvalidParameterSyntax;
    const name = pair[0..eq];
    const value = pair[eq + 1 ..];
    if (name.len == 0 or value.len == 0) return error.InvalidParameterSyntax;
    if (!isKnownParam(name)) return error.InvalidParameter;
    if (std.fmt.parseInt(i64, value, 10)) |int_value| {
        try env.setIntParam(name, int_value);
    } else |_| if (std.fmt.parseFloat(f64, value)) |dbl_value| {
        try env.setDblParam(name, dbl_value);
    } else |_| {
        try env.setStrParam(name, value);
    }
}

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;
    const io_context = init.io;

    var model_path: ?[]const u8 = null;
    var solution_path: ?[]const u8 = null;
    var param_pairs = std.ArrayList([]const u8).empty;
    defer param_pairs.deinit(allocator);

    const stdout = std.Io.File.stdout();

    var args = try init.minimal.args.iterateAllocator(allocator);
    defer args.deinit();
    _ = args.next();
    while (args.next()) |arg| {
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            try stdout.writeStreamingAll(io_context, usage_text);
            return;
        } else if (std.mem.eql(u8, arg, "--param")) {
            const pair = args.next() orelse
                fail(io_context, "--param requires a Key=Value argument", .{});
            try param_pairs.append(allocator, pair);
        } else if (std.mem.startsWith(u8, arg, "--param=")) {
            try param_pairs.append(allocator, arg["--param=".len..]);
        } else if (std.mem.eql(u8, arg, "--solution")) {
            if (solution_path != null) fail(io_context, "--solution given more than once", .{});
            solution_path = args.next() orelse
                fail(io_context, "--solution requires a file path", .{});
        } else if (std.mem.startsWith(u8, arg, "--solution=")) {
            if (solution_path != null) fail(io_context, "--solution given more than once", .{});
            solution_path = arg["--solution=".len..];
        } else if (std.mem.startsWith(u8, arg, "--")) {
            fail(io_context, "unknown option '{s}' (try --help)", .{arg});
        } else {
            if (model_path != null) fail(io_context, "multiple model files given; expected exactly one", .{});
            model_path = arg;
        }
    }
    const path = model_path orelse
        fail(io_context, "missing model file (try --help)", .{});

    var env = try api.Env.initSimple(allocator);
    defer env.deinit();
    for (param_pairs.items) |pair| {
        applyParam(&env, pair) catch
            fail(io_context, "invalid parameter '{s}'", .{pair});
    }

    var model = try api.Model.init(allocator, &env, path);
    defer model.deinit();

    model.read(path) catch
        fail(io_context, "cannot read model file '{s}'", .{path});
    const num_vars = try model.getIntAttr(.num_vars);
    const num_constrs = try model.getIntAttr(.num_constrs);
    const is_mip = (try model.getIntAttr(.is_mip)) != 0;

    var line = try std.fmt.allocPrint(allocator, "model {s}: {d} variables, {d} constraints, class={s}\n", .{
        path, num_vars, num_constrs, if (is_mip) "MILP" else "LP",
    });
    try stdout.writeStreamingAll(io_context, line);
    allocator.free(line);

    const started = nowNs();
    try model.optimize();
    const solve_ns: u64 = @intCast(nowNs() - started);

    const status = model.status;
    const obj_val = try model.getDblAttr(.obj_val);
    const sol_count = try model.getIntAttr(.sol_count);
    const have_solution = sol_count != 0;
    line = try std.fmt.allocPrint(allocator, "status: {s}\ntime: {d:.3} s\n", .{
        status.label(), @as(f64, @floatFromInt(solve_ns)) / std.time.ns_per_s,
    });
    try stdout.writeStreamingAll(io_context, line);
    allocator.free(line);
    if (have_solution) {
        line = try std.fmt.allocPrint(allocator, "objective: {d:.6}\n", .{obj_val});
        defer allocator.free(line);
        try stdout.writeStreamingAll(io_context, line);
    } else {
        try stdout.writeStreamingAll(io_context, "objective: unavailable\n");
    }

    if (is_mip) {
        const obj_bound = try model.getDblAttr(.obj_bound);
        const node_count = try model.getIntAttr(.node_count);
        const scale = @max(@abs(obj_val), 1e-10);
        const gap = if (have_solution and std.math.isFinite(obj_val) and std.math.isFinite(obj_bound))
            @abs(obj_val - obj_bound) / scale
        else
            std.math.inf(f64);
        line = try std.fmt.allocPrint(allocator, "bound: {d:.6}\ngap: {d:.4}\nnodes: {d}\nsolutions: {d}\n", .{
            obj_bound, gap, node_count, sol_count,
        });
        try stdout.writeStreamingAll(io_context, line);
        allocator.free(line);
    }

    if (!have_solution) return;

    try stdout.writeStreamingAll(io_context, "solution (non-zero):\n");
    var printed: usize = 0;
    for (0..@intCast(num_vars)) |i| {
        const value = model.getDblAttrElement(.x, i) catch continue;
        if (@abs(value) <= 1e-9) continue;
        const name = model.getStrAttrElement(.var_name, i) catch null;
        line = if (name) |n|
            try std.fmt.allocPrint(allocator, "  {s} = {d:.6}\n", .{ n, value })
        else
            try std.fmt.allocPrint(allocator, "  x[{d}] = {d:.6}\n", .{ i, value });
        try stdout.writeStreamingAll(io_context, line);
        allocator.free(line);
        printed += 1;
    }
    if (printed == 0) try stdout.writeStreamingAll(io_context, "  (all variables at zero)\n");

    if (solution_path) |solution_file| {
        const file = std.Io.Dir.cwd().createFile(io_context, solution_file, .{}) catch
            fail(io_context, "cannot write solution file '{s}'", .{solution_file});
        defer file.close(io_context);
        var write_buffer: [4096]u8 = undefined;
        var file_writer = file.writer(io_context, &write_buffer);
        for (0..@intCast(num_vars)) |i| {
            const raw_value = model.getDblAttrElement(.x, i) catch continue;
            const value = if (raw_value == 0.0) @as(f64, 0.0) else raw_value;
            const name = model.getStrAttrElement(.var_name, i) catch null;
            if (name) |n| {
                file_writer.interface.print("{s} {d:.10}\n", .{ n, value }) catch
                    fail(io_context, "cannot write solution file '{s}'", .{solution_file});
            } else {
                file_writer.interface.print("x[{d}] {d:.10}\n", .{ i, value }) catch
                    fail(io_context, "cannot write solution file '{s}'", .{solution_file});
            }
        }
        file_writer.interface.flush() catch
            fail(io_context, "cannot write solution file '{s}'", .{solution_file});
        line = try std.fmt.allocPrint(allocator, "solution written to {s}\n", .{solution_file});
        try stdout.writeStreamingAll(io_context, line);
        allocator.free(line);
    }
}
