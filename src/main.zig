//! Command-line entry point: read a model file (MPS/LP) and solve it.
//!
//! Usage:
//!   zhighs <model.mps|model.lp>
//!
//! Prints the solve status, objective value, MIP bound/gap/node count when
//! the model contains integer variables, wall-clock solve time, and the
//! non-zero entries of the primal solution.

const std = @import("std");
const zhighs = @import("zhighs");

const api = zhighs.api;

fn nowNs() i128 {
    const timestamp = std.Io.Clock.now(.awake, std.Io.Threaded.global_single_threaded.io());
    return timestamp.nanoseconds;
}

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;
    const io_context = init.io;

    var args = try init.minimal.args.iterateAllocator(allocator);
    defer args.deinit();
    _ = args.next();
    const path = args.next() orelse {
        try std.Io.File.stderr().writeStreamingAll(io_context, "usage: zhighs <model.mps|model.lp>\n");
        return error.MissingModelPath;
    };
    if (args.next() != null) return error.InvalidArguments;

    var env = try api.Env.initSimple(allocator);
    defer env.deinit();
    var model = try api.Model.init(allocator, &env, path);
    defer model.deinit();

    try model.read(path);
    const num_vars = try model.getIntAttr(.num_vars);
    const num_constrs = try model.getIntAttr(.num_constrs);
    const is_mip = (try model.getIntAttr(.is_mip)) != 0;

    const stdout = std.Io.File.stdout();
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
    line = try std.fmt.allocPrint(allocator, "status: {s}\nobjective: {d:.6}\ntime: {d:.3} s\n", .{
        status.label(), obj_val, @as(f64, @floatFromInt(solve_ns)) / std.time.ns_per_s,
    });
    try stdout.writeStreamingAll(io_context, line);
    allocator.free(line);

    if (is_mip) {
        const obj_bound = try model.getDblAttr(.obj_bound);
        const node_count = try model.getIntAttr(.node_count);
        const sol_count = try model.getIntAttr(.sol_count);
        const scale = @max(@abs(obj_val), 1e-10);
        const gap = if (std.math.isFinite(obj_val) and std.math.isFinite(obj_bound))
            @abs(obj_val - obj_bound) / scale
        else
            std.math.inf(f64);
        line = try std.fmt.allocPrint(allocator, "bound: {d:.6}\ngap: {d:.4}\nnodes: {d}\nsolutions: {d}\n", .{
            obj_bound, gap, node_count, sol_count,
        });
        try stdout.writeStreamingAll(io_context, line);
        allocator.free(line);
    }

    if (status != .optimal and sol_count_unavailable: {
        break :sol_count_unavailable (model.getIntAttr(.sol_count) catch 0) == 0;
    }) return;

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
}
