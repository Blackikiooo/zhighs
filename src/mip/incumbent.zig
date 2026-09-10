//! Validation and ownership transfer for primal incumbent candidates.

const std = @import("std");
const ProblemView = @import("types.zig").ProblemView;

pub const Candidate = struct { primal: []f64, objective: f64 };

pub fn copyIfFeasible(allocator: std.mem.Allocator, problem: ProblemView, values: []const f64, tolerance: f64) !?Candidate {
    if (values.len != problem.relaxation.num_cols) return null;
    const primal = try allocator.dupe(f64, values);
    errdefer allocator.free(primal);
    for (problem.integrality, primal, problem.relaxation.col_lower, problem.relaxation.col_upper) |kind, *value, lower, upper| {
        if (!std.math.isFinite(value.*)) {
            allocator.free(primal);
            return null;
        }
        switch (kind) {
            .continuous => {},
            .integer => {
                if (@abs(value.* - @round(value.*)) > tolerance) {
                    allocator.free(primal);
                    return null;
                }
                value.* = @round(value.*);
            },
            .binary => {
                if (@abs(value.* - @round(value.*)) > tolerance or value.* < -tolerance or value.* > 1 + tolerance) {
                    allocator.free(primal);
                    return null;
                }
                value.* = @round(value.*);
            },
            .semi_continuous, .semi_integer => {
                if (@abs(value.*) <= tolerance) {
                    value.* = 0;
                    continue;
                }
                if (kind == .semi_integer) {
                    if (@abs(value.* - @round(value.*)) > tolerance) {
                        allocator.free(primal);
                        return null;
                    }
                    value.* = @round(value.*);
                }
            },
        }
        if (value.* < lower - tolerance or value.* > upper + tolerance) {
            allocator.free(primal);
            return null;
        }
    }

    const activity = try allocator.alloc(f64, problem.relaxation.num_rows);
    defer allocator.free(activity);
    @memset(activity, 0);
    for (0..problem.relaxation.num_cols) |column| {
        const begin = problem.relaxation.matrix.col_starts[column];
        const end = problem.relaxation.matrix.col_starts[column + 1];
        for (problem.relaxation.matrix.row_indices[begin..end], problem.relaxation.matrix.values[begin..end]) |row, coefficient|
            activity[row.toUsize()] += coefficient * primal[column];
    }
    for (activity, problem.relaxation.row_lower, problem.relaxation.row_upper) |value, lower, upper| {
        if (value < lower - tolerance or value > upper + tolerance) {
            allocator.free(primal);
            return null;
        }
    }
    var objective = problem.relaxation.objective_offset;
    for (problem.relaxation.col_cost, primal) |cost, value| objective += cost * value;
    return .{ .primal = primal, .objective = objective };
}
