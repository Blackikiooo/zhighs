//! Matrix production-gate tests: structural properties and exhaustive OOM.

const std = @import("std");
const matrix = @import("matrix");
const foundation = @import("foundation");

fn allocationFailureScenario(allocator: std.mem.Allocator) !void {
    var builder = try matrix.MatrixBuilder.init(8, 8);
    defer builder.deinit(allocator);
    for (0..8) |col| {
        const col_id = try foundation.ColId.fromUsize(col);
        try builder.append(allocator, try foundation.RowId.fromUsize(col), col_id, 2.0);
        if (col + 1 < 8)
            try builder.append(allocator, try foundation.RowId.fromUsize(col + 1), col_id, -1.0);
    }
    var source = try builder.freeze(allocator, 0.0);
    defer source.deinit(allocator);

    const selected = [_]foundation.RowId{
        try foundation.RowId.init(0),
        try foundation.RowId.init(2),
        try foundation.RowId.init(4),
        try foundation.RowId.init(6),
    };
    var sliced = try matrix.extractRows(allocator, source, &selected);
    defer sliced.deinit(allocator);
    const row_map = [_]foundation.RowId{
        try foundation.RowId.init(3),
        try foundation.RowId.init(1),
        try foundation.RowId.init(0),
        try foundation.RowId.init(2),
    };
    const col_map = [_]foundation.ColId{
        try foundation.ColId.init(7),
        try foundation.ColId.init(6),
        try foundation.ColId.init(5),
        try foundation.ColId.init(4),
        try foundation.ColId.init(3),
        try foundation.ColId.init(2),
        try foundation.ColId.init(1),
        try foundation.ColId.init(0),
    };
    var permuted = try matrix.permute(allocator, sliced, &row_map, &col_map);
    defer permuted.deinit(allocator);
    try permuted.validate();
}

test "all structural owning allocations are leak-free under injected OOM" {
    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        allocationFailureScenario,
        .{},
    );
}

test "random structural transform chain remains canonical" {
    var prng = std.Random.DefaultPrng.init(0x6d61_7472_6978_6761);
    const random = prng.random();
    for (0..100) |_| {
        const rows_count = random.intRangeAtMost(usize, 1, 24);
        const cols_count = random.intRangeAtMost(usize, 1, 20);
        var builder = try matrix.MatrixBuilder.init(rows_count, cols_count);
        defer builder.deinit(std.testing.allocator);
        for (0..200) |_| {
            const row = random.intRangeLessThan(usize, 0, rows_count);
            const col = random.intRangeLessThan(usize, 0, cols_count);
            const value: f64 = @floatFromInt(random.intRangeAtMost(i8, -4, 4));
            try builder.append(
                std.testing.allocator,
                try foundation.RowId.fromUsize(row),
                try foundation.ColId.fromUsize(col),
                value,
            );
        }
        var source = try builder.freeze(std.testing.allocator, 0.0);
        defer source.deinit(std.testing.allocator);
        try source.validate();

        var transposed = try matrix.transpose(std.testing.allocator, source);
        defer transposed.deinit(std.testing.allocator);
        var restored = try matrix.transpose(std.testing.allocator, transposed);
        defer restored.deinit(std.testing.allocator);
        try restored.validate();
        try std.testing.expect(matrix.eql(source, restored));

        const appended_row_starts = [_]usize{ 0, 1, 1 };
        const appended_cols = [_]foundation.ColId{try foundation.ColId.fromUsize(random.intRangeLessThan(usize, 0, cols_count))};
        const appended_values = [_]f64{3.0};
        var appended = try matrix.appendRowsFromCsr(
            std.testing.allocator,
            source,
            &appended_row_starts,
            &appended_cols,
            &appended_values,
        );
        defer appended.deinit(std.testing.allocator);
        try appended.validate();
    }
}
