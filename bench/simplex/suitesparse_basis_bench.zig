//! Derive a deterministic, certified-invertible basis from Matrix Market data.
//!
//! For the leading square submatrix, off-diagonal row i is scaled by
//! `0.5 / sum(abs(a_i*))` and the diagonal is replaced by one. Hence every
//! generated basis is strictly row diagonally dominant and nonsingular.

const std = @import("std");
const zhighs = @import("zhighs");

const Header = struct { rows: usize, columns: usize, pattern: bool, symmetric: bool };

fn header(lines: *std.mem.TokenIterator(u8, .scalar)) !Header {
    var parts = std.mem.tokenizeAny(u8, std.mem.trim(u8, lines.next() orelse return error.InvalidMatrixMarket, " \t\r"), " \t\r");
    if (!std.mem.eql(u8, parts.next() orelse return error.InvalidMatrixMarket, "%%MatrixMarket")) return error.InvalidMatrixMarket;
    _ = parts.next() orelse return error.InvalidMatrixMarket;
    if (!std.ascii.eqlIgnoreCase(parts.next() orelse return error.InvalidMatrixMarket, "coordinate")) return error.UnsupportedMatrixMarket;
    const field = parts.next() orelse return error.InvalidMatrixMarket;
    const symmetry = parts.next() orelse return error.InvalidMatrixMarket;
    var dimension_line: ?[]const u8 = null;
    while (lines.next()) |raw| {
        const line = std.mem.trim(u8, raw, " \t\r");
        if (line.len != 0 and line[0] != '%') { dimension_line = line; break; }
    }
    var dims = std.mem.tokenizeAny(u8, dimension_line orelse return error.InvalidMatrixMarket, " \t\r");
    return .{
        .rows = try std.fmt.parseUnsigned(usize, dims.next() orelse return error.InvalidMatrixMarket, 10),
        .columns = try std.fmt.parseUnsigned(usize, dims.next() orelse return error.InvalidMatrixMarket, 10),
        .pattern = std.ascii.eqlIgnoreCase(field, "pattern"),
        .symmetric = std.ascii.eqlIgnoreCase(symmetry, "symmetric") or std.ascii.eqlIgnoreCase(symmetry, "hermitian"),
    };
}

fn consume(content: []const u8, n: usize, row_sums: []f64, builder: ?*zhighs.matrix.MatrixBuilder, allocator: std.mem.Allocator) !void {
    var lines = std.mem.tokenizeScalar(u8, content, '\n');
    const info = try header(&lines);
    while (lines.next()) |raw| {
        const line = std.mem.trim(u8, raw, " \t\r");
        if (line.len == 0 or line[0] == '%') continue;
        var parts = std.mem.tokenizeAny(u8, line, " \t\r");
        const row = (try std.fmt.parseUnsigned(usize, parts.next() orelse return error.InvalidMatrixMarket, 10)) - 1;
        const column = (try std.fmt.parseUnsigned(usize, parts.next() orelse return error.InvalidMatrixMarket, 10)) - 1;
        const value: f64 = if (info.pattern) 1 else try std.fmt.parseFloat(f64, parts.next() orelse return error.InvalidMatrixMarket);
        if (row < n and column < n and row != column) {
            if (builder) |out| try out.append(allocator, try zhighs.RowId.fromUsize(row), try zhighs.ColId.fromUsize(column), value * 0.5 / row_sums[row]) else row_sums[row] += @abs(value);
        }
        if (info.symmetric and row != column and column < n and row < n) {
            if (builder) |out| try out.append(allocator, try zhighs.RowId.fromUsize(column), try zhighs.ColId.fromUsize(row), value * 0.5 / row_sums[column]) else row_sums[column] += @abs(value);
        }
    }
}

pub fn main(init: std.process.Init) !void {
    var args = std.process.Args.iterate(init.minimal.args);
    _ = args.next();
    const path = args.next() orelse return error.InvalidArguments;
    const limit = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 256;
    if (limit == 0 or args.next() != null) return error.InvalidArguments;
    const allocator = std.heap.c_allocator;
    const file = try std.Io.Dir.cwd().openFile(init.io, path, .{});
    defer file.close(init.io);
    var buffer: [64 * 1024]u8 = undefined;
    var reader = file.reader(init.io, &buffer);
    const content = try reader.interface.allocRemaining(allocator, .limited(2 * 1024 * 1024 * 1024));
    defer allocator.free(content);
    var lines = std.mem.tokenizeScalar(u8, content, '\n');
    const info = try header(&lines);
    const n = @min(limit, @min(info.rows, info.columns));
    const row_sums = try allocator.alloc(f64, n);
    defer allocator.free(row_sums);
    @memset(row_sums, 0);
    try consume(content, n, row_sums, null, allocator);
    // Empty source rows remain the identity row and need no scaling divisor.
    for (row_sums) |*sum| if (sum.* == 0) { sum.* = 1; };
    var builder = try zhighs.matrix.MatrixBuilder.init(n, n);
    defer builder.deinit(allocator);
    try consume(content, n, row_sums, &builder, allocator);
    for (0..n) |diagonal| try builder.append(allocator, try zhighs.RowId.fromUsize(diagonal), try zhighs.ColId.fromUsize(diagonal), 1);
    var matrix = try builder.freeze(allocator, 0);
    defer matrix.deinit(allocator);
    const basis = zhighs.matrix.SparseBasisView{
        .dimension = n,
        .starts = matrix.compact_col_starts orelse return error.DimensionTooLarge,
        .rows = matrix.row_indices,
        .values = matrix.values,
    };
    // Verify the construction certificate on the canonical, duplicate-merged
    // matrix rather than trusting the source triplet stream.
    var minimum_margin = std.math.inf(f64);
    for (0..n) |row| {
        var diagonal: f64 = 0;
        var off_diagonal: f64 = 0;
        for (0..n) |column| for (matrix.col_starts[column]..matrix.col_starts[column + 1]) |entry| if (matrix.row_indices[entry].toUsize() == row) {
            if (row == column) diagonal = @abs(matrix.values[entry]) else off_diagonal += @abs(matrix.values[entry]);
        };
        const margin = diagonal - off_diagonal;
        if (!(margin > 0)) return error.DiagonalDominanceCertificateFailed;
        minimum_margin = @min(minimum_margin, margin);
    }
    var lu = zhighs.matrix.SparseLU.init(allocator);
    defer lu.deinit();
    try lu.factorize(basis);
    const rhs = try allocator.alloc(f64, n);
    defer allocator.free(rhs);
    for (rhs, 0..) |*value, index| value.* = 1 + @as(f64, @floatFromInt(index % 11)) * 0.25;
    const original = try allocator.dupe(f64, rhs);
    defer allocator.free(original);
    try lu.solve(rhs);
    var residual: f64 = 0;
    for (0..n) |row| {
        var product: f64 = 0;
        for (0..n) |column| for (matrix.col_starts[column]..matrix.col_starts[column + 1]) |entry|
            if (matrix.row_indices[entry].toUsize() == row) { product += matrix.values[entry] * rhs[column]; };
        residual = @max(residual, @abs(product - original[row]));
    }
    std.debug.print("suitesparse-derived,{s},{d},{d},{d},{d},{e},{e},{d}\n", .{
        path, n, basis.nnz(), lu.factorNonzeros(), lu.inserted_fill, minimum_margin, residual, lu.requestedBytes(),
    });
}
