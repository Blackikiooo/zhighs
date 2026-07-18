//! Deterministic sparse-LU acceptance corpus.
//!
//! Fixtures are built before timing. Each nonsingular family is strictly
//! row-diagonally dominant, so generation itself supplies an invertibility
//! certificate without relying on a dense reference factorization.

const std = @import("std");
const zhighs = @import("zhighs");

const Family = enum { high_fill, block_arrow, network, ill_conditioned, rank_deficient };

fn peakRssKb() usize {
    return @intCast(std.posix.getrusage(std.posix.rusage.SELF).maxrss);
}

fn nowNs() i128 {
    var ts: std.posix.timespec = undefined;
    switch (std.posix.errno(std.posix.system.clock_gettime(.MONOTONIC, &ts))) {
        .SUCCESS => {},
        else => unreachable,
    }
    return @as(i128, ts.sec) * std.time.ns_per_s + ts.nsec;
}

fn present(family: Family, row: usize, column: usize, n: usize) bool {
    const effective_column = if (family == .rank_deficient and column + 1 == n) n - 2 else column;
    if (row == effective_column) return true;
    return switch (family) {
        // A deterministic unsymmetric pattern whose fronts become much denser
        // than the input under elimination.
        .high_fill => (row > column and (row + 3 * column) % 7 == 0) or
            (row < column and (5 * row + column) % 11 == 0),
        .block_arrow => blk: {
            const block = @max(@as(usize, 4), n / 8);
            break :blk row / block == column / block or row + 1 == n or column + 1 == n;
        },
        // Grounded, diagonally shifted network adjacency with a few long arcs.
        .network => row + 1 == column or column + 1 == row or
            (row + 13) % n == column or (column + 13) % n == row,
        .ill_conditioned, .rank_deficient => row + 1 == effective_column or effective_column + 1 == row,
    };
}

fn coefficient(family: Family, row: usize, column: usize, n: usize) f64 {
    const effective_column = if (family == .rank_deficient and column + 1 == n) n - 2 else column;
    if (row != effective_column) return if ((row * 17 + effective_column * 29) & 1 == 0) 0.125 else -0.125;
    if (family == .ill_conditioned) {
        const exponent: i32 = @intCast((12 * column) / @max(n - 1, 1));
        return std.math.pow(f64, 10.0, -@as(f64, @floatFromInt(exponent)));
    }
    var off_diagonal: usize = 0;
    for (0..n) |other| if (other != row and present(family, row, other, n)) { off_diagonal += 1; };
    return 1.0 + @as(f64, @floatFromInt(off_diagonal)) * 0.125;
}

fn residual(basis: zhighs.matrix.SparseBasisView, x: []const f64, rhs: []const f64) f64 {
    var maximum: f64 = 0;
    for (0..basis.dimension) |row| {
        var product: f64 = 0;
        for (0..basis.dimension) |column|
            for (@as(usize, basis.starts[column])..@as(usize, basis.starts[column + 1])) |entry|
                if (basis.rows[entry].toUsize() == row) { product += basis.values[entry] * x[column]; };
        maximum = @max(maximum, @abs(product - rhs[row]));
    }
    return maximum;
}

pub fn main(init: std.process.Init) !void {
    var args = std.process.Args.iterate(init.minimal.args);
    _ = args.next();
    const n = if (args.next()) |text| try std.fmt.parseUnsigned(usize, text, 10) else 128;
    if (n < 16 or n > 2048 or args.next() != null) return error.InvalidArguments;
    const allocator = std.heap.c_allocator;
    const starts = try allocator.alloc(zhighs.foundation.HUInt, n + 1);
    defer allocator.free(starts);
    const rows = try allocator.alloc(zhighs.foundation.RowId, n * n);
    defer allocator.free(rows);
    const values = try allocator.alloc(f64, n * n);
    defer allocator.free(values);
    const rhs = try allocator.alloc(f64, n);
    defer allocator.free(rhs);
    const original = try allocator.alloc(f64, n);
    defer allocator.free(original);

    std.debug.print("family,dimension,basis_nnz,status,peeled,kernel_dim,kernel_nnz,max_row,max_col,ordering,factor_nnz,fill,invert_ns,max_residual,requested_bytes,peak_rss_kb\n", .{});
    inline for (std.meta.tags(Family)) |family| {
        var nnz: usize = 0;
        for (0..n) |column| {
            starts[column] = @intCast(nnz);
            for (0..n) |row| if (present(family, row, column, n)) {
                const value = coefficient(family, row, column, n);
                if (value != 0) {
                    rows[nnz] = zhighs.foundation.RowId.fromUsizeAssumeValid(row);
                    values[nnz] = value;
                    nnz += 1;
                }
            };
        }
        starts[n] = @intCast(nnz);
        const basis = zhighs.matrix.SparseBasisView{ .dimension = n, .starts = starts, .rows = rows[0..nnz], .values = values[0..nnz] };
        var lu = zhighs.matrix.SparseLU.init(allocator);
        defer lu.deinit();
        const began = nowNs();
        const result = lu.factorizeAssumeValid(basis);
        const elapsed: u64 = @intCast(nowNs() - began);
        if (result) |_| {
            if (family == .rank_deficient) return error.ExpectedSingular;
            for (rhs, original, 0..) |*value, *saved, index| {
                value.* = 1.0 + @as(f64, @floatFromInt(index % 17)) * 0.0625;
                saved.* = value.*;
            }
            try lu.solve(rhs);
            std.debug.print("{s},{d},{d},ok,{d},{d},{d},{d},{d},{s},{d},{d},{d},{e},{d},{d}\n", .{
                @tagName(family), n, nnz, lu.peeled_pivots, lu.kernel_dimension, lu.kernel_nonzeros,
                lu.kernel_maximum_row_count, lu.kernel_maximum_column_count, @tagName(lu.selected_ordering),
                lu.factorNonzeros(), lu.inserted_fill, elapsed,
                residual(basis, rhs, original), lu.requestedBytes(), peakRssKb(),
            });
        } else |err| {
            if (family != .rank_deficient or err != error.Singular) return err;
            std.debug.print("{s},{d},{d},expected-singular,{d},{d},{d},{d},{d},{s},0,0,{d},0,0,{d}\n", .{
                @tagName(family), n, nnz, lu.peeled_pivots, lu.kernel_dimension, lu.kernel_nonzeros,
                lu.kernel_maximum_row_count, lu.kernel_maximum_column_count, @tagName(lu.selected_ordering), elapsed, peakRssKb(),
            });
        }
    }
}
