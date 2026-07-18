//! Real Matrix Market validation and benchmark runner.
//!
//! Usage: matrix-dataset-runner DATASET_DIR REPORT.tsv

const std = @import("std");
const zhighs = @import("zhighs");

const Parsed = struct {
    matrix: zhighs.matrix.CscMatrix,
    parse_build_ns: u64,
};

const DatasetFile = struct { name: []const u8, size: u64 };

fn nowNs() i128 {
    var ts: std.posix.timespec = undefined;
    switch (std.posix.errno(std.posix.system.clock_gettime(.MONOTONIC, &ts))) {
        .SUCCESS => {},
        else => unreachable,
    }
    return @as(i128, ts.sec) * std.time.ns_per_s + ts.nsec;
}

inline fn clobber(pointer: anytype) void {
    asm volatile (""
        :
        : [pointer] "r" (pointer),
        : .{ .memory = true });
}

fn parseMatrixMarket(io: std.Io, allocator: std.mem.Allocator, path: []const u8) !Parsed {
    const started = nowNs();
    const file = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer file.close(io);
    var read_buffer: [64 * 1024]u8 = undefined;
    var reader = file.reader(io, &read_buffer);
    const content = try reader.interface.allocRemaining(allocator, .limited(2 * 1024 * 1024 * 1024));
    defer allocator.free(content);

    var lines = std.mem.tokenizeScalar(u8, content, '\n');
    const banner = std.mem.trim(u8, lines.next() orelse return error.InvalidMatrixMarket, " \t\r");
    var banner_parts = std.mem.tokenizeAny(u8, banner, " \t\r");
    if (!std.ascii.eqlIgnoreCase(banner_parts.next() orelse return error.InvalidMatrixMarket, "%%MatrixMarket")) return error.InvalidMatrixMarket;
    if (!std.ascii.eqlIgnoreCase(banner_parts.next() orelse return error.InvalidMatrixMarket, "matrix")) return error.InvalidMatrixMarket;
    if (!std.ascii.eqlIgnoreCase(banner_parts.next() orelse return error.InvalidMatrixMarket, "coordinate")) return error.UnsupportedMatrixMarket;
    const field = banner_parts.next() orelse return error.InvalidMatrixMarket;
    const symmetry = banner_parts.next() orelse return error.InvalidMatrixMarket;
    const is_pattern = std.ascii.eqlIgnoreCase(field, "pattern");
    if (!is_pattern and !std.ascii.eqlIgnoreCase(field, "real") and !std.ascii.eqlIgnoreCase(field, "integer")) return error.UnsupportedMatrixMarket;
    const is_symmetric = std.ascii.eqlIgnoreCase(symmetry, "symmetric") or std.ascii.eqlIgnoreCase(symmetry, "hermitian");
    if (!is_symmetric and !std.ascii.eqlIgnoreCase(symmetry, "general")) return error.UnsupportedMatrixMarket;

    var dimension_line: ?[]const u8 = null;
    while (lines.next()) |raw| {
        const line = std.mem.trim(u8, raw, " \t\r");
        if (line.len == 0 or line[0] == '%') continue;
        dimension_line = line;
        break;
    }
    var dims = std.mem.tokenizeAny(u8, dimension_line orelse return error.InvalidMatrixMarket, " \t\r");
    const rows = try std.fmt.parseUnsigned(usize, dims.next() orelse return error.InvalidMatrixMarket, 10);
    const cols = try std.fmt.parseUnsigned(usize, dims.next() orelse return error.InvalidMatrixMarket, 10);
    const stored_nnz = try std.fmt.parseUnsigned(usize, dims.next() orelse return error.InvalidMatrixMarket, 10);

    var builder = try zhighs.matrix.MatrixBuilder.init(rows, cols);
    defer builder.deinit(allocator);
    try builder.reserve(allocator, if (is_symmetric) stored_nnz * 2 else stored_nnz);
    var seen: usize = 0;
    while (lines.next()) |raw| {
        const line = std.mem.trim(u8, raw, " \t\r");
        if (line.len == 0 or line[0] == '%') continue;
        var parts = std.mem.tokenizeAny(u8, line, " \t\r");
        const row_one = try std.fmt.parseUnsigned(usize, parts.next() orelse return error.InvalidMatrixMarket, 10);
        const col_one = try std.fmt.parseUnsigned(usize, parts.next() orelse return error.InvalidMatrixMarket, 10);
        if (row_one == 0 or col_one == 0) return error.InvalidMatrixMarket;
        const value: f64 = if (is_pattern) 1.0 else try std.fmt.parseFloat(f64, parts.next() orelse return error.InvalidMatrixMarket);
        const row = try zhighs.RowId.fromUsize(row_one - 1);
        const col = try zhighs.ColId.fromUsize(col_one - 1);
        try builder.append(allocator, row, col, value);
        if (is_symmetric and row_one != col_one)
            try builder.append(allocator, try zhighs.RowId.fromUsize(col_one - 1), try zhighs.ColId.fromUsize(row_one - 1), value);
        seen += 1;
    }
    if (seen != stored_nnz) return error.InvalidMatrixMarket;
    const matrix = try builder.freeze(allocator, 0.0);
    return .{ .matrix = matrix, .parse_build_ns = @intCast(nowNs() - started) };
}

fn maxRelativeDifference(lhs: []const f64, rhs: []const f64) f64 {
    var result: f64 = 0.0;
    for (lhs, rhs) |a, b| result = @max(result, @abs(a - b) / @max(1.0, @max(@abs(a), @abs(b))));
    return result;
}

fn peakRssKb() usize {
    const value = std.posix.getrusage(std.posix.rusage.SELF).maxrss;
    return if (value > 0) @intCast(value) else 1;
}

fn currentRssKb(io: std.Io) !usize {
    const file = try std.Io.Dir.cwd().openFile(io, "/proc/self/status", .{});
    defer file.close(io);
    var read_buffer: [4096]u8 = undefined;
    var reader = file.reader(io, &read_buffer);
    var content: [4096]u8 = undefined;
    const length = try reader.interface.readSliceShort(&content);
    var lines = std.mem.tokenizeScalar(u8, content[0..length], '\n');
    while (lines.next()) |line| {
        if (!std.mem.startsWith(u8, line, "VmRSS:")) continue;
        var fields = std.mem.tokenizeAny(u8, line[6..], " \t");
        return try std.fmt.parseInt(usize, fields.next() orelse return error.InvalidProcStatus, 10);
    }
    return error.InvalidProcStatus;
}

fn writeLine(io: std.Io, file: std.Io.File, allocator: std.mem.Allocator, comptime format: []const u8, values: anytype) !void {
    const line = try std.fmt.allocPrint(allocator, format, values);
    defer allocator.free(line);
    try file.writeStreamingAll(io, line);
}

fn benchmarkOne(io: std.Io, allocator: std.mem.Allocator, dataset_dir: []const u8, name: []const u8, report: std.Io.File, details: std.Io.File) !void {
    const total_started = nowNs();
    const path = try std.fs.path.join(allocator, &.{ dataset_dir, name });
    defer allocator.free(path);
    var parsed = try parseMatrixMarket(io, allocator, path);
    defer parsed.matrix.deinit(allocator);
    try parsed.matrix.validate();

    const x = try allocator.alloc(f64, parsed.matrix.num_cols);
    defer allocator.free(x);
    const y_csc = try allocator.alloc(f64, parsed.matrix.num_rows);
    defer allocator.free(y_csc);
    const y_csr = try allocator.alloc(f64, parsed.matrix.num_rows);
    defer allocator.free(y_csr);
    for (x, 0..) |*value, index| value.* = 1.0 + @as(f64, @floatFromInt(index % 17)) * 0.03125;

    const repeats = std.math.clamp(100_000_000 / @max(parsed.matrix.nnz(), 1), 5, 50);
    parsed.matrix.multiplyAssumeValid(x, y_csc);
    var started = nowNs();
    for (0..repeats) |_| {
        parsed.matrix.multiplyAssumeValid(x, y_csc);
        clobber(y_csc.ptr);
    }
    const csc_spmv_ns: u64 = @intCast(@divTrunc(nowNs() - started, repeats));

    var csr_buffers = try zhighs.matrix.CsrBuffers.init(allocator, parsed.matrix.num_rows, parsed.matrix.nnz());
    defer csr_buffers.deinit(allocator);
    const transform_repeats: usize = 7;
    try zhighs.matrix.fillCsrFromCscAssumeValid(parsed.matrix, csr_buffers.row_starts, csr_buffers.col_indices, csr_buffers.values, csr_buffers.cursor);
    started = nowNs();
    for (0..transform_repeats) |_|
        try zhighs.matrix.fillCsrFromCscAssumeValid(parsed.matrix, csr_buffers.row_starts, csr_buffers.col_indices, csr_buffers.values, csr_buffers.cursor);
    const csr_build_ns: u64 = @intCast(@divTrunc(nowNs() - started, transform_repeats));
    const csr = zhighs.matrix.CsrView.initAssumeValid(parsed.matrix.num_rows, parsed.matrix.num_cols, csr_buffers.row_starts, csr_buffers.col_indices, csr_buffers.values);
    try csr.validate();
    csr.multiplyAssumeValid(x, y_csr);
    if (maxRelativeDifference(y_csc, y_csr) > 1e-12) return error.SemanticMismatch;
    started = nowNs();
    for (0..repeats) |_| {
        csr.multiplyAssumeValid(x, y_csr);
        clobber(y_csr.ptr);
    }
    const csr_spmv_ns: u64 = @intCast(@divTrunc(nowNs() - started, repeats));

    var transpose_buffers = try zhighs.matrix.TransposeBuffers.init(allocator, parsed.matrix.num_rows, parsed.matrix.nnz());
    defer transpose_buffers.deinit(allocator);
    try zhighs.matrix.transposeIntoAssumeValid(parsed.matrix, transpose_buffers.starts, transpose_buffers.rows, transpose_buffers.values, transpose_buffers.cursor);
    started = nowNs();
    for (0..transform_repeats) |_|
        try zhighs.matrix.transposeIntoAssumeValid(parsed.matrix, transpose_buffers.starts, transpose_buffers.rows, transpose_buffers.values, transpose_buffers.cursor);
    const transpose_ns: u64 = @intCast(@divTrunc(nowNs() - started, transform_repeats));
    const transposed = zhighs.matrix.CscMatrix.initBorrowedAssumeValid(parsed.matrix.num_cols, parsed.matrix.num_rows, transpose_buffers.starts, transpose_buffers.rows, transpose_buffers.values);
    try transposed.validate();

    const original_values = try allocator.dupe(f64, parsed.matrix.values);
    defer allocator.free(original_values);
    const row_factors = try allocator.alloc(f64, parsed.matrix.num_rows);
    defer allocator.free(row_factors);
    const col_factors = try allocator.alloc(f64, parsed.matrix.num_cols);
    defer allocator.free(col_factors);
    for (row_factors, 0..) |*factor, index| factor.* = if (index & 1 == 0) 0.5 else 2.0;
    for (col_factors, 0..) |*factor, index| factor.* = if (index & 1 == 0) 2.0 else 0.5;
    const scaling = zhighs.matrix.ScalingView{ .row = row_factors, .col = col_factors };
    try zhighs.matrix.applyScaling(&parsed.matrix, scaling);
    try zhighs.matrix.removeScaling(&parsed.matrix, scaling);
    if (!std.mem.eql(f64, original_values, parsed.matrix.values)) return error.ScalingRoundTripMismatch;

    const sample_cols = @min(parsed.matrix.num_cols, 4096);
    var sample = try zhighs.matrix.extractColumnRange(allocator, parsed.matrix, 0, sample_cols);
    defer sample.deinit(allocator);
    try sample.validate();
    const row_map = try allocator.alloc(zhighs.RowId, sample.num_rows);
    defer allocator.free(row_map);
    const col_map = try allocator.alloc(zhighs.ColId, sample.num_cols);
    defer allocator.free(col_map);
    for (row_map, 0..) |*id, index| id.* = zhighs.RowId.fromUsizeAssumeValid(sample.num_rows - 1 - index);
    for (col_map, 0..) |*id, index| id.* = zhighs.ColId.fromUsizeAssumeValid(sample.num_cols - 1 - index);
    var permuted = try zhighs.matrix.permuteAssumeValid(allocator, sample, row_map, col_map);
    defer permuted.deinit(allocator);
    try permuted.validate();
    var restored = try zhighs.matrix.permuteAssumeValid(allocator, permuted, row_map, col_map);
    defer restored.deinit(allocator);
    if (!zhighs.matrix.eql(sample, restored)) return error.PermutationRoundTripMismatch;

    const total_ns: u64 = @intCast(nowNs() - total_started);
    const rss = peakRssKb();
    try writeLine(io, report, allocator, "{s}\t{d}\t{d}\t{d}\t{d:.3}\t{d}\tPASS\n", .{
        name,                                    parsed.matrix.num_rows, parsed.matrix.num_cols, parsed.matrix.nnz(),
        @as(f64, @floatFromInt(total_ns)) / 1e6, rss,
    });
    try writeLine(io, details, allocator, "{s}\t{d}\t{d}\t{d}\t{d:.3}\t{d:.3}\t{d:.3}\t{d:.3}\t{d:.3}\t{d}\n", .{
        name,
        parsed.matrix.num_rows,
        parsed.matrix.num_cols,
        parsed.matrix.nnz(),
        @as(f64, @floatFromInt(parsed.parse_build_ns)) / 1e6,
        @as(f64, @floatFromInt(csc_spmv_ns)) / 1e6,
        @as(f64, @floatFromInt(csr_spmv_ns)) / 1e6,
        @as(f64, @floatFromInt(csr_build_ns)) / 1e6,
        @as(f64, @floatFromInt(transpose_ns)) / 1e6,
        rss,
    });
}

fn rssOne(io: std.Io, allocator: std.mem.Allocator, dataset_dir: []const u8, name: []const u8, report: std.Io.File) !void {
    const path = try std.fs.path.join(allocator, &.{ dataset_dir, name });
    defer allocator.free(path);
    var parsed = try parseMatrixMarket(io, allocator, path);
    defer parsed.matrix.deinit(allocator);
    try parsed.matrix.validate();
    const requested = if (parsed.matrix.storage) |storage| storage.len else parsed.matrix.col_starts.len * @sizeOf(usize) +
        parsed.matrix.row_indices.len * @sizeOf(zhighs.RowId) +
        parsed.matrix.values.len * @sizeOf(f64);
    try writeLine(io, report, allocator, "{s}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\n", .{
        name,      parsed.matrix.num_rows, parsed.matrix.num_cols, parsed.matrix.nnz(),
        requested, try currentRssKb(io),   peakRssKb(),
    });
}

pub fn main(init: std.process.Init) !void {
    const rss_only = init.environ_map.get("ZHIGHS_DATASET_RSS_ONLY") != null;
    const allocator = if (rss_only) std.heap.c_allocator else std.heap.smp_allocator;
    const io = init.io;
    var args = std.process.Args.iterate(init.minimal.args);
    _ = args.next();
    const dataset_dir = args.next() orelse return error.MissingDatasetDirectory;
    const report_path = args.next() orelse return error.MissingReportPath;
    if (args.next() != null) return error.TooManyArguments;

    const report = try std.Io.Dir.cwd().createFile(io, report_path, .{});
    defer report.close(io);
    if (rss_only)
        try report.writeStreamingAll(io, "dataset\trows\tcols\tnnz\trequested_bytes\tcurrent_rss_kb\tpeak_rss_kb\n")
    else
        try report.writeStreamingAll(io, "dataset\trows\tcols\tnnz\telapsed_ms\tpeak_rss_kb\tstatus\n");
    const details_path = try std.fmt.allocPrint(allocator, "{s}.details.tsv", .{report_path});
    defer allocator.free(details_path);
    const details = try std.Io.Dir.cwd().createFile(io, details_path, .{});
    defer details.close(io);
    if (!rss_only)
        try details.writeStreamingAll(io, "dataset\trows\tcols\tnnz\tparse_build_ms\tcsc_spmv_ms\tcsr_spmv_ms\tcsc_to_csr_ms\ttranspose_ms\tpeak_rss_kb\n");

    const filter = init.environ_map.get("ZHIGHS_DATASET_FILTER");

    var dir = try std.Io.Dir.cwd().openDir(io, dataset_dir, .{ .iterate = true });
    defer dir.close(io);
    var iterator = dir.iterate();
    var datasets: std.ArrayListUnmanaged(DatasetFile) = .empty;
    defer {
        for (datasets.items) |dataset| allocator.free(dataset.name);
        datasets.deinit(allocator);
    }
    while (try iterator.next(io)) |entry| {
        if (entry.kind != .file or !std.mem.endsWith(u8, entry.name, ".mtx")) continue;
        const stat = try dir.statFile(io, entry.name, .{});
        if (stat.size <= 1024 * 1024) continue;
        try datasets.append(allocator, .{ .name = try allocator.dupe(u8, entry.name), .size = stat.size });
    }
    std.mem.sort(DatasetFile, datasets.items, {}, struct {
        fn lessThan(_: void, lhs: DatasetFile, rhs: DatasetFile) bool {
            if (lhs.size != rhs.size) return lhs.size < rhs.size;
            return std.mem.lessThan(u8, lhs.name, rhs.name);
        }
    }.lessThan);
    if (datasets.items.len < 3) return error.InsufficientDatasets;
    for (datasets.items) |dataset| {
        if (filter) |selected| if (!std.mem.eql(u8, selected, dataset.name)) continue;
        std.debug.print("benchmarking {s}\n", .{dataset.name});
        if (rss_only)
            try rssOne(io, allocator, dataset_dir, dataset.name, report)
        else
            try benchmarkOne(io, allocator, dataset_dir, dataset.name, report, details);
    }
}
