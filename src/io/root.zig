//! High-throughput, model-independent optimization file I/O.
//!
//! Format frontends share ownership, diagnostics, limits, and canonical CSC
//! output, but retain independent grammars. Public model methods are adapters;
//! parsers in this module never mutate solver state.

const std = @import("std");

pub const types = @import("types.zig");
pub const format = @import("format.zig");
pub const ModelData = @import("model_data.zig").ModelData;
pub const lp = @import("lp/root.zig");
pub const mps = @import("mps/root.zig");
const output = @import("output.zig");

pub const IoError = types.IoError;
pub const Format = types.Format;
pub const Compression = types.Compression;
pub const FileKind = types.FileKind;
pub const ObjectiveSense = types.ObjectiveSense;
pub const RowSense = types.RowSense;
pub const VariableType = types.VariableType;
pub const Diagnostic = types.Diagnostic;
pub const ReadOptions = types.ReadOptions;
pub const WriteOptions = types.WriteOptions;
pub const ModelView = types.ModelView;

/// Read and parse a model file. The returned model owns canonical numeric data
/// and is independent of the input buffer.
pub fn readFile(io_context: std.Io, allocator: std.mem.Allocator, path: []const u8, options: ReadOptions) IoError!ModelData {
    const kind = try format.detect(path);
    if (kind.compression != .none) return error.UnsupportedCompression;
    const file = std.Io.Dir.cwd().openFile(io_context, path, .{}) catch |err| return mapOpenError(err);
    defer file.close(io_context);
    var buffer: [128 * 1024]u8 = undefined;
    var file_reader = file.reader(io_context, &buffer);
    const input = file_reader.interface.allocRemaining(allocator, .limited(options.max_file_bytes)) catch |err| return switch (err) {
        error.StreamTooLong => error.FileTooLarge,
        error.OutOfMemory => error.OutOfMemory,
        else => error.ReadFailed,
    };
    defer allocator.free(input);
    const name = modelName(path, kind);
    return switch (kind.format) {
        .lp, .rlp => lp.parse(allocator, input, name, options),
        .mps, .rew => mps.parse(allocator, input, name, options),
        else => error.UnsupportedFormat,
    };
}

/// Write a borrowed model using the representation selected by `path`.
pub fn writeFile(io_context: std.Io, allocator: std.mem.Allocator, path: []const u8, model: ModelView, options: WriteOptions) IoError!void {
    const kind = try format.detect(path);
    if (kind.compression != .none) return error.UnsupportedCompression;
    const file = std.Io.Dir.cwd().createFile(io_context, path, .{}) catch return error.WriteFailed;
    defer file.close(io_context);
    var buffer: [128 * 1024]u8 = undefined;
    var file_writer = file.writer(io_context, &buffer);
    var names = try output.Names.init(allocator, model);
    defer names.deinit();
    switch (kind.format) {
        .lp => try lp.write(&file_writer.interface, allocator, model, names, options),
        .mps => try mps.write(&file_writer.interface, allocator, model, names, options),
        else => return error.UnsupportedFormat,
    }
    file_writer.interface.flush() catch return error.WriteFailed;
}

fn mapOpenError(err: anyerror) IoError {
    return switch (err) {
        error.FileNotFound => error.FileNotFound,
        error.AccessDenied => error.PermissionDenied,
        else => error.ReadFailed,
    };
}

fn modelName(path: []const u8, kind: FileKind) []const u8 {
    var base = std.fs.path.basename(path);
    if (kind.compression != .none) base = base[0 .. base.len - std.fs.path.extension(base).len];
    const extension = std.fs.path.extension(base);
    return base[0 .. base.len - extension.len];
}

test {
    std.testing.refAllDecls(@This());
}

test "LP and MPS writers round trip the same canonical linear model" {
    const source =
        \\NAME ROUNDTRIP
        \\OBJSENSE
        \\ MAX
        \\ROWS
        \\ N OBJ
        \\ L CAP
        \\ E BAL
        \\COLUMNS
        \\ X OBJ -2 CAP 3
        \\ X BAL 1
        \\ Y OBJ 1 CAP -1
        \\ Y BAL 2
        \\RHS
        \\ RHS1 CAP 7 BAL 4
        \\ RHS1 OBJ -3
        \\RANGES
        \\ RNG1 CAP 2
        \\BOUNDS
        \\ UP BND X 4
        \\ BV BND Y
        \\ENDATA
    ;
    var original = try mps.parse(std.testing.allocator, source, "fallback", .{});
    defer original.deinit();
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const lp_path = try std.fmt.allocPrint(std.testing.allocator, ".zig-cache/tmp/{s}/roundtrip.lp", .{tmp.sub_path});
    defer std.testing.allocator.free(lp_path);
    const mps_path = try std.fmt.allocPrint(std.testing.allocator, ".zig-cache/tmp/{s}/roundtrip.mps", .{tmp.sub_path});
    defer std.testing.allocator.free(mps_path);
    try writeFile(std.testing.io, std.testing.allocator, lp_path, original.view(), .{});
    try writeFile(std.testing.io, std.testing.allocator, mps_path, original.view(), .{});
    var from_lp = try readFile(std.testing.io, std.testing.allocator, lp_path, .{});
    defer from_lp.deinit();
    var from_mps = try readFile(std.testing.io, std.testing.allocator, mps_path, .{});
    defer from_mps.deinit();
    try std.testing.expectEqual(original.col_cost.len, from_lp.col_cost.len);
    try std.testing.expectEqual(original.row_lower.len, from_lp.row_lower.len);
    try std.testing.expectEqual(original.matrix.nnz(), from_lp.matrix.nnz());
    try std.testing.expectEqual(original.col_cost.len, from_mps.col_cost.len);
    try std.testing.expectEqual(original.row_lower.len, from_mps.row_lower.len);
    try std.testing.expectEqual(original.matrix.nnz(), from_mps.matrix.nnz());
    try std.testing.expectEqualSlices(f64, original.col_cost, from_lp.col_cost);
    try std.testing.expectEqualSlices(f64, original.row_lower, from_mps.row_lower);
    try std.testing.expectEqualSlices(f64, original.row_upper, from_lp.row_upper);
    try std.testing.expectEqual(@as(f64, 5.0), from_lp.row_lower[0]);
    try std.testing.expectEqual(@as(f64, 7.0), from_mps.row_upper[0]);
    try std.testing.expectEqual(@as(f64, 3.0), from_lp.objective_offset);
    try std.testing.expectEqual(@as(f64, 3.0), from_mps.objective_offset);
}
