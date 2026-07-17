//! Owned lifetime wrapper for model-file input bytes.
//!
//! Format frontends consume only a borrowed `[]const u8`; this module decides
//! whether those bytes come from an allocator-owned buffer or a read-only file
//! mapping. Keeping that policy here lets LP and MPS share the same mmap path
//! without coupling their grammars to files, operating systems, or `std.Io`.

const builtin = @import("builtin");
const std = @import("std");
const types = @import("types.zig");

pub const FileInput = union(enum) {
    empty,
    buffered: []u8,
    memory_map: std.Io.File.MemoryMap,

    /// Load an already-open regular file according to `ReadOptions`.
    /// The file must remain open until `deinit` destroys a memory map.
    pub fn load(io: std.Io, allocator: std.mem.Allocator, file: std.Io.File, options: types.ReadOptions) types.IoError!FileInput {
        try options.checkCancelled();
        const stat = file.stat(io) catch |err| return mapStatError(err);
        if (stat.kind != .file) return error.ReadFailed;
        if (stat.size > options.max_file_bytes) return error.FileTooLarge;
        const size = std.math.cast(usize, stat.size) orelse return error.FileTooLarge;
        if (size == 0) return .empty;

        const use_memory_map = switch (options.input_mode) {
            .automatic => size >= options.memory_map_threshold_bytes,
            .buffered => false,
            .memory_map => true,
        };
        if (!use_memory_map) return loadBuffered(io, allocator, file, size, options);

        const mapping = file.createMemoryMap(io, .{
            .len = size,
            .protection = .{ .read = true },
            .populate = options.memory_map_populate,
        }) catch |err| return mapMemoryMapError(err);
        return .{ .memory_map = mapping };
    }

    /// Borrow the complete source. The slice becomes invalid after `deinit`.
    pub fn bytes(self: *const FileInput) []const u8 {
        return switch (self.*) {
            .empty => "",
            .buffered => |memory| memory,
            .memory_map => |mapping| mapping.memory,
        };
    }

    pub fn deinit(self: *FileInput, io: std.Io, allocator: std.mem.Allocator) void {
        switch (self.*) {
            .empty => {},
            .buffered => |memory| allocator.free(memory),
            .memory_map => |*mapping| mapping.destroy(io),
        }
        self.* = undefined;
    }

    /// True only when the backend obtained a native mapping rather than its
    /// aligned-read fallback. Intended for tests and input-path telemetry.
    pub fn isNativeMemoryMap(self: *const FileInput) bool {
        return switch (self.*) {
            .memory_map => |mapping| mapping.section != null,
            else => false,
        };
    }
};

fn loadBuffered(io: std.Io, allocator: std.mem.Allocator, file: std.Io.File, size: usize, options: types.ReadOptions) types.IoError!FileInput {
    var buffer: [128 * 1024]u8 = undefined;
    var file_reader = file.reader(io, &buffer);
    const memory = allocator.alloc(u8, size) catch return error.OutOfMemory;
    errdefer allocator.free(memory);
    var offset: usize = 0;
    while (offset < memory.len) {
        try options.checkCancelled();
        const end = @min(offset +| buffer.len, memory.len);
        const read = file_reader.interface.readSliceShort(memory[offset..end]) catch return error.ReadFailed;
        if (read == 0) return error.ReadFailed;
        offset += read;
    }
    try options.checkCancelled();
    return .{ .buffered = memory };
}

fn mapStatError(err: anyerror) types.IoError {
    return switch (err) {
        error.AccessDenied, error.PermissionDenied => error.PermissionDenied,
        else => error.ReadFailed,
    };
}

fn mapMemoryMapError(err: anyerror) types.IoError {
    return switch (err) {
        error.OutOfMemory => error.OutOfMemory,
        error.AccessDenied, error.PermissionDenied => error.PermissionDenied,
        else => error.ReadFailed,
    };
}

test "forced memory map borrows regular file bytes" {
    const source = "Minimize\n obj: x\nSubject To\n c: x <= 1\nEnd\n";
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const path = try std.fmt.allocPrint(std.testing.allocator, ".zig-cache/tmp/{s}/mapped.lp", .{tmp.sub_path});
    defer std.testing.allocator.free(path);

    const created = try std.Io.Dir.cwd().createFile(std.testing.io, path, .{});
    try created.writeStreamingAll(std.testing.io, source);
    created.close(std.testing.io);

    const file = try std.Io.Dir.cwd().openFile(std.testing.io, path, .{});
    defer file.close(std.testing.io);
    var input = try FileInput.load(std.testing.io, std.testing.allocator, file, .{ .input_mode = .memory_map });
    defer input.deinit(std.testing.io, std.testing.allocator);
    try std.testing.expectEqualStrings(source, input.bytes());
    if (builtin.os.tag == .linux) try std.testing.expect(input.isNativeMemoryMap());
}

test "buffered input is exact-size and honors cancellation" {
    const source = "Minimize\n obj: x\nEnd\n";
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const path = try std.fmt.allocPrint(std.testing.allocator, ".zig-cache/tmp/{s}/buffered.lp", .{tmp.sub_path});
    defer std.testing.allocator.free(path);
    const created = try std.Io.Dir.cwd().createFile(std.testing.io, path, .{});
    try created.writeStreamingAll(std.testing.io, source);
    created.close(std.testing.io);

    const file = try std.Io.Dir.cwd().openFile(std.testing.io, path, .{});
    defer file.close(std.testing.io);
    var input = try FileInput.load(std.testing.io, std.testing.allocator, file, .{ .input_mode = .buffered });
    defer input.deinit(std.testing.io, std.testing.allocator);
    try std.testing.expectEqualStrings(source, input.bytes());

    var interrupted = std.atomic.Value(bool).init(true);
    try std.testing.expectError(error.Cancelled, FileInput.load(std.testing.io, std.testing.allocator, file, .{
        .input_mode = .buffered,
        .interrupt_flag = &interrupted,
    }));
}

test "input checks size limit before mapping" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const path = try std.fmt.allocPrint(std.testing.allocator, ".zig-cache/tmp/{s}/limited.lp", .{tmp.sub_path});
    defer std.testing.allocator.free(path);

    const created = try std.Io.Dir.cwd().createFile(std.testing.io, path, .{});
    try created.writeStreamingAll(std.testing.io, "12345");
    created.close(std.testing.io);

    const file = try std.Io.Dir.cwd().openFile(std.testing.io, path, .{});
    defer file.close(std.testing.io);
    try std.testing.expectError(error.FileTooLarge, FileInput.load(std.testing.io, std.testing.allocator, file, .{
        .input_mode = .memory_map,
        .max_file_bytes = 4,
    }));
}

test "empty files do not request a zero-length mapping" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const path = try std.fmt.allocPrint(std.testing.allocator, ".zig-cache/tmp/{s}/empty.lp", .{tmp.sub_path});
    defer std.testing.allocator.free(path);

    const created = try std.Io.Dir.cwd().createFile(std.testing.io, path, .{});
    created.close(std.testing.io);
    const file = try std.Io.Dir.cwd().openFile(std.testing.io, path, .{});
    defer file.close(std.testing.io);
    var input = try FileInput.load(std.testing.io, std.testing.allocator, file, .{ .input_mode = .memory_map });
    defer input.deinit(std.testing.io, std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 0), input.bytes().len);
    try std.testing.expect(!input.isNativeMemoryMap());
}
