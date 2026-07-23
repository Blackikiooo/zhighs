//! Filename suffix recognition.  Dispatch is centralized so every frontend
//! agrees on aliases and compressed suffix handling.

const std = @import("std");
const types = @import("types.zig");

/// Detect the model grammar and optional compression wrapper from `path`.
///
/// Matching is ASCII case-insensitive. The compression suffix is removed
/// before inspecting the model extension. No file is opened and returned
/// slices never borrow from the path.
pub fn detect(path: []const u8) types.IoError!types.FileKind {
    if (path.len == 0) return error.UnsupportedFormat;
    var base = path;
    var compression: types.Compression = .none;
    const compressed = [_]struct { suffix: []const u8, kind: types.Compression }{
        .{ .suffix = ".gz", .kind = .gzip },
        .{ .suffix = ".bz2", .kind = .bzip2 },
        .{ .suffix = ".zip", .kind = .zip },
        .{ .suffix = ".7z", .kind = .seven_zip },
        .{ .suffix = ".xz", .kind = .xz },
    };
    for (compressed) |entry| {
        if (std.ascii.endsWithIgnoreCase(base, entry.suffix)) {
            compression = entry.kind;
            base = base[0 .. base.len - entry.suffix.len];
            break;
        }
    }
    const ext = std.fs.path.extension(base);
    const format: types.Format = if (std.ascii.eqlIgnoreCase(ext, ".lp"))
        .lp
    else if (std.ascii.eqlIgnoreCase(ext, ".rlp"))
        .rlp
    else if (std.ascii.eqlIgnoreCase(ext, ".mps"))
        .mps
    else if (std.ascii.eqlIgnoreCase(ext, ".rew"))
        .rew
    else if (std.ascii.eqlIgnoreCase(ext, ".dua"))
        .dua
    else if (std.ascii.eqlIgnoreCase(ext, ".dlp"))
        .dlp
    else if (std.ascii.eqlIgnoreCase(ext, ".ilp"))
        .ilp
    else if (std.ascii.eqlIgnoreCase(ext, ".opb"))
        .opb
    else
        return error.UnsupportedFormat;
    return .{ .format = format, .compression = compression };
}

test "detects model and compression suffix independently" {
    try std.testing.expectEqual(types.Format.mps, (try detect("large.MPS.gz")).format);
    try std.testing.expectEqual(types.Compression.gzip, (try detect("large.MPS.gz")).compression);
    try std.testing.expectEqual(types.Format.lp, (try detect("model.lp")).format);
}
