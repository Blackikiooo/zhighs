//! Stable, model-independent types shared by every file format frontend.

const std = @import("std");
const matrix = @import("matrix");

pub const IoError = error{
    OutOfMemory,
    FileNotFound,
    PermissionDenied,
    ReadFailed,
    WriteFailed,
    FileTooLarge,
    UnsupportedFormat,
    UnsupportedCompression,
    UnsupportedFeature,
    InvalidSyntax,
    InvalidNumber,
    InvalidName,
    DuplicateName,
    UnknownName,
    InvalidBounds,
    InvalidDimensions,
    NonFiniteValue,
};

pub const Format = enum {
    lp,
    rlp,
    mps,
    rew,
    dua,
    dlp,
    ilp,
    opb,
};

pub const Compression = enum { none, gzip, bzip2, zip, seven_zip, xz };

pub const FileKind = struct {
    format: Format,
    compression: Compression = .none,
};

pub const ObjectiveSense = enum(i8) { minimize = 1, maximize = -1 };
pub const RowSense = enum(u8) { less_equal = '<', equal = '=', greater_equal = '>' };
pub const VariableType = enum(u8) {
    continuous = 'C',
    binary = 'B',
    integer = 'I',
    semi_continuous = 'S',
    semi_integer = 'N',
};

pub const Diagnostic = struct {
    byte_offset: usize = 0,
    line: usize = 0,
    column: usize = 0,
    message: []const u8 = "",
};

/// Limits are part of the API so untrusted or unexpectedly large inputs cannot
/// consume unbounded resources.  A future streaming source can honor the same
/// contract without changing format parsers.
pub const ReadOptions = struct {
    max_file_bytes: usize = 8 * 1024 * 1024 * 1024,
    zero_tolerance: f64 = 0.0,
    keep_names: bool = true,
};

pub const WriteOptions = struct {
    emit_names: bool = true,
};

/// Borrowed, immutable model representation consumed by format writers.
/// It intentionally contains no dependency on the public `Model` type.
pub const ModelView = struct {
    name: []const u8,
    objective_sense: ObjectiveSense,
    objective_offset: f64,
    col_cost: []const f64,
    col_lower: []const f64,
    col_upper: []const f64,
    col_type: []const VariableType,
    col_names: []const ?[]const u8,
    row_lower: []const f64,
    row_upper: []const f64,
    row_names: []const ?[]const u8,
    matrix: matrix.CscView,

    pub inline fn numCols(self: ModelView) usize {
        return self.col_cost.len;
    }

    pub inline fn numRows(self: ModelView) usize {
        return self.row_lower.len;
    }
};

test {
    std.testing.refAllDecls(@This());
}
