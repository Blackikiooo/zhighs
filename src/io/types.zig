//! Stable, model-independent types shared by every file format frontend.

const std = @import("std");
const matrix = @import("matrix");

/// Closed error vocabulary exposed by all readers and writers.
///
/// Frontends translate allocator, filesystem and grammar-specific failures
/// into this set so callers do not need to depend on a particular format.
pub const IoError = error{
    OutOfMemory,
    FileNotFound,
    PermissionDenied,
    ReadFailed,
    WriteFailed,
    FileTooLarge,
    ResourceLimitExceeded,
    Cancelled,
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
    InvalidTolerance,
    NonFiniteValue,
};

/// Logical syntax selected after stripping an optional compression suffix.
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

/// Outer byte-stream encoding. Parsing the enclosed model is a separate step.
pub const Compression = enum { none, gzip, bzip2, zip, seven_zip, xz };

/// Result of filename detection: model grammar plus outer compression layer.
pub const FileKind = struct {
    /// Grammar used to interpret the decompressed bytes.
    format: Format,
    /// Compression wrapper; `.none` means the file is directly parseable.
    compression: Compression = .none,
};

/// Objective direction encoded as a multiplier for normalization to minimization.
pub const ObjectiveSense = enum(i8) { minimize = 1, maximize = -1 };
/// Relational operator used by row-oriented text formats.
pub const RowSense = enum(u8) { less_equal = '<', equal = '=', greater_equal = '>' };
/// Integrality/domain classification attached to a model column.
pub const VariableType = enum(u8) {
    continuous = 'C',
    binary = 'B',
    integer = 'I',
    semi_continuous = 'S',
    semi_integer = 'N',
};

/// Source location and human-readable explanation for a parse failure.
///
/// Offsets are byte based; line and column are parser-maintained display
/// coordinates and remain zero when the frontend cannot determine them.
pub const Diagnostic = struct {
    /// Zero-based byte offset in the uncompressed source.
    byte_offset: usize = 0,
    /// Line number reported by the parser, or zero when unavailable.
    line: usize = 0,
    /// Byte column within `line`, or zero when unavailable.
    column: usize = 0,
    /// Borrowed diagnostic text; its owner must outlive the consumer.
    message: []const u8 = "",
};

/// Selects how an uncompressed model file is exposed to a format parser.
pub const InputMode = enum {
    /// Buffer small files and memory-map files at or above the configured
    /// threshold. This is the recommended production setting.
    automatic,
    /// Always copy the file through `std.Io` into allocator-owned memory.
    buffered,
    /// Request a read-only memory map. The `std.Io` backend may transparently
    /// fall back to aligned buffered storage when mapping is unavailable.
    memory_map,
};

/// Limits are part of the API so untrusted or unexpectedly large inputs cannot
/// consume unbounded resources.  A future streaming source can honor the same
/// contract without changing format parsers.
pub const ReadOptions = struct {
    /// Maximum accepted uncompressed file length.
    max_file_bytes: usize = 8 * 1024 * 1024 * 1024,
    /// Storage policy used to expose file bytes to the parser.
    input_mode: InputMode = .automatic,
    /// Automatic mode maps files at or above this size. Small files stay on
    /// the buffered path to avoid mmap setup and page-fault overhead.
    memory_map_threshold_bytes: usize = 256 * 1024,
    /// Ask the OS to prefault mapped pages. Disabled by default so mapping a
    /// very large model does not immediately force the complete file into RSS.
    memory_map_populate: bool = false,
    /// Coefficients whose magnitude is at most this value may be discarded.
    zero_tolerance: f64 = 0.0,
    /// Retain model object names in the published result. When false, parsers
    /// still use temporary name indexes for resolution, then release them
    /// before CSC finalization; final row/column name tables are empty.
    keep_names: bool = true,
    /// Hard semantic limits are checked before growing parser-owned storage.
    max_rows: usize = std.math.maxInt(usize),
    max_columns: usize = std.math.maxInt(usize),
    /// Maximum temporary matrix entries. LP duplicate coordinates merged at
    /// append time count once; general/MPS triplets count as encountered.
    max_matrix_terms: usize = std.math.maxInt(usize),
    /// Sum of unique row and column name bytes retained by semantic metadata.
    max_name_bytes: usize = std.math.maxInt(usize),
    /// Reject adversarial physical records before MPS tokenization and LP
    /// records when their terminating newline/EOF is reached.
    max_line_bytes: usize = 16 * 1024 * 1024,
    /// Reject a single LP token before numeric conversion or name interning.
    max_token_bytes: usize = 1024 * 1024,
    /// Optional caller-owned atomic flag. It must outlive `readFile`/`parse`.
    interrupt_flag: ?*const std.atomic.Value(bool) = null,
    /// Number of parser work units between atomic flag loads. A value of zero
    /// is treated as one. Lower values reduce cancellation latency.
    interrupt_check_interval: usize = 4096,

    /// Poll the caller's cancellation flag immediately.
    ///
    /// This method performs no allocation and leaves all parser state intact
    /// when it returns `error.Cancelled`.
    pub inline fn checkCancelled(self: ReadOptions) IoError!void {
        if (self.interrupt_flag) |flag| if (flag.load(.acquire)) return error.Cancelled;
    }
};

/// Low-overhead cooperative cancellation poller shared by text frontends.
pub const ParseControl = struct {
    /// Borrowed cancellation flag, or null when polling is disabled.
    flag: ?*const std.atomic.Value(bool),
    /// Number of work units between atomic loads; always at least one.
    interval: usize,
    /// Work units left before the next atomic load.
    remaining: usize,

    /// Construct a poller from public read options and normalize zero interval.
    pub fn init(options: ReadOptions) ParseControl {
        const interval = @max(options.interrupt_check_interval, 1);
        return .{ .flag = options.interrupt_flag, .interval = interval, .remaining = interval };
    }

    /// Load the flag now and reset the countdown.
    pub inline fn checkNow(self: *ParseControl) IoError!void {
        if (self.flag) |flag| if (flag.load(.acquire)) return error.Cancelled;
        self.remaining = self.interval;
    }

    /// Account for one parser work unit and poll when the countdown expires.
    pub inline fn tick(self: *ParseControl) IoError!void {
        if (self.flag == null) return;
        if (self.remaining > 1) {
            self.remaining -= 1;
            return;
        }
        try self.checkNow();
    }
};

/// Formatting controls shared by model writers.
pub const WriteOptions = struct {
    /// Emit retained/generated row and column names when the format supports them.
    emit_names: bool = true,
};

/// Borrowed, immutable model representation consumed by format writers.
/// It intentionally contains no dependency on the public `Model` type.
pub const ModelView = struct {
    /// Borrowed model name; an empty slice represents an unnamed model.
    name: []const u8,
    /// Direction of the objective before any solver normalization.
    objective_sense: ObjectiveSense,
    /// Constant term added to the linear objective.
    objective_offset: f64,
    /// Objective coefficient for each structural column.
    col_cost: []const f64,
    /// Inclusive lower bound for each structural column.
    col_lower: []const f64,
    /// Inclusive upper bound for each structural column.
    col_upper: []const f64,
    /// Domain classification for each structural column.
    col_type: []const VariableType,
    /// Empty means names were intentionally discarded. Otherwise the slice
    /// length must equal `numCols()`.
    col_names: []const ?[]const u8,
    /// Inclusive lower activity bound for each row.
    row_lower: []const f64,
    /// Inclusive upper activity bound for each row.
    row_upper: []const f64,
    /// Empty means names were intentionally discarded. Otherwise the slice
    /// length must equal `numRows()`.
    row_names: []const ?[]const u8,
    /// Borrowed constraint matrix in canonical compressed-column form.
    matrix: matrix.CscView,

    /// Return the structural column count derived from the cost array.
    pub inline fn numCols(self: ModelView) usize {
        return self.col_cost.len;
    }

    /// Return the row count derived from the row-bound array.
    pub inline fn numRows(self: ModelView) usize {
        return self.row_lower.len;
    }
};

test {
    std.testing.refAllDecls(@This());
}

test "parse control observes atomic cancellation at configured interval" {
    var interrupted = std.atomic.Value(bool).init(false);
    var control = ParseControl.init(.{ .interrupt_flag = &interrupted, .interrupt_check_interval = 2 });
    try control.tick();
    interrupted.store(true, .release);
    try std.testing.expectError(error.Cancelled, control.tick());
}
