//! Free and traditional fixed-field MPS frontend.
//!
//! Both variants are token-compatible for ordinary names and numeric fields,
//! so the hot path scans whitespace-delimited fields without copying. COLUMNS
//! records are consumed in their native column order; the shared builder
//! merges duplicates in place and freezes directly into packed CSC storage.

const std = @import("std");
const types = @import("../types.zig");
const Builder = @import("../builder.zig").Builder;
const ModelData = @import("../model_data.zig").ModelData;
const output = @import("../output.zig");

const Section = enum { none, obj_sense, rows, columns, rhs, ranges, bounds, end };

const Parser = struct {
    allocator: std.mem.Allocator,
    builder: Builder,
    row_by_name: std.StringHashMap(usize),
    col_by_name: std.StringHashMap(usize),
    objective_row: ?[]const u8 = null,
    integer_mode: bool = false,
    rhs_set: ?[]const u8 = null,
    ranges_set: ?[]const u8 = null,
    bounds_set: ?[]const u8 = null,
    control: types.ParseControl,

    fn init(allocator: std.mem.Allocator, fallback_name: []const u8, options: types.ReadOptions) Parser {
        var builder = Builder.init(allocator);
        builder.configureLimits(options);
        builder.name = fallback_name;
        return .{
            .allocator = allocator,
            .builder = builder,
            .row_by_name = std.StringHashMap(usize).init(allocator),
            .col_by_name = std.StringHashMap(usize).init(allocator),
            .control = types.ParseControl.init(options),
        };
    }

    fn deinit(self: *Parser) void {
        self.row_by_name.deinit();
        self.col_by_name.deinit();
        self.builder.deinit();
    }

    /// ROWS/COLUMNS/RHS/RANGES/BOUNDS have all resolved their symbolic names
    /// before finalization, so neither hash table needs to overlap final CSC.
    fn releaseNameIndexes(self: *Parser) void {
        self.row_by_name.deinit();
        self.row_by_name = std.StringHashMap(usize).init(self.allocator);
        self.col_by_name.deinit();
        self.col_by_name = std.StringHashMap(usize).init(self.allocator);
    }

    fn column(self: *Parser, name: []const u8) types.IoError!usize {
        if (self.col_by_name.get(name)) |index| return index;
        const index = try self.builder.addColumn(.{ .name = name, .kind = if (self.integer_mode) .integer else .continuous });
        self.col_by_name.put(name, index) catch return error.OutOfMemory;
        return index;
    }
};

pub fn parse(allocator: std.mem.Allocator, input: []const u8, fallback_name: []const u8, options: types.ReadOptions) types.IoError!ModelData {
    if (input.len > options.max_file_bytes) return error.FileTooLarge;
    try options.checkCancelled();
    var parser = Parser.init(allocator, fallback_name, options);
    defer parser.deinit();
    var section: Section = .none;
    var saw_rows = false;
    var saw_columns = false;
    var saw_end = false;
    var lines = std.mem.tokenizeScalar(u8, input, '\n');
    while (lines.next()) |raw| {
        try parser.control.tick();
        if (raw.len > options.max_line_bytes) return error.ResourceLimitExceeded;
        const line = std.mem.trim(u8, raw, " \t\r");
        if (line.len == 0 or line[0] == '*') continue;
        var fields = try tokenizeRecord(allocator, raw, line, section);
        defer fields.deinit(allocator);
        if (fields.items.len == 0) continue;
        const first = fields.items[0];
        // Section names are legal set names in their own data sections. In
        // particular, `RHS ROW VALUE` is the conventional RHS-set spelling;
        // only a header-shaped record may change parser state.
        const header_shaped = fields.items.len == 1 or std.ascii.eqlIgnoreCase(first, "NAME") or
            std.ascii.eqlIgnoreCase(first, "OBJNAME");
        if (header_shaped) if (sectionKeyword(first)) |next| {
            section = next;
            switch (next) {
                .rows => saw_rows = true,
                .columns => saw_columns = true,
                .end => saw_end = true,
                else => {},
            }
            if (std.ascii.eqlIgnoreCase(first, "NAME") and fields.items.len >= 2) parser.builder.name = fields.items[1];
            continue;
        };
        switch (section) {
            .obj_sense => try parseObjectiveSense(&parser, fields.items),
            .rows => try parseRow(&parser, fields.items),
            .columns => try parseColumnRecord(&parser, fields.items),
            .rhs => try parseRhs(&parser, fields.items),
            .ranges => try parseRanges(&parser, fields.items),
            .bounds => try parseBound(&parser, fields.items),
            .end => break,
            else => return error.InvalidSyntax,
        }
    }
    if (!saw_rows or !saw_columns or !saw_end) return error.InvalidSyntax;
    parser.releaseNameIndexes();
    return parser.builder.finishColumnOrdered(options);
}

pub fn write(file: *std.Io.Writer, allocator: std.mem.Allocator, model: types.ModelView, names: output.Names, options: types.WriteOptions) types.IoError!void {
    _ = options;
    const model_name: []const u8 = if (model.name.len == 0) "MODEL" else model.name;
    if (!tokenName(model_name)) return error.InvalidName;
    for (names.columns) |name| if (!tokenName(name)) return error.InvalidName;
    for (names.rows) |name| if (!tokenName(name)) return error.InvalidName;
    var objective_name: ?[]const u8 = null;
    for ([_][]const u8{ "OBJ", "OBJROW", "_OBJROW" }) |candidate| {
        var collision = false;
        for (names.rows) |name| if (std.mem.eql(u8, name, candidate)) {
            collision = true;
            break;
        };
        if (!collision) {
            objective_name = candidate;
            break;
        }
    }
    const objective = objective_name orelse return error.DuplicateName;
    try output.print(file, allocator, "NAME {s}\nOBJSENSE\n {s}\nROWS\n N {s}\n", .{ model_name, if (model.objective_sense == .maximize) "MAX" else "MIN", objective });
    for (0..model.numRows()) |row| {
        const code: u8 = if (model.row_lower[row] == model.row_upper[row]) 'E' else if (std.math.isInf(model.row_lower[row])) 'L' else 'G';
        const name = names.rows[row];
        try output.print(file, allocator, " {c} {s}\n", .{ code, name });
    }
    try output.write(file, "COLUMNS\n");
    var integer_mode = false;
    var marker_index: usize = 0;
    for (0..model.numCols()) |column| {
        const integer = model.col_type[column] == .integer or model.col_type[column] == .binary or model.col_type[column] == .semi_integer;
        if (integer != integer_mode) {
            try output.print(file, allocator, " MARK{d} 'MARKER' '{s}'\n", .{ marker_index, if (integer) "INTORG" else "INTEND" });
            marker_index += 1;
            integer_mode = integer;
        }
        const begin = model.matrix.col_starts[column];
        const end = model.matrix.col_starts[column + 1];
        var wrote = false;
        if (model.col_cost[column] != 0.0) {
            try output.print(file, allocator, " {s} {s} {d}\n", .{ names.columns[column], objective, model.col_cost[column] });
            wrote = true;
        }
        for (begin..end) |position| {
            const row = model.matrix.row_indices[position].toUsize();
            try output.print(file, allocator, " {s} {s} {d}\n", .{ names.columns[column], names.rows[row], model.matrix.values[position] });
            wrote = true;
        }
        if (!wrote) try output.print(file, allocator, " {s} {s} 0\n", .{ names.columns[column], objective });
    }
    if (integer_mode) try output.print(file, allocator, " MARK{d} 'MARKER' 'INTEND'\n", .{marker_index});
    try output.write(file, "RHS\n");
    for (0..model.numRows()) |row| {
        const rhs = if (std.math.isInf(model.row_lower[row])) model.row_upper[row] else model.row_lower[row];
        if (rhs != 0.0) try output.print(file, allocator, " RHS1 {s} {d}\n", .{ names.rows[row], rhs });
    }
    if (model.objective_offset != 0.0) try output.print(file, allocator, " RHS1 {s} {d}\n", .{ objective, -model.objective_offset });
    var has_ranges = false;
    for (0..model.numRows()) |row| {
        if (!std.math.isInf(model.row_lower[row]) and !std.math.isInf(model.row_upper[row]) and model.row_lower[row] != model.row_upper[row]) {
            if (!has_ranges) {
                try output.write(file, "RANGES\n");
                has_ranges = true;
            }
            try output.print(file, allocator, " RNG1 {s} {d}\n", .{ names.rows[row], model.row_upper[row] - model.row_lower[row] });
        }
    }
    try output.write(file, "BOUNDS\n");
    for (0..model.numCols()) |column| try writeMpsBound(file, allocator, names.columns[column], model.col_type[column], model.col_lower[column], model.col_upper[column]);
    try output.write(file, "ENDATA\n");
}

fn tokenName(name: []const u8) bool {
    if (name.len == 0 or name.len > 255) return false;
    for (name) |char| if (std.ascii.isWhitespace(char)) return false;
    return true;
}

fn writeMpsBound(file: *std.Io.Writer, allocator: std.mem.Allocator, name: []const u8, kind: types.VariableType, lower: f64, upper: f64) types.IoError!void {
    if (lower > upper or std.math.isNan(lower) or std.math.isNan(upper)) return error.InvalidBounds;
    if (kind == .binary) {
        try output.print(file, allocator, " BV BND1 {s}\n", .{name});
        return;
    }
    if (std.math.isInf(lower) and lower < 0 and std.math.isInf(upper) and upper > 0) {
        try output.print(file, allocator, " FR BND1 {s}\n", .{name});
        return;
    }
    if (lower == upper) {
        try output.print(file, allocator, " FX BND1 {s} {d}\n", .{ name, lower });
        return;
    }
    const lower_code: []const u8 = if (kind == .integer) "LI" else "LO";
    const upper_code: []const u8 = if (kind == .integer) "UI" else if (kind == .semi_continuous) "SC" else if (kind == .semi_integer) "SI" else "UP";
    if (std.math.isInf(lower) and lower < 0) try output.print(file, allocator, " MI BND1 {s}\n", .{name}) else if (lower != 0.0) try output.print(file, allocator, " {s} BND1 {s} {d}\n", .{ lower_code, name, lower });
    if (!std.math.isInf(upper)) try output.print(file, allocator, " {s} BND1 {s} {d}\n", .{ upper_code, name, upper });
}

fn tokenize(allocator: std.mem.Allocator, line: []const u8) types.IoError!std.ArrayListUnmanaged([]const u8) {
    var fields: std.ArrayListUnmanaged([]const u8) = .empty;
    errdefer fields.deinit(allocator);
    var iterator = std.mem.tokenizeAny(u8, line, " \t\r");
    while (iterator.next()) |field| fields.append(allocator, field) catch return error.OutOfMemory;
    return fields;
}

/// Traditional fixed MPS permits blanks inside its eight-character row and
/// column names. Whitespace tokenization silently splits those names, so use
/// the standard field columns whenever the physical record has fixed-format
/// separator blanks. Free-format records fall back to the existing scanner.
fn tokenizeRecord(
    allocator: std.mem.Allocator,
    raw: []const u8,
    trimmed: []const u8,
    section: Section,
) types.IoError!std.ArrayListUnmanaged([]const u8) {
    if (!isFixedRecord(raw, section)) return tokenize(allocator, trimmed);
    var fields: std.ArrayListUnmanaged([]const u8) = .empty;
    errdefer fields.deinit(allocator);
    switch (section) {
        .rows => {
            try appendFixedField(allocator, &fields, raw, 1, 3);
            try appendFixedField(allocator, &fields, raw, 4, 12);
        },
        .columns, .rhs, .ranges => {
            try appendFixedField(allocator, &fields, raw, 4, 12);
            try appendFixedField(allocator, &fields, raw, 14, 22);
            try appendFixedField(allocator, &fields, raw, 24, 36);
            try appendFixedField(allocator, &fields, raw, 39, 47);
            try appendFixedField(allocator, &fields, raw, 49, 61);
        },
        .bounds => {
            try appendFixedField(allocator, &fields, raw, 1, 3);
            try appendFixedField(allocator, &fields, raw, 4, 12);
            try appendFixedField(allocator, &fields, raw, 14, 22);
            try appendFixedField(allocator, &fields, raw, 24, 36);
        },
        else => return tokenize(allocator, trimmed),
    }
    return fields;
}

fn isFixedRecord(raw: []const u8, section: Section) bool {
    if (raw.len < 5 or raw[0] != ' ' or raw[3] != ' ') return false;
    return switch (section) {
        .rows => raw.len <= 12 or raw[12] == ' ',
        .columns, .rhs, .ranges, .bounds => raw.len >= 14 and raw[12] == ' ' and raw[13] == ' ',
        else => false,
    };
}

fn appendFixedField(
    allocator: std.mem.Allocator,
    fields: *std.ArrayListUnmanaged([]const u8),
    raw: []const u8,
    begin: usize,
    end: usize,
) types.IoError!void {
    if (begin >= raw.len) return;
    const field = std.mem.trim(u8, raw[begin..@min(end, raw.len)], " \t\r");
    if (field.len != 0) fields.append(allocator, field) catch return error.OutOfMemory;
}

fn sectionKeyword(field: []const u8) ?Section {
    if (std.ascii.eqlIgnoreCase(field, "NAME")) return .none;
    if (std.ascii.eqlIgnoreCase(field, "OBJSENSE")) return .obj_sense;
    if (std.ascii.eqlIgnoreCase(field, "OBJNAME")) return .none;
    if (std.ascii.eqlIgnoreCase(field, "ROWS")) return .rows;
    if (std.ascii.eqlIgnoreCase(field, "COLUMNS")) return .columns;
    if (std.ascii.eqlIgnoreCase(field, "RHS")) return .rhs;
    if (std.ascii.eqlIgnoreCase(field, "RANGES")) return .ranges;
    if (std.ascii.eqlIgnoreCase(field, "BOUNDS")) return .bounds;
    if (std.ascii.eqlIgnoreCase(field, "ENDATA")) return .end;
    return null;
}

fn parseObjectiveSense(parser: *Parser, fields: []const []const u8) types.IoError!void {
    if (fields.len != 1) return error.InvalidSyntax;
    if (std.ascii.startsWithIgnoreCase(fields[0], "MAX")) parser.builder.objective_sense = .maximize else if (std.ascii.startsWithIgnoreCase(fields[0], "MIN")) parser.builder.objective_sense = .minimize else return error.InvalidSyntax;
}

fn parseRow(parser: *Parser, fields: []const []const u8) types.IoError!void {
    if (fields.len != 2 or fields[0].len != 1) return error.InvalidSyntax;
    const name = fields[1];
    switch (std.ascii.toUpper(fields[0][0])) {
        'N' => if (parser.objective_row == null) {
            parser.objective_row = name;
        },
        'L', 'G', 'E' => |code| {
            if (parser.row_by_name.contains(name)) return error.DuplicateName;
            const index = try parser.builder.addRow(if (code == 'L')
                .{ .name = name, .upper = 0.0 }
            else if (code == 'G')
                .{ .name = name, .lower = 0.0 }
            else
                .{ .name = name, .lower = 0.0, .upper = 0.0 });
            parser.row_by_name.put(name, index) catch return error.OutOfMemory;
        },
        else => return error.InvalidSyntax,
    }
}

fn unquoteMarker(value: []const u8) []const u8 {
    if (value.len >= 2 and ((value[0] == '\'' and value[value.len - 1] == '\'') or (value[0] == '"' and value[value.len - 1] == '"'))) return value[1 .. value.len - 1];
    return value;
}

fn parseColumnRecord(parser: *Parser, fields: []const []const u8) types.IoError!void {
    if (fields.len >= 3 and std.ascii.eqlIgnoreCase(unquoteMarker(fields[1]), "MARKER")) {
        const marker = unquoteMarker(fields[2]);
        if (std.ascii.eqlIgnoreCase(marker, "INTORG")) parser.integer_mode = true else if (std.ascii.eqlIgnoreCase(marker, "INTEND")) parser.integer_mode = false else return error.InvalidSyntax;
        return;
    }
    if (fields.len != 3 and fields.len != 5) return error.InvalidSyntax;
    const column = try parser.column(fields[0]);
    if (parser.integer_mode and parser.builder.columns.items[column].kind == .continuous) parser.builder.columns.items[column].kind = .integer;
    try addColumnPair(parser, column, fields[1], fields[2]);
    if (fields.len == 5) try addColumnPair(parser, column, fields[3], fields[4]);
}

fn addColumnPair(parser: *Parser, column: usize, row_name: []const u8, number: []const u8) types.IoError!void {
    const value = try parseFinite(number);
    if (parser.objective_row != null and std.mem.eql(u8, row_name, parser.objective_row.?)) {
        parser.builder.columns.items[column].cost += value;
    } else {
        const row = parser.row_by_name.get(row_name) orelse return error.UnknownName;
        try parser.builder.addTerm(row, column, value);
    }
}

fn parseRhs(parser: *Parser, fields: []const []const u8) types.IoError!void {
    const offset: usize = switch (fields.len) {
        2, 4 => 0,
        3, 5 => 1,
        else => return error.InvalidSyntax,
    };
    if (offset == 1) {
        if (parser.rhs_set == null) parser.rhs_set = fields[0];
        if (!std.mem.eql(u8, parser.rhs_set.?, fields[0])) return;
    }
    try setRhsPair(parser, fields[offset], fields[offset + 1]);
    if (fields.len - offset == 4) try setRhsPair(parser, fields[offset + 2], fields[offset + 3]);
}

fn setRhsPair(parser: *Parser, row_name: []const u8, number: []const u8) types.IoError!void {
    const value = try parseFinite(number);
    if (parser.objective_row != null and std.mem.eql(u8, row_name, parser.objective_row.?)) {
        parser.builder.objective_offset = -value;
    } else {
        const row = parser.row_by_name.get(row_name) orelse return error.UnknownName;
        const target = &parser.builder.rows.items[row];
        if (!std.math.isInf(target.lower) and !std.math.isInf(target.upper)) {
            target.lower = value;
            target.upper = value;
        } else if (std.math.isInf(target.lower)) {
            target.upper = value;
        } else {
            target.lower = value;
        }
    }
}

fn parseRanges(parser: *Parser, fields: []const []const u8) types.IoError!void {
    const offset: usize = switch (fields.len) {
        2, 4 => 0,
        3, 5 => 1,
        else => return error.InvalidSyntax,
    };
    if (offset == 1) {
        if (parser.ranges_set == null) parser.ranges_set = fields[0];
        if (!std.mem.eql(u8, parser.ranges_set.?, fields[0])) return;
    }
    try setRangePair(parser, fields[offset], fields[offset + 1]);
    if (fields.len - offset == 4) try setRangePair(parser, fields[offset + 2], fields[offset + 3]);
}

fn setRangePair(parser: *Parser, row_name: []const u8, number: []const u8) types.IoError!void {
    const signed_range = try parseFinite(number);
    const magnitude = @abs(signed_range);
    const row = parser.row_by_name.get(row_name) orelse return error.UnknownName;
    const target = &parser.builder.rows.items[row];
    if (std.math.isInf(target.lower)) {
        target.lower = target.upper - magnitude;
    } else if (std.math.isInf(target.upper)) {
        target.upper = target.lower + magnitude;
    } else if (target.lower == target.upper) {
        if (signed_range >= 0.0) target.upper += magnitude else target.lower -= magnitude;
    } else {
        return error.InvalidSyntax;
    }
}

test "MPS name indexes release all hash capacity before finalization" {
    var parser = Parser.init(std.testing.allocator, "names", .{});
    defer parser.deinit();
    parser.row_by_name.put("row", 0) catch return error.OutOfMemory;
    _ = try parser.column("column");
    try std.testing.expect(parser.row_by_name.capacity() > 0);
    try std.testing.expect(parser.col_by_name.capacity() > 0);

    parser.releaseNameIndexes();
    try std.testing.expectEqual(@as(u32, 0), parser.row_by_name.count());
    try std.testing.expectEqual(@as(u32, 0), parser.row_by_name.capacity());
    try std.testing.expectEqual(@as(u32, 0), parser.col_by_name.count());
    try std.testing.expectEqual(@as(u32, 0), parser.col_by_name.capacity());
}

test "MPS parser enforces record and semantic limits and cancellation" {
    const source =
        \\NAME LIMITS
        \\ROWS
        \\ N OBJ
        \\ L ROW
        \\COLUMNS
        \\ X ROW 1
        \\RHS
        \\ RHS1 ROW 1
        \\ENDATA
    ;
    try std.testing.expectError(error.ResourceLimitExceeded, parse(std.testing.allocator, source, "limits", .{ .max_line_bytes = 5 }));
    try std.testing.expectError(error.ResourceLimitExceeded, parse(std.testing.allocator, source, "limits", .{ .max_rows = 0 }));
    try std.testing.expectError(error.ResourceLimitExceeded, parse(std.testing.allocator, source, "limits", .{ .max_columns = 0 }));
    try std.testing.expectError(error.ResourceLimitExceeded, parse(std.testing.allocator, source, "limits", .{ .max_matrix_terms = 0 }));
    var interrupted = std.atomic.Value(bool).init(true);
    try std.testing.expectError(error.Cancelled, parse(std.testing.allocator, source, "limits", .{ .interrupt_flag = &interrupted }));
}

fn parseBound(parser: *Parser, fields: []const []const u8) types.IoError!void {
    const code = fields[0];
    const requires_value = std.ascii.eqlIgnoreCase(code, "LO") or std.ascii.eqlIgnoreCase(code, "UP") or
        std.ascii.eqlIgnoreCase(code, "FX") or std.ascii.eqlIgnoreCase(code, "LI") or
        std.ascii.eqlIgnoreCase(code, "UI") or std.ascii.eqlIgnoreCase(code, "SC") or
        std.ascii.eqlIgnoreCase(code, "SI");
    const omitted_set = (requires_value and fields.len == 3) or (!requires_value and fields.len == 2);
    if ((requires_value and fields.len != 3 and fields.len != 4) or
        (!requires_value and fields.len < 2 or fields.len > 4))
        return error.InvalidSyntax;
    const set_name: []const u8 = if (omitted_set) "" else fields[1];
    const column_name = fields[if (omitted_set) 1 else 2];
    const value_index: ?usize = if (requires_value)
        (if (omitted_set) 2 else 3)
    else if ((!omitted_set and fields.len == 4) or (omitted_set and fields.len == 3))
        fields.len - 1
    else
        null;
    if (parser.bounds_set == null) parser.bounds_set = set_name;
    if (!std.mem.eql(u8, parser.bounds_set.?, set_name)) return;
    const column = try parser.column(column_name);
    const value = if (value_index) |index| try parseFinite(fields[index]) else 0.0;
    const target = &parser.builder.columns.items[column];
    if (std.ascii.eqlIgnoreCase(code, "LO")) target.lower = value else if (std.ascii.eqlIgnoreCase(code, "UP")) target.upper = value else if (std.ascii.eqlIgnoreCase(code, "FX")) {
        target.lower = value;
        target.upper = value;
    } else if (std.ascii.eqlIgnoreCase(code, "FR")) {
        target.lower = -std.math.inf(f64);
        target.upper = std.math.inf(f64);
    } else if (std.ascii.eqlIgnoreCase(code, "MI")) target.lower = -std.math.inf(f64) else if (std.ascii.eqlIgnoreCase(code, "PL")) target.upper = std.math.inf(f64) else if (std.ascii.eqlIgnoreCase(code, "BV")) {
        target.lower = 0.0;
        target.upper = 1.0;
        target.kind = .binary;
    } else if (std.ascii.eqlIgnoreCase(code, "LI")) {
        target.lower = value;
        target.kind = .integer;
    } else if (std.ascii.eqlIgnoreCase(code, "UI")) {
        target.upper = value;
        target.kind = .integer;
    } else if (std.ascii.eqlIgnoreCase(code, "SC")) {
        target.upper = value;
        target.kind = .semi_continuous;
    } else if (std.ascii.eqlIgnoreCase(code, "SI")) {
        target.upper = value;
        target.kind = .semi_integer;
    } else return error.UnsupportedFeature;
    if (target.lower > target.upper) return error.InvalidBounds;
}

fn parseFinite(text: []const u8) types.IoError!f64 {
    const value = std.fmt.parseFloat(f64, text) catch return error.InvalidNumber;
    if (!std.math.isFinite(value)) return error.NonFiniteValue;
    return value;
}

test "parses free MPS objective rows integer markers RHS and bounds" {
    const input =
        \\NAME SAMPLE
        \\OBJSENSE
        \\ MAX
        \\ROWS
        \\ N OBJ
        \\ G DEMAND
        \\COLUMNS
        \\ X OBJ -2 DEMAND 3
        \\ MARK0000 'MARKER' 'INTORG'
        \\ Y OBJ 1 DEMAND -1
        \\ MARK0001 'MARKER' 'INTEND'
        \\RHS
        \\ RHS1 DEMAND -7
        \\BOUNDS
        \\ UP BND X 4
        \\ BV BND Y
        \\ENDATA
    ;
    var model = try parse(std.testing.allocator, input, "fallback", .{});
    defer model.deinit();
    try std.testing.expectEqualStrings("SAMPLE", model.name);
    try std.testing.expectEqual(types.ObjectiveSense.maximize, model.objective_sense);
    try std.testing.expectEqual(@as(usize, 2), model.col_cost.len);
    try std.testing.expectEqual(@as(usize, 2), model.matrix.nnz());
    try std.testing.expectEqual(types.VariableType.binary, model.col_type[1]);
    try std.testing.expectEqual(@as(f64, -7.0), model.row_lower[0]);
}

test "section keywords remain valid RHS and RANGES set names" {
    const input =
        \\NAME KEYWORD_SETS
        \\ROWS
        \\ N OBJ
        \\ G DEMAND
        \\COLUMNS
        \\ X OBJ 1 DEMAND 1
        \\RHS
        \\ RHS DEMAND 7
        \\RANGES
        \\ RANGES DEMAND 2
        \\ENDATA
    ;
    var model = try parse(std.testing.allocator, input, "fallback", .{});
    defer model.deinit();
    try std.testing.expectEqual(@as(f64, 7.0), model.row_lower[0]);
    try std.testing.expectEqual(@as(f64, 9.0), model.row_upper[0]);
}

test "fixed MPS accepts an omitted RHS set name" {
    const input =
        \\NAME FIXED_RHS
        \\ROWS
        \\ N OBJ
        \\ G FIRST
        \\ L SECOND
        \\COLUMNS
        \\ X OBJ 1 FIRST 1
        \\ X SECOND 1
        \\RHS
        \\ FIRST 7 SECOND 9
        \\ENDATA
    ;
    var model = try parse(std.testing.allocator, input, "fallback", .{});
    defer model.deinit();
    try std.testing.expectEqual(@as(f64, 7.0), model.row_lower[0]);
    try std.testing.expectEqual(@as(f64, 9.0), model.row_upper[1]);
}

test "fixed MPS preserves blanks inside eight-character names" {
    const input =
        \\NAME          FIXEDNAMES
        \\ROWS
        \\ N  OB1PNW20
        \\ E  DEDO3 1R
        \\COLUMNS
        \\    DEDO3 11  OB1PNW20        .02466   DEDO3 1R           -1.
        \\RHS
        \\    RHS 1     DEDO3 1R            3.
        \\BOUNDS
        \\ UP BND-1     DEDO3 11       200000.
        \\ENDATA
    ;
    var model = try parse(std.testing.allocator, input, "fixed-names", .{});
    defer model.deinit();
    try std.testing.expectEqual(@as(usize, 1), model.col_cost.len);
    try std.testing.expectEqual(@as(usize, 1), model.row_lower.len);
    try std.testing.expectEqual(@as(f64, 0.02466), model.col_cost[0]);
    try std.testing.expectEqual(@as(f64, 3.0), model.row_lower[0]);
    try std.testing.expectEqual(@as(f64, 3.0), model.row_upper[0]);
    try std.testing.expectEqual(@as(f64, 200000.0), model.col_upper[0]);
}

test "fixed MPS accepts an omitted bounds set name" {
    const input =
        \\NAME OMITTED_BOUND_SET
        \\ROWS
        \\ N OBJ
        \\ L ROW
        \\COLUMNS
        \\ X OBJ 1 ROW 1
        \\BOUNDS
        \\ UP           X                   7.
        \\ENDATA
    ;
    var model = try parse(std.testing.allocator, input, "omitted-bound-set", .{});
    defer model.deinit();
    try std.testing.expectEqual(@as(f64, 7.0), model.col_upper[0]);
}
