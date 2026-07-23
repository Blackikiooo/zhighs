//! CPLEX LP-compatible text frontend.
//!
//! The parser is allocation-conscious: tokens are borrowed slices into the
//! input buffer, names are interned by a hash table, and only semantic arrays
//! survive `finish`. Expressions are parsed directly without constructing an
//! AST.

const std = @import("std");
const types = @import("../types.zig");
const Builder = @import("../builder.zig").Builder;
const ModelData = @import("../model_data.zig").ModelData;
const output = @import("../output.zig");
const lexer_module = @import("lexer.zig");
const Token = lexer_module.Token;
const Tag = lexer_module.Tag;

pub const Lexer = lexer_module.Lexer;
pub const LexerError = lexer_module.Error;
pub const token = @import("token.zig");

const Section = enum { none, objective, constraints, bounds, binaries, generals, semicont, semiint, end };

const Parser = struct {
    /// Shared semantic output builder and temporary term storage.
    builder: Builder,
    /// Borrowed-name to temporary structural-column index.
    variables: std.StringHashMap(usize),
    /// Cooperative cancellation state.
    control: types.ParseControl,
    /// Objective scalar held until the next term/header disambiguates it.
    pending_objective_coefficient: ?f64 = null,
    /// Sign applied while parsing the current objective sense.
    objective_sign: f64 = 1.0,

    /// Initialize parser indexes and the LP column-chain builder fast path.
    fn init(allocator: std.mem.Allocator, name: []const u8, options: types.ReadOptions) Parser {
        var builder = Builder.init(allocator);
        builder.configureLimits(options);
        builder.enableColumnTermStorage();
        builder.name = name;
        return .{ .builder = builder, .variables = std.StringHashMap(usize).init(allocator), .control = types.ParseControl.init(options) };
    }

    /// Release name indexes and all unfinished semantic state.
    fn deinit(self: *Parser) void {
        self.variables.deinit();
        self.builder.deinit();
    }

    /// Name lookup ends with syntax parsing. Releasing the table before CSC
    /// finalization prevents it from overlapping the final matrix allocation.
    fn releaseVariableIndex(self: *Parser) void {
        self.variables.deinit();
        self.variables = std.StringHashMap(usize).init(self.builder.allocator);
    }

    /// Resolve or create one LP variable symbol.
    fn variable(self: *Parser, name: []const u8) types.IoError!usize {
        if (!validName(name)) return error.InvalidName;
        if (self.variables.get(name)) |index| return index;
        const index = try self.builder.addColumn(.{ .name = name });
        self.variables.put(name, index) catch return error.OutOfMemory;
        return index;
    }
};

/// Parse a complete borrowed CPLEX-LP source into owned canonical model data.
pub fn parse(allocator: std.mem.Allocator, input: []const u8, model_name: []const u8, options: types.ReadOptions) types.IoError!ModelData {
    if (input.len > options.max_file_bytes) return error.FileTooLarge;
    try options.checkCancelled();
    var parser = Parser.init(allocator, model_name, options);
    defer parser.deinit();
    var section: Section = .none;
    var saw_objective = false;
    var saw_end = false;
    var lexer = Lexer.initWithLimits(input, options.max_line_bytes, options.max_token_bytes);
    while (true) {
        try parser.control.tick();
        const line_start = lexer;
        const first = try nextToken(&lexer);
        switch (first.tag) {
            .newline => continue,
            .eof => break,
            else => {},
        }

        if (try consumeHeaderLine(&lexer, first)) |next| {
            if (section == .objective and next.section != .objective) {
                if (parser.pending_objective_coefficient) |constant| {
                    parser.builder.objective_offset += constant;
                    parser.pending_objective_coefficient = null;
                }
            }
            section = next.section;
            if (next.sense) |sense| {
                parser.builder.objective_sense = sense;
                saw_objective = true;
            }
            if (section == .end) saw_end = true;
            continue;
        }

        // A non-header line is parsed from its first token. Restoring this
        // cheap value checkpoint avoids a token buffer while header lines
        // retain the already-consumed first token.
        lexer = line_start;
        switch (section) {
            .objective => try parseObjective(&parser, &lexer),
            .constraints => try parseConstraint(&parser, &lexer),
            .bounds => try parseBound(&parser, &lexer),
            .binaries => try parseTypes(&parser, &lexer, .binary),
            .generals => try parseTypes(&parser, &lexer, .integer),
            .semicont => try parseTypes(&parser, &lexer, .semi_continuous),
            .semiint => try parseTypes(&parser, &lexer, .semi_integer),
            else => return error.InvalidSyntax,
        }
    }
    if (!saw_objective or !saw_end) return error.InvalidSyntax;
    parser.releaseVariableIndex();
    return parser.builder.finish(options);
}

/// Serialize a borrowed linear model in deterministic CPLEX-LP syntax.
pub fn write(file: *std.Io.Writer, allocator: std.mem.Allocator, model: types.ModelView, names: output.Names, options: types.WriteOptions) types.IoError!void {
    _ = options;
    for (names.columns) |name| if (!validName(name)) return error.InvalidName;
    for (names.rows) |name| if (!validName(name)) return error.InvalidName;
    try output.write(file, "\\ LP format generated by zhighs\n");
    try output.write(file, if (model.objective_sense == .minimize) "Minimize\n obj:" else "Maximize\n obj:");
    var first = true;
    for (model.col_cost, names.columns) |coefficient, name| try writeTerm(file, allocator, coefficient, name, &first);
    if (model.objective_offset != 0.0) {
        if (first)
            try output.print(file, allocator, " {d}", .{model.objective_offset})
        else if (model.objective_offset < 0.0)
            try output.print(file, allocator, " - {d}", .{@abs(model.objective_offset)})
        else
            try output.print(file, allocator, " + {d}", .{model.objective_offset});
        first = false;
    }
    if (first) try output.write(file, " 0");
    try output.write(file, "\nSubject To\n");

    const row_starts = allocator.alloc(usize, model.numRows() + 1) catch return error.OutOfMemory;
    defer allocator.free(row_starts);
    @memset(row_starts, 0);
    for (model.matrix.row_indices) |row| row_starts[row.toUsize() + 1] += 1;
    for (0..model.numRows()) |row| row_starts[row + 1] += row_starts[row];
    const row_columns = allocator.alloc(usize, model.matrix.nnz()) catch return error.OutOfMemory;
    defer allocator.free(row_columns);
    const row_values = allocator.alloc(f64, model.matrix.nnz()) catch return error.OutOfMemory;
    defer allocator.free(row_values);
    const cursor = allocator.dupe(usize, row_starts[0..model.numRows()]) catch return error.OutOfMemory;
    defer allocator.free(cursor);
    for (0..model.numCols()) |column| {
        for (model.matrix.col_starts[column]..model.matrix.col_starts[column + 1]) |position| {
            const row = model.matrix.row_indices[position].toUsize();
            const destination = cursor[row];
            cursor[row] += 1;
            row_columns[destination] = column;
            row_values[destination] = model.matrix.values[position];
        }
    }
    for (0..model.numRows()) |row| {
        try output.print(file, allocator, " {s}:", .{names.rows[row]});
        const lower = model.row_lower[row];
        const upper = model.row_upper[row];
        const ranged = !std.math.isInf(lower) and !std.math.isInf(upper) and lower != upper;
        if (ranged) try output.print(file, allocator, " {d} <=", .{lower});
        first = true;
        for (row_starts[row]..row_starts[row + 1]) |position| try writeTerm(file, allocator, row_values[position], names.columns[row_columns[position]], &first);
        if (first) try output.write(file, " 0");
        if (lower == upper)
            try output.print(file, allocator, " = {d}\n", .{lower})
        else if (std.math.isInf(lower))
            try output.print(file, allocator, " <= {d}\n", .{upper})
        else if (std.math.isInf(upper))
            try output.print(file, allocator, " >= {d}\n", .{lower})
        else
            try output.print(file, allocator, " <= {d}\n", .{upper});
    }
    try output.write(file, "Bounds\n");
    for (0..model.numCols()) |column| try writeBound(file, allocator, names.columns[column], model.col_lower[column], model.col_upper[column]);
    try writeTypes(file, allocator, model, names.columns, .binary, "Binaries");
    try writeTypes(file, allocator, model, names.columns, .integer, "Generals");
    try writeTypes(file, allocator, model, names.columns, .semi_continuous, "Semi-Continuous");
    try writeTypes(file, allocator, model, names.columns, .semi_integer, "Semi-Integer");
    try output.write(file, "End\n");
}

/// Emit one signed nonzero linear-expression term with stable spacing.
fn writeTerm(file: *std.Io.Writer, allocator: std.mem.Allocator, coefficient: f64, name: []const u8, first: *bool) types.IoError!void {
    if (coefficient == 0.0) return;
    if (!std.math.isFinite(coefficient)) return error.NonFiniteValue;
    if (first.*) try output.write(file, if (coefficient < 0.0) " - " else " ") else try output.write(file, if (coefficient < 0.0) " - " else " + ");
    const magnitude = @abs(coefficient);
    if (magnitude != 1.0) try output.print(file, allocator, "{d} ", .{magnitude});
    try output.write(file, name);
    first.* = false;
}

/// Emit the shortest LP bound statement representing one variable interval.
fn writeBound(file: *std.Io.Writer, allocator: std.mem.Allocator, name: []const u8, lower: f64, upper: f64) types.IoError!void {
    if (std.math.isNan(lower) or std.math.isNan(upper) or lower > upper) return error.InvalidBounds;
    if (std.math.isInf(lower) and lower < 0 and std.math.isInf(upper) and upper > 0)
        try output.print(file, allocator, " {s} free\n", .{name})
    else if (lower == upper)
        try output.print(file, allocator, " {s} = {d}\n", .{ name, lower })
    else if (!std.math.isInf(lower) and !std.math.isInf(upper))
        try output.print(file, allocator, " {d} <= {s} <= {d}\n", .{ lower, name, upper })
    else if (!std.math.isInf(lower))
        try output.print(file, allocator, " {s} >= {d}\n", .{ name, lower })
    else
        try output.print(file, allocator, " {s} <= {d}\n", .{ name, upper });
}

/// Emit one optional variable-domain section for matching columns.
fn writeTypes(file: *std.Io.Writer, allocator: std.mem.Allocator, model: types.ModelView, names: []const []u8, kind: types.VariableType, title: []const u8) types.IoError!void {
    var any = false;
    for (model.col_type) |actual| if (actual == kind) {
        any = true;
        break;
    };
    if (!any) return;
    try output.print(file, allocator, "{s}\n", .{title});
    for (model.col_type, names) |actual, name| if (actual == kind) try output.print(file, allocator, " {s}\n", .{name});
}

/// Result of recognizing one LP section header.
const Header = struct {
    /// Parser state entered after consuming the header.
    section: Section,
    /// Objective direction carried only by minimize/maximize headers.
    sense: ?types.ObjectiveSense = null,
};

/// Recognize a section header only at the parser's current line start. The
/// checkpoint is committed only when the complete line matches a header, so a
/// normal expression beginning with a keyword-like identifier is untouched.
fn consumeHeaderLine(lexer: *Lexer, first: Token) types.IoError!?Header {
    var probe = lexer.*;
    if (first.tag != .identifier) return null;

    var result: ?Header = null;
    if (equalsAny(first.lexeme, &.{ "minimize", "minimum", "min" })) {
        result = .{ .section = .objective, .sense = .minimize };
    } else if (equalsAny(first.lexeme, &.{ "maximize", "maximum", "max" })) {
        result = .{ .section = .objective, .sense = .maximize };
    } else if (std.ascii.eqlIgnoreCase(first.lexeme, "subject")) {
        const second = try nextToken(&probe);
        if (second.tag == .identifier and std.ascii.eqlIgnoreCase(second.lexeme, "to")) result = .{ .section = .constraints };
    } else if (std.ascii.eqlIgnoreCase(first.lexeme, "such")) {
        const second = try nextToken(&probe);
        if (second.tag == .identifier and std.ascii.eqlIgnoreCase(second.lexeme, "that")) result = .{ .section = .constraints };
    } else if (equalsAny(first.lexeme, &.{ "st", "s.t." })) {
        result = .{ .section = .constraints };
    } else if (std.ascii.eqlIgnoreCase(first.lexeme, "bounds")) {
        result = .{ .section = .bounds };
    } else if (equalsAny(first.lexeme, &.{ "binary", "binaries", "bin" })) {
        result = .{ .section = .binaries };
    } else if (equalsAny(first.lexeme, &.{ "general", "generals", "gen" })) {
        result = .{ .section = .generals };
    } else if (std.ascii.eqlIgnoreCase(first.lexeme, "semis")) {
        result = .{ .section = .semicont };
    } else if (std.ascii.eqlIgnoreCase(first.lexeme, "semi")) {
        // `Semi` is itself a common alias, while `Semi-Continuous` and
        // `Semi-Integer` are tokenized as identifier/minus/identifier. Probe
        // the suffix independently so a failed suffix leaves the lexer just
        // after `Semi` for the normal line-terminator check below.
        var suffix = probe;
        if (try consumeHyphenatedSemiHeader(&suffix)) |header_value| {
            probe = suffix;
            result = header_value;
        } else {
            result = .{ .section = .semicont };
        }
    } else if (std.ascii.eqlIgnoreCase(first.lexeme, "end")) {
        result = .{ .section = .end };
    }
    const header_value = result orelse return null;
    const terminator = try nextToken(&probe);
    if (terminator.tag != .newline and terminator.tag != .eof) return null;
    lexer.* = probe;
    return header_value;
}

/// Recognize the two-token `semi-continuous`/`semi-integer` header spellings.
fn consumeHyphenatedSemiHeader(lexer: *Lexer) types.IoError!?Header {
    const hyphen = try nextToken(lexer);
    if (hyphen.tag != .minus) return null;
    const kind = try nextToken(lexer);
    if (kind.tag != .identifier) return null;
    if (std.ascii.eqlIgnoreCase(kind.lexeme, "continuous")) return .{ .section = .semicont };
    if (std.ascii.eqlIgnoreCase(kind.lexeme, "integer")) return .{ .section = .semiint };
    return null;
}

/// Case-insensitively compare one spelling against a keyword set.
fn equalsAny(value: []const u8, candidates: []const []const u8) bool {
    for (candidates) |candidate| if (std.ascii.eqlIgnoreCase(value, candidate)) return true;
    return false;
}

/// Return whether a symbol can be emitted unquoted by this LP frontend.
fn validName(name: []const u8) bool {
    if (name.len == 0 or name.len > 255 or std.ascii.isDigit(name[0])) return false;
    if (std.mem.indexOfScalar(u8, "<>=()[],", name[0]) != null) return false;
    for (name) |char| if (std.ascii.isWhitespace(char) or std.mem.indexOfScalar(u8, "+-*^:\\", char) != null) return false;
    const reserved = [_][]const u8{ "min", "minimum", "minimize", "max", "maximum", "maximize", "st", "bounds", "binary", "binaries", "general", "generals", "end", "free" };
    for (reserved) |word| if (std.ascii.eqlIgnoreCase(name, word)) return false;
    return true;
}

/// Convert a numeric token and reject NaN or infinity.
fn parseFinite(text: []const u8) types.IoError!f64 {
    const value = std.fmt.parseFloat(f64, text) catch return error.InvalidNumber;
    if (!std.math.isFinite(value)) return error.NonFiniteValue;
    return value;
}

/// Read one token while mapping lexer errors to stable public I/O errors.
fn nextToken(lexer: *Lexer) types.IoError!Token {
    return lexer.next() catch |err| return switch (err) {
        error.ResourceLimitExceeded => error.ResourceLimitExceeded,
        error.InvalidCharacter => error.InvalidSyntax,
    };
}

/// Inspect one token without advancing parser state.
fn peekToken(lexer: *Lexer) types.IoError!Token {
    return lexer.peek() catch |err| return switch (err) {
        error.ResourceLimitExceeded => error.ResourceLimitExceeded,
        error.InvalidCharacter => error.InvalidSyntax,
    };
}

/// Consume one token and require the requested lexical tag.
fn expect(lexer: *Lexer, expected: Tag) types.IoError!Token {
    const found = try nextToken(lexer);
    if (found.tag != expected) return error.InvalidSyntax;
    return found;
}

/// Require newline or EOF after a complete LP statement/header.
fn expectLineEnd(lexer: *Lexer) types.IoError!void {
    const terminator = try nextToken(lexer);
    if (terminator.tag != .newline and terminator.tag != .eof) return error.InvalidSyntax;
}

/// Consume an optional `identifier:` prefix and return its borrowed name.
fn consumeOptionalLabel(lexer: *Lexer) types.IoError!?[]const u8 {
    var probe = lexer.*;
    const name = try nextToken(&probe);
    if (name.tag != .identifier) return null;
    const separator = try nextToken(&probe);
    if (separator.tag != .colon) return null;
    if (!validName(name.lexeme)) return error.InvalidName;
    lexer.* = probe;
    return name.lexeme;
}

/// Parse one objective line directly into column costs and constant offset.
fn parseObjective(parser: *Parser, lexer: *Lexer) types.IoError!void {
    _ = try consumeOptionalLabel(lexer);
    if (parser.pending_objective_coefficient) |coefficient| {
        const name = try nextToken(lexer);
        if (name.tag == .newline or name.tag == .eof) return;
        if (name.tag != .identifier) return error.InvalidSyntax;
        const column = try parser.variable(name.lexeme);
        parser.builder.columns.items[column].cost += coefficient;
        parser.pending_objective_coefficient = null;
    }
    while (true) {
        try parser.control.tick();
        const current = try nextToken(lexer);
        switch (current.tag) {
            .plus => parser.objective_sign = 1.0,
            .minus => parser.objective_sign = -1.0,
            .identifier => {
                const column = try parser.variable(current.lexeme);
                parser.builder.columns.items[column].cost += parser.objective_sign;
                parser.objective_sign = 1.0;
            },
            .number => {
                const coefficient = parser.objective_sign * try parseFinite(current.lexeme);
                parser.objective_sign = 1.0;
                const following = try nextToken(lexer);
                if (following.tag == .newline or following.tag == .eof) {
                    parser.pending_objective_coefficient = coefficient;
                    return;
                }
                if (following.tag != .identifier) return error.InvalidSyntax;
                const column = try parser.variable(following.lexeme);
                parser.builder.columns.items[column].cost += coefficient;
            },
            .newline, .eof => return,
            .star, .caret, .left_bracket, .right_bracket, .left_paren, .right_paren => return error.UnsupportedFeature,
            else => return error.InvalidSyntax,
        }
    }
}

/// Relation and constant separated from a parsed linear row expression.
const ParsedExpression = struct {
    /// One of less-equal, equal, or greater-equal.
    relation: Tag,
    /// Right-hand constant after moving expression constants.
    constant: f64,
};

/// Parse a row expression into column-chain terms and return relation/constant.
fn parseLinearExpression(parser: *Parser, lexer: *Lexer, row: usize) types.IoError!ParsedExpression {
    var sign: f64 = 1.0;
    var constant: f64 = 0.0;
    while (true) {
        try parser.control.tick();
        const current = try nextToken(lexer);
        switch (current.tag) {
            .plus => sign = 1.0,
            .minus => sign = -1.0,
            .less_equal, .greater_equal, .equal => return .{ .relation = current.tag, .constant = constant },
            .identifier => {
                try parser.builder.addColumnTerm(row, try parser.variable(current.lexeme), sign);
                sign = 1.0;
            },
            .number => {
                const value = try parseFinite(current.lexeme);
                const following = try peekToken(lexer);
                if (following.tag == .identifier) {
                    _ = try nextToken(lexer);
                    try parser.builder.addColumnTerm(row, try parser.variable(following.lexeme), sign * value);
                } else if (following.tag == .plus or following.tag == .minus or following.tag == .less_equal or following.tag == .greater_equal or following.tag == .equal or following.tag == .newline or following.tag == .eof) {
                    constant += sign * value;
                } else {
                    return error.UnsupportedFeature;
                }
                sign = 1.0;
            },
            .star, .caret, .left_bracket, .right_bracket, .left_paren, .right_paren => return error.UnsupportedFeature,
            else => return error.InvalidSyntax,
        }
    }
}

/// Parse one labelled or anonymous linear constraint and install row bounds.
fn parseConstraint(parser: *Parser, lexer: *Lexer) types.IoError!void {
    const row_name = try consumeOptionalLabel(lexer);
    var lower_prefix: ?f64 = null;
    var probe = lexer.*;
    if (isScalarStart(try nextToken(&probe))) {
        probe = lexer.*;
        const candidate = parseScalar(&probe, false) catch null;
        if (candidate) |value| {
            const relation = try nextToken(&probe);
            if (relation.tag == .less_equal) {
                lower_prefix = value;
                lexer.* = probe;
            }
        }
    }
    const row = try parser.builder.addRow(.{ .name = row_name });
    const expression = try parseLinearExpression(parser, lexer, row);
    const rhs = try parseScalar(lexer, false);
    try expectLineEnd(lexer);
    const adjusted = rhs - expression.constant;
    if (lower_prefix) |lower| {
        if (expression.relation != .less_equal) return error.InvalidSyntax;
        const adjusted_lower = lower - expression.constant;
        if (adjusted_lower > adjusted) return error.InvalidBounds;
        parser.builder.rows.items[row].lower = adjusted_lower;
        parser.builder.rows.items[row].upper = adjusted;
    } else switch (expression.relation) {
        .less_equal => parser.builder.rows.items[row].upper = adjusted,
        .greater_equal => parser.builder.rows.items[row].lower = adjusted,
        .equal => {
            parser.builder.rows.items[row].lower = adjusted;
            parser.builder.rows.items[row].upper = adjusted;
        },
        else => unreachable,
    }
}

/// Return whether a token may begin a signed numeric bound scalar.
fn isScalarStart(value: Token) bool {
    return value.tag == .number or value.tag == .plus or value.tag == .minus or value.tag == .identifier;
}

/// Parse a signed finite scalar or permitted infinity spelling.
fn parseScalar(lexer: *Lexer, allow_infinity: bool) types.IoError!f64 {
    var sign: f64 = 1.0;
    var value = try nextToken(lexer);
    if (value.tag == .plus or value.tag == .minus) {
        if (value.tag == .minus) sign = -1.0;
        value = try nextToken(lexer);
    }
    if (value.tag == .number) return sign * try parseFinite(value.lexeme);
    if (allow_infinity and value.tag == .identifier and (std.ascii.eqlIgnoreCase(value.lexeme, "inf") or std.ascii.eqlIgnoreCase(value.lexeme, "infinity"))) return sign * std.math.inf(f64);
    return error.InvalidNumber;
}

/// Parse one LP bound form and update the referenced column interval.
fn parseBound(parser: *Parser, lexer: *Lexer) types.IoError!void {
    var probe = lexer.*;
    const first = try nextToken(&probe);
    if (first.tag == .identifier) {
        lexer.* = probe;
        const column = try parser.variable(first.lexeme);
        const operation = try nextToken(lexer);
        if (operation.tag == .identifier and std.ascii.eqlIgnoreCase(operation.lexeme, "free")) {
            parser.builder.columns.items[column].lower = -std.math.inf(f64);
            parser.builder.columns.items[column].upper = std.math.inf(f64);
            return expectLineEnd(lexer);
        }
        if (operation.tag != .equal and operation.tag != .less_equal and operation.tag != .greater_equal) return error.InvalidSyntax;
        const value = try parseScalar(lexer, true);
        try expectLineEnd(lexer);
        switch (operation.tag) {
            .equal => {
                parser.builder.columns.items[column].lower = value;
                parser.builder.columns.items[column].upper = value;
            },
            .less_equal => parser.builder.columns.items[column].upper = value,
            .greater_equal => parser.builder.columns.items[column].lower = value,
            else => unreachable,
        }
    } else {
        const lower = try parseScalar(lexer, true);
        _ = try expect(lexer, .less_equal);
        const name = try expect(lexer, .identifier);
        _ = try expect(lexer, .less_equal);
        const upper = try parseScalar(lexer, true);
        try expectLineEnd(lexer);
        if (lower > upper) return error.InvalidBounds;
        const column = try parser.variable(name.lexeme);
        parser.builder.columns.items[column].lower = lower;
        parser.builder.columns.items[column].upper = upper;
    }
}

/// Parse one variable-name list and assign the requested domain type.
fn parseTypes(parser: *Parser, lexer: *Lexer, kind: types.VariableType) types.IoError!void {
    while (true) {
        try parser.control.tick();
        const name = try nextToken(lexer);
        if (name.tag == .newline or name.tag == .eof) return;
        if (name.tag != .identifier) return error.InvalidSyntax;
        parser.builder.columns.items[try parser.variable(name.lexeme)].kind = kind;
    }
}

test "parses a linear mixed integer LP into canonical CSC" {
    const input =
        \\Maximize
        \\ obj: - 2 x + y
        \\Subject To
        \\ demand: 3 x - y >= -7
        \\Bounds
        \\ x <= 4
        \\ 0 <= y <= 1
        \\Binaries
        \\ y
        \\End
    ;
    var model = try parse(std.testing.allocator, input, "sample", .{});
    defer model.deinit();
    try std.testing.expectEqual(@as(usize, 2), model.col_cost.len);
    try std.testing.expectEqual(@as(usize, 1), model.row_lower.len);
    try std.testing.expectEqual(@as(usize, 2), model.matrix.nnz());
    try std.testing.expectEqual(types.VariableType.binary, model.col_type[1]);
    try std.testing.expectEqual(@as(f64, -2.0), model.col_cost[0]);
}

test "lexer driven parser accepts compact expressions constants and infinite bounds" {
    const input =
        \\Minimize
        \\ obj:2x-3y+5
        \\Subject To
        \\ ranged:1<=x+2y+3<=9
        \\ balance:x-y=-2
        \\Bounds
        \\ -inf<=x<=4
        \\ y free
        \\Generals
        \\ y
        \\End
    ;
    var model = try parse(std.testing.allocator, input, "compact", .{});
    defer model.deinit();

    try std.testing.expectEqualSlices(f64, &.{ 2.0, -3.0 }, model.col_cost);
    try std.testing.expectEqual(@as(f64, 5.0), model.objective_offset);
    try std.testing.expectEqualSlices(f64, &.{ -2.0, -2.0 }, model.row_lower);
    try std.testing.expectEqualSlices(f64, &.{ 6.0, -2.0 }, model.row_upper);
    try std.testing.expect(std.math.isNegativeInf(model.col_lower[0]));
    try std.testing.expectEqual(@as(f64, 4.0), model.col_upper[0]);
    try std.testing.expect(std.math.isNegativeInf(model.col_lower[1]));
    try std.testing.expect(std.math.isPositiveInf(model.col_upper[1]));
    try std.testing.expectEqual(types.VariableType.integer, model.col_type[1]);
    try std.testing.expectEqual(@as(usize, 4), model.matrix.nnz());
}

test "single lexer streams comments blank lines and compound section headers" {
    const input =
        "\\ leading comment\n" ++
        "Maximize \\ header comment\n" ++
        " obj: x\n" ++
        "      + 2y \\ expression comment\n" ++
        "\n" ++
        "Such That\n" ++
        " c0: x + y <= 3\n" ++
        "Bounds\n" ++
        " x free\n" ++
        "Semi-Continuous\n" ++
        " x\n" ++
        "Semi-Integer\n" ++
        " y\n" ++
        "End";

    var model = try parse(std.testing.allocator, input, "stream", .{});
    defer model.deinit();

    try std.testing.expectEqualSlices(f64, &.{ 1.0, 2.0 }, model.col_cost);
    try std.testing.expectEqual(@as(usize, 1), model.row_lower.len);
    try std.testing.expectEqual(@as(usize, 2), model.matrix.nnz());
    try std.testing.expectEqual(types.VariableType.semi_continuous, model.col_type[0]);
    try std.testing.expectEqual(types.VariableType.semi_integer, model.col_type[1]);
}

test "LP name index releases all hash capacity before finalization" {
    var parser = Parser.init(std.testing.allocator, "names", .{});
    defer parser.deinit();
    _ = try parser.variable("alpha");
    _ = try parser.variable("beta");
    try std.testing.expect(parser.variables.capacity() > 0);

    parser.releaseVariableIndex();
    try std.testing.expectEqual(@as(u32, 0), parser.variables.count());
    try std.testing.expectEqual(@as(u32, 0), parser.variables.capacity());
}

test "LP parser enforces file line token and semantic resource limits" {
    const source =
        \\Minimize
        \\ obj: alpha + beta
        \\Subject To
        \\ row: alpha + beta <= 1
        \\End
    ;
    try std.testing.expectError(error.FileTooLarge, parse(std.testing.allocator, source, "limits", .{ .max_file_bytes = source.len - 1 }));
    try std.testing.expectError(error.ResourceLimitExceeded, parse(std.testing.allocator, source, "limits", .{ .max_line_bytes = 7 }));
    try std.testing.expectError(error.ResourceLimitExceeded, parse(std.testing.allocator, source, "limits", .{ .max_token_bytes = 4 }));
    try std.testing.expectError(error.ResourceLimitExceeded, parse(std.testing.allocator, source, "limits", .{ .max_columns = 1 }));
    try std.testing.expectError(error.ResourceLimitExceeded, parse(std.testing.allocator, source, "limits", .{ .max_rows = 0 }));
    try std.testing.expectError(error.ResourceLimitExceeded, parse(std.testing.allocator, source, "limits", .{ .max_matrix_terms = 1 }));
    try std.testing.expectError(error.ResourceLimitExceeded, parse(std.testing.allocator, source, "limits", .{ .max_name_bytes = 4 }));
}

test "LP parser honors a pre-set cancellation flag" {
    var interrupted = std.atomic.Value(bool).init(true);
    try std.testing.expectError(error.Cancelled, parse(std.testing.allocator, "Minimize\n obj: x\nEnd\n", "cancel", .{ .interrupt_flag = &interrupted }));
}
