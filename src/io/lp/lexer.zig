//! Allocation-free, zero-copy lexer for LP model files.
//!
//! `Lexer` is intentionally a small value type. Copying it creates a cheap
//! checkpoint, which lets the parser probe ambiguous prefixes without token
//! buffering or heap allocation.

const std = @import("std");
const token = @import("token.zig");

pub const Token = token.Token;
pub const Tag = token.Tag;
pub const Location = token.Location;

/// Lexical errors are deliberately narrow. Grammar and semantic failures are
/// owned by the parser and mapped to the richer `io.IoError` error set there.
pub const Error = error{ InvalidCharacter, ResourceLimitExceeded };

pub const Lexer = struct {
    /// Source bytes are borrowed. Every token lexeme points into this slice,
    /// so the caller must keep it alive while tokens are being consumed.
    input: []const u8,
    /// Cursor relative to `input`; `base_offset + index` is the absolute byte
    /// position in the complete model file.
    index: usize = 0,
    /// Source locations are one-based for diagnostics and editor integration.
    line: usize = 1,
    column: usize = 1,
    /// Absolute byte offset of `input[0]`. It is zero for a whole-file lexer.
    base_offset: usize = 0,
    line_start_index: usize = 0,
    max_line_bytes: usize = std.math.maxInt(usize),
    max_token_bytes: usize = std.math.maxInt(usize),

    /// Start lexing a complete in-memory LP source.
    pub fn init(input: []const u8) Lexer {
        return .{ .input = input };
    }

    pub fn initWithLimits(input: []const u8, max_line_bytes: usize, max_token_bytes: usize) Lexer {
        return .{ .input = input, .max_line_bytes = max_line_bytes, .max_token_bytes = max_token_bytes };
    }

    /// Construct a lexer over a source window while retaining locations in the
    /// complete file. This is reserved for future chunked/windowed inputs; the
    /// current LP parser uses one whole-file stream over a buffer or mmap.
    pub fn initWindow(input: []const u8, byte_offset: usize, line: usize, column: usize) Lexer {
        return .{ .input = input, .base_offset = byte_offset, .line = line, .column = column };
    }

    /// Return the next significant token.
    ///
    /// Horizontal whitespace and comments are skipped. Newlines remain tokens
    /// because the parser consumes the complete file as one stream and uses
    /// them as statement and section-header boundaries. No allocation occurs.
    pub fn next(self: *Lexer) Error!Token {
        while (self.index < self.input.len) {
            const char = self.input[self.index];
            switch (char) {
                ' ', '\t', '\r' => self.advanceByte(),
                '\\' => self.skipComment(),
                '\n' => return self.newlineToken(),
                '+' => return self.single(.plus),
                '-' => return self.single(.minus),
                '*' => return self.single(.star),
                '^' => return self.single(.caret),
                ':' => return self.single(.colon),
                ',' => return self.single(.comma),
                '(' => return self.single(.left_paren),
                ')' => return self.single(.right_paren),
                '[' => return self.single(.left_bracket),
                ']' => return self.single(.right_bracket),
                '=' => return self.single(.equal),
                '<' => return self.relation(.less_equal, '='),
                '>' => return self.relation(.greater_equal, '='),
                '.', '0'...'9' => return self.number(),
                else => return self.identifier(),
            }
        }
        if (self.index - self.line_start_index > self.max_line_bytes) return error.ResourceLimitExceeded;
        return .{ .tag = .eof, .lexeme = self.input[self.input.len..], .location = self.location() };
    }

    /// Inspect the next token without advancing this lexer. Since `Lexer`
    /// contains only a borrowed slice and scalar cursor state, a checkpoint is
    /// just a small value copy and requires neither allocation nor rollback.
    pub fn peek(self: *Lexer) Error!Token {
        var checkpoint = self.*;
        return checkpoint.next();
    }

    fn location(self: Lexer) Location {
        return .{ .byte_offset = self.base_offset + self.index, .line = self.line, .column = self.column };
    }

    fn advanceByte(self: *Lexer) void {
        self.index += 1;
        self.column += 1;
    }

    fn skipComment(self: *Lexer) void {
        // In LP files a backslash starts a comment extending to end-of-line.
        // Do not consume '\n': it still carries grammar and location meaning.
        while (self.index < self.input.len and self.input[self.index] != '\n') self.advanceByte();
    }

    fn single(self: *Lexer, tag: Tag) Token {
        const start = self.index;
        const source_location = self.location();
        self.advanceByte();
        return .{ .tag = tag, .lexeme = self.input[start..self.index], .location = source_location };
    }

    fn relation(self: *Lexer, tag: Tag, expected: u8) Error!Token {
        const start = self.index;
        const source_location = self.location();
        self.advanceByte();
        // The supported LP grammar uses <= and >=, never bare < or >. Failing
        // here gives the parser a precise lexical error instead of two tokens.
        if (self.index >= self.input.len or self.input[self.index] != expected) return error.InvalidCharacter;
        self.advanceByte();
        return .{ .tag = tag, .lexeme = self.input[start..self.index], .location = source_location };
    }

    fn newlineToken(self: *Lexer) Error!Token {
        if (self.index - self.line_start_index > self.max_line_bytes) return error.ResourceLimitExceeded;
        const start = self.index;
        const source_location = self.location();
        self.index += 1;
        self.line += 1;
        self.column = 1;
        self.line_start_index = self.index;
        return .{ .tag = .newline, .lexeme = self.input[start..self.index], .location = source_location };
    }

    fn number(self: *Lexer) Error!Token {
        // Sign characters are separate tokens so the parser can distinguish a
        // unary sign from the operator between two linear terms.
        const start = self.index;
        const source_location = self.location();
        var digits_before: usize = 0;
        while (self.index < self.input.len and std.ascii.isDigit(self.input[self.index])) {
            self.advanceByte();
            digits_before += 1;
        }
        var digits_after: usize = 0;
        if (self.index < self.input.len and self.input[self.index] == '.') {
            self.advanceByte();
            while (self.index < self.input.len and std.ascii.isDigit(self.input[self.index])) {
                self.advanceByte();
                digits_after += 1;
            }
        }
        if (digits_before == 0 and digits_after == 0) return error.InvalidCharacter;
        if (self.index < self.input.len and (self.input[self.index] == 'e' or self.input[self.index] == 'E')) {
            // Once e/E is consumed it must form a complete exponent. Accepting
            // a partial exponent would silently split malformed numeric input.
            self.advanceByte();
            if (self.index < self.input.len and (self.input[self.index] == '+' or self.input[self.index] == '-')) self.advanceByte();
            const exponent_start = self.index;
            while (self.index < self.input.len and std.ascii.isDigit(self.input[self.index])) self.advanceByte();
            if (self.index == exponent_start) return error.InvalidCharacter;
        }
        if (self.index - start > self.max_token_bytes) return error.ResourceLimitExceeded;
        return .{ .tag = .number, .lexeme = self.input[start..self.index], .location = source_location };
    }

    fn identifier(self: *Lexer) Error!Token {
        const start = self.index;
        const source_location = self.location();
        while (self.index < self.input.len and !isDelimiter(self.input[self.index])) self.advanceByte();
        if (self.index == start) return error.InvalidCharacter;
        if (self.index - start > self.max_token_bytes) return error.ResourceLimitExceeded;
        return .{ .tag = .identifier, .lexeme = self.input[start..self.index], .location = source_location };
    }

    fn isDelimiter(char: u8) bool {
        // Keywords are intentionally returned as identifiers. Section and
        // context-specific keyword recognition belongs to the parser, keeping
        // this scanner compact and making names zero-copy.
        return std.ascii.isWhitespace(char) or switch (char) {
            '+', '-', '*', '^', ':', ',', '(', ')', '[', ']', '=', '<', '>', '\\' => true,
            else => false,
        };
    }
};

test "lexer preserves exact lexemes in compact linear syntax" {
    const source = "c0:x+1.25e-2y<=3\\ comment\nBounds";
    var lexer = Lexer.init(source);
    const Expected = struct { tag: Tag, lexeme: []const u8 };
    const expected = [_]Expected{
        .{ .tag = .identifier, .lexeme = "c0" },
        .{ .tag = .colon, .lexeme = ":" },
        .{ .tag = .identifier, .lexeme = "x" },
        .{ .tag = .plus, .lexeme = "+" },
        .{ .tag = .number, .lexeme = "1.25e-2" },
        .{ .tag = .identifier, .lexeme = "y" },
        .{ .tag = .less_equal, .lexeme = "<=" },
        .{ .tag = .number, .lexeme = "3" },
        .{ .tag = .newline, .lexeme = "\n" },
        .{ .tag = .identifier, .lexeme = "Bounds" },
        .{ .tag = .eof, .lexeme = "" },
    };
    for (expected) |wanted| {
        const actual = try lexer.next();
        try std.testing.expectEqual(wanted.tag, actual.tag);
        try std.testing.expectEqualStrings(wanted.lexeme, actual.lexeme);
    }
}

test "lexer recognizes punctuation reserved for future nonlinear grammar" {
    var lexer = Lexer.init("x*y^2,[z]");
    const expected = [_]Tag{
        .identifier,    .star,  .identifier,   .caret,
        .number,        .comma, .left_bracket, .identifier,
        .right_bracket, .eof,
    };
    for (expected) |tag| try std.testing.expectEqual(tag, (try lexer.next()).tag);
}

test "lexer reports stable byte line and column locations" {
    var lexer = Lexer.init("x\n  y");
    const x = try lexer.next();
    const newline = try lexer.next();
    const y = try lexer.next();
    try std.testing.expectEqual(@as(usize, 0), x.location.byte_offset);
    try std.testing.expectEqual(@as(usize, 1), newline.location.line);
    try std.testing.expectEqual(@as(usize, 2), y.location.line);
    try std.testing.expectEqual(@as(usize, 3), y.location.column);
    try std.testing.expectEqualStrings("y", y.lexeme);
}

test "copying lexer creates an independent cursor checkpoint" {
    var lexer = Lexer.init("1 <= x");
    var probe = lexer;
    try std.testing.expectEqual(Tag.number, (try probe.next()).tag);
    try std.testing.expectEqual(Tag.less_equal, (try probe.next()).tag);
    try std.testing.expectEqual(Tag.number, (try lexer.next()).tag);
}

test "window lexer preserves absolute source location" {
    var lexer = Lexer.initWindow("  x", 128, 17, 9);
    const name = try lexer.next();
    try std.testing.expectEqual(Tag.identifier, name.tag);
    try std.testing.expectEqual(@as(usize, 130), name.location.byte_offset);
    try std.testing.expectEqual(@as(usize, 17), name.location.line);
    try std.testing.expectEqual(@as(usize, 11), name.location.column);
}

test "lexer rejects incomplete relations and malformed exponents" {
    var relation = Lexer.init("x < 1");
    _ = try relation.next();
    try std.testing.expectError(error.InvalidCharacter, relation.next());

    var exponent = Lexer.init("1e+");
    try std.testing.expectError(error.InvalidCharacter, exponent.next());
}
