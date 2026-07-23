//! Tokens produced by the zero-copy LP lexer.

/// Lexical category emitted by the zero-copy LP lexer.
pub const Tag = enum(u8) {
    // Payload tokens: `lexeme` contains the source spelling.
    identifier,
    number,

    // Linear-expression punctuation and relations.
    plus,
    minus,
    star,
    caret,
    colon,
    comma,
    left_paren,
    right_paren,
    left_bracket,
    right_bracket,
    less_equal,
    greater_equal,
    equal,

    // Stream boundaries. Horizontal whitespace and comments are not emitted.
    newline,
    eof,
};

/// Position of a token's first byte in the original source buffer.
pub const Location = struct {
    /// Zero-based byte offset in the complete source.
    byte_offset: usize,
    /// One-based source line.
    line: usize,
    /// One-based byte column. LP model syntax is ASCII-oriented.
    column: usize,
};

/// One lexer item; both spelling and location refer to the original source.
pub const Token = struct {
    /// Category used by the parser's grammar dispatch.
    tag: Tag,
    /// Borrowed directly from the lexer's input buffer.
    lexeme: []const u8,
    /// Location of the first byte in `lexeme` (or the end cursor for EOF).
    location: Location,
};
