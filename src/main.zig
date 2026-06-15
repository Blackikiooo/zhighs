const std = @import("std");
const Io = std.Io;

const zhighs = @import("zhighs");

pub fn main(init: std.process.Init) !void {
    _ = init; // Prints to stderr, unbuffered, ignoring potential errors.
}
