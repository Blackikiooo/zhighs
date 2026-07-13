//! Fine-grained model revision counters used for solver-cache invalidation.
//!
//! Each counter changes independently so persistent solver sessions can retain
//! basis factorizations across objective, bound, and RHS edits while reliably
//! invalidating them after structural or matrix-value changes.

const std = @import("std");

pub const RevisionKind = enum { structure, matrix_values, bounds, objective };

pub const RevisionSet = struct {
    structure: u64 = 0,
    matrix_values: u64 = 0,
    bounds: u64 = 0,
    objective: u64 = 0,

    pub fn bump(self: *RevisionSet, kind: RevisionKind) error{RevisionOverflow}!void {
        const target = switch (kind) {
            .structure => &self.structure,
            .matrix_values => &self.matrix_values,
            .bounds => &self.bounds,
            .objective => &self.objective,
        };
        target.* = std.math.add(u64, target.*, 1) catch return error.RevisionOverflow;
    }
};

test "fine-grained revisions change independently" {
    var revisions = RevisionSet{};
    try revisions.bump(.bounds);
    try std.testing.expectEqual(@as(u64, 1), revisions.bounds);
    try std.testing.expectEqual(@as(u64, 0), revisions.structure);
}
