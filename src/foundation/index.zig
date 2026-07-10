const std = @import("std");
const Int = @import("int.zig");

const HInt = Int.HInt;
const HUInt = Int.HUInt;

const needs_huint_to_usize_check = @bitSizeOf(HUInt) > @bitSizeOf(usize);
const needs_usize_to_huint_check = @bitSizeOf(usize) > @bitSizeOf(HUInt);
const hint_fits_usize = @bitSizeOf(HInt) <= @bitSizeOf(usize);

pub const IndexError = error{
    NegativeIndex,
    IndexOverflow,
    ReservedIndex,
};

const none_raw = std.math.maxInt(HUInt);

const IndexKind = enum(u1) {
    row,
    column,
};

/// Creates a strongly typed, non-negative model identifier.
///
/// Each `kind` instantiation produces a distinct enum type, so row and column
/// identifiers cannot be mixed accidentally even though both use `HUInt` as
/// their runtime representation.
fn IndexId(comptime kind_value: IndexKind) type {
    return enum(HUInt) {
        _,

        const Self = @This();

        pub const is_row = kind_value == .row;

        /// Creates an identifier from the native unsigned HiGHS type.
        ///
        /// The cast also establishes the invariant required by `toUsize`.
        /// On targets where every HUInt fits in usize, the check is removed at
        /// compile time.
        pub inline fn init(value: HUInt) IndexError!Self {
            if (value == none_raw) return IndexError.ReservedIndex;

            if (comptime needs_huint_to_usize_check) {
                _ = std.math.cast(usize, value) orelse
                    return IndexError.IndexOverflow;
            }

            return @enumFromInt(value);
        }

        /// Creates an identifier without an error union for trusted hot paths.
        ///
        /// The caller guarantees that `value` is not reserved and fits in
        /// `usize`. Debug builds assert these preconditions; ReleaseFast emits
        /// only the enum conversion.
        pub inline fn initAssumeValid(value: HUInt) Self {
            std.debug.assert(value != none_raw);

            if (comptime needs_huint_to_usize_check) {
                std.debug.assert(std.math.cast(usize, value) != null);
            }

            return @enumFromInt(value);
        }

        /// Converts an external, signed HiGHS index at an API boundary.
        pub inline fn fromHInt(value: HInt) IndexError!Self {
            if (value < 0) return IndexError.NegativeIndex;

            return init(@intCast(value));
        }

        /// Converts a previously validated signed index without an error union.
        pub inline fn fromHIntAssumeValid(value: HInt) Self {
            std.debug.assert(value >= 0);
            return initAssumeValid(@intCast(value));
        }

        /// Converts a Zig memory index to the configured HiGHS index width.
        pub inline fn fromUsize(value: usize) IndexError!Self {
            const raw_value: HUInt =
                if (comptime needs_usize_to_huint_check)
                    std.math.cast(HUInt, value) orelse
                        return IndexError.IndexOverflow
                else
                    @intCast(value);

            return init(raw_value);
        }

        /// Converts a previously validated slice index without an error union.
        pub inline fn fromUsizeAssumeValid(value: usize) Self {
            if (comptime needs_usize_to_huint_check) {
                std.debug.assert(std.math.cast(HUInt, value) != null);
            }

            return initAssumeValid(@intCast(value));
        }

        pub inline fn raw(self: Self) HUInt {
            return @intFromEnum(self);
        }

        /// Returns the slice index established as valid by the constructors.
        pub inline fn toUsize(self: Self) usize {
            return @intCast(@intFromEnum(self));
        }
    };
}

pub const RowId = IndexId(.row);
pub const ColId = IndexId(.column);

/// Creates a compact optional representation for a strongly typed ID.
///
/// `none` uses the reserved maximum HUInt value, so this type has exactly the
/// same size and alignment as its non-optional ID. Calling `get` produces a
/// temporary Zig optional for ergonomic branching without changing the stored
/// representation.
fn OptionalId(comptime Id: type) type {
    return enum(HUInt) {
        none = none_raw,
        _,

        const Self = @This();

        pub const Value = Id;

        pub inline fn some(id: Id) Self {
            std.debug.assert(id.raw() != none_raw);
            return @enumFromInt(id.raw());
        }

        pub inline fn fromNullable(id: ?Id) Self {
            return if (id) |value| some(value) else .none;
        }

        pub inline fn get(self: Self) ?Id {
            const raw_value = @intFromEnum(self);
            if (raw_value == none_raw) return null;

            return @enumFromInt(raw_value);
        }

        pub inline fn isSome(self: Self) bool {
            return self != .none;
        }

        pub inline fn isNone(self: Self) bool {
            return self == .none;
        }

        pub inline fn raw(self: Self) HUInt {
            return @intFromEnum(self);
        }
    };
}

pub const OptionalRowId = OptionalId(RowId);
pub const OptionalColId = OptionalId(ColId);

test "RowId and ColId are distinct types" {
    try std.testing.expect(RowId != ColId);
    try std.testing.expect(RowId.is_row);
    try std.testing.expect(!ColId.is_row);
}

test "compact optional identifier types remain distinct" {
    try std.testing.expect(OptionalRowId != OptionalColId);
    try std.testing.expect(OptionalRowId.Value == RowId);
    try std.testing.expect(OptionalColId.Value == ColId);
}

test "unsigned identifiers round trip" {
    const row = try RowId.init(12);
    const col = try ColId.init(23);

    try std.testing.expectEqual(@as(HUInt, 12), row.raw());
    try std.testing.expectEqual(@as(usize, 12), row.toUsize());
    try std.testing.expectEqual(@as(HUInt, 23), col.raw());
    try std.testing.expectEqual(@as(usize, 23), col.toUsize());
}

test "signed boundary conversion rejects negative identifiers" {
    try std.testing.expectError(IndexError.NegativeIndex, RowId.fromHInt(-1));
    try std.testing.expectError(IndexError.NegativeIndex, ColId.fromHInt(-1));
}

test "signed and usize conversions round trip" {
    const row = try RowId.fromHInt(7);
    const col = try ColId.fromUsize(9);

    try std.testing.expectEqual(@as(HUInt, 7), row.raw());
    try std.testing.expectEqual(@as(usize, 9), col.toUsize());
}

test "trusted constructors return bare strongly typed identifiers" {
    const direct: RowId = RowId.initAssumeValid(11);
    const from_hint: RowId = RowId.fromHIntAssumeValid(12);
    const from_usize: ColId = ColId.fromUsizeAssumeValid(13);

    try std.testing.expectEqual(@as(HUInt, 11), direct.raw());
    try std.testing.expectEqual(@as(HUInt, 12), from_hint.raw());
    try std.testing.expectEqual(@as(HUInt, 13), from_usize.raw());
    try std.testing.expect(@TypeOf(RowId.initAssumeValid(0)) == RowId);
    try std.testing.expect(@TypeOf(ColId.fromUsizeAssumeValid(0)) == ColId);
}

test "maximum HUInt value is reserved for compact none" {
    try std.testing.expectError(
        IndexError.ReservedIndex,
        RowId.init(none_raw),
    );
    try std.testing.expectError(
        IndexError.ReservedIndex,
        ColId.init(none_raw),
    );
}

test "compact optional identifiers store some and none" {
    const row = try RowId.init(17);
    const some_row = OptionalRowId.some(row);
    const no_row: OptionalRowId = .none;

    try std.testing.expect(some_row.isSome());
    try std.testing.expect(!some_row.isNone());
    try std.testing.expectEqual(row, some_row.get().?);

    try std.testing.expect(no_row.isNone());
    try std.testing.expect(!no_row.isSome());
    try std.testing.expect(no_row.get() == null);
    try std.testing.expectEqual(none_raw, no_row.raw());
}

test "compact optional identifiers convert from Zig optional" {
    const col = try ColId.init(19);

    try std.testing.expectEqual(
        col,
        OptionalColId.fromNullable(col).get().?,
    );
    try std.testing.expect(
        OptionalColId.fromNullable(null).isNone(),
    );
}

test "compact optional identifiers have no storage overhead" {
    try std.testing.expectEqual(@sizeOf(RowId), @sizeOf(OptionalRowId));
    try std.testing.expectEqual(@alignOf(RowId), @alignOf(OptionalRowId));
    try std.testing.expectEqual(@sizeOf(ColId), @sizeOf(OptionalColId));
    try std.testing.expectEqual(@alignOf(ColId), @alignOf(OptionalColId));
}

test "usize conversion checks configured integer width" {
    if (comptime needs_usize_to_huint_check) {
        const too_large = @as(usize, std.math.maxInt(HUInt)) + 1;
        try std.testing.expectError(
            IndexError.IndexOverflow,
            RowId.fromUsize(too_large),
        );
    }
}

test "zero round trips through all constructors" {
    const direct = try RowId.init(0);
    try std.testing.expectEqual(@as(HUInt, 0), direct.raw());
    try std.testing.expectEqual(@as(usize, 0), direct.toUsize());

    const from_hint = try RowId.fromHInt(0);
    try std.testing.expectEqual(@as(HUInt, 0), from_hint.raw());

    const from_usize = try RowId.fromUsize(@as(usize, 0));
    try std.testing.expectEqual(@as(HUInt, 0), from_usize.raw());
    try std.testing.expectEqual(@as(usize, 0), from_usize.toUsize());
}

test "maximum non-reserved value is accepted by init" {
    if (comptime !needs_huint_to_usize_check) {
        const max_valid = none_raw - 1;
        const row = try RowId.init(max_valid);
        try std.testing.expectEqual(max_valid, row.raw());
        try std.testing.expectEqual(@as(usize, max_valid), row.toUsize());
    }
}

test "fromHInt maximum positive value" {
    if (comptime hint_fits_usize) {
        const max_hint = std.math.maxInt(HInt);
        const row = try RowId.fromHInt(max_hint);
        try std.testing.expectEqual(@as(HUInt, @intCast(max_hint)), row.raw());
    }
}

test "fromUsize rejects reserved maximum value" {
    const max_huint = std.math.maxInt(HUInt);
    const bounded = std.math.cast(usize, max_huint) orelse return;
    try std.testing.expectError(
        IndexError.ReservedIndex,
        RowId.fromUsize(bounded),
    );
}

test "compact optional round trips at maximum valid value" {
    const max_valid = try RowId.init(none_raw - 1);
    const wrapped = OptionalRowId.some(max_valid);

    try std.testing.expect(wrapped.isSome());
    try std.testing.expect(!wrapped.isNone());
    try std.testing.expectEqual(max_valid, wrapped.get().?);
    try std.testing.expectEqual(none_raw - 1, wrapped.raw());
}
