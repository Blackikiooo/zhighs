//! Print the concrete Zig ABI layout used by matrix control blocks and views.
//! This complements per-allocation bytes/nnz accounting: slice headers live in
//! these structs, while matrix elements live in separately allocated SoA data.

const std = @import("std");
const zhighs = @import("zhighs");

const BuilderTripletLayout = struct {
    row: zhighs.RowId,
    col: zhighs.ColId,
    value: f64,
    sequence: usize,
};

fn printStruct(comptime T: type) void {
    std.debug.print("type\t{s}\tsize\t{d}\talign\t{d}\n", .{ @typeName(T), @sizeOf(T), @alignOf(T) });
    inline for (@typeInfo(T).@"struct".fields) |field| {
        std.debug.print("field\t{s}\t{s}\toffset\t{d}\tsize\t{d}\talign\t{d}\n", .{
            field.name,
            @typeName(field.type),
            @offsetOf(T, field.name),
            @sizeOf(field.type),
            @alignOf(field.type),
        });
    }
}

pub fn main() void {
    std.debug.print("scalar\tusize\tsize\t{d}\talign\t{d}\n", .{ @sizeOf(usize), @alignOf(usize) });
    std.debug.print("scalar\tHUInt\tsize\t{d}\talign\t{d}\n", .{ @sizeOf(zhighs.HUInt), @alignOf(zhighs.HUInt) });
    std.debug.print("scalar\tRowId\tsize\t{d}\talign\t{d}\n", .{ @sizeOf(zhighs.RowId), @alignOf(zhighs.RowId) });
    std.debug.print("scalar\tColId\tsize\t{d}\talign\t{d}\n", .{ @sizeOf(zhighs.ColId), @alignOf(zhighs.ColId) });
    std.debug.print("scalar\t[]f64\tsize\t{d}\talign\t{d}\n", .{ @sizeOf([]f64), @alignOf([]f64) });
    std.debug.print("scalar\t?[]HUInt\tsize\t{d}\talign\t{d}\n", .{ @sizeOf(?[]zhighs.HUInt), @alignOf(?[]zhighs.HUInt) });

    printStruct(zhighs.matrix.CscMatrix);
    printStruct(zhighs.matrix.CscView);
    printStruct(zhighs.matrix.CsrView);
    printStruct(zhighs.matrix.CsrBuffers);
    printStruct(zhighs.matrix.TransposeBuffers);
    printStruct(zhighs.matrix.CscBuildBuffers);
    printStruct(zhighs.matrix.CscTransformBuffers);
    printStruct(zhighs.matrix.MatrixBuilder);
    printStruct(zhighs.matrix.SparseVectorView(zhighs.RowId));
    printStruct(zhighs.matrix.SparseVector(zhighs.RowId));
    printStruct(BuilderTripletLayout);
}
