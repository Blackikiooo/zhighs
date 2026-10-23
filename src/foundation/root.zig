const std = @import("std");
const Int = @import("int.zig");
const Double = @import("double.zig");
const Index = @import("index.zig");
const tolerance = @import("tolerance.zig");

pub const HInt = Int.HInt;
pub const HUInt = Int.HUInt;
pub const HD = Double.HD;
pub const HCD = Double.HCD;
pub const IndexError = Index.IndexError;
pub const RowId = Index.RowId;
pub const ColId = Index.ColId;
pub const OptionalRowId = Index.OptionalRowId;
pub const OptionalColId = Index.OptionalColId;
/// Model-facing coefficient floor; see `tolerance.zig`.
pub const model_coefficient_tolerance = tolerance.model_coefficient_tolerance;
/// Whether a model coefficient survives the shared dropping threshold.
pub const retainsModelCoefficient = tolerance.retainsModelCoefficient;

test {
    std.testing.refAllDecls(@This());
}
