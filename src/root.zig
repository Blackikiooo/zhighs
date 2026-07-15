pub const std = @import("std");
pub const foundation = @import("foundation");
/// Compatibility alias for callers using the pre-refactor module name.
pub const types = foundation;
pub const matrix = @import("matrix");
pub const model = @import("model");
pub const nla = @import("nla/root.zig");
pub const lp = @import("lp");
pub const presolve = @import("presolve/root.zig");
pub const analysis = @import("analysis/root.zig");
pub const qp = @import("qp/root.zig");
pub const ipm = @import("ipm/root.zig");
pub const pdlp = @import("pdlp/root.zig");
pub const mip = @import("mip/root.zig");
pub const framework = @import("framework/root.zig");
pub const plugin = @import("plugin/root.zig");
pub const plugins_builtin = @import("plugins_builtin/root.zig");
pub const solver = @import("solver");
pub const api = @import("api/root.zig");
pub const io = @import("io");
pub const diagnostics = @import("diagnostics/root.zig");
pub const parallel = @import("parallel/root.zig");
pub const bindings = @import("bindings/root.zig");
pub const HInt = foundation.HInt;
pub const HUInt = foundation.HUInt;
pub const HD = foundation.HD;
pub const HCD = foundation.HCD;
pub const IndexError = foundation.IndexError;
pub const RowId = foundation.RowId;
pub const ColId = foundation.ColId;
pub const OptionalRowId = foundation.OptionalRowId;
pub const OptionalColId = foundation.OptionalColId;
test {
    std.testing.refAllDecls(@This());
}
