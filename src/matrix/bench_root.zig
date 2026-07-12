//! Matrix-only package surface for benchmarks and profiling tools.
//!
//! Keeping the benchmark executable independent of model/API modules lets
//! matrix performance work continue while higher layers are being edited, and
//! prevents unrelated code from affecting benchmark compilation or code size.

const foundation = @import("foundation");

pub const matrix = @import("root.zig");
pub const HUInt = foundation.HUInt;
pub const HCD = foundation.HCD;
pub const RowId = foundation.RowId;
pub const ColId = foundation.ColId;
