//! Matrix benchmark package adapter.
//!
//! This file intentionally lives with benchmark sources. It exposes the small
//! compatibility surface used by the existing benchmark programs without
//! pulling model or API modules into performance binaries.

const foundation = @import("foundation");

pub const matrix = @import("matrix");
pub const HUInt = foundation.HUInt;
pub const HCD = foundation.HCD;
pub const RowId = foundation.RowId;
pub const ColId = foundation.ColId;
