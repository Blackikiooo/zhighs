//! Model I/O methods.
//!
//! Handles reading and writing model files in LP/MPS/etc. formats.
//!
//! ## Responsibility
//!
//! Defines the model-file boundary, including format dispatch and the update
//! semantics applied before reading or writing.  `read` is only a compatibility
//! alias for the canonical `readModel` implementation.  Parameter files are
//! handled by `model_params.zig`.

const types = @import("types.zig");
const Model = @import("model.zig").Model;

const ModelError = types.ModelError;

/// Write the model to a file
/// The format is inferred from the file extension (.mps, .lp, .mps.gz, ...).
pub fn writeModel(self: *Model, filename: []const u8) ModelError!void {
    // Flush first so the written model is current.
    try self.updateModel();
    _ = filename;
    // Future: delegate to io module for LP/MPS writing.
    return error.FeatureNotAvailable;
}

/// Read model data from a file into the current model.
///
/// This is the canonical implementation for both the `readModel` and `read`
/// public entry points.  Keeping the compatibility name as an alias prevents
/// the two APIs from acquiring different update/error semantics later.
pub fn readModel(self: *Model, filename: []const u8) ModelError!void {
    // Reading will eventually replace/merge committed state.  Flush pending
    // changes first, matching the rest of the model I/O boundary.
    try self.updateModel();
    _ = filename;
    return error.FeatureNotAvailable;
}

/// Short compatibility spelling for `readModel`.
pub const read = readModel;
