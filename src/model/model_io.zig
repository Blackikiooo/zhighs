//! Model I/O methods.
//!
//! Handles reading and writing model files in LP/MPS/etc. formats.

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

/// Read model data from a file into the current model
/// This can be used to read LP/MPS files.
pub fn read(self: *Model, filename: []const u8) ModelError!void {
    _ = self;
    _ = filename;
    return error.FeatureNotAvailable;
}
