//! Solver callback registration and callback-context methods for `Model`.
//!
//! ## Responsibility
//!
//! Owns callback registration, interruption, callback-time queries/actions,
//! and callback-settable parameter entry points.  Environment logging callbacks
//! remain in `model_params.zig`; normal model edits remain in their domain
//! modules.

const types = @import("types.zig");
const Model = @import("model.zig").Model;
const ParamValue = @import("env.zig").ParamValue;

const ModelError = types.ModelError;
const CallbackFunc = types.CallbackFunc;

pub fn setCallbackFunc(self: *Model, callback: CallbackFunc, usrstate: ?*anyopaque) void {
    self.env.setCallbackFunc(callback, usrstate);
}

pub fn getCallbackFunc(self: Model, usrstate: *?*anyopaque) ?CallbackFunc {
    usrstate.* = self.env.usrstate;
    return self.env.callback;
}

pub fn terminate(self: *Model) void {
    self.interrupted.store(true, .release);
}

pub fn cbGet(self: Model, what: i32, result: *anyopaque) ModelError!void {
    _ = self;
    _ = what;
    _ = result;
    return error.FeatureNotAvailable;
}

pub fn cbCut(self: *Model, cutlen: usize, cutind: []const i32, cutval: []const f64) ModelError!void {
    _ = self;
    _ = cutlen;
    _ = cutind;
    _ = cutval;
    return error.FeatureNotAvailable;
}

pub fn cbLazy(self: *Model, lazylen: usize, lazyind: []const i32, lazyval: []const f64) ModelError!void {
    _ = self;
    _ = lazylen;
    _ = lazyind;
    _ = lazyval;
    return error.FeatureNotAvailable;
}

pub fn cbSolution(self: *Model, solution: []const f64, obj: *f64) ModelError!void {
    _ = self;
    _ = solution;
    _ = obj;
    return error.FeatureNotAvailable;
}

pub fn cbSetDblParam(self: *Model, name: []const u8, value: f64) ModelError!void {
    _ = self;
    _ = name;
    _ = value;
    return error.FeatureNotAvailable;
}

pub fn cbSetIntParam(self: *Model, name: []const u8, value: i64) ModelError!void {
    _ = self;
    _ = name;
    _ = value;
    return error.FeatureNotAvailable;
}

pub fn cbSetStrParam(self: *Model, name: []const u8, value: []const u8) ModelError!void {
    _ = self;
    _ = name;
    _ = value;
    return error.FeatureNotAvailable;
}

pub fn cbSetParam(self: *Model, name: []const u8, value: ParamValue) ModelError!void {
    _ = self;
    _ = name;
    _ = value;
    return error.FeatureNotAvailable;
}

pub fn cbProceed(self: *Model) ModelError!void {
    _ = self;
    return error.FeatureNotAvailable;
}

pub fn cbStopOneMultiObj(self: *Model, multiobjnum: i32) ModelError!void {
    _ = self;
    _ = multiobjnum;
    return error.FeatureNotAvailable;
}

pub fn setCallbackFuncAdv(self: *Model, callback: CallbackFunc, usrstate: ?*anyopaque, wheres: u32) void {
    _ = wheres;
    self.env.setCallbackFunc(callback, usrstate);
}
