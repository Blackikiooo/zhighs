//! Parameter, environment, metadata, and logging methods for `Model`.
//!
//! ## Responsibility
//!
//! Delegates model-scoped parameter operations to `Env` and exposes parameter
//! metadata, environment access, library version, error messages, and ordinary
//! logging callbacks.  Solver callbacks and callback-context parameter changes
//! belong to `model_callback.zig`.

const types = @import("types.zig");
const Model = @import("model.zig").Model;
const Env = @import("env.zig").Env;
const ParamValue = @import("env.zig").ParamValue;
const AttributeInfo = @import("attrs.zig").AttributeInfo;

const ModelError = types.ModelError;
const CallbackFunc = types.CallbackFunc;
const Version = types.Version;

pub fn version(self: Model) Version {
    _ = self;
    return Env.version();
}

pub fn tuneModel(self: *Model) ModelError!void {
    _ = self;
    return error.FeatureNotAvailable;
}

pub fn getTuneResult(self: Model, n: usize) ModelError!void {
    _ = self;
    _ = n;
    return error.FeatureNotAvailable;
}

pub fn getDblParamInfo(self: Model, name: []const u8, min_val: *f64, max_val: *f64, default_val: *f64) ModelError!void {
    _ = self;
    _ = name;
    _ = min_val;
    _ = max_val;
    _ = default_val;
    return error.FeatureNotAvailable;
}

pub fn getIntParamInfo(self: Model, name: []const u8, min_val: *i32, max_val: *i32, default_val: *i32) ModelError!void {
    _ = self;
    _ = name;
    _ = min_val;
    _ = max_val;
    _ = default_val;
    return error.FeatureNotAvailable;
}

pub fn getStrParamInfo(self: Model, name: []const u8) ModelError!struct { default: []const u8 } {
    _ = self;
    _ = name;
    return error.FeatureNotAvailable;
}

pub fn setIntParam(self: *Model, name: []const u8, value: i64) ModelError!void {
    try self.env.setIntParam(name, value);
}

pub fn getIntParam(self: Model, name: []const u8) ModelError!i64 {
    return self.env.getIntParam(name);
}

pub fn setDblParam(self: *Model, name: []const u8, value: f64) ModelError!void {
    try self.env.setDblParam(name, value);
}

pub fn getDblParam(self: Model, name: []const u8) ModelError!f64 {
    return self.env.getDblParam(name);
}

pub fn setStrParam(self: *Model, name: []const u8, value: []const u8) ModelError!void {
    try self.env.setStrParam(name, value);
}

pub fn getStrParam(self: Model, name: []const u8) ModelError![]const u8 {
    return self.env.getStrParam(name);
}

pub fn setParam(self: *Model, name: []const u8, value: ParamValue) ModelError!void {
    try self.env.setParam(name, value);
}

pub fn writeParams(self: *Model, filename: []const u8) ModelError!void {
    try self.env.writeParams(filename);
}

pub fn readParams(self: *Model, filename: []const u8) ModelError!void {
    try self.env.readParams(filename);
}

pub fn resetParams(self: *Model) ModelError!void {
    try self.env.resetParams();
}

pub fn getErrormsg(self: Model) []const u8 {
    return self.env.getErrorMessage();
}

pub fn getEnv(self: Model) *Env {
    return self.env;
}

pub fn getAttrInfo(self: Model, name: []const u8) ?AttributeInfo {
    _ = self;
    return @import("attrs.zig").lookup(name);
}

pub fn msg(self: *Model, comptime fmt: []const u8, args: anytype) void {
    self.env.log(fmt, args);
}

pub fn setLogCallbackFunc(self: *Model, callback: CallbackFunc, usrstate: ?*anyopaque) void {
    self.env.setLogCallbackFunc(callback, usrstate);
}

pub fn getLogCallbackFunc(self: Model, usrstate: *?*anyopaque) ?CallbackFunc {
    usrstate.* = self.env.log_usrstate;
    return self.env.log_callback;
}
