//! Model environment.
//!
//! The environment owns configuration shared across all models that live
//! within it — most notably parameter defaults and the log/output channel.
//! Callers create a single `Env` (via `Env.init`) and then create one or more
//! `Model` handles against it.

const std = @import("std");
const types = @import("types.zig");

const ModelError = types.ModelError;
const INFINITY = types.INFINITY;
const Version = types.Version;
const CallbackFunc = types.CallbackFunc;

// ── Parameter storage ───────────────────────────────────────────────────

/// Value stored per parameter; matched to the parameter's declared type.
pub const ParamValue = union(enum) {
    int: i64,
    double: f64,
    string: []const u8,
};

// ── Environment ─────────────────────────────────────────────────────────

pub const Env = struct {
    allocator: std.mem.Allocator,

    /// Named parameter table.  Parameters are resolved by string name.
    params: std.StringHashMap(ParamValue),

    /// If non-null, solver output is also written to this file.
    log_filename: ?[]const u8 = null,

    /// Last error message (cleared on each new operation).
    last_error: ?[]const u8 = null,

    /// User-registered callback for solver progress.
    callback: ?CallbackFunc = null,
    usrstate: ?*anyopaque = null,

    /// User-registered log callback.
    log_callback: ?CallbackFunc = null,
    log_usrstate: ?*anyopaque = null,

    const Self = @This();

    /// Create an environment with the default parameter set.
    ///
    /// `log_filename` may be omitted (or empty) to suppress file logging.
    /// Pass a file path to write a solver log alongside stderr output.
    pub fn init(allocator: std.mem.Allocator, log_filename: ?[]const u8) ModelError!Self {
        var self = Self{
            .allocator = allocator,
            .params = std.StringHashMap(ParamValue).init(allocator),
            .log_filename = if (log_filename) |n| if (n.len > 0) try allocator.dupe(u8, n) else null else null,
        };
        try self.setDefaults();
        return self;
    }

    /// Create an environment without file logging.
    pub fn initSimple(allocator: std.mem.Allocator) ModelError!Self {
        var self = Self{
            .allocator = allocator,
            .params = std.StringHashMap(ParamValue).init(allocator),
        };
        try self.setDefaults();
        return self;
    }

    /// Release all resources held by the environment.
    pub fn deinit(self: *Self) void {
        if (self.log_filename) |n| self.allocator.free(n);
        self.params.deinit();
        self.* = undefined;
    }

    // ── Parameter get/set ───────────────────────────────────────────────

    /// Set an integer-valued parameter.
    pub fn setIntParam(self: *Self, name: []const u8, value: i64) ModelError!void {
        try validateIntParam(name, value);
        self.params.put(name, .{ .int = value }) catch return error.OutOfMemory;
    }

    /// Set a double-valued parameter.
    pub fn setDblParam(self: *Self, name: []const u8, value: f64) ModelError!void {
        try validateDblParam(name, value);
        self.params.put(name, .{ .double = value }) catch return error.OutOfMemory;
    }

    /// Set a string-valued parameter.
    pub fn setStrParam(self: *Self, name: []const u8, value: []const u8) ModelError!void {
        try validateStrParam(name, value);
        self.params.put(name, .{ .string = value }) catch return error.OutOfMemory;
    }

    /// Get an integer parameter value.  Returns `error.InvalidParameter`
    /// if the parameter does not exist or is not of integer type.
    pub fn getIntParam(self: Self, name: []const u8) ModelError!i64 {
        const v = self.params.get(name) orelse return error.InvalidParameter;
        return switch (v) {
            .int => |i| i,
            else => error.InvalidParameter,
        };
    }

    /// Get a double parameter value.
    pub fn getDblParam(self: Self, name: []const u8) ModelError!f64 {
        const v = self.params.get(name) orelse return error.InvalidParameter;
        return switch (v) {
            .double => |d| d,
            else => error.InvalidParameter,
        };
    }

    /// Get a string parameter value.
    pub fn getStrParam(self: Self, name: []const u8) ModelError![]const u8 {
        const v = self.params.get(name) orelse return error.InvalidParameter;
        return switch (v) {
            .string => |s| s,
            else => error.InvalidParameter,
        };
    }

    // ── Logging helper ─────────────────────────────────────────────────

    pub fn log(self: Self, comptime fmt: []const u8, args: anytype) void {
        _ = self;
        std.debug.print(fmt ++ "\n", args);
    }

    // ── Version ────────────────────────────────────────────────────────

    /// Return the version of the zhighs library.
    pub fn version() Version {
        return Version{ .major = 0, .minor = 1, .technical = 0 };
    }

    // ── Error handling ─────────────────────────────────────────────────

    /// Retrieve the last error message. Returns an empty string when no
    /// error has been recorded.
    pub fn getErrorMessage(self: Self) []const u8 {
        return if (self.last_error) |e| e else "";
    }

    // ── Parameter file I/O ─────────────────────────────────────────────

    /// Write all current parameters to a file.
    /// Format: one `param_name value` per line.
    pub fn writeParams(self: Self, filename: []const u8) ModelError!void {
        const file = std.fs.cwd().createFile(filename, .{}) catch return error.IoError;
        defer file.close();
        var it = self.params.iterator();
        while (it.next()) |entry| {
            const line = switch (entry.value_ptr.*) {
                .int => |v| std.fmt.allocPrint(self.allocator, "{s} {d}\n", .{ entry.key_ptr.*, v }) catch return error.OutOfMemory,
                .double => |v| std.fmt.allocPrint(self.allocator, "{s} {e}\n", .{ entry.key_ptr.*, v }) catch return error.OutOfMemory,
                .string => |v| std.fmt.allocPrint(self.allocator, "{s} {s}\n", .{ entry.key_ptr.*, v }) catch return error.OutOfMemory,
            };
            defer self.allocator.free(line);
            file.writeAll(line) catch return error.IoError;
        }
    }

    /// Read parameters from a file.
    /// Format: one `param_name value` per line.
    pub fn readParams(self: *Self, filename: []const u8) ModelError!void {
        const file = std.fs.cwd().openFile(filename, .{}) catch return error.IoError;
        defer file.close();
        const content = file.readToEndAlloc(self.allocator, 1 << 20) catch return error.IoError;
        defer self.allocator.free(content);

        var lines = std.mem.tokenizeScalar(u8, content, '\n');
        while (lines.next()) |line| {
            var parts = std.mem.tokenizeScalar(u8, line, ' ');
            const name = parts.next() orelse continue;
            const val_str = parts.next() orelse continue;

            // Determine parameter type by looking it up in current params.
            if (self.params.get(name)) |existing| {
                switch (existing) {
                    .int => {
                        const v = std.fmt.parseInt(i64, val_str, 10) catch return error.InvalidParameter;
                        try self.setIntParam(name, v);
                    },
                    .double => {
                        const v = std.fmt.parseFloat(f64, val_str) catch return error.InvalidParameter;
                        try self.setDblParam(name, v);
                    },
                    .string => try self.setStrParam(name, val_str),
                }
            }
        }
    }

    /// Reset all parameters to their default values.
    pub fn resetParams(self: *Self) ModelError!void {
        self.params.clearRetainingCapacity();
        try self.setDefaults();
    }

    // ── Type-dispatching param get/set ──────────────────────────────────

    /// Set a parameter whose type is determined at runtime.
    pub fn setParam(self: *Self, name: []const u8, value: ParamValue) ModelError!void {
        switch (value) {
            .int => |v| try self.setIntParam(name, v),
            .double => |v| try self.setDblParam(name, v),
            .string => |v| try self.setStrParam(name, v),
        }
    }

    // ── Callback support ───────────────────────────────────────────────

    /// Register a callback function.
    pub fn setCallbackFunc(self: *Self, callback: CallbackFunc, usrstate: ?*anyopaque) void {
        self.callback = callback;
        self.usrstate = usrstate;
    }

    /// Register a log callback function.
    pub fn setLogCallbackFunc(self: *Self, callback: CallbackFunc, usrstate: ?*anyopaque) void {
        self.log_callback = callback;
        self.log_usrstate = usrstate;
    }

    // ── Internal helpers ───────────────────────────────────────────────

    /// Populate the environment with the default parameter values.
    fn setDefaults(self: *Self) ModelError!void {
        const defaults = [_]struct { name: []const u8, value: ParamValue }{
            // ── Algorithm control ──────────────────────────────────────
            .{ .name = "Method",           .value = .{ .int = -1 } },     // -1=auto
            .{ .name = "Threads",          .value = .{ .int = 0 } },      // 0=auto
            .{ .name = "Presolve",         .value = .{ .int = -1 } },     // -1=auto
            .{ .name = "DualReductions",   .value = .{ .int = 1 } },
            .{ .name = "InfUnbdInfo",      .value = .{ .int = 0 } },
            .{ .name = "NormAdjust",       .value = .{ .int = -1 } },
            .{ .name = "Crossover",        .value = .{ .int = 0 } },     // 0=auto
            .{ .name = "BarHomogeneous",   .value = .{ .int = 1 } },
            .{ .name = "BarOrder",         .value = .{ .int = -1 } },    // -1=auto
            .{ .name = "SimplexPricing",   .value = .{ .int = -1 } },

            // ── Tolerances ─────────────────────────────────────────────
            .{ .name = "FeasibilityTol",   .value = .{ .double = 1e-6 } },
            .{ .name = "OptimalityTol",    .value = .{ .double = 1e-6 } },
            .{ .name = "BarQCPConvTol",    .value = .{ .double = 1e-6 } },
            .{ .name = "MarkowitzTol",     .value = .{ .double = 0.5 } },
            .{ .name = "PSDTol",           .value = .{ .double = 1e-6 } },

            // ── Limits ─────────────────────────────────────────────────
            .{ .name = "TimeLimit",        .value = .{ .double = INFINITY } },
            .{ .name = "IterationLimit",   .value = .{ .int = std.math.maxInt(i64) } },
            .{ .name = "NodeLimit",        .value = .{ .int = std.math.maxInt(i64) } },
            .{ .name = "SolutionLimit",    .value = .{ .int = std.math.maxInt(i64) } },
            .{ .name = "GapLimit",         .value = .{ .double = 1e-10 } },
            .{ .name = "BarIterLimit",     .value = .{ .int = std.math.maxInt(i64) } },

            // ── Output ─────────────────────────────────────────────────
            .{ .name = "OutputFlag",       .value = .{ .int = 1 } },
            .{ .name = "LogToConsole",     .value = .{ .int = 1 } },
            .{ .name = "DisplayInterval",  .value = .{ .int = 5 } },
            .{ .name = "MIPGap",           .value = .{ .double = 1e-4 } },
            .{ .name = "MIPGapAbs",        .value = .{ .double = 1e-10 } },

            // ── LP-specific ────────────────────────────────────────────
            .{ .name = "LPWarmStart",      .value = .{ .int = 1 } },
            .{ .name = "ObjScale",         .value = .{ .int = 0 } },

            // ── QP-specific ────────────────────────────────────────────
            .{ .name = "QCPDual",          .value = .{ .int = 0 } },
        };

        for (&defaults) |d| {
            try self.params.put(d.name, d.value);
        }
    }
};

// ── Parameter validation helpers ────────────────────────────────────────

fn validateIntParam(name: []const u8, value: i64) ModelError!void {
    _ = name;
    _ = value;
    // Future: range-check per parameter name.
}

fn validateDblParam(name: []const u8, value: f64) ModelError!void {
    _ = name;
    if (!std.math.isFinite(value)) return error.InvalidArgument;
}

fn validateStrParam(name: []const u8, value: []const u8) ModelError!void {
    _ = name;
    _ = value;
}

// ── Tests ───────────────────────────────────────────────────────────────

test "Env.initSimple creates a usable environment" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();

    try std.testing.expectEqual(@as(i64, -1), try env.getIntParam("Method"));
    try std.testing.expect(std.math.isInf(try env.getDblParam("TimeLimit")));
}

test "Env.setIntParam round-trips a parameter" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();

    try env.setIntParam("Threads", 4);
    try std.testing.expectEqual(@as(i64, 4), try env.getIntParam("Threads"));
}

test "Env.setDblParam round-trips a parameter" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();

    try env.setDblParam("FeasibilityTol", 1e-8);
    try std.testing.expectEqual(@as(f64, 1e-8), try env.getDblParam("FeasibilityTol"));
}

test "Env.getIntParam returns error on unknown parameter" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();

    try std.testing.expectError(error.InvalidParameter, env.getIntParam("NonExistentParam"));
}

test "Env.setDblParam rejects non-finite values" {
    var env = try Env.initSimple(std.testing.allocator);
    defer env.deinit();

    try std.testing.expectError(error.InvalidArgument, env.setDblParam("TimeLimit", std.math.nan(f64)));
    try std.testing.expectError(error.InvalidArgument, env.setDblParam("TimeLimit", std.math.inf(f64)));
}
