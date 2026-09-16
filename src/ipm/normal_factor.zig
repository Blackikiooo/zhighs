const std = @import("std");
const Dense = @import("cholesky.zig").DenseCholesky;
const Profile = @import("profile_cholesky.zig").ProfileCholesky;
const StandardForm = @import("standard_form.zig").StandardForm;

pub const NormalFactor = union(enum) {
    dense: Dense,
    profile: Profile,

    pub fn init(a: std.mem.Allocator, sf: *const StandardForm, dense: bool) !NormalFactor {
        return if (dense) .{ .dense = try Dense.init(a, sf.num_rows) } else .{ .profile = try Profile.init(a, sf) };
    }
    pub fn deinit(self: *NormalFactor) void {
        switch (self.*) {
            inline else => |*value| value.deinit(),
        }
    }
    pub fn storageLen(self: *const NormalFactor) usize {
        return switch (self.*) {
            inline else => |value| value.data.len,
        };
    }
    pub fn index(self: *const NormalFactor, row: usize, col: usize) usize {
        return switch (self.*) {
            .dense => |value| row * value.n + col,
            .profile => |*value| value.index(row, col),
        };
    }
    pub fn factor(self: *NormalFactor, values: []const f64) error{NotPositiveDefinite}!void {
        switch (self.*) {
            inline else => |*value| try value.factor(values),
        }
    }
    pub fn solve(self: *const NormalFactor, rhs: []f64) void {
        switch (self.*) {
            inline else => |*value| value.solve(rhs),
        }
    }
};
