//! Primal/dual ratio-test policies.

const std = @import("std");

pub const RatioRule = enum { standard, harris_two_pass, bound_flipping };
pub const LeavingChoice = struct { row: ?u32 = null, step: f64 = 0.0 };
pub const RatioTest = struct {
    rule: RatioRule = .harris_two_pass,
    tolerance: f64 = 1e-9,

    pub fn chooseLeaving(self: *const RatioTest, direction: []const f64, rhs: []const f64) LeavingChoice {
        var choice = LeavingChoice{ .step = std.math.inf(f64) };
        for (direction, rhs, 0..) |coefficient, value, i| {
            if (coefficient > self.tolerance) {
                const step = value / coefficient;
                if (step >= -self.tolerance and step < choice.step) {
                    choice = .{ .row = @intCast(i), .step = @max(step, 0.0) };
                }
            }
        }
        return choice;
    }
};

test {
    std.testing.refAllDecls(@This());
}
