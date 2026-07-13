//! Entering-variable pricing policies.

const std = @import("std");

pub const PricingRule = enum { dantzig, devex, steepest_edge, partial, hyper_sparse };
pub const Pricing = struct {
    rule: PricingRule = .devex,
    devex_reset_period: usize = 100,
    iterations: usize = 0,

    pub fn chooseEntering(self: *Pricing, reduced_cost: []const f64, tolerance: f64) ?u32 {
        self.iterations += 1;
        var best: ?u32 = null;
        var best_value = tolerance;
        for (reduced_cost, 0..) |value, i| {
            const magnitude = @abs(value);
            if (magnitude > best_value) {
                best_value = magnitude;
                best = @intCast(i);
            }
        }
        return best;
    }
};

test {
    std.testing.refAllDecls(@This());
}
