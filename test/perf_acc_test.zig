const std = @import("std");

const sentinel: f64 = std.math.floatMin(f64);

const Accum = struct {
    dimension: usize,
    dense_values: []f64,
    active_ptr: [*]u32,
    active_len: usize,
    active_cap: usize,

    fn init(alloc: std.mem.Allocator, dim: usize) !Accum {
        const dv = try alloc.alloc(f64, dim);
        @memset(dv, 0);
        return .{ .dimension = dim, .dense_values = dv, .active_ptr = undefined, .active_len = 0, .active_cap = 0 };
    }
    fn deinit(self: *Accum, alloc: std.mem.Allocator) void {
        alloc.free(self.dense_values);
        if (self.active_cap > 0) alloc.free(self.active_ptr[0..self.active_cap]);
    }
    fn reserve(self: *Accum, alloc: std.mem.Allocator, cap: usize) !void {
        if (cap > self.active_cap) {
            const buf = try alloc.alloc(u32, cap);
            if (self.active_cap > 0) { @memcpy(buf[0..self.active_len], self.active_ptr[0..self.active_len]); alloc.free(self.active_ptr[0..self.active_cap]); }
            self.active_ptr = buf.ptr; self.active_cap = cap;
        }
    }
    fn clear(self: *Accum) void {
        if (10 * self.active_len < 3 * self.dimension) { var i: usize = 0; while (i < self.active_len) : (i += 1) self.dense_values[self.active_ptr[i]] = 0.0; }
        else { @memset(self.dense_values, 0); }
        self.active_len = 0;
    }
    fn addAssumeValid(self: *Accum, id: usize, value: f64) void {
        if (value == 0.0) return;
        const dv = self.dense_values;
        const curr = dv[id];
        if (curr != 0.0) {
            const sum = curr + value;
            dv[id] = if (sum == 0.0) sentinel else sum;
        } else {
            dv[id] = value;
            self.active_ptr[self.active_len] = @intCast(id);
            self.active_len += 1;
        }
    }
};

pub fn main() !void {
    const alloc = std.heap.page_allocator;
    const dim: usize = 50000;
    const repeats: usize = 3000;
    var acc = try Accum.init(alloc, dim);
    defer acc.deinit(alloc);
    try acc.reserve(alloc, dim);
    for (0..repeats) |_| {
        acc.clear();
        var i: usize = 0;
        while (i < dim) : (i += 1) {
            acc.addAssumeValid(i, 1.0);
            acc.addAssumeValid(i, -0.5);
        }
    }
    std.debug.print("{d}\n", .{acc.dense_values[dim / 2]});
}
