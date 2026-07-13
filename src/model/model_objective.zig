//! Objective-specific methods for `Model`.
//!
//! ## Responsibility
//!
//! Owns objective representations that need dedicated model storage, currently
//! piecewise-linear objectives.  Linear coefficients are handled through
//! `model_linear.zig`, quadratic terms through `model_constraints.zig`, and
//! solver-wide multi-objective workflows through `model_advanced.zig`.

const types = @import("types.zig");
const Model = @import("model.zig").Model;

const ModelError = types.ModelError;

// ══════════════════════════════════════════════════════════════════════════
//  Piecewise-linear objective
// ══════════════════════════════════════════════════════════════════════════

/// Set a piecewise-linear objective for a variable.
/// `(x[i], y[i])` are the breakpoints; NUM_PTS must be ≥ 2.
pub fn setPWLObj(self: *Model, var_idx: usize, num_pts: usize, x: []const f64, y: []const f64) ModelError!void {
    const alloc = self.allocator;
    if (x.len < num_pts or y.len < num_pts or num_pts < 2) return error.InvalidArgument;

    // Check if this variable already has a PWL objective — replace it.
    for (self.pwlobj_var, 0..) |v, i| {
        if (v == var_idx) {
            // Replace existing entry.
            const old_npts = self.pwlobj_npts[i];
            // Compute offsets into packed data.
            var x_off: usize = 0;
            var y_off: usize = 0;
            for (0..i) |j| {
                x_off += self.pwlobj_npts[j];
                y_off += self.pwlobj_npts[j];
            }
            // Replace the x/y data in-place (may change length).
            const diff = @as(isize, @intCast(num_pts)) - @as(isize, @intCast(old_npts));
            if (diff > 0) {
                const grow = @as(usize, @intCast(diff));
                const old_xlen = self.pwlobj_xdata.len;
                self.pwlobj_xdata = try alloc.realloc(self.pwlobj_xdata, old_xlen + grow);
                self.pwlobj_ydata = try alloc.realloc(self.pwlobj_ydata, old_xlen + grow);
                // Shift trailing entries.
                if (y_off + old_npts < old_xlen) {
                    const tail = old_xlen - (y_off + old_npts);
                    @memcpy(self.pwlobj_xdata[y_off + num_pts ..][0..tail], self.pwlobj_xdata[y_off + old_npts ..][0..tail]);
                    @memcpy(self.pwlobj_ydata[y_off + num_pts ..][0..tail], self.pwlobj_ydata[y_off + old_npts ..][0..tail]);
                }
            } else if (diff < 0) {
                const shrink = @as(usize, @intCast(-diff));
                const old_xlen = self.pwlobj_xdata.len;
                if (y_off + old_npts < old_xlen) {
                    const tail = old_xlen - (y_off + old_npts);
                    @memcpy(self.pwlobj_xdata[y_off + num_pts ..][0..tail], self.pwlobj_xdata[y_off + old_npts ..][0..tail]);
                    @memcpy(self.pwlobj_ydata[y_off + num_pts ..][0..tail], self.pwlobj_ydata[y_off + old_npts ..][0..tail]);
                }
                self.pwlobj_xdata = alloc.realloc(self.pwlobj_xdata, old_xlen - shrink) catch unreachable;
                self.pwlobj_ydata = alloc.realloc(self.pwlobj_ydata, old_xlen - shrink) catch unreachable;
            }
            @memcpy(self.pwlobj_xdata[x_off..][0..num_pts], x[0..num_pts]);
            @memcpy(self.pwlobj_ydata[y_off..][0..num_pts], y[0..num_pts]);
            self.pwlobj_npts[i] = num_pts;
            self.revision += 1;
            return;
        }
    }

    // New entry.
    const old = self.pwlobj_count;
    const new = old + 1;
    self.pwlobj_var = try alloc.realloc(self.pwlobj_var, new);
    self.pwlobj_var[old] = var_idx;
    self.pwlobj_npts = try alloc.realloc(self.pwlobj_npts, new);
    self.pwlobj_npts[old] = num_pts;
    const old_xlen = self.pwlobj_xdata.len;
    const new_xlen = old_xlen + num_pts;
    self.pwlobj_xdata = try alloc.realloc(self.pwlobj_xdata, new_xlen);
    self.pwlobj_ydata = try alloc.realloc(self.pwlobj_ydata, new_xlen);
    @memcpy(self.pwlobj_xdata[old_xlen..new_xlen], x[0..num_pts]);
    @memcpy(self.pwlobj_ydata[old_xlen..new_xlen], y[0..num_pts]);
    self.pwlobj_count = new;
    self.revision += 1;
}

/// Retrieve piecewise-linear objective data for a variable.
/// Returns the number of points actually available for this variable.
pub fn getPWLObj(self: Model, var_idx: usize, num_pts: *usize, x: []f64, y: []f64) ModelError!usize {
    for (self.pwlobj_var, 0..) |v, i| {
        if (v == var_idx) {
            const np = self.pwlobj_npts[i];
            num_pts.* = np;
            // Compute offset into packed data.
            var off: usize = 0;
            for (0..i) |j| off += self.pwlobj_npts[j];
            if (x.len >= np and y.len >= np) {
                @memcpy(x[0..np], self.pwlobj_xdata[off..][0..np]);
                @memcpy(y[0..np], self.pwlobj_ydata[off..][0..np]);
            }
            return np;
        }
    }
    return error.NotInModel;
}
