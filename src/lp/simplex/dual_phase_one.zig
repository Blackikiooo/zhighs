//! Persistent SoA work state for the serial dual simplex.
//!
//! BasisState owns the live work bounds/value/dual arrays. This structure owns
//! the complementary cost, shift, range and move arrays plus immutable bound
//! snapshots. It survives Phase I/II transitions and never copies the model
//! matrix.

const std = @import("std");
const basis_module = @import("basis.zig");
const highs_random = @import("highs_random.zig");

/// Changes made while correcting nonbasic dual infeasibilities at rebuild.
pub const DualCorrection = struct {
    /// Free variables that cannot be corrected by flipping or shifting.
    free_infeasibility_count: usize = 0,
    /// Fixed/boxed nonbasic variables moved to their opposite bound.
    flip_count: usize = 0,
    /// One-sided reduced costs corrected by temporary cost shifts.
    shift_count: usize = 0,
};

/// Aggregate count, maximum and sum of one feasibility measure.
pub const InfeasibilityStats = struct {
    /// Entries at or above the supplied feasibility tolerance.
    count: usize = 0,
    /// Largest positive infeasibility, including values below tolerance.
    maximum: f64 = 0.0,
    /// Sum of all positive infeasibilities.
    sum: f64 = 0.0,
};

/// Non-owning binding of dual work arrays used by iteration kernels.
pub const DualWorkView = struct {
    /// Working objective coefficients before temporary `shift`.
    cost: []f64,
    /// Temporary per-column objective correction.
    shift: []f64,
    /// Live working lower bounds owned by `BasisState`.
    lower: []f64,
    /// Live working upper bounds owned by `BasisState`.
    upper: []f64,
    /// Cached `upper-lower` range used by BFRT.
    range: []f64,
    /// Live primal values owned by `BasisState`.
    value: []f64,
    /// Live reduced costs/work duals owned by `BasisState`.
    dual: []f64,
    /// Persistent nonbasic move in {-1,0,+1}.
    move: []i8,
};

/// Owning persistent SoA workspace shared by serial dual Phase I and II.
pub const DualPhaseOneWorkspace = struct {
    /// Allocator owning every retained buffer below.
    allocator: std.mem.Allocator,
    /// Original/Phase-II lower bounds saved before installing Phase-I bounds.
    saved_lower: []f64 = &.{},
    /// Original/Phase-II upper bounds saved before installing Phase-I bounds.
    saved_upper: []f64 = &.{},
    /// Working objective coefficient for every non-artificial column.
    work_cost: []f64 = &.{},
    /// Temporary objective shifts kept separate from `work_cost`.
    work_shift: []f64 = &.{},
    /// Cached working bound width for BFRT capacity calculations.
    work_range: []f64 = &.{},
    /// Per-column scratch for dual infeasibility construction.
    dual_infeasibility: []f64 = &.{},
    /// Deterministic random values used by cost perturbation.
    perturbation: []f64 = &.{},
    /// HiGHS random permutation used for CHUZR/CHUZC tie order.
    highs_permutation: []u32 = &.{},
    /// Explicit Phase-I nonbasic direction: +1 moves up from the lower
    /// bound, -1 moves down from the upper bound, and zero is basic/fixed/free.
    nonbasic_move: []i8 = &.{},
    /// Residual basic-bound violation after the current ratio-test flip set.
    remaining_violation: []f64 = &.{},
    /// `HEkkDualRHS` state. Infeasibilities are squared, matching HiGHS when
    /// dual steepest-edge pricing is active. The index/mark pair deliberately
    /// retains feasible rows until the next rebuild: insertion order is the
    /// CHUZR tie-break state, not merely a cache of the current candidates.
    primal_infeasibility: []f64 = &.{},
    /// Sparse row order used by serial CHUZR.
    infeasibility_index: []u32 = &.{},
    /// Whether a row has ever been inserted since the last rebuild.
    infeasibility_mark: []u8 = &.{},
    /// HVector-equivalent row order changed by the latest FTRAN.
    changed_row_index: []u32 = &.{},
    /// Active prefix length of `infeasibility_index`.
    infeasibility_count: usize = 0,
    /// Active prefix length of `changed_row_index`.
    changed_row_count: usize = 0,
    /// Hyper-sparse merit cutoff; zero until the >500 candidate path exists.
    infeasibility_cutoff: f64 = 0.0,
    /// True when CHUZR scans every row instead of the sparse index.
    infeasibility_dense: bool = false,
    /// Whether the persistent RHS list is initialized for the current basis.
    infeasibility_list_active: bool = false,
    /// Reusable transactional checkpoint for the basis entering the shifted
    /// cold-dual path. Status is one contiguous SoA block; the basis head is
    /// dense. These buffers grow with the existing workspace and never
    /// allocate on a steady-state solve.
    checkpoint_status: []basis_module.BasisStatus = &.{},
    /// Basis head paired with `checkpoint_status`.
    checkpoint_basic_index: []u32 = &.{},
    /// Whether checkpoint prefixes contain a complete restorable basis.
    checkpoint_valid: bool = false,
    /// Monotonic identifier incremented when a Phase-I/cost epoch begins.
    basis_epoch: u64 = 0,
    /// HiGHS-compatible RNG stream consumed by dual correction and CHUZR.
    correction_random: highs_random.RandomStream = highs_random.RandomStream.init(0),
    /// Basic rows infeasible at the most recently recorded rebuild.
    rebuild_primal_infeasibility_count: usize = 0,
    /// Largest basic-bound violation at the most recent rebuild.
    rebuild_primal_infeasibility_max: f64 = 0.0,
    /// Sum of basic-bound violations at the most recent rebuild.
    rebuild_primal_infeasibility_sum: f64 = 0.0,
    /// Working dual objective at the most recent rebuild.
    rebuild_dual_objective: f64 = 0.0,
    /// Monotonic number of recorded rebuild solution states.
    rebuild_epoch: u64 = 0,
    /// Basic-coordinate envelope in the engine's scaled coordinates when the
    /// current Phase-I epoch was installed.
    working_radius: f64 = 1.0,
    /// Whether Phase-I bounds are currently installed and require restoration.
    active: bool = false,

    /// Construct an empty workspace whose arrays grow on first binding.
    pub fn init(allocator: std.mem.Allocator) DualPhaseOneWorkspace {
        return .{ .allocator = allocator };
    }

    /// Release every retained SoA buffer.
    pub fn deinit(self: *DualPhaseOneWorkspace) void {
        self.allocator.free(self.saved_lower);
        self.allocator.free(self.saved_upper);
        self.allocator.free(self.work_cost);
        self.allocator.free(self.work_shift);
        self.allocator.free(self.work_range);
        self.allocator.free(self.dual_infeasibility);
        self.allocator.free(self.perturbation);
        self.allocator.free(self.highs_permutation);
        self.allocator.free(self.nonbasic_move);
        self.allocator.free(self.remaining_violation);
        self.allocator.free(self.primal_infeasibility);
        self.allocator.free(self.infeasibility_index);
        self.allocator.free(self.infeasibility_mark);
        self.allocator.free(self.changed_row_index);
        self.allocator.free(self.checkpoint_status);
        self.allocator.free(self.checkpoint_basic_index);
        self.* = .{ .allocator = self.allocator };
    }

    /// Grow all SoA arrays together. Existing capacity is retained across
    /// solves, so phase transitions allocate only when the model grows.
    pub fn ensureCapacity(self: *DualPhaseOneWorkspace, count: usize) !void {
        if (self.saved_lower.len >= count) return;
        const lower = try self.allocator.alloc(f64, count);
        errdefer self.allocator.free(lower);
        const upper = try self.allocator.alloc(f64, count);
        errdefer self.allocator.free(upper);
        const cost = try self.allocator.alloc(f64, count);
        errdefer self.allocator.free(cost);
        const shift = try self.allocator.alloc(f64, count);
        errdefer self.allocator.free(shift);
        const range = try self.allocator.alloc(f64, count);
        errdefer self.allocator.free(range);
        const infeasibility = try self.allocator.alloc(f64, count);
        errdefer self.allocator.free(infeasibility);
        const perturbation = try self.allocator.alloc(f64, count);
        errdefer self.allocator.free(perturbation);
        const highs_permutation = try self.allocator.alloc(u32, count);
        errdefer self.allocator.free(highs_permutation);
        const move = try self.allocator.alloc(i8, count);
        errdefer self.allocator.free(move);
        const remaining = try self.allocator.alloc(f64, count);
        errdefer self.allocator.free(remaining);
        const primal_infeasibility = try self.allocator.alloc(f64, count);
        errdefer self.allocator.free(primal_infeasibility);
        const infeasibility_index = try self.allocator.alloc(u32, count);
        errdefer self.allocator.free(infeasibility_index);
        const infeasibility_mark = try self.allocator.alloc(u8, count);
        errdefer self.allocator.free(infeasibility_mark);
        const changed_row_index = try self.allocator.alloc(u32, count);
        errdefer self.allocator.free(changed_row_index);
        const checkpoint_status = try self.allocator.alloc(basis_module.BasisStatus, count);
        errdefer self.allocator.free(checkpoint_status);
        const checkpoint_basic_index = try self.allocator.alloc(u32, count);
        errdefer self.allocator.free(checkpoint_basic_index);

        self.allocator.free(self.saved_lower);
        self.allocator.free(self.saved_upper);
        self.allocator.free(self.work_cost);
        self.allocator.free(self.work_shift);
        self.allocator.free(self.work_range);
        self.allocator.free(self.dual_infeasibility);
        self.allocator.free(self.perturbation);
        self.allocator.free(self.highs_permutation);
        self.allocator.free(self.nonbasic_move);
        self.allocator.free(self.remaining_violation);
        self.allocator.free(self.primal_infeasibility);
        self.allocator.free(self.infeasibility_index);
        self.allocator.free(self.infeasibility_mark);
        self.allocator.free(self.changed_row_index);
        self.allocator.free(self.checkpoint_status);
        self.allocator.free(self.checkpoint_basic_index);
        self.saved_lower = lower;
        self.saved_upper = upper;
        self.work_cost = cost;
        self.work_shift = shift;
        self.work_range = range;
        self.dual_infeasibility = infeasibility;
        self.perturbation = perturbation;
        self.highs_permutation = highs_permutation;
        self.nonbasic_move = move;
        self.remaining_violation = remaining;
        self.primal_infeasibility = primal_infeasibility;
        self.infeasibility_index = infeasibility_index;
        self.infeasibility_mark = infeasibility_mark;
        self.changed_row_index = changed_row_index;
        self.checkpoint_status = checkpoint_status;
        self.checkpoint_basic_index = checkpoint_basic_index;
        self.checkpoint_valid = false;
    }

    /// Return squared basic-bound violation, or zero inside tolerance.
    fn squaredPrimalInfeasibility(value: f64, lower: f64, upper: f64, tolerance: f64) f64 {
        const violation = if (value < lower - tolerance)
            lower - value
        else if (value > upper + tolerance)
            value - upper
        else
            0.0;
        return violation * violation;
    }

    /// Rebuild `HEkkDualRHS::work_infeasibility/workIndex/workMark` in row
    /// order. The hyper-sparse cutoff is intentionally left at zero until its
    /// >500-candidate selection workspace is migrated; dense-vs-sparse mode
    /// and the 20% boundary are exact for the serial corpus path.
    pub fn createPrimalInfeasibilityList(
        self: *DualPhaseOneWorkspace,
        basis: *const basis_module.BasisState,
        num_rows: usize,
        tolerance: f64,
    ) bool {
        if (num_rows > basis.basic_value.len or num_rows > self.primal_infeasibility.len or
            num_rows > self.infeasibility_index.len or num_rows > self.infeasibility_mark.len) return false;
        @memset(self.infeasibility_mark[0..num_rows], 0);
        self.infeasibility_count = 0;
        self.infeasibility_cutoff = 0.0;
        self.infeasibility_dense = false;
        for (0..num_rows) |row| {
            const infeasibility = squaredPrimalInfeasibility(
                basis.basic_value[row],
                basis.basic_lower[row],
                basis.basic_upper[row],
                tolerance,
            );
            self.primal_infeasibility[row] = infeasibility;
            if (infeasibility == 0.0) continue;
            self.infeasibility_mark[row] = 1;
            self.infeasibility_index[self.infeasibility_count] = @intCast(row);
            self.infeasibility_count += 1;
        }
        if (self.infeasibility_count * 5 > num_rows) self.infeasibility_dense = true;
        self.infeasibility_list_active = true;
        return true;
    }

    /// Update values for the FTRAN output rows, then append newly infeasible
    /// rows in the solver-provided HVector-equivalent order. Existing marks
    /// are never cleared here, matching `HEkkDualRHS::updateInfeasList`.
    pub fn updatePrimalInfeasibilityList(
        self: *DualPhaseOneWorkspace,
        basis: *const basis_module.BasisState,
        rows: []const u32,
        tolerance: f64,
    ) void {
        if (!self.infeasibility_list_active) return;
        for (rows) |row_u32| {
            const row: usize = @intCast(row_u32);
            if (row >= basis.basic_value.len or row >= self.primal_infeasibility.len) continue;
            const infeasibility = squaredPrimalInfeasibility(
                basis.basic_value[row],
                basis.basic_lower[row],
                basis.basic_upper[row],
                tolerance,
            );
            self.primal_infeasibility[row] = infeasibility;
            if (self.infeasibility_dense or infeasibility == 0.0 or self.infeasibility_mark[row] != 0) continue;
            if (self.infeasibility_count >= self.infeasibility_index.len) continue;
            self.infeasibility_index[self.infeasibility_count] = row_u32;
            self.infeasibility_mark[row] = 1;
            self.infeasibility_count += 1;
        }
    }

    /// Choose a CHUZR row using HiGHS list order and weighted squared merit.
    pub fn choosePrimalInfeasibleRow(
        self: *DualPhaseOneWorkspace,
        basis: *const basis_module.BasisState,
        num_rows: usize,
    ) ?u32 {
        if (!self.infeasibility_list_active or num_rows == 0) return null;
        const scan_count = if (self.infeasibility_dense) num_rows else self.infeasibility_count;
        if (scan_count == 0) return null;
        const random_start = self.correction_random.integer(scan_count);
        var best_row: ?u32 = null;
        var best_merit: f64 = 0.0;
        for (0..scan_count) |offset| {
            const position = (random_start + offset) % scan_count;
            const row_u32: u32 = if (self.infeasibility_dense) @intCast(position) else self.infeasibility_index[position];
            const row: usize = @intCast(row_u32);
            if (row >= basis.row_edge_weight.len) continue;
            const infeasibility = self.primal_infeasibility[row];
            if (infeasibility <= 0.0) continue;
            const weight = basis.row_edge_weight[row];
            if (best_merit * weight < infeasibility) {
                best_merit = infeasibility / weight;
                best_row = row_u32;
            }
        }
        // A positive cutoff can require a recursive full-list rebuild in
        // HiGHS. The cutoff path is not enabled until its >500-row selection
        // workspace is present, so cutoff is currently always zero.
        return best_row;
    }

    /// Save status and basis-head prefixes for transactional fallback.
    pub fn captureBasisCheckpoint(
        self: *DualPhaseOneWorkspace,
        basis: *const basis_module.BasisState,
        original_count: usize,
        num_rows: usize,
    ) !void {
        try self.ensureCapacity(original_count);
        if (num_rows > self.checkpoint_basic_index.len) return error.OutOfMemory;
        @memcpy(self.checkpoint_status[0..original_count], basis.col_status[0..original_count]);
        @memcpy(self.checkpoint_basic_index[0..num_rows], basis.basic_index[0..num_rows]);
        self.checkpoint_valid = true;
    }

    /// Bind the complementary SoA storage to BasisState's live work arrays.
    /// No allocation or copy occurs: lower/upper/value/dual remain the single
    /// arrays used by primal and dual kernels.
    pub fn view(self: *DualPhaseOneWorkspace, basis: *basis_module.BasisState, count: usize) ?DualWorkView {
        if (self.work_cost.len < count or self.work_shift.len < count or self.work_range.len < count or
            self.nonbasic_move.len < count or basis.col_lower.len < count or basis.col_upper.len < count or
            basis.primal.len < count or basis.reduced_cost.len < count) return null;
        return .{
            .cost = self.work_cost[0..count],
            .shift = self.work_shift[0..count],
            .lower = basis.col_lower[0..count],
            .upper = basis.col_upper[0..count],
            .range = self.work_range[0..count],
            .value = basis.primal[0..count],
            .dual = basis.reduced_cost[0..count],
            .move = self.nonbasic_move[0..count],
        };
    }

    /// Rebuild explicit nonbasic moves and snap values to compatible bounds.
    fn initializeNonbasicMove(self: *DualPhaseOneWorkspace, basis: *basis_module.BasisState, count: usize) void {
        for (0..count) |column| {
            const lower = basis.col_lower[column];
            const upper = basis.col_upper[column];
            if (basis.col_status[column] == .basic) {
                self.nonbasic_move[column] = 0;
                continue;
            }
            // HiGHS derives fixed/free movement from work bounds, not from a
            // possibly stale imported status label.
            if (lower == upper) {
                self.nonbasic_move[column] = 0;
                basis.col_status[column] = .fixed;
                basis.primal[column] = lower;
                continue;
            }
            if (!std.math.isFinite(lower) and !std.math.isFinite(upper)) {
                self.nonbasic_move[column] = 0;
                basis.col_status[column] = .free;
                basis.primal[column] = 0.0;
                continue;
            }
            if (std.math.isFinite(lower)) {
                if (std.math.isFinite(upper) and basis.col_status[column] == .at_upper) {
                    self.nonbasic_move[column] = -1;
                    basis.col_status[column] = .at_upper;
                    basis.primal[column] = upper;
                } else {
                    // Lower or boxed with invalid original move.
                    self.nonbasic_move[column] = 1;
                    basis.col_status[column] = .at_lower;
                    basis.primal[column] = lower;
                }
            } else {
                self.nonbasic_move[column] = -1;
                basis.col_status[column] = .at_upper;
                basis.primal[column] = upper;
            }
        }
    }

    /// Bind original Phase-II bounds and initialize complementary work arrays.
    pub fn bindPhaseTwo(self: *DualPhaseOneWorkspace, basis: *basis_module.BasisState, count: usize) !DualWorkView {
        try self.ensureCapacity(count);
        self.infeasibility_list_active = false;
        @memcpy(self.saved_lower[0..count], basis.col_lower[0..count]);
        @memcpy(self.saved_upper[0..count], basis.col_upper[0..count]);
        for (self.work_range[0..count], basis.col_lower[0..count], basis.col_upper[0..count]) |*range_value, lower, upper|
            range_value.* = upper - lower;
        @memset(self.work_shift[0..count], 0.0);
        self.initializeNonbasicMove(basis, count);
        return self.view(basis, count) orelse error.DimensionMismatch;
    }

    /// Position the correction RNG after both cold-solve random-vector initializations.
    pub fn resetCorrectionRandom(self: *DualPhaseOneWorkspace, num_cols: usize, count: usize) void {
        self.correction_random = highs_random.RandomStream.afterVectorInitialization(num_cols, count);
    }

    /// Match `computeDualInfeasibilitiesWithFixedVariableFlips`: fixed
    /// variables have move zero and are ignored; free variables use |dual|.
    pub fn dualInfeasibilityStats(
        self: *const DualPhaseOneWorkspace,
        basis: *const basis_module.BasisState,
        count: usize,
        dual_tolerance: f64,
    ) InfeasibilityStats {
        var stats = InfeasibilityStats{};
        for (0..count) |column| {
            if (basis.col_status[column] == .basic) continue;
            const lower = basis.col_lower[column];
            const upper = basis.col_upper[column];
            const dual = basis.reduced_cost[column];
            const infeasibility = if (!std.math.isFinite(lower) and !std.math.isFinite(upper))
                @abs(dual)
            else
                -@as(f64, @floatFromInt(self.nonbasic_move[column])) * dual;
            if (infeasibility <= 0.0) continue;
            if (infeasibility >= dual_tolerance) stats.count += 1;
            stats.maximum = @max(stats.maximum, infeasibility);
            stats.sum += infeasibility;
        }
        return stats;
    }

    /// Source-equivalent core of `HEkkDual::correctDualInfeasibilities`.
    /// Free infeasibilities are reported, fixed/boxed variables are flipped
    /// when allowed, and all remaining one-sided infeasibilities are removed
    /// by changing workCost and the matching live workDual.
    pub fn correctDualInfeasibilities(
        self: *DualPhaseOneWorkspace,
        basis: *basis_module.BasisState,
        count: usize,
        dual_tolerance: f64,
        force_phase_two: bool,
    ) DualCorrection {
        var result = DualCorrection{};
        for (0..count) |column| {
            if (basis.col_status[column] == .basic) continue;
            const lower = basis.col_lower[column];
            const upper = basis.col_upper[column];
            const current_dual = basis.reduced_cost[column];
            const move = self.nonbasic_move[column];
            const free = !std.math.isFinite(lower) and !std.math.isFinite(upper);
            if (free) {
                if (@abs(current_dual) >= dual_tolerance) result.free_infeasibility_count += 1;
                continue;
            }
            const dual_infeasibility = -@as(f64, @floatFromInt(move)) * current_dual;
            if (dual_infeasibility < dual_tolerance) continue;
            const fixed = lower == upper;
            const boxed = std.math.isFinite(lower) and std.math.isFinite(upper);
            if (fixed or (boxed and !force_phase_two)) {
                self.nonbasic_move[column] = -move;
                if (self.nonbasic_move[column] > 0) {
                    basis.col_status[column] = .at_lower;
                    basis.primal[column] = lower;
                } else if (self.nonbasic_move[column] < 0) {
                    basis.col_status[column] = .at_upper;
                    basis.primal[column] = upper;
                } else {
                    basis.col_status[column] = .fixed;
                    basis.primal[column] = lower;
                }
                result.flip_count += 1;
                continue;
            }
            const new_dual = if (move > 0)
                (1.0 + self.correction_random.fraction()) * dual_tolerance
            else
                -(1.0 + self.correction_random.fraction()) * dual_tolerance;
            const shift = new_dual - current_dual;
            self.work_cost[column] += shift;
            basis.reduced_cost[column] = new_dual;
            result.shift_count += 1;
        }
        return result;
    }

    /// Materialize the state produced by HiGHS' rebuild after computePrimal:
    /// the row-indexed primal infeasibilities and the current dual objective.
    /// `remaining_violation` is reused as the persistent dense backing array;
    /// CHUZR may still maintain a sparse candidate structure independently.
    pub fn recordRebuildState(
        self: *DualPhaseOneWorkspace,
        basis: *const basis_module.BasisState,
        count: usize,
        primal_tolerance: f64,
    ) bool {
        if (self.remaining_violation.len < basis.basic_value.len or count > basis.col_status.len) return false;
        self.rebuild_primal_infeasibility_count = 0;
        self.rebuild_primal_infeasibility_max = 0.0;
        self.rebuild_primal_infeasibility_sum = 0.0;
        @memset(self.remaining_violation[0..basis.basic_value.len], 0.0);
        for (basis.basic_value, basis.basic_lower, basis.basic_upper, 0..) |value, lower, upper, row| {
            const violation = @max(@max(lower - value, value - upper), 0.0);
            self.remaining_violation[row] = violation;
            if (violation < primal_tolerance) continue;
            self.rebuild_primal_infeasibility_count += 1;
            self.rebuild_primal_infeasibility_max = @max(self.rebuild_primal_infeasibility_max, violation);
            self.rebuild_primal_infeasibility_sum += violation;
        }
        self.rebuild_dual_objective = self.dualObjective(basis, count);
        self.rebuild_epoch +%= 1;
        return true;
    }

    /// Install one temporary cost shift; returns false for invalid/double-shifted columns.
    pub fn shiftCost(self: *DualPhaseOneWorkspace, column: usize, amount: f64) bool {
        if (column >= self.work_shift.len or self.work_shift[column] != 0.0) return false;
        self.work_shift[column] = amount;
        return true;
    }

    /// Remove a column's temporary shift from its live reduced cost.
    pub fn shiftBack(self: *DualPhaseOneWorkspace, dual: []f64, column: usize) bool {
        if (column >= self.work_shift.len or column >= dual.len) return false;
        dual[column] -= self.work_shift[column];
        self.work_shift[column] = 0.0;
        return true;
    }

    /// Start a dual Phase-I epoch and snapshot original bounds for restoration.
    pub fn begin(self: *DualPhaseOneWorkspace, basis: *basis_module.BasisState, original_count: usize) !void {
        try self.ensureCapacity(original_count);
        self.infeasibility_list_active = false;
        @memcpy(self.saved_lower[0..original_count], basis.col_lower[0..original_count]);
        @memcpy(self.saved_upper[0..original_count], basis.col_upper[0..original_count]);
        for (self.work_range[0..original_count], basis.col_lower[0..original_count], basis.col_upper[0..original_count]) |*range_value, lower, upper|
            range_value.* = upper - lower;
        @memset(self.work_shift[0..original_count], 0.0);
        @memset(self.dual_infeasibility[0..original_count], 0.0);
        @memset(self.perturbation[0..original_count], 0.0);
        self.initializeNonbasicMove(basis, original_count);
        @memset(self.remaining_violation[0..original_count], 0.0);
        self.basis_epoch +%= 1;
        self.active = true;
    }

    /// Start a work-cost-only epoch. Bounds are snapshotted for perturbation
    /// policy, but no Phase-I bounds are installed and no deferred restore is
    /// required.
    pub fn beginCostEpoch(self: *DualPhaseOneWorkspace, basis: *basis_module.BasisState, original_count: usize) !void {
        try self.ensureCapacity(original_count);
        self.infeasibility_list_active = false;
        @memcpy(self.saved_lower[0..original_count], basis.col_lower[0..original_count]);
        @memcpy(self.saved_upper[0..original_count], basis.col_upper[0..original_count]);
        for (self.work_range[0..original_count], basis.col_lower[0..original_count], basis.col_upper[0..original_count]) |*range_value, lower, upper|
            range_value.* = upper - lower;
        @memset(self.work_shift[0..original_count], 0.0);
        @memset(self.dual_infeasibility[0..original_count], 0.0);
        @memset(self.perturbation[0..original_count], 0.0);
        self.initializeNonbasicMove(basis, original_count);
        @memset(self.remaining_violation[0..original_count], 0.0);
        self.basis_epoch +%= 1;
        self.active = false;
    }

    /// Map original bounds to the dual Phase-I subproblem (HiGHS-style).
    ///
    /// Working bounds:
    ///   FREE           → [-1000, 1000]
    ///   lower-unbounded → [-1, 0]
    ///   upper-unbounded → [0, 1]
    ///   boxed / fixed  → [0, 0]    (collapsed)
    ///
    /// `initialiseNonbasicValueAndMove` is then applied to these Phase-I
    /// bounds using the original nonbasic move. Thus ordinary lower/upper
    /// variables start at the zero endpoint; they are not pre-shifted to
    /// ±1 from their reduced-cost sign.
    ///
    /// Boxed columns collapse to [0, 0] and have zero Phase-I move, exactly
    /// like HiGHS' initialiseNonbasicValueAndMove. Their original side is
    /// selected again only after the Phase-II bounds have been restored.
    pub fn installWorkingBounds(
        self: *DualPhaseOneWorkspace,
        basis: *basis_module.BasisState,
        original_count: usize,
    ) void {
        self.working_radius = 1000.0;

        for (0..original_count) |column| {
            const lower = self.saved_lower[column];
            const upper = self.saved_upper[column];
            const lower_finite = std.math.isFinite(lower);
            const upper_finite = std.math.isFinite(upper);

            // ── Set working bounds ──
            if (!lower_finite and !upper_finite) {
                basis.col_lower[column] = -1000.0;
                basis.col_upper[column] = 1000.0;
            } else if (!lower_finite) {
                basis.col_lower[column] = -1.0;
                basis.col_upper[column] = 0.0;
            } else if (!upper_finite) {
                basis.col_lower[column] = 0.0;
                basis.col_upper[column] = 1.0;
            } else {
                // Both finite: collapse to [0, 0] regardless of equality.
                // nonbasic_move records the original side.
                basis.col_lower[column] = 0.0;
                basis.col_upper[column] = 0.0;
            }
            self.work_range[column] = basis.col_upper[column] - basis.col_lower[column];

            if (basis.col_status[column] == .basic) continue;

            // ── Infeasibility from original reduced cost ──
            const reduced = basis.reduced_cost[column];
            const infeasibility = if (!lower_finite and !upper_finite)
                @abs(reduced)
            else if (!lower_finite)
                @max(reduced, 0.0)
            else if (!upper_finite)
                @max(-reduced, 0.0)
            else switch (basis.col_status[column]) {
                .at_lower => @max(-reduced, 0.0),
                .at_upper => @max(reduced, 0.0),
                else => 0.0,
            };
            self.dual_infeasibility[column] = infeasibility;

            // ── Pinned HiGHS initialiseNonbasicValueAndMove ──
            if (basis.col_lower[column] == basis.col_upper[column]) {
                // Fixed or collapsed boxed variables do not participate in
                // the Phase-I ratio test or bound flips.
                self.nonbasic_move[column] = 0;
                basis.col_status[column] = .fixed;
                basis.primal[column] = 0.0;
            } else if (!lower_finite and !upper_finite) {
                // Original free move is zero. Once represented by finite
                // Phase-I bounds HiGHS corrects that invalid boxed move to
                // the lower endpoint.
                basis.col_status[column] = .at_lower;
                basis.primal[column] = -1000.0;
                self.nonbasic_move[column] = 1;
            } else if (!lower_finite) {
                // Original upper variable keeps its downward move and is
                // placed at the Phase-I upper endpoint, zero.
                basis.col_status[column] = .at_upper;
                basis.primal[column] = 0.0;
                self.nonbasic_move[column] = -1;
            } else {
                // Original lower variable keeps its upward move and is
                // placed at the Phase-I lower endpoint, zero.
                basis.col_status[column] = .at_lower;
                basis.primal[column] = 0.0;
                self.nonbasic_move[column] = 1;
            }
        }
        for (basis.basic_index, 0..) |column, row| {
            basis.basic_lower[row] = basis.col_lower[column];
            basis.basic_upper[row] = basis.col_upper[column];
        }
    }

    /// Mirror the first `HEkkDual::rebuild -> correctDualInfeasibilities`
    /// after Phase-I bounds have been installed. Every non-fixed Phase-I
    /// variable is boxed, so a dual-infeasible variable is corrected by
    /// flipping to the opposite endpoint, not by shifting its cost. The
    /// resulting nonbasic value (±1, or ±1000 for an originally free
    /// variable) is what encodes the Phase-I dual objective.
    pub fn correctInitialDualInfeasibilities(
        self: *DualPhaseOneWorkspace,
        basis: *basis_module.BasisState,
        original_count: usize,
        dual_tolerance: f64,
    ) usize {
        var flip_count: usize = 0;
        for (0..original_count) |column| {
            if (basis.col_status[column] == .basic) continue;
            const move = self.nonbasic_move[column];
            if (move == 0 or basis.col_lower[column] == basis.col_upper[column]) continue;
            const dual_infeasibility = -@as(f64, @floatFromInt(move)) * basis.reduced_cost[column];
            self.dual_infeasibility[column] = @max(dual_infeasibility, 0.0);
            if (dual_infeasibility < dual_tolerance) continue;
            if (move > 0) {
                basis.col_status[column] = .at_upper;
                basis.primal[column] = basis.col_upper[column];
                self.nonbasic_move[column] = -1;
            } else {
                basis.col_status[column] = .at_lower;
                basis.primal[column] = basis.col_lower[column];
                self.nonbasic_move[column] = 1;
            }
            flip_count += 1;
        }
        return flip_count;
    }

    /// Compute the nonbasic contribution to the current working dual objective.
    pub fn dualObjective(
        self: *const DualPhaseOneWorkspace,
        basis: *const basis_module.BasisState,
        original_count: usize,
    ) f64 {
        _ = self;
        var objective: f64 = 0.0;
        for (0..original_count) |column| {
            if (basis.col_status[column] == .basic) continue;
            objective += basis.primal[column] * basis.reduced_cost[column];
        }
        return objective;
    }

    /// Restore the original model bounds. nonbasic_move (tracked through
    /// the Phase-I pivot path) determines the side for boxed columns that
    /// were collapsed to [0, 0]; for other columns the fresh reduced-cost
    /// sign is used. Basic membership is unchanged.
    pub fn restoreOriginalBounds(self: *DualPhaseOneWorkspace, basis: *basis_module.BasisState, original_count: usize) void {
        @memcpy(basis.col_lower[0..original_count], self.saved_lower[0..original_count]);
        @memcpy(basis.col_upper[0..original_count], self.saved_upper[0..original_count]);
        for (self.work_range[0..original_count], basis.col_lower[0..original_count], basis.col_upper[0..original_count]) |*range_value, lower, upper|
            range_value.* = upper - lower;
        for (0..original_count) |column| {
            if (basis.col_status[column] == .basic) continue;
            const lower = basis.col_lower[column];
            const upper = basis.col_upper[column];
            const lower_finite = std.math.isFinite(lower);
            const upper_finite = std.math.isFinite(upper);
            const move = self.nonbasic_move[column];
            if (lower == upper) {
                basis.col_status[column] = .fixed;
                basis.primal[column] = lower;
            } else if (!lower_finite and !upper_finite) {
                basis.col_status[column] = .free;
                basis.primal[column] = 0.0;
            } else if (!lower_finite) {
                basis.col_status[column] = .at_upper;
                basis.primal[column] = upper;
            } else if (!upper_finite) {
                basis.col_status[column] = .at_lower;
                basis.primal[column] = lower;
            } else if (move != 0) {
                // Boxed column that was collapsed: use the Phase-I-tracked
                // direction. move > 0 → at_lower, move < 0 → at_upper.
                if (move > 0) {
                    basis.col_status[column] = .at_lower;
                    basis.primal[column] = lower;
                } else {
                    basis.col_status[column] = .at_upper;
                    basis.primal[column] = upper;
                }
            } else if (basis.reduced_cost[column] >= 0.0) {
                basis.col_status[column] = .at_lower;
                basis.primal[column] = lower;
            } else {
                basis.col_status[column] = .at_upper;
                basis.primal[column] = upper;
            }
        }
        for (basis.basic_index, 0..) |column, row| {
            basis.basic_lower[row] = basis.col_lower[column];
            basis.basic_upper[row] = basis.col_upper[column];
        }
        self.active = false;
    }

    /// Mirror a committed nonbasic bound flip in the persistent move array.
    pub fn noteBoundFlip(self: *DualPhaseOneWorkspace, column: usize) void {
        if (column >= self.nonbasic_move.len) return;
        self.nonbasic_move[column] = -self.nonbasic_move[column];
    }

    /// Update explicit moves after a committed basis exchange.
    pub fn notePivot(
        self: *DualPhaseOneWorkspace,
        entering_column: usize,
        leaving_column: usize,
        leaving_bound: basis_module.BasisStatus,
    ) void {
        if (entering_column >= self.nonbasic_move.len or leaving_column >= self.nonbasic_move.len) return;
        self.nonbasic_move[entering_column] = 0;
        self.nonbasic_move[leaving_column] = switch (leaving_bound) {
            .at_lower => 1,
            .at_upper => -1,
            .fixed => 0,
            else => 0,
        };
    }

    /// Record leaving violation remaining after the proposed BFRT flip set.
    pub fn recordRemainingViolation(
        self: *DualPhaseOneWorkspace,
        row: usize,
        violation: f64,
        tableau: []const f64,
        flip_columns: []const u32,
        lower: []const f64,
        upper: []const f64,
    ) void {
        if (!self.active or row >= self.remaining_violation.len) return;
        var corrected: f64 = 0.0;
        for (flip_columns) |column_u32| {
            const column: usize = @intCast(column_u32);
            if (column >= tableau.len or column >= lower.len or column >= upper.len) continue;
            corrected += @abs(tableau[column]) * (upper[column] - lower[column]);
        }
        self.remaining_violation[row] = @max(violation - corrected, 0.0);
    }

    /// Return bytes requested by every retained workspace allocation.
    pub fn requestedBytes(self: *const DualPhaseOneWorkspace) usize {
        return std.mem.sliceAsBytes(self.saved_lower).len +
            std.mem.sliceAsBytes(self.saved_upper).len +
            std.mem.sliceAsBytes(self.work_cost).len +
            std.mem.sliceAsBytes(self.work_shift).len +
            std.mem.sliceAsBytes(self.work_range).len +
            std.mem.sliceAsBytes(self.dual_infeasibility).len +
            std.mem.sliceAsBytes(self.perturbation).len +
            std.mem.sliceAsBytes(self.highs_permutation).len +
            std.mem.sliceAsBytes(self.nonbasic_move).len +
            std.mem.sliceAsBytes(self.remaining_violation).len +
            std.mem.sliceAsBytes(self.primal_infeasibility).len +
            std.mem.sliceAsBytes(self.infeasibility_index).len +
            std.mem.sliceAsBytes(self.infeasibility_mark).len +
            std.mem.sliceAsBytes(self.changed_row_index).len +
            std.mem.sliceAsBytes(self.checkpoint_status).len +
            std.mem.sliceAsBytes(self.checkpoint_basic_index).len;
    }
};

/// Compatibility alias while Phase-I-specific helpers are migrated into the
/// unified serial dual controller.
pub const DualWorkState = DualPhaseOneWorkspace;

test "dual work view aliases live basis arrays and owns cost lifecycle" {
    var basis = try basis_module.BasisState.init(std.testing.allocator, 1, 1);
    defer basis.deinit();
    var work = DualWorkState.init(std.testing.allocator);
    defer work.deinit();
    basis.col_status[0..2].* = .{ .at_lower, .at_upper };
    basis.col_lower[0..2].* = .{ 0.0, -2.0 };
    basis.col_upper[0..2].* = .{ 3.0, 4.0 };
    try work.beginCostEpoch(&basis, 2);
    const view = work.view(&basis, 2).?;
    try std.testing.expectEqualSlices(i8, &.{ 1, -1 }, view.move);
    try std.testing.expectEqualSlices(f64, &.{ 3.0, 6.0 }, view.range);
    view.value[0] = 7.0;
    try std.testing.expectEqual(@as(f64, 7.0), basis.primal[0]);
    try std.testing.expect(work.shiftCost(0, 0.25));
    view.dual[0] = 1.5;
    try std.testing.expect(work.shiftBack(view.dual, 0));
    try std.testing.expectEqual(@as(f64, 1.25), basis.reduced_cost[0]);
    try std.testing.expectEqual(@as(f64, 0.0), view.shift[0]);
}

test "dual correction flips boxed and shifts one-sided work costs" {
    var basis = try basis_module.BasisState.init(std.testing.allocator, 1, 4);
    defer basis.deinit();
    var work = DualWorkState.init(std.testing.allocator);
    defer work.deinit();
    basis.col_lower[0..4].* = .{ 0.0, 0.0, -std.math.inf(f64), -std.math.inf(f64) };
    basis.col_upper[0..4].* = .{ 2.0, std.math.inf(f64), 3.0, std.math.inf(f64) };
    basis.col_status[0..4].* = .{ .at_lower, .at_lower, .at_upper, .free };
    basis.primal[0..4].* = .{ 0.0, 0.0, 3.0, 0.0 };
    basis.reduced_cost[0..4].* = .{ -4.0, -5.0, 6.0, 7.0 };
    try work.beginCostEpoch(&basis, 4);
    @memset(work.work_cost[0..4], 0.0);
    work.resetCorrectionRandom(2, 4);
    const correction = work.correctDualInfeasibilities(&basis, 4, 1e-7, false);
    try std.testing.expectEqual(@as(usize, 1), correction.flip_count);
    try std.testing.expectEqual(@as(usize, 2), correction.shift_count);
    try std.testing.expectEqual(@as(usize, 1), correction.free_infeasibility_count);
    try std.testing.expectEqual(basis_module.BasisStatus.at_upper, basis.col_status[0]);
    try std.testing.expectEqual(@as(f64, 2.0), basis.primal[0]);
    try std.testing.expect(basis.reduced_cost[1] > 1e-7 and basis.reduced_cost[1] < 2e-7);
    try std.testing.expect(basis.reduced_cost[2] < -1e-7 and basis.reduced_cost[2] > -2e-7);
    try std.testing.expectEqual(basis.reduced_cost[1] + 5.0, work.work_cost[1]);
    try std.testing.expectEqual(basis.reduced_cost[2] - 6.0, work.work_cost[2]);
}

test "dual RHS sparse infeasibility list retains insertion order" {
    var basis = try basis_module.BasisState.init(std.testing.allocator, 5, 0);
    defer basis.deinit();
    var work = DualWorkState.init(std.testing.allocator);
    defer work.deinit();
    try work.ensureCapacity(5);
    @memset(basis.basic_value, 0.0);
    @memset(basis.basic_lower, 0.0);
    @memset(basis.basic_upper, 10.0);
    @memset(basis.row_edge_weight, 1.0);
    basis.basic_value[2] = -2.0;

    try std.testing.expect(work.createPrimalInfeasibilityList(&basis, 5, 1e-7));
    try std.testing.expect(!work.infeasibility_dense);
    try std.testing.expectEqual(@as(usize, 1), work.infeasibility_count);
    try std.testing.expectEqual(@as(u32, 2), work.infeasibility_index[0]);

    // A feasible marked row remains in place, while a newly infeasible row is
    // appended according to the FTRAN index order.
    basis.basic_value[2] = 0.0;
    basis.basic_value[4] = -3.0;
    work.updatePrimalInfeasibilityList(&basis, &.{ 2, 4 }, 1e-7);
    try std.testing.expectEqual(@as(usize, 2), work.infeasibility_count);
    try std.testing.expectEqualSlices(u32, &.{ 2, 4 }, work.infeasibility_index[0..2]);
    try std.testing.expectEqual(@as(?u32, 4), work.choosePrimalInfeasibleRow(&basis, 5));
}

test "dual infeasibility statistics use explicit move and ignore fixed" {
    var basis = try basis_module.BasisState.init(std.testing.allocator, 1, 4);
    defer basis.deinit();
    var work = DualWorkState.init(std.testing.allocator);
    defer work.deinit();
    basis.col_lower[0..4].* = .{ 0.0, -2.0, 3.0, -std.math.inf(f64) };
    basis.col_upper[0..4].* = .{ 5.0, 7.0, 3.0, std.math.inf(f64) };
    basis.col_status[0..4].* = .{ .at_lower, .at_upper, .fixed, .free };
    basis.reduced_cost[0..4].* = .{ -2.0, 4.0, -100.0, -3.0 };
    try work.beginCostEpoch(&basis, 4);
    const stats = work.dualInfeasibilityStats(&basis, 4, 1e-7);
    try std.testing.expectEqual(@as(usize, 3), stats.count);
    try std.testing.expectEqual(@as(f64, 4.0), stats.maximum);
    try std.testing.expectEqual(@as(f64, 9.0), stats.sum);
}

test "dual rebuild state records primal infeasibilities and objective" {
    var basis = try basis_module.BasisState.init(std.testing.allocator, 2, 2);
    defer basis.deinit();
    var work = DualWorkState.init(std.testing.allocator);
    defer work.deinit();
    try work.beginCostEpoch(&basis, 4);
    basis.basic_value[0..2].* = .{ -3.0, 7.0 };
    basis.basic_lower[0..2].* = .{ -1.0, 0.0 };
    basis.basic_upper[0..2].* = .{ 4.0, 5.0 };
    basis.col_status[0..4].* = .{ .at_lower, .at_upper, .basic, .basic };
    basis.primal[0..4].* = .{ 2.0, 3.0, -3.0, 7.0 };
    basis.reduced_cost[0..4].* = .{ 5.0, -2.0, 0.0, 0.0 };
    try std.testing.expect(work.recordRebuildState(&basis, 4, 1e-7));
    try std.testing.expectEqual(@as(usize, 2), work.rebuild_primal_infeasibility_count);
    try std.testing.expectEqual(@as(f64, 2.0), work.rebuild_primal_infeasibility_max);
    try std.testing.expectEqual(@as(f64, 4.0), work.rebuild_primal_infeasibility_sum);
    try std.testing.expectEqualSlices(f64, &.{ 2.0, 2.0 }, work.remaining_violation[0..2]);
    try std.testing.expectEqual(@as(f64, 4.0), work.rebuild_dual_objective);
    try std.testing.expectEqual(@as(u64, 1), work.rebuild_epoch);
}

test "dual Phase-I initial correction flips infeasibilities to objective endpoints" {
    var basis = try basis_module.BasisState.init(std.testing.allocator, 1, 4);
    defer basis.deinit();
    var workspace = DualPhaseOneWorkspace.init(std.testing.allocator);
    defer workspace.deinit();
    // col0: LOWER  [0, inf]   rc=-2 at_lower   → infeasible → value= 1
    // col1: UPPER  [-inf, 3]  rc= 3 at_upper   → infeasible → value=-1
    // col2: FREE   [-inf, inf] rc=-4 free       → infeasible → value=-1
    // col3: BOXED  [-2, 2]    rc=-5 at_lower   → collapsed [0,0], move=1
    basis.col_lower[0..4].* = .{ 0.0, -std.math.inf(f64), -std.math.inf(f64), -2.0 };
    basis.col_upper[0..4].* = .{ std.math.inf(f64), 3.0, std.math.inf(f64), 2.0 };
    basis.col_status[0..4].* = .{ .at_lower, .at_upper, .free, .at_lower };
    basis.reduced_cost[0..4].* = .{ -2.0, 3.0, -4.0, -5.0 };
    try workspace.begin(&basis, 4);
    workspace.installWorkingBounds(&basis, 4);
    // Working bounds: LOWER [0,1]  UPPER [-1,0]  FREE [-1000,1000]  BOXED [0,0]
    try std.testing.expectEqualSlices(f64, &.{ 0.0, -1.0, -1000.0, 0.0 }, basis.col_lower[0..4]);
    try std.testing.expectEqualSlices(f64, &.{ 1.0, 0.0, 1000.0, 0.0 }, basis.col_upper[0..4]);
    // initialiseNonbasicValueAndMove first uses the original move and zero
    // endpoint (free is corrected to the finite Phase-I lower endpoint).
    try std.testing.expectEqualSlices(f64, &.{ 0.0, 0.0, -1000.0, 0.0 }, basis.primal[0..4]);
    try std.testing.expectEqualSlices(i8, &.{ 1, -1, 1, 0 }, workspace.nonbasic_move[0..4]);
    // The first rebuild flips each Phase-I boxed dual infeasibility. This is
    // where ±1/±1000 values and the negative Phase-I objective arise.
    try std.testing.expectEqual(@as(usize, 3), workspace.correctInitialDualInfeasibilities(&basis, 4, 1e-7));
    try std.testing.expectEqualSlices(f64, &.{ 1.0, -1.0, 1000.0, 0.0 }, basis.primal[0..4]);
    try std.testing.expectEqualSlices(i8, &.{ -1, 1, -1, 0 }, workspace.nonbasic_move[0..4]);
    try std.testing.expectApproxEqAbs(@as(f64, -4005.0), workspace.dualObjective(&basis, 4), 1e-12);
    workspace.noteBoundFlip(0);
    try std.testing.expectEqual(@as(i8, 1), workspace.nonbasic_move[0]);
    workspace.notePivot(1, 3, .at_upper);
    try std.testing.expectEqual(@as(i8, 0), workspace.nonbasic_move[1]);
    try std.testing.expectEqual(@as(i8, -1), workspace.nonbasic_move[3]);
    workspace.restoreOriginalBounds(&basis, 4);
    try std.testing.expectEqualSlices(f64, &.{ 0.0, -std.math.inf(f64), -std.math.inf(f64), -2.0 }, basis.col_lower[0..4]);
}

test "dual Phase-I records violation left after bound flips" {
    var workspace = DualPhaseOneWorkspace.init(std.testing.allocator);
    defer workspace.deinit();
    var basis = try basis_module.BasisState.init(std.testing.allocator, 1, 1);
    defer basis.deinit();
    try workspace.begin(&basis, 2);

    workspace.recordRemainingViolation(
        0,
        10,
        &.{ 2, -1 },
        &.{0},
        &.{ 0, 0 },
        &.{ 3, 4 },
    );
    try std.testing.expectApproxEqAbs(@as(f64, 4), workspace.remaining_violation[0], 1e-12);
}

test "dual Phase-I HiGHS bounds: boxed collapsed, basic bounds follow" {
    var basis = try basis_module.BasisState.init(std.testing.allocator, 1, 0);
    defer basis.deinit();
    var workspace = DualPhaseOneWorkspace.init(std.testing.allocator);
    defer workspace.deinit();
    basis.basic_index[0] = 0;
    basis.col_status[0] = .basic;
    basis.col_lower[0] = 0;
    basis.col_upper[0] = 1;
    basis.basic_value[0] = 12; // irrelevant: bounds are hardcoded

    try workspace.begin(&basis, 1);
    workspace.installWorkingBounds(&basis, 1);

    // BOXED column → collapsed [0, 0]; basic bounds follow column bounds
    try std.testing.expectApproxEqAbs(@as(f64, 1000), workspace.working_radius, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0), basis.col_lower[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0), basis.col_upper[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0), basis.basic_lower[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0), basis.basic_upper[0], 1e-12);
}
