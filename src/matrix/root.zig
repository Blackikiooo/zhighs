//! Sparse and dense matrix data structures.
//!
//! Construction uses mutable SoA builders; solver kernels consume canonical,
//! owning sparse vectors and CSC matrices.
//!
//! API stability policy:
//! - A declaration documented as `Stable API` is intended for production use.
//! - Every other public declaration is currently an `Experimental API`. Its
//!   signature, ownership model, validation contract, or implementation may be
//!   corrected as the LP/presolve/simplex layers expose new requirements.
//!
//! This policy prevents today's matrix prototypes from accidentally becoming a
//! permanent compatibility promise. Experimental APIs remain fully testable and
//! usable inside this repository, but callers should isolate them behind their
//! own adapter.

const std = @import("std");
const sparse_vector = @import("sparse_vector.zig");
const sparse_vector_builder = @import("sparse_vector_builder.zig");
const csc = @import("csc.zig");
const builder = @import("builder.zig");
const csr_view = @import("csr_view.zig");
const ops = @import("ops.zig");
const transpose_module = @import("transpose.zig");
const slice_module = @import("slice.zig");
const scaling = @import("scaling.zig");
const permutation = @import("permutation.zig");
const edit = @import("edit.zig");
const dynamic_rows = @import("dynamic_rows.zig");
const store = @import("store.zig");
const sparse_sum = @import("sparse_sum.zig");
const memory = @import("memory.zig");

pub const SparseVectorError = sparse_vector.SparseVectorError;
pub const SparseVectorView = sparse_vector.SparseVectorView;
pub const SparseVector = sparse_vector.SparseVector;
pub const SparseVectorBuilder = sparse_vector_builder.SparseVectorBuilder;
pub const MatrixError = csc.MatrixError;
pub const CscMatrix = csc.CscMatrix;
/// Stable API: read-only borrowed view over canonical CSC storage.
pub const CscView = csc.CscView;
pub const MatrixBuilder = builder.MatrixBuilder;
/// Stable API: reusable caller-owned storage for canonical CSC construction.
pub const CscBuildBuffers = builder.CscBuildBuffers;
/// Experimental API: may change as sorted-input construction is integrated.
pub const freezeFromSortedArraysAssumeValid = builder.freezeFromSortedArraysAssumeValid;
/// Experimental API: may change or be replaced by the reusable-buffer path.
pub const freezeFromCanonicalArraysAssumeValid = builder.freezeFromCanonicalArraysAssumeValid;
/// Stable API: allocation-free canonical CSC construction into caller buffers.
pub const freezeCanonicalIntoAssumeValid = builder.freezeCanonicalIntoAssumeValid;
pub const CsrView = csr_view.CsrView;
pub const CsrCache = csr_view.CsrCache;
pub const CsrBuffers = csr_view.CsrBuffers;
pub const fillCsrFromCsc = csr_view.fillFromCsc;
pub const fillCsrFromCscAssumeValid = csr_view.fillFromCscAssumeValid;
pub const maxAbs = ops.maxAbs;
pub const AbsoluteRange = ops.AbsoluteRange;
pub const ValueAssessment = ops.ValueAssessment;
pub const eql = ops.eql;
pub const absoluteRange = ops.absoluteRange;
pub const assessValues = ops.assessValues;
pub const hasLargeValue = ops.hasLargeValue;
pub const columnOneNorms = ops.columnOneNorms;
pub const columnOneNormsAssumeValid = ops.columnOneNormsAssumeValid;
pub const rowOneNorms = ops.rowOneNorms;
pub const rowOneNormsAssumeValid = ops.rowOneNormsAssumeValid;
pub const oneNorm = ops.oneNorm;
pub const infinityNorm = ops.infinityNorm;
pub const frobeniusNorm = ops.frobeniusNorm;
pub const addProduct = ops.addProduct;
pub const addProductAssumeValid = ops.addProductAssumeValid;
pub const addProductSkippingZeros = ops.addProductSkippingZeros;
pub const addProductSkippingZerosAssumeValid = ops.addProductSkippingZerosAssumeValid;
pub const addTransposeProduct = ops.addTransposeProduct;
pub const addTransposeProductAssumeValid = ops.addTransposeProductAssumeValid;
pub const multiplyHighPrecision = ops.multiplyHighPrecision;
pub const multiplyHighPrecisionAssumeValid = ops.multiplyHighPrecisionAssumeValid;
pub const multiplyHighPrecisionFastAssumeValid = ops.multiplyHighPrecisionFastAssumeValid;
pub const multiplyCompensatedAssumeValid = ops.multiplyCompensatedAssumeValid;
pub const transposeMultiplyHighPrecision = ops.transposeMultiplyHighPrecision;
pub const transposeMultiplyHighPrecisionAssumeValid = ops.transposeMultiplyHighPrecisionAssumeValid;
pub const transposeMultiplyHighPrecisionFastAssumeValid = ops.transposeMultiplyHighPrecisionFastAssumeValid;
pub const transposeMultiplyCompensatedAssumeValid = ops.transposeMultiplyCompensatedAssumeValid;
pub const transpose = transpose_module.transpose;
pub const transposeAssumeValid = transpose_module.transposeAssumeValid;
/// Experimental API: owning lean layout and allocation policy may change.
pub const transposeLeanAssumeValid = transpose_module.transposeLeanAssumeValid;
/// Experimental API: compact-offset representation may change after solver profiling.
pub const transposeLeanAssumeValidCompact = transpose_module.transposeLeanAssumeValidCompact;
pub const transposeInto = transpose_module.transposeInto;
pub const transposeIntoAssumeValid = transpose_module.transposeIntoAssumeValid;
/// Experimental API: compact scratch/output contract may change after integration.
pub const transposeIntoAssumeValidCompact = transpose_module.transposeIntoAssumeValidCompact;
pub const TransposeBuffers = transpose_module.TransposeBuffers;
pub const extractColumns = slice_module.extractColumns;
pub const extractColumnsAssumeValid = slice_module.extractColumnsAssumeValid;
pub const extractRows = slice_module.extractRows;
pub const extractRowsAssumeValid = slice_module.extractRowsAssumeValid;
pub const extractColumnRange = slice_module.extractColumnRange;
pub const ScalingView = scaling.ScalingView;
pub const applyScaling = scaling.apply;
pub const applyScalingAssumeValid = scaling.applyAssumeValid;
pub const removeScaling = scaling.remove;
pub const removeScalingAssumeValid = scaling.removeAssumeValid;
pub const computeMaxEquilibration = scaling.computeMaxEquilibration;
pub const scaleColumn = scaling.scaleColumn;
pub const scaleRow = scaling.scaleRow;
pub const applyColumnFactors = scaling.applyColumnFactors;
pub const applyRowFactors = scaling.applyRowFactors;
pub const computePowerOfTwoColumnFactors = scaling.computePowerOfTwoColumnFactors;
pub const computePowerOfTwoRowFactors = scaling.computePowerOfTwoRowFactors;
pub const permute = permutation.permute;
pub const permuteAssumeValid = permutation.permuteAssumeValid;
pub const appendColumns = edit.appendColumns;
pub const appendRows = edit.appendRows;
pub const deleteColumns = edit.deleteColumns;
pub const deleteRows = edit.deleteRows;
pub const DynamicRowMatrix = dynamic_rows.DynamicRowMatrix;
pub const MatrixStoreError = store.MatrixStoreError;
pub const MatrixStore = store.MatrixStore;
pub const SparseAccumulator = sparse_sum.SparseAccumulator;
pub const clearF64 = memory.clearF64;

test {
    std.testing.refAllDecls(@This());
}
