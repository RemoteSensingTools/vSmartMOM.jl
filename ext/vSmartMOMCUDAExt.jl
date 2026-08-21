"""
    vSmartMOMCUDAExt

Package extension that loads CUDA support when CUDA.jl is available.
This allows vSmartMOM to work without CUDA on systems without compatible GPUs.

Requires Julia 1.9+ with package extensions feature.
"""
module vSmartMOMCUDAExt

using vSmartMOM
using vSmartMOM.Architectures
using vSmartMOM.CoreRT
using vSmartMOM.Scattering
using CUDA
using KernelAbstractions

# Re-export CUDAKernels for code that needs CUDADevice
import CUDA.CUDAKernels

# Extend GPU device to use CUDA backend
Architectures.devi(::vSmartMOM.Architectures.GPU) = CUDA.CUDABackend(; always_inline=true)

# Extend array_type to return CuArray for GPU
Architectures.array_type(::vSmartMOM.Architectures.GPU) = CuArray

# Extend architecture detection for CuArrays
Architectures.architecture(::CuArray) = vSmartMOM.Architectures.GPU()

# Memory-aware interaction staging (v0.7+). Non-staged SI_11 interaction is ~6×
# faster than CPU-staged. Choose it only when the FULL GPU footprint it triggers
# (all allocated AFTER this call) will fit free VRAM; else fall back to staging
# (memory-frugal, PCIe-bound). That footprint is:
#   InteractionWorkspace 4-D : gpu_ie_mat_A/B (2·mat4d) + gpu_ie_src (1·src4d)
#   InteractionWorkspace 3-D : tmp_inv, tmpR⁻⁺, tmpR⁺⁻, tmpT⁻⁻, tmpT⁺⁺ (5·mat3d)
#                              + tmpJ₀⁻, tmpJ₀⁺ (2·src3d, negligible)
#   per-call non-staged temps: tmpieT⁻⁻ + tmpieT⁺⁺ (2·mat4d) + tmpieJ₀⁺ (1·src4d)
# The 3-D terms (mat3d = mat4d/nRaman) matter at small nRaman (codex P2). +25% margin
# so we err toward staging when uncertain — an OOM is worse than a slow run.
function CoreRT._use_staged_interaction(::vSmartMOM.Architectures.GPU, composite_layer)
    N, _, nSpec, nRaman = size(composite_layer.ieR⁻⁺)
    FT = eltype(composite_layer.ieR⁻⁺)
    mat4d = N * N * nSpec * nRaman * sizeof(FT)
    src4d = N * nSpec * nRaman * sizeof(FT)
    mat3d = N * N * nSpec * sizeof(FT)
    src3d = N * nSpec * sizeof(FT)
    footprint = 4 * mat4d + 2 * src4d + 5 * mat3d + 2 * src3d
    return CUDA.available_memory() < 1.25 * footprint
end

# Architecture → KernelAbstractions backend mapping for the GPU Mie pipeline.
# `always_inline=false` keeps the (large) Mie kernels from blowing up
# compilation; the batched RT path uses its own backend via `devi`.
Architectures.ka_backend(::vSmartMOM.Architectures.GPU) = CUDA.CUDABackend()

# Backend policy: on CUDA the ScatteringInterface_11 interaction updates are
# plain batched GEMMs at N≈36 — cuBLAS's best case — and the fused KA kernels
# measured SLOWER (3.27 s vs 2.89 s, production satellite solve, A100). Keep
# the `_bmm!` ladder under `:auto` here; `_FUSED_INTERACTION_MODE[] = :on`
# still forces the fused path (needed for cuBLAS-free CUDA-graph capture).
# Non-CUDA GPU backends keep CoreRT's generic `true` default.
CoreRT._fused_interaction_default(::CUDA.CUDABackend) = false

# Default GPU Mie precision policy, FT-aware:
#   Float64 → NativeFloat64 (hardware FP64 Dₙ recursion) — right for datacenter
#             GPUs (A100/V100) and bit-for-bit with the CPU reference.
#   Float32 → DSEmulated (Float32 double-single Dₙ recursion). A Float32 model
#             CANNOT use NativeFloat64: compute_aerosol_optical_properties_gpu
#             asserts FT === Float64 for that policy. Existing GPU YAMLs built
#             with float_type: Float32 would otherwise fail during aerosol-optics
#             construction. DSEmulated is the Float32-native path; its host-side
#             bulk reduction is already Float64-widened + Neumaier-compensated, so
#             accuracy sits at the Float32 representational floor.
Architectures.default_mie_precision_policy(::vSmartMOM.Architectures.GPU,
                                           ::Type{Float64}) =
    Scattering.NativeFloat64()
Architectures.default_mie_precision_policy(::vSmartMOM.Architectures.GPU,
                                           ::Type{Float32}) =
    Scattering.DSEmulated()
# Any other FT (e.g. a transient Dual that should never reach the GPU router):
# fall back to NativeFloat64 to preserve the historical default.
Architectures.default_mie_precision_policy(::vSmartMOM.Architectures.GPU,
                                           ::Type) =
    Scattering.NativeFloat64()

# CUDA GPU has a full KernelAbstractions Mie pipeline → opt into the GPU Mie
# path. (Metal defines no Mie path, so MetalGPU keeps the base `false`.)
Architectures.has_gpu_mie(::vSmartMOM.Architectures.GPU) = true

# No need to override synchronize_if_gpu - it checks has_cuda() and calls CUDA.synchronize()
# which will be available from this extension

# Include GPU-specific batched operations
include("gpu_batched_cuda.jl")

# Module initialization - called when extension is loaded
function __init__()
    # CUDA.functional() checks:
    # 1. CUDA libraries are available
    # 2. A compatible GPU is present
    # 3. The CUDA toolkit version matches Julia's requirements
    # 4. The driver version is compatible

    if CUDA.functional()
        try
            # Make CUBLAS available to CoreRT for batched GPU operations (make_added_layer, doubling)
            CoreRT.CUBLAS_ref[] = CUDA.CUBLAS
            # Additional safety check: try to actually use CUDA
            # This catches edge cases where CUDA.functional() returns true
            # but operations fail due to version mismatches
            test_arr = CUDA.CuArray([1.0f0])
            CUDA.synchronize()

            # If we got here, CUDA really works!
            Architectures._has_cuda[] = true
            Architectures._sync_gpu[] = CUDA.synchronize

            CUDA.allowscalar(false)

        catch e
            @warn "vSmartMOM GPU initialization failed, falling back to CPU" exception=e
            Architectures._has_cuda[] = false
        end
    else
        Architectures._has_cuda[] = false
    end
end

end # module vSmartMOMCUDAExt
