"""
    vSmartMOMMetalExt

Package extension that loads Apple Metal support when Metal.jl is available.
Metal is an optional dependency; vSmartMOM continues to load without it on
non-Mac systems.
"""
module vSmartMOMMetalExt

using vSmartMOM
using vSmartMOM.Architectures
using vSmartMOM.CoreRT
using vSmartMOM.Scattering
using Metal

const _MetalArray3{FT} = Metal.MtlArray{FT,3}
const _MetalArrayOrView3{FT} = Union{Metal.MtlArray{FT,3}, SubArray{FT,3,<:Metal.MtlArray}}
const METAL_BATCH_INV_LOCALMEM_LIMIT_BYTES = 32 * 1024

Architectures.devi(::vSmartMOM.Architectures.MetalGPU) = Metal.MetalBackend()
Architectures.array_type(::vSmartMOM.Architectures.MetalGPU) = Metal.MtlArray
Architectures.architecture(::Metal.MtlArray) = vSmartMOM.Architectures.MetalGPU()
Architectures.architecture(::Type{<:Metal.MtlArray}) = vSmartMOM.Architectures.MetalGPU()

# ============================================================================
# Metal Mie route
# ============================================================================
# Architecture → KernelAbstractions backend mapping for the GPU Mie pipeline
# (mirrors `Architectures.ka_backend(::GPU) = CUDA.CUDABackend()` in
# vSmartMOMCUDAExt.jl). Needed by BOTH the log-normal `MieModel` GPU router
# (`_dispatch_aerosol_optics` in phase_function_autodiff.jl) and the
# caller-node GPU functions (`compute_aerosol_optical_properties_nodes_gpu` /
# `compute_aerosol_extinction_nodes_gpu`).
Architectures.ka_backend(::vSmartMOM.Architectures.MetalGPU) = Metal.MetalBackend()

# Real-dispatch refinement of Scattering._is_metal_backend (defaults to
# `false` on `Any` in gpu_precision.jl) -- the same weak-dependency pattern as
# every other Metal hook in this file: a concrete method on the actual
# `Metal.MetalBackend` type, not name/string matching. Used by the caller-node
# GPU preamble (`Scattering._prepare_node_mie_gpu`) to refuse (with a clear
# `ArgumentError`) any policy/element-type combination that would need a
# Float64 device array on Metal.
Scattering._is_metal_backend(::Metal.MetalBackend) = true

# Metal has NO Float64 hardware support at all -- there is no "NativeFloat64
# on Metal", regardless of caller. Fail fast with a clear, actionable error
# instead of letting a Float64 auto-select attempt reach a device-array
# allocation deep inside KernelAbstractions/Metal.jl.
Architectures.default_mie_precision_policy(::vSmartMOM.Architectures.MetalGPU,
                                           ::Type{Float64}) =
    throw(ArgumentError(
        "MetalGPU() has no Float64 hardware support: there is no NativeFloat64 " *
        "policy on Metal. Use Float32 node/model inputs (which auto-select " *
        "Scattering.NativeFloat32() on Metal)."))

# Float32 on Metal auto-selects NativeFloat32 -- the pure-Float32,
# never-Float64-widened policy (Dₙ recursion AND size-distribution reduction
# AND Greek projection all in native Float32 + Neumaier compensation). This is
# the ONLY policy registered for Metal this round: `DSEmulated`'s historical
# reduction discipline widens to Float64 for a Float32 kernel (see
# `Scattering.MiePrecisionPolicy`'s docstring), which is illegal on Metal
# device arrays; making DSEmulated-on-Metal correct would need a
# Float32-COMPENSATED (not widened) reduction variant, which is NOT
# implemented this round (tracked as a follow-up). The caller-node GPU
# preamble (`Scattering._prepare_node_mie_gpu`) additionally guards this at
# the call site -- an explicit `precision_policy = DSEmulated()` passed
# directly to the `_gpu` functions on a Metal backend throws a clear
# `ArgumentError` there too, rather than crashing inside `KernelAbstractions.allocate`.
Architectures.default_mie_precision_policy(::vSmartMOM.Architectures.MetalGPU,
                                           ::Type{Float32}) =
    Scattering.NativeFloat32()
# Any other FT (e.g. a transient Dual that should never reach this router, or
# Float16): fall back to NativeFloat32, the only policy Metal supports at all
# (mirrors the CUDA extension's "any other FT" catch-all, but NativeFloat64
# would be nonsensical here since Metal has no Float64 hardware whatsoever).
Architectures.default_mie_precision_policy(::vSmartMOM.Architectures.MetalGPU,
                                           ::Type) =
    Scattering.NativeFloat32()

# Metal now has a GPU Mie pipeline (NativeFloat32 only -- see above), unlike
# before this trait was refined: `has_gpu_mie` is a single architecture-level
# capability trait shared by the log-normal `MieModel` GPU router AND the
# caller-node GPU functions, so this opts BOTH in. This is safe for BOTH:
# `compute_aerosol_optical_properties_gpu` (compute_NAI2_gpu.jl, the
# log-normal path) now dispatches Kernel 1 via the SAME shared
# `Scattering._mie_kernel1(policy, backend)` helper the caller-node path uses
# (NativeFloat64/NativeFloat32 -> native kernel, DSEmulated -> double-single
# kernel) and carries its own Metal-only Float64 device-array guard
# (`Scattering._is_metal_backend`) -- its host-side size-distribution
# reduction is plain `Array`-copied-back CPU arithmetic regardless of policy,
# so (unlike the caller-node path's device-resident reduction) it was never at
# risk of an illegal Float64 DEVICE array in the first place. The caller-node
# GPU path is fully `NativeFloat32`-aware end to end (see
# `Scattering._prepare_node_mie_gpu`).
Architectures.has_gpu_mie(::vSmartMOM.Architectures.MetalGPU) = true

"Return the backing Metal array used to allocate outputs for arrays or views."
@inline _metal_storage(A::Metal.MtlArray) = A
@inline _metal_storage(A::SubArray{FT,3,<:Metal.MtlArray}) where {FT} = parent(A)

"Batched matrix multiply for 3D Metal arrays or views using a portable KA kernel."
function vSmartMOM.CoreRT.batched_mul(A::_MetalArrayOrView3{FT},
                                      B::_MetalArrayOrView3{FT}) where {FT}
    C = similar(_metal_storage(A), FT, (size(A, 1), size(B, 2), size(A, 3)))
    vSmartMOM.CoreRT.ka_batched_mul!(C, A, B, Metal.MetalBackend())
end

"Given 3D Metal arrays A and B, fill in X[:,:,k] = A[:,:,k] \\ B[:,:,k]."
function vSmartMOM.CoreRT.batch_solve!(X::_MetalArray3{FT},
                                       A::_MetalArray3{FT},
                                       B::_MetalArray3{FT}) where {FT}
    temp = similar(A)
    vSmartMOM.CoreRT.batch_inv!(temp, A)
    X .= vSmartMOM.CoreRT.batched_mul(temp, B)
    Metal.synchronize()
    return X
end

"Given 3D Metal array A, fill in X[:,:,k] = inv(A[:,:,k])."
function vSmartMOM.CoreRT.batch_inv!(X::_MetalArray3{FT}, A::_MetalArray3{FT}) where {FT}
    vSmartMOM.CoreRT.ka_batch_inv_lu!(
        X, A, Metal.MetalBackend();
        max_localmem_bytes=METAL_BATCH_INV_LOCALMEM_LIMIT_BYTES,
    )
    Metal.synchronize()
    return X
end

function vSmartMOM.CoreRT.batch_inv!(X::_MetalArray3{FT},
                                     A::_MetalArray3{FT},
                                     ::Nothing,
                                     ::Nothing) where {FT}
    vSmartMOM.CoreRT.batch_inv!(X, A)
end

function vSmartMOM.CoreRT.batch_inv!(X::_MetalArray3{FT},
                                     A::_MetalArray3{FT},
                                     ws::vSmartMOM.CoreRT.RTWorkspace) where {FT}
    vSmartMOM.CoreRT.batch_inv!(X, A)
end

function __init__()
    if !Sys.isapple()
        Architectures._has_metal[] = false
        return
    end

    metal_functional = try
        Metal.functional()
    catch e
        @warn "vSmartMOM Metal availability check failed, falling back to CPU" exception=e
        false
    end

    if !metal_functional
        Architectures._has_metal[] = false
        @info "Metal.jl is loaded but no functional Metal device is available; default_architecture() will not select MetalGPU()."
        return
    end

    try
        test_arr = Metal.MtlArray([1.0f0])
        Metal.synchronize()
        test_arr = nothing

        Architectures._has_metal[] = true
        if !Architectures.has_cuda()
            Architectures._sync_gpu[] = Metal.synchronize
        end
    catch e
        @warn "vSmartMOM Metal initialization failed, falling back to CPU" exception=e
        Architectures._has_metal[] = false
    end
end

end # module vSmartMOMMetalExt
