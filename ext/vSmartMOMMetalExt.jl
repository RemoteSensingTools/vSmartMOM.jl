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
using Metal

const _MetalArray3{FT} = Metal.MtlArray{FT,3}
const _MetalArrayOrView3{FT} = Union{Metal.MtlArray{FT,3}, SubArray{FT,3,<:Metal.MtlArray}}
const METAL_BATCH_INV_LOCALMEM_LIMIT_BYTES = 32 * 1024

Architectures.devi(::vSmartMOM.Architectures.MetalGPU) = Metal.MetalBackend()
Architectures.array_type(::vSmartMOM.Architectures.MetalGPU) = Metal.MtlArray
Architectures.architecture(::Metal.MtlArray) = vSmartMOM.Architectures.MetalGPU()
Architectures.architecture(::Type{<:Metal.MtlArray}) = vSmartMOM.Architectures.MetalGPU()

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
