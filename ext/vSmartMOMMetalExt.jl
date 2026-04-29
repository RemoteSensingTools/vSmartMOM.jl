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

Architectures.devi(::vSmartMOM.Architectures.MetalGPU) = Metal.MetalBackend()
Architectures.array_type(::vSmartMOM.Architectures.MetalGPU) = Metal.MtlArray
Architectures.architecture(::Metal.MtlArray) = vSmartMOM.Architectures.MetalGPU()
Architectures.architecture(::Type{<:Metal.MtlArray}) = vSmartMOM.Architectures.MetalGPU()

"Batched matrix multiply for 3D Metal arrays using a portable KA kernel."
function vSmartMOM.CoreRT.batched_mul(A::_MetalArray3{FT}, B::_MetalArray3{FT}) where {FT}
    vSmartMOM.CoreRT.ka_batched_mul(A, B, Metal.MetalBackend())
end

"Batched multiply for 3D Metal array views."
@inline _as_mtlarray3(A::SubArray{FT,3,<:Metal.MtlArray}) where {FT} = copy(A)

function vSmartMOM.CoreRT.batched_mul(A::SubArray{FT,3,<:Metal.MtlArray}, B::_MetalArray3{FT}) where {FT}
    vSmartMOM.CoreRT.batched_mul(_as_mtlarray3(A), B)
end

function vSmartMOM.CoreRT.batched_mul(A::_MetalArray3{FT}, B::SubArray{FT,3,<:Metal.MtlArray}) where {FT}
    vSmartMOM.CoreRT.batched_mul(A, _as_mtlarray3(B))
end

function vSmartMOM.CoreRT.batched_mul(A::SubArray{FT,3,<:Metal.MtlArray},
                                      B::SubArray{FT,3,<:Metal.MtlArray}) where {FT}
    vSmartMOM.CoreRT.batched_mul(_as_mtlarray3(A), _as_mtlarray3(B))
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
    vSmartMOM.CoreRT.ka_batch_inv_lu!(X, A, Metal.MetalBackend())
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

"""
    CoreRT.make_gpu_rt_workspace(FT, NquadN, nSpec, ::MetalGPU)

Create an RT workspace backed by Metal `MtlArray`s.
"""
function vSmartMOM.CoreRT.make_gpu_rt_workspace(FT::Type,
                                                NquadN::Int,
                                                nSpec::Int,
                                                ::vSmartMOM.Architectures.MetalGPU)
    dims3 = (NquadN, NquadN, nSpec)
    dims_J = (NquadN, 1, nSpec)

    vSmartMOM.CoreRT.RTWorkspace(
        Metal.MtlArray(zeros(FT, dims3)),
        Metal.MtlArray(zeros(FT, dims3)),
        Metal.MtlArray(zeros(FT, dims_J)),
        Metal.MtlArray(zeros(FT, dims_J)),
        zeros(Cint, NquadN, nSpec),
        zeros(Cint, nSpec),
        Metal.MtlArray(zeros(FT, dims3)),
        Metal.MtlArray(zeros(FT, dims3)),
        Metal.MtlArray(zeros(FT, dims3)),
        Metal.MtlArray(zeros(FT, dims3)),
    )
end

function __init__()
    if Sys.isapple() && Metal.functional()
        try
            test_arr = Metal.MtlArray([1.0f0])
            Metal.synchronize()
            test_arr = nothing

            Architectures._has_metal[] = true
            Architectures._sync_gpu[] = Metal.synchronize
        catch e
            @warn "vSmartMOM Metal initialization failed, falling back to CPU" exception=e
            Architectures._has_metal[] = false
        end
    else
        Architectures._has_metal[] = false
    end
end

end # module vSmartMOMMetalExt
