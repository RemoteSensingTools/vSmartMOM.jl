#=

This file contains CPU implementations of batched linear algebra operations.

These are fallback methods for AbstractArray types (CPU arrays).
GPU-accelerated versions using CuArray are defined in vSmartMOMCUDAExt 
when CUDA is loaded.

=#

# Default synchronization (no-op for CPU, overridden in CUDAExt for GPU)
@inline synchronize() = nothing

"""
    batched_pointer_cache(A)

Return backend-specific batched pointer metadata for `A`, or `nothing` when
the backend does not need pointer arrays. CUDA overrides this for `CuArray`
so the existing CUBLAS batched inverse path can use `unsafe_strided_batch`.
Portable KernelAbstractions backends, including Metal, use `nothing`.
"""
@inline batched_pointer_cache(::AbstractArray) = nothing

"Given 3D Julia Arrays A and B, fill in X[:,:,k] = A[:,:,k] \\ B[:,:,k]"
function batch_solve!(X::AbstractArray{FT,3}, A::AbstractArray{FT,3}, B::AbstractArray{FT,3}) where {FT}
    Threads.@threads for i = 1:size(A, 3)
        @views ldiv!(X[:,:,i], qr!(A[:,:,i]), B[:,:,i])
    end
end

"Given 3D Julia Array A, fill in X[:,:,k] = A[:,:,k] \\ I"
function batch_inv!(X::AbstractArray{FT,3}, A::AbstractArray{FT,3}) where {FT}
    Threads.@threads for i = 1:size(A, 3)
        @views X[:,:,i] = A[:,:,i]\I
    end
    return X
end

"Given 3D Julia Array A, fill in X[:,:,k] = A[:,:,k] \\ I (pointer version, ignores pointers for CPU)" 
function batch_inv!(X::AbstractArray{FT,3}, A::AbstractArray{FT,3}, xx::Nothing, yy::Nothing) where {FT}
    batch_inv!(X, A)
end

"Workspace-aware batch_inv! (CPU fallback: workspace ignored, just calls plain version)"
function batch_inv!(X::AbstractArray{FT,3}, A::AbstractArray{FT,3}, ws::RTWorkspace) where {FT}
    batch_inv!(X, A)
end

"Batched multiplication for 3D Matrix with 1D vector (i.e. repeated)"
function batched_mul(A::AbstractArray{FT,3}, B::AbstractArray{FT,1}) where {FT}
    return batched_mul(A, reshape(B, (size(B,1), 1)))
end

"""
    batched_mul(A, B)

Threaded CPU override for `NNlib.batched_mul` on dense `Array`s of `BlasFloat`
element type. Parallelizes over the batch (third) axis with one BLAS `gemm!`
per slice. Pin BLAS to one thread (`LinearAlgebra.BLAS.set_num_threads(1)`)
when running with `Threads.nthreads() > 1` to avoid nested BLAS×Julia
parallelism on the small per-slice matrices used in vSmartMOM kernels. The
`ForwardDiff.Dual` element type is *not* `BlasFloat` and falls through to
NNlib's existing generic path unchanged.
"""
function batched_mul(A::Array{T,3}, B::Array{T,3}) where {T<:LinearAlgebra.BLAS.BlasFloat}
    @assert size(A,3) == size(B,3) "batch dim mismatch: $(size(A,3)) vs $(size(B,3))"
    @assert size(A,2) == size(B,1) "inner dim mismatch: $(size(A,2)) vs $(size(B,1))"
    # Singleton-batch defers to NNlib's generic path. The threaded `mul!` loop
    # over `1:1` was empirically producing wrong (≈ 1/50) intensities in
    # multi-layer Stokes_I + per-layer-injection elastic RT — surfaced by the
    # VLIDORT baseline Case B (see dev_notes/nspec1_multilayer_stokesI_bug.md).
    if size(A, 3) == 1
        return invoke(NNlib.batched_mul,
                      Tuple{AbstractArray{T,3}, AbstractArray{T,3}},
                      A, B)
    end
    C = Array{T,3}(undef, size(A,1), size(B,2), size(A,3))
    Threads.@threads for k in 1:size(C,3)
        @views mul!(C[:,:,k], A[:,:,k], B[:,:,k])
    end
    return C
end
