#=

This file contains CPU implementations of batched linear algebra operations.

These are fallback methods for AbstractArray types (CPU arrays).
GPU-accelerated versions using CuArray are defined in vSmartMOMCUDAExt 
when CUDA is loaded.

=#

# Default synchronization (no-op for CPU, overridden in CUDAExt for GPU)
@inline synchronize() = nothing

"Given 3D Julia Arrays A and B, fill in X[:,:,k] = A[:,:,k] \\ B[:,:,k]" 
function batch_solve!(X::AbstractArray{FT,3}, A::AbstractArray{FT,3}, B::AbstractArray{FT,3}) where {FT}
    for i = 1:size(A, 3)
        @views ldiv!(X[:,:,i], qr!(A[:,:,i]), B[:,:,i])
    end
end

"Given 3D Julia Array A, fill in X[:,:,k] = A[:,:,k] \\ I" 
function batch_inv!(X::AbstractArray{FT,3}, A::AbstractArray{FT,3}) where {FT}
    for i = 1:size(A, 3)
        @views X[:,:,i] = A[:,:,i]\I;
    end
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
