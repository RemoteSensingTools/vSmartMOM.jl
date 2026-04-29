#=
 
Define all the scientific constants needed in this module here
 
=#

# Avogadro's number:
Nₐ = 6.0221415e23;

# Optional CUDA batch pointer provider. Set by vSmartMOMCUDAExt when CUDA is loaded.
# Non-CUDA backends use `batched_pointer_cache(A) === nothing` and dispatch to
# portable batched kernels instead of pointer-array CUBLAS routines.
const CUBLAS_ref = Ref{Any}(nothing)
