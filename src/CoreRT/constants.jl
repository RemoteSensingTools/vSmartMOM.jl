#=
 
Define all the scientific constants needed in this module here
 
=#

# Avogadro's number:
Nₐ = 6.0221415e23;

# Optional GPU batch pointer provider. Set by vSmartMOMCUDAExt when CUDA is loaded.
# When nothing, GPU arrays (CuArray) cannot be used for RT layers that need batched pointers.
const CUBLAS_ref = Ref{Any}(nothing)
