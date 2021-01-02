using LinearAlgebra
using StaticArrays
using BenchmarkTools
using KernelAbstractions
using CUDA

n = 300
FT = Float32

@kernel function rt_interaction_kernel!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, M1, M2, M3)
    N = @index(Global)
    rt_interaction!(R⁻⁺[N], T⁺⁺[N], R⁺⁻[N], T⁻⁻[N], r⁻⁺[N], t⁺⁺[N], r⁺⁻[N], t⁻⁻[N], M1[N], M2[N], M3[N])
    @synchronize
end

function rt_interaction_stupid!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, M1, M2, M3)
    for N = 1:length(R⁻⁺)
        rt_interaction!(R⁻⁺[N], T⁺⁺[N], R⁺⁻[N], T⁻⁻[N], r⁻⁺[N], t⁺⁺[N], r⁺⁻[N], t⁻⁻[N], M1[N], M2[N], M3[N])
    end
    # @synchronize
end


function rt_interaction!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, M1, M2)
    # ToDo: Important output from this routine is R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻ (can be renamed to 𝐓⁻⁻, etc later)
    # Need to check with paper nomenclature. This is basically eqs. 23-28 in vSmartMOM)
    # for i = 1:1000
    M1[:] = (I - R⁺⁻ * r⁻⁺) \ T⁺⁺
    M2[:] = (I - r⁻⁺ * R⁺⁻) \ t⁻⁻

        # repeating for mirror-reflected directions
        
    R⁺⁻[:] = r⁺⁻ + t⁺⁺ * R⁺⁻ * M2 
    T⁻⁻[:] = T⁺⁺ * M2
    T⁺⁺[:] = t⁺⁺ * M1
    R⁻⁺[:] = R⁻⁺ + T⁻⁻ * r⁻⁺ * M1
    return nothing
    # end 
        # return t_R⁻⁺, t_T⁺⁺, t_R⁺⁻, t_T⁻⁻
end
function rt_interaction!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, aux1, aux2, aux3)
    mul!(aux1, R⁺⁻, r⁻⁺)        # R⁺⁻ * r⁻⁺
    aux1[:] = I - aux1   # (I - R⁺⁻ * r⁻⁺)
    aux2[:] = aux1 \ T⁺⁺
    # ldiv!(aux2, qr!(aux1), T⁺⁺) # M1 = (I - R⁺⁻ * r⁻⁺) \ T⁺⁺

        # t_R⁻⁺ = R⁻⁺ + T⁻⁻ * r⁻⁺ * M1
    mul!(aux1, r⁻⁺, aux2)   # r⁻⁺ * M1
    mul!(aux3, T⁻⁻, aux1)   # T⁻⁻ * r⁻⁺ * M1
    R⁻⁺[:] = R⁻⁺ + aux3     # t_R⁻⁺ = R⁻⁺ + T⁻⁻ * r⁻⁺ * M1
        
        # t_T⁺⁺ = t⁺⁺ * M1
    mul!(T⁺⁺, t⁺⁺, aux2)

        # Repeating for mirror-reflected directions

        # M1 = (I - r⁻⁺ * R⁺⁻) \ t⁻⁻
    mul!(aux1, r⁻⁺, R⁺⁻)        # r⁻⁺ * R⁺⁻
    aux1[:] = I - aux1   # (I - r⁻⁺ * R⁺⁻)
    aux2[:] = aux1 \ t⁻⁻
    # ldiv!(aux2, qr!(aux1), t⁻⁻) # M1 = (I - r⁻⁺ * R⁺⁻) \ t⁻⁻

        # t_R⁺⁻ = r⁺⁻ + t⁺⁺ * R⁺⁻ * M1
    mul!(aux3, R⁺⁻, aux2)   # R⁺⁻ * M1
    mul!(aux1, t⁺⁺, aux3)   # t⁺⁺ * R⁺⁻ * M1
    R⁺⁻[:] = r⁺⁻ + aux1     # t_R⁺⁻ = r⁺⁻ + t⁺⁺ * R⁺⁻ * M1

        # t_T⁻⁻ = T⁺⁺ * M1
    mul!(T⁻⁻, T⁺⁺, aux2)
    return nothing
end


r⁻⁺ = CuArray(randn(FT, n, n));
t⁺⁺ = CuArray(randn(FT, n, n));
r⁺⁻ = CuArray(randn(FT, n, n));
t⁻⁻ = CuArray(randn(FT, n, n));
R⁻⁺ = CuArray(randn(FT, n, n));
T⁺⁺ = CuArray(randn(FT, n, n));
R⁺⁻ = CuArray(randn(FT, n, n));
T⁻⁻ = CuArray(randn(FT, n, n));
M1 = similar(r⁺⁻);
M2 = similar(r⁺⁻);
M3 = similar(r⁺⁻);

@btime rt_interaction!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, M1, M2, M3)

nSpec = 1000
AA = [CuArray(randn(FT, n, n)) for i = 1:nSpec]
BB = [CuArray(randn(FT, n, n)) for i = 1:nSpec]
CC = [CuArray(randn(FT, n, n)) for i = 1:nSpec]


_r⁻⁺ = [CuArray(randn(FT, n, n)) for i = 1:nSpec]
_t⁺⁺ = [CuArray(randn(FT, n, n)) for i = 1:nSpec]
_r⁺⁻ = [CuArray(randn(FT, n, n)) for i = 1:nSpec]
_t⁻⁻ = [CuArray(randn(FT, n, n)) for i = 1:nSpec]
_R⁻⁺ = [CuArray(randn(FT, n, n)) for i = 1:nSpec]
_T⁺⁺ = [CuArray(randn(FT, n, n)) for i = 1:nSpec]
_R⁺⁻ = [CuArray(randn(FT, n, n)) for i = 1:nSpec]
_T⁻⁻ = [CuArray(randn(FT, n, n)) for i = 1:nSpec]
_M1 = [CuArray(randn(FT, n, n)) for i = 1:nSpec]
_M2 = [CuArray(randn(FT, n, n)) for i = 1:nSpec]
_M3 = [CuArray(randn(FT, n, n)) for i = 1:nSpec]
kernel_inter! = rt_interaction_kernel!(CUDADevice())
kernel_inter!(_R⁻⁺, _T⁺⁺, _R⁺⁻, _T⁻⁻, _r⁻⁺, _t⁺⁺, _r⁺⁻, _t⁻⁻, _M1, _M2,_M3, ndrange=nSpec)

@kernel function mulCu!(A, B, C)
    N = @index(Global)
    rt_interaction!(R⁻⁺[N], T⁺⁺[N], R⁺⁻[N], T⁻⁻[N], r⁻⁺[N], t⁺⁺[N], r⁺⁻[N], t⁻⁻[N], M1[N], M2[N], M3[N])
    @synchronize
end

@kernel function mat_test!(A, B, C)
    N = @index(Global)
    @tensor A[i, j,N] = B[i, k, N] * C[k, j,N]
end

function testGPU!(A, B, C)
    event = ker!(A, B, C, ndrange=nSpec)
    wait(CUDADevice(), event)
    return nothing
end
function matCPU_test!(A, B, C)
    d = size(C, 3)
    @inbounds for N = 1:d
        A[:, :,N] = B[:, :, N] * C[:, :,N]
    end
end
