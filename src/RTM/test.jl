using LinearAlgebra
using StaticArrays
using BenchmarkTools
using KernelAbstractions
using CUDA

n = 30
FT = Float32

@kernel function rt_interaction_kernel!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, M1, M2, M3)
    N = @index(Global)
    rt_interaction!(R⁻⁺[N], T⁺⁺[N], R⁺⁻[N], T⁻⁻[N], r⁻⁺[N], t⁺⁺[N], r⁺⁻[N], t⁻⁻[N], M1[N], M2[N], M3[N])
    @synchronize
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
dτ = FT(0.001);
ndoubl = 12;
τ_total = FT(2^ndoubl * dτ);

if n < 20
    r⁻⁺ = @SMatrix randn(Float32, n, n);
    t⁺⁺ = @SMatrix randn(Float32, n, n);
    r⁺⁻ = @SMatrix randn(Float32, n, n);
    t⁻⁻ = @SMatrix randn(Float32, n, n);



    @btime rt_doubling!(dτ, τ_total, ndoubl, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
end

r⁻⁺ = randn(FT, n, n);
t⁺⁺ = randn(FT, n, n);
r⁺⁻ = randn(FT, n, n);
t⁻⁻ = randn(FT, n, n);
R⁻⁺ = randn(FT, n, n);
T⁺⁺ = randn(FT, n, n);
R⁺⁻ = randn(FT, n, n);
T⁻⁻ = randn(FT, n, n);
M1 = similar(r⁺⁻);
M2 = similar(r⁺⁻);
M3 = similar(r⁺⁻);

@btime rt_doubling!(dτ, τ_total, ndoubl, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)

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

nSpec = 1000
_r⁻⁺ = [r⁻⁺ for i = 1:nSpec]
_t⁺⁺ = [t⁺⁺ for i = 1:nSpec]
_r⁺⁻ = [r⁺⁻ for i = 1:nSpec]
_t⁻⁻ = [t⁻⁻ for i = 1:nSpec]
_R⁻⁺ = [R⁻⁺ for i = 1:nSpec]
_T⁺⁺ = [T⁺⁺ for i = 1:nSpec]
_R⁺⁻ = [R⁺⁻ for i = 1:nSpec]
_T⁻⁻ = [T⁻⁻ for i = 1:nSpec]
_M1 = [similar(T⁻⁻) for i = 1:nSpec]
_M2 = [similar(T⁻⁻) for i = 1:nSpec]
_M3 = [similar(T⁻⁻) for i = 1:nSpec]
kernel_inter! = rt_interaction_kernel!(CUDADevice())
kernel_inter!(_R⁻⁺, _T⁺⁺, _R⁺⁻, _T⁻⁻, _r⁻⁺, _t⁺⁺, _r⁺⁻, _t⁻⁻, _M1, _M2,_M3, ndrange=nSpec)

@btime rt_doubling!(dτ, τ_total, ndoubl, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
