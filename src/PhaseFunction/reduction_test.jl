using CUDA 
using KernelAbstractions

NMax = 370;

an = (zeros(Complex{Float32},NMax)) 
an .+= 1.3-0.01im
bn = (zeros(Complex{Float32},NMax))
bn .+= 1.2-0.1im
wignerA = (zeros(Complex{Float32},2NMax+1,NMax,NMax)); fill!(wignerA,1.0-0im);
wignerB = (zeros(Complex{Float32},2NMax+1,NMax,NMax)); fill!(wignerB,1.0-0im);
Sνν = zeros(2NMax)

@kernel function compute_Sl_νν!(@Const(wignerA),@Const(wignerB), an, bn, Sνν, @Const(lMax))
    # Get length of an:
    san = size(an)[1];
    # Indices over n and m
    n, m = @index(Global, NTuple)
    
    # Outer loop over l
    for l = 1:lMax
        if max(l-n,n) <= m <= min(san,n+l)
            #println((an[n]' + bn[n]') * (an[m] + bn[m]) * wignerA[l,n,m]^2)
            Sνν[l] += real((an[n]' + bn[n]') * (an[m] + bn[m]) * wignerA[l,n,m]^2)
            Sνν[l] += 1/2 *  real(abs2(an[n] + bn[n]) * wignerA[l,n,m]^2)
        end
    end
    
end

kernel! = compute_Sl_νν!(CPU(),4);
kernelGPU! = compute_Sl_νν!(CUDADevice());
@time event = kernel!(wignerA, wignerB, an, bn, Sνν, 2NMax, ndrange=(NMax,NMax)); wait(event)

wignerA_CU = CuArray(wignerA);
wignerB_CU = CuArray(wignerB);
an_CU      = CuArray(an);
bn_CU      = CuArray(bn);
Sνν_CU     = CuArray(Sνν);
@time event = kernelGPU!(wignerA_CU, wignerB_CU, an_CU, bn_CU, Sνν_CU, 2NMax, ndrange=(NMax,NMax)); wait(event)
