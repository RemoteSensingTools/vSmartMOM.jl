using Test, Random, LinearAlgebra, vSmartMOM
using vSmartMOM.CoreRT
using KernelAbstractions

function check_blocked_jacobians(AT)
    rng = Xoshiro(129)
    for FT in (Float32,Float64), n in (33,48,57,64)
        ns,np = 3,2
        A,B = randn(rng,FT,n,n,ns),randn(rng,FT,n,n,ns)
        dA,dB = randn(rng,FT,n,n,ns,np),randn(rng,FT,n,n,ns,np)
        output = AT(zeros(FT,n,n,ns,np))
        backend = KernelAbstractions.get_backend(output)
        side = 16cld(n,16)
        range = (side,side,ns*np)
        for (left,right) in ((A,dB),(dA,B),(dA,dB))
            fill!(output,FT(0.7))
            CoreRT._jac_mul_blocked!(backend,(16,16,1))(output,AT(left),AT(right),FT(0.3),Val(n);ndrange=range)
            result = Array(output)
            for p in 1:np, s in 1:ns
                l = ndims(left)==3 ? left[:,:,s] : left[:,:,s,p]
                r = ndims(right)==3 ? right[:,:,s] : right[:,:,s,p]
                @test result[:,:,s,p] ≈ l*r .+ FT(0.3)*FT(0.7) rtol=100eps(FT) atol=100eps(FT)
            end
        end
        CoreRT._jac_product_blocked!(backend,(16,16,1))(output,AT(A),AT(dA),AT(B),AT(dB),Val(n);ndrange=range)
        result = Array(output)
        for p in 1:np, s in 1:ns
            @test result[:,:,s,p] ≈ dA[:,:,s,p]*B[:,:,s]+A[:,:,s]*dB[:,:,s,p] rtol=100eps(FT) atol=100eps(FT)
        end
    end
end
@testset "Blocked tangent products CPU" begin
    check_blocked_jacobians(Array)
end
if get(ENV,"VSMARTMOM_JACOBIAN_GPU_TEST","false") == "true"
    using CUDA
    CUDA.allowscalar(false)
    @testset "Blocked tangent products CUDA" begin
        check_blocked_jacobians(CUDA.CuArray)
    end
end
