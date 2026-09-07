using Test, Random, LinearAlgebra, CUDA, vSmartMOM
const _JT = vSmartMOM.CoreRT
CUDA.functional() || error("CUDA required for tiled Jacobian tests")
CUDA.allowscalar(false)
@testset "Shared-memory Jacobian products" begin
    old = _JT._TILED_JACOBIANS_ENABLED[]
    try
        for FT in (Float32,Float64), n in (6,9,12,15,16,18,32)
            ns, np = 512, 3
            A = CUDA.rand(FT,n,n,ns); B = CUDA.rand(FT,n,n,ns)
            da = CUDA.rand(FT,n,n,ns,np+1); db = CUDA.rand(FT,n,n,ns,np+1)
            dA = @view da[:,:,:,1:np]; dB = @view db[:,:,:,1:np]
            out, reference = similar(dA), similar(dA)
            tol = FT === Float32 ? 2e-5 : 1e-12
            _JT._TILED_JACOBIANS_ENABLED[] = false
            _JT._jprod!(reference,A,dA,B,dB)
            _JT._TILED_JACOBIANS_ENABLED[] = true
            @test _JT._use_jacobian_tiles(CUDA.CUDABackend(),out,A,B)
            _JT._jprod!(out,A,dA,B,dB)
            @test out ≈ reference rtol=tol
            for (left,right) in ((A,dB),(dA,B)), β in (zero(FT),FT(0.5))
                out .= FT(0.3); reference .= FT(0.3)
                _JT._TILED_JACOBIANS_ENABLED[] = false
                _JT._jmul!(reference,left,right,β)
                _JT._TILED_JACOBIANS_ENABLED[] = true
                _JT._jmul!(out,left,right,β)
                @test out ≈ reference rtol=tol
            end
        end
    finally
        _JT._TILED_JACOBIANS_ENABLED[] = old
    end
end

@testset "Fused Jacobian geometric inverse" begin
    old = _JT._JACOBIAN_FUSED_INVERSE_ENABLED[]
    try
        for FT in (Float32,Float64), n in (16,18,32)
            ns = 512
            A = CUDA.rand(FT,n,n,ns) .* FT(0.02)
            B = CUDA.rand(FT,n,n,ns) .* FT(0.03)
            ident = Diagonal(CUDA.ones(FT,n))
            tmp, G, reference = similar(A), similar(A), similar(A)
            _JT._JACOBIAN_FUSED_INVERSE_ENABLED[] = false
            _JT._jacobian_geometric_inverse!(reference,tmp,A,B,ident)
            _JT._JACOBIAN_FUSED_INVERSE_ENABLED[] = true
            _JT._jacobian_geometric_inverse!(G,tmp,A,B,ident)
            tol = FT === Float32 ? FT(2e-5) : FT(1e-12)
            @test G ≈ reference rtol=tol
            _JT._bmm!(tmp,A,B); tmp .= ident .- tmp
            _JT._bmm!(reference,G,tmp)
            @test maximum(abs,reference .- ident) < tol
        end
    finally
        _JT._JACOBIAN_FUSED_INVERSE_ENABLED[] = old
    end
end
