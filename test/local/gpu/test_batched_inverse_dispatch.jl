using Test, LinearAlgebra, CUDA
using vSmartMOM: CoreRT

@testset "CUDA inverse without pointer caches" begin
    CUDA.allowscalar(false)
    for FT in (Float32, Float64), batches in (1, 4)
        host = zeros(FT, 3, 3, batches)
        for k in 1:batches
            host[:, :, k] .= FT[4 1 0; 1 3 1; 0 1 2] .+ FT(k) .* Matrix{FT}(I, 3, 3)
        end
        A = CuArray(host)
        X = similar(A)
        CoreRT.batch_inv!(X, A, nothing, nothing)
        result = Array(X)
        for k in 1:batches
            @test host[:, :, k] * result[:, :, k] ≈ Matrix{FT}(I, 3, 3) rtol=50eps(FT)
        end
    end
end
