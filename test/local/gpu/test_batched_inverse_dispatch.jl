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
        @test CoreRT.batch_inv!(X, A, nothing, nothing) === X
        result = Array(X)
        for k in 1:batches
            @test host[:, :, k] * result[:, :, k] ≈ Matrix{FT}(I, 3, 3) rtol=50eps(FT)
        end

        rhs_host = reshape(FT.(1:6batches), 3, 2, batches)
        A = CuArray(host)
        rhs = CuArray(rhs_host)
        solution = similar(rhs)
        @test CoreRT.batch_solve!(solution, A, rhs) === solution
        @test Array(rhs) == rhs_host
        solution_host = Array(solution)
        for k in 1:batches
            @test solution_host[:, :, k] ≈ host[:, :, k] \ rhs_host[:, :, k] rtol=50eps(FT)
        end
    end
end
