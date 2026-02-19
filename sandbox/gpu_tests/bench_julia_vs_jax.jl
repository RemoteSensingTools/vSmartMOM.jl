using CUDA
using LinearAlgebra
using BenchmarkTools
using Printf

function batch_inv_nosync!(X::CuArray{FT,3}, A::CuArray{FT,3}) where {FT}
    pivot, info = CUDA.CUBLAS.getrf_strided_batched!(A, true)
    CUDA.CUBLAS.getri_strided_batched!(A, X, pivot)
    CUDA.synchronize()
end

matrix_sizes = [12, 24, 36]
batch_sizes  = [6837, 10000, 50000]
num_iters    = 20

println("Julia CUDA batched operations benchmark")
println("GPU: ", CUDA.name(CUDA.device()))
println("="^75)

for FT in [Float64, Float32]
    println("\n── $FT ──")
    println()
    @printf("  %-8s  %8s  │ %12s  %10s  │ %12s  %10s\n",
            "size", "batch", "matmul (ms)", "μs/mat", "inv (ms)", "μs/mat")
    println("  ", "─"^71)

    for n in matrix_sizes
        for nBatch in batch_sizes
            # --- Batched multiply ---
            A = CUDA.randn(FT, n, n, nBatch)
            B = CUDA.randn(FT, n, n, nBatch)
            C = similar(A)
            CUDA.CUBLAS.gemm_strided_batched!('N', 'N', one(FT), A, B, zero(FT), C)
            CUDA.synchronize()

            mul_times = Float64[]
            for _ in 1:num_iters
                CUDA.synchronize()
                t = @elapsed begin
                    CUDA.CUBLAS.gemm_strided_batched!('N', 'N', one(FT), A, B, zero(FT), C)
                    CUDA.synchronize()
                end
                push!(mul_times, t)
            end
            mul_med = sort(mul_times)[num_iters ÷ 2] * 1000  # ms

            A = nothing; B = nothing; C = nothing
            GC.gc(); CUDA.reclaim()

            # --- Batched inversion ---
            A_cpu = randn(FT, n, n, nBatch)
            for k in 1:nBatch
                A_cpu[:,:,k] .+= FT(n) * FT.(I(n))
            end
            A_gpu = CuArray(A_cpu)
            X_gpu = similar(A_gpu)

            A_tmp = copy(A_gpu)
            batch_inv_nosync!(X_gpu, A_tmp)

            inv_times = Float64[]
            for _ in 1:num_iters
                A_tmp = copy(A_gpu)
                CUDA.synchronize()
                t = @elapsed begin
                    batch_inv_nosync!(X_gpu, A_tmp)
                end
                push!(inv_times, t)
            end
            inv_med = sort(inv_times)[num_iters ÷ 2] * 1000  # ms

            @printf("  %2d×%-4d  %8d  │ %10.3f ms  %8.3f μs  │ %10.3f ms  %8.3f μs\n",
                    n, n, nBatch,
                    mul_med, mul_med * 1000 / nBatch,
                    inv_med, inv_med * 1000 / nBatch)

            A_gpu = nothing; X_gpu = nothing
            GC.gc(); CUDA.reclaim()
        end
    end
end
