using CUDA
using LinearAlgebra
using BenchmarkTools

function batch_inv_strided_nosync!(X::CuArray{FT,3}, A::CuArray{FT,3}) where {FT}
    pivot, info = CUDA.CUBLAS.getrf_strided_batched!(A, true)
    CUDA.CUBLAS.getri_strided_batched!(A, X, pivot)
    CUDA.synchronize()
end

matrix_sizes = [12, 24, 36]
batch_sizes = [6837, 50000]

test_batches = [6837, 50000]

println("="^70)
println("1) Batch INVERSION (getrf + getri)")
println("="^70)

for FT in [Float64, Float32]
    println("\n  $FT:")
    for n in matrix_sizes
        for nBatch in batch_sizes
            A_cpu = randn(FT, n, n, nBatch)
            for k in 1:nBatch
                A_cpu[:,:,k] .+= FT(n) * FT.(I(n))
            end
            A_gpu = CuArray(A_cpu)
            X_gpu = similar(A_gpu)

            A_tmp = copy(A_gpu)
            batch_inv_strided_nosync!(X_gpu, A_tmp)

            b = @benchmark begin
                A_tmp = copy($A_gpu)
                batch_inv_strided_nosync!($X_gpu, A_tmp)
            end

            med = BenchmarkTools.median(b).time / 1e6
            us_per = med * 1000 / nBatch
            println("    $(n)×$(n) × $(nBatch)  →  $(round(med, digits=2)) ms  ($(round(us_per, digits=3)) μs/matrix)")
            
            A_gpu = nothing; X_gpu = nothing
            GC.gc(); CUDA.reclaim()
        end
    end
end

println("\n\n")
println("="^70)
println("2) Batch MULTIPLICATION (gemm_strided_batched)")
println("="^70)

for FT in [Float64, Float32]
    println("\n  $FT:")
    for n in matrix_sizes
        for nBatch in batch_sizes
            A_gpu = CUDA.randn(FT, n, n, nBatch)
            B_gpu = CUDA.randn(FT, n, n, nBatch)
            C_gpu = similar(A_gpu)

            # Warmup
            CUDA.CUBLAS.gemm_strided_batched!('N', 'N', one(FT), A_gpu, B_gpu, zero(FT), C_gpu)
            CUDA.synchronize()

            b = @benchmark begin
                CUDA.CUBLAS.gemm_strided_batched!('N', 'N', one($FT), $A_gpu, $B_gpu, zero($FT), $C_gpu)
                CUDA.synchronize()
            end

            med = BenchmarkTools.median(b).time / 1e6
            us_per = med * 1000 / nBatch
            println("    $(n)×$(n) × $(nBatch)  →  $(round(med, digits=3)) ms  ($(round(us_per, digits=3)) μs/matrix)")

            A_gpu = nothing; B_gpu = nothing; C_gpu = nothing
            GC.gc(); CUDA.reclaim()
        end
    end
end
