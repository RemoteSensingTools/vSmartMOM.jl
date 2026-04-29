using KernelAbstractions

@testset "Portable KA batched kernels" begin
    FT = Float32
    n = 4
    nbatch = 6
    A = rand(FT, n, n, nbatch)
    B = rand(FT, n, 3, nbatch)

    for k in axes(A, 3), i in 1:n
        A[i, i, k] += one(FT)
    end

    backend = KernelAbstractions.CPU()

    X = similar(A)
    CoreRT.ka_batch_inv_lu!(X, A, backend)
    X_ref = cat([inv(A[:, :, k]) for k in axes(A, 3)]...; dims = 3)
    @test maximum(abs.(X .- X_ref)) < 50eps(FT)

    C = CoreRT.ka_batched_mul(A, B, backend)
    C_ref = cat([A[:, :, k] * B[:, :, k] for k in axes(A, 3)]...; dims = 3)
    @test maximum(abs.(C .- C_ref)) < 50eps(FT)
end

@testset "Batched pointer metadata fallback" begin
    A = zeros(Float32, 2, 2, 3)
    @test CoreRT.batched_pointer_cache(A) === nothing
end
