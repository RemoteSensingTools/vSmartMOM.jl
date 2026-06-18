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
    # Cross-platform Float32 tolerance: these compare the KA kernels against a
    # different reference (LAPACK `inv` / naive matmul). Agreement is bounded by
    # Float32 precision, the (unseeded) random matrices' conditioning, and the
    # platform BLAS/LAPACK — Windows OpenBLAS `inv` ran ~6e-5 off the LU kernel,
    # so the old 50*eps(FT) ≈ 6e-6 bound was far too tight. Scale by |ref|.
    @test maximum(abs.(X .- X_ref)) < 1f-3 * maximum(abs.(X_ref))

    C = CoreRT.ka_batched_mul(A, B, backend)
    C_ref = cat([A[:, :, k] * B[:, :, k] for k in axes(A, 3)]...; dims = 3)
    @test maximum(abs.(C .- C_ref)) < 1f-4 * maximum(abs.(C_ref))

    A_pad = rand(FT, n + 1, n + 1, nbatch)
    B_pad = rand(FT, n + 1, 4, nbatch)
    A_view = @view A_pad[2:n+1, 1:n, :]
    B_view = @view B_pad[1:n, 2:4, :]
    C_view = CoreRT.ka_batched_mul(A_view, B_view, backend)
    C_view_ref = cat([A_view[:, :, k] * B_view[:, :, k] for k in axes(A_view, 3)]...; dims = 3)
    @test maximum(abs.(C_view .- C_view_ref)) < 1f-4 * maximum(abs.(C_view_ref))

    @test CoreRT.ka_batch_inv_localmem_bytes(FT, n) == 2 * n * n * sizeof(FT) + n * sizeof(Int32)
    @test_throws ArgumentError CoreRT.ka_batch_inv_lu!(
        similar(A), A, backend; max_localmem_bytes=CoreRT.ka_batch_inv_localmem_bytes(FT, n) - 1
    )
end

@testset "Batched pointer metadata fallback" begin
    A = zeros(Float32, 2, 2, 3)
    @test CoreRT.batched_pointer_cache(A) === nothing
end

@testset "Metal batched kernels local smoke" begin
    if Sys.isapple() && Base.find_package("Metal") !== nothing
        @eval import Metal
        if Metal.functional()
            FT = Float32
            n = 4
            nbatch = 3
            A_cpu = rand(FT, n, n, nbatch)
            B_cpu = rand(FT, n, 2, nbatch)
            for k in axes(A_cpu, 3), i in 1:n
                A_cpu[i, i, k] += one(FT)
            end

            A = Metal.MtlArray(A_cpu)
            B = Metal.MtlArray(B_cpu)

            X = similar(A)
            CoreRT.batch_inv!(X, A)
            X_ref = cat([inv(A_cpu[:, :, k]) for k in axes(A_cpu, 3)]...; dims = 3)
            @test maximum(abs.(Array(X) .- X_ref)) < 1f-3 * maximum(abs.(X_ref))

            C = CoreRT.batched_mul(A, B)
            C_ref = cat([A_cpu[:, :, k] * B_cpu[:, :, k] for k in axes(A_cpu, 3)]...; dims = 3)
            @test maximum(abs.(Array(C) .- C_ref)) < 1f-4 * maximum(abs.(C_ref))
        else
            @info "Skipping local Metal batched-kernel smoke because Metal.functional() is false."
            @test true
        end
    else
        @info "Skipping local Metal batched-kernel smoke on non-Mac or without Metal.jl in the active environment."
        @test true
    end
end
