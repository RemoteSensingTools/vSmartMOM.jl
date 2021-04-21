# Is included in latest CUBLAS version, can be removed in the future!

## getriBatched - performs batched matrix inversion

# Loop over all these CUBLAS functions and create custom overloads of getri_batched!
for (fname, elty) in
    ((:cublasDgetriBatched, :Float64),
     (:cublasSgetriBatched, :Float32),
     (:cublasZgetriBatched, :ComplexF64),
     (:cublasCgetriBatched, :ComplexF32))
    @eval begin
    # cublasStatus_t cublasDgetriBatched(
    #   cublasHandle_t handle, int n, double **A,
    #   int lda, int *PivotArray, double **C,
    #   int ldc, int *info, int batchSize)
        function getri_batched!(n, Aptrs::CuVector{CuPtr{$elty}},
                          lda, Cptrs::CuVector{CuPtr{$elty}},ldc,
                          pivotArray::CuArray{Cint})
            batchSize = length(Aptrs)
            info = CuArray{Cint}(undef, batchSize)
            CUBLAS.$fname(CUBLAS.handle(), n, Aptrs, lda, pivotArray, Cptrs, ldc, info, batchSize)
            CUBLAS.unsafe_free!(Cptrs)
            CUBLAS.unsafe_free!(Aptrs)
            return info
        end
    end
end

# Loop over all these CUBLAS functions and create custom overloads of getri_batched!
for (fname, elty) in
        ((:cublasDgetriBatched, :Float64),
         (:cublasSgetriBatched, :Float32),
         (:cublasZgetriBatched, :ComplexF64),
         (:cublasCgetriBatched, :ComplexF32))
    @eval begin
        # cublasStatus_t cublasDgetriBatched(
        #   cublasHandle_t handle, int n, double **A,
        #   int lda, int *PivotArray, double **C,
        #   int ldc, int *info, int batchSize)
        function getri_batched!(A::Vector{<:CuMatrix{$elty}},
                              C::Vector{<:CuMatrix{$elty}},
                              pivotArray::CuMatrix{Cint})
            n = size(A[1])[1]
            lda = max(1, stride(A[1], 2))
            ldc = max(1, stride(C[1], 2))
            Aptrs = CUBLAS.unsafe_batch(A)
            Cptrs = CUBLAS.unsafe_batch(C)
            info = CuArrays.zeros(Cint, length(A))
            CUBLAS.$fname(CUBLAS.handle(), n, Aptrs, lda, pivotArray, Cptrs, ldc, info, length(A))
            CUBLAS.unsafe_free!(Cptrs)
            CUBLAS.unsafe_free!(Aptrs)

            return info
        end
    end
end

# CUDA has no strided batched getri, but we can at least avoid constructing costly views (copied this over from gertf)
function getri_strided_batched!(A::CuArray{<:Any,3}, C::CuArray{<:Any,3}, pivot::CuArray{Cint})
    m, n = size(A, 1), size(A, 2)
    if m != n
        throw(DimensionMismatch("All matrices must be square!"))
    end
    ldc = max(1, stride(C, 2))
    lda = max(1, stride(A, 2))
    Cptrs = CUBLAS.unsafe_strided_batch(C)
    Aptrs = CUBLAS.unsafe_strided_batch(A)
    return getri_batched!(n, Aptrs, lda, Cptrs, ldc, pivot)
end
