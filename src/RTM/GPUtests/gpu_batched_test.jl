using CUDA

elty = Float32
FT = Float32
n = 100
nSpec = 10000


C_ = (randn(FT, n, n, nSpec));
B_ = (randn(FT, n, n, nSpec));
A_ = (randn(FT, n, n, nSpec));

C = CuArray(C_);
B = CuArray(B_);
A = CuArray(A_);



function test(d_A, d_B, d_C)
    pivot, info      = CUBLAS.getrf_strided_batched!(d_A, true)
    getri_strided_batched(d_A, d_C, pivot)
    d_C = CUBLAS.gemm_strided_batched('N','N',d_A,d_B)
    #CUBLAS.gemm_strided_batched('N','N',1.0, d_A,d_B, )
    synchronize()
end
@time test(A,B,C)


A = CuArray(randn(FT, m, m, nSpec));
function test2(d_A)
    pivot, info = CUBLAS.getrf_strided_batched!(A, false);
    synchronize()
end
@time test2(A)

function testLU(d_A)
    pivot, info      = CUBLAS.getrf_strided_batched!(d_A, true)
    synchronize()
end

# CUDA has no strided batched getrf, but we can at least avoid constructing costly views
function getri_strided_batched(A::AbstractArray{Float32, 3},C::AbstractArray{Float32, 3},pivotArray::CuMatrix{Cint})
    m,n = size(A,1), size(A,2)
    if m != n
        throw(DimensionMismatch("All matrices must be square!"))
    end
    n = size(A,1)
    ldc = max(1,stride(C,2))
    lda = max(1,stride(A,2))
    info = CUDA.zeros(Cint,size(A,3))
    Cptrs = CUBLAS.unsafe_strided_batch(C)
    Aptrs = CUBLAS.unsafe_strided_batch(A)
    CUBLAS.cublasSgetriBatched(CUBLAS.handle(), n, Aptrs, lda, pivotArray, Cptrs, ldc, info, size(A,3))
    return nothing
end