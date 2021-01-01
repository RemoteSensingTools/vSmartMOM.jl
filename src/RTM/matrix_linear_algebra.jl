# LU decomposition with partial pivoting
# LU decomposed Matrix of A, change in place (A will be destroyed)
function LU_decomposition!(A::AbstractArray{FT,2}, pivot::AbstractArray{Int,1}, n) where {FT}
    local imax = 0;
    @inbounds begin
        for i = 1:n
            pivot[i] = i;
        end
        @unroll for i = 1:n
            maxA = FT(0);
            imax = i;
            @unroll for k = i:n
                if abs(A[k,i]) > maxA
                    maxA = abs(A[k,i]);
                    imax = k;
                end
            end
            # Ensure Matrix is not singular: 
            # @assert maxA > 2eps(FT) "Matrix in LU decompose is singular"

            # interchange the two rows.
            if imax != i
                iP = pivot[i];
                pivot[i] = pivot[imax];
                pivot[imax] = iP
                for j = 1:n 
                    tmp = A[i,j];
                    A[i,j] = A[imax,j];
                    A[imax,j] = tmp;
                end
            end

            # Find the upper triangular matrix elements for row k. 
            @unroll for j = (i + 1):n
                A[j,i] /= A[i,i];
                @unroll for k = (i + 1):n
                    A[j,k] -= A[j,i] * A[i,k];
                end   
            end
        end
    end
    return nothing
end



# Input LU decomposed Matrix of A
# Output: Solution to AX=B (all matrices), Y is just internal storage
function LU_solve!(LU::AbstractArray{FT,2}, B::AbstractArray{FT,2}, X::AbstractArray{FT,2}, Y::AbstractArray{FT,2}, pivot::AbstractArray{Int,1}, n) where {FT}
    fill!(X, FT(0));
    fill!(Y, FT(0));
    
    @inbounds begin
        @unroll for col = 1:n
            # find solution of Ly = b   
            @unroll for i = 1:n
                alpha = FT(0);
                for k = 1:i
                    alpha += LU[i, k] * Y[k,col];
                end
                Y[i, col] = B[pivot[i],col] - alpha;
            end
            # find solution of Ux = y
            @unroll for i = n:(-1):1
                alpha = FT(0);
                for k = (i + 1):n
                    alpha += LU[i, k] * X[k,col];
                end
                X[i,col] = (Y[i,col] - alpha) / LU[i, i];
            end
        end
    end
end

@kernel function LU_solve_kernel!(LU, B, X, Y, pivot, n) where {FT}
    col,batch, = @index(Global,NTuple)
    #@show batch,col
    fill!(X, 0);
    fill!(Y, 0);
    @inbounds begin
        # find solution of Ly = b   
        @unroll for i = 1:n
            for k = 1:i
                Y[i, col,batch] -= LU[i, k, batch] * Y[k,col,batch];
            end
            Y[i, col,batch] += B[pivot[i,batch],col] 
        end
        # find solution of Ux = y
        @unroll for i = n:(-1):1
            for k = (i + 1):n
                X[i,col,batch] -= LU[i, k,batch] * X[k,col,batch];
            end
            X[i,col,batch] = (Y[i,col,batch] + X[i,col,batch]) / LU[i, i,batch];
        end
    end
    nothing
end


"Allocation free matrix multiplication"
function mat_multiply!(out::AbstractArray{FT,2}, in1::AbstractArray{FT,2}, in2::AbstractArray{FT,2}, n) where {FT}
    @inbounds begin
        @unroll for i = 1:n # for each column but the last
            @unroll for j = 1:n
                out[i,j] = in1[i,1] * in2[1,j];
                @unroll for k = 2:n
                    out[i, j] += in1[i, k] * in2[k, j];
                end
            end
        end
    end
end

"Allocation free addition"
function mat_add!(out::AbstractArray{FT,2}, in::AbstractArray{FT,2}, n) where {FT}
    @inbounds begin
        for i = 1:n # for each column but the last
            for j = 1:n
                out[j,i] += in[j,i] ;
            end
        end
    end
end

function mat_add!(out::AbstractArray{FT,2}, in::AbstractArray{FT,2},in2::AbstractArray{FT,2}, n) where {FT}
    @inbounds begin
        for i = 1:n # for each column but the last
            for j = 1:n
                out[j,i] = in[j,i] + in2[j,i] ;
            end
        end
    end
end

"Allocation free matrix multiplication"
function eye_minus_matrix!(out::AbstractArray{FT,3},n) where {FT}
    @inbounds begin
        for i = 1:n # for each column but the last
            for j=1:n
                out[j,i] = -out[j,i]
            end
            out[i,i] += 1 
        end
    end
end

function mat_multiply!(out::AbstractArray{FT,3}, in1::AbstractArray{FT,3}, in2::AbstractArray{FT,3}) where {FT}
    @tensor out[i,k,N] = in1[i,j,N]*in2[j,k,N]
end

# Batched Matrix Multiplication
@macroexpand @kernel function matmul_kernel!(a, b, c, tmp)
    i, j, n = @index(Global, NTuple)
    # creating a temporary sum variable for matrix multiplication
    @unroll for k = 1:size(a,2)
        tmp[n] += a[i,k,n] * b[k, j,n]
    end
    c[i,j,n] = tmp[n]
end

function matmul!(a, b, c)
    @assert size(a,2) == size(b,1) "Matrix size mismatch!"
    for i=1:size(a,3)
        kernel!(view(a,:,:,i), view(b,:,:,i), view(c,:,:,i), ndrange=size(c[:,:,1]))
    end
    synchronize() 
end

@kernel function matmul_kernel!(a, b, c)
    i, j = @index(Global, NTuple)

    # creating a temporary sum variable for matrix multiplication
    tmp_sum = zero(eltype(c))
    for k = 1:size(a)[2]
        tmp_sum += a[i,k] * b[k, j]
    end

    c[i,j] = tmp_sum
end