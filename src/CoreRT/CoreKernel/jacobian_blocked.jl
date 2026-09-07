# Blocked product rules for medium polarized operators. Unlike the small
# whole-matrix tiles, these use 16×16 workgroups regardless of matrix size;
# N>32 therefore does not exceed CUDA's 1024-thread workgroup limit.
# The third launch axis combines wavelength and tangent direction. Matrix
# operands shared across directions are read without replication.

@kernel function _jac_mul_blocked!(C,@Const(A),@Const(B),β,::Val{N}) where {N}
    li,lj,_ = @index(Local,NTuple)
    gi,gj,batch = @index(Group,NTuple)
    a = @localmem eltype(C) (17,16)
    b = @localmem eltype(C) (17,16)
    value = @private eltype(C) 1
    value[1] = zero(eltype(C))
    for tile in 0:cld(N,16)-1
        # Reconstruct derived indices inside each barrier phase for KA's
        # CPU lowering; only @private accumulators span tile iterations.
        i,j = 16(gi-1)+li,16(gj-1)+lj
        s,p = mod1(batch,size(C,3)),(batch-1) ÷ size(C,3)+1
        ka,kb = 16tile+lj,16tile+li
        @inbounds begin
            a[li,lj] = i <= N && ka <= N ? _jac_entry(A,i,ka,s,p) : zero(eltype(C))
            b[li,lj] = kb <= N && j <= N ? _jac_entry(B,kb,j,s,p) : zero(eltype(C))
        end
        @synchronize
        @inbounds for k in 1:16
            value[1] += a[li,k]*b[k,lj]
        end
        @synchronize
    end
    i,j = 16(gi-1)+li,16(gj-1)+lj
    s,p = mod1(batch,size(C,3)),(batch-1) ÷ size(C,3)+1
    if i <= N && j <= N
        @inbounds C[i,j,s,p] = iszero(β) ? value[1] : value[1]+β*C[i,j,s,p]
    end
end

# S2014 (C.6): d(AB)=dA B+A dB. Keep the two sums separate so their
# association agrees with the small-operator and vendor-BLAS reference paths.
@kernel function _jac_product_blocked!(C,@Const(A),@Const(dA),@Const(B),@Const(dB),
                                      ::Val{N}) where {N}
    li,lj,_ = @index(Local,NTuple)
    gi,gj,batch = @index(Group,NTuple)
    a = @localmem eltype(C) (17,16)
    da = @localmem eltype(C) (17,16)
    b = @localmem eltype(C) (17,16)
    db = @localmem eltype(C) (17,16)
    value = @private eltype(C) 2
    value[1] = zero(eltype(C)); value[2] = zero(eltype(C))
    for tile in 0:cld(N,16)-1
        # Reconstruct derived indices inside each barrier phase for KA's
        # CPU lowering; only @private accumulators span tile iterations.
        i,j = 16(gi-1)+li,16(gj-1)+lj
        s,p = mod1(batch,size(C,3)),(batch-1) ÷ size(C,3)+1
        ka,kb = 16tile+lj,16tile+li
        @inbounds begin
            a[li,lj] = i <= N && ka <= N ? _jac_entry(A,i,ka,s,p) : zero(eltype(C))
            da[li,lj] = i <= N && ka <= N ? _jac_entry(dA,i,ka,s,p) : zero(eltype(C))
            b[li,lj] = kb <= N && j <= N ? _jac_entry(B,kb,j,s,p) : zero(eltype(C))
            db[li,lj] = kb <= N && j <= N ? _jac_entry(dB,kb,j,s,p) : zero(eltype(C))
        end
        @synchronize
        @inbounds for k in 1:16
            value[1] += da[li,k]*b[k,lj]
            value[2] += a[li,k]*db[k,lj]
        end
        @synchronize
    end
    i,j = 16(gi-1)+li,16(gj-1)+lj
    s,p = mod1(batch,size(C,3)),(batch-1) ÷ size(C,3)+1
    if i <= N && j <= N
        @inbounds C[i,j,s,p] = value[1]+value[2]
    end
end
