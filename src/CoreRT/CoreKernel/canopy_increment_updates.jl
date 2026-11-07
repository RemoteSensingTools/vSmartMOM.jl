"""
Update solar sources in one launch, using C = T*(I-RR)^(-1)-I and the beam
increment beta. The source update reads the old R and sources throughout.
One workgroup per spectral point, with four shared vectors; matrices are read
from global memory. No near-identity transmission or beam factor is formed.
"""
function ka_canopy_increment_sources!(jd,ju,R,C,beta,backend)
    N = size(R,1)
    kernel! = canopy_increment_sources_kernel!(backend,N)
    kernel!(jd,ju,R,C,beta,Val(N);ndrange=(N*size(R,3),))
    backend isa KernelAbstractions.CPU && KernelAbstractions.synchronize(backend)
    nothing
end

@kernel function canopy_increment_sources_kernel!(jd,ju,@Const(R),@Const(C),
        @Const(beta),::Val{N}) where {N}
    k = @index(Group,Linear)
    i = @index(Local,Linear)
    oldd = @localmem eltype(jd) (N,)
    oldu = @localmem eltype(jd) (N,)
    vu = @localmem eltype(jd) (N,)
    vd = @localmem eltype(jd) (N,)
    @inbounds begin
        oldd[i] = jd[i,1,k]
        oldu[i] = ju[i,1,k]
    end
    @synchronize()
    @inbounds begin
        bd = beta[k]
        reflected_down,reflected_up = zero(eltype(jd)),zero(eltype(jd))
        for j in 1:N
            r = R[i,j,k]
            reflected_down += r*oldd[j]
            reflected_up += r*(oldu[j]+bd*oldu[j])
        end
        vu[i] = reflected_down+(oldu[i]+bd*oldu[i])
        vd[i] = reflected_up+oldd[i]
    end
    @synchronize()
    @inbounds begin
        up,down = zero(eltype(jd)),zero(eltype(jd))
        for j in 1:N
            c = C[i,j,k]
            up += c*vu[j]
            down += c*vd[j]
        end
        ju[i,1,k] = oldu[i]+(vu[i]+up)
        jd[i,1,k] = (oldd[i]+beta[k]*oldd[i])+vd[i]+down
    end
end

"""
Update R, Delta = T-I and beta in one launch. Each workitem owns an output
column. Three shared N×N tiles retain old R, C and Delta; a private vector
retains Q = R[:,j]+R*Delta[:,j]. The pre-write
barrier ensures every thread has captured its old inputs before any mutation.
"""
function ka_canopy_increment_operators!(R,delta,C,beta,backend)
    N = size(R,1)
    kernel! = canopy_increment_operators_kernel!(backend,N)
    kernel!(R,delta,C,beta,Val(N);ndrange=(N*size(R,3),))
    backend isa KernelAbstractions.CPU && KernelAbstractions.synchronize(backend)
    nothing
end

@kernel function canopy_increment_operators_kernel!(R,delta,@Const(C),beta,
        ::Val{N}) where {N}
    k = @index(Group,Linear)
    column = @index(Local,Linear)
    oldR = @localmem eltype(R) (N,N)
    oldC = @localmem eltype(R) (N,N)
    oldDelta = @localmem eltype(R) (N,N)
    q = @private eltype(R) (N,)
    @inbounds for j in 1:N
        oldR[column,j] = R[column,j,k]
        oldC[column,j] = C[column,j,k]
        oldDelta[column,j] = delta[column,j,k]
    end
    @synchronize()
    @inbounds begin
        for i in 1:N
            product = zero(eltype(R))
            for j in 1:N
                product += oldR[i,j]*oldDelta[j,column]
            end
            q[i] = oldR[i,column]+product
        end
        for i in 1:N
            cq,cd = zero(eltype(R)),zero(eltype(R))
            for j in 1:N
                cq += oldC[i,j]*q[j]
                cd += oldC[i,j]*oldDelta[j,column]
            end
            R[i,column,k] = oldR[i,column]+(q[i]+cq)
            delta[i,column,k] = oldC[i,column]+oldDelta[i,column]+cd
        end
        if column==1
            beta[k] = beta[k]*(eltype(beta)(2)+beta[k])
        end
    end
end
