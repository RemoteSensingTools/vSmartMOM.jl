# Keep the reference propagation available for numerical/performance A/B checks.
const _BATCHED_JACOBIANS_ENABLED = Ref(true)

function make_jacobian_workspace(A::AbstractArray{FT,3}, nparams) where {FT}
    FT <: Union{Float32,Float64} || return nothing
    # The portable GPU products are intended for small diffuse operators.
    # Larger operators retain vendor-BLAS propagation until benchmarked.
    backend = KernelAbstractions.get_backend(A)
    backend isa KernelAbstractions.CPU || size(A, 1) <= 32 || return nothing
    n, _, ns = size(A)
    matrix() = similar(A)
    vector() = similar(A, n, 1, ns)
    tangent() = similar(A, n, n, ns, nparams)
    source_tangent() = similar(A, n, 1, ns, nparams)
    return JacobianPropagationWorkspace(
        ntuple(_ -> matrix(), 6)..., ntuple(_ -> vector(), 4)...,
        ntuple(_ -> tangent(), 6)..., ntuple(_ -> source_tangent(), 4)...)
end

@inline _use_batched_jacobians(a::AddedLayerLin) =
    _BATCHED_JACOBIANS_ENABLED[] && a.propagation_workspace !== nothing

@inline _jac_active(A, n) = size(A, 4) == n ? A : view(A, :, :, :, 1:n)
@inline _jac_entry(A::AbstractArray{T,3}, i, j, s, p) where {T} = @inbounds A[i,j,s]
@inline _jac_entry(A::AbstractArray{T,4}, i, j, s, p) where {T} = @inbounds A[i,j,s,p]
@inline _jac_slice(A::AbstractArray{T,3}, s, p) where {T} = view(A, :, :, s)
@inline _jac_slice(A::AbstractArray{T,4}, s, p) where {T} = view(A, :, :, s, p)

# Fold (wavelength, parameter) into one launch axis. A 3D operand is shared
# across parameters without repeat(), pointer arrays, or materialized views.
@kernel function _jac_mul_kernel!(C, @Const(A), @Const(B), β, ::Val{K}) where {K}
    i, j, batch = @index(Global, NTuple)
    s = mod1(batch, size(C,3))
    p = (batch - 1) ÷ size(C,3) + 1
    value = zero(eltype(C))
    @inbounds for k in 1:K
        value += _jac_entry(A,i,k,s,p) * _jac_entry(B,k,j,s,p)
    end
    @inbounds C[i,j,s,p] = iszero(β) ? value : value + β*C[i,j,s,p]
end

# Fusing both product-rule terms saves one launch and one full tangent temp.
@kernel function _jac_product_kernel!(C, @Const(A), @Const(dA), @Const(B), @Const(dB),
                                      ::Val{K}) where {K}
    i, j, batch = @index(Global, NTuple)
    s = mod1(batch, size(C,3))
    p = (batch - 1) ÷ size(C,3) + 1
    left = zero(eltype(C))
    right = zero(eltype(C))
    @inbounds for k in 1:K
        left += dA[i,k,s,p] * B[k,j,s]
        right += A[i,k,s] * dB[k,j,s,p]
    end
    @inbounds C[i,j,s,p] = left + right
end

function _jac_cpu_batches!(f, C, k)
    nbatch = size(C,3) * size(C,4)
    # Avoid a thread barrier for every tiny product. Parallel work is grouped
    # over both spectral and parameter axes, with BLAS pinned by the RT driver.
    if Threads.nthreads() > 1 && length(C)*k >= 2_000_000
        Threads.@threads for b in 1:nbatch
            f(mod1(b, size(C,3)), (b-1) ÷ size(C,3) + 1)
        end
    else
        for b in 1:nbatch
            f(mod1(b, size(C,3)), (b-1) ÷ size(C,3) + 1)
        end
    end
    return C
end

"C = A*B + β*C over wavelength × parameter; C must not alias either operand."
function _jmul!(C, A, B, β=zero(eltype(C)))
    isempty(C) && return C
    backend = KernelAbstractions.get_backend(C)
    if backend isa KernelAbstractions.CPU
        _jac_cpu_batches!(C, size(A,2)) do s, p
            mul!(_jac_slice(C,s,p), _jac_slice(A,s,p), _jac_slice(B,s,p), one(eltype(C)), β)
        end
    else
        _jac_mul_kernel!(backend)(C, A, B, β, Val(size(A,2));
            ndrange=(size(C,1),size(C,2),size(C,3)*size(C,4)))
    end
    return C
end

"dC = dA*B + A*dB for every physical parameter; dC must not alias inputs."
function _jprod!(dC, A, dA, B, dB)
    isempty(dC) && return dC
    backend = KernelAbstractions.get_backend(dC)
    if backend isa KernelAbstractions.CPU
        _jac_cpu_batches!(dC, size(A,2)) do s, p
            C = _jac_slice(dC,s,p)
            mul!(C, _jac_slice(dA,s,p), _jac_slice(B,s,p))
            mul!(C, _jac_slice(A,s,p), _jac_slice(dB,s,p), one(eltype(C)), one(eltype(C)))
        end
    else
        _jac_product_kernel!(backend)(dC, A, dA, B, dB, Val(size(A,2));
            ndrange=(size(dC,1),size(dC,2),size(dC,3)*size(dC,4)))
    end
    return dC
end

function doubling_batched_lin!(pol_type, expk, ndoubl, a, da, I_static, dτ, μ₀;
                               N_active=0)
    ndoubl == 0 && return nothing
    w = da.propagation_workspace
    n = N_active > 0 ? N_active : size(da.ap_ṙ⁻⁺,4)
    r, t, jp, jm = a.r⁻⁺, a.t⁺⁺, a.j₀⁺, a.j₀⁻
    dr, dt = _jac_active(da.ap_ṙ⁻⁺,n), _jac_active(da.ap_ṫ⁺⁺,n)
    djp, djm = _jac_active(da.ap_J̇₀⁺,n), _jac_active(da.ap_J̇₀⁻,n)
    dG, dH = _jac_active(w.dG,n), _jac_active(w.dH,n)
    dm1, dm2 = _jac_active(w.dm1,n), _jac_active(w.dm2,n)
    dR, dT = _jac_active(w.dR,n), _jac_active(w.dT,n)
    dv1, dv2 = _jac_active(w.dv1,n), _jac_active(w.dv2,n)
    dJm, dJp = _jac_active(w.dJminus,n), _jac_active(w.dJplus,n)
    dJ1m = _jac_active(da.dbl_ap_J̇₁⁻,n)
    de = view(da.dbl_ap_expk_lin, :, 1:n)
    de .= -reshape(expk,:,1) .* view(dτ,:,1:n) ./ μ₀
    e3, e4 = reshape(expk,1,1,:), reshape(expk,1,1,:,1)
    de4 = reshape(de,1,1,size(de,1),n)
    for _ in 1:ndoubl
        _bmm!(w.m1, r, r)
        w.m1 .= I_static .- w.m1
        batch_inv!(w.G, w.m1)
        _bmm!(w.H, t, w.G)
        _jprod!(dm1, r, dr, r, dr)
        _jmul!(dm2, w.G, dm1)
        _jmul!(dG, dm2, w.G)
        _jprod!(dH, t, dt, w.G, dG)

        # Preserve both old sources and their tangents until the complete
        # source update is formed. In-place overwrites here change the physics.
        w.Jplus .= jp .* e3
        w.Jminus .= jm .* e3
        dJp .= djp .* e4 .+ jp .* de4
        dJ1m .= djm .* e4 .+ jm .* de4
        _bmm!(w.v1, r, jp); w.v1 .+= w.Jminus
        _bmm!(w.v2, r, w.Jminus); w.v2 .+= jp
        _jprod!(dv1, r, dr, jp, djp); dv1 .+= dJ1m
        _jprod!(dv2, w.H, dH, w.v1, dv1)
        dJm .= djm .+ dv2
        _jprod!(dv1, r, dr, w.Jminus, dJ1m)
        dv1 .+= djp
        _jprod!(dv2, w.H, dH, w.v2, dv1)
        dJp .+= dv2
        _bmm!(w.Jminus, w.H, w.v1); w.Jminus .+= jm
        _bmm!(w.Jplus, w.H, w.v2); w.Jplus .+= jp .* e3

        _bmm!(w.m1, r, t)
        _jprod!(dm1, r, dr, t, dt)
        _jprod!(dR, w.H, dH, w.m1, dm1); dR .+= dr
        _jprod!(dT, w.H, dH, t, dt)
        _bmm!(w.R, w.H, w.m1); w.R .+= r
        _bmm!(w.T, w.H, t)
        r .= w.R; t .= w.T
        dr .= dR; dt .= dT
        jp .= w.Jplus; jm .= w.Jminus
        djp .= dJp; djm .= dJm
        de .*= 2 .* reshape(expk,:,1)
        expk .*= expk
    end
    apply_D_matrix!(pol_type.n, r, t, a.r⁺⁻, a.t⁻⁻,
        da.ap_ṙ⁻⁺, da.ap_ṫ⁺⁺, da.ap_ṙ⁺⁻, da.ap_ṫ⁻⁻)
    apply_D_matrix_SFI!(pol_type.n, jm, da.ap_J̇₀⁻)
    return nothing
end

# Compute one direction of the general adding equations. The caller defers
# committing the first direction until the second has consumed all old state.
function _interaction_direction_lin!(w, R, dR, Rout, dRout, Tleft, dTleft,
                                     Tright, dTright, Tother, dTother,
                                     Jbase, dJbase, Jin, dJin, j, dj, I_static)
    _bmm!(w.m1, R, Rout)
    w.m1 .= I_static .- w.m1
    batch_inv!(w.G, w.m1)
    _jprod!(w.dm1, R, dR, Rout, dRout)
    _jmul!(w.dm2, w.G, w.dm1)
    _jmul!(w.dG, w.dm2, w.G)
    _bmm!(w.H, Tleft, w.G)
    _jprod!(w.dH, Tleft, dTleft, w.G, w.dG)
    _bmm!(w.m1, R, Tright)
    _jprod!(w.dm1, R, dR, Tright, dTright)
    _jprod!(w.dm2, w.H, w.dH, w.m1, w.dm1)
    _bmm!(w.m2, w.H, w.m1)
    _jprod!(w.dm1, w.H, w.dH, Tother, dTother)
    _bmm!(w.T, w.H, Tother)
    _bmm!(w.v1, R, Jin); w.v1 .+= j
    _jprod!(w.dv1, R, dR, Jin, dJin); w.dv1 .+= dj
    _jprod!(w.dv2, w.H, w.dH, w.v1, w.dv1); w.dv2 .+= dJbase
    _bmm!(w.v2, w.H, w.v1); w.v2 .+= Jbase
    return nothing
end

function interaction_batched_lin!(c, dc, a, da, I_static)
    w = da.propagation_workspace
    # Upward output: composite reflection and downward-facing transmission.
    _interaction_direction_lin!(w, a.r⁻⁺, da.ap_ṙ⁻⁺, c.R⁺⁻, dc.Ṙ⁺⁻,
        c.T⁻⁻, dc.Ṫ⁻⁻, c.T⁺⁺, dc.Ṫ⁺⁺, a.t⁻⁻, da.ap_ṫ⁻⁻,
        c.J₀⁻, dc.J̇₀⁻, c.J₀⁺, dc.J̇₀⁺, a.j₀⁻, da.ap_J̇₀⁻, I_static)
    w.R .= c.R⁻⁺ .+ w.m2; w.dR .= dc.Ṙ⁻⁺ .+ w.dm2
    # Save the first direction outside scratch reused for the second.
    w.dT .= w.dm1
    w.Jminus .= w.v2; w.dJminus .= w.dv2
    # w.m1 is scratch in the next call, so the first forward T is saved in
    # added.temp1, which is dead between doubling calls.
    a.temp1 .= w.T

    _interaction_direction_lin!(w, c.R⁺⁻, dc.Ṙ⁺⁻, a.r⁻⁺, da.ap_ṙ⁻⁺,
        a.t⁺⁺, da.ap_ṫ⁺⁺, a.t⁻⁻, da.ap_ṫ⁻⁻, c.T⁺⁺, dc.Ṫ⁺⁺,
        a.j₀⁺, da.ap_J̇₀⁺, a.j₀⁻, da.ap_J̇₀⁻, c.J₀⁺, dc.J̇₀⁺, I_static)
    c.R⁺⁻ .= a.r⁺⁻ .+ w.m2; dc.Ṙ⁺⁻ .= da.ap_ṙ⁺⁻ .+ w.dm2
    c.T⁺⁺ .= w.T; dc.Ṫ⁺⁺ .= w.dm1
    c.J₀⁺ .= w.v2; dc.J̇₀⁺ .= w.dv2
    c.R⁻⁺ .= w.R; dc.Ṙ⁻⁺ .= w.dR
    c.T⁻⁻ .= a.temp1; dc.Ṫ⁻⁻ .= w.dT
    c.J₀⁻ .= w.Jminus; dc.J̇₀⁻ .= w.dJminus
    return nothing
end
