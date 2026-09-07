# Scientific notation used below (at fixed wavelength and Fourier order):
# dA means ∂A/∂p for one supplied tangent direction p. Matrix multiplication
# contracts only the quadrature × Stokes index, never wavelength or parameter.
# A[i,j,s] is shared by all p; dA[i,j,s,p] holds independent tangent matrices.
# These may be physical-parameter tangents or supplied core-optics directions.
#
# References: S2014 = Sanghavi, Davis & Eldering, JQSRT 133 (2014), 412–433;
# SF2023-II = Sanghavi & Frankenberg, JQSRT 311 (2023), 108791.
# S2014 (C.6): d(AB) = dA B + A dB.
# S2014 (C.7): d(A⁻¹) = -A⁻¹ dA A⁻¹ (order cannot be interchanged).
# See docs/dev_notes/theory_references.md for the verified paper/code map.

# Keep the reference propagation available for numerical/performance A/B checks.
const _BATCHED_JACOBIANS_ENABLED = Ref(true)
const _TILED_JACOBIANS_ENABLED = Ref(true)
const _MEDIUM_JACOBIANS_ENABLED = Ref(true)
const _JACOBIAN_FUSED_INVERSE_ENABLED = Ref(true)
_jacobian_tiles_supported(::Any) = false

@inline _use_jacobian_tiles(backend, C, A, B) =
    _TILED_JACOBIANS_ENABLED[] && _jacobian_tiles_supported(backend) &&
    16 <= size(C,1) <= 32 && size(C,3) >= 512 &&
    size(C,1) == size(C,2) == size(A,1) == size(A,2) == size(B,1) == size(B,2)

@inline _use_blocked_jacobians(backend,C,A,B) =
    _MEDIUM_JACOBIANS_ENABLED[] && _jacobian_tiles_supported(backend) &&
    32 < size(C,1) <= 64 && size(C,3) >= 512 &&
    size(C,1) == size(C,2) == size(A,1) == size(A,2) == size(B,1) == size(B,2)

function make_jacobian_workspace(A::AbstractArray{FT,3}, nparams) where {FT}
    FT <: Union{Float32,Float64} || return nothing
    # CUDA tiles are validated for operators up to 64 with at least 512
    # spectral points. Smaller batches above 32 and larger operators keep
    # vendor-BLAS propagation; Metal retains the existing 32-operator limit.
    backend = KernelAbstractions.get_backend(A)
    gpu_limit = _MEDIUM_JACOBIANS_ENABLED[] && _jacobian_tiles_supported(backend) &&
                size(A,3) >= 512 ? 64 : 32
    backend isa KernelAbstractions.CPU || size(A, 1) <= gpu_limit || return nothing
    n, _, ns = size(A)
    matrix() = similar(A)
    vector() = similar(A, n, 1, ns)
    tangent() = similar(A, n, n, ns, nparams)
    source_tangent() = similar(A, n, 1, ns, nparams)
    return JacobianPropagationWorkspace(
        ntuple(_ -> matrix(), 6)..., ntuple(_ -> vector(), 4)...,
        ntuple(_ -> tangent(), 5)..., ntuple(_ -> source_tangent(), 4)...)
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

# One wavelength/parameter matrix per workgroup. The Stokes-expanded square
# products reuse each loaded operand across all output columns/rows instead
# of reloading it inside every dot product. Source-vector products retain the
# simpler kernel because almost all tile threads would be idle there.
@kernel function _jac_mul_tiled!(C, @Const(A), @Const(B), β, ::Val{N}) where {N}
    batch = @index(Group, Linear)
    tid = @index(Local, Linear)
    s, p = mod1(batch,size(C,3)), (batch-1) ÷ size(C,3) + 1
    i, j = mod1(tid,N), (tid-1) ÷ N + 1
    a = @localmem eltype(C) (N,N)
    b = @localmem eltype(C) (N,N)
    @inbounds a[i,j] = _jac_entry(A,i,j,s,p)
    @inbounds b[i,j] = _jac_entry(B,i,j,s,p)
    @synchronize
    value = zero(eltype(C))
    @inbounds for k in 1:N
        value += a[i,k] * b[k,j]
    end
    @inbounds C[i,j,s,p] = iszero(β) ? value : value + β*C[i,j,s,p]
end

@kernel function _jac_product_tiled!(C, @Const(A), @Const(dA), @Const(B), @Const(dB),
                                     ::Val{N}) where {N}
    batch = @index(Group, Linear)
    tid = @index(Local, Linear)
    s, p = mod1(batch,size(C,3)), (batch-1) ÷ size(C,3) + 1
    i, j = mod1(tid,N), (tid-1) ÷ N + 1
    a = @localmem eltype(C) (N,N)
    da = @localmem eltype(C) (N,N)
    b = @localmem eltype(C) (N,N)
    db = @localmem eltype(C) (N,N)
    @inbounds begin
        a[i,j] = A[i,j,s]; da[i,j] = dA[i,j,s,p]
        b[i,j] = B[i,j,s]; db[i,j] = dB[i,j,s,p]
    end
    @synchronize
    left, right = zero(eltype(C)), zero(eltype(C))
    @inbounds for k in 1:N
        left += da[i,k] * b[k,j]
        right += a[i,k] * db[k,j]
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
    elseif _use_blocked_jacobians(backend,C,A,B)
        n = size(C,1)
        side = 16cld(n,16)
        _jac_mul_blocked!(backend,(16,16,1))(C,A,B,β,Val(n);
            ndrange=(side,side,size(C,3)*size(C,4)))
    elseif _use_jacobian_tiles(backend,C,A,B)
        n = size(C,1)
        _jac_mul_tiled!(backend,n*n)(C,A,B,β,Val(n);
            ndrange=n*n*size(C,3)*size(C,4))
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
    elseif _use_blocked_jacobians(backend,dC,A,B)
        n = size(dC,1)
        side = 16cld(n,16)
        _jac_product_blocked!(backend,(16,16,1))(dC,A,dA,B,dB,Val(n);
            ndrange=(side,side,size(dC,3)*size(dC,4)))
    elseif _use_jacobian_tiles(backend,dC,A,B)
        n = size(dC,1)
        _jac_product_tiled!(backend,n*n)(dC,A,dA,B,dB,Val(n);
            ndrange=n*n*size(dC,3)*size(dC,4))
    else
        _jac_product_kernel!(backend)(dC, A, dA, B, dB, Val(size(A,2));
            ndrange=(size(dC,1),size(dC,2),size(dC,3)*size(dC,4)))
    end
    return dC
end

# G = (I-R₁R₂)⁻¹. The forward fused right-solve computes X(I-R₁R₂)=B;
# choosing B=I yields exactly G, which all supplied tangent directions share.
# Use the same pivoted LU kernel as the forward solve instead of returning
# to cuBLAS getrf/getri for small polarized CUDA matrices. The temporary
# identity is separate from G: the fused kernel must read the old RHS.
function _jacobian_geometric_inverse!(G, tmp, R₁, R₂, I_static)
    backend = KernelAbstractions.get_backend(G)
    if _JACOBIAN_FUSED_INVERSE_ENABLED[] && _jacobian_tiles_supported(backend) &&
       size(G,1) <= 32 && _use_fused_solve(G)
        tmp .= I_static
        ka_fused_solve!(G, R₁, R₂, tmp, backend)
    else
        _bmm!(tmp, R₁, R₂)
        tmp .= I_static .- tmp
        batch_inv!(G, tmp)
    end
    return G
end

function doubling_batched_lin!(pol_type, expk, ndoubl, a, da, I_static, dτ, μ₀;
                               N_active=0)
    ndoubl == 0 && return nothing
    w = da.propagation_workspace
    n = N_active > 0 ? N_active : size(da.ap_ṙ⁻⁺,4)
    r, t, jp, jm = a.r⁻⁺, a.t⁺⁺, a.j₀⁺, a.j₀⁻
    dr, dt = _jac_active(da.ap_ṙ⁻⁺,n), _jac_active(da.ap_ṫ⁺⁺,n)
    djp, djm = _jac_active(da.ap_J̇₀⁺,n), _jac_active(da.ap_J̇₀⁻,n)
    dH = _jac_active(w.dH,n)
    dm1, dm2 = _jac_active(w.dm1,n), _jac_active(w.dm2,n)
    dR, dT = _jac_active(w.dR,n), _jac_active(w.dT,n)
    dv1, dv2 = _jac_active(w.dv1,n), _jac_active(w.dv2,n)
    dJm, dJp = _jac_active(w.dJminus,n), _jac_active(w.dJplus,n)
    dJ1m = _jac_active(da.dbl_ap_J̇₁⁻,n)
    de = view(da.dbl_ap_expk_lin, :, 1:n)
    # Direct-beam attenuation: e = exp(-dτ/μ₀), de = -e d(dτ)/μ₀.
    # dτ here already contains the supplied elemental-thickness tangents;
    # geometry and the integer doubling count are held fixed.
    de .= -reshape(expk,:,1) .* view(dτ,:,1:n) ./ μ₀
    e3, e4 = reshape(expk,1,1,:), reshape(expk,1,1,:,1)
    de4 = reshape(de,1,1,size(de,1),n)
    for _ in 1:ndoubl
        # D-transformed homogeneous layer: r is the starred reflection of
        # S2014 (31), so G = (I-r r)⁻¹ and H = t G. From (C.6)–(C.7):
        # dG = G (dr r + r dr) G; dH = dt G + t dG.
        # Substitute H=tG to evaluate dH = [dt + H(dr r + r dr)]G.
        # This avoids materializing dG and saves two tangent products while
        # keeping matrix order intact. The two inverse-rule signs cancel.
        _jacobian_geometric_inverse!(w.G, w.m1, r, r, I_static)
        _bmm!(w.H, t, w.G)
        _jprod!(dm1, r, dr, r, dr)
        _jmul!(dm2, w.H, dm1)
        dm2 .+= dt
        _jmul!(dH, dm2, w.G)

        # Solar-source adding, SF2023-II (12), with identical half layers:
        # v₁ = r j⁺ + e j⁻,  v₂ = r(e j⁻) + j⁺,
        # j⁻new = j⁻ + H v₁, j⁺new = e j⁺ + H v₂.
        # Differentiate each product using S2014 (C.6), including de j±.
        # The attenuation factors belong to the solar source; the thermal
        # source expressions in S2014 App. C alone do not supply these terms.
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

        # S2014 (23)–(24) in the starred doubling basis:
        # r_new = r + H r t; t_new = H t.
        # dr_new = dr + dH(r t) + H(dr t + r dt),
        # dt_new = dH t + H dt; cf. (C.11)–(C.12).
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
        # Doubling the thickness squares e: d(e²) = 2 e de, using old e.
        de .*= 2 .* reshape(expk,:,1)
        expk .*= expk
    end
    # Recover physical Stokes signs and opposite-direction operators.
    # D is independent of p: d(D r D) = D dr D; S2014 (C.17)–(C.19).
    apply_D_matrix!(pol_type.n, r, t, a.r⁺⁻, a.t⁻⁻,
        da.ap_ṙ⁻⁺, da.ap_ṫ⁺⁺, da.ap_ṙ⁺⁻, da.ap_ṫ⁻⁻)
    apply_D_matrix_SFI!(pol_type.n, jm, da.ap_J̇₀⁻)
    return nothing
end

# Compute one direction of the general adding equations. The caller defers
# committing the first direction until the second has consumed all old state.
# S2014 (23)–(28), equivalently SF2023-II (12), in local argument names:
# G = (I-R Rout)⁻¹, H = Tleft G,
# ΔR = H R Tright, Tnew = H Tother, Jnew = Jbase + H(R Jin + j).
# The caller adds the original reflection to ΔR. The two directions require
# different inverses: (I-r⁻⁺ R⁺⁻)⁻¹ above, (I-R⁺⁻ r⁻⁺)⁻¹ below.
# S2014 (C.11)–(C.16) follow by repeated (C.6)–(C.7):
# dG = G(dR Rout + R dRout)G; dH = dTleft G + Tleft dG;
# substitute H=Tleft G: dH = [dTleft + H(dR Rout + R dRout)]G.
# Only dH is used downstream, so compute it without storing dG.
# dΔR = dH(R Tright) + H(dR Tright + R dTright),
# dTnew = dH Tother + H dTother,
# dJnew = dJbase + dH(R Jin+j) + H(dR Jin+R dJin+dj).
function _interaction_direction_lin!(w, R, dR, Rout, dRout, Tleft, dTleft,
                                     Tright, dTright, Tother, dTother,
                                     Jbase, dJbase, Jin, dJin, j, dj, I_static)
    _jacobian_geometric_inverse!(w.G, w.m1, R, Rout, I_static)
    _jprod!(w.dm1, R, dR, Rout, dRout)
    _bmm!(w.H, Tleft, w.G)
    _jmul!(w.dm2, w.H, w.dm1)
    w.dm2 .+= dTleft
    _jmul!(w.dH, w.dm2, w.G)
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
    # Upward output: R⁻⁺, T⁻⁻, J⁻; G = (I-r⁻⁺ R⁺⁻)⁻¹.
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

    # Downward output: R⁺⁻, T⁺⁺, J⁺; G = (I-R⁺⁻ r⁻⁺)⁻¹.
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
