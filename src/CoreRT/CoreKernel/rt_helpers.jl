#=
Shared helper functions for the RT kernel (elemental, doubling, interaction).
Extracted from common patterns across forward, linearized, and inelastic paths.
=#

using NNlib: batched_mul as ⊠

"""
    fourier_weight(m, FT)

Azimuthal integration weight for Fourier moment `m`:
  m=0 → 0.5  (from ∫₀²π cos²(0·ϕ) dϕ / 4π)
  m>0 → 0.25 (from ∫₀²π cos²(mϕ) dϕ / 4π)
"""
@inline fourier_weight(m::Int, ::Type{FT}) where {FT} = m == 0 ? FT(0.50) : FT(0.25)

"""
    scaled_weights(m, wt_μN)

Quadrature weights scaled by the azimuthal Fourier factor:
  m=0 → wt_μN / 2
  m>0 → wt_μN / 4
"""
@inline scaled_weights(m::Int, wt_μN) = m == 0 ? wt_μN / 2 : wt_μN / 4

"""
    expdiff_neg(a, b)

Stably compute `exp(-a) - exp(-b)` for nearby positive optical-depth
arguments without losing the small difference to cancellation.
"""
@inline function expdiff_neg(a, b)
    if a == b
        return zero(a - b)
    elseif a < b
        return exp(-a) * (-expm1(-(b - a)))
    else
        return -exp(-b) * (-expm1(-(a - b)))
    end
end

"""
    rt_tol(FT, x)

Cast a scalar tolerance to the kernel element type. Use this in GPU-facing
kernels instead of Float64 literals so Float32 kernels stay Float32.
"""
@inline rt_tol(::Type{FT}, x) where {FT} = FT(x)

"""
    rt_weight_tol(FT)

Return the same-type cutoff for ignoring effectively zero quadrature weights
inside elemental kernels.
"""
@inline rt_weight_tol(::Type{FT}) where {FT} = rt_tol(FT, 1e-8)
@inline rt_weight_tol(::Type{Float32}) = 1f-8
@inline rt_weight_tol(::Type{Float64}) = 1e-8

"""
    rt_close_tol(FT)

Return the same-type near-singularity tolerance used by elemental Raman branch
checks that historically used `1e-8`.
"""
@inline rt_close_tol(::Type{FT}) where {FT} = rt_tol(FT, 1e-8)
@inline rt_close_tol(::Type{Float32}) = 1f-8
@inline rt_close_tol(::Type{Float64}) = 1e-8

"""
    rt_loose_tol(FT)

Return the same-type near-singularity tolerance used by elemental Raman branch
checks that historically used `1e-6`.
"""
@inline rt_loose_tol(::Type{FT}) where {FT} = rt_tol(FT, 1e-6)
@inline rt_loose_tol(::Type{Float32}) = 1f-6
@inline rt_loose_tol(::Type{Float64}) = 1e-6

"""
    _bmm!(C, A, B)

In-place batched matrix multiply `C = A ⊠ B` (NNlib dispatch: CUBLAS on
CUDA, threaded BLAS on CPU, generic fallback elsewhere). `C` must not
alias `A` or `B`. Phase-0 preallocation seam: the Phase-1 fused KA
kernels replace calls to this helper wholesale, so keep hot-path batched
products routed through it.
"""
@inline _bmm!(C, A, B) = NNlib.batched_mul!(C, A, B)

"""
    _interaction_scratch(added_layer) -> (mA, v1, v2)

Phase-0 scratch for the in-place interaction step: one matrix buffer and
two J-shaped buffers. Reuses the doubling workspace (`dbl_gp_refl`,
`dbl_v1`, `dbl_v2`) — dead by the time interaction runs — and falls back
to `similar` for AddedLayers built without it (e.g. the lin path).
"""
@inline function _interaction_scratch(added_layer)
    mA = added_layer.dbl_gp_refl === nothing ? similar(added_layer.t⁺⁺) : added_layer.dbl_gp_refl
    v1 = added_layer.dbl_v1 === nothing ? similar(added_layer.j₀⁺) : added_layer.dbl_v1
    v2 = added_layer.dbl_v2 === nothing ? similar(added_layer.j₀⁺) : added_layer.dbl_v2
    return mA, v1, v2
end

"""
    compute_geometric_progression!(gp_refl, tt_gp, r⁻⁺, t⁺⁺, I_static, temp2, temp1_ptr, temp2_ptr)

Compute the matrix geometric-series factor `(E − R·R)⁻¹` that captures all
internal reflections within a homogeneous layer being doubled, then
pre-multiply by `T⁺⁺` to form the helper `tt_gp = T·(E − R·R)⁻¹` reused by
[`doubling_source_update!`](@ref) and [`doubling_rt_update!`](@ref).

This is the inner factor of Sanghavi et al. (2014, JQSRT 133:412–433),
Eqs. (23)–(28) — the matrix-operator-method *adding equations*. For a layer
combined with an identical copy of itself, the only place an inverse matrix
appears is `(E − R·R)⁻¹`; representing it as the geometric series
`I + R·R + (R·R)² + …` makes clear that this term sums up an infinite series
of internal reflections.

The batched matrix inversion `batch_inv!` is dispatched by array type to
threaded BLAS (CPU), CUBLAS `getrf_strided_batched! + getri_strided_batched!`
(CUDA), or the portable KernelAbstractions LU kernel (Metal). One batched
call covers all spectral points; see [Concepts/07](../../docs/src/pages/concepts/07_architecture.md).

Mutates `gp_refl` and `tt_gp` in place.
"""

@inline function compute_geometric_progression!(gp_refl, tt_gp, r⁻⁺, t⁺⁺, I_static, temp2, temp1_ptr, temp2_ptr)
    if _use_fused_gp(tt_gp)
        # Fused single-kernel path: build (E − R·R), LU-factorise it, and
        # right-solve `tt_gp · (E − R·R) = T` directly — no explicit inverse and
        # no separate batched GEMM. For the small RT matrices this is several×
        # faster than the inverse-then-multiply path below (see
        # test/benchmarks/batched_fused_v2_benchmark.jl). The fallback runs when
        # the 3N² local-memory tile exceeds the device budget (large N) or on CPU.
        ka_fused_gp_solve!(tt_gp, r⁻⁺, t⁺⁺, KernelAbstractions.get_backend(tt_gp))
    else
        # In-place (Phase 0): tt_gp doubles as scratch for R·R before its
        # final assignment — no allocation, no extra buffer.
        _bmm!(tt_gp, r⁻⁺, r⁻⁺)                           # R·R (scratch use)
        temp2 .= I_static .- tt_gp                       # (E − R·R)
        batch_inv!(gp_refl, temp2, temp1_ptr, temp2_ptr) # (E − R·R)⁻¹
        _bmm!(tt_gp, t⁺⁺, gp_refl)                       # T · (E − R·R)⁻¹
    end
    return nothing
end

"""
    doubling_source_update!(j₀⁺, j₀⁻, j₁⁺, j₁⁻, r⁻⁺, tt_gp, expk)

Update the source-function vectors during a single doubling step.

Applies the adding equations for `J₂₀⁻` and `J₀₂⁺` from Sanghavi et al. (2014),
Eqs. (27)–(28), restated for two identical layers as Eqs. (8) of Sanghavi &
Frankenberg (2023):

    j₁⁺  =  j₀⁺ · expk         (direct-beam attenuation across the lower copy)
    j₁⁻  =  j₀⁻ · expk
    j₀⁻  ←  j₀⁻ + tt_gp · (j₁⁻ + r⁻⁺ · j₀⁺)
    j₀⁺  ←  j₁⁺ + tt_gp · (j₀⁺ + r⁻⁺ · j₁⁻)

`expk` is the scalar attenuation factor `exp(−δτ/μ₀)` for the elemental
optical thickness; it is *squared* by [`doubling_rt_update!`](@ref) at the
end of each iteration so the layer thickness doubles. `tt_gp` is the helper
`T · (E − R·R)⁻¹` produced by [`compute_geometric_progression!`](@ref).

Common to the forward, linearized, and inelastic paths.
"""
@inline function doubling_source_update!(j₀⁺, j₀⁻, j₁⁺, j₁⁻, r⁻⁺, tt_gp, expk, v1, v2)
    @inbounds @views j₁⁺[:,1,:] .= j₀⁺[:,1,:] .* expk'
    @inbounds @views j₁⁻[:,1,:] .= j₀⁻[:,1,:] .* expk'
    # In-place (Phase 0) restatement of
    #   j₀⁻ .= j₀⁻ .+ (tt_gp ⊠ (j₁⁻ .+ r⁻⁺ ⊠ j₀⁺))
    #   j₀⁺ .= j₁⁺ .+ (tt_gp ⊠ (j₀⁺ .+ r⁻⁺ ⊠ j₁⁻))
    # Identical operation order (each ⊠ replaced by the same batched GEMM,
    # each dotted chain by the same fused broadcast), so results are
    # bit-identical; only the ⊠ output allocations are gone. Note the j₀⁺
    # update reads the PRE-update j₀⁺ — v1 snapshots the RHS before the
    # assignment, exactly as the allocating form did.
    _bmm!(v1, r⁻⁺, j₀⁺);  v1 .= j₁⁻ .+ v1
    _bmm!(v2, tt_gp, v1); j₀⁻ .= j₀⁻ .+ v2
    _bmm!(v1, r⁻⁺, j₁⁻);  v1 .= j₀⁺ .+ v1
    _bmm!(v2, tt_gp, v1); j₀⁺ .= j₁⁺ .+ v2
    return nothing
end

# Back-compat allocating wrapper (non-hot callers / external code).
@inline doubling_source_update!(j₀⁺, j₀⁻, j₁⁺, j₁⁻, r⁻⁺, tt_gp, expk) =
    doubling_source_update!(j₀⁺, j₀⁻, j₁⁺, j₁⁻, r⁻⁺, tt_gp, expk,
                            similar(j₀⁺), similar(j₀⁺))

"""
    doubling_rt_update!(r⁻⁺, t⁺⁺, tt_gp, expk)

Update the reflection and transmission supermatrices during a single
doubling step.

Applies the adding equations from Sanghavi et al. (2014), Eqs. (23)–(24),
restated for a homogeneous layer combined with an identical copy of itself:

    R₂₀  =  R₁₀ + T₀₁ · (E − R₂₁·R₀₁)⁻¹ · R₂₁ · T₁₀
    T₂₀  =  T₂₁ · (E − R₀₁·R₂₁)⁻¹ · T₁₀

In the symmetric case (both layers identical), the geometric-series helper
`tt_gp = T · (E − R·R)⁻¹` from [`compute_geometric_progression!`](@ref)
collapses both updates to:

    r⁻⁺  ←  r⁻⁺  +  tt_gp · r⁻⁺ · t⁺⁺
    t⁺⁺  ←  tt_gp · t⁺⁺

`expk` (= `exp(−δτ/μ₀)`) is squared so the elemental thickness doubles
between iterations: `n` doublings give layer thickness `2ⁿ · δτ`. That's
the **logarithmic-in-τ** scaling that makes MOM cheap for thick atmospheres.
"""
@inline function doubling_rt_update!(r⁻⁺, t⁺⁺, tt_gp, expk, m1, m2)
    # In-place (Phase 0) restatement of
    #   r⁻⁺ .= r⁻⁺ .+ (tt_gp ⊠ r⁻⁺ ⊠ t⁺⁺)   (⊠ is left-associative)
    #   t⁺⁺ .= tt_gp ⊠ t⁺⁺
    # Same GEMMs in the same order — bit-identical, allocation-free.
    _bmm!(m1, tt_gp, r⁻⁺)
    _bmm!(m2, m1, t⁺⁺);  r⁻⁺ .= r⁻⁺ .+ m2
    _bmm!(m1, tt_gp, t⁺⁺); t⁺⁺ .= m1
    expk .= expk .^ 2
    return nothing
end

# Back-compat allocating wrapper (non-hot callers / external code).
@inline doubling_rt_update!(r⁻⁺, t⁺⁺, tt_gp, expk) =
    doubling_rt_update!(r⁻⁺, t⁺⁺, tt_gp, expk, similar(t⁺⁺), similar(t⁺⁺))

"""
    zero_added_noscat!(added_layer, τ_λ, qp_μN)

Zero out reflectance/source fields and set Beer-law transmission for a
non-scattering layer.  Replaces the repeated pattern in `rt_kernel!` methods.
"""
@inline function zero_added_noscat!(added_layer, τ_λ, qp_μN)
    added_layer.r⁻⁺[:] .= 0
    added_layer.r⁺⁻[:] .= 0
    added_layer.j₀⁻[:] .= 0
    _set_transmission_noscat!(added_layer.t⁺⁺, added_layer.t⁻⁻, τ_λ, qp_μN)
    return nothing
end

"""
    zero_added_noscat_ie!(added_layer, τ_λ, qp_μN)

Like `zero_added_noscat!` but also zeros the inelastic fields
(`ier⁻⁺`, `ier⁺⁻`, `ieJ₀⁻`, `iet⁻⁻`, `iet⁺⁺`, `ieJ₀⁺`).
"""
@inline function zero_added_noscat_ie!(added_layer, τ_λ, qp_μN)
    zero_added_noscat!(added_layer, τ_λ, qp_μN)
    added_layer.ier⁻⁺[:] .= 0
    added_layer.ier⁺⁻[:] .= 0
    added_layer.ieJ₀⁻[:] .= 0
    added_layer.iet⁻⁻[:] .= 0
    added_layer.iet⁺⁺[:] .= 0
    added_layer.ieJ₀⁺[:] .= 0
    return nothing
end

"""
    copy_added_to_composite!(composite_layer, added_layer)

Copy all fields from an `AddedLayer` into the `CompositeLayer` (for TOA, iz==1).
"""
@inline function copy_added_to_composite!(composite_layer, added_layer)
    composite_layer.T⁺⁺[:], composite_layer.T⁻⁻[:] = (added_layer.t⁺⁺, added_layer.t⁻⁻)
    composite_layer.R⁻⁺[:], composite_layer.R⁺⁻[:] = (added_layer.r⁻⁺, added_layer.r⁺⁻)
    composite_layer.J₀⁺[:], composite_layer.J₀⁻[:] = (added_layer.j₀⁺, added_layer.j₀⁻)
    # v0.7 Phase A.2a — copy per-source slots too (if any). AddedLayerRS lacks
    # these slots; Raman surface-emission sources are injected directly into
    # the layer's legacy j₀ arrays before this copy.
    if hasproperty(added_layer, :j₀_by_src)
        for (key, slot) in pairs(added_layer.j₀_by_src)
            cslot = composite_layer.J₀_by_src[key]
            cslot.J₀⁺[:] = slot.j₀⁺
            cslot.J₀⁻[:] = slot.j₀⁻
        end
    end
    return nothing
end

"""
    copy_added_to_composite_ie!(composite_layer, added_layer)

Like `copy_added_to_composite!` but also copies inelastic fields.
"""
@inline function copy_added_to_composite_ie!(composite_layer, added_layer)
    copy_added_to_composite!(composite_layer, added_layer)
    composite_layer.ieT⁺⁺[:], composite_layer.ieT⁻⁻[:] = (added_layer.iet⁺⁺, added_layer.iet⁻⁻)
    composite_layer.ieR⁻⁺[:], composite_layer.ieR⁺⁻[:] = (added_layer.ier⁻⁺, added_layer.ier⁺⁻)
    composite_layer.ieJ₀⁺[:], composite_layer.ieJ₀⁻[:] = (added_layer.ieJ₀⁺, added_layer.ieJ₀⁻)
    return nothing
end
