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
    compute_geometric_progression!(gp_refl, tt_gp, r⁻⁺, t⁺⁺, I_static, temp1_ptr, temp2_ptr)

Compute the geometric-progression factor `(I - R·R)⁻¹` and pre-multiply by `T⁺⁺`.
Used in doubling and interaction steps across all RT paths.

Mutates `gp_refl` (temp1) and `tt_gp` in place.
"""
@inline function compute_geometric_progression!(gp_refl, tt_gp, r⁻⁺, t⁺⁺, I_static, temp2, temp1_ptr, temp2_ptr)
    temp2 .= I_static .- r⁻⁺ ⊠ r⁻⁺
    batch_inv!(gp_refl, temp2, temp1_ptr, temp2_ptr)
    tt_gp .= t⁺⁺ ⊠ gp_refl
    return nothing
end

"""
    doubling_source_update!(j₀⁺, j₀⁻, j₁⁺, j₁⁻, r⁻⁺, tt_gp, expk)

Update source functions during a single doubling step.
Applies the adding equations for `J⁻₀₂` and `J⁺₂₀` (Eqs.8 in Raman paper).
Common to forward, linearized, and inelastic paths.
"""
@inline function doubling_source_update!(j₀⁺, j₀⁻, j₁⁺, j₁⁻, r⁻⁺, tt_gp, expk)
    @inbounds @views j₁⁺[:,1,:] .= j₀⁺[:,1,:] .* expk'
    @inbounds @views j₁⁻[:,1,:] .= j₀⁻[:,1,:] .* expk'
    j₀⁻ .= j₀⁻ .+ (tt_gp ⊠ (j₁⁻ .+ r⁻⁺ ⊠ j₀⁺))
    j₀⁺ .= j₁⁺ .+ (tt_gp ⊠ (j₀⁺ .+ r⁻⁺ ⊠ j₁⁻))
    return nothing
end

"""
    doubling_rt_update!(r⁻⁺, t⁺⁺, tt_gp, expk)

Update reflection and transmission matrices during a single doubling step.
Applies `R₂₀ = R₁₀ + T₀₁·(I-R₂₁R₀₁)⁻¹·R₂₁·T₁₀` and
       `T₂₀ = T₂₁·(I-R₀₁R₂₁)⁻¹·T₁₀`.
"""
@inline function doubling_rt_update!(r⁻⁺, t⁺⁺, tt_gp, expk)
    r⁻⁺ .= r⁻⁺ .+ (tt_gp ⊠ r⁻⁺ ⊠ t⁺⁺)
    t⁺⁺ .= tt_gp ⊠ t⁺⁺
    expk .= expk .^ 2
    return nothing
end

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
