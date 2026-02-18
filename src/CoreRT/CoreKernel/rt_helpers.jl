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
