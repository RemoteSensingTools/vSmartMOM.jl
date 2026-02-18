#=
Shared helper functions for the RT kernel (elemental, doubling, interaction).
=#

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
