# Mode A — sectional (bins as they are).

"""
    Sectional <: AerosolRepresentation

Use the 15 TOMAS bins directly. In-bin integration policy:
  - `:direct`   — evaluate Mie at the bin center (fast; diagnostic).
  - `:constant` — Gauss–Legendre in log-size, constant dN/dlogD per bin.
  - `:linear`   — GL with piecewise-linear (monotone-limited) dN/dlogD.

PORT TARGET: the `{DirectBinSum, ConstantIntegrationPerBin, LinearIntegrationPerBin}`
policies and `_integrated_bin_extinction` in gchp-io:src/Aerosols/aod_diagnostics.jl.
"""
struct Sectional <: AerosolRepresentation
    integration::Symbol   # :direct | :constant | :linear
    nquad::Int            # GL nodes for :constant / :linear (gchp-io uses 160)
end
Sectional(; integration=:constant, nquad=160) = Sectional(integration, nquad)

"""
    SplineSectional <: AerosolRepresentation   # Mode A2 (NEW)

Fit a smooth spline through the 15 bin values of `dN/dlogD` (per layer), then quadrature on
the smooth curve instead of the bin staircase. Reduces discretization artifacts while
staying non-parametric. No reference impl yet — this is new work.
"""
struct SplineSectional <: AerosolRepresentation
    nquad::Int
    boundary::Symbol      # spline end condition, e.g. :natural | :flat
end
SplineSectional(; nquad=160, boundary=:natural) = SplineSectional(nquad, boundary)

# function layer_distributions(rep::Sectional, aer, ilev) ... end        # TODO (port)
# function layer_distributions(rep::SplineSectional, aer, ilev) ... end  # TODO (new)
