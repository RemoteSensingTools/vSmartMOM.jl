# Mode B — parametric fit. This is the piece that maps GCHP onto the classic vSmartMOM
# parameter file (median radius + width + refractive index), and it is currently a STUB on
# gchp-io (`LogNormalFit()` type exists, fitting not implemented).

"""
    LogNormalFit <: AerosolRepresentation   # Mode B (NEW — implement the fit)

Fit a continuous distribution to the binned `dN/dlogD` per layer:
  - `n_modes == 1` → single `LogNormal(μ, σ)` (the classic parameter-file shape, maps to
    `vSmartMOM.Aerosol(LogNormal(μ, σ), nᵣ, nᵢ)`).
  - `n_modes > 1`  → sum of lognormals; natural place to attach per-mode speciation.

THE HARD PART (see brief §4): speciation in a parametric fit. Options:
  (i)  one mode per dominant species (multi-modal), each with its own RI;
  (ii) one size mode + a single effective RI from the bulk per-layer composition (hybrid:
       parametric size, sectional composition).
Record which you choose; start with (ii) for simplicity.
"""
struct LogNormalFit <: AerosolRepresentation
    n_modes::Int
    speciation::Symbol    # :bulk_effective_ri (ii) | :per_species_mode (i)
end
LogNormalFit(; n_modes=1, speciation=:bulk_effective_ri) = LogNormalFit(n_modes, speciation)

"""
    fit_lognormal(diameters, dNdlogD; n_modes=1) -> Vector{LogNormal}

Fit `n_modes` lognormals to a binned distribution. Suggested: weighted least squares in
log-size, or moment matching (M0/M2/M3 from the bins → μ, σ) as a robust initial guess.
"""
function fit_lognormal end   # TODO (new)

# function layer_distributions(rep::LogNormalFit, aer, ilev) ... end     # TODO (new)
