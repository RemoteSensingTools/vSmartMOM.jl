# =============================================================================
# sif_raman.jl — SIF retrieval demo via Legendre-polynomial fit (Phase 6 port
# from sanghavi).
#
# Reads grid spectra produced by `creategrid_O2Aband_RamanSIF.jl` (or similar)
# for each combination of (SZA, surface albedo, surface pressure) and fits a
# simple Legendre-polynomial × solar-reference + constant-SIF model to the
# 758-759.2 nm window. Returns the relative SIF (SIF / max-radiance) as a
# (nρ × nSZA) matrix, with and without Raman scattering in the forward model.
#
# The original sanghavi script was a notebook-style post-processor with
# hard-coded `/home/sanghavi/RamanSIFgrid/` paths and inline Plots calls. This
# port:
#   1. Wraps the fitting math in the reusable `fit_sif_legendre(ν, rad, ref;
#      fit_window, npoly)` function.
#   2. Wraps the grid scan in `run_sif_grid(griddir; ...)` which checks for
#      file existence (non-public data path in the original).
#   3. Drops Plots — the fitted spectrum + residuals are returned to the
#      caller, who can choose to plot.
# =============================================================================

using LegendrePolynomials
using DelimitedFiles

"""
    fit_sif_legendre(wl, rad, solar_ref; fit_window, npoly = 3)

Fit a Legendre-polynomial spectral shape plus a constant SIF offset to
`rad`, using `solar_ref` as the spectral reference.

Returns a NamedTuple with
  `fit`      — fitted spectrum over the window indices
  `SIF`      — absolute SIF offset (last coefficient in the fit)
  `rel_SIF`  — SIF / maximum(rad_in_window)
  `ind`      — indices of the fit window into `wl`
  `poly`     — Legendre polynomial basis × solar_ref
"""
function fit_sif_legendre(wl::AbstractVector, rad::AbstractVector,
                          solar_ref::AbstractVector;
                          fit_window::Tuple{<:Real, <:Real},
                          npoly::Integer = 3)
    ind = findall(x -> x > fit_window[1] && x < fit_window[2], wl)
    isempty(ind) && return (fit = Float64[], SIF = 0.0, rel_SIF = 0.0,
                            ind = ind, poly = zeros(0, 0))
    iLeg = range(-1, 1, length(ind))
    poly_basis = Pl.(iLeg, (0:npoly)')
    Kref = solar_ref[ind] .* poly_basis
    K    = [Kref ones(length(ind))]
    x    = K \ rad[ind]
    fit  = K * x
    SIF  = x[end]
    rel_SIF = SIF / maximum(rad[ind])
    return (; fit, SIF, rel_SIF, ind, poly = Kref)
end

"""
    run_sif_grid(griddir; fit_window=(758.0, 759.2), sza_list, ρ_list, psurf_list)

Scan a RamanSIF grid directory produced by `creategrid_O2Aband_RamanSIF.jl`
and return `(xSIF, xSIF_noRS)` where each is a (nρ × nSZA) matrix of
relative SIF per (surface-pressure, albedo, SZA) scenario.

File names follow sanghavi convention
  `rayl(SIF)_sza{SZA}_alb{ALB}_psurf{PSURF}hpa_{nors|rrs}_ABO2.dat`

Returns `nothing` if `griddir` is not a directory — the function is
idempotent when the data bundle is not present on disk.
"""
function run_sif_grid(griddir::AbstractString;
                      fit_window::Tuple{<:Real,<:Real} = (758.0, 759.2),
                      sza_list::AbstractVector{<:Integer},
                      ρ_list::AbstractVector{<:Real},
                      psurf_list::AbstractVector{<:Integer} = [1000, 750, 500],
                      Tfit::Type = Float64)
    isdir(griddir) || return nothing

    ρ_str = [replace(string(round(r, digits = 2)), "." => "p") for r in ρ_list]

    xSIF      = zeros(Tfit, length(ρ_list), length(sza_list))
    xSIF_noRS = zeros(Tfit, length(ρ_list), length(sza_list))

    for isurf in eachindex(psurf_list)
        for iρ in eachindex(ρ_list)
            for iA in eachindex(sza_list)
                sza_s   = string(sza_list[iA])
                alb_s   = ρ_str[iρ]
                psurf_s = string(psurf_list[isurf])

                f_rrs   = joinpath(griddir, "raylSIF_sza$(sza_s)_alb$(alb_s)_psurf$(psurf_s)hpa_rrs_ABO2.dat")
                f_nors  = joinpath(griddir, "raylSIF_sza$(sza_s)_alb$(alb_s)_psurf$(psurf_s)hpa_nors_ABO2.dat")
                (isfile(f_rrs) && isfile(f_nors)) || continue

                specRRS  = readdlm(f_rrs)
                specNoRS = readdlm(f_nors)

                wl       = 1e7 ./ specNoRS[:, 1]
                rad_rrs  = specRRS[:, 2] .+ specRRS[:, 5]   # elastic + inelastic
                rad_nors = specNoRS[:, 2]
                ref      = specNoRS[:, 5]                    # solar reference

                fit_rrs  = fit_sif_legendre(wl, rad_rrs, ref; fit_window)
                fit_nors = fit_sif_legendre(wl, rad_nors, ref; fit_window)

                xSIF[iρ, iA]      = fit_rrs.rel_SIF
                xSIF_noRS[iρ, iA] = fit_nors.rel_SIF
            end
        end
    end

    return (; xSIF, xSIF_noRS, sza_list, ρ_list, psurf_list, fit_window)
end

# Optional entry point: when run as a script against a known grid dir.
if abspath(PROGRAM_FILE) == @__FILE__
    griddir = get(ENV, "RAMANSIF_GRID_DIR", "/home/sanghavi/RamanSIFgrid")
    sza_acos = [acosd(i / 20) for i in 0:20]
    sza      = reverse(Int.(ceil.(sza_acos[8:21])))
    ρ        = collect(0:14) .* 0.05
    out = run_sif_grid(griddir; sza_list = sza, ρ_list = ρ)
    if isnothing(out)
        @info "sif_raman: no data at $(griddir) — skipping grid scan" env="RAMANSIF_GRID_DIR"
    else
        @info "sif_raman grid complete" xSIF_extrema = extrema(out.xSIF)
    end
end
