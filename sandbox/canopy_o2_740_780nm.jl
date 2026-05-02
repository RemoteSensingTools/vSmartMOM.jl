# =================================================================
# Sample run: homogeneous canopy + within-canopy O2 absorption,
# isolated (no real atmosphere above), with PROSPECT-derived
# spectral leaf optics.
#
#   - canopy depth ≈ 30 m  (canopy_dp = 3 hPa, O2 mixed in at 21%)
#   - LAI           = 4.0  (4 sub-layers, leaves interleaved with O2)
#   - LAD           = spherical leaves (G(mu) = 0.5)
#   - typical green-leaf chlorophyll (Ccab = 40 μg/cm²)
#   - 740–780 nm at 0.2 cm⁻¹ resolution (full O2 A-band core).  The
#     v5.2 O2 ABSCO LUT only spans ≈ 754.6–784.6 nm; outside that
#     range σ_O2 is extrapolated to zero (correct physically — A-band
#     is the only O2 feature in this window).
#
# The YAML config sets the column above the canopy to a vanishingly
# thin layer (Δp = 0.01 hPa) so the TOA reflectance matrix returned
# by rt_run is essentially the canopy's own R⁻⁺ — the O2 imprint you
# see comes only from the in-canopy 3 hPa column.
# =================================================================

using vSmartMOM
using vSmartMOM.CoreRT
using CanopyOptics
using CairoMakie

# Run with `julia --project=docs sandbox/canopy_o2_740_780nm.jl`.
if !haskey(ENV, "VSMARTMOM_O2_ABSCO")
    default_o2_lut = "/net/squid/data1/ABSCO_CS_Database/v5.2_final/o2_v52_v2.jld2"
    isfile(default_o2_lut) || error(
        "Set ENV[\"VSMARTMOM_O2_ABSCO\"] to the O2 ABSCO JLD2 LUT path.")
    ENV["VSMARTMOM_O2_ABSCO"] = default_o2_lut
end

# ── 1) Atmosphere + O2 absorption + spectral grid (from YAML) ────────
yaml_path = joinpath(@__DIR__, "canopy_o2_740_780nm.yaml")
params    = parameters_from_yaml(yaml_path)

# ── 2) Build a typical green-leaf PROSPECT spectrum ──────────────────
leaf = CanopyOptics.LeafProspectProProperties(
    N      = 1.5,    # leaf structure
    Ccab   = 40.0,   # chlorophyll a+b  [μg/cm²]
    Ccar   = 8.0,    # carotenoids      [μg/cm²]
    Canth  = 0.0,
    Cbrown = 0.0,
    Cw     = 0.012,  # equivalent water thickness [cm]
    Cm     = 0.009,  # dry matter content        [g/cm²]
    Cprot  = 0.0,
    Ccbc   = 0.0,
)

# Spectral leaf R/T on a 1 nm grid; CanopySurface interpolates onto
# the RT wavenumber grid internally.
LAD = CanopyOptics.spherical_leaves(Float64)
canopy = CanopySurface_from_prospect(
    leaf, 400.0:1.0:2500.0;
    soil        = LambertianSurfaceScalar(0.10),
    LAI         = 4.0,
    n_layers    = 4,        # need ≥ 2 sublayers to interleave atm
    LAD         = LAD,
    include_atm = true,
    canopy_dp   = 3.0,      # ~30 m forest ≈ 3 hPa
)

# Replace the per-band soil BRDF with the canopy (with O2 inside).
for ib in eachindex(params.brdf)
    params.brdf[ib] = canopy
end

# ── 3) Build model and run forward RT ────────────────────────────────
model = model_from_parameters(params)
R, T  = rt_run(model)
println("R size (vza, stokes, λ): ", size(R))

# ── 4) Pull nadir Stokes I and plot vs wavelength ────────────────────
νs   = collect(model.atmosphere.spec_bands[1])     # cm⁻¹
λ_nm = 1e7 ./ νs

vza_idx = findfirst(==(0.0), model.obs_geom.vza)
@assert vza_idx !== nothing "Nadir VZA=0 must be in the geometry."

R_nadir = R[vza_idx, 1, :]
T_nadir = T[vza_idx, 1, :]
println("R(I, nadir) min/max: ", extrema(R_nadir))
println("T(I, nadir) min/max: ", extrema(T_nadir))

fig = Figure(size = (900, 760))
ax_R = Axis(fig[1, 1];
    ylabel = "Nadir TOA reflectance  (Stokes I)",
    title  = "Spherical-LAD canopy (LAI=4, Ccab=40 μg/cm²) with within-canopy O₂ — 740–780 nm",
)
lines!(ax_R, λ_nm, R_nadir; color = :darkgreen, linewidth = 1.0)
xlims!(ax_R, 740, 780)
hidexdecorations!(ax_R; grid = false)

ax_T = Axis(fig[2, 1];
    xlabel = "Wavelength (nm)",
    ylabel = "Nadir BOA transmittance  (Stokes I)",
)
lines!(ax_T, λ_nm, T_nadir; color = :darkorange, linewidth = 1.0)
xlims!(ax_T, 740, 780)

linkxaxes!(ax_R, ax_T)

out_png = joinpath(@__DIR__, "canopy_o2_740_780nm.png")
save(out_png, fig)
println("Saved nadir R and T plot → ", out_png)
fig
