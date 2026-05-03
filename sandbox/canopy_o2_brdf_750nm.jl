# =================================================================
# Canopy BRDF phase plot at 750 nm
#
# Uses the same leaf/soil setup as the spectral
# script (sandbox/canopy_o2_740_780nm.jl) but reduces to a single
# wavelength (750 nm — outside the O2 A-band, so we see the canopy
# scattering structure cleanly) and sweeps a dense (VZA, VAZ) grid
# at fixed SZA = 30°.
#
# More Fourier moments are needed than for nadir-only spectra to
# resolve the azimuthal BRDF structure (hot spot, principal-plane
# asymmetry); we use max_m = 8 here.
#
# The canopy uses one layer and a planophile leaf-angle distribution:
# LAD = CanopyOptics.planophile_leaves2(Float64).
# =================================================================

using vSmartMOM
using vSmartMOM.CoreRT
using CanopyOptics
using CairoMakie

# Run with `julia --project=docs sandbox/canopy_o2_brdf_750nm.jl`.
if !haskey(ENV, "VSMARTMOM_O2_ABSCO")
    default_o2_lut = "/net/squid/data1/ABSCO_CS_Database/v5.2_final/o2_v52_v2.jld2"
    isfile(default_o2_lut) || error(
        "Set ENV[\"VSMARTMOM_O2_ABSCO\"] to the O2 ABSCO JLD2 LUT path.")
    ENV["VSMARTMOM_O2_ABSCO"] = default_o2_lut
end

# ── 1) Base config (atmosphere + O2 absorption + canopy floor) ──────
yaml_path = joinpath(@__DIR__, "canopy_o2_740_780nm.yaml")
params    = parameters_from_yaml(yaml_path)

# ── 2) Single wavelength at 750 nm ──────────────────────────────────
λ_target = 750.0                     # nm
ν_target = 1e7 / λ_target            # ≈ 13333.33 cm⁻¹  (just outside the O2 A-band)
params.spec_bands = [Float64[ν_target]]

# ── 3) Dense angular grid: SZA = 30°, VZA × VAZ ─────────────────────
vzas_deg = collect(0.0:5.0:80.0)     # 17 polar zenith angles
vazs_deg = collect(0.0:5.0:355.0)    # 72 relative azimuth angles (180 = backscatter)
n_vza, n_vaz = length(vzas_deg), length(vazs_deg)

vza_all = repeat(vzas_deg; inner=n_vaz)
vaz_all = repeat(vazs_deg; outer=n_vza)
params.vza = vza_all
params.vaz = vaz_all

# ── 4) Fourier moments + streams ────────────────────────────────────
# Canopy uses compute_Z_matrices_aniso (no δ-M truncation) so l_trunc here
# only buys quadrature streams via rt_set_streams. BiLambertianCanopyScattering
# is azimuthally smooth, so the Fourier series in φ converges by m ≈ 8.
params.max_m   = 8
params.l_trunc = 64

# ── 5) One-layer PROSPECT canopy with planophile LAD ───────────────
leaf = CanopyOptics.LeafProspectProProperties(
    N      = 1.5,
    Ccab   = 40.0,
    Ccar   = 8.0,
    Canth  = 0.0,
    Cbrown = 0.0,
    Cw     = 0.012,
    Cm     = 0.009,
    Cprot  = 0.0,
    Ccbc   = 0.0,
)
LAD = CanopyOptics.planophile_leaves2(Float64)
canopy = CanopySurface_from_prospect(
    leaf, 400.0:1.0:2500.0;
    soil        = LambertianSurfaceScalar(0.10),
    LAI         = 4.0,
    n_layers    = 1,
    LAD         = LAD,
)
for ib in eachindex(params.brdf)
    params.brdf[ib] = canopy
end

# ── 6) Build model and run RT at the (VZA, VAZ) grid ────────────────
model = model_from_parameters(params)
R, _  = rt_run(model)
println("R size (dirs, stokes, λ): ", size(R))

# Unpack into (vza, vaz) grid following the iteration order above.
R_grid = permutedims(reshape(R[:, 1, 1], n_vaz, n_vza))
println("R_grid extrema: ", extrema(R_grid))

# ── 7) Polar phase plot ─────────────────────────────────────────────
fig = Figure(size = (820, 720))
ax  = PolarAxis(fig[1, 1];
    title       = "Planophile-LAD one-layer canopy BRDF at 750 nm  (SZA = 30°, LAI = 4, Ccab = 40)",
    theta_0     = π / 2,    # 0° azimuth at the top
    direction   = -1,        # azimuth increases clockwise (visual convention)
    rticks      = 0:20:80,
    rticklabelsize = 11,
    thetaticks  = (deg2rad.(0:30:330), string.(0:30:330) .* "°"),
)

θ_rad = deg2rad.(vcat(vazs_deg, 360.0))
R_plot = hcat(R_grid, R_grid[:, 1])
hm = surface!(ax, θ_rad, vzas_deg, R_plot';
              colormap = :viridis, shading = NoShading)

Colorbar(fig[1, 2], hm; label = "Reflectance  (Stokes I)")

# Mark the solar-opposition direction. In vSmartMOM's principal-plane
# convention, the hotspot/backscatter side is at VAZ=180 and VZA=SZA.
sca_θ = deg2rad(180.0)
sca_r = 30.0
scatter!(ax, [sca_θ], [sca_r];
         color = :red, markersize = 12, marker = :star5,
         strokecolor = :white, strokewidth = 1.0)

out_png = joinpath(@__DIR__, "canopy_o2_brdf_750nm.png")
save(out_png, fig)
println("Saved BRDF plot → ", out_png)
fig
