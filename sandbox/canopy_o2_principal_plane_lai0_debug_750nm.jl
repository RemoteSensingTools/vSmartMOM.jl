# =================================================================
# LAI=0 debug principal-plane reflectance slice at 750 nm
#
# Same plot style as sandbox/canopy_o2_principal_plane_750nm.jl, but with
# the one-layer planophile canopy set to LAI=0. This is a diagnostic for the
# zero-canopy limit over the same Lambertian soil.
#
# Signed-VZA convention follows docs/make.jl:
#   negative VZA = VAZ 180 deg, backscatter / hotspot side
#   positive VZA = VAZ   0 deg, forward-scatter side
# =================================================================

using vSmartMOM
using vSmartMOM.CoreRT
using CanopyOptics
using CairoMakie

# Run with `julia --project=docs sandbox/canopy_o2_principal_plane_lai0_debug_750nm.jl`.
if !haskey(ENV, "VSMARTMOM_O2_ABSCO")
    default_o2_lut = "/net/squid/data1/ABSCO_CS_Database/v5.2_final/o2_v52_v2.jld2"
    isfile(default_o2_lut) || error(
        "Set ENV[\"VSMARTMOM_O2_ABSCO\"] to the O2 ABSCO JLD2 LUT path.")
    ENV["VSMARTMOM_O2_ABSCO"] = default_o2_lut
end

# -- 1) Base config and single wavelength --------------------------
yaml_path = joinpath(@__DIR__, "canopy_o2_740_780nm.yaml")
params = parameters_from_yaml(yaml_path)

lambda_target = 750.0
nu_target = 1e7 / lambda_target
params.spec_bands = [Float64[nu_target]]

# -- 2) Principal-plane geometry -----------------------------------
signed_vza = collect(-80.0:5.0:80.0)
params.vza = abs.(signed_vza)
params.vaz = ifelse.(signed_vza .< 0, 180.0, 0.0)

params.max_m = 8
params.l_trunc = 64

# -- 3) One-layer planophile canopy, zero LAI -----------------------
leaf = CanopyOptics.LeafProspectProProperties(
    N = 1.5,
    Ccab = 40.0,
    Ccar = 8.0,
    Canth = 0.0,
    Cbrown = 0.0,
    Cw = 0.012,
    Cm = 0.009,
    Cprot = 0.0,
    Ccbc = 0.0,
)
LAD = CanopyOptics.planophile_leaves2(Float64)
canopy = CanopySurface_from_prospect(
    leaf, 400.0:1.0:2500.0;
    soil = LambertianSurfaceScalar(0.10),
    LAI = 0.0,
    n_layers = 1,
    LAD = LAD,
)
for ib in eachindex(params.brdf)
    params.brdf[ib] = canopy
end

# -- 4) Run RT ------------------------------------------------------
model = model_from_parameters(params)
R, _ = rt_run(model)
R_principal = R[:, 1, 1]
nadir_idx = findfirst(==(0.0), signed_vza)
sec_vza_scaled = R_principal[nadir_idx] ./ cosd.(abs.(signed_vza))

println("R size (dirs, stokes, lambda): ", size(R))
println("LAI=0 principal-plane R(I) min/max: ", extrema(R_principal))
for x in (-params.sza, 0.0, params.sza)
    idx = findfirst(==(x), signed_vza)
    idx === nothing || println("R(I), signed VZA = ", x, " deg: ", R_principal[idx])
end

# -- 5) Plot --------------------------------------------------------
fig = Figure(size = (820, 520))
ax = Axis(fig[1, 1];
    xlabel = "Signed VZA in solar principal plane (deg)",
    ylabel = "TOA reflectance (Stokes I)",
    title = "LAI=0 debug: planophile one-layer canopy at 750 nm (SZA = 30 deg)",
)
lines!(ax, signed_vza, R_principal; color = :darkgreen, linewidth = 2.5,
       label = "LAI=0 reflectance")
scatter!(ax, signed_vza, R_principal; color = :darkgreen, markersize = 7)
lines!(ax, signed_vza, sec_vza_scaled; color = :black, linestyle = :dash,
       linewidth = 2.0, label = "R(nadir) / cos(|VZA|)")

vlines!(ax, [-params.sza]; color = :red, linestyle = :dash, linewidth = 1.8,
        label = "Backscatter side")
vlines!(ax, [0.0]; color = :gray45, linestyle = :dot, linewidth = 1.4,
        label = "Nadir")

text!(ax, -78, minimum(R_principal);
      text = "negative: VAZ = 180 deg\npositive: VAZ = 0 deg",
      align = (:left, :bottom), fontsize = 13, color = :gray30)

axislegend(ax; position = :lt, framevisible = false)
xlims!(ax, -82, 82)

out_png = joinpath(@__DIR__, "canopy_o2_principal_plane_lai0_debug_750nm.png")
save(out_png, fig)
println("Saved LAI=0 debug principal-plane plot -> ", out_png)
fig
