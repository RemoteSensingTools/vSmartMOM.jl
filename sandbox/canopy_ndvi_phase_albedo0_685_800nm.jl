# =================================================================
# Full phase-space canopy NDVI BRDF effect at 685 and 800 nm
#
# This is a no-gas, no-LUT sandbox. It runs 685 and 800 nm together over a
# dense (VZA, VAZ) grid, computes NDVI, and plots the angular structure for a
# black Lambertian soil (albedo = 0).
#
# Run with:
#   julia --project=docs sandbox/canopy_ndvi_phase_albedo0_685_800nm.jl
# =================================================================

using vSmartMOM
using vSmartMOM.CoreRT
using CanopyOptics
using CairoMakie
using Statistics

const λ_nm = [685.0, 800.0]
const ν_cm⁻¹ = 1e7 ./ λ_nm
const sza_deg = 30.0

vzas_deg = collect(0.0:5.0:80.0)
vazs_deg = collect(0.0:5.0:355.0)
n_vza, n_vaz = length(vzas_deg), length(vazs_deg)

vza_all = repeat(vzas_deg; inner=n_vaz)
vaz_all = repeat(vazs_deg; outer=n_vza)

function no_absorption_phase_params()
    # Minimal, no-absorption scene. The pressure thickness is 0.01 hPa, so
    # Rayleigh scattering above the canopy is negligible but the usual RT
    # machinery still has a valid atmospheric layer.
    return parameters_from_dict(Dict{String,Any}(
        "radiative_transfer" => Dict{String,Any}(
            "spec_bands" => [string("[", join(ν_cm⁻¹, " "), "]")],
            "surface" => ["LambertianSurfaceScalar(0.0)"],
            "quadrature_type" => "RadauQuad()",
            "polarization_type" => "Stokes_I()",
            "max_m" => 8,
            "Δ_angle" => 2.0,
            "l_trunc" => 64,
            "depol" => -1.0,
            "float_type" => "Float64",
            "architecture" => "CPU()",
        ),
        "geometry" => Dict{String,Any}(
            "sza" => sza_deg,
            "vza" => vza_all,
            "vaz" => vaz_all,
            "obs_alt" => 1000.0,
        ),
        "atmospheric_profile" => Dict{String,Any}(
            "T" => [285.0],
            "p" => [1012.99, 1013.0],
            "profile_reduction" => -1,
        ),
    ))
end

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

opti = CanopyOptics.createLeafOpticalStruct(400.0:1.0:2500.0)
T_leaf_grid, R_leaf_grid = CanopyOptics.prospect(leaf, opti)
λ_grid = [Float64(v.val) for v in opti.λ]
grid_idx = [findmin(abs.(λ_grid .- λ))[2] for λ in λ_nm]
R_leaf = R_leaf_grid[grid_idx]
T_leaf = T_leaf_grid[grid_idx]

params = no_absorption_phase_params()
params.brdf[1] = CanopySurface(;
    soil = LambertianSurfaceScalar(0.0),
    LAI = 4.0,
    n_layers = 1,
    LAD = CanopyOptics.planophile_leaves2(Float64),
    leaf_reflectance = R_leaf,
    leaf_transmittance = T_leaf,
    leaf_optics_grid = λ_nm,
    grid_unit = :nm,
)

model = model_from_parameters(params)
R, _ = rt_run(model)
println("R size (dirs, stokes, λ): ", size(R))

R_red_grid = permutedims(reshape(R[:, 1, 1], n_vaz, n_vza))
R_nir_grid = permutedims(reshape(R[:, 1, 2], n_vaz, n_vza))
NDVI_grid = (R_nir_grid .- R_red_grid) ./ (R_nir_grid .+ R_red_grid)

nadir_ndvi = mean(NDVI_grid[1, :])
ΔNDVI_grid = NDVI_grid .- nadir_ndvi

println("685 nm reflectance extrema: ", extrema(R_red_grid))
println("800 nm reflectance extrema: ", extrema(R_nir_grid))
println("NDVI extrema: ", extrema(NDVI_grid))
println("NDVI at nadir: ", nadir_ndvi)
println("ΔNDVI from nadir extrema: ", extrema(ΔNDVI_grid))

θ_rad = deg2rad.(vcat(vazs_deg, 360.0))
NDVI_plot = hcat(NDVI_grid, NDVI_grid[:, 1])
ΔNDVI_plot = hcat(ΔNDVI_grid, ΔNDVI_grid[:, 1])
Δlim = maximum(abs, ΔNDVI_grid)

fig = Figure(size = (1180, 640))

ax_ndvi = PolarAxis(fig[1, 1];
    title = "NDVI phase space",
    theta_0 = π / 2,
    direction = -1,
    rticks = 0:20:80,
    rticklabelsize = 11,
    thetaticks = (deg2rad.(0:30:330), string.(0:30:330) .* "°"),
)
hm_ndvi = surface!(ax_ndvi, θ_rad, vzas_deg, NDVI_plot';
                   colormap = :viridis, shading = NoShading)
Colorbar(fig[1, 2], hm_ndvi; label = "NDVI")

ax_delta = PolarAxis(fig[1, 3];
    title = "NDVI - NDVI(nadir)",
    theta_0 = π / 2,
    direction = -1,
    rticks = 0:20:80,
    rticklabelsize = 11,
    thetaticks = (deg2rad.(0:30:330), string.(0:30:330) .* "°"),
)
hm_delta = surface!(ax_delta, θ_rad, vzas_deg, ΔNDVI_plot';
                    colormap = :balance, colorrange = (-Δlim, Δlim),
                    shading = NoShading)
Colorbar(fig[1, 4], hm_delta; label = "ΔNDVI")

# In vSmartMOM's relative-azimuth convention, the hotspot/backscatter side is
# at VAZ = 180 deg and VZA = SZA.
for ax in (ax_ndvi, ax_delta)
    scatter!(ax, [deg2rad(180.0)], [sza_deg];
             color = :red, markersize = 12, marker = :star5,
             strokecolor = :white, strokewidth = 1.0)
end

Label(fig[0, 1:4],
      "Canopy NDVI BRDF effect, soil albedo = 0.0 (SZA = 30 deg, LAI = 4, planophile LAD)";
      fontsize = 18)

out_png = joinpath(@__DIR__, "canopy_ndvi_phase_albedo0_685_800nm.png")
save(out_png, fig)
println("Saved NDVI phase-space plot -> ", out_png)
fig
