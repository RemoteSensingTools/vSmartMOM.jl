# =================================================================
# Fast canopy BRDF and NDVI slice at two wavelengths: red (685 nm) and NIR
# (800 nm)
#
# This is a no-gas, no-LUT sandbox. It constructs the RT parameters in Julia,
# runs 685 and 800 nm together as one two-point spectral band, and uses a
# vanishingly thin atmosphere so the plotted signal is essentially the
# canopy+soil BRDF.
#
# Run with:
#   julia --project=docs sandbox/canopy_red_nir_brdf_685_800nm.jl
# =================================================================

using vSmartMOM
using vSmartMOM.CoreRT
using CanopyOptics
using CairoMakie

const λ_nm = [685.0, 800.0]
const ν_cm⁻¹ = 1e7 ./ λ_nm
const soil_cases = [
    ("soil albedo = 0.0", 0.0),
    ("soil albedo = 1.0", 1.0),
]
const sza_deg = 30.0

# Principal-plane geometry:
#   negative VZA = VAZ 180 deg, backscatter / hotspot side
#   positive VZA = VAZ   0 deg, forward-scatter side
signed_vza = collect(-80.0:5.0:80.0)
vza = abs.(signed_vza)
vaz = ifelse.(signed_vza .< 0, 180.0, 0.0)

function no_absorption_params()
    # Minimal, no-absorption scene. The pressure thickness is 0.01 hPa, so
    # Rayleigh scattering above the canopy is negligible but the usual RT
    # machinery still has a valid atmospheric layer.
    return parameters_from_dict(Dict{String,Any}(
        "radiative_transfer" => Dict{String,Any}(
            "spec_bands" => [string("[", join(ν_cm⁻¹, " "), "]")],
            "surface" => ["LambertianSurfaceScalar(0.0)"],
            "quadrature_type" => "RadauQuad()",
            "polarization_type" => "Stokes_I()",
            "max_m" => 6,
            "Δ_angle" => 2.0,
            "l_trunc" => 32,
            "depol" => -1.0,
            "float_type" => "Float64",
            "architecture" => "CPU()",
        ),
        "geometry" => Dict{String,Any}(
            "sza" => sza_deg,
            "vza" => vza,
            "vaz" => vaz,
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

# PROSPECT is evaluated on a 1 nm grid once, then the two requested leaf R/T
# points are handed to the spectral canopy surface.
opti = CanopyOptics.createLeafOpticalStruct(400.0:1.0:2500.0)
T_leaf_grid, R_leaf_grid = CanopyOptics.prospect(leaf, opti)
λ_grid = [Float64(v.val) for v in opti.λ]
grid_idx = [findmin(abs.(λ_grid .- λ))[2] for λ in λ_nm]
R_leaf = R_leaf_grid[grid_idx]
T_leaf = T_leaf_grid[grid_idx]

R_principal = zeros(Float64, length(signed_vza), length(λ_nm), length(soil_cases))
for (isoil, (soil_label, soil_input)) in pairs(soil_cases)
    params = no_absorption_params()
    params.brdf[1] = CanopySurface(;
        soil = LambertianSurfaceScalar(soil_input),
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
    R_principal[:, :, isoil] .= R[:, 1, :]

    println(soil_label, " (input ", soil_input, "), R size: ", size(R))
    for (iλ, λ) in pairs(λ_nm)
        println("  λ = ", λ, " nm, leaf R/T = ", (R_leaf[iλ], T_leaf[iλ]),
                ", principal-plane R(I) min/max: ",
                extrema(R_principal[:, iλ, isoil]))
        for x in (-sza_deg, 0.0, sza_deg)
            idx = findfirst(==(x), signed_vza)
            idx === nothing || println("    R(I), signed VZA = ", x, " deg: ",
                                       R_principal[idx, iλ, isoil])
        end
    end
end

NDVI_principal = (R_principal[:, 2, :] .- R_principal[:, 1, :]) ./
                 (R_principal[:, 2, :] .+ R_principal[:, 1, :])
nadir_idx = findfirst(==(0.0), signed_vza)
NDVI_nadir = NDVI_principal[nadir_idx, :]
NDVI_delta_from_nadir = NDVI_principal .- reshape(NDVI_nadir, 1, :)

for (isoil, (soil_label, _)) in pairs(soil_cases)
    println(soil_label, ", NDVI min/max: ", extrema(NDVI_principal[:, isoil]))
    println("  NDVI at nadir: ", NDVI_principal[nadir_idx, isoil])
    for x in (-sza_deg, sza_deg)
        idx = findfirst(==(x), signed_vza)
        idx === nothing || println("  NDVI, signed VZA = ", x, " deg: ",
                                   NDVI_principal[idx, isoil],
                                   ", ΔNDVI from nadir: ",
                                   NDVI_delta_from_nadir[idx, isoil])
    end
end

fig = Figure(size = (1080, 520))
colors = [:firebrick3, :darkgreen]
labels = ["685 nm (red)", "800 nm (NIR)"]

for (isoil, (soil_label, _)) in pairs(soil_cases)
    ax = Axis(fig[1, isoil];
        xlabel = "Signed VZA in solar principal plane (deg)",
        ylabel = isoil == 1 ? "TOA reflectance (Stokes I)" : "",
        title = soil_label,
    )

    for iλ in eachindex(λ_nm)
        lines!(ax, signed_vza, R_principal[:, iλ, isoil];
               color = colors[iλ], linewidth = 2.5, label = labels[iλ])
        scatter!(ax, signed_vza, R_principal[:, iλ, isoil];
                 color = colors[iλ], markersize = 6)
    end

    vlines!(ax, [-sza_deg]; color = :gray35, linestyle = :dash, linewidth = 1.5,
            label = "Backscatter / hotspot")
    vlines!(ax, [0.0]; color = :gray55, linestyle = :dot, linewidth = 1.3,
            label = "Nadir")

    text!(ax, -78, minimum(R_principal[:, :, isoil]);
          text = "negative: VAZ = 180 deg\npositive: VAZ = 0 deg",
          align = (:left, :bottom), fontsize = 12, color = :gray30)

    axislegend(ax; position = :lt, framevisible = false)
    xlims!(ax, -82, 82)
end

Label(fig[0, 1:2],
      "Planophile-LAD one-layer canopy BRDF slice (SZA = 30 deg, LAI = 4)";
      fontsize = 18)

out_png = joinpath(@__DIR__, "canopy_red_nir_brdf_685_800nm.png")
save(out_png, fig)
println("Saved BRDF slice plot -> ", out_png)

fig_ndvi = Figure(size = (1080, 520))
soil_colors = [:teal, :purple4]

ax_ndvi = Axis(fig_ndvi[1, 1];
    xlabel = "Signed VZA in solar principal plane (deg)",
    ylabel = "NDVI",
    title = "NDVI from 800 and 685 nm reflectance",
)
ax_delta = Axis(fig_ndvi[1, 2];
    xlabel = "Signed VZA in solar principal plane (deg)",
    ylabel = "NDVI - NDVI(nadir)",
    title = "Angular NDVI BRDF effect",
)

for (isoil, (soil_label, _)) in pairs(soil_cases)
    lines!(ax_ndvi, signed_vza, NDVI_principal[:, isoil];
           color = soil_colors[isoil], linewidth = 2.5, label = soil_label)
    scatter!(ax_ndvi, signed_vza, NDVI_principal[:, isoil];
             color = soil_colors[isoil], markersize = 6)

    lines!(ax_delta, signed_vza, NDVI_delta_from_nadir[:, isoil];
           color = soil_colors[isoil], linewidth = 2.5, label = soil_label)
    scatter!(ax_delta, signed_vza, NDVI_delta_from_nadir[:, isoil];
             color = soil_colors[isoil], markersize = 6)
end

for ax in (ax_ndvi, ax_delta)
    hlines!(ax, [0.0]; color = :gray70, linewidth = 1.0)
    vlines!(ax, [-sza_deg]; color = :gray35, linestyle = :dash, linewidth = 1.5,
            label = "Backscatter / hotspot")
    vlines!(ax, [0.0]; color = :gray55, linestyle = :dot, linewidth = 1.3,
            label = "Nadir")
    xlims!(ax, -82, 82)
end

axislegend(ax_ndvi; position = :rb, framevisible = false)
axislegend(ax_delta; position = :lb, framevisible = false)

Label(fig_ndvi[0, 1:2],
      "Canopy NDVI angular dependence (SZA = 30 deg, LAI = 4, planophile LAD)";
      fontsize = 18)

out_ndvi_png = joinpath(@__DIR__, "canopy_ndvi_brdf_effect_685_800nm.png")
save(out_ndvi_png, fig_ndvi)
println("Saved NDVI BRDF-effect plot -> ", out_ndvi_png)
fig_ndvi
