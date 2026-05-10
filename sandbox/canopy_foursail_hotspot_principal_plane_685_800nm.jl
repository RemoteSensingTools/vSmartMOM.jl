# =================================================================
# Principal-plane 4SAIL hotspot diagnostic against coupled vSmartMOM
#
# This is a no-gas red/NIR canopy comparison:
#   - wavelengths: 685 and 800 nm
#   - SZA = 30 deg
#   - LAI = 4
#   - one-layer planophile canopy
#   - black Lambertian soil
#   - 4SAIL hotspot h = 0.05
#
# vSmartMOM currently has no canopy-hotspot source-term hook here, so its curve
# is the coupled multiple-scattering canopy BRF without hotspot. 4SAIL is shown
# both with and without hotspot, plus the embedded direct/direct term `rsost`.
#
# Principal-plane convention:
#   signed VZA < 0 => VAZ = 180 deg, backscatter / hotspot side
#   signed VZA > 0 => VAZ =   0 deg, forward-scatter side
#
# Run with:
#   julia --project=docs sandbox/canopy_foursail_hotspot_principal_plane_685_800nm.jl
# =================================================================

using vSmartMOM
using vSmartMOM.CoreRT
using CanopyOptics
using CairoMakie
using Statistics

const lambda_nm = [685.0, 800.0]
const nu_cm_inv = 1e7 ./ lambda_nm
const sza_deg = 30.0
const μ0 = cosd(sza_deg)
const lai = 4.0
const soil_albedo = 0.0
const hotspot_h = 0.05

signed_vza = collect(-80.0:2.0:80.0)
vza = abs.(signed_vza)
vaz = ifelse.(signed_vza .< 0, 180.0, 0.0)
raa_4sail = abs.(vaz .- 180.0)

function no_absorption_params()
    return parameters_from_dict(Dict{String,Any}(
        "radiative_transfer" => Dict{String,Any}(
            "spec_bands" => [string("[", join(nu_cm_inv, " "), "]")],
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

opti = CanopyOptics.createLeafOpticalStruct(400.0:1.0:2500.0)
T_leaf_grid, R_leaf_grid = CanopyOptics.prospect(leaf, opti)
lambda_grid = [Float64(v.val) for v in opti.λ]
grid_idx = [findmin(abs.(lambda_grid .- lambda))[2] for lambda in lambda_nm]
R_leaf = R_leaf_grid[grid_idx]
T_leaf = T_leaf_grid[grid_idx]

lad = CanopyOptics.planophile_leaves2(Float64)

params = no_absorption_params()
params.brdf[1] = CanopySurface(;
    soil = LambertianSurfaceScalar(soil_albedo),
    LAI = lai,
    n_layers = 1,
    LAD = lad,
    leaf_reflectance = R_leaf,
    leaf_transmittance = T_leaf,
    leaf_optics_grid = lambda_nm,
    grid_unit = :nm,
)

model = model_from_parameters(params)
R_vmom, _ = rt_run(model)
vmom_brf = (π / μ0) .* R_vmom[:, 1, :]

geoms = CanopyOptics.FourSAILGeometrySet(lad;
    sza_deg = sza_deg,
    vza_deg = vza,
    raa_deg = raa_4sail,
    quadrature = CanopyOptics.CanopyQuadrature(n_leaf = 64, n_azimuth = 16),
)

sail_nohot = CanopyOptics.foursail(R_leaf, T_leaf, soil_albedo, geoms, lai;
                                   hotspot = 0.0, mode = :full)
sail_hot = CanopyOptics.foursail(R_leaf, T_leaf, soil_albedo, geoms, lai;
                                 hotspot = hotspot_h, mode = :full)

println("vSmartMOM raw SFI converted to BRF with π / μ0 = ", π / μ0)
println("4SAIL hotspot parameter h = ", hotspot_h)

hotspot_idx = findfirst(==(-sza_deg), signed_vza)
nadir_idx = findfirst(==(0.0), signed_vza)
forward_idx = findfirst(==(sza_deg), signed_vza)

for (ilambda, lambda) in pairs(lambda_nm)
    println("lambda = ", lambda, " nm")
    println("  leaf R/T: ", (R_leaf[ilambda], T_leaf[ilambda]))
    println("  vSmartMOM BRF extrema: ", extrema(vmom_brf[:, ilambda]))
    println("  4SAIL no-hotspot full extrema: ", extrema(sail_nohot.rsot[ilambda, :]))
    println("  4SAIL hotspot full extrema: ", extrema(sail_hot.rsot[ilambda, :]))
    println("  4SAIL hotspot increment extrema: ",
            extrema(sail_hot.rsot[ilambda, :] .- sail_nohot.rsot[ilambda, :]))
    for (label, idx) in (("hotspot", hotspot_idx),
                         ("nadir", nadir_idx),
                         ("forward", forward_idx))
        idx === nothing && continue
        println("  ", label, " signed VZA = ", signed_vza[idx], " deg:")
        println("    vSmartMOM = ", vmom_brf[idx, ilambda])
        println("    4SAIL no-hotspot full = ", sail_nohot.rsot[ilambda, idx])
        println("    4SAIL hotspot full = ", sail_hot.rsot[ilambda, idx])
        println("    4SAIL hotspot direct rsost = ", sail_hot.rsost[ilambda, idx])
    end
end

fig = Figure(size = (1250, 760))
colors = (
    vmom = :black,
    sail_nohot = :dodgerblue3,
    sail_hot = :seagreen4,
    sail_direct = :darkorange3,
)

for (ilambda, lambda) in pairs(lambda_nm)
    ax = Axis(fig[ilambda, 1];
        xlabel = ilambda == length(lambda_nm) ? "Signed VZA in solar principal plane (deg)" : "",
        ylabel = "BRF",
        title = "$(Int(round(lambda))) nm principal-plane BRF",
    )
    lines!(ax, signed_vza, vmom_brf[:, ilambda];
           color = colors.vmom, linewidth = 2.8, label = "vSmartMOM coupled")
    lines!(ax, signed_vza, sail_nohot.rsot[ilambda, :];
           color = colors.sail_nohot, linewidth = 2.2, linestyle = :dash,
           label = "4SAIL full, h = 0")
    lines!(ax, signed_vza, sail_hot.rsot[ilambda, :];
           color = colors.sail_hot, linewidth = 2.5,
           label = "4SAIL full, h = $(hotspot_h)")
    lines!(ax, signed_vza, sail_hot.rsost[ilambda, :];
           color = colors.sail_direct, linewidth = 2.0, linestyle = :dot,
           label = "4SAIL direct rsost, h = $(hotspot_h)")

    axd = Axis(fig[ilambda, 2];
        xlabel = ilambda == length(lambda_nm) ? "Signed VZA in solar principal plane (deg)" : "",
        ylabel = "BRF difference",
        title = "$(Int(round(lambda))) nm 4SAIL - vSmartMOM",
    )
    lines!(axd, signed_vza, sail_nohot.rsot[ilambda, :] .- vmom_brf[:, ilambda];
           color = colors.sail_nohot, linewidth = 2.2, linestyle = :dash,
           label = "full, h = 0")
    lines!(axd, signed_vza, sail_hot.rsot[ilambda, :] .- vmom_brf[:, ilambda];
           color = colors.sail_hot, linewidth = 2.5,
           label = "full, h = $(hotspot_h)")
    lines!(axd, signed_vza, sail_hot.rsot[ilambda, :] .- sail_nohot.rsot[ilambda, :];
           color = :purple4, linewidth = 2.0, linestyle = :dot,
           label = "hotspot increment")
    hlines!(axd, [0.0]; color = :gray70, linewidth = 1.0)

    for a in (ax, axd)
        vlines!(a, [-sza_deg]; color = :gray30, linestyle = :dashdot,
                linewidth = 1.4, label = "backscatter/hotspot")
        vlines!(a, [0.0]; color = :gray55, linestyle = :dot,
                linewidth = 1.2, label = "nadir")
        xlims!(a, -82, 82)
    end

    axislegend(ax; position = :lt, framevisible = false)
    axislegend(axd; position = :lb, framevisible = false)
end

Label(fig[0, 1:2],
      "4SAIL hotspot principal-plane comparison against coupled vSmartMOM (soil albedo = 0, SZA = 30 deg, LAI = 4)";
      fontsize = 18)

out_png = joinpath(@__DIR__, "canopy_foursail_hotspot_principal_plane_685_800nm.png")
save(out_png, fig)
println("Saved hotspot principal-plane plot -> ", out_png)
fig
