# =================================================================
# Compare 4SAIL single-scattering BRDF against vSmartMOM canopy-only RT
#
# This uses a canopy-only red/NIR setup:
#   - wavelengths: 685 and 800 nm
#   - SZA = 30 deg
#   - LAI = 4
#   - one-layer planophile canopy
#   - black Lambertian soil
#   - atmospheric Rayleigh, aerosol, and gas absorption optical depths are
#     explicitly set to zero after model construction
#   - Δ_angle = 0.0 so none of the usual δ-truncation/exclusion-angle machinery
#     is requested; canopy Fourier moments are carried through m = 0:11
#
# vSmartMOM uses VAZ = 180 deg for the backscatter/hotspot side. The 4SAIL
# relative azimuth convention has that same geometry at RAA = 0 deg, so this
# script maps RAA_4SAIL = abs(VAZ_vSmartMOM - 180 deg).
#
# rt_run returns radiance normalized by unit solar irradiance. For a BRF/HDRF
# comparison against 4SAIL rsot/rsost, convert vSmartMOM by π / μ0. The 4SAIL
# direct single-scattering panel uses the embedded `rsost` component from the
# same full solve (`rsot = rsost + rsodt`), not a separate reduced-mode run.
#
# Run with:
#   julia --project=docs sandbox/canopy_foursail_ss_vs_vsmartmom_brdf_685_800nm.jl
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

vzas_deg = collect(0.0:5.0:80.0)
vazs_deg = collect(0.0:5.0:355.0)
n_vza, n_vaz = length(vzas_deg), length(vazs_deg)

vza_all = repeat(vzas_deg; inner=n_vaz)
vaz_all = repeat(vazs_deg; outer=n_vza)
raa_4sail = abs.(vaz_all .- 180.0)

function canopy_only_phase_params()
    return parameters_from_dict(Dict{String,Any}(
        "radiative_transfer" => Dict{String,Any}(
            "spec_bands" => [string("[", join(nu_cm_inv, " "), "]")],
            "surface" => ["LambertianSurfaceScalar(0.0)"],
            "quadrature_type" => "RadauQuad()",
            "polarization_type" => "Stokes_I()",
            "max_m" => 12,
            "Δ_angle" => 0.0,
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
            # One valid layer is still needed by the generic RT driver. Its
            # optical depths are zeroed below, so it is a transparent identity
            # layer rather than a Rayleigh/absorbing atmosphere.
            "T" => [285.0],
            "p" => [1012.99, 1013.0],
            "profile_reduction" => -1,
        ),
    ))
end

function zero_atmosphere!(model)
    for τ in model.τ_rayl
        fill!(τ, 0)
    end
    for τ in model.τ_abs
        fill!(τ, 0)
    end
    for τ in model.τ_aer
        fill!(τ, 0)
    end
    return model
end

max_or_zero(arrays) = maximum([isempty(A) ? 0.0 : maximum(A) for A in arrays])

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

params = canopy_only_phase_params()
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
zero_atmosphere!(model)
println("Atmosphere zeroed: max τ_rayl = ", max_or_zero(model.τ_rayl),
        ", max τ_abs = ", max_or_zero(model.τ_abs),
        ", max τ_aer = ", max_or_zero(model.τ_aer))
R_vmom, _ = rt_run(model)
println("vSmartMOM R size (dirs, stokes, wavelength): ", size(R_vmom))
println("Converting raw vSmartMOM SFI output to BRF/HDRF with π / μ0 = ",
        π / μ0)

vmom_grid = Array{Float64}(undef, length(lambda_nm), n_vza, n_vaz)
for ilambda in eachindex(lambda_nm)
    vmom_grid[ilambda, :, :] .= (π / μ0) .*
                                 permutedims(reshape(R_vmom[:, 1, ilambda],
                                                     n_vaz, n_vza))
end

geoms = CanopyOptics.FourSAILGeometrySet(lad;
    sza_deg = sza_deg,
    vza_deg = vza_all,
    raa_deg = raa_4sail,
    quadrature = CanopyOptics.CanopyQuadrature(n_leaf = 64, n_azimuth = 16),
)

R_4sail = CanopyOptics.foursail(R_leaf, T_leaf, soil_albedo, geoms, lai;
                                hotspot = 0.0, mode = :full)

sail_full_grid = Array{Float64}(undef, length(lambda_nm), n_vza, n_vaz)
sail_ss_grid = Array{Float64}(undef, length(lambda_nm), n_vza, n_vaz)
for ilambda in eachindex(lambda_nm)
    sail_full_grid[ilambda, :, :] .= permutedims(reshape(R_4sail.rsot[ilambda, :],
                                                         n_vaz, n_vza))
    sail_ss_grid[ilambda, :, :] .= permutedims(reshape(R_4sail.rsost[ilambda, :],
                                                       n_vaz, n_vza))
end

diff_full_grid = sail_full_grid .- vmom_grid
diff_ss_grid = sail_ss_grid .- vmom_grid

for (ilambda, lambda) in pairs(lambda_nm)
    diff_full = diff_full_grid[ilambda, :, :]
    diff_ss = diff_ss_grid[ilambda, :, :]
    println("lambda = ", lambda, " nm")
    println("  leaf R/T: ", (R_leaf[ilambda], T_leaf[ilambda]))
    println("  vSmartMOM canopy-only BRF extrema: ", extrema(vmom_grid[ilambda, :, :]))
    println("  4SAIL full extrema: ", extrema(sail_full_grid[ilambda, :, :]))
    println("  4SAIL SS extrema: ", extrema(sail_ss_grid[ilambda, :, :]))
    println("  4SAIL full - vSmartMOM extrema: ", extrema(diff_full))
    println("  4SAIL SS - vSmartMOM extrema: ", extrema(diff_ss))
    println("  4SAIL full mean absolute difference: ", mean(abs.(diff_full)))
    println("  4SAIL SS mean absolute difference: ", mean(abs.(diff_ss)))
end

theta_rad = deg2rad.(vcat(vazs_deg, 360.0))

fig = Figure(size = (2040, 920))

for (ilambda, lambda) in pairs(lambda_nm)
    vmom_plot = hcat(vmom_grid[ilambda, :, :], vmom_grid[ilambda, :, 1])
    sail_full_plot = hcat(sail_full_grid[ilambda, :, :], sail_full_grid[ilambda, :, 1])
    sail_ss_plot = hcat(sail_ss_grid[ilambda, :, :], sail_ss_grid[ilambda, :, 1])
    diff_ss_plot = hcat(diff_ss_grid[ilambda, :, :], diff_ss_grid[ilambda, :, 1])

    refl_range = extrema(vcat(vec(vmom_grid[ilambda, :, :]),
                              vec(sail_full_grid[ilambda, :, :]),
                              vec(sail_ss_grid[ilambda, :, :])))
    diff_lim = maximum(abs, diff_ss_grid[ilambda, :, :])

    ax_vmom = PolarAxis(fig[ilambda, 1];
        title = "$(Int(round(lambda))) nm vSmartMOM canopy-only BRF",
        theta_0 = pi / 2,
        direction = -1,
        rticks = 0:20:80,
        rticklabelsize = 10,
        thetaticks = (deg2rad.(0:45:315), string.(0:45:315) .* " deg"),
    )
    hm_vmom = surface!(ax_vmom, theta_rad, vzas_deg, vmom_plot';
                       colormap = :viridis, colorrange = refl_range,
                       shading = NoShading)

    ax_sail = PolarAxis(fig[ilambda, 3];
        title = "$(Int(round(lambda))) nm 4SAIL full",
        theta_0 = pi / 2,
        direction = -1,
        rticks = 0:20:80,
        rticklabelsize = 10,
        thetaticks = (deg2rad.(0:45:315), string.(0:45:315) .* " deg"),
    )
    hm_sail = surface!(ax_sail, theta_rad, vzas_deg, sail_full_plot';
                       colormap = :viridis, colorrange = refl_range,
                       shading = NoShading)

    ax_ss = PolarAxis(fig[ilambda, 5];
        title = "$(Int(round(lambda))) nm 4SAIL direct SS",
        theta_0 = pi / 2,
        direction = -1,
        rticks = 0:20:80,
        rticklabelsize = 10,
        thetaticks = (deg2rad.(0:45:315), string.(0:45:315) .* " deg"),
    )
    hm_ss = surface!(ax_ss, theta_rad, vzas_deg, sail_ss_plot';
                     colormap = :viridis, colorrange = refl_range,
                     shading = NoShading)

    ax_diff = PolarAxis(fig[ilambda, 7];
        title = "$(Int(round(lambda))) nm 4SAIL SS - vSmartMOM",
        theta_0 = pi / 2,
        direction = -1,
        rticks = 0:20:80,
        rticklabelsize = 10,
        thetaticks = (deg2rad.(0:45:315), string.(0:45:315) .* " deg"),
    )
    hm_diff = surface!(ax_diff, theta_rad, vzas_deg, diff_ss_plot';
                       colormap = :balance, colorrange = (-diff_lim, diff_lim),
                       shading = NoShading)

    for ax in (ax_vmom, ax_sail, ax_ss, ax_diff)
        scatter!(ax, [deg2rad(180.0)], [sza_deg];
                 color = :red, markersize = 11, marker = :star5,
                 strokecolor = :white, strokewidth = 1.0)
    end

    Colorbar(fig[ilambda, 2], hm_vmom; label = "BRF")
    Colorbar(fig[ilambda, 4], hm_sail; label = "BRF")
    Colorbar(fig[ilambda, 6], hm_ss; label = "BRF")
    Colorbar(fig[ilambda, 8], hm_diff; label = "Difference")
end

Label(fig[0, 1:8],
      "4SAIL full and direct single-scattering BRF vs vSmartMOM canopy-only BRF (soil albedo = 0, SZA = 30 deg, LAI = 4)";
      fontsize = 18)

out_png = joinpath(@__DIR__, "canopy_foursail_ss_vs_vsmartmom_vacuum_brdf_685_800nm.png")
save(out_png, fig)
println("Saved comparison plot -> ", out_png)
fig
