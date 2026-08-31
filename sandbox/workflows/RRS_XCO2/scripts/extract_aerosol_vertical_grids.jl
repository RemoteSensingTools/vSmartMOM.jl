#!/usr/bin/env julia

"""Export exact layer geometry and model-used aerosol AOD at 760 nm."""

include(joinpath(@__DIR__, "generate_truth_map.jl"))
using DelimitedFiles

const VERTICAL_OUT = joinpath(ROOT, "truth_map", "aerosol_vertical_grids.dat")

function model_for_layers(state, nlayer, output_ν, solve_ν)
    p = parameters_from_yaml(RAMAN_CONFIG)
    current = load_parameters()
    p.scattering_params = scattering_as_ft32(current.scattering_params)
    p.p, p.T, p.q = FT.(current.p), FT.(current.T), FT.(current.q)
    set_common!(p, state)
    p.profile_reduction_n = nlayer
    surface = transformed_surface(state.coeff[1], output_ν, solve_ν)
    select_band!(p, 1, solve_ν, surface)
    return model_from_parameters(p; external_solar=false)
end

function main_vertical()
    state = read_states()[9]
    output_ν = output_grids()[1]
    solve_ν, _ = o2_solve_grid(output_ν)
    open(VERTICAL_OUT, "w") do io
        println(io, "# Exact model geometry and untruncated aerosol AOD at the grid point nearest 760 nm")
        println(io, "# nlayer layer z_top_km z_bottom_km z_center_km dz_km sulfate_aod760 organic_aod760 stratospheric_aod760")
        for nlayer in (12, 16)
            model = model_for_layers(state, nlayer, output_ν, solve_ν)
            z_half = CoreRT.half_level_altitudes(model.profile)
            j760 = argmin(abs.(solve_ν .- FT(1e7 / 760)))
            for iz in eachindex(model.profile.Δz)
                aod = [Array(model.τ_aer[1])[ia, j760, iz] for ia in 1:3]
                @printf(io, "%d %d %.9f %.9f %.9f %.9f %.12e %.12e %.12e\n",
                        nlayer, iz, z_half[iz], z_half[iz+1],
                        (z_half[iz] + z_half[iz+1]) / 2,
                        model.profile.Δz[iz] / 1000, aod...)
            end
        end
    end
    println(VERTICAL_OUT)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main_vertical()
