# =============================================================================
# regression_newLUT.jl — Quick regression check of the rebuilt HITRAN LUTs
#
# Builds a vSmartMOM model with the 7-species new LUTs, extracts total column
# τ_abs + τ_rayl, prints them at diagnostic wavelengths alongside the matching
# values from MODTRAN_out/MODTRAN_H2O0.0500_AOT0.0001_GNDALT0.000_TSZ168.0.dat,
# and checks for the pathologies we saw in the old cubic LUT (negative τ,
# pointwise spikes inside the O2 A-band).
#
# Run: julia --project=test test/benchmarks/regression_newLUT.jl
# =============================================================================

using vSmartMOM
using vSmartMOM.CoreRT
using Printf
using Statistics

const YAML_PATH = joinpath(pkgdir(vSmartMOM), "test", "test_parameters",
                           "ParamsEMIT_MODTRANcomp_newLUT.yaml")
const MODTRAN_REF = expanduser("~/data/EMIT_MODTRANcomp/MODTRAN_out/MODTRAN_H2O0.0500_AOT0.0001_GNDALT0.000_TSZ168.0.dat")

function load_modtran_table()
    data = Dict("λ"=>Float64[], "td_dir"=>Float64[], "tu_dir"=>Float64[])
    open(MODTRAN_REF) do io
        readline(io); readline(io)
        for line in eachline(io)
            vals = parse.(Float64, split(line))
            push!(data["λ"], vals[1])
            push!(data["td_dir"], vals[4])
            push!(data["tu_dir"], vals[6])
        end
    end
    return data
end

"""Find nearest index to `val` in sorted vector `xs`."""
function nearest_idx(xs, val)
    _, i = findmin(abs.(xs .- val))
    return i
end

function diag()
    @info "Loading YAML + building model" YAML_PATH
    params = parameters_from_yaml(YAML_PATH)
    model  = model_from_parameters(params)

    # τ_abs[band] is (nSpec, nLayers), so sum over layers (dim 2).
    τ_abs_col  = vec(sum(model.τ_abs[1],  dims = 2))
    τ_rayl_col = vec(sum(model.τ_rayl[1], dims = 2))
    τ_tot_col  = τ_abs_col .+ τ_rayl_col

    ν_axis = collect(params.spec_bands[1])
    λ_nm   = 1e7 ./ ν_axis
    μ_s    = cos(deg2rad(Float64(params.sza)))

    MOD = load_modtran_table()

    # Sanity flags on the whole spectrum.
    n_neg_abs    = count(<(-1e-3), τ_abs_col)
    n_spike_abs  = count(>(10.0),  τ_abs_col)
    n_nan_abs    = count(isnan,    τ_abs_col)
    @info "Sanity of the new τ_abs column" nλ = length(τ_abs_col) n_neg = n_neg_abs n_spikes_gt10 = n_spike_abs n_nan = n_nan_abs

    # Point-check at diagnostic wavelengths.
    check_λ = [380.0, 400.0, 450.0, 500.0, 550.0, 630.0,
               760.0, 761.0, 762.5, 765.0,          # O2 A-band
               940.0, 1270.0, 1580.0,
               1930.0, 1940.0, 1950.0,              # H2O 1.94 μm band
               2000.0, 2100.0, 2160.0, 2200.0, 2300.0, 2500.0]

    @printf("\n  %-8s  %-11s  %-11s  %-11s  %-11s  %-11s\n",
            "λ_nm", "τ_abs", "τ_rayl", "τ_vert_tot", "T_dir_us", "T_dir_MODTRAN")
    println("-"^80)
    for λt in check_λ
        i  = nearest_idx(λ_nm, λt)
        iM = nearest_idx(MOD["λ"], λt)
        Tus  = exp(-τ_tot_col[i] / μ_s)
        Tmod = MOD["td_dir"][iM]
        @printf("  %-8.1f  %-11.4e  %-11.4e  %-11.4e  %-11.6f  %-11.6f\n",
                λt, τ_abs_col[i], τ_rayl_col[i], τ_tot_col[i], Tus, Tmod)
    end

    println("\n--- O2 A-band fine scan (759.5–762 nm) ---")
    abands = findall(j -> 759.5 ≤ λ_nm[j] ≤ 762.0, eachindex(λ_nm))
    τslice = τ_tot_col[abands]
    @printf("  n pts in window: %d   min τ: %.3e   max τ: %.3e   median τ: %.3e\n",
            length(τslice), minimum(τslice), maximum(τslice),
            sort(τslice)[(length(τslice)+1)÷2])

    println("\n--- Summary ---")
    println("  any τ_abs < -1e-3?  → ", n_neg_abs, " points")
    println("  any τ_abs > 10?     → ", n_spike_abs, " points")
    println("  any τ_abs = NaN?    → ", n_nan_abs, " points")
    return (; τ_abs_col, τ_rayl_col, ν_axis, λ_nm, MOD)
end

if abspath(PROGRAM_FILE) == @__FILE__
    diag()
end
