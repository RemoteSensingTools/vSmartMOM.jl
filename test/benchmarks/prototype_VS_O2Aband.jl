# =============================================================================
# prototype_VS_O2Aband.jl — Vibrational Scattering (VS) smoke over the O2
# A-band (Phase 6 port from sanghavi).
#
# The sanghavi script was a 442-line notebook that (a) ran the VS forward
# model over a grid of source-wavenumber sections (1e7/725→1e7/625 cm⁻¹)
# writing each section's output to `/home/sanghavi/julia_output/`, and
# (b) post-processed the grid into convolved IQ plots. This port drops
# the grid-generation scaffold (trivial to add once a target output dir
# is chosen) and provides the core single-section VS + RRS + noRS
# compute path as a reproducible driver.
#
# Usage:
#   # default: one narrow ν₀ section (faster smoke), section #1 of 500
#   julia --project=test test/benchmarks/prototype_VS_O2Aband.jl
#
#   # a specific section (1..500)
#   VS_SECTION=5 julia --project=test test/benchmarks/prototype_VS_O2Aband.jl
#
#   # full span (expensive — runs over all 500 sections)
#   VS_SECTION=all julia --project=test test/benchmarks/prototype_VS_O2Aband.jl
#
#   # custom output dir (default: test/benchmarks/VS_O2Aband_output/)
#   VS_OUT_DIR=/tmp/vs_out julia --project=test test/benchmarks/prototype_VS_O2Aband.jl
# =============================================================================

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using vSmartMOM.SolarModel
using Statistics
using DelimitedFiles

const FT = Float64

const VS_YAML_VS  = get(ENV, "VS_YAML_VS",
                        joinpath(pkgdir(vSmartMOM), "test", "test_parameters",
                                 "O2ABandParamsVS.yaml"))
const VS_YAML_RRS = get(ENV, "VS_YAML_RRS",
                        joinpath(pkgdir(vSmartMOM), "test", "test_parameters",
                                 "O2ABandParamsRRS.yaml"))
const VS_OUT_DIR  = get(ENV, "VS_OUT_DIR",
                        joinpath(pkgdir(vSmartMOM), "test", "benchmarks",
                                 "VS_O2Aband_output"))

const VS_NSECT_TOTAL = 500
const VS_SECTION_SPEC = get(ENV, "VS_SECTION", "1")

# Optional solar-transmission table. Users who have it should set
# VS_SOLAR_FILE; otherwise we use a flat solar spectrum (transmission = 1).
const VS_SOLAR_FILE = get(ENV, "VS_SOLAR_FILE", "")

mkpath(VS_OUT_DIR)

"""
    run_vs_section(nsect; ν₀_full, νₐ, yaml_vs, yaml_rrs, out_dir)

Compute VS + RRS + noRS over a single source-wavenumber section `nsect`
of 1..VS_NSECT_TOTAL partitions of the full `ν₀_full` range. Writes two
`.dat` files to `out_dir`:
  - `out_src_O2Aband_625nm_725nm_{nsect}.dat` — source spectrum R₀, R₀_Q
  - `out_VS_625nm_725nm_{nsect}.dat`          — target O2-A band VS Rₐ, ieRₐ

Returns a NamedTuple with the computed arrays for the caller.
"""
function run_vs_section(nsect::Integer;
                        ν₀_full = (1e7/725.):0.1:(1e7/625.),
                        νₐ      = (1e7/780.):0.1:(1e7/750.),
                        yaml_vs::AbstractString = VS_YAML_VS,
                        yaml_rrs::AbstractString = VS_YAML_RRS,
                        out_dir::AbstractString = VS_OUT_DIR,
                        solar_file::AbstractString = VS_SOLAR_FILE,
                        N_sect::Integer = VS_NSECT_TOTAL)
    # Slice ν₀_full into the nsect-th section
    l_orig = length(ν₀_full)
    Δl     = Int(ceil(l_orig / N_sect))
    i_lo   = (nsect - 1) * Δl + 1
    i_hi   = min(nsect * Δl, l_orig)
    i_lo > l_orig && error("nsect $(nsect) out of range for $(N_sect) sections")
    ν₀ = ν₀_full[i_lo:i_hi]

    T_sun = 5777.0
    P₀    = planck_spectrum_wn(T_sun, collect(ν₀))
    Pₐ    = planck_spectrum_wn(T_sun, collect(νₐ))

    # Optional solar transmission. Fall back to unity when file absent.
    Tsolar_interp = if isfile(solar_file)
        t = solar_transmission_from_file(solar_file)
        using_interp(t[4:end, 1], t[4:end, 2])
    else
        ν -> 1.0
    end

    R₀    = zeros(FT, length(ν₀))
    R₀_Q  = zeros(FT, length(ν₀))
    Rₐ    = zeros(FT, length(νₐ))
    Rₐ_Q  = zeros(FT, length(νₐ))
    ieRₐ  = zeros(FT, length(νₐ))
    ieRₐ_Q = zeros(FT, length(νₐ))

    parameters_vs = parameters_from_yaml(yaml_vs)

    for (i, ν0_i) in enumerate(ν₀)
        λ₀ = 1e7 / ν0_i
        n2, o2 = InelasticScattering.getRamanAtmoConstants(1.7 / λ₀, FT(300.0))
        RS_type = InelasticScattering.VS_0to1_plus{FT}(n2 = n2, o2 = o2)

        model = model_from_parameters(RS_type, λ₀, parameters_vs)

        nSpec_total = sum(length(RS_type.grid_in[j])
                           for j in 1:length(RS_type.iBand))
        nPol = parameters_vs.polarization_type.n
        RS_type.F₀ = zeros(FT, nPol, nSpec_total)

        t_offset = 0
        corridx  = Int[]
        bandidx  = Int[]
        for iB in 1:length(RS_type.iBand)
            νiB = collect(RS_type.grid_in[iB])
            for iii in 1:length(νiB)
                ν_here = νiB[iii]
                soltmp = if iB == 1
                    Tsolar_interp(ν_here) * P₀[i]
                else
                    j = findmin(abs.(νₐ .- ν_here))[2]
                    if abs(νₐ[j] - ν_here) < 1
                        push!(corridx, j); push!(bandidx, iii + t_offset)
                        Tsolar_interp(ν_here) * Pₐ[j]
                    else
                        0.0
                    end
                end
                RS_type.F₀[1, iii + t_offset] = soltmp
            end
            t_offset += length(νiB)
        end

        R_rt, _, ieR_rt, _ = CoreRT.rt_run_test(RS_type, model, RS_type.iBand)

        R₀[i]   = R_rt[1, 1, 1]
        R₀_Q[i] = R_rt[1, 2, 1]
        for k in eachindex(corridx)
            Rₐ[corridx[k]]    = R_rt[1, 1, bandidx[k]]
            Rₐ_Q[corridx[k]]  = R_rt[1, 2, bandidx[k]]
            ieRₐ[corridx[k]]  += ieR_rt[1, 1, bandidx[k]]
            ieRₐ_Q[corridx[k]] += ieR_rt[1, 2, bandidx[k]]
        end
        @info "VS section" nsect=nsect i=i of=length(ν₀) ν0=ν0_i
    end

    mkpath(out_dir)
    writedlm(joinpath(out_dir, "out_src_O2Aband_625nm_725nm_$(nsect).dat"),
             [ν₀ R₀ R₀_Q])
    writedlm(joinpath(out_dir, "out_VS_625nm_725nm_$(nsect).dat"),
             [νₐ Rₐ Rₐ_Q ieRₐ ieRₐ_Q])

    return (; ν₀, R₀, R₀_Q, νₐ, Rₐ, Rₐ_Q, ieRₐ, ieRₐ_Q)
end

# Use a local helper for optional solar-transmission interpolation; avoids
# `using Interpolations` hitting users who don't have the solar data.
function using_interp(xs, ys)
    itp = sort(collect(zip(xs, ys)))
    sx = getindex.(itp, 1); sy = getindex.(itp, 2)
    function (ν)
        ν <= sx[1]   && return sy[1]
        ν >= sx[end] && return sy[end]
        j = searchsortedfirst(sx, ν)
        t = (ν - sx[j - 1]) / (sx[j] - sx[j - 1])
        sy[j - 1] + t * (sy[j] - sy[j - 1])
    end
end

"""
    run_rrs_o2aband(; yaml_rrs, iBand=1, out_dir=VS_OUT_DIR)

Compute RRS (rotational Raman) and noRS forward outputs over the O₂ A-band
and write `.dat` files consumable by the VS post-processing. Returns the
R/ieR/RnoRS arrays.
"""
function run_rrs_o2aband(; yaml_rrs::AbstractString = VS_YAML_RRS,
                           iBand::Integer = 1,
                           out_dir::AbstractString = VS_OUT_DIR,
                           solar_file::AbstractString = VS_SOLAR_FILE)
    parameters = parameters_from_yaml(yaml_rrs)
    model = model_from_parameters(parameters)

    ν     = model.atmosphere.spec_bands[iBand]
    ν̃    = mean(ν)
    i_ref = argmin(abs.(ν .- ν̃))
    effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
    n2, o2 = InelasticScattering.getRamanAtmoConstants(ν̃, effT)

    nPol  = parameters.polarization_type.n
    nSpec = length(ν)
    F₀   = zeros(FT, nPol, nSpec)
    SIF₀ = zeros(FT, nPol, nSpec)

    Tsolar_interp = if isfile(solar_file)
        t = solar_transmission_from_file(solar_file)
        using_interp(t[4:end, 1], t[4:end, 2])
    else
        ν -> 1.0
    end
    T_sun = 5777.0
    P = planck_spectrum_wn(T_sun, ν)
    for i in 1:length(P)
        F₀[1, i] = Tsolar_interp(ν[i]) * P[i]
    end

    RS_rrs = InelasticScattering.RRS(
        n2 = n2, o2 = o2,
        greek_raman = InelasticScattering.GreekCoefs([FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)]),
        fscattRayl  = [FT(1)], ϖ_Cabannes = [FT(1)],
        ϖ_λ₁λ₀      = zeros(FT, 1),
        i_λ₁λ₀      = zeros(Int, 1),
        Z⁻⁺_λ₁λ₀    = zeros(FT, 1, 1), Z⁺⁺_λ₁λ₀ = zeros(FT, 1, 1),
        i_ref       = i_ref, n_Raman = 0,
        F₀ = F₀, SIF₀ = SIF₀,
    )
    CoreRT.getRamanSSProp!(RS_rrs, 1e7 / ν̃, ν)

    R, T, ieR, ieT = CoreRT.rt_run_test(RS_rrs, model, iBand)

    RS_noRS = InelasticScattering.noRS{FT}(F₀ = F₀, SIF₀ = SIF₀)
    RnoRS, TnoRS, _, _ = CoreRT.rt_run_test(RS_noRS, model, iBand)

    mkpath(out_dir)
    writedlm(joinpath(out_dir, "out_RRS_630nm_730nm.dat"),
             [ν R[1,1,:] R[1,2,:] ieR[1,1,:] ieR[1,2,:]])
    writedlm(joinpath(out_dir, "out_noRS_630nm_730nm.dat"),
             [ν RnoRS[1,1,:] RnoRS[1,2,:]])

    return (; ν, R, T, ieR, ieT, RnoRS, TnoRS)
end

# --- Script entry point ----------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__
    @info "prototype_VS_O2Aband: starting" VS_YAML_VS VS_YAML_RRS VS_OUT_DIR VS_SECTION_SPEC
    if VS_SECTION_SPEC == "all"
        for ns in 1:VS_NSECT_TOTAL
            run_vs_section(ns)
        end
    else
        ns = parse(Int, VS_SECTION_SPEC)
        run_vs_section(ns)
    end
    run_rrs_o2aband()
    @info "prototype_VS_O2Aband: done"
end
