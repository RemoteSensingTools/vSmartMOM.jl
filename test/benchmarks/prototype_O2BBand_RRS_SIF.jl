# =============================================================================
# prototype_O2BBand_RRS_SIF.jl — O₂ B-band RRS + noRS with SIF emission
# (Phase 6 port from sanghavi).
#
# Same structure as prototype_O2ABand_RRS_SIF.jl but uses
# O2Parameters1_SIF.yaml (O₂ B-band region ~686 nm). Outputs have `_BBO2`
# suffix for easy distinguishing from A-band results. The sanghavi
# original dedicated ~700 lines to inline Plots analysis; those are
# dropped in favor of data-file output.
#
# ENV overrides:
#   O2B_SIF_YAML       — test YAML (default O2Parameters1_SIF.yaml)
#   O2B_SIF_OUT_DIR    — output dir (default test/benchmarks/O2Bband_RRS_SIF_output)
#   O2B_SIF_SOLAR_FILE — solar transmission file; "" → unity
# =============================================================================

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using vSmartMOM.SolarModel
using Statistics
using DelimitedFiles

include(joinpath(pkgdir(vSmartMOM), "src", "SIF_emission", "sif_loader.jl"))

const O2B_SIF_YAML        = get(ENV, "O2B_SIF_YAML",
                                 joinpath(pkgdir(vSmartMOM), "test", "test_parameters",
                                          "O2Parameters1_SIF.yaml"))
const O2B_SIF_OUT_DIR     = get(ENV, "O2B_SIF_OUT_DIR",
                                 joinpath(pkgdir(vSmartMOM), "test", "benchmarks",
                                          "O2Bband_RRS_SIF_output"))
const O2B_SIF_SOLAR_FILE  = get(ENV, "O2B_SIF_SOLAR_FILE", "")

function _o2b_solar_transm(ν, file)
    isfile(file) || return ones(length(ν))
    t = solar_transmission_from_file(file)
    xs = t[4:end, 1]; ys = t[4:end, 2]
    out = similar(ν, Float64)
    @inbounds for (k, νk) in enumerate(ν)
        if νk ≤ xs[1] || νk ≥ xs[end]
            out[k] = 1.0
        else
            j = searchsortedfirst(xs, νk)
            tt = (νk - xs[j - 1]) / (xs[j] - xs[j - 1])
            out[k] = ys[j - 1] + tt * (ys[j] - ys[j - 1])
        end
    end
    return out
end

function _o2b_sif_spectrum_interp(ν)
    ν_sif, jSIF = load_sif_spectrum()
    xs = collect(ν_sif); ys = collect(jSIF)
    out = similar(ν, Float64)
    @inbounds for (k, νk) in enumerate(ν)
        if νk ≤ xs[1] || νk ≥ xs[end]
            out[k] = 0.0
        else
            j = searchsortedfirst(xs, νk)
            tt = (νk - xs[j - 1]) / (xs[j] - xs[j - 1])
            out[k] = ys[j - 1] + tt * (ys[j] - ys[j - 1])
        end
    end
    return out
end

function run_prototype_O2Bband_RRS_SIF(; yaml       = O2B_SIF_YAML,
                                          solar_file = O2B_SIF_SOLAR_FILE,
                                          out_dir    = O2B_SIF_OUT_DIR,
                                          iBand      = 1,
                                          FT         = Float64,
                                          Tsun       = 5777.0)
    parameters = parameters_from_yaml(yaml)
    model      = model_from_parameters(parameters)
    mkpath(out_dir)

    ν     = model.atmosphere.spec_bands[iBand]
    nSpec = length(ν)
    ν̃     = mean(ν)
    i_ref = argmin(abs.(ν .- ν̃))
    effT  = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
    n2, o2 = InelasticScattering.getRamanAtmoConstants(ν̃, effT)

    nPol = parameters.polarization_type.n

    P       = FT.(planck_spectrum_wn(Tsun, ν) .* 2.1629e-5 * π)
    Ts      = FT.(_o2b_solar_transm(ν, solar_file))
    F₀mat   = zeros(FT, nPol, nSpec); F₀mat[1, :] = Ts .* P

    SIF_values = FT.(_o2b_sif_spectrum_interp(ν))
    SIF_mat    = zeros(FT, nPol, nSpec); SIF_mat[1, :] = SIF_values

    RS_rrs = InelasticScattering.RRS(
        n2 = n2, o2 = o2,
        greek_raman = InelasticScattering.GreekCoefs([FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)]),
        fscattRayl = [FT(1)], ϖ_Cabannes = [FT(1)],
        ϖ_λ₁λ₀ = zeros(FT, 1), i_λ₁λ₀ = zeros(Int, 1),
        Z⁻⁺_λ₁λ₀ = zeros(FT, 1, 1), Z⁺⁺_λ₁λ₀ = zeros(FT, 1, 1),
        i_ref = i_ref, n_Raman = 0,
        F₀ = F₀mat, SIF₀ = zeros(FT, nPol, nSpec),
    )
    CoreRT.getRamanSSProp!(RS_rrs, 1e7 / ν̃, ν)
    R_rrs, T_rrs, ieR, ieT = CoreRT.rt_run_test(RS_rrs, model, iBand)

    RS_noRS = InelasticScattering.noRS{FT}(F₀ = F₀mat, SIF₀ = zeros(FT, nPol, nSpec))
    R_nors, T_nors, _, _ = CoreRT.rt_run_test(RS_noRS, model, iBand)

    RS_noRS_sif = InelasticScattering.noRS{FT}(F₀ = F₀mat, SIF₀ = SIF_mat)
    R_nors_sif, T_nors_sif, _, _ = CoreRT.rt_run_test(RS_noRS_sif, model, iBand)

    sza_i = Int(round(parameters.sza))
    vza_i = Int(round(first(parameters.vza)))
    tag   = "sza$(sza_i)_vza$(vza_i)"

    writedlm(joinpath(out_dir, "rayl_$(tag)_rrs_BBO2.dat"),
             [ν R_rrs[1, 1, :] R_rrs[1, 2, :] R_rrs[1, 3, :] ieR[1, 1, :] ieR[1, 2, :] ieR[1, 3, :]])
    writedlm(joinpath(out_dir, "rayl_$(tag)_nors_BBO2.dat"),
             [ν R_nors[1, 1, :] R_nors[1, 2, :] R_nors[1, 3, :]])
    writedlm(joinpath(out_dir, "rayl_$(tag)_nors_BBO2_SIF.dat"),
             [ν R_nors_sif[1, 1, :] R_nors_sif[1, 2, :] R_nors_sif[1, 3, :]])

    @info "prototype_O2BBand_RRS_SIF: wrote" out_dir tag

    return (; ν, R_rrs, T_rrs, ieR, ieT, R_nors, T_nors, R_nors_sif, T_nors_sif)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_prototype_O2Bband_RRS_SIF()
end
