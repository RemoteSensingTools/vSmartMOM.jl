# =============================================================================
# prototype_O2ABand_RRS.jl — O₂ A-band RRS + noRS forward driver (Phase 6
# port from sanghavi).
#
# Runs RRS (rotational-Raman scattering) and noRS (pure-elastic) forward
# models over the O₂ A-band using `O2Parameters2.yaml` and writes
# `out_ray_sza*_vza*_{rrs,nors}_MOM.dat` files suitable for comparison
# against literature reference spectra. The sanghavi original also had
# ~500 lines of inline Plots.jl analysis; those are dropped — data-only
# deliverable, plot downstream if wanted.
#
# Path / config ENV overrides:
#   O2A_YAML        (default test_parameters/O2Parameters2.yaml)
#   O2A_SOLAR_FILE  (default "" → unit solar transmission)
#   O2A_OUT_DIR     (default test/benchmarks/O2Aband_RRS_output)
# =============================================================================

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using vSmartMOM.SolarModel
using Statistics
using DelimitedFiles

const O2A_YAML       = get(ENV, "O2A_YAML",
                            joinpath(pkgdir(vSmartMOM), "test", "test_parameters",
                                     "O2Parameters2.yaml"))
const O2A_SOLAR_FILE = get(ENV, "O2A_SOLAR_FILE", "")
const O2A_OUT_DIR    = get(ENV, "O2A_OUT_DIR",
                            joinpath(pkgdir(vSmartMOM), "test", "benchmarks",
                                     "O2Aband_RRS_output"))

function _solar_interp(ν, file)
    isfile(file) || return ones(length(ν))
    t = solar_transmission_from_file(file)
    xs = t[4:end, 1]; ys = t[4:end, 2]
    # Linear interp with unity outside the table range.
    out = similar(ν, Float64)
    for (k, νk) in enumerate(ν)
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

function run_prototype_O2Aband_RRS(; yaml        = O2A_YAML,
                                      solar_file = O2A_SOLAR_FILE,
                                      out_dir    = O2A_OUT_DIR,
                                      iBand      = 1,
                                      Tsun       = 5777.0,
                                      FT         = Float64)
    parameters = parameters_from_yaml(yaml)
    model      = model_from_parameters(parameters)
    mkpath(out_dir)

    ν     = model.atmosphere.spec_bands[iBand]
    nSpec = length(ν)
    ν̃     = mean(ν)
    i_ref = argmin(abs.(ν .- ν̃))
    effT  = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
    n2, o2 = InelasticScattering.getRamanAtmoConstants(ν̃, effT)

    nPol  = parameters.polarization_type.n
    F₀   = zeros(FT, nPol, nSpec)
    SIF₀ = zeros(FT, nPol, nSpec)

    P       = planck_spectrum_wn(Tsun, ν)
    Tsolar  = _solar_interp(ν, solar_file)
    F₀[1, :] = Tsolar .* P

    # RRS forward
    RS_rrs = InelasticScattering.RRS(
        n2 = n2, o2 = o2,
        greek_raman = InelasticScattering.GreekCoefs([FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)]),
        fscattRayl = [FT(1)], ϖ_Cabannes = [FT(1)],
        ϖ_λ₁λ₀ = zeros(FT, 1), i_λ₁λ₀ = zeros(Int, 1),
        Z⁻⁺_λ₁λ₀ = zeros(FT, 1, 1), Z⁺⁺_λ₁λ₀ = zeros(FT, 1, 1),
        i_ref = i_ref, n_Raman = 0,
        F₀ = F₀, SIF₀ = SIF₀,
    )
    CoreRT.getRamanSSProp!(RS_rrs, 1e7 / ν̃, ν)
    R, T, ieR, ieT = CoreRT.rt_run_test(RS_rrs, model, iBand)

    # noRS forward (uses same F₀)
    RS_noRS = InelasticScattering.noRS{FT}(F₀ = F₀, SIF₀ = SIF₀)
    RnoRS, TnoRS, _, _ = CoreRT.rt_run_test(RS_noRS, model, iBand)

    # Infer label from parameters.sza / parameters.vza for filenames
    sza_i = Int(round(parameters.sza))
    vza_i = Int(round(first(parameters.vza)))
    tag   = "sza$(sza_i)_vza$(vza_i)"

    rrs_cols  = hcat(ν, R[1, 1, :], R[1, 2, :], R[1, 3, :],
                        ieR[1, 1, :], ieR[1, 2, :], ieR[1, 3, :])
    nors_cols = hcat(ν, RnoRS[1, 1, :], RnoRS[1, 2, :], RnoRS[1, 3, :])

    writedlm(joinpath(out_dir, "out_ray_$(tag)_rrs_MOM.dat"),  rrs_cols)
    writedlm(joinpath(out_dir, "out_ray_$(tag)_nors_MOM.dat"), nors_cols)

    @info "prototype_O2ABand_RRS: wrote" rrs = joinpath(out_dir, "out_ray_$(tag)_rrs_MOM.dat") nors = joinpath(out_dir, "out_ray_$(tag)_nors_MOM.dat")

    return (; ν, R, T, ieR, ieT, RnoRS, TnoRS, parameters, model)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_prototype_O2ABand_RRS()
end
