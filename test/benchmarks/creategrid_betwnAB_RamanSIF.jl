# =========================================================================
# Between O₂ A/B bands RRS + SIF grid driver (sanghavi Phase 6 port).
# =========================================================================
#
# Mirror of `creategrid_O2Aband_RamanSIF.jl` but targets the spectral
# region BETWEEN the O₂ A and B bands (~700-740 nm) via
# `O2_parameters1p5_SIF_grid.yaml`. Useful for retrieval tests that need
# aerosol/Rayleigh info from a minimal-absorption window adjacent to the
# SIF emission peak.
#
# Usage:
#   julia --project=test test/benchmarks/creategrid_betwnAB_RamanSIF.jl
#   SZA=45 ALBEDO=0.3 SIF_ON=0 julia --project=test test/benchmarks/creategrid_betwnAB_RamanSIF.jl
# =========================================================================

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using vSmartMOM.SolarModel
using JLD2
using Statistics

include(joinpath(pkgdir(vSmartMOM), "src", "SIF_emission", "sif_loader.jl"))

const GRID_YAML_BTWN = get(ENV, "GRID_YAML_BTWN",
    joinpath(pkgdir(vSmartMOM), "test", "test_parameters",
             "O2_parameters1p5_SIF_grid.yaml"))
const GRID_OUTDIR_BTWN = get(ENV, "GRID_OUTDIR_BTWN",
    joinpath(pkgdir(vSmartMOM), "test", "benchmarks", "aerosol_exploration_output",
             "grid_betwnAB_RamanSIF"))
const SZA_OVERRIDE    = tryparse(Float64, get(ENV, "SZA", ""))
const ALBEDO_OVERRIDE = tryparse(Float64, get(ENV, "ALBEDO", ""))
const SIF_ON          = get(ENV, "SIF_ON", "1") != "0"

mkpath(GRID_OUTDIR_BTWN)
@info "Between-AB RRS+SIF driver" GRID_YAML_BTWN GRID_OUTDIR_BTWN SIF_ON

parameters = parameters_from_yaml(GRID_YAML_BTWN)
if SZA_OVERRIDE !== nothing
    parameters.sza = SZA_OVERRIDE
end
if ALBEDO_OVERRIDE !== nothing
    for iB in 1:length(parameters.brdf)
        parameters.brdf[iB].albedo = parameters.float_type(ALBEDO_OVERRIDE)
    end
end
FT    = parameters.float_type
model = model_from_parameters(parameters)
nPol  = parameters.polarization_type.n

iBand = 1
ν     = parameters.spec_bands[iBand]
nSpec = length(ν)
ν̃     = mean(ν)
i_ref = argmin(abs.(ν .- ν̃))
effT  = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
n2, o2 = InelasticScattering.getRamanAtmoConstants(ν̃, effT)

T_sun = FT(5777)
P     = planck_spectrum_wn(T_sun, ν) .* FT(2.1629e-05 * π)
F₀    = zeros(FT, nPol, nSpec); F₀[1, :] .= FT.(P)

rs_rrs = InelasticScattering.RRS(
    n2 = n2, o2 = o2,
    greek_raman = InelasticScattering.GreekCoefs(
        [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)]),
    fscattRayl  = [FT(1)], ϖ_Cabannes = [FT(1)],
    ϖ_λ₁λ₀ = zeros(FT, 1), i_λ₁λ₀ = zeros(Int, 1),
    Z⁻⁺_λ₁λ₀ = zeros(FT, 1, 1), Z⁺⁺_λ₁λ₀ = zeros(FT, 1, 1),
    i_ref = i_ref, n_Raman = 0,
    F₀ = F₀, SIF₀ = zeros(FT, nPol, nSpec))
CoreRT.getRamanSSProp!(rs_rrs, 1e7 / ν̃, ν)

if SIF_ON
    ν_sif, jSIF = load_sif_spectrum()
    build_sif_source(rs_rrs, collect(ν), ν_sif, jSIF; pol_component = 1)
end

rs_norrs = InelasticScattering.noRS{FT}(
    fscattRayl = [FT(1)], ϖ_Cabannes = [FT(1)],
    bandSpecLim = UnitRange{Int}[], iBand = [iBand],
    F₀ = copy(F₀),
    SIF₀ = SIF_ON ? copy(rs_rrs.SIF₀) : zeros(FT, nPol, nSpec))

R_rrs, T_rrs, ieR, ieT = CoreRT.rt_run_test(rs_rrs, model, iBand)
R_nors, T_nors, _, _   = CoreRT.rt_run_test(rs_norrs, model, iBand)

sza_str = replace(string(round(parameters.sza, digits = 1)), "." => "p")
alb     = parameters.brdf[iBand].albedo
alb_str = replace(string(round(alb, digits = 2)), "." => "p")
sif_tag = SIF_ON ? "sifon" : "sifoff"
outpath = joinpath(GRID_OUTDIR_BTWN, "rrs_grid_betwnAB_sza$(sza_str)_alb$(alb_str)_$(sif_tag).jld2")

jldopen(outpath, "w") do f
    f["ν"]       = collect(ν)
    f["F₀"]      = collect(F₀)
    f["SIF₀"]    = collect(rs_rrs.SIF₀)
    f["R_rrs"]   = collect(R_rrs)
    f["T_rrs"]   = collect(T_rrs)
    f["ieR_rrs"] = collect(ieR)
    f["ieT_rrs"] = collect(ieT)
    f["R_nors"]  = collect(R_nors)
    f["T_nors"]  = collect(T_nors)
    f["sza"]     = parameters.sza
    f["albedo"]  = alb
    f["sif_on"]  = SIF_ON
    f["band"]    = "btwnAB"
end

@info "Between-AB RRS+SIF driver: wrote $outpath"
