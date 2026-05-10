# =========================================================================
# O2 A-band RRS + SIF grid driver (sanghavi Phase 3 acceptance script)
# =========================================================================
#
# Ported from sanghavi `test/benchmarks/creategrid_O2Aband_RamanSIF.jl`
# (~930 lines) in Phase 3 of the sanghavi-unified merge. The original
# script generated a (isurf × ρ × sza × iBand) grid of forward-model
# outputs, each written to `/home/sanghavi/data/RamanSIFgrid_new/`-style
# hardcoded paths. This port is a minimal reproducible driver that
# exercises the RRS + SIF path once, with portable paths and env-var
# overrides, and dumps R / T / ieR / ieT / F₀ / ν to JLD2.
#
# Grid generation (looping over many scenarios) remains a research-facing
# concern and can be re-landed later on top of this driver as needed.
#
# Usage:
#   # default config: Phase1b narrow O2 A-band RRS, SIF on, 1 (sza, ρ)
#   julia --project=test test/benchmarks/creategrid_O2Aband_RamanSIF.jl
#
#   # different albedo / sza
#   SZA=45 ALBEDO=0.3 \
#     julia --project=test test/benchmarks/creategrid_O2Aband_RamanSIF.jl
#
#   # turn SIF off
#   SIF_ON=0 \
#     julia --project=test test/benchmarks/creategrid_O2Aband_RamanSIF.jl

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using vSmartMOM.SolarModel
using JLD2
using Statistics

# TODO: The (0.5π / maximum(J_SIF)) rescaling inside `load_sif_spectrum`
# is an intentional hack to make SIF magnitude data-independent for grid
# generation. This normalizes shape but discards absolute physical
# magnitude. Revisit: confirm the downstream physics depends only on SIF
# shape (not absolute flux), or replace with physical units
# (mW/m²/cm⁻¹). Not a merge blocker. (Per plan §Phase 3 amendments §2.5.)

# --- Configuration ---------------------------------------------------------

const GRID_YAML = get(ENV, "GRID_YAML",
    joinpath(pkgdir(vSmartMOM), "test", "test_parameters",
             "Phase1b_RRS_761-764nm.yaml"))
const GRID_OUTDIR = get(ENV, "GRID_OUTDIR",
    joinpath(pkgdir(vSmartMOM), "test", "benchmarks", "aerosol_exploration_output",
             "grid_O2Aband_RamanSIF"))
const SZA_OVERRIDE = tryparse(Float64, get(ENV, "SZA", ""))
const ALBEDO_OVERRIDE = tryparse(Float64, get(ENV, "ALBEDO", ""))
const SIF_ON = get(ENV, "SIF_ON", "1") != "0"

mkpath(GRID_OUTDIR)
@info "O2 A-band RRS+SIF driver: YAML $(GRID_YAML)"
@info "O2 A-band RRS+SIF driver: outdir $(GRID_OUTDIR)"
@info "O2 A-band RRS+SIF driver: SIF $(SIF_ON ? "on" : "off")"

# --- Model build -----------------------------------------------------------

parameters = parameters_from_yaml(GRID_YAML)
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

# --- Per-band RRS + SIF run -----------------------------------------------

iBand = 1
ν     = parameters.spec_bands[iBand]
nSpec = length(ν)
ν̃     = mean(ν)
i_ref = argmin(abs.(ν .- ν̃))
effT  = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
n2, o2 = InelasticScattering.getRamanAtmoConstants(ν̃, effT)

# Solar spectrum (Planck at 5777 K, no solar-line absorption since solar.out
# isn't bundled — see prototype_EMIT_aer_ht.jl for the SOLAR_OUT override).
T_sun = FT(5777)
P     = planck_spectrum_wn(T_sun, ν) .* FT(2.1629e-05 * π)  # → mW/m²/cm⁻¹
F₀    = zeros(FT, nPol, nSpec); F₀[1, :] .= FT.(P)

# RRS type with Raman cross-wavelength operators
rs_rrs = InelasticScattering.RRS(
    n2 = n2, o2 = o2,
    greek_raman = InelasticScattering.GreekCoefs(
        [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)]),
    fscattRayl  = [FT(1)], ϖ_Cabannes = [FT(1)],
    ϖ_λ₁λ₀ = zeros(FT, 1), i_λ₁λ₀ = zeros(Int, 1),
    Z⁻⁺_λ₁λ₀ = zeros(FT, 1, 1), Z⁺⁺_λ₁λ₀ = zeros(FT, 1, 1),
    i_ref = i_ref, n_Raman = 0,
    F₀ = F₀, SIF₀ = zeros(FT, nPol, nSpec))
CoreRT.getRamanSSProp!(rs_rrs, 1e7/ν̃, ν)

# Populate SIF₀ via the loader (normalized to peak = 0.5π by default).
if SIF_ON
    ν_sif, jSIF = load_sif_spectrum()
    build_sif_source(rs_rrs, collect(ν), ν_sif, jSIF; pol_component=1)
    @info "SIF₀ loaded: peak = $(maximum(rs_rrs.SIF₀[1, :]))  length = $(size(rs_rrs.SIF₀, 2))"
end

# Parallel noRS run for delta-vs-SIF comparison
rs_norrs = InelasticScattering.noRS{FT}(
    fscattRayl  = [FT(1)], ϖ_Cabannes = [FT(1)],
    bandSpecLim = UnitRange{Int}[], iBand = [iBand],
    F₀ = copy(F₀),
    SIF₀ = SIF_ON ? copy(rs_rrs.SIF₀) : zeros(FT, nPol, nSpec))

R_rrs, T_rrs, ieR, ieT = CoreRT.rt_run_test(rs_rrs, model, iBand)
R_nors, T_nors, _, _    = CoreRT.rt_run_test(rs_norrs, model, iBand)

# --- Save outputs ----------------------------------------------------------

sza_str = replace(string(round(parameters.sza, digits=1)), "." => "p")
alb = parameters.brdf[iBand].albedo
alb_str = replace(string(round(alb, digits=2)), "." => "p")
sif_tag = SIF_ON ? "sifon" : "sifoff"
outpath = joinpath(GRID_OUTDIR, "rrs_grid_sza$(sza_str)_alb$(alb_str)_$(sif_tag).jld2")

jldopen(outpath, "w") do f
    f["ν"]        = collect(ν)
    f["F₀"]       = collect(F₀)
    f["SIF₀"]     = collect(rs_rrs.SIF₀)
    f["R_rrs"]    = collect(R_rrs)
    f["T_rrs"]    = collect(T_rrs)
    f["ieR_rrs"]  = collect(ieR)
    f["ieT_rrs"]  = collect(ieT)
    f["R_nors"]   = collect(R_nors)
    f["T_nors"]   = collect(T_nors)
    f["sza"]      = parameters.sza
    f["albedo"]   = alb
    f["sif_on"]   = SIF_ON
end

@info "O2 A-band RRS+SIF driver: wrote $outpath"
