# =========================================================================
# EMIT-style noRS forward + linearized RT acceptance script
# =========================================================================
#
# Ported from sanghavi `test/benchmarks/prototype_EMIT_aer_ht.jl` (Phase 3
# of the sanghavi-unified merge). Exercises the full forward + linearized
# rt_run path with aerosols + Lambertian surface + O2/CH4 line-by-line
# absorption.
#
# The hundreds of lines of post-processing/visualization from the
# sanghavi version (writing to `/home/sanghavi/EMIT/...`, multi-scenario
# cross-comparison plots) were research detritus and are intentionally
# omitted. This script is now a minimal, reproducible forward+Jacobian
# driver that writes R, T, Ṙ, Ṫ, λ_EMIT outputs to JLD2 for downstream
# analysis.
#
# Usage:
#   # default reduced-range config (<1 min on CPU)
#   julia --project=test test/benchmarks/prototype_EMIT_aer_ht.jl
#
#   # full-range config (requires network-mounted O2 LUT)
#   EMIT_YAML=test_parameters/ParamsEMIT.yaml \
#     julia --project=test test/benchmarks/prototype_EMIT_aer_ht.jl
#
#   # write outputs under a custom dir
#   EMIT_OUTDIR=/tmp/emit_run \
#     julia --project=test test/benchmarks/prototype_EMIT_aer_ht.jl

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.Architectures
using vSmartMOM.SolarModel
using vSmartMOM.InelasticScattering
using DataInterpolations
using DelimitedFiles
using JLD2
using Statistics

# Perturb-parameter helper for finite-difference Jacobian comparisons
include(joinpath(pkgdir(vSmartMOM), "src", "Testing", "perturb_parameters.jl"))

# --- Configuration ---------------------------------------------------------

const EMIT_YAML = get(ENV, "EMIT_YAML",
    joinpath(pkgdir(vSmartMOM), "test", "test_parameters", "ParamsEMIT_fast.yaml"))
const EMIT_OUTDIR = get(ENV, "EMIT_OUTDIR",
    joinpath(pkgdir(vSmartMOM), "test", "benchmarks", "aerosol_exploration_output",
             "emit_aer_ht_run"))

mkpath(EMIT_OUTDIR)
@info "EMIT acceptance: loading YAML $(EMIT_YAML)"
@info "EMIT acceptance: writing outputs to $(EMIT_OUTDIR)"

# --- Forward + linearized model build -------------------------------------

parameters = parameters_from_yaml(EMIT_YAML)
pert_pct   = parameters.float_type(0.01)
pert_parameters = perturb_parameters(parameters, pert_pct)
FT = parameters.float_type

model, lin_model = model_from_parameters(LinMode(), parameters)

# --- Solar spectrum prep --------------------------------------------------

T_sun     = FT(5777)

# Solar transmission file path: defaults to the 32 MB `solar.out` if present
# under `src/SolarModel/` (sanghavi convention), else falls back to a unit
# transmission (i.e. pure Planck-scaled F₀). The large file is not bundled
# in the unified repo to keep it lean; set SOLAR_OUT to point at a local
# copy if needed.
const _solar_out_path = get(ENV, "SOLAR_OUT",
    joinpath(pkgdir(vSmartMOM), "src", "SolarModel", "solar.out"))
Tsolar_interp = if isfile(_solar_out_path)
    solar_dat = solar_transmission_from_file(_solar_out_path)
    LinearInterpolation(solar_dat[4:end, 2], solar_dat[4:end, 1];
                        extrapolation=ExtrapolationType.Linear)
else
    @warn "solar.out not found at $(_solar_out_path); using unit solar transmission"
    ν -> one(FT)
end

# --- Per-band forward + linearized run ------------------------------------

spec_bands = parameters.spec_bands
n_bands    = length(spec_bands)
nPol       = parameters.polarization_type.n

# Band range — defaults to all bands declared in the YAML. Override via
# `EMIT_BANDS=1:1` (or any range) to restrict.
const EMIT_BANDS = let s = get(ENV, "EMIT_BANDS", "1:$(n_bands)")
    rng_parts = parse.(Int, split(s, ':'))
    length(rng_parts) == 2 ? (rng_parts[1]:rng_parts[2]) :
    length(rng_parts) == 1 ? (rng_parts[1]:rng_parts[1]) :
    error("EMIT_BANDS must be \"a\" or \"a:b\", got $(s)")
end

allband_R = Dict{Int, Any}()
allband_T = Dict{Int, Any}()
allband_Ṙ = Dict{Int, Any}()
allband_Ṫ = Dict{Int, Any}()
allband_ν = Dict{Int, Any}()

for iBand in EMIT_BANDS
    ν     = spec_bands[iBand]
    nSpec = length(ν)

    # Per-band solar irradiance scaled by Planck spectrum at T_sun.
    P  = planck_spectrum_wn(T_sun, ν)
    F₀ = zeros(FT, nPol, nSpec)
    for i in 1:nSpec
        F₀[1, i] = FT(Tsolar_interp(ν[i]) * P[i])
    end

    rs = InelasticScattering.noRS{FT}(
        fscattRayl  = [FT(1)],
        ϖ_Cabannes  = [FT(1)],
        bandSpecLim = UnitRange{Int}[],
        iBand       = [iBand],
        F₀          = F₀,
        SIF₀        = zeros(FT, nPol, nSpec))

    NAer  = length(parameters.scattering_params.rt_aerosols)
    NGas  = 1 + length(parameters.absorption_params.variable_molecules[iBand])
    NSurf = 1   # one Lambertian albedo parameter per band

    @info "Band $iBand: nSpec=$nSpec  NAer=$NAer  NGas=$NGas  NSurf=$NSurf"
    R, T, Ṙ, Ṫ = CoreRT.rt_run_test(rs, model, lin_model, NAer, NGas, NSurf, iBand)

    allband_ν[iBand] = collect(ν)
    allband_R[iBand] = R
    allband_T[iBand] = T
    allband_Ṙ[iBand] = Ṙ
    allband_Ṫ[iBand] = Ṫ
end

# --- Save outputs ----------------------------------------------------------

jldopen(joinpath(EMIT_OUTDIR, "rt_outputs.jld2"), "w") do f
    for iBand in sort(collect(keys(allband_R)))
        band_grp = "band$(iBand)"
        f["$(band_grp)/ν"] = allband_ν[iBand]
        f["$(band_grp)/R"] = allband_R[iBand]
        f["$(band_grp)/T"] = allband_T[iBand]
        f["$(band_grp)/Rdot"] = allband_Ṙ[iBand]
        f["$(band_grp)/Tdot"] = allband_Ṫ[iBand]
    end
end

@info "EMIT acceptance: wrote $(joinpath(EMIT_OUTDIR, "rt_outputs.jld2"))"
