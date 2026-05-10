# =================================================================
# EMIT Aerosol Height Sensitivity Prototype — Linearized RT
# =================================================================
# Ported from the sanghavi branch:
#   test/benchmarks/prototype_EMIT_aer_ht.jl
#
# This script runs the linearized radiative transfer model across
# EMIT spectral bands (O2 UV/B/A, CO2, CH4), producing both
# reflectance/transmittance (R, T) and their analytical Jacobians
# (Ṙ, Ṫ) with respect to aerosol and gas properties.
#
# The high-res spectra are then convolved to EMIT's instrument
# resolution (~7 nm sampling) using a Gaussian kernel.
#
# Every change from the original script is marked with [CHANGED]
# to help the original author (VS) follow the migration.
# =================================================================

##
# [CHANGED] Old: using CUDA; device!(1)
# New: CUDA is optional. Uncomment the two lines below for GPU runs.
#      The YAML sets `architecture: default_architecture` which auto-detects GPU.
# using CUDA
# CUDA.device!(0)  # Select GPU device (0-indexed in CUDA.jl)

using Revise
using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.SolarModel
using vSmartMOM.InelasticScattering
using Statistics
using Interpolations
using Distributions
using DelimitedFiles
using LinearAlgebra

# [CHANGED] Old: include("/home/sanghavi/code/github/vSmartMOM.jl/src/Testing/perturb_parameters.jl")
# New: Locate perturb_parameters.jl portably via pathof(vSmartMOM).
#      This file defines `perturb_parameters(params, pert_pct)` and
#      `compute_FD_modelJacobian(params, pert_pct)` for finite-difference
#      Jacobian validation.
include(joinpath(dirname(pathof(vSmartMOM)), "Testing", "perturb_parameters.jl"))

# RT mode types — these are defined in src/vSmartMOM.jl
fwd_mode = FwdMode()
lin_mode = LinMode()
FT = Float64

## =================================================================
# Load parameters and build models
# =================================================================

# [CHANGED] Old: parameters_from_yaml("/home/sanghavi/code/github/vSmartMOM.jl/test/test_parameters/ParamsEMIT.yaml")
# New: Use @__DIR__ to locate the YAML relative to this script.
parameters = parameters_from_yaml(joinpath(@__DIR__, "ParamsEMIT.yaml"))

pert_pct = FT(0.01)
pert_parameters = perturb_parameters(parameters, pert_pct)

# [UNCHANGED] Build both the forward model and the linearized model.
# model_from_parameters(lin_mode, params) returns (model, lin_model):
#   model    — standard vSmartMOM_Model with precomputed optical properties
#   lin_model — vSmartMOM_Lin with derivative fields (τ̇_abs, τ̇_aer, lin_aerosol_optics)
model, lin_model = model_from_parameters(lin_mode, parameters)

@show "Model built successfully"

## =================================================================
# Solar irradiance setup
# =================================================================

# [CHANGED] Old: Tsolar = solar_transmission_from_file("/home/sanghavi/.../solar.out")
# New: Use default_solar_transmission() which auto-downloads the Toon solar
#      line-list if not already cached, then returns [ν, transmission].
#      We only need the raw file for manual interpolation (matching the
#      original approach of skipping the first 3 header rows).
#
# Alternative (simpler): use default_solar_transmission(ν) directly per band.
# We keep the manual interpolation approach here to match the original workflow.
solar_file = joinpath(dirname(pathof(vSmartMOM)), "SolarModel", "solar.out")
if !isfile(solar_file)
    download("http://web.gps.caltech.edu/~cfranken/hitran_2016/solar_merged_20160127_600_26316_100.out",
             solar_file)
end
Tsolar = solar_transmission_from_file(solar_file)
Tsolar_interp = LinearInterpolation(Tsolar[4:end, 1], Tsolar[4:end, 2])

T_sun = FT(5777.0)  # Solar effective temperature (K)

## =================================================================
# Instrument convolution kernel (EMIT-like Gaussian)
# =================================================================
# [CHANGED] Old: used InstrumentOperator.create_instrument_kernel + ImageFiltering.imfilter
# New: Simple 1D Gaussian convolution using base Julia (no extra dependencies).
#      The kernel σ = 12.5 cm⁻¹ approximates EMIT's spectral response.

x = collect(-40:0.05:40)
kernel_raw = exp.(-x.^2 / (2 * 12.5^2))
kernel = kernel_raw / sum(kernel_raw)

"""
    gauss_conv(signal, kernel)

1D convolution of `signal` with `kernel`, returning same-length output (centered).
Replaces ImageFiltering.imfilter for this use case.
"""
function gauss_conv(signal::AbstractVector, kernel::AbstractVector)
    n = length(signal)
    m = length(kernel)
    hm = m ÷ 2
    out = similar(signal)
    for i in 1:n
        s = zero(eltype(signal))
        w = zero(eltype(kernel))
        for j in 1:m
            idx = i - hm + j - 1
            if 1 <= idx <= n
                s += signal[idx] * kernel[j]
                w += kernel[j]
            end
        end
        out[i] = s / w
    end
    return out
end

## =================================================================
# Main loop: per-band linearized RT + convolution
# =================================================================

I_conv = []
λEMIT  = []
I_EMIT = []
λ_all  = []
Δλ = FT(7.0)  # EMIT sampling interval (nm)

for iBand in 1:length(model.params.spec_bands)
    @show iBand

    # Wavenumber grid for this band
    ν = model.params.spec_bands[iBand]
    ν̃ = mean(ν)

    # EMIT wavelength grid for this band (nm)
    push!(λEMIT, 1e7/ν[end]:Δλ:1e7/ν[1])

    # Reference index for Raman (not used for noRS, but kept for compatibility)
    i_ref = argmin(abs.(ν .- ν̃))

    # Effective temperature for Raman calculations (unused for noRS)
    effT = FT(300.0)

    # [UNCHANGED] Construct a noRS (no Raman scattering) object.
    # The noRS struct now includes a SIF₀ field (added in unified branch)
    # but the constructor defaults work fine.
    RS_type = InelasticScattering.noRS(
        fscattRayl  = [FT(1)],
        ϖ_Cabannes  = [FT(1)],
        bandSpecLim = [],
        iBand       = [1],
        F₀          = zeros(FT, 1, 1),
        SIF₀        = zeros(FT, 1, 1))

    # [UNCHANGED] Compute Planck blackbody spectrum at solar temperature
    P = planck_spectrum_wn(T_sun, ν)

    # Build the solar irradiance Stokes vector: F₀ = Planck × solar_transmission
    F₀ = zeros(FT, length(P))
    RS_type.F₀ = zeros(FT, model.params.polarization_type.n, length(P))
    for i in 1:length(P)
        sol_trans = Tsolar_interp(ν[i])
        F₀[i] = FT(sol_trans * P[i])
        RS_type.F₀[1, i] = F₀[i]
    end

    # [UNCHANGED] Count Jacobian parameters:
    #   NAer = number of aerosol types (each contributes 7 params: τ,nᵣ,nᵢ,μ,σ,p₀,σp)
    #   NGas = 1 (water vapor q) + number of variable gas molecules
    #   NSurf = number of surface albedo parameters (one per band)
    NAer  = length(model.params.scattering_params.rt_aerosols)
    NGas  = 1 + length(model.params.absorption_params.variable_molecules)
    NSurf = length(model.params.spec_bands)
    Nparams = NAer * 7 + NGas + NSurf

    # [UNCHANGED] Run the linearized RT.
    # Returns: R (reflectance), T (transmittance), Ṙ (dR/dx), Ṫ (dT/dx)
    # R dimensions: (nVZA, nStokes, nSpec)
    # Ṙ dimensions: (nVZA, nStokes, nSpec, Nparams)
    R, T, Ṙ, Ṫ = CoreRT.rt_run_test(RS_type, model, lin_model, NAer, NGas, NSurf, iBand)

    # === Convolution to EMIT instrument grid ===
    # Convert radiance from mW/m²-str-cm⁻¹ to mW/m²-str-nm
    convfct = 1e7 ./ ν.^2
    push!(I_conv, gauss_conv(R[1, 1, :] .* convfct, kernel))

    # Sample convolved spectrum at EMIT wavelength grid
    tmpI_EMIT = []
    for j in 1:length(λEMIT[iBand])
        ii = findall(x -> 1e7/(λEMIT[iBand][j] + Δλ/2) < x < 1e7/(λEMIT[iBand][j] - Δλ/2), ν)
        tmpI = 0.0
        δλ   = 0.0
        for k in 1:length(ii)
            if ii[k] == 1
                dλ = (1e7/ν[ii[k]] - 1e7/ν[ii[k]+1])
            elseif ii[k] == length(ν)
                dλ = (1e7/ν[ii[k]-1] - 1e7/ν[ii[k]])
            else
                dλ = 0.5 * (1e7/ν[ii[k]-1] - 1e7/ν[ii[k]+1])
            end
            tmpI += I_conv[iBand][ii[k]] * dλ
            δλ   += dλ
        end
        push!(tmpI_EMIT, tmpI / δλ)
    end
    push!(I_EMIT, tmpI_EMIT)
    push!(λ_all, 1e7 ./ model.params.spec_bands[iBand])
end

@show "All bands completed"

## =================================================================
# Results summary
# =================================================================
# At this point:
#   I_conv[iBand]  — convolved hi-res radiance per band (mW/m²-str-nm)
#   I_EMIT[iBand]  — sampled EMIT-resolution radiance per band
#   λEMIT[iBand]   — EMIT wavelength grid (nm)
#   λ_all[iBand]   — hi-res wavelength grid (nm)
#   R, T, Ṙ, Ṫ    — last band's full RT output + Jacobians
#
# The Jacobian Ṙ has shape (nVZA, nStokes, nSpec, Nparams) where
# Nparams = NAer*7 + NGas + NSurf. The parameter ordering is:
#   [q, var_mol_1, ..., var_mol_N,
#    aer1_τ, aer1_nᵣ, aer1_nᵢ, aer1_μ, aer1_σ, aer1_p₀, aer1_σp,
#    ...,
#    surf_albedo_band1, ..., surf_albedo_bandN]

## =================================================================
# Post-analysis / Plotting (stub)
# =================================================================
# The original script (lines ~100-826) reads pre-computed results from
# /home/sanghavi/EMIT/ for multiple scenarios:
#   - 3 surface types: arid, ocean, semiarid
#   - AOD values: 0, 0.5, 5.0
#   - Aerosol heights: 500 hPa, 900 hPa
#   - CH4 levels: high, low
#
# These files are not available in this repository. To reproduce:
# 1. Run this script with different YAML configs (varying τ_ref, p₀, CH4 VMR)
# 2. Save results with writedlm() to local output directory
# 3. Use Plots.jl to create comparison figures
#
# Example of original data-reading pattern:
#   tmp = readdlm("path/to/I_EMIT_aer0p0_hiCH4_band1.dat")
#
# See the full original script at:
#   https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/sanghavi/test/benchmarks/prototype_EMIT_aer_ht.jl
