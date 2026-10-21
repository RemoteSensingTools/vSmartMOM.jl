## =============================================================================
## lake_mary_vsmartmom.jl  (uni_vSmartMOM, v2.1.0)
##
## Reconstruct EMIT TOA radiance over Lake Mary, CA (snow scene, 27 Mar 2025)
## via DIRECT vSmartMOM forward simulation.
##
## Base config: test/test_parameters/ParamsEMIT_MODTRANcomp.yaml
##   - Lake Mary geometry: SZA=39.27, SAA=212.97, VZA=5.23, VAA=48.49
##   - Surface:  LambertianSurfaceSpectrum from Niklas's in-situ ρ_s
##   - Aerosol:  τ_ref = 0.04 (AOD550)
##   - H₂O:      scale q to column ≈ 0.3 g/cm²
##
## Run:
##   VSMARTMOM_EMIT_LUT_DIR=/home/sanghavi/data/EMIT_MODTRANcomp/LUTs \
##     JULIA_LOAD_PATH="@:/home/sanghavi/.julia/environments/v1.11:@stdlib" \
##     julia --project=/home/sanghavi/code/github/uni_vSmartMOM \
##           test/benchmarks/lake_mary_vsmartmom.jl
## =============================================================================

using CUDA
CUDA.device!(1)

using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.InelasticScattering
using NCDatasets
using DelimitedFiles
using Interpolations
using Printf
using Dates
using Statistics

# ---- Configuration ----------------------------------------------------------

const YAML_PATH     = joinpath(@__DIR__, "..", "test_parameters",
                              "ParamsEMIT_MODTRANcomp_newLUT.yaml")
const LAKE_MARY_DIR = expanduser("~/EMIT/Lake_Mary_Example")
const OUT_PATH      = joinpath(LAKE_MARY_DIR, "output", "lake_mary_TOA_vsmartmom.dat")
const LUT_PATH      = expanduser("~/data/EMIT_MODTRANcomp/emit_20250404.nc")

# Geometry (Lake Mary, 27 Mar 2025, fid emit20250327T212148)
const SZA    = 39.27
const SAA    = 212.97
const VZA    = 5.23
const VAA    = 48.49
# vSmartMOM internal convention: RAA = 180 - (VAA - SAA)
# 180 - (48.49 - 212.97) = 344.48° = -15.52° (wrapped to (-180, 180])
const RAA    = mod(180.0 - (VAA - SAA) + 360.0, 360.0)   # 344.48° (= -15.52°)
const μ₀     = cosd(SZA)

# Scene atmosphere
const GNDALT = 2.7    # km
const AOT550 = 0.04
const H2OCOL = 0.3    # g/cm² (column water vapor target)

# Ensure HITRAN LUT dir env var is set (YAML uses ${ENV:VSMARTMOM_HITRAN_LUT_DIR})
if !haskey(ENV, "VSMARTMOM_HITRAN_LUT_DIR")
    ENV["VSMARTMOM_HITRAN_LUT_DIR"] = expanduser("~/data/HITRAN_LUTs")
end
@info "VSMARTMOM_HITRAN_LUT_DIR = $(ENV["VSMARTMOM_HITRAN_LUT_DIR"])"

# ---- Read EMIT inputs -------------------------------------------------------

@info "Reading Lake Mary inputs"

emit_wl_data = readdlm(joinpath(LAKE_MARY_DIR, "emit_wavelengths.txt"), comments=true)
λ_emit_μm    = Float64.(emit_wl_data[:, 2])
fwhm_emit_μm = Float64.(emit_wl_data[:, 3])
λ_emit_nm    = λ_emit_μm .* 1000
fwhm_emit_nm = fwhm_emit_μm .* 1000
n_emit       = length(λ_emit_nm)

surf_data = readdlm(joinpath(LAKE_MARY_DIR, "surface_reflectance.txt"), ','; skipstart=1)
ρ_surf_mean = Float64.(surf_data[:, 2])

emit_rad_data = readdlm(joinpath(LAKE_MARY_DIR, "emit_toa_radiance.txt"), comments=true)
L_emit_meas = Float64.(emit_rad_data[:])

@info "  EMIT bands: $n_emit, λ range = [$(λ_emit_nm[1]), $(λ_emit_nm[end])] nm"

ρ_surf_itp = linear_interpolation(λ_emit_nm, ρ_surf_mean, extrapolation_bc=Flat())

# ---- Solar irradiance from MODTRAN LUT -------------------------------------

@info "Loading solar_irr from MODTRAN LUT"
ds_lut = NCDataset(LUT_PATH, "r")
wl_lut_nm   = Float64.(ds_lut["wl"][:])
solar_irr_lut = Float64.(ds_lut["solar_irr"][:])
close(ds_lut)
solar_itp = linear_interpolation(wl_lut_nm, solar_irr_lut, extrapolation_bc=Flat())

# ---- Load YAML and apply Lake Mary overrides --------------------------------

@info "Loading vSmartMOM parameters from $YAML_PATH"
params = parameters_from_yaml(YAML_PATH)
FT = params.float_type

# Geometry override
params.sza = FT(SZA)
params.vza = FT[VZA]
params.vaz = FT[RAA]

# Override aerosol τ_ref to Lake Mary's AOD-550
params.scattering_params.rt_aerosols[1].τ_ref = FT(AOT550)

# Override surface with LambertianSurfaceSpline (spectrum from Niklas's in-situ ρ_s).
spec_band_cm = params.spec_bands[1]
spec_band_nm = FT.(1e7 ./ spec_band_cm)
n_spec = length(spec_band_nm)
# Build interpolator over EMIT λ grid (typed as FT), and pass spec_band_nm as wlGrid.
ρ_surf_FT = FT.(ρ_surf_mean)
λ_emit_nm_FT = FT.(λ_emit_nm)
surf_interp = linear_interpolation(λ_emit_nm_FT, ρ_surf_FT, extrapolation_bc=Flat())
params.brdf = CoreRT.AbstractSurfaceType[CoreRT.LambertianSurfaceSpline{FT}(surf_interp, spec_band_nm)]
ρ_s_native = FT.(clamp.([ρ_surf_itp(λ) for λ in spec_band_nm], FT(0), FT(1)))   # for diagnostics

@info "  Geometry: SZA=$(params.sza), vza=$(params.vza), vaz=$(params.vaz)"
@info "  Aerosol τ_ref = $(params.scattering_params.rt_aerosols[1].τ_ref)"
@info "  vSmartMOM spec_bands: $n_spec points  range_nm = [$(minimum(spec_band_nm)), $(maximum(spec_band_nm))]"
@info "  Surface ρ extrema (native grid): $(round.(extrema(ρ_s_native); sigdigits=4))"

# H2O column scaling — uni_vSmartMOM newLUT YAML drives H2O via q-profile.
if hasproperty(params, :q) && params.q !== nothing && !isempty(params.q)
    @info "  Scaling q-profile by $H2OCOL for column H2O target = $H2OCOL g/cm²"
    params.q .= FT.(params.q .* H2OCOL)
else
    @warn "  No q profile found in params — H2O column unchanged"
end

@info "Building vSmartMOM model"
t0 = time()
model = model_from_parameters(params)
elapsed_build = (time() - t0) / 60.0
@info "  Model built in $(round(elapsed_build; digits=1)) min"

# ---- RT run -----------------------------------------------------------------

n_pol = CoreRT.polarization_type(model).n
RS_type = InelasticScattering.noRS(
    fscattRayl=[FT(1)], ϖ_Cabannes=[FT(1)],
    bandSpecLim=[], iBand=[1],
    F₀=zeros(FT, n_pol, n_spec),
    SIF₀=zeros(FT, n_pol, n_spec))
RS_type.F₀[1, :] .= FT(1)

@info "Running rt_run (multiple scattering)"
t0 = time()
result = CoreRT.rt_run(RS_type, model, 1)
elapsed_ms = (time() - t0) / 60.0
@info "  MS done in $(round(elapsed_ms; digits=1)) min, returned $(length(result))-tuple"

R_SFI_tot = result[1]
T_SFI_tot = result[2]
bhr_uw_vec = length(result) >= 6 ? Array(result[6]) : zeros(FT, n_spec)
bhr_dw_vec = length(result) >= 7 ? Array(result[7]) : zeros(FT, n_spec)

R_native = Array(R_SFI_tot)[1, 1, :]
T_native = Array(T_SFI_tot)[1, 1, :]

@info "  R_SFI extrema (native): $(round.(extrema(R_native); sigdigits=4))"
@info "  bhr_uw extrema: $(round.(extrema(bhr_uw_vec); sigdigits=4))"

# ---- Convolve to EMIT bands -------------------------------------------------

@info "Convolving native spectrum to EMIT bands"

perm = sortperm(spec_band_nm)
λ_native = spec_band_nm[perm]
R_sorted = R_native[perm]
T_sorted = T_native[perm]

function trapz_weights(λ)
    n = length(λ)
    w = zeros(eltype(λ), n)
    w[1] = (λ[2] - λ[1]) / 2
    w[end] = (λ[end] - λ[end-1]) / 2
    @inbounds for i in 2:n-1
        w[i] = (λ[i+1] - λ[i-1]) / 2
    end
    return w
end
dλ_native = trapz_weights(λ_native)

function convolve_to_emit(spec_native)
    out = zeros(Float64, n_emit)
    for b in 1:n_emit
        λ_c   = λ_emit_nm[b]
        σ     = fwhm_emit_nm[b] / 2.3548200450309493
        half_w = 3σ
        i_lo  = searchsortedfirst(λ_native, λ_c - half_w)
        i_hi  = searchsortedlast( λ_native, λ_c + half_w)
        i_lo  = max(i_lo, 1); i_hi = min(i_hi, length(λ_native))
        if i_hi < i_lo
            out[b] = 0.0
            continue
        end
        num = 0.0
        den = 0.0
        @inbounds for i in i_lo:i_hi
            k = exp(-((λ_native[i] - λ_c)^2) / (2σ^2))
            num += spec_native[i] * k * dλ_native[i]
            den += k * dλ_native[i]
        end
        out[b] = den > 0 ? num / den : 0.0
    end
    return out
end

R_emit_band = convolve_to_emit(Float64.(R_sorted))
T_emit_band = convolve_to_emit(Float64.(T_sorted))

ρ_s_emit       = ρ_surf_mean
solar_irr_emit = [solar_itp(λ) for λ in λ_emit_nm]
L_TOA_modeled  = R_emit_band .* solar_irr_emit
T_BOA_modeled  = T_emit_band .* solar_irr_emit

# ---- Save output ------------------------------------------------------------

@info "Writing output to $OUT_PATH"

open(OUT_PATH, "w") do io
    @printf(io, "# Lake Mary, CA — EMIT scene emit20250327T212148  direct vSmartMOM reconstruction\n")
    @printf(io, "# Generated: %s\n", string(now()))
    @printf(io, "# Geometry: SZA=%.3f  SAA=%.3f  VZA=%.3f  VAA=%.3f  RAA=%.3f (deg)\n", SZA, SAA, VZA, VAA, RAA)
    @printf(io, "# Scene:    GNDALT=%.2f km  AOT550=%.4f  H2O column scale=%.3f\n", GNDALT, AOT550, H2OCOL)
    @printf(io, "# Model:    YAML=%s\n", YAML_PATH)
    @printf(io, "#           FT=%s, n_native_spec=%d  Lambertian SurfaceSpectrum (Niklas in-situ)\n", FT, n_spec)
    @printf(io, "# Times:    model_build=%.1f min   MS_run=%.1f min\n", elapsed_build, elapsed_ms)
    @printf(io, "# Convolve: Gaussian, per-band FWHM from emit_wavelengths.txt\n")
    @printf(io, "# Radiance: L_TOA(μW/cm²/nm/sr) = R_SFI · solar_irr_LUT(λ)\n")
    @printf(io, "#\n")
    @printf(io, "# %-12s %-12s %-12s %-12s %-14s %-14s %-14s\n",
            "wl_nm", "fwhm_nm", "rho_surf", "R_SFI_MS", "L_TOA_MS", "T_BOA_MS", "L_TOA_meas")
    for i in 1:n_emit
        @printf(io, "  %12.6f %12.6f %12.6e %12.6e %14.6e %14.6e %14.6e\n",
                λ_emit_nm[i], fwhm_emit_nm[i],
                ρ_s_emit[i], R_emit_band[i],
                L_TOA_modeled[i], T_BOA_modeled[i],
                L_emit_meas[i])
    end
end

@info "Summary"
@info "  R_SFI (EMIT convolved) extrema:           $(round.(extrema(R_emit_band); sigdigits=4))"
@info "  L_TOA_modeled (MS) extrema (μW/cm²/nm/sr): $(round.(extrema(L_TOA_modeled); sigdigits=4))"
@info "  L_TOA_emit_meas extrema:                  $(round.(extrema(L_emit_meas); sigdigits=4))"
@info "Done — $OUT_PATH"
