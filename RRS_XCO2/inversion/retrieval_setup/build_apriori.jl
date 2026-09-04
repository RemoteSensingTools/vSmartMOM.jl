#!/usr/bin/env julia

"""
Build the four surface-specific RRS-XCO2 a priori states and covariances.

Usage:

    julia --project=. RRS_XCO2/inversion/retrieval_setup/build_apriori.jl \
        [SIF_REFERENCE_RADIANCE_MW_NM]

The reference radiance defaults to the retrieval-suite value
`0.1 mW m^-2 sr^-1 nm^-1`. It is both the small nonzero SIF760 a priori and
the scale that turns the fractional slope `g=-0.035 nm^-1` into the absolute
slope `dL_lambda/dlambda=g*L_scale`. Supplying an argument remains useful for
controlled sensitivity tests.

Environment variables:

- `APRIORI_OUTPUT`: NetCDF output path (default `apriori_states.nc` beside
  this script).
- `APRIORI_SUMMARY`: ASCII summary path (default `apriori_states.dat` beside
  this script).
- `AEROSOL_LN_AOD_SIGMA`: common one-sigma prior width for all three
  `ln(AOD760)` retrieval coordinates (default `0.75`). Override it only for
  controlled sensitivity studies; use a distinct `APRIORI_OUTPUT` for such
  tests so the authoritative prior is not replaced.
- `SURFACE_P1_SIGMA` and `SURFACE_P2_SIGMA`: common one-sigma surface slope
  and curvature priors for all three bands (defaults `1e-3` and `1e-4`). The
  bottom-layer CO2 campaign uses `2e-3` for both.
- `O2A_SURFACE_P{1,2}_SIGMA`, `WEAK_CO2_SURFACE_P{1,2}_SIGMA`, and
  `STRONG_CO2_SURFACE_P{1,2}_SIGMA`: optional per-band overrides of those
  common values for controlled sensitivity studies.
- `SIF_SLOPE_MEAN_MW_NM2`: wavelength-space absolute SIF-slope prior mean.
  Its historical default is `-0.0035`; the no-SIF-first bottom-layer campaign
  uses zero.
- `SIF_SLOPE_SIGMA_MW_NM2`: absolute wavelength-space SIF-slope 1-sigma
  width. Its historical default is `0.00875*0.1 = 8.75e-4`; the bottom-layer
  campaign uses three times that value, `2.625e-3`.
"""

using Dates
using LinearAlgebra
using NCDatasets
using Printf

const OUTPUT = get(
    ENV, "APRIORI_OUTPUT", joinpath(@__DIR__, "apriori_states.nc"))
const SUMMARY = get(
    ENV, "APRIORI_SUMMARY", joinpath(@__DIR__, "apriori_states.dat"))
const N_PARAMETER = 34
const LAMBDA_SIF_NM = 760.0
const WAVENUMBER_CONVERSION = 1.0e7
const DEFAULT_SIF_REFERENCE_RADIANCE_MW_NM = 0.1
const DEFAULT_SIF_FRACTIONAL_SLOPE_PER_NM = -0.035
const DEFAULT_SIF_SLOPE_MEAN_MW_NM2 =
    DEFAULT_SIF_FRACTIONAL_SLOPE_PER_NM *
    DEFAULT_SIF_REFERENCE_RADIANCE_MW_NM
const DEFAULT_SIF_SLOPE_SIGMA_MW_NM2 =
    0.25 * abs(DEFAULT_SIF_FRACTIONAL_SLOPE_PER_NM) *
    DEFAULT_SIF_REFERENCE_RADIANCE_MW_NM
const CO2_COVARIANCE_FILE = joinpath(@__DIR__, "co2_prior_covariances.dat")
const CO2_COVARIANCE_MODEL = "acos_mapped"
const DEFAULT_AEROSOL_LN_AOD_SIGMA = 0.75
const DEFAULT_SURFACE_P1_SIGMA = 1e-3
const DEFAULT_SURFACE_P2_SIGMA = 1e-4

# Exact 16-layer centers from truth_map/aerosol_vertical_grids.dat, in the
# internal TOA-to-BOA order.
const Z_CENTER_KM = [
    40.792617798, 16.678894043, 13.178297043, 10.971210480,
     9.325729370,  7.991297245,  6.854076385,  5.850268841,
     4.944913864,  4.126883507,  3.370160103,  2.663361073,
     2.015735626,  1.404733658,  0.820331514,  0.268311203,
]

# Exact layer-center pressures after reducing the common profile to 16 layers
# at psurf=1000 hPa. These are the target nodes used to map the 20-level ACOS
# covariance in normalized pressure.
const P_CENTER_HPA = [
     31.385625,  93.876875, 156.368125, 218.859375,
    281.350625, 343.841875, 406.333125, 468.824375,
    531.315625, 593.806875, 656.298125, 718.789375,
    781.280625, 843.771875, 906.263125, 968.754375,
]

# Dry-column weights from that same reduced 1000 hPa profile. They are used
# only to validate and document the implied CO2-only prior uncertainty in
# XCO2; the retrieval recomputes its pressure-dependent diagnostic weights.
const DRY_COLUMN_WEIGHTS = [
    0.06271458366696699, 0.06267999123579230,
    0.06264827279449964, 0.06261831034665068,
    0.06258999065374835, 0.06255915670625688,
    0.06253127935109044, 0.06250621030757231,
    0.06248208275638748, 0.06245969764265345,
    0.06243817296080489, 0.06241161154141432,
    0.06238468284683572, 0.06235227383324737,
    0.06232590314728601, 0.06229778020879323,
]

const AEROSOL_Z0_KM = [1.525651538, 2.112319568, 12.120602005]
const SURFACES = (:urban, :rural, :desert, :forest)
const BANDS = (:o2a, :weak_co2, :strong_co2)
const SURFACE_P0 = Dict(
    :urban  => [0.2715374708583, 0.2522486169619, 0.2176911056292],
    :rural  => [0.4316337681355, 0.2649756999501, 0.1335417267332],
    :desert => [0.4186071552534, 0.4972482231007, 0.4821355300133],
    :forest => [0.4682630263673, 0.2557722421504, 0.1100470188428],
)

const PARAMETER_NAMES = vcat(
    ["psurf"],
    [@sprintf("co2_vmr_layer%02d", layer) for layer in 1:16],
    ["ln_sulfate_aod760", "ln_organic_carbon_aod760", "ln_utls_sulfate_aod760"],
    ["ln_sulfate_z0", "ln_organic_carbon_z0", "ln_utls_sulfate_z0"],
    ["o2a_surface_P0", "o2a_surface_P1", "o2a_surface_P2"],
    ["weak_co2_surface_P0", "weak_co2_surface_P1", "weak_co2_surface_P2"],
    ["strong_co2_surface_P0", "strong_co2_surface_P1", "strong_co2_surface_P2"],
    ["SIF760", "mSIF"],
)

const PARAMETER_UNITS = vcat(
    ["hPa"],
    fill("mol mol-1", 16),
    fill("1", 6),
    fill("1", 9),
    ["mW m-2 sr-1 (cm-1)-1", "mW m-2 sr-1 (cm-1)-2"],
)

length(PARAMETER_NAMES) == N_PARAMETER || error("parameter-name layout is inconsistent")
length(PARAMETER_UNITS) == N_PARAMETER || error("parameter-unit layout is inconsistent")

"""Read one symmetric 16-layer CO2 covariance from the comparison table."""
function read_co2_covariance(model::AbstractString;
                             path::AbstractString=CO2_COVARIANCE_FILE)
    isfile(path) || throw(ArgumentError("missing CO2 covariance table: $path"))
    covariance_ppm2 = zeros(Float64, 16, 16)
    seen = Set{Tuple{Int,Int}}()
    for line in eachline(path)
        stripped = strip(line)
        (isempty(stripped) || startswith(stripped, '#')) && continue
        fields = split(stripped)
        length(fields) == 4 || error(
            "invalid CO2 covariance row in $path: $stripped")
        fields[1] == model || continue
        row = parse(Int, fields[2])
        column = parse(Int, fields[3])
        1 <= column <= row <= 16 || error(
            "CO2 covariance must tabulate the lower 16-layer triangle")
        entry = (row, column)
        entry in seen && error("duplicate $model CO2 covariance entry $entry")
        push!(seen, entry)
        value = parse(Float64, fields[4])
        isfinite(value) || error("non-finite $model CO2 covariance entry $entry")
        covariance_ppm2[row, column] = value
        covariance_ppm2[column, row] = value
    end
    isempty(seen) && error("CO2 covariance model '$model' is absent from $path")
    all(iszero, covariance_ppm2[1:4, :]) || error(
        "$model covariance must keep CO2 layers 1:4 fixed")
    active = Symmetric(covariance_ppm2[5:16, 5:16])
    isposdef(active) || error("$model active CO2 covariance is not positive definite")
    return covariance_ppm2 .* 1e-12 # ppm^2 -> (mol mol^-1)^2
end

function co2_xco2_sigma_ppm(covariance_vmr2::AbstractMatrix)
    size(covariance_vmr2) == (16, 16) || throw(DimensionMismatch(
        "CO2 covariance must be 16 by 16"))
    variance = dot(DRY_COLUMN_WEIGHTS,
                   covariance_vmr2 * DRY_COLUMN_WEIGHTS)
    variance >= 0 || error("CO2 covariance gives a negative XCO2 variance")
    return 1e6 * sqrt(variance)
end

"""
Transform local wavelength-density SIF coefficients into model-native
wavenumber-density coefficients at 760 nm.

For `u=[a,b]=[L_lambda, dL_lambda/dlambda]` and
`v=[SIF760,mSIF]=[L_nu,dL_nu/dnu]`, `v=T*u`. This transformation includes the
spectral-density Jacobian, not only the independent-variable chain rule.
"""
function sif_transform()
    lambda = LAMBDA_SIF_NM
    C = WAVENUMBER_CONVERSION
    return [lambda^2 / C  0.0;
            -2lambda^3 / C^2  -lambda^4 / C^2]
end

function build_prior(surface::Symbol, slope_reference_radiance_mw_nm::Real;
                     sif_wavelength_slope_mean::Real=
                         DEFAULT_SIF_FRACTIONAL_SLOPE_PER_NM *
                         slope_reference_radiance_mw_nm,
                     sif_wavelength_slope_sigma::Real=
                         0.25 * abs(DEFAULT_SIF_FRACTIONAL_SLOPE_PER_NM) *
                         slope_reference_radiance_mw_nm,
                     aerosol_ln_aod_sigma::Real=
                         DEFAULT_AEROSOL_LN_AOD_SIGMA,
                     surface_p1_sigmas::NTuple{3,<:Real}=
                         ntuple(_ -> DEFAULT_SURFACE_P1_SIGMA, 3),
                     surface_p2_sigmas::NTuple{3,<:Real}=
                         ntuple(_ -> DEFAULT_SURFACE_P2_SIGMA, 3))
    surface in SURFACES || throw(ArgumentError("unknown surface: $surface"))
    isfinite(slope_reference_radiance_mw_nm) || throw(ArgumentError(
        "SIF slope reference radiance must be finite"))
    slope_reference_radiance_mw_nm > 0 || throw(ArgumentError(
        "SIF slope reference radiance must be positive and nonzero"))
    isfinite(sif_wavelength_slope_mean) || throw(ArgumentError(
        "SIF wavelength-space slope mean must be finite"))
    isfinite(sif_wavelength_slope_sigma) && sif_wavelength_slope_sigma > 0 ||
        throw(ArgumentError(
            "SIF wavelength-space slope sigma must be finite and positive"))
    isfinite(aerosol_ln_aod_sigma) && aerosol_ln_aod_sigma > 0 ||
        throw(ArgumentError(
            "aerosol ln(AOD) sigma must be finite and positive"))
    all(sigma -> isfinite(sigma) && sigma > 0, surface_p1_sigmas) ||
        throw(ArgumentError("surface P1 sigmas must be finite and positive"))
    all(sigma -> isfinite(sigma) && sigma > 0, surface_p2_sigmas) ||
        throw(ArgumentError("surface P2 sigmas must be finite and positive"))

    xa = zeros(Float64, N_PARAMETER)
    Sa = zeros(Float64, N_PARAMETER, N_PARAMETER)

    xa[1] = 1000.0
    Sa[1, 1] = 50.0^2

    xa[2:17] .= 400e-6
    co2_covariance = read_co2_covariance(CO2_COVARIANCE_MODEL)
    Sa[2:17, 2:17] .= co2_covariance
    isapprox(co2_xco2_sigma_ppm(co2_covariance), 13.71563311;
             atol=1e-7, rtol=0) || error(
        "mapped ACOS covariance failed its XCO2-uncertainty regression")

    # Positivity-preserving retrieval coordinates u=ln(q). Earlier widths of
    # sigma_u=5 and then 2 allowed numerically extreme aerosol trial loadings.
    # A noiseless bottom-layer sensitivity comparison on clear state 001 and
    # aerosol state 013 supported the tighter sigma_u=0.75 production prior.
    xa[18:20] .= log(0.02)
    for index in 18:20
        Sa[index, index] = Float64(aerosol_ln_aod_sigma)^2
    end

    xa[21:23] .= log.(AEROSOL_Z0_KM)
    for index in 21:23
        Sa[index, index] = 0.10^2
    end

    for (iband, p0) in enumerate(SURFACE_P0[surface])
        first_index = 24 + 3 * (iband - 1)
        xa[first_index:first_index + 2] .= (p0, 0.0, 0.0)
        Sa[first_index, first_index] = (0.10 * p0)^2
        Sa[first_index + 1, first_index + 1] =
            Float64(surface_p1_sigmas[iband])^2
        p2_sigma = Float64(surface_p2_sigmas[iband])
        Sa[first_index + 2, first_index + 2] = p2_sigma^2
    end

    # User-space wavelength coefficients and independent 1-sigma errors.
    a_lambda = Float64(slope_reference_radiance_mw_nm)
    sigma_a_lambda = 0.25
    b_lambda = Float64(sif_wavelength_slope_mean)
    sigma_b_lambda = Float64(sif_wavelength_slope_sigma)

    T = sif_transform()
    sif_native = T * [a_lambda, b_lambda]
    sif_cov_native = T * Diagonal([sigma_a_lambda^2, sigma_b_lambda^2]) * T'
    xa[33:34] .= sif_native
    Sa[33:34, 33:34] .= sif_cov_native

    active = diag(Sa) .> 0
    count(active) == 30 || error("expected 30 active parameters, got $(count(active))")
    isposdef(Symmetric(Sa[active, active])) || error("active covariance is not positive definite")
    return (; xa, Sa, active, active_indices=findall(active),
            sif_wavelength_state=[a_lambda, b_lambda],
            sif_wavelength_sigma=[sigma_a_lambda, sigma_b_lambda],
            slope_reference_radiance_mw_nm=Float64(slope_reference_radiance_mw_nm),
            sif_wavelength_slope_mean=Float64(sif_wavelength_slope_mean),
            sif_wavelength_slope_sigma=Float64(sif_wavelength_slope_sigma),
            aerosol_ln_aod_sigma=Float64(aerosol_ln_aod_sigma),
            surface_p1_sigmas=Tuple(Float64.(surface_p1_sigmas)),
            surface_p2_sigmas=Tuple(Float64.(surface_p2_sigmas)))
end

function write_netcdf(priors; output_path=OUTPUT)
    mkpath(dirname(output_path))
    isfile(output_path) && rm(output_path)
    n_active = count(first(values(priors)).active)
    NCDataset(output_path, "c") do output
        defDim(output, "parameter", N_PARAMETER)
        defDim(output, "parameter_2", N_PARAMETER)
        defDim(output, "active_parameter", n_active)
        defDim(output, "active_parameter_2", n_active)
        defDim(output, "surface", length(SURFACES))
        defDim(output, "sif_wavelength_parameter", 2)

        xa = defVar(output, "xa", Float64, ("parameter", "surface"))
        xa.attrib["long_name"] = "surface-specific retrieval-coordinate a priori state"
        covariance = defVar(
            output, "Sa", Float64, ("parameter", "parameter_2", "surface"))
        covariance.attrib["long_name"] = "full retrieval-coordinate a priori covariance"
        sigma = defVar(output, "prior_sigma", Float64, ("parameter", "surface"))
        sigma.attrib["long_name"] = "square root of the full covariance diagonal"
        active_covariance = defVar(
            output, "Sa_active", Float64,
            ("active_parameter", "active_parameter_2", "surface"))
        active_covariance.attrib["long_name"] =
            "positive-definite covariance after removing fixed CO2 entries"
        active_index = defVar(output, "active_parameter_index", Int16, ("active_parameter",))
        active_index.attrib["index_convention"] = "one-based index into full 34-element state"
        active_index[:] = Int16.(first(values(priors)).active_indices)
        fixed_mask = defVar(output, "active_mask", Int8, ("parameter",))
        fixed_mask.attrib["flag_values"] = Int8[0, 1]
        fixed_mask.attrib["flag_meanings"] = "fixed active"
        fixed_mask[:] = Int8.(first(values(priors)).active)
        centers = defVar(output, "co2_layer_center_height", Float64, ("parameter",))
        centers.attrib["units"] = "km"
        centers.attrib["comment"] = "NaN for parameters other than the 16 CO2 layers"
        center_values = fill(NaN, N_PARAMETER)
        center_values[2:17] .= Z_CENTER_KM
        centers[:] = center_values
        pressures = defVar(output, "co2_layer_center_pressure", Float64,
                           ("parameter",))
        pressures.attrib["units"] = "hPa"
        pressures.attrib["comment"] =
            "NaN for parameters other than the 16 CO2 layers; psurf=1000 hPa grid"
        pressure_values = fill(NaN, N_PARAMETER)
        pressure_values[2:17] .= P_CENTER_HPA
        pressures[:] = pressure_values
        sif_wavelength = defVar(
            output, "sif_wavelength_state", Float64,
            ("sif_wavelength_parameter", "surface"))
        sif_wavelength.attrib["units"] =
            "mixed: mW m-2 sr-1 nm-1; mW m-2 sr-1 nm-2"
        sif_wavelength.attrib["order"] = "Llambda_760 dLlambda_dlambda_760"
        sif_wavelength_sigma = defVar(
            output, "sif_wavelength_sigma", Float64,
            ("sif_wavelength_parameter", "surface"))
        sif_wavelength_sigma.attrib["order"] = "sigma_Llambda_760 sigma_dLlambda_dlambda_760"

        for (isurface, surface) in enumerate(SURFACES)
            prior = priors[surface]
            xa[:, isurface] = prior.xa
            covariance[:, :, isurface] = prior.Sa
            sigma[:, isurface] = sqrt.(diag(prior.Sa))
            active_covariance[:, :, isurface] =
                prior.Sa[prior.active, prior.active]
            sif_wavelength[:, isurface] = prior.sif_wavelength_state
            sif_wavelength_sigma[:, isurface] = prior.sif_wavelength_sigma
        end

        output.attrib["surface_order"] = join(string.(SURFACES), " ")
        output.attrib["parameter_names"] = join(PARAMETER_NAMES, " ")
        output.attrib["parameter_units"] = join(PARAMETER_UNITS, " | ")
        output.attrib["state_order"] =
            "psurf; 16 CO2 TOA-to-BOA; 3 ln(AOD760); 3 ln(aerosol z0/km); 9 surface; SIF760; mSIF"
        output.attrib["aerosol_coordinate_transform"] =
            "u=ln(q); q=exp(u); K_u=K_q*q at the current state; Sa_u=diag(1/q_a)*Sa_q*diag(1/q_a)"
        output.attrib["aerosol_ln_aod_sigma"] =
            first(values(priors)).aerosol_ln_aod_sigma
        output.attrib["uncertainty_convention"] = "quoted +/- interpreted as one standard deviation"
        output.attrib["correlation_convention"] =
            "mapped ACOS covariance among active CO2 layers; exact 2x2 SIF coefficient transform; other blocks independent"
        output.attrib["co2_covariance_model"] = CO2_COVARIANCE_MODEL
        output.attrib["co2_covariance_source_file"] = abspath(CO2_COVARIANCE_FILE)
        output.attrib["co2_covariance_source_dataset"] =
            "RetrievalResults/apriori_covariance_matrix[1:20,1:20] from the four oco_gain.jl orbit files"
        output.attrib["co2_covariance_mapping"] =
            "linear interpolation H in normalized pressure; S16=H*S20*transpose(H); active marginal block layers 5:16"
        output.attrib["co2_prior_xco2_sigma_ppm_at_1000hpa"] =
            co2_xco2_sigma_ppm(first(values(priors)).Sa[2:17, 2:17])
        output.attrib["sif_slope_reference_radiance_mw_m2_sr_nm"] =
            first(values(priors)).slope_reference_radiance_mw_nm
        output.attrib["sif_wavelength_slope_prior_mw_m2_sr_nm2"] =
            first(values(priors)).sif_wavelength_slope_mean
        output.attrib["sif_fractional_slope_per_nm"] =
            DEFAULT_SIF_FRACTIONAL_SLOPE_PER_NM
        output.attrib["sif_wavelength_slope_sigma_mw_m2_sr_nm2"] =
            first(values(priors)).sif_wavelength_slope_sigma
        output.attrib["sif_fractional_slope_sigma_per_nm"] =
            first(values(priors)).sif_wavelength_slope_sigma /
            first(values(priors)).slope_reference_radiance_mw_nm
        output.attrib["sif_slope_prior_interpretation"] =
            "absolute wavelength-space slope mean is independent of the 0.1 radiance scale used to set its sigma"
        output.attrib["sif_reference_wavelength_nm"] = LAMBDA_SIF_NM
        output.attrib["surface_p1_sigma_o2a"] =
            first(values(priors)).surface_p1_sigmas[1]
        output.attrib["surface_p1_sigma_weak_co2"] =
            first(values(priors)).surface_p1_sigmas[2]
        output.attrib["surface_p1_sigma_strong_co2"] =
            first(values(priors)).surface_p1_sigmas[3]
        output.attrib["surface_p2_sigma_o2a"] =
            first(values(priors)).surface_p2_sigmas[1]
        output.attrib["surface_p2_sigma_weak_co2"] =
            first(values(priors)).surface_p2_sigmas[2]
        output.attrib["surface_p2_sigma_strong_co2"] =
            first(values(priors)).surface_p2_sigmas[3]
        output.attrib["full_covariance_status"] =
            "positive semidefinite and singular because four CO2 entries are fixed"
        output.attrib["created_utc"] = string(now(UTC))
        output.attrib["apriori_complete"] = 1
    end
end

function write_summary(priors; output_path=SUMMARY)
    mkpath(dirname(output_path))
    open(output_path, "w") do io
        println(io, "# Retrieval-coordinate 34-element a priori states and marginal sigmas.")
        println(io, "# Aerosol AOD and z0 entries are natural logarithms of physical values.")
        println(io, "# CO2 layers 5:16 use the mapped ACOS off-diagonal covariance.")
        println(io, "# The SIF block is also correlated; use the NetCDF or covariance table for full matrices.")
        @printf(io, "# CO2-only sigma(XCO2) at 1000 hPa = %.9f ppm\n",
                co2_xco2_sigma_ppm(first(values(priors)).Sa[2:17, 2:17]))
        @printf(io, "# Aerosol ln(AOD760) sigma: %.12e\n",
                first(values(priors)).aerosol_ln_aod_sigma)
        @printf(io, "# Surface P1 sigmas: O2A=%.12e weak_CO2=%.12e strong_CO2=%.12e\n",
                first(values(priors)).surface_p1_sigmas...)
        @printf(io, "# Surface P2 sigmas: O2A=%.12e weak_CO2=%.12e strong_CO2=%.12e\n",
                first(values(priors)).surface_p2_sigmas...)
        @printf(io, "# SIF wavelength-space absolute slope: mean=% .12e sigma=%.12e mW m-2 sr-1 nm-2\n",
                first(values(priors)).sif_wavelength_slope_mean,
                first(values(priors)).sif_wavelength_sigma[2])
        println(io, "# surface index parameter xa sigma units active")
        for surface in SURFACES
            prior = priors[surface]
            marginal_sigma = sqrt.(diag(prior.Sa))
            for index in 1:N_PARAMETER
                @printf(io, "%-7s %2d %-29s % .12e %.12e %-29s %d\n",
                        String(surface), index, PARAMETER_NAMES[index], prior.xa[index],
                        marginal_sigma[index], PARAMETER_UNITS[index], prior.active[index])
            end
        end
        println(io, "# SIF covariance blocks (surface values are identical)")
        for surface in SURFACES
            block = priors[surface].Sa[33:34, 33:34]
            @printf(io, "# %-7s Sa33_33=%.12e Sa33_34=% .12e Sa34_34=%.12e\n",
                    String(surface), block[1, 1], block[1, 2], block[2, 2])
        end
    end
end

function main(args=ARGS)
    length(args) <= 1 || error(
        "usage: build_apriori.jl [SIF_REFERENCE_RADIANCE_MW_NM]")
    slope_reference = isempty(args) ?
        DEFAULT_SIF_REFERENCE_RADIANCE_MW_NM : parse(Float64, only(args))
    sif_wavelength_slope_mean = parse(
        Float64, get(ENV, "SIF_SLOPE_MEAN_MW_NM2",
                     string(DEFAULT_SIF_FRACTIONAL_SLOPE_PER_NM *
                            slope_reference)))
    sif_wavelength_slope_sigma = parse(
        Float64, get(ENV, "SIF_SLOPE_SIGMA_MW_NM2",
                     string(0.25 * abs(DEFAULT_SIF_FRACTIONAL_SLOPE_PER_NM) *
                            slope_reference)))
    aerosol_ln_aod_sigma = parse(
        Float64, get(ENV, "AEROSOL_LN_AOD_SIGMA",
                     string(DEFAULT_AEROSOL_LN_AOD_SIGMA)))
    common_surface_p1_sigma = parse(
        Float64, get(ENV, "SURFACE_P1_SIGMA",
                     string(DEFAULT_SURFACE_P1_SIGMA)))
    surface_p1_sigmas = (
        parse(Float64, get(ENV, "O2A_SURFACE_P1_SIGMA",
                           string(common_surface_p1_sigma))),
        parse(Float64, get(ENV, "WEAK_CO2_SURFACE_P1_SIGMA",
                           string(common_surface_p1_sigma))),
        parse(Float64, get(ENV, "STRONG_CO2_SURFACE_P1_SIGMA",
                           string(common_surface_p1_sigma))),
    )
    common_surface_p2_sigma = parse(
        Float64, get(ENV, "SURFACE_P2_SIGMA",
                     string(DEFAULT_SURFACE_P2_SIGMA)))
    surface_p2_sigmas = (
        parse(Float64, get(ENV, "O2A_SURFACE_P2_SIGMA",
                           string(common_surface_p2_sigma))),
        parse(Float64, get(ENV, "WEAK_CO2_SURFACE_P2_SIGMA",
                           string(common_surface_p2_sigma))),
        parse(Float64, get(ENV, "STRONG_CO2_SURFACE_P2_SIGMA",
                           string(common_surface_p2_sigma))),
    )
    priors = Dict(
        surface => build_prior(
            surface, slope_reference;
            sif_wavelength_slope_mean=sif_wavelength_slope_mean,
            sif_wavelength_slope_sigma=sif_wavelength_slope_sigma,
            aerosol_ln_aod_sigma=aerosol_ln_aod_sigma,
            surface_p1_sigmas=surface_p1_sigmas,
            surface_p2_sigmas=surface_p2_sigmas)
        for surface in SURFACES)
    write_netcdf(priors)
    write_summary(priors)
    println(OUTPUT)
    println(SUMMARY)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
