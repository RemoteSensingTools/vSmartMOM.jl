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
"""

using Dates
using LinearAlgebra
using NCDatasets
using Printf

const OUTPUT = joinpath(@__DIR__, "apriori_states.nc")
const SUMMARY = joinpath(@__DIR__, "apriori_states.dat")
const N_PARAMETER = 34
const LAMBDA_SIF_NM = 760.0
const WAVENUMBER_CONVERSION = 1.0e7
const DEFAULT_SIF_REFERENCE_RADIANCE_MW_NM = 0.1

# Exact 16-layer centers from truth_map/aerosol_vertical_grids.dat, in the
# internal TOA-to-BOA order.
const Z_CENTER_KM = [
    40.792617798, 16.678894043, 13.178297043, 10.971210480,
     9.325729370,  7.991297245,  6.854076385,  5.850268841,
     4.944913864,  4.126883507,  3.370160103,  2.663361073,
     2.015735626,  1.404733658,  0.820331514,  0.268311203,
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

"One-sigma CO2 uncertainty in VMR at a layer-center altitude."
function co2_sigma_vmr(z_km)
    z_km < 1 && return 0.10 * 400e-6
    z_km < 3 && return 0.02 * 400e-6
    z_km < 10 && return 0.01 * 400e-6
    return 0.0
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

function build_prior(surface::Symbol, slope_reference_radiance_mw_nm::Real)
    surface in SURFACES || throw(ArgumentError("unknown surface: $surface"))
    isfinite(slope_reference_radiance_mw_nm) || throw(ArgumentError(
        "SIF slope reference radiance must be finite"))
    slope_reference_radiance_mw_nm > 0 || throw(ArgumentError(
        "SIF slope reference radiance must be positive and nonzero"))

    xa = zeros(Float64, N_PARAMETER)
    Sa = zeros(Float64, N_PARAMETER, N_PARAMETER)

    xa[1] = 1000.0
    Sa[1, 1] = 50.0^2

    xa[2:17] .= 400e-6
    for (index, z_km) in enumerate(Z_CENTER_KM)
        sigma = co2_sigma_vmr(z_km)
        Sa[index + 1, index + 1] = sigma^2
    end

    # Positivity-preserving retrieval coordinates u=ln(q). The first retrieval
    # suite deliberately uses sigma_u=2 rather than the earlier first-order
    # mapping of the broad physical-space +/-0.10 AOD statement, which gave 5
    # and allowed numerically extreme trial loadings.
    xa[18:20] .= log(0.02)
    for index in 18:20
        Sa[index, index] = 2.0^2
    end

    xa[21:23] .= log.(AEROSOL_Z0_KM)
    for index in 21:23
        Sa[index, index] = 0.10^2
    end

    for (iband, p0) in enumerate(SURFACE_P0[surface])
        first_index = 24 + 3 * (iband - 1)
        xa[first_index:first_index + 2] .= (p0, 0.0, 0.0)
        Sa[first_index, first_index] = (0.10 * p0)^2
        Sa[first_index + 1, first_index + 1] = 1e-3^2
        Sa[first_index + 2, first_index + 2] = 1e-4^2
    end

    # User-space wavelength coefficients and independent 1-sigma errors.
    a_lambda = Float64(slope_reference_radiance_mw_nm)
    sigma_a_lambda = 0.25
    fractional_slope = -0.035
    sigma_fractional_slope = 0.25 * abs(fractional_slope)
    b_lambda = fractional_slope * slope_reference_radiance_mw_nm
    sigma_b_lambda = sigma_fractional_slope * slope_reference_radiance_mw_nm

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
            slope_reference_radiance_mw_nm=Float64(slope_reference_radiance_mw_nm))
end

function write_netcdf(priors; output_path=OUTPUT)
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
        output.attrib["uncertainty_convention"] = "quoted +/- interpreted as one standard deviation"
        output.attrib["correlation_convention"] =
            "independent user-space parameters; exact 2x2 SIF coefficient transform retained"
        output.attrib["sif_slope_reference_radiance_mw_m2_sr_nm"] =
            first(values(priors)).slope_reference_radiance_mw_nm
        output.attrib["sif_fractional_slope_per_nm"] = -0.035
        output.attrib["sif_fractional_slope_sigma_per_nm"] = 0.00875
        output.attrib["sif_reference_wavelength_nm"] = LAMBDA_SIF_NM
        output.attrib["full_covariance_status"] =
            "positive semidefinite and singular because four CO2 entries are fixed"
        output.attrib["created_utc"] = string(now(UTC))
        output.attrib["apriori_complete"] = 1
    end
end

function write_summary(priors; output_path=SUMMARY)
    open(output_path, "w") do io
        println(io, "# Retrieval-coordinate 34-element a priori states and marginal sigmas.")
        println(io, "# Aerosol AOD and z0 entries are natural logarithms of physical values.")
        println(io, "# Covariance entry Sa[33,34] is nonzero; use the NetCDF for the full matrix.")
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
    priors = Dict(surface => build_prior(surface, slope_reference) for surface in SURFACES)
    write_netcdf(priors)
    write_summary(priors)
    println(OUTPUT)
    println(SUMMARY)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
