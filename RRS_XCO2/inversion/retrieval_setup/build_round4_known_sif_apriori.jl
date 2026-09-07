#!/usr/bin/env julia

module Round4KnownSIFPriorBuilder

"""
Production prior builder for the round-4 known-SIF retrieval campaign.

This is the sole supported production path for round-4 prior files. It keeps
the approved round-3 non-SIF blocks exactly, fixes the SIF radiance at 759 nm,
and retains the existing physical wavelength-slope constraint at 760 nm.
Gaussian conditioning of the old uncertain-amplitude two-coordinate SIF block
is scientifically different: it shifts and narrows the physical slope prior.
It must not be substituted for this construction in round 4.
"""

using Dates
using LinearAlgebra
using NCDatasets
using Printf
using SHA

include(joinpath(@__DIR__, "..", "Round4SIFTruthConvention.jl"))
using .Round4SIFTruthConvention

export BASE_ACTIVE_TO_FULL,
       ROUND4_MODEL,
       ROUND4_OFF_ACTIVE_TO_FULL,
       ROUND4_ON_ACTIVE_TO_FULL,
       default_output_directory,
       file_sha256,
       native_slope_prior,
       output_filename,
       build_round4_prior,
       write_round4_prior,
       write_round4_summary,
       main

const FULL_STATE_COUNT = 34
const CORE_REFERENCE_WAVELENGTH_NM = 760.0
const KNOWN_SIF_WAVELENGTH_NM = 759.0
const WAVENUMBER_CONVERSION_NM_CM1 = 1.0e7
const CORE_SIF760_FULL_INDEX = 33
const CORE_MSIF_FULL_INDEX = 34
const BASE_ACTIVE_TO_FULL = vcat(1, collect(6:34))
const ROUND4_OFF_ACTIVE_TO_FULL = vcat(1, collect(6:32))
const ROUND4_ON_ACTIVE_TO_FULL = vcat(ROUND4_OFF_ACTIVE_TO_FULL, 34)
const ROUND4_OFF_ACTIVE_TO_CORE = collect(1:28)
const ROUND4_ON_ACTIVE_TO_CORE = vcat(collect(1:28), 30)
const ROUND4_MODEL =
    "round4_known_sif759_acos_mapped_tapered_vertical_correlation"
const REQUIRED_CO2_MODEL = "acos_mapped_tapered_vertical_correlation"
const DEFAULT_SOURCE_PRIOR = normpath(joinpath(
    @__DIR__, "..", "..", "bottom_layer_XCO2_retrievals",
    "retrieval_setup",
    "apriori_states_acos_mapped_tapered_vertical_correlation.nc"))
const DEFAULT_OUTPUT_DIRECTORY = normpath(joinpath(
    @__DIR__, "..", "..", "bottom_layer_XCO2_retrievals",
    "round4_known_sif759", "retrieval_setup"))

default_output_directory() = DEFAULT_OUTPUT_DIRECTORY

function _check_mode(mode::Symbol)
    mode in (:on, :off) || throw(ArgumentError(
        "round-4 SIF mode must be :on or :off, got :$mode"))
    return mode
end

function output_filename(mode::Symbol; extension::AbstractString="nc")
    _check_mode(mode)
    extension in ("nc", "dat") || throw(ArgumentError(
        "round-4 prior extension must be nc or dat"))
    return "apriori_states_round4_known_sif759_$(mode)_" *
           "acos_mapped_tapered_vertical_correlation.$extension"
end

function file_sha256(path::AbstractString)
    isfile(path) || throw(ArgumentError("missing file for SHA-256: $path"))
    return open(path, "r") do io
        bytes2hex(sha256(io))
    end
end

function _attributes(attributes)
    return Dict{String,Any}(String(key) => value
                            for (key, value) in pairs(attributes))
end

function _variable_attributes(variable)
    return _attributes(variable.attrib)
end

function _read_source_prior(path::AbstractString)
    isfile(path) || throw(ArgumentError(
        "missing round-3 tapered source prior: $path"))
    source_sha256 = file_sha256(path)
    return NCDataset(path, "r") do dataset
        get(dataset.attrib, "apriori_complete", 0) == 1 || error(
            "source prior is not marked complete: $path")
        get(dataset.attrib, "co2_covariance_model", "") ==
            REQUIRED_CO2_MODEL || error(
            "round 4 requires the tapered mapped-ACOS source prior")
        active = Int.(dataset["active_parameter_index"][:])
        active == BASE_ACTIVE_TO_FULL || error(
            "source prior does not have the approved 30-parameter mapping")
        size(dataset["xa"]) == (FULL_STATE_COUNT, 4) || error(
            "source prior must contain four 34-element means")
        size(dataset["Sa"]) == (FULL_STATE_COUNT, FULL_STATE_COUNT, 4) ||
            error("source prior must contain four 34x34 covariances")

        slope_mean = Float64(get(
            dataset.attrib, "sif_wavelength_slope_prior_mw_m2_sr_nm2", NaN))
        slope_sigma = Float64(get(
            dataset.attrib, "sif_wavelength_slope_sigma_mw_m2_sr_nm2", NaN))
        iszero(slope_mean) || error(
            "round 4 expects the current zero-centered wavelength-slope prior")
        isapprox(slope_sigma, 0.002625; atol=0, rtol=8eps(Float64)) || error(
            "round 4 expects the current wavelength-slope sigma of 0.002625")

        source_xa = Float64.(dataset["xa"][:, :])
        source_Sa = Float64.(dataset["Sa"][:, :, :])
        all(iszero, source_Sa[33:34, 1:32, :]) || error(
            "source prior couples SIF to non-SIF coordinates; refusing to discard that covariance")
        all(iszero, source_Sa[1:32, 33:34, :]) || error(
            "source prior couples non-SIF to SIF coordinates; refusing to discard that covariance")

        variable_names = (
            "co2_layer_center_height",
            "co2_layer_center_pressure",
            "co2_correlation_taper_adjacent_retention",
            "co2_base_adjacent_correlation",
            "co2_selected_adjacent_correlation",
        )
        auxiliary = Dict{String,Any}()
        auxiliary_attributes = Dict{String,Dict{String,Any}}()
        for name in variable_names
            haskey(dataset, name) || error(
                "source prior is missing required variable $name")
            auxiliary[name] = Array(dataset[name][:])
            auxiliary_attributes[name] = _variable_attributes(dataset[name])
        end

        return (;
            path=abspath(path),
            sha256=source_sha256,
            attributes=_attributes(dataset.attrib),
            xa=source_xa,
            Sa=source_Sa,
            surface_order=split(String(dataset.attrib["surface_order"])),
            parameter_names=split(String(dataset.attrib["parameter_names"])),
            parameter_units=split(String(dataset.attrib["parameter_units"]), " | "),
            slope_mean,
            slope_sigma,
            auxiliary,
            auxiliary_attributes,
        )
    end
end

"""
    native_slope_prior(Lnu759, wavelength_slope_mean, wavelength_slope_sigma)

Convert the retained round-3 physical slope constraint at 760 nm into the
native free coordinate `mSIF=dLnu/dnu`, while enforcing exact `Lnu(759)`.

The retrieval source is linear in wavenumber,

`Lnu(nu) = Lnu759 + mSIF*(nu-nu759)`.

The spectral-density Jacobian gives an affine relationship
`dLlambda/dlambda|760 = alpha*mSIF + beta`.  This function transforms the
mean and one-sigma width of that wavelength-space slope without reintroducing
an uncertain SIF amplitude.
"""
function native_slope_prior(Lnu759::Real,
                            wavelength_slope_mean::Real,
                            wavelength_slope_sigma::Real)
    all(isfinite, (Lnu759, wavelength_slope_mean,
                   wavelength_slope_sigma)) || throw(ArgumentError(
        "round-4 SIF prior inputs must be finite"))
    Lnu759 >= 0 || throw(ArgumentError("known Lnu759 must be nonnegative"))
    wavelength_slope_sigma >= 0 || throw(ArgumentError(
        "wavelength-space slope sigma must be nonnegative"))

    C = WAVENUMBER_CONVERSION_NM_CM1
    lambda = CORE_REFERENCE_WAVELENGTH_NM
    nu_anchor = C / KNOWN_SIF_WAVELENGTH_NM
    nu_reference = C / lambda
    delta_nu = nu_reference - nu_anchor
    beta = -2C * Float64(Lnu759) / lambda^3
    alpha = -C^2 / lambda^4 - 2C * delta_nu / lambda^3
    iszero(alpha) && error("degenerate SIF coordinate transformation")

    mean = (Float64(wavelength_slope_mean) - beta) / alpha
    sigma = Float64(wavelength_slope_sigma) / abs(alpha)
    SIF760_mean = Float64(Lnu759) + delta_nu * mean
    SIF760_sigma = abs(delta_nu) * sigma
    return (;
        mean,
        sigma,
        variance=sigma^2,
        SIF760_mean,
        SIF760_sigma,
        delta_nu_760_minus_759=delta_nu,
        wavelength_slope_alpha=alpha,
        wavelength_slope_beta=beta,
    )
end

function _active_indices(mode::Symbol)
    _check_mode(mode)
    return mode == :on ? copy(ROUND4_ON_ACTIVE_TO_FULL) :
                         copy(ROUND4_OFF_ACTIVE_TO_FULL)
end

function _active_core_indices(mode::Symbol)
    _check_mode(mode)
    return mode == :on ? copy(ROUND4_ON_ACTIVE_TO_CORE) :
                         copy(ROUND4_OFF_ACTIVE_TO_CORE)
end

"""
    build_round4_prior(source_path, mode)

Build one four-surface prior in memory. Every mean and covariance entry in
full-state coordinates 1:32 is copied exactly from the approved round-3
tapered prior. Only its SIF block (33:34) and explicit active mapping change.
"""
function build_round4_prior(source_path::AbstractString, mode::Symbol)
    _check_mode(mode)
    # Do not replace this construction with Gaussian conditioning of the old
    # 2-D SIF block. That block encoded uncertain amplitude at 760 nm, whereas
    # round 4 declares the 759-nm amplitude exact and explicitly retains the
    # physical wavelength-slope prior below.
    source = _read_source_prior(source_path)
    convention = validate_round4_sif_truth_convention()

    known_Lnu759 = mode == :on ? Float64(convention.diagnostic.Lnu) : 0.0
    known_Llambda759 = mode == :on ?
        Float64(convention.diagnostic.Llambda) : 0.0
    slope_mean = mode == :on ? source.slope_mean : 0.0
    slope_sigma = mode == :on ? source.slope_sigma : 0.0
    slope = native_slope_prior(
        known_Lnu759, slope_mean, slope_sigma)

    xa = copy(source.xa)
    Sa = copy(source.Sa)
    xa[CORE_SIF760_FULL_INDEX:CORE_MSIF_FULL_INDEX, :] .= 0.0
    Sa[CORE_SIF760_FULL_INDEX:CORE_MSIF_FULL_INDEX, :, :] .= 0.0
    Sa[:, CORE_SIF760_FULL_INDEX:CORE_MSIF_FULL_INDEX, :] .= 0.0

    if mode == :on
        delta_nu = slope.delta_nu_760_minus_759
        sif_covariance = slope.variance .* [
            delta_nu^2 delta_nu
            delta_nu   1.0
        ]
        for surface_index in axes(xa, 2)
            xa[CORE_SIF760_FULL_INDEX, surface_index] = slope.SIF760_mean
            xa[CORE_MSIF_FULL_INDEX, surface_index] = slope.mean
            Sa[CORE_SIF760_FULL_INDEX:CORE_MSIF_FULL_INDEX,
               CORE_SIF760_FULL_INDEX:CORE_MSIF_FULL_INDEX,
               surface_index] .= sif_covariance
        end
    end

    active_to_full = _active_indices(mode)
    active_to_core = _active_core_indices(mode)
    active_mask = falses(FULL_STATE_COUNT)
    active_mask[active_to_full] .= true
    Sa_active = Array{Float64}(undef,
        length(active_to_full), length(active_to_full), size(Sa, 3))
    for surface_index in axes(Sa, 3)
        Sa_active[:, :, surface_index] .=
            Sa[active_to_full, active_to_full, surface_index]
        isposdef(Symmetric(Sa_active[:, :, surface_index])) || error(
            "round-4 $mode active covariance is not positive definite for " *
            "surface $(source.surface_order[surface_index])")

        # Round 4 is allowed to change only the two SIF coordinates.
        xa[1:32, surface_index] == source.xa[1:32, surface_index] || error(
            "round-4 builder changed a non-SIF prior mean")
        Sa[1:32, 1:32, surface_index] ==
            source.Sa[1:32, 1:32, surface_index] || error(
            "round-4 builder changed a non-SIF prior covariance")
    end

    if mode == :on
        delta_nu = slope.delta_nu_760_minus_759
        for surface_index in axes(xa, 2)
            reconstructed = xa[CORE_SIF760_FULL_INDEX, surface_index] -
                delta_nu * xa[CORE_MSIF_FULL_INDEX, surface_index]
            isapprox(reconstructed, known_Lnu759;
                     atol=4eps(known_Lnu759), rtol=4eps(Float64)) || error(
                "round-4 prior mean violates its exact 759-nm SIF anchor")
        end
    else
        all(iszero, xa[33:34, :]) || error("SIF-off mean is not exact zero")
        all(iszero, Sa[33:34, :, :]) || error(
            "SIF-off covariance rows are not exact zero")
        all(iszero, Sa[:, 33:34, :]) || error(
            "SIF-off covariance columns are not exact zero")
    end

    return (;
        mode,
        source,
        xa,
        Sa,
        Sa_active,
        active_to_full,
        active_to_core,
        active_mask,
        known_Lnu759,
        known_Llambda759,
        wavelength_slope_mean=slope_mean,
        wavelength_slope_sigma=slope_sigma,
        native_slope=slope,
    )
end

function _copy_attributes!(target, attributes)
    for (name, value) in attributes
        target.attrib[name] = value
    end
    return target
end

function _define_and_copy_auxiliary!(output, prior, name::String,
                                     dimensions)
    values = prior.source.auxiliary[name]
    variable = defVar(output, name, eltype(values), dimensions)
    _copy_attributes!(variable,
        prior.source.auxiliary_attributes[name])
    variable[:] = values
    return variable
end

"""Write one round-4 prior NetCDF, refusing replacement by default."""
function write_round4_prior(prior;
                            output_path::AbstractString,
                            overwrite::Bool=false)
    isfile(output_path) && !overwrite && throw(ArgumentError(
        "round-4 prior already exists: $output_path"))
    mkpath(dirname(output_path))
    isfile(output_path) && rm(output_path)

    nactive = length(prior.active_to_full)
    nsurface = length(prior.source.surface_order)
    NCDataset(output_path, "c") do output
        defDim(output, "parameter", FULL_STATE_COUNT)
        defDim(output, "parameter_2", FULL_STATE_COUNT)
        defDim(output, "active_parameter", nactive)
        defDim(output, "active_parameter_2", nactive)
        defDim(output, "surface", nsurface)
        defDim(output, "sif_wavelength_parameter", 2)
        defDim(output, "co2_adjacent_layer_pair", 11)

        xa = defVar(output, "xa", Float64, ("parameter", "surface"))
        xa.attrib["long_name"] =
            "surface-specific full-state prior mean; fixed/derived entries retained"
        xa[:, :] = prior.xa

        covariance = defVar(output, "Sa", Float64,
            ("parameter", "parameter_2", "surface"))
        covariance.attrib["long_name"] =
            "full retrieval-coordinate covariance; singular at exact constraints"
        covariance[:, :, :] = prior.Sa

        sigma = defVar(output, "prior_sigma", Float64,
            ("parameter", "surface"))
        sigma.attrib["long_name"] =
            "square root of the full covariance diagonal"
        for surface_index in 1:nsurface
            sigma[:, surface_index] =
                sqrt.(max.(diag(prior.Sa[:, :, surface_index]), 0.0))
        end

        active_covariance = defVar(output, "Sa_active", Float64,
            ("active_parameter", "active_parameter_2", "surface"))
        active_covariance.attrib["long_name"] =
            "positive-definite covariance in the round-4 numerical solve"
        active_covariance[:, :, :] = prior.Sa_active

        active_index = defVar(output, "active_parameter_index", Int16,
            ("active_parameter",))
        active_index.attrib["index_convention"] =
            "one-based index into the canonical full 34-element state"
        active_index[:] = Int16.(prior.active_to_full)

        active_core_index = defVar(output, "active_core_parameter_index", Int16,
            ("active_parameter",))
        active_core_index.attrib["index_convention"] =
            "one-based index into the 30-column OCO_RRS_synth core state"
        active_core_index[:] = Int16.(prior.active_to_core)

        active_mask = defVar(output, "active_mask", Int8, ("parameter",))
        active_mask.attrib["flag_values"] = Int8[0, 1]
        active_mask.attrib["flag_meanings"] = "fixed_or_derived active"
        active_mask[:] = Int8.(prior.active_mask)

        _define_and_copy_auxiliary!(output, prior,
            "co2_layer_center_height", ("parameter",))
        _define_and_copy_auxiliary!(output, prior,
            "co2_layer_center_pressure", ("parameter",))
        _define_and_copy_auxiliary!(output, prior,
            "co2_correlation_taper_adjacent_retention",
            ("co2_adjacent_layer_pair",))
        _define_and_copy_auxiliary!(output, prior,
            "co2_base_adjacent_correlation",
            ("co2_adjacent_layer_pair",))
        _define_and_copy_auxiliary!(output, prior,
            "co2_selected_adjacent_correlation",
            ("co2_adjacent_layer_pair",))

        wavelength_state = defVar(output, "sif_wavelength_state", Float64,
            ("sif_wavelength_parameter", "surface"))
        wavelength_state.attrib["units"] =
            "mixed: mW m-2 sr-1 nm-1; mW m-2 sr-1 nm-2"
        wavelength_state.attrib["order"] =
            "known_Llambda_759 dLlambda_dlambda_760"
        wavelength_sigma = defVar(output, "sif_wavelength_sigma", Float64,
            ("sif_wavelength_parameter", "surface"))
        wavelength_sigma.attrib["order"] =
            "sigma_known_Llambda_759 sigma_dLlambda_dlambda_760"
        for surface_index in 1:nsurface
            wavelength_state[:, surface_index] =
                [prior.known_Llambda759, prior.wavelength_slope_mean]
            wavelength_sigma[:, surface_index] =
                [0.0, prior.wavelength_slope_sigma]
        end

        # Carry all round-3 physical/covariance provenance forward, then make
        # every round-4 change explicit below.
        _copy_attributes!(output, prior.source.attributes)
        output.attrib["retrieval_state_model"] = "round4_known_sif759"
        output.attrib["round4_prior_model"] = ROUND4_MODEL
        output.attrib["round4_prior_definition_version"] = Int32(1)
        output.attrib["round4_sif_case"] = String(prior.mode)
        output.attrib["known_sif_wavelength_nm"] = KNOWN_SIF_WAVELENGTH_NM
        output.attrib["known_sif_Lnu_mW_m-2_sr-1_per_cm-1"] =
            prior.known_Lnu759
        output.attrib["known_sif_Llambda_mW_m-2_sr-1_nm-1"] =
            prior.known_Llambda759
        output.attrib["known_sif_source"] = prior.mode == :on ?
            "corrected-v2 full truth-template interpolation at 759 nm" :
            "exact zero for the no-SIF truth class"
        output.attrib["core_sif_reference_wavelength_nm"] =
            CORE_REFERENCE_WAVELENGTH_NM
        output.attrib["round4_sif_state_mapping"] =
            "SIF760=Lnu759+mSIF*(nu760-nu759); Lnu759 exact; mSIF active only for SIF-on"
        output.attrib["round4_sif_jacobian_mapping"] =
            "K_mSIF_round4=K_SIF760*(nu760-nu759)+K_mSIF_core"
        output.attrib["round4_slope_constraint_coordinate"] =
            "dLlambda/dlambda at 760 nm"
        output.attrib["round4_slope_prior_mean_mW_m-2_sr-1_nm-2"] =
            prior.wavelength_slope_mean
        output.attrib["round4_slope_prior_sigma_mW_m-2_sr-1_nm-2"] =
            prior.wavelength_slope_sigma
        # Override the legacy generic fields as well, so a SIF-off file does
        # not simultaneously advertise the round-3 nonzero slope variance.
        output.attrib["sif_wavelength_slope_prior_mw_m2_sr_nm2"] =
            prior.wavelength_slope_mean
        output.attrib["sif_wavelength_slope_sigma_mw_m2_sr_nm2"] =
            prior.wavelength_slope_sigma
        output.attrib["sif_fractional_slope_sigma_per_nm"] =
            prior.mode == :on ? Float64(prior.source.attributes[
                "sif_fractional_slope_sigma_per_nm"]) : 0.0
        output.attrib["sif_slope_prior_interpretation"] =
            prior.mode == :on ?
            "absolute wavelength-space slope at 760 nm retained from round 3 while Lnu759 is exact" :
            "SIF-off slope is canonically zero with zero variance and is excluded from the solve"
        output.attrib["round4_native_mSIF_prior_mean_mW_m-2_sr-1_per_cm-2"] =
            prior.native_slope.mean
        output.attrib["round4_native_mSIF_prior_sigma_mW_m-2_sr-1_per_cm-2"] =
            prior.native_slope.sigma
        output.attrib["round4_native_SIF760_prior_mean_mW_m-2_sr-1_per_cm-1"] =
            prior.native_slope.SIF760_mean
        output.attrib["round4_native_SIF760_prior_sigma_mW_m-2_sr-1_per_cm-1"] =
            prior.native_slope.SIF760_sigma
        output.attrib["round4_delta_nu_760_minus_759_cm-1"] =
            prior.native_slope.delta_nu_760_minus_759
        output.attrib["source_prior_path"] = prior.source.path
        output.attrib["source_prior_sha256"] = prior.source.sha256
        output.attrib["source_prior_co2_covariance_model"] =
            String(prior.source.attributes["co2_covariance_model"])
        output.attrib["non_sif_prior_invariant"] =
            "xa[1:32] and Sa[1:32,1:32] copied exactly from source_prior_sha256"
        output.attrib["round4_production_prior_policy"] =
            "use this dedicated builder; do not Gaussian-condition the old uncertain-amplitude SIF block because that changes the retained physical slope prior"
        output.attrib["active_state_count"] = Int32(nactive)
        output.attrib["round4_active_to_full"] =
            join(prior.active_to_full, " ")
        output.attrib["round4_active_to_core"] =
            join(prior.active_to_core, " ")
        output.attrib["correlation_convention"] =
            "source non-SIF blocks unchanged; exact known-SIF constraint; one-dimensional SIF slope block for SIF-on"
        output.attrib["full_covariance_status"] = prior.mode == :on ?
            "positive semidefinite: four upper CO2 entries fixed and SIF760 derived from one active slope" :
            "positive semidefinite: four upper CO2 entries and both SIF coordinates fixed"
        output.attrib["created_utc"] = string(now(UTC))
        output.attrib["apriori_complete"] = Int32(1)
    end
    return abspath(output_path)
end

"""Write a human-readable audit table paired with one round-4 prior."""
function write_round4_summary(prior;
                              output_path::AbstractString,
                              overwrite::Bool=false)
    isfile(output_path) && !overwrite && throw(ArgumentError(
        "round-4 prior summary already exists: $output_path"))
    mkpath(dirname(output_path))
    isfile(output_path) && rm(output_path)

    active_position = Dict(index => position for
        (position, index) in enumerate(prior.active_to_full))
    core_position = Dict(index => core for
        (index, core) in zip(prior.active_to_full, prior.active_to_core))
    open(output_path, "w") do io
        println(io, "# Round-4 known-SIF a priori: ", ROUND4_MODEL)
        println(io, "# SIF case: ", prior.mode)
        println(io, "# Source tapered prior: ", prior.source.path)
        println(io, "# Source prior SHA-256: ", prior.source.sha256)
        println(io, "# Production policy: generated only by the dedicated round-4 builder; Gaussian conditioning of the old uncertain-amplitude SIF block is not equivalent and must not be substituted.")
        @printf(io, "# Exact SIF anchor: lambda=%.1f nm Lnu=%.15e mW m-2 sr-1 (cm-1)-1 Llambda=%.15e mW m-2 sr-1 nm-1\n",
                KNOWN_SIF_WAVELENGTH_NM, prior.known_Lnu759,
                prior.known_Llambda759)
        @printf(io, "# Retained wavelength-space slope at 760 nm: mean=% .15e sigma=%.15e mW m-2 sr-1 nm-2\n",
                prior.wavelength_slope_mean, prior.wavelength_slope_sigma)
        @printf(io, "# Native mSIF prior: mean=% .15e sigma=%.15e mW m-2 sr-1 (cm-1)-2\n",
                prior.native_slope.mean, prior.native_slope.sigma)
        @printf(io, "# Derived SIF760 prior: mean=% .15e sigma=%.15e mW m-2 sr-1 (cm-1)-1\n",
                prior.native_slope.SIF760_mean,
                prior.native_slope.SIF760_sigma)
        println(io, "# Active state count: ", length(prior.active_to_full))
        println(io, "# Active full-state indices: ",
                join(prior.active_to_full, " "))
        println(io, "# Active OCO_RRS_synth core indices: ",
                join(prior.active_to_core, " "))
        println(io, "# surface full_index active_index core_index parameter mean sigma units role")
        for (surface_index, surface) in enumerate(prior.source.surface_order)
            for full_index in 1:FULL_STATE_COUNT
                active_index = get(active_position, full_index, 0)
                core_index = get(core_position, full_index, 0)
                role = full_index == CORE_SIF760_FULL_INDEX ?
                    (prior.mode == :on ? "derived_from_anchor_and_slope" : "fixed_zero") :
                    full_index == CORE_MSIF_FULL_INDEX ?
                    (prior.mode == :on ? "active_slope" : "fixed_zero") :
                    (active_index == 0 ? "fixed" : "active")
                @printf(io, "%-7s %2d %2d %2d %-29s % .15e %.15e %-29s %s\n",
                        surface, full_index, active_index, core_index,
                        prior.source.parameter_names[full_index],
                        prior.xa[full_index, surface_index],
                        sqrt(max(prior.Sa[full_index, full_index,
                                          surface_index], 0.0)),
                        prior.source.parameter_units[full_index], role)
            end
        end
        println(io, "# Native SIF covariance blocks")
        for (surface_index, surface) in enumerate(prior.source.surface_order)
            block = prior.Sa[33:34, 33:34, surface_index]
            @printf(io, "# %-7s Sa33_33=%.15e Sa33_34=% .15e Sa34_34=%.15e\n",
                    surface, block[1, 1], block[1, 2], block[2, 2])
        end
    end
    return abspath(output_path)
end

function main(args=ARGS)
    isempty(args) || error(
        "usage: build_round4_known_sif_apriori.jl (configure paths with environment variables)")
    source = get(ENV, "ROUND4_SOURCE_PRIOR", DEFAULT_SOURCE_PRIOR)
    output_directory = get(
        ENV, "ROUND4_PRIOR_OUTPUT_DIRECTORY", DEFAULT_OUTPUT_DIRECTORY)
    overwrite = get(ENV, "ROUND4_PRIOR_OVERWRITE", "0") == "1"

    for mode in (:off, :on)
        prior = build_round4_prior(source, mode)
        netcdf_path = joinpath(output_directory,
            output_filename(mode; extension="nc"))
        summary_path = joinpath(output_directory,
            output_filename(mode; extension="dat"))
        write_round4_prior(prior;
            output_path=netcdf_path, overwrite=overwrite)
        write_round4_summary(prior;
            output_path=summary_path, overwrite=overwrite)
        println(netcdf_path)
        println(summary_path)
        println(file_sha256(netcdf_path), "  ", netcdf_path)
    end
end

end # module Round4KnownSIFPriorBuilder

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Round4KnownSIFPriorBuilder.main()
end
