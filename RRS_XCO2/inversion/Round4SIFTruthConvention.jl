module Round4SIFTruthConvention

include(joinpath(@__DIR__, "..", "scripts", "common.jl"))
using .RRSXCO2Common

export ROUND4_SIF_DIAGNOSTIC_WAVELENGTH_NM,
       lnu_to_llambda, llambda_to_lnu,
       truth_sif_at_wavelength_nm, round4_sif_truth_convention,
       expected_round4_sif_provenance,
       validate_round4_sif_truth_convention,
       validate_round4_sif_provenance

const ROUND4_SIF_DIAGNOSTIC_WAVELENGTH_NM = 759.0
const WAVENUMBER_PER_WAVELENGTH = 1.0e7

"""
    lnu_to_llambda(Lnu, wavelength_nm)

Convert spectral radiance density per cm⁻¹ to spectral radiance density per
nm.  With `nu = 1e7 / wavelength_nm`, the density Jacobian is
`|dnu/dlambda| = 1e7 / wavelength_nm^2`.
"""
function lnu_to_llambda(Lnu::Real, wavelength_nm::Real)
    wavelength_nm > 0 || throw(ArgumentError(
        "wavelength_nm must be positive"))
    return Lnu * WAVENUMBER_PER_WAVELENGTH / wavelength_nm^2
end

"""Inverse of [`lnu_to_llambda`](@ref)."""
function llambda_to_lnu(Llambda::Real, wavelength_nm::Real)
    wavelength_nm > 0 || throw(ArgumentError(
        "wavelength_nm must be positive"))
    return Llambda * wavelength_nm^2 / WAVENUMBER_PER_WAVELENGTH
end

function _interpolate_sorted(x::AbstractVector, y::AbstractVector, query::Real)
    length(x) == length(y) || throw(DimensionMismatch(
        "spectral coordinate and radiance arrays have different lengths"))
    length(x) >= 2 || throw(ArgumentError(
        "at least two SIF template samples are required"))
    issorted(x) || throw(ArgumentError(
        "SIF template wavenumbers must be monotonically increasing"))
    first(x) <= query <= last(x) || throw(ArgumentError(
        "requested wavenumber $query lies outside the SIF template"))

    upper = searchsortedfirst(x, query)
    upper <= length(x) && x[upper] == query && return y[upper]
    upper > 1 || return y[1]
    lower = upper - 1
    fraction = (query - x[lower]) / (x[upper] - x[lower])
    return y[lower] + fraction * (y[upper] - y[lower])
end

"""
    truth_sif_at_wavelength_nm(wavelength_nm=759.0; state=campaign_sif_state())

Evaluate the *full corrected-v2 truth template* at a wavelength.  `Lnu` is in
`mW m^-2 sr^-1 (cm^-1)^-1` and `Llambda` is in
`mW m^-2 sr^-1 nm^-1`.

The additional `linear_Lnu` and `linear_Llambda` fields evaluate the retrieval
forward model's local two-parameter tangent at the same wavelength.  They are
reported to make the deliberate truth-template/retrieval-model distinction
explicit; they must not replace `Lnu`/`Llambda` when specifying truth.
"""
function truth_sif_at_wavelength_nm(
        wavelength_nm::Real=ROUND4_SIF_DIAGNOSTIC_WAVELENGTH_NM;
        state=RRSXCO2Common.campaign_sif_state())
    wavelength_nm > 0 || throw(ArgumentError(
        "wavelength_nm must be positive"))
    nu = WAVENUMBER_PER_WAVELENGTH / wavelength_nm
    Lnu = _interpolate_sorted(state.ν, state.spectrum, nu)
    linear_Lnu = state.SIF760 + state.mSIF * (nu - state.ν_ref)
    return (;
        wavelength_nm=Float64(wavelength_nm),
        wavenumber_cm1=nu,
        Lnu,
        Llambda=lnu_to_llambda(Lnu, wavelength_nm),
        linear_Lnu,
        linear_Llambda=lnu_to_llambda(linear_Lnu, wavelength_nm),
    )
end

"""
    round4_sif_truth_convention(; state=campaign_sif_state())

Return the immutable scalar facts needed to pin the corrected-v2 SIF truth
convention in a round-4 campaign.  The requested `0.5` is the unweighted
upward-solid-angle integral `2pi*Llambda(760 nm)`, not `SIF760` and not the
cosine-weighted hemispheric irradiance.
"""
function round4_sif_truth_convention(;
        state=RRSXCO2Common.campaign_sif_state())
    diagnostic = truth_sif_at_wavelength_nm(
        ROUND4_SIF_DIAGNOSTIC_WAVELENGTH_NM; state)
    return (;
        definition_version=RRSXCO2Common.SIF_DEFINITION_VERSION,
        case_on=RRSXCO2Common.SIF_CASE_ON,
        reference_wavelength_nm=
            RRSXCO2Common.SIF_REFERENCE_WAVELENGTH_NM,
        upwelling_solid_angle_sr=
            RRSXCO2Common.SIF_UPWELLING_SOLID_ANGLE_SR,
        angular_integral_760=
            RRSXCO2Common.SIF_ANGULAR_INTEGRAL_760,
        radiance_Llambda_760=state.radiance_760,
        cosine_weighted_irradiance_760=pi * state.radiance_760,
        SIF760=state.SIF760,
        mSIF=state.mSIF,
        template_wavelength_integral=state.wavelength_integral,
        diagnostic,
    )
end

"""
    validate_round4_sif_truth_convention([convention])

Fail closed unless the convention has the exact corrected-v2 campaign
identity and is internally consistent in wavelength and wavenumber units.
This function performs no I/O and does not modify truth products.
"""
function validate_round4_sif_truth_convention(
        convention=round4_sif_truth_convention())
    convention.definition_version == 2 || error(
        "round-4 SIF requires definition version 2")
    convention.case_on == "angular_integral760_0p5" || error(
        "round-4 SIF uses a stale case label")
    convention.reference_wavelength_nm == 760.0 || error(
        "round-4 SIF reference wavelength must be 760 nm")
    isapprox(convention.upwelling_solid_angle_sr, 2pi;
             atol=0, rtol=8eps(Float64)) || error(
        "round-4 SIF has an inconsistent upward solid angle")
    isapprox(convention.angular_integral_760, 0.5;
             atol=0, rtol=0) || error(
        "round-4 SIF angular integral at 760 nm must be 0.5")
    isapprox(convention.radiance_Llambda_760, 0.5 / (2pi);
             atol=2e-16, rtol=0) || error(
        "round-4 SIF stream radiance at 760 nm is inconsistent")
    isapprox(convention.SIF760,
             llambda_to_lnu(convention.radiance_Llambda_760, 760.0);
             atol=2e-18, rtol=0) || error(
        "round-4 SIF760 has inconsistent spectral-density units")
    isapprox(convention.SIF760, 0.004596394756493938;
             atol=2e-18, rtol=0) || error(
        "round-4 SIF760 does not match the corrected-v2 truth")
    isapprox(convention.mSIF, 1.2291230681458325e-5;
             atol=2e-19, rtol=0) || error(
        "round-4 mSIF does not match the corrected-v2 truth template")
    isapprox(convention.cosine_weighted_irradiance_760,
             pi * convention.radiance_Llambda_760;
             atol=2e-16, rtol=0) || error(
        "round-4 cosine-weighted SIF irradiance is inconsistent")

    diagnostic = convention.diagnostic
    diagnostic.wavelength_nm == 759.0 || error(
        "round-4 SIF diagnostic must be evaluated at 759 nm")
    isapprox(diagnostic.Lnu, 0.004818031987713776;
             atol=3e-18, rtol=0) || error(
        "759-nm truth-template Lnu has drifted")
    isapprox(diagnostic.Llambda, 0.0836346275560863;
             atol=3e-17, rtol=0) || error(
        "759-nm truth-template Llambda has drifted")
    isapprox(llambda_to_lnu(diagnostic.Llambda, 759.0), diagnostic.Lnu;
             atol=3e-18, rtol=0) || error(
        "759-nm wavelength/wavenumber density conversion is inconsistent")
    return convention
end

"""
    expected_round4_sif_provenance(; enabled=true)

Construct the versioned NetCDF provenance record expected for a corrected-v2
truth product.  This only builds an in-memory dictionary.
"""
function expected_round4_sif_provenance(; enabled::Bool=true)
    attributes = Dict{String,Any}()
    RRSXCO2Common.write_sif_provenance!(attributes, enabled)
    return attributes
end

"""
    validate_round4_sif_provenance(attributes; enabled=true, source="input")

Validate all corrected-v2 SIF provenance fields and their unit identities.
SIF-off records use the same versioned definition with zero-valued amplitudes.
"""
function validate_round4_sif_provenance(attributes;
                                        enabled::Bool=true,
                                        source::AbstractString="input")
    expected = expected_round4_sif_provenance(; enabled)
    for (key, wanted) in expected
        haskey(attributes, key) || error(
            "$source is missing corrected-v2 SIF provenance '$key'")
        actual = attributes[key]
        matches = wanted isa Real && actual isa Real ?
            isapprox(Float64(actual), Float64(wanted);
                     atol=1e-14, rtol=1e-12) :
            string(actual) == string(wanted)
        matches || error(
            "$source has incorrect corrected-v2 SIF provenance '$key'")
    end

    active = enabled ? 1.0 : 0.0
    Llambda = Float64(
        attributes["sif_radiance_760_mW_m-2_sr-1_nm-1"])
    angular = Float64(
        attributes["sif_angular_integral_760_mW_m-2_nm-1"])
    irradiance = Float64(
        attributes["sif_cosine_weighted_irradiance_760_mW_m-2_nm-1"])
    Lnu = Float64(
        attributes["sif_SIF760_mW_m-2_sr-1_per_cm-1"])
    isapprox(angular, active * 2pi * Llambda;
             atol=1e-14, rtol=1e-12) || error(
        "$source has inconsistent angular-integral SIF provenance")
    isapprox(irradiance, pi * Llambda;
             atol=1e-14, rtol=1e-12) || error(
        "$source has inconsistent cosine-weighted SIF provenance")
    isapprox(Lnu, llambda_to_lnu(Llambda, 760.0);
             atol=1e-14, rtol=1e-12) || error(
        "$source has inconsistent wavelength/wavenumber SIF provenance")
    return Dict{String,Any}(key => attributes[key] for key in keys(expected))
end

end # module
