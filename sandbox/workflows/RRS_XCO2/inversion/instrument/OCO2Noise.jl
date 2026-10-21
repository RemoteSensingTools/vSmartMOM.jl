module OCO2Noise

using NCDatasets

export NoiseBandSpec,
       NOISE_BAND_SPECS,
       noise_band_spec,
       energy_per_photon,
       energy_to_photon_radiance,
       photon_to_energy_radiance,
       noise_equivalent_radiance_photon,
       noise_statistics,
       read_representative_snr_coefficients

const PLANCK = 6.62607015e-34
const LIGHT_SPEED = 299792458.0

"""OCO-2 dynamic-range constants from L1B ATBD Tables 3-5 and 3-6."""
struct NoiseBandSpec{FT<:AbstractFloat}
    name::Symbol
    max_ms::FT
    min_ms::FT
end

const NOISE_BAND_SPECS = (
    NoiseBandSpec(:o2a, 7.00e20, 7.50e16),
    NoiseBandSpec(:weak_co2, 2.45e20, 2.15e16),
    NoiseBandSpec(:strong_co2, 1.25e20, 2.15e16),
)

function noise_band_spec(name::Symbol)
    index = findfirst(spec -> spec.name == name, NOISE_BAND_SPECS)
    isnothing(index) && throw(ArgumentError("unknown OCO-2 noise band: $name"))
    return NOISE_BAND_SPECS[index]
end

"""Photon energy in joules at wavelengths expressed in nanometres."""
function energy_per_photon(wavelength_nm::AbstractVector)
    all(isfinite, wavelength_nm) ||
        throw(ArgumentError("wavelength contains non-finite values"))
    all(wavelength_nm .> 0) ||
        throw(ArgumentError("wavelength must be positive"))
    return PLANCK * LIGHT_SPEED ./ (Float64.(wavelength_nm) .* 1e-9)
end

"""
Convert mW m^-2 sr^-1 nm^-1 to photons s^-1 m^-2 sr^-1 um^-1.

The numerical value in mW/nm is identical to W/um, so division by the photon
energy is the only numerical conversion required.
"""
function energy_to_photon_radiance(wavelength_nm::AbstractVector,
                                   radiance_mw_per_nm::AbstractVector)
    length(wavelength_nm) == length(radiance_mw_per_nm) ||
        throw(DimensionMismatch("wavelength and radiance lengths differ"))
    all(isfinite, radiance_mw_per_nm) ||
        throw(ArgumentError("radiance contains non-finite values"))
    return Float64.(radiance_mw_per_nm) ./ energy_per_photon(wavelength_nm)
end

"""Inverse of [`energy_to_photon_radiance`](@ref)."""
function photon_to_energy_radiance(wavelength_nm::AbstractVector,
                                   radiance_photon::AbstractVector)
    length(wavelength_nm) == length(radiance_photon) ||
        throw(DimensionMismatch("wavelength and radiance lengths differ"))
    all(isfinite, radiance_photon) ||
        throw(ArgumentError("photon radiance contains non-finite values"))
    return Float64.(radiance_photon) .* energy_per_photon(wavelength_nm)
end

"""
Evaluate OCO-2 L1B ATBD Eq. (3-8) in photon-radiance units.

`signal_photon` and `max_ms` have units photons s^-1 m^-2 sr^-1 um^-1.
`c_photon` and `c_background` are the two wavelength-dependent coefficients
from `InstrumentHeader/snr_coef`. The absolute value inside Eq. (3-8) is
retained exactly.
"""
function noise_equivalent_radiance_photon(signal_photon::AbstractVector,
                                          max_ms::Real,
                                          c_photon::AbstractVector,
                                          c_background::AbstractVector)
    length(signal_photon) == length(c_photon) == length(c_background) ||
        throw(DimensionMismatch("signal and SNR coefficient lengths differ"))
    max_ms > 0 || throw(ArgumentError("MaxMS must be positive"))
    all(isfinite, signal_photon) ||
        throw(ArgumentError("signal contains non-finite values"))
    all(isfinite, c_photon) && all(isfinite, c_background) ||
        throw(ArgumentError("SNR coefficients contain non-finite values"))
    all(c_photon .>= 0) && all(c_background .>= 0) ||
        throw(ArgumentError("SNR coefficients must be non-negative"))

    normalized_signal = abs.(100 .* Float64.(signal_photon) ./ Float64(max_ms))
    return Float64(max_ms) / 100 .* sqrt.(
        normalized_signal .* Float64.(c_photon) .^ 2 .+
        Float64.(c_background) .^ 2)
end

"""
Return per-sample OCO-2 noise diagnostics and the diagonal covariance.

The requested diagonal measurement-error covariance is represented exactly by
`variance_energy`; a dense matrix is `Diagonal(variance_energy)`. MinMS is
used only to flag samples outside the Table 3-6 requirement. It is not a
clamping value and does not enter ATBD Eq. (3-8).
"""
function noise_statistics(wavelength_nm::AbstractVector,
                          radiance_mw_per_nm::AbstractVector,
                          c_photon::AbstractVector,
                          c_background::AbstractVector,
                          spec::NoiseBandSpec)
    signal_photon = energy_to_photon_radiance(
        wavelength_nm, radiance_mw_per_nm)
    nen_photon = noise_equivalent_radiance_photon(
        signal_photon, spec.max_ms, c_photon, c_background)
    nen_energy = photon_to_energy_radiance(wavelength_nm, nen_photon)
    variance_energy = nen_energy .^ 2
    snr = abs.(signal_photon) ./ nen_photon
    return (
        signal_photon=signal_photon,
        nen_photon=nen_photon,
        nen_energy=nen_energy,
        variance_energy=variance_energy,
        snr=snr,
        below_min_ms=abs.(signal_photon) .< spec.min_ms,
        above_max_ms=abs.(signal_photon) .> spec.max_ms,
    )
end

"""Read representative wavelength-dependent OCO-2 SNR coefficients."""
function read_representative_snr_coefficients(path::AbstractString)
    isfile(path) || throw(ArgumentError(
        "missing representative SNR coefficient file: $path"))
    return NCDataset(path) do dataset
        Dict(spec.name => begin
            name = String(spec.name)
            wavelength = Float64.(nomissing(
                dataset["$(name)_wavelength"][:], NaN))
            c_photon = Float64.(nomissing(
                dataset["c_photon_$(name)"][:], NaN))
            c_background = Float64.(nomissing(
                dataset["c_background_$(name)"][:], NaN))
            extrapolated_source_count = Int16.(nomissing(
                dataset["extrapolated_source_count_$(name)"][:], typemax(Int16)))
            all(isfinite, wavelength) && all(isfinite, c_photon) &&
                all(isfinite, c_background) || error(
                    "representative SNR file contains non-finite $name values")
            (
                wavelength=wavelength,
                c_photon=c_photon,
                c_background=c_background,
                extrapolated_source_count=extrapolated_source_count,
            )
        end for spec in NOISE_BAND_SPECS)
    end
end

end # module OCO2Noise
