module SyntheticOCO2

using NCDatasets
using Statistics

export BandSpec,
       BAND_SPECS,
       band_spec,
       fwhm_to_sigma,
       synthetic_grid,
       project_oco_analyzer,
       per_wavenumber_to_per_wavelength,
       required_wavenumber_extensions,
       gaussian_convolve_resample,
       process_stokes_spectrum,
       process_stokes_jacobian,
       read_representative_coefficients

"""
Synthetic OCO-2 band definition in wavelength space.

The limits, sampling intervals, and FWHM values are the OCO-2 entries in
Table 1 of `RRS_XCO2/resources/The_Cross-Calibration_of_Spectral_Radiances_and_Cr.pdf`.
The real OCO-2 ILS is deliberately replaced by a wavelength-space Gaussian.
"""
struct BandSpec{FT<:AbstractFloat}
    name::Symbol
    long_name::String
    wavelength_start_nm::FT
    wavelength_stop_nm::FT
    sampling_interval_nm::FT
    fwhm_nm::FT
end

const BAND_SPECS = (
    BandSpec(:o2a, "O2 A-band", 758.0, 772.0, 0.015, 0.04),
    BandSpec(:weak_co2, "weak CO2 band", 1594.0, 1619.0, 0.031, 0.08),
    BandSpec(:strong_co2, "strong CO2 band", 2042.0, 2082.0, 0.04, 0.10),
)

function band_spec(name::Symbol)
    index = findfirst(spec -> spec.name == name, BAND_SPECS)
    isnothing(index) && throw(ArgumentError("unknown synthetic OCO-2 band: $name"))
    return BAND_SPECS[index]
end

fwhm_to_sigma(fwhm) = fwhm / (2 * sqrt(2 * log(2)))

"""Return the requested `start:sampling:stop` wavelength grid for one band."""
synthetic_grid(spec::BandSpec) = collect(
    spec.wavelength_start_nm:spec.sampling_interval_nm:spec.wavelength_stop_nm)

"""
    project_oco_analyzer(stokes, coefficients; normalize=:none)

Apply the OCO L1B Stokes analyzer coefficients to a high-resolution Stokes
array whose layout is `(stokes, wavelength)`. The current truth simulations
contain `(I,Q,U)`. Following `OCORaman/src/OCOPlots/oco_gain.jl`, the scalar
detector response is

```
M11*I - M12*Q + M13*U.
```

With `normalize=:none` (the production default), return the raw detector
response exactly as in `oco_gain.jl`. `normalize=:m11` remains available as a
diagnostic that divides by the unpolarized throughput, but is not used for the
retrieval measurements.

The OCO L1B field contains a four-element analyzer row, not a complete 4x4
Mueller matrix. A determinant is therefore undefined. Moreover, the ideal
linear-polarizer matrix implied by these rows is singular. Consequently,
`:determinant` is rejected explicitly instead of silently inventing a
normalization.
"""
function project_oco_analyzer(stokes::AbstractMatrix,
                              coefficients::AbstractVector;
                              normalize::Symbol=:none)
    size(stokes, 1) == 3 || throw(DimensionMismatch(
        "truth-map analyzer projection expects three Stokes rows (I,Q,U); " *
        "received $(size(stokes, 1))"))
    length(coefficients) >= 3 || throw(DimensionMismatch(
        "OCO analyzer requires at least the I, Q, and U coefficients"))

    m11, m12, m13 = coefficients[1], coefficients[2], coefficients[3]
    isfinite(m11) && !iszero(m11) || throw(ArgumentError(
        "the representative OCO M11 coefficient must be finite and nonzero"))

    response = m11 .* view(stokes, 1, :) .-
               m12 .* view(stokes, 2, :) .+
               m13 .* view(stokes, 3, :)

    if normalize == :m11
        return response ./ m11
    elseif normalize == :none
        return response
    elseif normalize == :determinant
        throw(ArgumentError(
            "determinant normalization is undefined: OCO L1B stores only " *
            "a Stokes analyzer row, and its ideal-polarizer completion has det(M)=0"))
    else
        throw(ArgumentError("unsupported analyzer normalization: $normalize"))
    end
end

"""
Convert spectral radiance per cm^-1 to spectral radiance per nm.

For `nu = 1e7/lambda_nm`, the positive density Jacobian is
`abs(dnu/dlambda_nm) = 1e7/lambda_nm^2`.
"""
function per_wavenumber_to_per_wavelength(wavelength_nm::AbstractVector,
                                          radiance_per_cm::AbstractVector)
    length(wavelength_nm) == length(radiance_per_cm) ||
        throw(DimensionMismatch("wavelength and radiance lengths differ"))
    return radiance_per_cm .* (1e7 ./ wavelength_nm .^ 2)
end

function _ascending_spectrum(wavelength_nm::AbstractVector,
                             spectrum::AbstractVector)
    length(wavelength_nm) == length(spectrum) ||
        throw(DimensionMismatch("wavelength and spectrum lengths differ"))
    all(isfinite, wavelength_nm) ||
        throw(ArgumentError("source wavelengths contain non-finite values"))
    all(isfinite, spectrum) ||
        throw(ArgumentError("source spectrum contains non-finite/fill values"))

    order = sortperm(wavelength_nm)
    wavelength = Float64.(wavelength_nm[order])
    values = Float64.(spectrum[order])
    all(diff(wavelength) .> 0) ||
        throw(ArgumentError("source wavelengths must be unique"))
    return wavelength, values
end

function _trapezoid_weights(x::AbstractVector)
    length(x) >= 2 || throw(ArgumentError("at least two source wavelengths are required"))
    weights = similar(x, Float64)
    weights[1] = (x[2] - x[1]) / 2
    @views weights[2:end-1] .= (x[3:end] .- x[1:end-2]) ./ 2
    weights[end] = (x[end] - x[end-1]) / 2
    return weights
end

function _check_gaussian_support(source_wavelength_nm,
                                 target_wavelength_nm,
                                 sigma_nm,
                                 support_sigma)
    radius = support_sigma * sigma_nm
    required_min = minimum(target_wavelength_nm) - radius
    required_max = maximum(target_wavelength_nm) + radius
    available_min, available_max = extrema(source_wavelength_nm)
    if available_min > required_min || available_max < required_max
        throw(ArgumentError(
            "insufficient high-resolution wavelength support for a " *
            "$(support_sigma)-sigma Gaussian: available " *
            "[$available_min, $available_max] nm, required " *
            "[$required_min, $required_max] nm"))
    end
end

"""
    required_wavenumber_extensions(base_wavenumber, spec;
                                   step_cm=0.1, support_sigma=6)

Return the 0.1 cm^-1 nodes required beyond either end of an existing ascending
wavenumber grid to cover every synthetic sample through +/-`support_sigma`.
The `short` field contains higher-wavenumber/shorter-wavelength nodes and
`long` contains lower-wavenumber/longer-wavelength nodes.
"""
function required_wavenumber_extensions(base_wavenumber::AbstractVector,
                                        spec::BandSpec;
                                        step_cm::Real=0.1,
                                        support_sigma::Real=6)
    length(base_wavenumber) >= 2 ||
        throw(ArgumentError("base wavenumber grid needs at least two points"))
    all(diff(base_wavenumber) .> 0) ||
        throw(ArgumentError("base wavenumber grid must be strictly ascending"))
    step_cm > 0 || throw(ArgumentError("wavenumber step must be positive"))
    support_sigma > 0 || throw(ArgumentError("support_sigma must be positive"))

    target = synthetic_grid(spec)
    sigma_nm = fwhm_to_sigma(spec.fwhm_nm)
    required_short_nm = minimum(target) - support_sigma * sigma_nm
    required_long_nm = maximum(target) + support_sigma * sigma_nm
    available_wavelength = 1e7 ./ Float64.(base_wavenumber)

    short_nodes = Float64[]
    if minimum(available_wavelength) > required_short_nm
        required_max_wavenumber = 1e7 / required_short_nm
        count = ceil(Int, (required_max_wavenumber - Float64(last(base_wavenumber))) /
                          Float64(step_cm))
        count > 0 && (short_nodes = Float64(last(base_wavenumber)) .+
                                    Float64(step_cm) .* collect(1:count))
    end

    long_nodes = Float64[]
    if maximum(available_wavelength) < required_long_nm
        required_min_wavenumber = 1e7 / required_long_nm
        count = ceil(Int, (Float64(first(base_wavenumber)) - required_min_wavenumber) /
                          Float64(step_cm))
        count > 0 && (long_nodes = Float64(first(base_wavenumber)) .-
                                   Float64(step_cm) .* collect(reverse(1:count)))
    end

    return (short=short_nodes,
            long=long_nodes,
            required_short_nm=required_short_nm,
            required_long_nm=required_long_nm)
end

"""
    gaussian_convolve_resample(source_wavelength_nm, source_spectrum,
                               target_wavelength_nm, fwhm_nm;
                               support_sigma=6, require_full_support=true)

Convolve a spectral-density spectrum with a normalized wavelength-space
Gaussian and evaluate it directly at target sample centers. The source grid
may be nonuniform and may run in either direction. Trapezoidal wavelength
weights are included in both numerator and normalization, which preserves a
constant spectrum and correctly handles the truth map's uniform-wavenumber
(therefore nonuniform-wavelength) basis.

By default, every target center must have source coverage through +/-6 sigma.
This prevents a clipped edge kernel from being silently renormalized into a
biased measurement.
"""
function gaussian_convolve_resample(source_wavelength_nm::AbstractVector,
                                    source_spectrum::AbstractVector,
                                    target_wavelength_nm::AbstractVector,
                                    fwhm_nm;
                                    support_sigma::Real=6,
                                    require_full_support::Bool=true)
    fwhm_nm > 0 || throw(ArgumentError("Gaussian FWHM must be positive"))
    support_sigma > 0 || throw(ArgumentError("support_sigma must be positive"))
    all(isfinite, target_wavelength_nm) ||
        throw(ArgumentError("target wavelengths contain non-finite values"))

    wavelength, spectrum = _ascending_spectrum(source_wavelength_nm, source_spectrum)
    sigma_nm = fwhm_to_sigma(float(fwhm_nm))
    require_full_support && _check_gaussian_support(
        wavelength, target_wavelength_nm, sigma_nm, support_sigma)

    quadrature_weights = _trapezoid_weights(wavelength)
    radius = support_sigma * sigma_nm
    output = Vector{Float64}(undef, length(target_wavelength_nm))

    for (index, center) in pairs(target_wavelength_nm)
        first_index = searchsortedfirst(wavelength, center - radius)
        last_index = searchsortedlast(wavelength, center + radius)
        first_index <= last_index || throw(ArgumentError(
            "no source samples overlap the Gaussian centered at $center nm"))

        source_range = first_index:last_index
        delta = @view(wavelength[source_range]) .- center
        kernel = exp.(-0.5 .* (delta ./ sigma_nm) .^ 2)
        weights = kernel .* @view(quadrature_weights[source_range])
        normalization = sum(weights)
        normalization > 0 || throw(ArgumentError(
            "Gaussian normalization vanished at $center nm"))
        output[index] = sum(weights .* @view(spectrum[source_range])) / normalization
    end
    return output
end

"""
Apply the complete synthetic OCO-2 measurement operator to one Stokes
spectral component: analyzer projection, density conversion, Gaussian
convolution, then sampling at the synthetic OCO wavelengths.
"""
function process_stokes_spectrum(source_wavelength_nm::AbstractVector,
                                 stokes_per_cm::AbstractMatrix,
                                 coefficients::AbstractVector,
                                 spec::BandSpec;
                                 normalize::Symbol=:none,
                                 support_sigma::Real=6,
                                 require_full_support::Bool=true)
    size(stokes_per_cm, 2) == length(source_wavelength_nm) ||
        throw(DimensionMismatch("Stokes spectral dimension does not match wavelength"))
    detector_per_cm = project_oco_analyzer(
        stokes_per_cm, coefficients; normalize=normalize)
    detector_per_nm = per_wavenumber_to_per_wavelength(
        source_wavelength_nm, detector_per_cm)
    return gaussian_convolve_resample(
        source_wavelength_nm,
        detector_per_nm,
        synthetic_grid(spec),
        spec.fwhm_nm;
        support_sigma=support_sigma,
        require_full_support=require_full_support,
    )
end

"""
    process_stokes_jacobian(source_wavelength_nm, jacobian_per_cm,
                            coefficients, spec; ...)

Apply the same fixed measurement operator to every parameter column of a
canonical `(stokes, wavelength, parameter)` high-resolution Jacobian. The
returned layout is `(sample, parameter)`, matching the usual optimal-
estimation convention. No derivative of the instrument operator is needed
because the synthetic grid, Gaussian ILS, and representative analyzer row are
held fixed.
"""
function process_stokes_jacobian(source_wavelength_nm::AbstractVector,
                                 jacobian_per_cm::AbstractArray{<:Real,3},
                                 coefficients::AbstractVector,
                                 spec::BandSpec;
                                 normalize::Symbol=:none,
                                 support_sigma::Real=6,
                                 require_full_support::Bool=true)
    size(jacobian_per_cm, 1) == 3 || throw(DimensionMismatch(
        "Jacobian layout must be (3 Stokes, wavelength, parameter)"))
    size(jacobian_per_cm, 2) == length(source_wavelength_nm) ||
        throw(DimensionMismatch("Jacobian spectral dimension does not match wavelength"))

    output = Matrix{Float64}(
        undef, length(synthetic_grid(spec)), size(jacobian_per_cm, 3))
    for parameter_index in axes(jacobian_per_cm, 3)
        output[:, parameter_index] = process_stokes_spectrum(
            source_wavelength_nm,
            @view(jacobian_per_cm[:, :, parameter_index]),
            coefficients,
            spec;
            normalize=normalize,
            support_sigma=support_sigma,
            require_full_support=require_full_support,
        )
    end
    return output
end

"""Read the `(coefficient, band)` representative analyzer array."""
function read_representative_coefficients(path::AbstractString)
    isfile(path) || throw(ArgumentError("missing representative analyzer file: $path"))
    return NCDataset(path) do dataset
        variable = dataset["representative_stokes_coefficients"]
        coefficients = Array(variable[:, :])
        size(coefficients) == (4, length(BAND_SPECS)) ||
            throw(DimensionMismatch(
                "expected representative coefficients with size (4,3), got " *
                "$(size(coefficients))"))
        Dict(spec.name => Float64.(coefficients[:, index])
             for (index, spec) in pairs(BAND_SPECS))
    end
end

end # module SyntheticOCO2
