#=
Loaders for the solar-induced fluorescence (SIF) emission spectra and
leaf-reflectance samples bundled under `src/SIF_emission/`.
=#

using DelimitedFiles: readdlm
using DataInterpolations: LinearInterpolation

export load_sif_spectrum, sif_reference_state, load_ficus_reflectance,
       sif_data_path, build_sif_source

"""
    sif_data_path(filename) -> String

Absolute path to `filename` inside `src/SIF_emission/`. Benchmark scripts
should use this helper rather than hardcoded absolute paths so that the
package remains relocatable via `pkgdir(vSmartMOM)`.
"""
sif_data_path(filename::AbstractString) = joinpath(@__DIR__, filename)

"""
    load_sif_spectrum(path=sif_data_path("sif-spectra.csv"); column=:SIF_OLD,
                      rescale_to_peak=true) -> (ν::Vector, jSIF::Vector)

Load a SIF emission spectrum from a CSV file and return
`(ν, jSIF)` sorted by increasing wavenumber (cm⁻¹).

- `ν` is the wavenumber grid in cm⁻¹ (converted from the file's nm column).
- `jSIF` is the SIF flux in mW/m²/cm⁻¹ (converted from mW/m²/nm via the
  dλ/dν = 1e7/ν² Jacobian).

`column` selects which SIF spectrum to read — the `sif-spectra.csv`
file ships three: `:SIF_OLD`, `:SIF_NEW`, `:SIF_DEF`. Defaults to
`:SIF_OLD` to match the benchmark scripts.

TODO: The (0.5π / maximum(J_SIF)) rescaling is an intentional hack to
make SIF magnitude data-independent for grid generation. This normalizes
shape but discards absolute physical magnitude. Revisit: confirm the
downstream physics depends only on SIF shape (not absolute flux), or
replace with physical units (mW/m²/cm⁻¹). Not a merge blocker.
"""
function load_sif_spectrum(path::AbstractString = sif_data_path("sif-spectra.csv");
                           column::Symbol = :SIF_OLD,
                           rescale_to_peak::Bool = true)
    raw = readdlm(path, ',')
    header = Symbol.(raw[1, :])
    col_idx = findfirst(==(column), header)
    col_idx === nothing && throw(ArgumentError(
        "column $column not found in $path (available: $(header[2:end]))"))

    wl_nm = Float64.(raw[2:end, 1])
    flux_nm = Float64.(raw[2:end, col_idx])

    if rescale_to_peak
        peak = maximum(flux_nm)
        peak > 0 || throw(ArgumentError("SIF spectrum peak is zero; cannot rescale"))
        flux_nm .*= (0.5π / peak)
    end

    ν = reverse(1e7 ./ wl_nm)
    jSIF = reverse(flux_nm)
    jSIF .*= 1e7 ./ (ν .^ 2)

    return ν, jSIF
end

"""
    sif_reference_state(; total_sif=0.5, reference_wavelength_nm=760,
                        column=:SIF_OLD, path=sif_data_path("sif-spectra.csv"))

Normalize the supplied SIF shape so that its wavelength-space integral is
`total_sif` (mW m⁻² sr⁻¹), convert it to spectral radiance per cm⁻¹, and
return the local wavenumber-linear state at `reference_wavelength_nm`:

`(SIF760=SIF(ν₇₆₀), mSIF=dSIF/dν|₇₆₀, ν_ref=ν₇₆₀, ν, spectrum)`.

The `SIF760` property retains that name for the retrieval convention even if
a different diagnostic reference wavelength is requested. The slope is a
centered secant through the nearest tabulated samples bracketing the reference.
"""
function sif_reference_state(; total_sif::Real=0.5,
                             reference_wavelength_nm::Real=760,
                             column::Symbol=:SIF_OLD,
                             path::AbstractString=sif_data_path("sif-spectra.csv"))
    total_sif >= 0 || throw(ArgumentError("total_sif must be nonnegative"))
    ν, shape_ν = load_sif_spectrum(path; column, rescale_to_peak=false)
    area = sum(((@view shape_ν[1:end-1]) .+ (@view shape_ν[2:end])) .* diff(ν)) / 2
    area > 0 || throw(ArgumentError("SIF spectrum has nonpositive integrated area"))
    spectrum = shape_ν .* (total_sif / area)
    ν_ref = 1e7 / reference_wavelength_nm
    ihi = searchsortedfirst(ν, ν_ref)
    1 < ihi <= length(ν) || throw(ArgumentError(
        "reference wavelength $reference_wavelength_nm nm lies outside the SIF spectrum"))
    ilo = ihi - 1
    SIF760 = spectrum[ilo] + (spectrum[ihi] - spectrum[ilo]) *
             (ν_ref - ν[ilo]) / (ν[ihi] - ν[ilo])
    # Use one sample on either side of an exact tabulated reference so the
    # derivative is centered rather than selecting one piecewise-linear side.
    if isapprox(ν[ihi], ν_ref; rtol=0, atol=8eps(ν_ref)) && ihi < length(ν)
        ilo = ihi - 1
        ihi += 1
    end
    mSIF = (spectrum[ihi] - spectrum[ilo]) / (ν[ihi] - ν[ilo])
    return (; SIF760, mSIF, ν_ref, ν, spectrum)
end

"""
    load_ficus_reflectance(path=sif_data_path("ficus_refl_600to800nm.dat")) ->
        (λ_μm::Vector, R::Vector)

Load a leaf-reflectance sample from a simple two-column file (λ in μm,
reflectance in %). The 600-800 nm file has no header; the full-range
`ficus_refl.dat` file has a metadata header and the numeric block starts
after the first blank line. Both are handled. Returned reflectance is
dimensionless (the % is divided out).
"""
function load_ficus_reflectance(path::AbstractString = sif_data_path("ficus_refl_600to800nm.dat"))
    lines = readlines(path)
    start = something(findfirst(isempty ∘ strip, lines), 0) + 1
    data = readdlm(IOBuffer(join(lines[start:end], '\n')))
    λ_μm = Float64.(data[:, 1])
    R = Float64.(data[:, 2]) ./ 100
    return λ_μm, R
end

"""
    build_sif_source(SIF₀, ν_model, ν_sif, jSIF; pol_component=1)
    build_sif_source(RS_type, ν_model, ν_sif, jSIF; pol_component=1)

Interpolate `(ν_sif, jSIF)` onto `ν_model` (both in cm⁻¹) and write the
result into `SIF₀[pol_component, :]`. Other Stokes components are left
untouched (default-zero for Lambertian/unpolarized SIF).

The matrix method is the preferred path: pass the returned matrix into
`SurfaceSIF(SIF₀=SIF₀)`. The `RS_type` method remains only for older scripts
that still carry a legacy `SIF₀` field.
"""
function build_sif_source(SIF₀::AbstractMatrix, ν_model::AbstractVector,
                          ν_sif::AbstractVector, jSIF::AbstractVector;
                          pol_component::Integer = 1)
    1 ≤ pol_component ≤ size(SIF₀, 1) ||
        throw(ArgumentError("pol_component=$pol_component is outside SIF₀ rows 1:$(size(SIF₀, 1))"))
    size(SIF₀, 2) == length(ν_model) ||
        throw(DimensionMismatch("SIF₀ has size $(size(SIF₀)) but ν_model length $(length(ν_model))"))

    FT = eltype(SIF₀)
    interp = LinearInterpolation(jSIF, ν_sif; extrapolation = ExtrapolationType.Linear)
    @inbounds for (i, ν) in enumerate(ν_model)
        SIF₀[pol_component, i] = FT(interp(ν))
    end
    return SIF₀
end

function build_sif_source(RS_type, ν_model::AbstractVector, ν_sif::AbstractVector,
                          jSIF::AbstractVector; pol_component::Integer = 1)
    hasproperty(RS_type, :SIF₀) ||
        throw(ArgumentError("RS_type $(typeof(RS_type)) has no SIF₀ field"))
    return build_sif_source(RS_type.SIF₀, ν_model, ν_sif, jSIF; pol_component)
end
