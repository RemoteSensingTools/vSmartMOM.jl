#=
Loaders for the solar-induced fluorescence (SIF) emission spectra and
leaf-reflectance samples bundled under `src/SIF_emission/`.
=#

using DelimitedFiles: readdlm
using DataInterpolations: LinearInterpolation

export load_sif_spectrum, load_ficus_reflectance, sif_data_path

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
    build_sif_source(RS_type, ν_model, ν_sif, jSIF; pol_component=1)

Interpolate `(ν_sif, jSIF)` onto `ν_model` (both in cm⁻¹) and write the
result into `RS_type.SIF₀[pol_component, :]`. Other Stokes components are
left untouched (default-zero for Lambertian/unpolarized SIF).

Requires `RS_type` to have a `SIF₀::Array{FT,2}` field already sized
to `(pol_type.n, length(ν_model))`. `rt_run` and `rt_run_ss` resize
`RS_type.SIF₀` before first use, so call `build_sif_source` *after*
one of them has run or after a manual resize.
"""
function build_sif_source(RS_type, ν_model::AbstractVector, ν_sif::AbstractVector,
                          jSIF::AbstractVector; pol_component::Integer = 1)
    hasproperty(RS_type, :SIF₀) ||
        throw(ArgumentError("RS_type $(typeof(RS_type)) has no SIF₀ field"))
    SIF₀ = RS_type.SIF₀
    size(SIF₀, 2) == length(ν_model) ||
        throw(DimensionMismatch("RS_type.SIF₀ has size $(size(SIF₀)) but ν_model length $(length(ν_model))"))
    FT = eltype(SIF₀)
    interp = LinearInterpolation(jSIF, ν_sif; extrapolation = ExtrapolationType.Linear)
    for (i, ν) in enumerate(ν_model)
        SIF₀[pol_component, i] = FT(interp(ν))
    end
    return SIF₀
end
