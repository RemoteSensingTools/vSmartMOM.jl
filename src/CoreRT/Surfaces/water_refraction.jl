#=
Built-in complex refractive index of liquid water, Segelstein (1981).

The table (`data/water_refractive_index_segelstein1981.txt`, next to this
file) is the complete Segelstein compilation as redistributed by the
refractiveindex.info database (CC0); its header records the source, retrieval
date, checksum, terms and the cross-check against a second redistribution.

`water_refractive_index(λ_nm)` returns `Complex(n, k)` with k ≥ 0, the
convention of `fresnel_coefficients`. Between table rows, n and k are each
interpolated linearly in wavelength. The rows are 0.23–2.6 % apart in λ
between 0.2 and 4 µm, so this differs from interpolating in log λ or log k
only at second order. Within one table interval the interpolated index moves
on a straight segment of the complex plane, which is why the spectral
Cox–Munk operator needs few knots (`_n_water_knots` in coxmunk_surface.jl).
=#

"Path of the Segelstein (1981) water refractive-index table."
const _SEGELSTEIN_1981_PATH =
    joinpath(@__DIR__, "data", "water_refractive_index_segelstein1981.txt")
include_dependency(_SEGELSTEIN_1981_PATH)

function _read_water_ri_table(path::AbstractString)
    λ_μm = Float64[]; n = Float64[]; k = Float64[]
    for line in eachline(path)
        (isempty(strip(line)) || startswith(line, '#')) && continue
        fields = split(line)
        length(fields) == 3 || error("malformed row in $path: $line")
        push!(λ_μm, parse(Float64, fields[1]))
        push!(n, parse(Float64, fields[2]))
        push!(k, parse(Float64, fields[3]))
    end
    issorted(λ_μm; lt = <=) || error("$path: wavelengths must increase strictly")
    return λ_μm .* 1000, n, k   # wavelength in nm
end

const _WATER_RI_TABLE = _read_water_ri_table(_SEGELSTEIN_1981_PATH)
const _WATER_λ_NM  = _WATER_RI_TABLE[1]
const _WATER_N_REAL = _WATER_RI_TABLE[2]
const _WATER_K_IMAG = _WATER_RI_TABLE[3]

"""
    water_refractive_index(λ_nm)

Complex refractive index `n + i k` (k ≥ 0) of liquid water at 25 °C and
vacuum wavelength `λ_nm` (nanometres), from Segelstein (1981).

Between the tabulated wavelengths, n and k are interpolated linearly in
wavelength. The table spans 33.96 nm to 10 m; a wavelength outside that range
throws an `ArgumentError`.

Returns a `Complex{Float64}`. `CoxMunkSurface(n_water = nothing)` evaluates it
at the wavelength of every spectral point.

Reference: D. J. Segelstein, "The complex refractive index of water", M.S.
thesis, University of Missouri–Kansas City (1981).
"""
function water_refractive_index(λ_nm::Real)
    λ = Float64(λ_nm)
    λ_lo, λ_hi = first(_WATER_λ_NM), last(_WATER_λ_NM)
    λ_lo <= λ <= λ_hi || throw(ArgumentError(
        "water_refractive_index: λ = $(λ_nm) nm is outside the Segelstein (1981) " *
        "table ($(λ_lo)–$(λ_hi) nm); pass `n_water` explicitly"))
    # Last row with λ_table ≤ λ, kept inside 1:N-1 so that i+1 exists.
    i = clamp(searchsortedlast(_WATER_λ_NM, λ), 1, length(_WATER_λ_NM) - 1)
    t = (λ - _WATER_λ_NM[i]) / (_WATER_λ_NM[i+1] - _WATER_λ_NM[i])
    n = _WATER_N_REAL[i] + t * (_WATER_N_REAL[i+1] - _WATER_N_REAL[i])
    k = _WATER_K_IMAG[i] + t * (_WATER_K_IMAG[i+1] - _WATER_K_IMAG[i])
    return Complex{Float64}(n, k)
end
