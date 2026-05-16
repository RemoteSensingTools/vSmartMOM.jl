module AerosolMieExtinctionLUT

using Interpolations
using JLD2
using vSmartMOM.Scattering

export MieExtinctionLUT, build_mie_extinction_lut, evaluate_extinction_lut,
       save_mie_extinction_lut, load_mie_extinction_lut, single_particle_qext

"""
    MieExtinctionLUT

Small diagnostic lookup table for aerosol AOD maps. The table stores
dimensionless Mie extinction efficiency `Q_ext(x, n_real, n_imag)` on a
3-D gridded interpolator. Runtime code converts to cross-section via
`C_ext = πr² Q_ext`.

This LUT is intentionally separate from the full Greek-coefficient LUT:
it is the fast path for extinction/AOD diagnostics only.
"""
struct MieExtinctionLUT{FT, IT}
    x_grid::Vector{FT}
    n_real_grid::Vector{FT}
    n_imag_grid::Vector{FT}
    q_ext::IT
    metadata::Dict{String, Any}
end

function single_particle_qext(x::FT, n_real::FT, n_imag::FT) where {FT<:AbstractFloat}
    x > zero(FT) || throw(ArgumentError("size parameter x must be positive"))
    n_real > zero(FT) || throw(ArgumentError("n_real must be positive"))
    n_imag >= zero(FT) || throw(ArgumentError("n_imag must be nonnegative"))

    n_mie = Scattering.get_n_max(x)
    m_ref = Complex{FT}(n_real, -n_imag)
    y = x * m_ref
    nmx = round(Int, max(FT(n_mie), FT(abs(y))) + FT(51))

    an = zeros(Complex{FT}, n_mie)
    bn = zeros(Complex{FT}, n_mie)
    Dn = zeros(Complex{FT}, nmx)
    Scattering.compute_mie_ab!(x, m_ref, an, bn, Dn)

    series = zero(FT)
    @inbounds for n in 1:n_mie
        series += FT(2n + 1) * real(an[n] + bn[n])
    end
    return max(zero(FT), FT(2) * series / x^2)
end

function build_mie_extinction_lut(x_grid, n_real_grid, n_imag_grid;
                                  FT::Type{<:AbstractFloat} = Float32,
                                  metadata = Dict{String, Any}())
    xs = collect(FT, x_grid)
    nrs = collect(FT, n_real_grid)
    nis = collect(FT, n_imag_grid)
    issorted(xs) || throw(ArgumentError("x_grid must be sorted ascending"))
    issorted(nrs) || throw(ArgumentError("n_real_grid must be sorted ascending"))
    issorted(nis) || throw(ArgumentError("n_imag_grid must be sorted ascending"))
    all(>(zero(FT)), xs) || throw(ArgumentError("x_grid values must be positive"))
    all(>(zero(FT)), nrs) || throw(ArgumentError("n_real_grid values must be positive"))
    all(>=(zero(FT)), nis) || throw(ArgumentError("n_imag_grid values must be nonnegative"))

    q_ext = zeros(FT, length(xs), length(nrs), length(nis))
    @inbounds for ix in eachindex(xs), ir in eachindex(nrs), ii in eachindex(nis)
        q_ext[ix, ir, ii] = single_particle_qext(xs[ix], nrs[ir], nis[ii])
    end

    md = Dict{String, Any}(metadata)
    md["schema"] = "vSmartMOM AerosolMieExtinctionLUT"
    md["axes"] = ["size_parameter", "n_real", "n_imag"]
    md["quantity"] = "Q_ext"
    md["note"] = "AOD diagnostic only; does not contain scattering phase information."

    itp = interpolate((xs, nrs, nis), q_ext, Gridded(Linear()))
    return MieExtinctionLUT(xs, nrs, nis, itp, md)
end

evaluate_extinction_lut(lut::MieExtinctionLUT, x::Real, n_real::Real, n_imag::Real) =
    max(zero(eltype(lut.x_grid)), lut.q_ext(x, n_real, n_imag))

function save_mie_extinction_lut(filepath::AbstractString, lut::MieExtinctionLUT)
    JLD2.jldsave(filepath; lut)
    return filepath
end

function load_mie_extinction_lut(filepath::AbstractString)
    return JLD2.load(filepath, "lut")
end

end # module
