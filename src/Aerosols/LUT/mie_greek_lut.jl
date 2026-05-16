module AerosolMieGreekLUT

using Interpolations
using JLD2
using FastGaussQuadrature: gausslegendre
using LinearAlgebra: dot
using vSmartMOM.Scattering

export MieGreekLUT, build_mie_greek_lut, evaluate_lut,
       save_mie_greek_lut, load_mie_greek_lut,
       single_particle_mie_greek

"""
    MieGreekLUT

Development lookup table for homogeneous-sphere Mie quantities on
`(size_parameter, n_real, n_imag)`. Each Greek coefficient is represented by a
vector of 3-D `Interpolations.jl` interpolators, one interpolator per Legendre
order. This avoids interpolating across the discrete Legendre index while still
letting JLD2 serialize the interpolation objects, matching the existing
absorption-LUT pattern.
"""
struct MieGreekLUT{FT, IT}
    x_grid::Vector{FT}
    n_real_grid::Vector{FT}
    n_imag_grid::Vector{FT}
    l_values::Vector{Int}
    q_ext::IT
    q_sca::IT
    alpha::Vector{IT}
    beta::Vector{IT}
    gamma::Vector{IT}
    delta::Vector{IT}
    epsilon::Vector{IT}
    zeta::Vector{IT}
    metadata::Dict{String, Any}
end

function _grid_interpolator(xs, nrs, nis, values)
    return interpolate((xs, nrs, nis), values, Gridded(Linear()))
end

function _coefficient_interpolators(xs, nrs, nis, values::AbstractArray)
    return [_grid_interpolator(xs, nrs, nis, Array(view(values, il, :, :, :)))
            for il in axes(values, 1)]
end

"""
    single_particle_mie_greek(x, n_real, n_imag; lmax=nothing)

Compute Mie efficiencies and normalized Greek coefficients for one size
parameter `x = 2*pi*r/lambda` and refractive index `n_real - im*n_imag`.
Returned Greek coefficients follow the same convention as
`Scattering.compute_aerosol_optical_properties`.
"""
function single_particle_mie_greek(x::FT, n_real::FT, n_imag::FT;
                                   lmax::Union{Nothing, Integer} = nothing) where {FT<:AbstractFloat}
    x > zero(FT) || throw(ArgumentError("size parameter x must be positive"))
    n_real > zero(FT) || throw(ArgumentError("n_real must be positive"))
    n_imag >= zero(FT) || throw(ArgumentError("n_imag must be nonnegative"))

    n_mie = max(Scattering.get_n_max(x), 2)
    lmax_eff = lmax === nothing ? n_mie - 1 : Int(lmax)
    lmax_eff >= 1 || throw(ArgumentError("lmax must be at least 1"))

    m_ref = Complex{FT}(n_real, -n_imag)
    y = x * m_ref
    nmx = round(Int, max(FT(n_mie), FT(abs(y))) + FT(51))

    an = zeros(Complex{FT}, n_mie)
    bn = zeros(Complex{FT}, n_mie)
    Dn = zeros(Complex{FT}, nmx)
    Scattering.compute_mie_ab!(x, m_ref, an, bn, Dn)

    q_ext = zero(FT)
    q_sca = zero(FT)
    @inbounds for n in 1:n_mie
        weight = FT(2n + 1)
        q_ext += weight * real(an[n] + bn[n])
        q_sca += weight * (abs2(an[n]) + abs2(bn[n]))
    end
    q_ext *= FT(2) / x^2
    q_sca *= FT(2) / x^2

    l_count = lmax_eff + 1
    n_mu = max(2n_mie - 1, 2l_count + 1)
    mu, w_mu = gausslegendre(n_mu)
    mu = FT.(mu)
    w_mu = FT.(w_mu)

    leg_pi, leg_tau = Scattering.compute_mie_π_τ(mu, n_mie)
    S1 = zeros(Complex{FT}, n_mu)
    S2 = zeros(Complex{FT}, n_mu)
    Scattering.compute_mie_S₁S₂!(an, bn, leg_pi, leg_tau, S1, S2)

    f11 = zeros(FT, n_mu)
    f33 = zeros(FT, n_mu)
    f12 = zeros(FT, n_mu)
    f34 = zeros(FT, n_mu)
    inv_x2 = FT(0.5) / x^2
    @inbounds for imu in 1:n_mu
        s1 = S1[imu]
        s2 = S2[imu]
        f11[imu] = inv_x2 * real(abs2(s1) + abs2(s2))
        f33[imu] = inv_x2 * real(s1 * conj(s2) + s2 * conj(s1))
        f12[imu] = -inv_x2 * real(abs2(s1) - abs2(s2))
        f34[imu] = -inv_x2 * imag(s1 * conj(s2) - s2 * conj(s1))
    end

    # Match compute_NAI2.jl normalization for a monodisperse particle.
    if q_sca > zero(FT)
        f_scale = FT(4) / q_sca
        f11 .*= f_scale
        f33 .*= f_scale
        f12 .*= f_scale
        f34 .*= f_scale
    else
        fill!(f11, zero(FT))
        fill!(f33, zero(FT))
        fill!(f12, zero(FT))
        fill!(f34, zero(FT))
    end

    P, P2, R2, T2 = Scattering.compute_legendre_poly(mu, l_count)
    alpha = zeros(FT, l_count)
    beta = zeros(FT, l_count)
    gamma = zeros(FT, l_count)
    delta = zeros(FT, l_count)
    epsilon = zeros(FT, l_count)
    zeta = zeros(FT, l_count)

    @inbounds for l in 0:lmax_eff
        il = l + 1
        fac = l >= 2 ? FT(2l + 1) / FT(2) *
              sqrt(one(FT) / (FT(l - 1) * FT(l) * FT(l + 1) * FT(l + 2))) :
              zero(FT)
        beta[il] = FT(2l + 1) / FT(2) * dot(w_mu, f11 .* P[:, il])
        delta[il] = FT(2l + 1) / FT(2) * dot(w_mu, f33 .* P[:, il])
        gamma[il] = fac * dot(w_mu, f12 .* P2[:, il])
        epsilon[il] = fac * dot(w_mu, f34 .* P2[:, il])
        zeta[il] = fac * dot(w_mu, f33 .* R2[:, il] .+ f11 .* T2[:, il])
        alpha[il] = fac * dot(w_mu, f11 .* R2[:, il] .+ f33 .* T2[:, il])
    end

    greek = Scattering.GreekCoefs(alpha, beta, gamma, delta, epsilon, zeta)
    return (; q_ext = max(zero(FT), q_ext),
             q_sca = max(zero(FT), q_sca),
             greek_coefs = greek)
end

"""
    build_mie_greek_lut(x_grid, n_real_grid, n_imag_grid; lmax, FT=Float64)

Build a JLD2-serializable LUT for Mie efficiencies and Greek coefficients.
The axes are size parameter and complex refractive index; wavelength dependence
enters later through `x = 2*pi*r/lambda` and the refractive-index database.
"""
function build_mie_greek_lut(x_grid, n_real_grid, n_imag_grid;
                             lmax::Integer,
                             FT::Type{<:AbstractFloat} = Float64,
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

    lmax_i = Int(lmax)
    l_values = collect(0:lmax_i)
    shape = (length(xs), length(nrs), length(nis))
    coeff_shape = (length(l_values), shape...)

    q_ext = zeros(FT, shape)
    q_sca = zeros(FT, shape)
    alpha = zeros(FT, coeff_shape)
    beta = zeros(FT, coeff_shape)
    gamma = zeros(FT, coeff_shape)
    delta = zeros(FT, coeff_shape)
    epsilon = zeros(FT, coeff_shape)
    zeta = zeros(FT, coeff_shape)

    @inbounds for ix in eachindex(xs), ir in eachindex(nrs), ii in eachindex(nis)
        result = single_particle_mie_greek(xs[ix], nrs[ir], nis[ii]; lmax = lmax_i)
        q_ext[ix, ir, ii] = result.q_ext
        q_sca[ix, ir, ii] = result.q_sca
        greek = result.greek_coefs
        alpha[:, ix, ir, ii] .= greek.α
        beta[:, ix, ir, ii] .= greek.β
        gamma[:, ix, ir, ii] .= greek.γ
        delta[:, ix, ir, ii] .= greek.δ
        epsilon[:, ix, ir, ii] .= greek.ϵ
        zeta[:, ix, ir, ii] .= greek.ζ
    end

    md = Dict{String, Any}(metadata)
    md["schema"] = "vSmartMOM AerosolMieGreekLUT"
    md["axes"] = ["size_parameter", "n_real", "n_imag"]
    md["lmax"] = lmax_i

    q_ext_itp = _grid_interpolator(xs, nrs, nis, q_ext)
    q_sca_itp = _grid_interpolator(xs, nrs, nis, q_sca)
    alpha_itp = _coefficient_interpolators(xs, nrs, nis, alpha)
    beta_itp = _coefficient_interpolators(xs, nrs, nis, beta)
    gamma_itp = _coefficient_interpolators(xs, nrs, nis, gamma)
    delta_itp = _coefficient_interpolators(xs, nrs, nis, delta)
    epsilon_itp = _coefficient_interpolators(xs, nrs, nis, epsilon)
    zeta_itp = _coefficient_interpolators(xs, nrs, nis, zeta)

    return MieGreekLUT(xs, nrs, nis, l_values, q_ext_itp, q_sca_itp,
                       alpha_itp, beta_itp, gamma_itp, delta_itp,
                       epsilon_itp, zeta_itp, md)
end

function _eval_coeffs(interpolators, x, n_real, n_imag, nterms::Int)
    return [interpolators[i](x, n_real, n_imag) for i in 1:nterms]
end

"""
    evaluate_lut(lut, x, n_real, n_imag; lmax=maximum(lut.l_values))

Interpolate Mie efficiencies and Greek coefficients from a LUT.
"""
function evaluate_lut(lut::MieGreekLUT, x::Real, n_real::Real, n_imag::Real;
                      lmax::Integer = maximum(lut.l_values))
    nterms = Int(lmax) + 1
    nterms <= length(lut.l_values) ||
        throw(ArgumentError("requested lmax=$(lmax), LUT only contains $(maximum(lut.l_values))"))

    q_ext = lut.q_ext(x, n_real, n_imag)
    q_sca = lut.q_sca(x, n_real, n_imag)
    greek = Scattering.GreekCoefs(
        _eval_coeffs(lut.alpha, x, n_real, n_imag, nterms),
        _eval_coeffs(lut.beta, x, n_real, n_imag, nterms),
        _eval_coeffs(lut.gamma, x, n_real, n_imag, nterms),
        _eval_coeffs(lut.delta, x, n_real, n_imag, nterms),
        _eval_coeffs(lut.epsilon, x, n_real, n_imag, nterms),
        _eval_coeffs(lut.zeta, x, n_real, n_imag, nterms),
    )
    return (; q_ext, q_sca, greek_coefs = greek)
end

function save_mie_greek_lut(filepath::AbstractString, lut::MieGreekLUT)
    JLD2.jldsave(filepath; lut)
    return filepath
end

function load_mie_greek_lut(filepath::AbstractString)
    return JLD2.load(filepath, "lut")
end

end # module
