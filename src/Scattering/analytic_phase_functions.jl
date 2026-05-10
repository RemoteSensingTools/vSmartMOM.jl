_phase_float_type_of(x::Real) = typeof(float(x))
_phase_float_type_of(::Type{FT}) where {FT<:Real} = float(FT)
_phase_float_type(xs...) = promote_type(map(_phase_float_type_of, xs)...)

function _validate_hg_asymmetry(g)
    abs(g) < one(g) ||
        throw(ArgumentError("Henyey-Greenstein asymmetry parameter must satisfy abs(g) < 1"))
    return nothing
end

function _validate_synthetic_polarization(p)
    abs(p) <= one(p) ||
        throw(ArgumentError("synthetic polarization_fraction must satisfy abs(p) <= 1"))
    return nothing
end

function phase_function(phase::HenyeyGreensteinPhaseFunction, cosΘ::Real)
    FT = _phase_float_type(phase.g, cosΘ)
    g = convert(FT, phase.g)
    _validate_hg_asymmetry(g)
    x = convert(FT, cosΘ)
    return (one(FT) - g^2) / (one(FT) + g^2 - FT(2) * g * x)^(FT(1.5))
end

phase_function(phase::HenyeyGreensteinPhaseFunction,
               cosΘ::AbstractArray{<:Real}) =
    phase_function.(Ref(phase), cosΘ)

phase_function(phase::SyntheticPolarizedHenyeyGreensteinPhaseFunction,
               cosΘ::Real) =
    phase_function(HenyeyGreensteinPhaseFunction(phase.g), cosΘ)

phase_function(phase::SyntheticPolarizedHenyeyGreensteinPhaseFunction,
               cosΘ::AbstractArray{<:Real}) =
    phase_function.(Ref(phase), cosΘ)

function scattering_matrix(phase::HenyeyGreensteinPhaseFunction,
                           μ::AbstractVector{<:Real})
    f11 = phase_function(phase, μ)
    z = zero.(f11)
    return ScatteringMatrix(f11, z, copy(f11), copy(f11), z, copy(f11))
end

function scattering_matrix(phase::SyntheticPolarizedHenyeyGreensteinPhaseFunction,
                           μ::AbstractVector{<:Real})
    FT = _phase_float_type(phase.g, phase.polarization_fraction, eltype(μ))
    p = convert(FT, phase.polarization_fraction)
    _validate_synthetic_polarization(p)

    μ_ft = FT.(μ)
    f11 = FT.(phase_function(phase, μ_ft))
    f12 = similar(f11)
    @inbounds for i in eachindex(μ_ft)
        x² = μ_ft[i]^2
        f12[i] = p * f11[i] * (one(FT) - x²) / (one(FT) + x²)
    end

    z = zero.(f11)
    return ScatteringMatrix(f11, f12, copy(f11), copy(f11), z, copy(f11))
end

"""
    greek_coefficients_from_scattering_matrix(μ, w, matrix)

Convert angle-space scattering-matrix elements sampled at quadrature nodes
`μ` with weights `w` into Sanghavi Greek coefficients. This is the same
coefficient projection used by the Mie/NAI2 path, factored for analytic phase
functions.
"""
function greek_coefficients_from_scattering_matrix(
    μ::AbstractVector{<:Real}, w::AbstractVector{<:Real},
    matrix::ScatteringMatrix; l_max::Int = length(μ))

    length(μ) == length(w) ||
        throw(ArgumentError("μ and w must have the same length"))
    l_max > 0 || throw(ArgumentError("l_max must be positive"))
    l_max <= length(μ) ||
        throw(ArgumentError("l_max must be no larger than the number of quadrature nodes"))
    P, P², R², T² = compute_legendre_poly(μ, l_max)
    FT = promote_type(eltype(μ), eltype(w), eltype(matrix.f₁₁),
                      eltype(matrix.f₁₂), eltype(matrix.f₃₃),
                      eltype(matrix.f₃₄))

    α = zeros(FT, l_max)
    β = zeros(FT, l_max)
    γ = zeros(FT, l_max)
    δ = zeros(FT, l_max)
    ϵ = zeros(FT, l_max)
    ζ = zeros(FT, l_max)

    f11 = FT.(matrix.f₁₁)
    f12 = FT.(matrix.f₁₂)
    f33 = FT.(matrix.f₃₃)
    f34 = FT.(matrix.f₃₄)
    w_ft = FT.(w)

    @inbounds for l in 0:(l_max - 1)
        half_factor = FT(2l + 1) / FT(2)
        fac = l >= 2 ?
              half_factor * sqrt(inv(FT((l - 1) * l * (l + 1) * (l + 2)))) :
              zero(FT)
        il = l + 1
        β[il] = half_factor * dot(w_ft, f11 .* view(P, :, il))
        δ[il] = half_factor * dot(w_ft, f33 .* view(P, :, il))
        γ[il] = fac * dot(w_ft, f12 .* view(P², :, il))
        ϵ[il] = fac * dot(w_ft, f34 .* view(P², :, il))
        ζ[il] = fac * dot(w_ft, f33 .* view(R², :, il) .+
                                f11 .* view(T², :, il))
        α[il] = fac * dot(w_ft, f11 .* view(R², :, il) .+
                                f33 .* view(T², :, il))
    end

    return GreekCoefs(α, β, γ, δ, ϵ, ζ)
end

"""
    greek_coefficients(phase; l_max=64, nquad=max(2l_max+1, 64))

Project an analytic phase/scattering matrix into Greek coefficients so it can
enter the same CoreRT MOM path as Mie-derived aerosol optics.
"""
function greek_coefficients(phase::AbstractAnalyticPhaseFunction;
                            l_max::Int = 64,
                            nquad::Int = max(2l_max + 1, 64))
    l_max > 0 || throw(ArgumentError("l_max must be positive"))
    nquad >= l_max ||
        throw(ArgumentError("nquad must be at least l_max"))
    μ64, w64 = gausslegendre(nquad)
    FT = phase isa HenyeyGreensteinPhaseFunction ?
         _phase_float_type(phase.g) :
         _phase_float_type(phase.g, phase.polarization_fraction)
    μ = FT.(μ64)
    w = FT.(w64)
    return greek_coefficients_from_scattering_matrix(
        μ, w, scattering_matrix(phase, μ); l_max)
end

"""
    analytic_aerosol_optics(phase; single_scattering_albedo=1, extinction_cross_section=1, l_max=64)

Create an `AerosolOptics` object from an analytic phase/scattering matrix.
The returned object is directly consumable by the existing CoreRT aerosol
mixing and Fourier-moment code.
"""
function analytic_aerosol_optics(
    phase::AbstractAnalyticPhaseFunction;
    single_scattering_albedo::Real = 1,
    extinction_cross_section::Real = 1,
    l_max::Int = 64,
    nquad::Int = max(2l_max + 1, 64))

    greek = greek_coefficients(phase; l_max, nquad)
    FT = promote_type(eltype(greek.β), typeof(float(single_scattering_albedo)),
                      typeof(float(extinction_cross_section)))
    return AerosolOptics(greek_coefs=greek,
                         ω̃=convert(FT, single_scattering_albedo),
                         k=convert(FT, extinction_cross_section),
                         fᵗ=zero(FT))
end

compute_aerosol_optical_properties(
    phase::AbstractAnalyticPhaseFunction; kwargs...) =
    analytic_aerosol_optics(phase; kwargs...)

_scattering_angle_cosine(μ₀, μv, Δϕ) =
    -μ₀ * μv + sqrt(max(zero(μ₀), one(μ₀) - μ₀^2)) *
               sqrt(max(zero(μv), one(μv) - μv^2)) * cos(Δϕ)

function _rotation_from_scattering_plane(μ₀::FT, μv::FT, Δϕ::FT,
                                         cosΘ::FT) where {FT}
    sinΘ² = max(zero(FT), one(FT) - cosΘ^2)
    if sinΘ² <= eps(FT)
        return one(FT), zero(FT)
    end

    sinΘ = sqrt(sinΘ²)
    sin₀ = sqrt(max(zero(FT), one(FT) - μ₀^2))
    sinv = sqrt(max(zero(FT), one(FT) - μv^2))
    cosχ = (μ₀ * sinv + μv * sin₀ * cos(Δϕ)) / sinΘ
    sinχ = sin₀ * sin(Δϕ) / sinΘ
    return cosχ^2 - sinχ^2, FT(2) * sinχ * cosχ
end

"""
    phase_matrix_first_column(greek, μ₀, μv, Δϕ, Val(N))

Evaluate the first column of the scattering phase matrix for an unpolarized
incident beam at one plane-parallel sun-view geometry. `μ₀` and `μv` are
positive cosines and `Δϕ` is the vSmartMOM relative azimuth in radians.
"""
function phase_matrix_first_column(greek::GreekCoefs, μ₀::Real, μv::Real,
                                   Δϕ::Real, ::Val{N}) where {N}
    FT = promote_type(eltype(greek.β), typeof(float(μ₀)),
                      typeof(float(μv)), typeof(float(Δϕ)))
    μ₀_ft = convert(FT, μ₀)
    μv_ft = convert(FT, μv)
    Δϕ_ft = convert(FT, Δϕ)
    cosΘ = _scattering_angle_cosine(μ₀_ft, μv_ft, Δϕ_ft)
    matrix = reconstruct_phase(greek, [cosΘ])
    f11 = convert(FT, matrix.f₁₁[1])
    f12 = convert(FT, matrix.f₁₂[1])
    cos2χ, sin2χ = _rotation_from_scattering_plane(μ₀_ft, μv_ft, Δϕ_ft, cosΘ)
    q = f12 * cos2χ
    u = f12 * sin2χ
    return ntuple(i -> i == 1 ? f11 :
                       i == 2 ? q :
                       i == 3 ? u : zero(FT), Val(N))
end
