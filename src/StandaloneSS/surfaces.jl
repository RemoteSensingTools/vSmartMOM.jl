_real_numeric_type(x::Real) = typeof(x)
_real_numeric_type(x::Complex) = typeof(real(x))
_real_numeric_type(x::AbstractArray) = _real_numeric_type(eltype(x))
_real_numeric_type(::Type{Complex{FT}}) where {FT} = FT
_real_numeric_type(::Type{FT}) where {FT<:Real} = FT

_surface_numeric_type(surface::LambertianSSSurface) =
    _real_numeric_type(surface.albedo)

_surface_numeric_type(::CoxMunkSSSurface{W, Nothing, A}) where {W, A} =
    promote_type(W, A)

_surface_numeric_type(surface::CoxMunkSSSurface{W, N, A}) where {W, N, A} =
    promote_type(W, A, _real_numeric_type(surface.n_water))

_contributor_numeric_type(c::RayleighSSContributor) =
    promote_type(_real_numeric_type(c.τ), typeof(c.depol))
_contributor_numeric_type(c::AbsorptionSSContributor) =
    _real_numeric_type(c.τ)
_contributor_numeric_type(c::HGAerosolSSContributor) =
    promote_type(_real_numeric_type(c.τ), typeof(c.g), typeof(c.ϖ))
_contributor_numeric_type(c::GreekCoefsSSContributor) =
    promote_type(_real_numeric_type(c.τ), typeof(c.ϖ),
                 _real_numeric_type(c.greek_coefs.β))

_contributors_numeric_type(::Tuple{}) = Bool
_contributors_numeric_type(contributors::Tuple) =
    promote_type(_contributor_numeric_type(first(contributors)),
                 _contributors_numeric_type(Base.tail(contributors)))
_contributors_numeric_type(contributors) =
    mapreduce(_contributor_numeric_type, promote_type, contributors; init=Bool)

function _config_numeric_type(config::ExactSSConfig)
    return promote_type(typeof(config.geometry.μ₀),
                        eltype(config.geometry.μv),
                        eltype(config.geometry.Δϕ),
                        _real_numeric_type(config.I0),
                        _surface_numeric_type(config.surface),
                        _contributors_numeric_type(config.contributors))
end

_precompute_surface_brdf(surface::AbstractSSSurface, geometry::SSGeometry,
                         n_spec::Int, ::Type{FT}) where {FT} =
    _precompute_surface_brdf(Val(1), surface, geometry, n_spec, FT)

function _precompute_surface_brdf(::Val{1}, surface::LambertianSSSurface,
                                  geometry::SSGeometry, n_spec::Int,
                                  ::Type{FT}) where {FT}
    albedo = _vectorize_albedo(surface, n_spec, FT)
    n_geom = length(geometry.μv)
    ρ = zeros(FT, n_geom, n_spec)
    @inbounds for ispec in 1:n_spec, iv in 1:n_geom
        ρ[iv, ispec] = albedo[ispec] / FT(pi)
    end
    return ρ
end

function _precompute_surface_brdf(::Val{N}, surface::LambertianSSSurface,
                                  geometry::SSGeometry, n_spec::Int,
                                  ::Type{FT}) where {N,FT}
    albedo = _vectorize_albedo(surface, n_spec, FT)
    n_geom = length(geometry.μv)
    ρ = zeros(FT, n_geom, N, n_spec)
    @inbounds for ispec in 1:n_spec, iv in 1:n_geom
        ρ[iv, 1, ispec] = albedo[ispec] / FT(pi)
    end
    return ρ
end

function _precompute_surface_brdf(::Val{1}, surface::CoxMunkSSSurface,
                                  geometry::SSGeometry, n_spec::Int,
                                  ::Type{FT}) where {FT}
    n_geom = length(geometry.μv)
    ρ = zeros(FT, n_geom, n_spec)
    μ₀ = convert(FT, geometry.μ₀)
    @inbounds for ispec in 1:n_spec, iv in 1:n_geom
        μᵥ = convert(FT, geometry.μv[iv])
        Δϕ = convert(FT, geometry.Δϕ[iv])
        ρ[iv, ispec] = _coxmunk_brdf_scalar(surface, μᵥ, μ₀, Δϕ,
                                            ispec, n_spec, FT)
    end
    return ρ
end

function _precompute_surface_brdf(::Val{N}, surface::CoxMunkSSSurface,
                                  geometry::SSGeometry, n_spec::Int,
                                  ::Type{FT}) where {N,FT}
    n_geom = length(geometry.μv)
    ρ = zeros(FT, n_geom, N, n_spec)
    μ₀ = convert(FT, geometry.μ₀)
    core_surface = CoxMunkSurface{FT}(
        wind_speed=convert(FT, surface.wind_speed),
        n_water=nothing,
        whitecap_albedo=convert(FT, surface.whitecap_albedo),
        include_whitecaps=surface.include_whitecaps,
        shadowing=surface.shadowing)
    @inbounds for ispec in 1:n_spec, iv in 1:n_geom
        μᵥ = convert(FT, geometry.μv[iv])
        Δϕ = convert(FT, geometry.Δϕ[iv])
        n_water = _coxmunk_n_water(surface, ispec, n_spec, FT)
        M = coxmunk_brdf_mueller(core_surface, N, μᵥ, μ₀, Δϕ;
                                  n_water)
        for istokes in 1:N
            ρ[iv, istokes, ispec] = M[istokes, 1]
        end
    end
    return ρ
end

function _coxmunk_n_water(surface::CoxMunkSSSurface, ispec::Int, n_spec::Int,
                          ::Type{FT}) where {FT}
    n_water = surface.n_water
    n_water === nothing && return Complex{FT}(FT(1.34), zero(FT))
    if n_water isa Complex
        return Complex{FT}(convert(FT, real(n_water)),
                           convert(FT, imag(n_water)))
    end
    length(n_water) == n_spec ||
        throw(ArgumentError("CoxMunkSSSurface.n_water vector length must match nSpec"))
    n = n_water[ispec]
    return Complex{FT}(convert(FT, real(n)), convert(FT, imag(n)))
end

@inline _wind_to_sigma2(U::FT) where {FT} = FT(0.003) + FT(0.00512) * U

@inline function _coxmunk_pdf(zx::FT, zy::FT, σ²::FT) where {FT}
    return exp(-(zx^2 + zy^2) / (FT(2) * σ²)) / (FT(2π) * σ²)
end

@inline function _whitecap_fraction(U::FT) where {FT}
    U <= zero(FT) && return zero(FT)
    return FT(2.95e-6) * U^FT(3.52)
end

function _smith_Λ(μ::FT, σ²::FT) where {FT}
    μ <= zero(FT) && return FT(1e10)
    σ = sqrt(σ²)
    cotθ = μ / sqrt(max(FT(1e-30), one(FT) - μ^2))
    ν = cotθ / (sqrt(FT(2)) * σ)
    Λ = (exp(-ν^2) / (sqrt(FT(2π)) * ν) - erfc(ν)) / FT(2)
    return max(zero(FT), Λ)
end

@inline function _shadow_factor(μᵢ::FT, μᵣ::FT, σ²::FT) where {FT}
    return one(FT) / (one(FT) + _smith_Λ(μᵢ, σ²) + _smith_Λ(μᵣ, σ²))
end

function _coxmunk_geometry(μᵢ::FT, μᵣ::FT, Δϕ::FT) where {FT}
    sinᵢ = sqrt(max(zero(FT), one(FT) - μᵢ^2))
    sinᵣ = sqrt(max(zero(FT), one(FT) - μᵣ^2))
    cosΔ = cos(Δϕ)
    sinΔ = sin(Δϕ)

    nx = -sinᵢ + sinᵣ * cosΔ
    ny = sinᵣ * sinΔ
    nz = μᵢ + μᵣ
    norm_n = sqrt(nx^2 + ny^2 + nz^2)
    if norm_n < FT(1e-15)
        return (; cosβ=one(FT), cosθ_local=one(FT),
                zx=zero(FT), zy=zero(FT))
    end

    nx /= norm_n
    ny /= norm_n
    nz /= norm_n

    cosβ = max(FT(1e-10), nz)
    cosθ_local = clamp((μᵢ + μᵣ) / (FT(2) * cosβ), zero(FT), one(FT))
    zx = -nx / cosβ
    zy = -ny / cosβ
    return (; cosβ, cosθ_local, zx, zy)
end

function _fresnel_reflectance(n_rel::Complex{FT}, cosθ::FT) where {FT}
    sinθ² = max(zero(FT), one(FT) - cosθ^2)
    cosθt = sqrt(one(Complex{FT}) - sinθ² / n_rel^2)
    rs = (cosθ - n_rel * cosθt) / (cosθ + n_rel * cosθt)
    rp = (n_rel * cosθ - cosθt) / (n_rel * cosθ + cosθt)
    return (abs2(rs) + abs2(rp)) / FT(2)
end

function _coxmunk_brdf_scalar(surface::CoxMunkSSSurface, μᵢ::FT, μᵣ::FT,
                              Δϕ::FT, ispec::Int, n_spec::Int,
                              ::Type{FT}) where {FT}
    σ² = _wind_to_sigma2(convert(FT, surface.wind_speed))
    geom = _coxmunk_geometry(μᵢ, μᵣ, Δϕ)
    (; cosβ, cosθ_local, zx, zy) = geom

    P = _coxmunk_pdf(zx, zy, σ²)
    n_water = _coxmunk_n_water(surface, ispec, n_spec, FT)
    fresnel = _fresnel_reflectance(n_water, cosθ_local)
    prefactor = P / (FT(4) * μᵢ * μᵣ * cosβ^4)
    if surface.shadowing
        prefactor *= _shadow_factor(μᵢ, μᵣ, σ²)
    end

    ρ_glint = prefactor * fresnel
    if surface.include_whitecaps
        f_wc = _whitecap_fraction(convert(FT, surface.wind_speed))
        ρ_wc = convert(FT, surface.whitecap_albedo) / FT(pi)
        return (one(FT) - f_wc) * ρ_glint + f_wc * ρ_wc
    end
    return ρ_glint
end
