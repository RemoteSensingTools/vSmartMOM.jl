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
    core_surface = _coxmunk_core_surface(surface, FT)
    @inbounds for ispec in 1:n_spec, iv in 1:n_geom
        μᵥ = convert(FT, geometry.μv[iv])
        Δϕ = convert(FT, geometry.Δϕ[iv])
        n_water = _coxmunk_n_water(surface, ispec, n_spec, FT)
        M = coxmunk_brdf_mueller(core_surface, 1, μᵥ, μ₀, Δϕ;
                                  n_water)
        ρ[iv, ispec] = M[1, 1]
    end
    return ρ
end

function _precompute_surface_brdf(::Val{N}, surface::CoxMunkSSSurface,
                                  geometry::SSGeometry, n_spec::Int,
                                  ::Type{FT}) where {N,FT}
    n_geom = length(geometry.μv)
    ρ = zeros(FT, n_geom, N, n_spec)
    μ₀ = convert(FT, geometry.μ₀)
    core_surface = _coxmunk_core_surface(surface, FT)
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

function _coxmunk_core_surface(surface::CoxMunkSSSurface, ::Type{FT}) where {FT}
    return CoxMunkSurface{FT}(
        wind_speed=convert(FT, surface.wind_speed),
        n_water=nothing,
        whitecap_albedo=convert(FT, surface.whitecap_albedo),
        include_whitecaps=surface.include_whitecaps,
        shadowing=surface.shadowing)
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
