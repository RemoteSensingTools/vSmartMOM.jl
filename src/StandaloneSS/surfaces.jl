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

"""
    surface_brdf_wind_jacobian(config)
    surface_brdf_wind_jacobian(surface, geometry, n_spec; polarization_type=nothing)

Return the StandaloneSS seam derivative of Cox-Munk direct-beam surface BRDF
with respect to wind speed. For scalar Stokes it returns an array with shape
`(nGeom, nSpec, 1)`. For vector Stokes it returns
`(nGeom, nStokes, nSpec, 1)`, matching
[`chain_rule_combine_surface_brdf`](@ref)'s `dρ_dp` input.
"""
function surface_brdf_wind_jacobian(config::ExactSSConfig)
    _, n_spec = _dims(config.contributors)
    return surface_brdf_wind_jacobian(
        config.surface, config.geometry, n_spec;
        polarization_type=config.polarization_type)
end

function surface_brdf_wind_jacobian(surface::CoxMunkSSSurface,
                                    geometry::SSGeometry,
                                    n_spec::Integer;
                                    polarization_type = nothing)
    n_spec > 0 ||
        throw(ArgumentError("surface_brdf_wind_jacobian requires n_spec > 0"))
    FT = promote_type(typeof(geometry.μ₀), eltype(geometry.μv),
                      eltype(geometry.Δϕ), _surface_numeric_type(surface))
    stokes_dispatch = _surface_jacobian_stokes_dispatch(polarization_type)
    return _precompute_surface_brdf_wind_jacobian(
        stokes_dispatch, surface, geometry, Int(n_spec), FT)
end

surface_brdf_wind_jacobian(surface::AbstractSSSurface, args...; kwargs...) =
    throw(ArgumentError("surface_brdf_wind_jacobian currently supports " *
                        "CoxMunkSSSurface only; got $(typeof(surface))"))

_surface_jacobian_stokes_dispatch(::Nothing) = Val(1)
_surface_jacobian_stokes_dispatch(::Scattering.Stokes_I) = Val(1)
_surface_jacobian_stokes_dispatch(::Scattering.Stokes_IQ) = Val(2)
_surface_jacobian_stokes_dispatch(::Scattering.Stokes_IQU) = Val(3)
_surface_jacobian_stokes_dispatch(::Scattering.Stokes_IQUV) = Val(4)

function _surface_jacobian_stokes_dispatch(polarization_type)
    hasproperty(polarization_type, :n) ||
        throw(ArgumentError("polarization_type must expose an `n` field"))
    return Val(Int(getproperty(polarization_type, :n)))
end

function _precompute_surface_brdf_wind_jacobian(
    surface::CoxMunkSSSurface,
    geometry::SSGeometry,
    n_spec::Int,
    ::Type{FT},
) where {FT}
    return _precompute_surface_brdf_wind_jacobian(
        Val(1), surface, geometry, n_spec, FT)
end

function _precompute_surface_brdf_wind_jacobian(
    ::Val{1},
    surface::CoxMunkSSSurface,
    geometry::SSGeometry,
    n_spec::Int,
    ::Type{FT},
) where {FT}
    n_geom = length(geometry.μv)
    dρ = zeros(FT, n_geom, n_spec, 1)
    μ₀ = convert(FT, geometry.μ₀)
    core_surface = _coxmunk_core_surface(surface, FT)
    @inbounds for ispec in 1:n_spec, iv in 1:n_geom
        μᵥ = convert(FT, geometry.μv[iv])
        Δϕ = convert(FT, geometry.Δϕ[iv])
        n_water = _coxmunk_n_water(surface, ispec, n_spec, FT)
        _, dM = coxmunk_brdf_mueller_and_deriv(
            core_surface, 1, μᵥ, μ₀, Δϕ; n_water)
        dρ[iv, ispec, 1] = dM[1, 1]
    end
    return dρ
end

function _precompute_surface_brdf_wind_jacobian(
    ::Val{N},
    surface::CoxMunkSSSurface,
    geometry::SSGeometry,
    n_spec::Int,
    ::Type{FT},
) where {N,FT}
    n_geom = length(geometry.μv)
    dρ = zeros(FT, n_geom, N, n_spec, 1)
    μ₀ = convert(FT, geometry.μ₀)
    core_surface = _coxmunk_core_surface(surface, FT)
    @inbounds for ispec in 1:n_spec, iv in 1:n_geom
        μᵥ = convert(FT, geometry.μv[iv])
        Δϕ = convert(FT, geometry.Δϕ[iv])
        n_water = _coxmunk_n_water(surface, ispec, n_spec, FT)
        _, dM = coxmunk_brdf_mueller_and_deriv(
            core_surface, N, μᵥ, μ₀, Δϕ; n_water)
        for istokes in 1:N
            dρ[iv, istokes, ispec, 1] = dM[istokes, 1]
        end
    end
    return dρ
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
