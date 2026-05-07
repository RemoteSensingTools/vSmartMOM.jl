_validate_l_trunc(l_trunc::Integer) =
    l_trunc >= 2 ? Int(l_trunc) :
    throw(ArgumentError("truncated_ss_path1 requires l_trunc >= 2"))

_validate_l_trunc(l_trunc) =
    throw(ArgumentError("truncated_ss_path1 requires an integer l_trunc; got $(typeof(l_trunc))"))

_validate_max_m(max_m::Integer) =
    max_m > 0 ? Int(max_m) :
    throw(ArgumentError("truncated_ss_path2 requires max_m > 0"))

_validate_max_m(max_m) =
    throw(ArgumentError("truncated_ss_path2 requires an integer max_m; got $(typeof(max_m))"))

function _truncate_greek(greek::Scattering.GreekCoefs, n_moments::Int)
    n = minimum((n_moments, length(greek.α), length(greek.β),
                 length(greek.γ), length(greek.δ), length(greek.ϵ),
                 length(greek.ζ)))
    n >= 2 ||
        throw(ArgumentError("Greek coefficient truncation requires at least two retained moments"))
    return Scattering.GreekCoefs(
        greek.α[1:n], greek.β[1:n], greek.γ[1:n],
        greek.δ[1:n], greek.ϵ[1:n], greek.ζ[1:n])
end

function _hg_truncated_greek(c::HGAerosolSSContributor, n_moments::Int)
    FT = promote_type(typeof(c.g), typeof(c.ϖ), _real_numeric_type(c.τ))
    g = convert(FT, c.g)
    abs(g) < one(FT) ||
        throw(ArgumentError("Henyey-Greenstein asymmetry parameter must satisfy abs(g) < 1"))
    n = _validate_l_trunc(n_moments)
    α = zeros(FT, n)
    β = zeros(FT, n)
    γ = zeros(FT, n)
    δ = zeros(FT, n)
    ϵ = zeros(FT, n)
    ζ = zeros(FT, n)
    @inbounds for l in 0:(n - 1)
        β[l + 1] = FT(2l + 1) * g^l
    end
    return Scattering.GreekCoefs(α, β, γ, δ, ϵ, ζ)
end

_moment_limited_contributor(c::AbsorptionSSContributor, ::Int) = c

function _moment_limited_contributor(c::RayleighSSContributor, n_moments::Int)
    FT = promote_type(typeof(c.depol), _real_numeric_type(c.τ))
    greek = Scattering.get_greek_rayleigh(convert(FT, c.depol))
    return GreekCoefsSSContributor(
        greek_coefs = _truncate_greek(greek, n_moments),
        ϖ = one(FT),
        τ = c.τ)
end

function _moment_limited_contributor(c::HGAerosolSSContributor, n_moments::Int)
    return GreekCoefsSSContributor(
        greek_coefs = _hg_truncated_greek(c, n_moments),
        ϖ = c.ϖ,
        τ = c.τ)
end

function _moment_limited_contributor(c::GreekCoefsSSContributor,
                                     n_moments::Int)
    return GreekCoefsSSContributor(
        greek_coefs = _truncate_greek(c.greek_coefs, n_moments),
        ϖ = c.ϖ,
        τ = c.τ)
end

_moment_limited_contributors(contributors, n_moments::Int) =
    map(c -> _moment_limited_contributor(c, n_moments), contributors)

function _moment_limited_config(config::ExactSSConfig, n_moments::Int)
    return ExactSSConfig(
        geometry = config.geometry,
        surface = config.surface,
        contributors = _moment_limited_contributors(config.contributors,
                                                    n_moments),
        I0 = config.I0,
        n_stokes = config.n_stokes,
        polarization_type = config.polarization_type,
        inner_nquad = config.inner_nquad,
        azimuth_nquad = config.azimuth_nquad,
        architecture = config.architecture)
end

"""
    truncated_ss_path1(config, l_trunc)

Run standalone atmospheric single-scatter path 1 with scattering phase
functions reconstructed from only the first `l_trunc` Greek/Legendre moments.
Optical depths and single-scattering albedos are unchanged. This is the
post-hoc first-order analogue of the MOM angular truncation, not a full δ-M or
δ-BGE optical-depth rescaling.
"""
function truncated_ss_path1(config::ExactSSConfig, l_trunc)
    n_moments = _validate_l_trunc(l_trunc)
    truncated_config = _moment_limited_config(config, n_moments)
    return run_exact_ss(truncated_config; paths = :path1).path1
end

@inline function _ss_azimuthal_kernel(si::Int, sj::Int, m::Int,
                                      dϕ::FT) where {FT}
    row_is_iq = si <= 2
    col_is_iq = sj <= 2
    return row_is_iq == col_is_iq ? cos(m * dϕ) : sin(m * dϕ)
end

function _gauss_legendre_0π(n::Int, ::Type{FT}) where {FT}
    n > 0 || throw(ArgumentError("azimuth quadrature must be positive"))
    x, w = gausslegendre(n)
    scale = FT(pi) / FT(2)
    return scale .* (FT.(x) .+ one(FT)), scale .* FT.(w)
end

function _coxmunk_fourier_coeff_element(core_surface::CoxMunkSurface{FT},
                                        ϕ, w, n_stokes::Int,
                                        si::Int, sj::Int,
                                        μᵥ::FT, μ₀::FT, m::Int,
                                        n_water::Complex{FT}) where {FT}
    result = zero(FT)
    @inbounds for iϕ in eachindex(ϕ)
        M = coxmunk_brdf_mueller(core_surface, n_stokes, μᵥ, μ₀, ϕ[iϕ];
                                  n_water)
        result += w[iϕ] * M[si, sj] *
                  _ss_azimuthal_kernel(si, sj, m, ϕ[iϕ])
    end
    ff = m == 0 ? one(FT) : FT(2)
    return ff * result / FT(pi)
end

function _precompute_surface_brdf_fourier(::Val{1},
                                          surface::CoxMunkSSSurface,
                                          geometry::SSGeometry,
                                          n_spec::Int, ::Type{FT},
                                          max_m::Int) where {FT}
    n_geom = length(geometry.μv)
    ρ = zeros(FT, n_geom, n_spec)
    μ₀ = convert(FT, geometry.μ₀)
    core_surface = _coxmunk_core_surface(surface, FT)
    ϕ, w = _gauss_legendre_0π(100, FT)
    @inbounds for ispec in 1:n_spec, iv in 1:n_geom
        μᵥ = convert(FT, geometry.μv[iv])
        dϕ = convert(FT, geometry.Δϕ[iv])
        n_water = _coxmunk_n_water(surface, ispec, n_spec, FT)
        value = zero(FT)
        for m in 0:(max_m - 1)
            weight_m = m == 0 ? FT(0.5) : one(FT)
            coeff = _coxmunk_fourier_coeff_element(
                core_surface, ϕ, w, 1, 1, 1, μᵥ, μ₀, m, n_water)
            value += weight_m * _ss_azimuthal_kernel(1, 1, m, dϕ) * coeff
        end
        ρ[iv, ispec] = value
    end
    return ρ
end

function _precompute_surface_brdf_fourier(::Val{N},
                                          surface::CoxMunkSSSurface,
                                          geometry::SSGeometry,
                                          n_spec::Int, ::Type{FT},
                                          max_m::Int) where {N,FT}
    n_geom = length(geometry.μv)
    ρ = zeros(FT, n_geom, N, n_spec)
    μ₀ = convert(FT, geometry.μ₀)
    core_surface = _coxmunk_core_surface(surface, FT)
    ϕ, w = _gauss_legendre_0π(100, FT)
    @inbounds for ispec in 1:n_spec, iv in 1:n_geom
        μᵥ = convert(FT, geometry.μv[iv])
        dϕ = convert(FT, geometry.Δϕ[iv])
        n_water = _coxmunk_n_water(surface, ispec, n_spec, FT)
        for si in 1:N
            value = zero(FT)
            for m in 0:(max_m - 1)
                weight_m = m == 0 ? FT(0.5) : one(FT)
                coeff = _coxmunk_fourier_coeff_element(
                    core_surface, ϕ, w, N, si, 1, μᵥ, μ₀, m, n_water)
                value += weight_m * _ss_azimuthal_kernel(si, 1, m, dϕ) * coeff
            end
            ρ[iv, si, ispec] = value
        end
    end
    return ρ
end

function _run_path2_with_surface_brdf(config::ExactSSConfig, surface_brdf)
    FT = _config_numeric_type(config)
    _validate_config(config, :path2)
    architecture = config.architecture
    backend = _architecture_backend(architecture)
    n_geom = length(config.geometry.μv)
    n_stokes = _config_n_stokes(config)
    stokes_dispatch = Val(n_stokes)
    optics = _precompute_optics(Val(1), config)
    n_spec = size(optics.τ_cum, 2)
    path2 = _path_zeros(architecture, FT, n_geom, n_stokes, n_spec)
    _run_path2_kernel!(
        stokes_dispatch,
        path2,
        _to_arch(architecture, optics.τ_total_column),
        config.geometry.μ₀,
        _to_arch(architecture, config.geometry.μv),
        _to_arch(architecture, surface_brdf),
        _to_arch(architecture, _vectorize_I0(config.I0, n_spec, FT)),
        backend)
    return path2
end

_precompute_truncated_surface_brdf(stokes_dispatch, surface::LambertianSSSurface,
                                   geometry, n_spec, FT, ::Int) =
    _precompute_surface_brdf(stokes_dispatch, surface, geometry, n_spec, FT)

_precompute_truncated_surface_brdf(stokes_dispatch, surface::CoxMunkSSSurface,
                                   geometry, n_spec, FT, max_m::Int) =
    _precompute_surface_brdf_fourier(stokes_dispatch, surface, geometry, n_spec,
                                     FT, max_m)

_precompute_truncated_surface_brdf(_, surface::AbstractSSSurface, args...) =
    throw(ArgumentError("truncated_ss_path2 supports LambertianSSSurface and " *
                        "CoxMunkSSSurface only; got $(typeof(surface))"))

"""
    truncated_ss_path2(config, max_m)

Run standalone direct-beam surface path 2 with the surface BRDF reconstructed
from Fourier moments `m = 0:max_m-1`. Lambertian surfaces are exactly
represented by `m = 0`; Cox-Munk surfaces use the same direct-beam Fourier
coefficient convention as the CoreRT surface correction.
"""
function truncated_ss_path2(config::ExactSSConfig, max_m)
    m_count = _validate_max_m(max_m)
    FT = _config_numeric_type(config)
    _validate_config(config, :path2)
    _, n_spec = _dims(config.contributors)
    n_stokes = _config_n_stokes(config)
    stokes_dispatch = Val(n_stokes)
    surface_brdf = _precompute_truncated_surface_brdf(
        stokes_dispatch, config.surface, config.geometry, n_spec, FT, m_count)
    return _run_path2_with_surface_brdf(config, surface_brdf)
end

function _ss_back_correction(config::ExactSSConfig, l_trunc, max_m)
    exact = run_exact_ss(config; paths = :paths_1_2).total
    truncated = truncated_ss_path1(config, l_trunc)
    truncated .+= truncated_ss_path2(config, max_m)
    return exact .- truncated
end

"""
    apply_back_correction!(R_SFI, config; l_trunc, max_m)
    apply_back_correction!(R_SFI, model; i_band=1, l_trunc=nothing, max_m=nothing, I0=nothing)

Apply the standalone first-order back-correction in place:
`exact(path1+path2) - truncated(path1+path2)`. This is a diagnostic/post-hoc
single-scattering correction and does not propagate a corrected source through
the multiple-scattering operator.
"""
function apply_back_correction!(R_SFI::AbstractArray, config::ExactSSConfig;
                                l_trunc, max_m)
    correction = _ss_back_correction(config, l_trunc, max_m)
    size(R_SFI) == size(correction) ||
        throw(DimensionMismatch("R_SFI has size $(size(R_SFI)); expected $(size(correction))"))
    R_SFI .+= correction
    return R_SFI
end

function apply_back_correction!(R_SFI::AbstractArray, model::RTModel;
                                i_band::Integer = 1, l_trunc = nothing,
                                max_m = nothing, I0 = nothing)
    ib = Int(i_band)
    l = l_trunc === nothing ? model.solver.l_max[ib] : l_trunc
    # The public `max_m` kwarg here is a count of Fourier moments; the
    # solver stores order in `m_max_bands`, so convert via +1 (Phase B).
    m = max_m === nothing ? n_fourier_moments_bands(model)[ib] : max_m
    config = exact_ss_config_from_model(model; i_band = ib, I0)
    return apply_back_correction!(R_SFI, config; l_trunc = l, max_m = m)
end
