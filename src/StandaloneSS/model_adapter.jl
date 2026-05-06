"""
    exact_ss_config_from_model(model::RTModel; i_band=1, I0=nothing)

Build a standalone exact single-scattering configuration from a precomputed
CoreRT model. The adapter is intentionally conservative: it reuses the model's
Rayleigh Greek coefficients, absorption optical depths, supported surface
types, and any aerosol optics whose truncation factor is zero.

If aerosol optics have already been δ-M/δ-BGE truncated (`fᵗ != 0`), the
adapter throws instead of silently mixing truncated phase matrices with
unmodified optical-depth state.
"""
function exact_ss_config_from_model(model::RTModel; i_band::Integer = 1,
                                    I0 = nothing)
    FT = float_type(model)
    n_spec = length(model.atmosphere.spec_bands[i_band])
    I0_vec = I0 === nothing ? ones(FT, n_spec) : I0
    return ExactSSConfig(
        geometry = _ss_geometry_from_model(model, FT),
        surface = _ss_surface_from_core(model.surfaces[i_band], FT),
        contributors = _ss_contributors_from_model(model, i_band, n_spec, FT),
        I0 = I0_vec,
        polarization_type = model.solver.polarization_type,
        architecture = model.architecture)
end

function _ss_geometry_from_model(model::RTModel, ::Type{FT}) where {FT}
    obs = model.obs_geom
    return SSGeometry(
        μ₀ = cosd(convert(FT, obs.sza)),
        μv = cosd.(FT.(obs.vza)),
        Δϕ = deg2rad.(FT.(obs.vaz)))
end

_ss_surface_from_core(surface::LambertianSurfaceScalar, ::Type{FT}) where {FT} =
    LambertianSSSurface(albedo = convert(FT, surface.albedo))

_ss_surface_from_core(surface::LambertianSurfaceSpectrum, ::Type{FT}) where {FT} =
    LambertianSSSurface(albedo = FT.(surface.albedo))

function _ss_surface_from_core(surface::CoxMunkSurface, ::Type{FT}) where {FT}
    return CoxMunkSSSurface(
        wind_speed = convert(FT, surface.wind_speed),
        n_water = _convert_n_water(surface.n_water, FT),
        whitecap_albedo = convert(FT, surface.whitecap_albedo),
        include_whitecaps = surface.include_whitecaps,
        shadowing = surface.shadowing)
end

_ss_surface_from_core(surface, ::Type) =
    throw(ArgumentError("exact_ss_config_from_model does not yet support surface $(typeof(surface))"))

_convert_n_water(::Nothing, ::Type{FT}) where {FT} = nothing
_convert_n_water(n::Complex, ::Type{FT}) where {FT} =
    Complex{FT}(convert(FT, real(n)), convert(FT, imag(n)))
_convert_n_water(n::AbstractVector, ::Type{FT}) where {FT} =
    Complex{FT}.(FT.(real.(n)), FT.(imag.(n)))

function _layer_spec_matrix(τ_spec_layer::AbstractMatrix, ::Type{FT}) where {FT}
    return permutedims(FT.(τ_spec_layer), (2, 1))
end

function _aerosol_layer_spec_matrix(τ_aer_band::AbstractMatrix, i_aer::Int,
                                    n_spec::Int, ::Type{FT}) where {FT}
    τ_layer = FT.(view(τ_aer_band, i_aer, :))
    return repeat(reshape(τ_layer, :, 1), 1, n_spec)
end

function _aerosol_layer_spec_matrix(τ_aer_band::AbstractArray{<:Real,3},
                                    i_aer::Int, n_spec::Int,
                                    ::Type{FT}) where {FT}
    return permutedims(FT.(view(τ_aer_band, i_aer, :, :)), (2, 1))
end

_nonzero_optical_depth(τ) = any(!iszero, τ)
_zero_truncation_factor(f::Real) = iszero(f)
_zero_truncation_factor(f) = all(iszero, f)

function _ss_contributors_from_model(model::RTModel, i_band::Integer,
                                     n_spec::Int, ::Type{FT}) where {FT}
    contributors = AbstractSSContributor[]

    τ_rayl = _layer_spec_matrix(model.τ_rayl[i_band], FT)
    if _nonzero_optical_depth(τ_rayl)
        push!(contributors, GreekCoefsSSContributor(
            greek_coefs = model.greek_rayleigh[i_band],
            ϖ = one(FT),
            τ = τ_rayl))
    end

    τ_abs = _layer_spec_matrix(model.τ_abs[i_band], FT)
    if _nonzero_optical_depth(τ_abs)
        push!(contributors, AbsorptionSSContributor(τ = τ_abs))
    end

    τ_aer_band = model.τ_aer[i_band]
    for i_aer in axes(τ_aer_band, 1)
        τ_aer = _aerosol_layer_spec_matrix(τ_aer_band, i_aer, n_spec, FT)
        _nonzero_optical_depth(τ_aer) || continue

        optics = model.aerosol_optics[i_band][i_aer]
        _zero_truncation_factor(optics.fᵗ) ||
            throw(ArgumentError(
                "exact_ss_config_from_model currently supports only untruncated aerosol optics (fᵗ == 0)."))
        optics.ω̃ isa Real ||
            throw(ArgumentError(
                "exact_ss_config_from_model currently supports scalar aerosol single-scattering albedo only."))

        push!(contributors, GreekCoefsSSContributor(
            greek_coefs = optics.greek_coefs,
            ϖ = convert(FT, optics.ω̃),
            τ = τ_aer))
    end

    isempty(contributors) &&
        throw(ArgumentError("model has no nonzero Rayleigh, aerosol, or absorption optical depth in band $i_band"))
    return Tuple(contributors)
end
