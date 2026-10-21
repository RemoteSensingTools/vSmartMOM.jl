"""
    supports_surface_sif(surface) -> Bool

Whether the surface implements fluorescence injection through
`surface_source_contribute!` and `surface_source_contribute_lin!`.
The scalar, spectral, Legendre, and spline Lambertian surfaces support SIF.
An extension implementing both contributions for another surface should also
specialize this trait. The default is `false`.
"""
supports_surface_sif(::AbstractSurfaceType) = false
supports_surface_sif(::Union{LambertianSurfaceScalar, LambertianSurfaceSpectrum,
                            LambertianSurfaceLegendre, LambertianSurfaceSpline}) = true

# A single direct-beam carrier is shared by elemental and surface kernels.
# BlackbodySource returns SolarBeam, so it uses the same counting rule.
_solar_beam_count(::AbstractSource) = 0
_solar_beam_count(::Union{SolarBeam,PreparedSolarBeam}) = 1
_solar_beam_count(s::SourceSet) = sum(_solar_beam_count, s.sources; init=0)

_validate_source_geometry(::AbstractSource, sza) = nothing
function _validate_source_geometry(s::SolarBeam, sza)
    s.sza === nothing && return nothing
    # Compare at the narrower precision, in either direction. Exact equality
    # there admits representation rounding, not a different physical angle.
    source_angle, model_angle = float(s.sza), float(sza)
    FT = precision(source_angle) <= precision(model_angle) ?
         typeof(source_angle) : typeof(model_angle)
    isfinite(source_angle) && isfinite(model_angle) &&
        FT(source_angle) == FT(model_angle) || throw(ArgumentError(
        "SolarBeam(sza=$(s.sza)) does not match the model solar zenith angle " *
        "$sza degrees. Set params.sza before model_from_parameters, or use " *
        "remake_geometry; source-specific beam geometries are unsupported."))
    return nothing
end
function _validate_source_geometry(s::SourceSet, sza)
    foreach(src -> _validate_source_geometry(src, sza), s.sources)
    return nothing
end

# A prescribed zero is a harmless placeholder. Retrievable SIF requires the
# injection path even at zero radiance, because its tangent is nonzero.
_requires_surface_sif(s::SurfaceSIF) = surface_sif_parameter_count(s) > 0 ||
    (s.SIF₀ !== nothing && !iszero(s.SIF₀))
_requires_surface_sif(s::PreparedSurfaceSIF) = s.n_parameters > 0 || !iszero(s.SIF₀)

validate_source_surface(::Union{AbstractSource,AbstractPreparedSource}, ::AbstractSurfaceType) = nothing
function validate_source_surface(s::Union{SurfaceSIF,PreparedSurfaceSIF}, surface::AbstractSurfaceType)
    if _requires_surface_sif(s) && !supports_surface_sif(surface)
        throw(ArgumentError(
            "SurfaceSIF with nonzero emission or retrievable coefficients is " *
            "unsupported for $(typeof(surface)). Use a scalar, spectral, " *
            "Legendre, or spline Lambertian surface with SIF injection support."))
    end
    return nothing
end
function validate_source_surface(s::SourceSet, surface::AbstractSurfaceType)
    foreach(src -> validate_source_surface(src, surface), s.sources)
    return nothing
end

"""
    validate_source_requests(sources, sza, surface)

Check the source contract at the solve boundary, before preparing source arrays.
The current solver carries at most one collimated beam, whose optional SZA must
agree with the model geometry. Surface fluorescence requires an implemented
source/surface contribution. Split surface replay also calls
`validate_source_surface` on its cached prepared sources and new surface.
"""
function validate_source_requests(sources::AbstractSource, sza, surface::AbstractSurfaceType)
    _solar_beam_count(sources) <= 1 || throw(ArgumentError(
        "Multiple SolarBeam sources are unsupported: the solver carries one " *
        "direct beam. For a shared geometry, sum their irradiance spectra " *
        "into one SolarBeam(F₀=F₀₁ + F₀₂). BlackbodySource also counts as a beam."))
    _validate_source_geometry(sources, sza)
    validate_source_surface(sources, surface)
    validate_sif_solar_spectrum(sources)
    return nothing
end
