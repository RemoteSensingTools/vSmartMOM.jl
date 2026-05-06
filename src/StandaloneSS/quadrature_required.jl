"""
    determine_required_l_from_moments(β; target_relative_error=1e-4)

Estimate the highest Greek/Legendre moment needed from a supplied coefficient
series by comparing the squared tail energy to the requested relative error.
The input is indexed as Julia arrays are: `β[1]` is the `l=0` coefficient.
"""
function determine_required_l_from_moments(β::AbstractVector{<:Real};
                                           target_relative_error::Real = 1e-4)
    isempty(β) && throw(ArgumentError("moment coefficient vector must not be empty"))
    target_relative_error > 0 ||
        throw(ArgumentError("target_relative_error must be positive"))
    length(β) == 1 && return 0

    β² = abs2.(β)
    tail_cumulative = reverse(cumsum(reverse(β²)))
    total_energy = tail_cumulative[1]
    total_energy == 0 && return 0
    threshold = total_energy * target_relative_error^2

    for l in 0:(length(β) - 2)
        tail_after_l = tail_cumulative[l + 2]
        tail_after_l <= threshold && return l
    end
    return length(β) - 1
end

"""
    determine_required_l_aerosol(contributor; target_relative_error=1e-4)

Return the required phase-function moment order for a Phase 1 contributor.
Rayleigh is exact through `l=2`; pure absorption contributes no phase series.
For the Henyey-Greenstein contributor, the estimate uses the closed-form
asymptotic `g^l` decay.
"""
determine_required_l_aerosol(::RayleighSSContributor;
                             target_relative_error::Real = 1e-4) = 2

determine_required_l_aerosol(::AbsorptionSSContributor;
                             target_relative_error::Real = 1e-4) = 0

function determine_required_l_aerosol(c::HGAerosolSSContributor;
                                      target_relative_error::Real = 1e-4)
    target_relative_error > 0 ||
        throw(ArgumentError("target_relative_error must be positive"))
    g = abs(c.g)
    g < 1 || throw(ArgumentError("Henyey-Greenstein |g| must be < 1"))
    g == 0 && return 0

    # If coefficients decay like g^l, the omitted tail is O(g^L).
    return max(0, ceil(Int, log(target_relative_error) / log(g)))
end

determine_required_l_aerosol(c::GreekCoefsSSContributor;
                             target_relative_error::Real = 1e-4) =
    determine_required_l_from_moments(c.greek_coefs.β; target_relative_error)

"""
    determine_required_nbrdf(surface; target_relative_error=1e-4)

Return the required BRDF Fourier moment count for a surface. Lambertian
surfaces only need `m=0`, represented as a count of 1.
"""
determine_required_nbrdf(::LambertianSSSurface;
                         target_relative_error::Real = 1e-4) = 1

determine_required_nbrdf(surface::CoxMunkSSSurface;
                         target_relative_error::Real = 1e-4) =
    determine_required_nbrdf_coxmunk(surface.wind_speed; target_relative_error)

"""
    determine_required_nbrdf_coxmunk(wind_speed_ms; target_relative_error=1e-4)

Heuristic Cox-Munk BRDF Fourier moment count from the capillary-slope width.
This helper is standalone and does not depend on the CoreRT CoxMunk type.
"""
function determine_required_nbrdf_coxmunk(wind_speed_ms::Real;
                                          target_relative_error::Real = 1e-4)
    wind_speed_ms >= 0 || throw(ArgumentError("wind_speed_ms must be non-negative"))
    target_relative_error > 0 ||
        throw(ArgumentError("target_relative_error must be positive"))
    σ² = 0.003 + 0.00512 * wind_speed_ms
    m_max = ceil(Int, log(1 / target_relative_error) / sqrt(σ²))
    return clamp(m_max, 16, 200)
end

function _max_required_l(contributors; target_relative_error::Real)
    maximum((determine_required_l_aerosol(c; target_relative_error)
             for c in contributors); init=0)
end

"""
    determine_required_nstreams(contributors, surface; target_relative_error=1e-4)

Return a diagnostic record for the discrete-ordinate stream count implied by
Phase 1 contributors and the surface BRDF.
"""
function determine_required_nstreams(contributors, surface;
                                     target_relative_error::Real = 1e-4)
    L_aerosol = _max_required_l(contributors; target_relative_error)
    Nstreams_aero = ceil(Int, (L_aerosol + 1) / 2)
    Nbrdf = determine_required_nbrdf(surface; target_relative_error)
    return (; Nstreams=max(Nstreams_aero, 1),
            NSTREAMS_BRDF=Nbrdf,
            L_aerosol,
            surface_kind=typeof(surface).name.name,
            target_relative_error)
end

_residual_asymmetry(::RayleighSSContributor) = 0
_residual_asymmetry(::AbsorptionSSContributor) = 0
_residual_asymmetry(c::HGAerosolSSContributor) = abs(c.g)
_residual_asymmetry(::GreekCoefsSSContributor) = 0

"""
    determine_required_nquad_inner(τ_total, contributors; target_relative_error=1e-4)

Estimate the smooth inner quadrature order needed by later paths 3+4. Phase 1
does not use inner quadrature, but this diagnostic is part of the standalone
solver API for the next path slice.
"""
function determine_required_nquad_inner(τ_total, contributors;
                                        target_relative_error::Real = 1e-4)
    target_relative_error > 0 ||
        throw(ArgumentError("target_relative_error must be positive"))
    τ_max = τ_total isa Number ? τ_total : maximum(τ_total)
    τ_max >= 0 || throw(ArgumentError("τ_total must be non-negative"))
    g_residual = maximum((_residual_asymmetry(c) for c in contributors); init=0)
    τ_factor = log(1 + τ_max)
    N = ceil(Int, 8 + 4 * g_residual + 6 * τ_factor)
    return clamp(N, 8, 64)
end

"""
    determine_required_nquad(config; target_relative_error=1e-4)

Combined Phase 2 diagnostic for the standalone solver: stream count,
surface Fourier count, and later paths-3+4 inner quadrature estimate.
"""
function determine_required_nquad(config::ExactSSConfig;
                                  target_relative_error::Real = 1e-4)
    nstream_diag = determine_required_nstreams(config.contributors, config.surface;
                                               target_relative_error)
    optics = _precompute_optics(config)
    inner = determine_required_nquad_inner(optics.τ_total_column, config.contributors;
                                           target_relative_error)
    return merge(nstream_diag, (; inner_nquad=inner))
end
