# Diagnostic bridge from GCHP/TOMAS microphysics to Mie extinction only.
# This does not produce phase matrices or wire aerosols into
# `model_from_parameters`.

const _DEFAULT_RI_DATABASES = Dict{DataType, Any}()
const _GAUSSLEGENDRE_CACHE = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}()
const _DEFAULT_RI_PATH = joinpath(@__DIR__, "..", "..", "data",
                                  "refractive_indices_database.yaml")

function _default_refractive_index_database(::Type{FT}) where {FT}
    db = get!(_DEFAULT_RI_DATABASES, FT) do
        load_refractive_index_database(_DEFAULT_RI_PATH, FT)
    end
    return db::RefractiveIndexDatabase{FT}
end

function _grid_radius_um(grid::AerosolSizeGrid{FT}, ibin::Integer) where {FT}
    value = grid.centers[ibin]
    value_um = if grid.units in ("μm", "um", "micrometer", "micrometers")
        value
    elseif grid.units in ("nm", "nanometer", "nanometers")
        value / FT(1000)
    elseif grid.units in ("m", "meter", "meters")
        value * FT(1e6)
    else
        throw(ArgumentError("Cannot convert aerosol size-grid units '$(grid.units)' to μm"))
    end

    if grid.coordinate === :diameter
        return value_um / FT(2)
    elseif grid.coordinate === :radius
        return value_um
    else
        throw(ArgumentError("AOD diagnostic needs radius/diameter size grid; got coordinate :$(grid.coordinate)"))
    end
end

function _size_to_um(value::FT, units::AbstractString) where {FT}
    if units in ("μm", "um", "micrometer", "micrometers")
        return value
    elseif units in ("nm", "nanometer", "nanometers")
        return value / FT(1000)
    elseif units in ("m", "meter", "meters")
        return value * FT(1e6)
    else
        throw(ArgumentError("Cannot convert aerosol size-grid units '$units' to μm"))
    end
end

function _bin_log_bounds(grid::AerosolSizeGrid{FT}, ibin::Integer) where {FT}
    grid.edges === nothing &&
        throw(ArgumentError("bin integration requires explicit size-grid edges"))
    lo = _size_to_um(grid.edges[ibin], grid.units)
    hi = _size_to_um(grid.edges[ibin + 1], grid.units)
    lo > zero(FT) && hi > lo ||
        throw(ArgumentError("invalid aerosol bin edges for bin $ibin: ($lo, $hi) $(grid.units)"))
    return log10(lo), log10(hi)
end

function _radius_from_log_size(grid::AerosolSizeGrid{FT}, log_size::FT) where {FT}
    value_um = FT(10)^log_size
    if grid.coordinate === :diameter
        return value_um / FT(2)
    elseif grid.coordinate === :radius
        return value_um
    else
        throw(ArgumentError("AOD diagnostic needs radius/diameter size grid; got coordinate :$(grid.coordinate)"))
    end
end

function _bin_log_center(grid::AerosolSizeGrid{FT}, ibin::Integer) where {FT}
    lo, hi = _bin_log_bounds(grid, ibin)
    return (lo + hi) / FT(2)
end

function _bin_log_width(grid::AerosolSizeGrid{FT}, ibin::Integer) where {FT}
    lo, hi = _bin_log_bounds(grid, ibin)
    return hi - lo
end

function _bin_mean_dndlog(aerosols::SectionalAerosolData{FT},
                          ibin::Integer, ilev::Integer) where {FT}
    width = _bin_log_width(aerosols.size_grid, ibin)
    width > zero(FT) || return zero(FT)
    return aerosols.number[ibin, ilev] * FT(1e6) / width
end

function _limited_linear_slope(aerosols::SectionalAerosolData{FT},
                               ibin::Integer, ilev::Integer,
                               mean_i::FT,
                               xlo::FT, xmid::FT, xhi::FT) where {FT}
    if ibin == 1
        mean_i <= zero(FT) && return zero(FT)
        return mean_i / (xmid - xlo)
    elseif ibin == nbins(aerosols)
        mean_i <= zero(FT) && return zero(FT)
        return -mean_i / (xhi - xmid)
    end

    xprev = _bin_log_center(aerosols.size_grid, ibin - 1)
    xnext = _bin_log_center(aerosols.size_grid, ibin + 1)
    mean_prev = _bin_mean_dndlog(aerosols, ibin - 1, ilev)
    mean_next = _bin_mean_dndlog(aerosols, ibin + 1, ilev)
    s_left = (mean_i - mean_prev) / (xmid - xprev)
    s_right = (mean_next - mean_i) / (xnext - xmid)

    slope = if s_left * s_right <= zero(FT)
        zero(FT)
    else
        sign(s_left) * min(abs(s_left), abs(s_right))
    end

    if mean_i > zero(FT)
        max_abs_slope = mean_i / max(abs(xlo - xmid), abs(xhi - xmid))
        slope = clamp(slope, -max_abs_slope, max_abs_slope)
    else
        slope = zero(FT)
    end
    return slope
end

function _gausslegendre_nodes_weights(nquad::Integer)
    return get!(_GAUSSLEGENDRE_CACHE, Int(nquad)) do
        gausslegendre(Int(nquad))
    end
end

function _layer_thickness_m(p_top_hpa::FT, p_bot_hpa::FT,
                            T::FT, q::FT) where {FT}
    Rd = FT(287.05)
    g = FT(9.80665)
    p_top = max(p_top_hpa * FT(100), FT(1e-3))
    p_bot = max(p_bot_hpa * FT(100), p_top * (one(FT) + eps(FT)))
    Tv = T * (one(FT) + FT(0.61) * q)
    return max(zero(FT), Rd * Tv / g * log(p_bot / p_top))
end

function _mie_extinction_cross_section(radius_um::FT,
                                       wavelength_um::FT,
                                       n_eff::Complex{FT}) where {FT}
    radius_um > zero(FT) || return zero(FT)
    wavelength_um > zero(FT) || throw(ArgumentError("wavelength must be positive"))

    x = FT(2π) * radius_um / wavelength_um
    x > zero(FT) || return zero(FT)

    n_max = Scattering.get_n_max(x)
    m_ref = Complex{FT}(real(n_eff), -imag(n_eff))
    y = x * m_ref
    nmx = round(Int, max(FT(n_max), FT(abs(y))) + FT(51))

    an = zeros(Complex{FT}, n_max)
    bn = zeros(Complex{FT}, n_max)
    Dn = zeros(Complex{FT}, nmx)
    Scattering.compute_mie_ab!(x, m_ref, an, bn, Dn)

    k = FT(2π) / wavelength_um
    series = zero(FT)
    @inbounds for n in 1:n_max
        series += FT(2n + 1) * real(an[n] + bn[n])
    end
    c_ext_um2 = FT(2π) / k^2 * series
    return max(zero(FT), c_ext_um2) * FT(1e-12)
end

"""
    compute_column_aod(aerosols, p_half, T, q, wavelengths_um; ...)

Compute column aerosol optical depth for a sectional aerosol column without
entering the CoreRT solver. `p_half`, `T`, and `q` are TOA-to-BOA. The
returned vector follows `wavelengths_um`.
"""
function compute_column_aod(aerosols::SectionalAerosolData{FT},
                            p_half::AbstractVector,
                            T::AbstractVector,
                            q::AbstractVector,
                            wavelengths_um::AbstractVector;
                            ri_database::Union{Nothing, RefractiveIndexDatabase{FT}} = aerosols.ri_database,
                            mixing_rule::AbstractMixingRule = aerosols.mixing_rule,
                            integration::AbstractBinIntegration = aerosols.integration,
                            mie_cache = nothing,
                            mie_lut = nothing,
                            ri_round_digits::Union{Nothing, Integer} = nothing) where {FT}
    nlay = nlayers(aerosols)
    length(p_half) == nlay + 1 ||
        throw(ArgumentError("p_half length must be nlayers + 1"))
    length(T) == nlay ||
        throw(ArgumentError("T length must match aerosol nlayers"))
    length(q) == nlay ||
        throw(ArgumentError("q length must match aerosol nlayers"))

    db = ri_database === nothing ? _default_refractive_index_database(FT) : ri_database
    λs = FT.(wavelengths_um)
    out = zeros(FT, length(λs))

    dz = [_layer_thickness_m(FT(p_half[ilev]), FT(p_half[ilev + 1]),
                             FT(T[ilev]), FT(q[ilev]))
          for ilev in 1:nlay]

    @inbounds for (iλ, λ) in enumerate(λs)
        τ = zero(FT)
        for ilev in 1:nlay
            dz_layer = dz[ilev]
            dz_layer > zero(FT) || continue
            for ibin in 1:nbins(aerosols)
                aerosols.number[ibin, ilev] > zero(FT) || continue

                comp = bin_composition(aerosols, ibin, ilev; basis=:volume)
                sum(comp.fractions) > zero(FT) || continue
                n_eff = effective_ri(comp, db, λ, mixing_rule, aerosols.scheme)
                τ += _integrated_bin_extinction(aerosols, ibin, ilev, λ,
                                                n_eff, integration,
                                                mie_cache, mie_lut,
                                                ri_round_digits) *
                     dz_layer
            end
        end
        out[iλ] = max(zero(FT), τ)
    end

    return out
end

function _integrated_bin_extinction(aerosols::SectionalAerosolData{FT},
                                    ibin::Integer, ilev::Integer,
                                    λ::FT, n_eff::Complex{FT},
                                    ::DirectBinSum,
                                    mie_cache, mie_lut, ri_round_digits) where {FT}
    N_m3 = aerosols.number[ibin, ilev] * FT(1e6)
    r_um = _grid_radius_um(aerosols.size_grid, ibin)
    c_ext = _cached_mie_extinction(r_um, λ, n_eff, mie_cache, mie_lut,
                                   ri_round_digits)
    return N_m3 * c_ext
end

function _integrated_bin_extinction(aerosols::SectionalAerosolData{FT},
                                    ibin::Integer, ilev::Integer,
                                    λ::FT, n_eff::Complex{FT},
                                    integration::ConstantIntegrationPerBin,
                                    mie_cache, mie_lut, ri_round_digits) where {FT}
    xlo, xhi = _bin_log_bounds(aerosols.size_grid, ibin)
    xmid = (xlo + xhi) / FT(2)
    half_width = (xhi - xlo) / FT(2)
    mean_dndlog = _bin_mean_dndlog(aerosols, ibin, ilev)
    nodes, weights = _gausslegendre_nodes_weights(integration.nquad)

    β = zero(FT)
    @inbounds for i in eachindex(nodes)
        x = xmid + half_width * FT(nodes[i])
        r_um = _radius_from_log_size(aerosols.size_grid, x)
        c_ext = _cached_mie_extinction(r_um, λ, n_eff, mie_cache, mie_lut,
                                       ri_round_digits)
        β += FT(weights[i]) * mean_dndlog * c_ext
    end
    return half_width * β
end

function _integrated_bin_extinction(aerosols::SectionalAerosolData{FT},
                                    ibin::Integer, ilev::Integer,
                                    λ::FT, n_eff::Complex{FT},
                                    integration::LinearIntegrationPerBin,
                                    mie_cache, mie_lut, ri_round_digits) where {FT}
    xlo, xhi = _bin_log_bounds(aerosols.size_grid, ibin)
    xmid = (xlo + xhi) / FT(2)
    half_width = (xhi - xlo) / FT(2)
    mean_dndlog = _bin_mean_dndlog(aerosols, ibin, ilev)
    slope = _limited_linear_slope(aerosols, ibin, ilev, mean_dndlog,
                                  xlo, xmid, xhi)
    nodes, weights = _gausslegendre_nodes_weights(integration.nquad)

    β = zero(FT)
    @inbounds for i in eachindex(nodes)
        x = xmid + half_width * FT(nodes[i])
        dndlog = mean_dndlog + slope * (x - xmid)
        dndlog = max(zero(FT), dndlog)
        r_um = _radius_from_log_size(aerosols.size_grid, x)
        c_ext = _cached_mie_extinction(r_um, λ, n_eff, mie_cache, mie_lut,
                                       ri_round_digits)
        β += FT(weights[i]) * dndlog * c_ext
    end
    return half_width * β
end

_integrated_bin_extinction(::SectionalAerosolData, _ibin, _ilev, _λ, _n_eff,
                           ::LogNormalFit, _mie_cache, _mie_lut, _ri_round_digits) =
    error("LogNormalFit AOD integration is not implemented for sectional data yet")

function _lut_mie_extinction(radius_um::FT,
                             wavelength_um::FT,
                             n_eff::Complex{FT},
                             mie_lut) where {FT}
    x = FT(2π) * radius_um / wavelength_um
    n_real = real(n_eff)
    n_imag = max(zero(FT), imag(n_eff))
    q_ext = mie_lut.q_ext(x, n_real, n_imag)
    return max(zero(FT), FT(q_ext)) * FT(π) * radius_um^2 * FT(1e-12)
end

function _cached_mie_extinction(radius_um::FT,
                                wavelength_um::FT,
                                n_eff::Complex{FT},
                                cache,
                                mie_lut,
                                ri_round_digits) where {FT}
    if mie_lut !== nothing
        return _lut_mie_extinction(radius_um, wavelength_um, n_eff, mie_lut)
    end

    cache === nothing && return _mie_extinction_cross_section(radius_um, wavelength_um, n_eff)

    key = if ri_round_digits === nothing
        (radius_um, wavelength_um, real(n_eff), imag(n_eff))
    else
        digits = Int(ri_round_digits)
        (radius_um, wavelength_um,
         round(real(n_eff); digits=digits),
         round(imag(n_eff); digits=digits))
    end

    return get!(cache, key) do
        n_cached = Complex{FT}(FT(key[3]), FT(key[4]))
        _mie_extinction_cross_section(radius_um, wavelength_um, n_cached)
    end
end
