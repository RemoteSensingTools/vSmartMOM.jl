# Azimuthal Fourier-series convergence helpers.

"Highest Fourier order completed by the most recent forward/linearized solve."
const _LAST_FOURIER_M_USED = Ref{Int}(-1)

_fourier_convergence_active(::AllFourierMoments) = false
_fourier_convergence_active(::IntensityConvergence) = true
_fourier_convergence_active(::StokesConvergence) = true

"""
Highest order protected from convergence testing for this solve.

`convergence.min_m - 1` is the configured low-order safety guard. Clamping it
to `m_max` means a physically shorter series is evaluated in full rather than
being extended solely to exercise the runtime test.
"""
@inline _fourier_guard_through_m(convergence, m_max::Int) =
    min(convergence.min_m - 1, m_max)

"Drop unavailable diagnostics while retaining a stable tuple of output arrays."
_fourier_outputs(outputs...) = Tuple(output for output in outputs if output !== nothing)

"Create previous-partial-sum snapshots used to isolate the next m term."
_fourier_snapshots(outputs::Tuple) = map(copy, outputs)

"Copy current accumulated outputs into reusable previous-moment snapshots."
function _update_fourier_snapshots!(snapshots::Tuple, outputs::Tuple)
    length(snapshots) == length(outputs) || throw(DimensionMismatch(
        "Fourier convergence snapshot/output counts differ"))
    foreach(copyto!, snapshots, outputs)
    return snapshots
end

@inline _fourier_component_indices(::IntensityConvergence, output) = 1:1
@inline _fourier_component_indices(::StokesConvergence, output) = axes(output, 2)

"""
Return `true` when the latest contribution passes for every component selected
by `convergence`, at every view and wavelength in every supplied output.
Directional radiance arrays use the package layout `(view, stokes, spectral)`.
"""
function _fourier_outputs_pass(convergence::Union{IntensityConvergence,StokesConvergence},
                               outputs::Tuple,
                               snapshots::Tuple)
    length(outputs) == length(snapshots) || throw(DimensionMismatch(
        "Fourier convergence snapshot/output counts differ"))
    tolerance = convergence.tolerance
    for (output, previous) in zip(outputs, snapshots)
        ndims(output) == 3 || throw(DimensionMismatch(
            "Fourier convergence expects (view, stokes, spectral) arrays; " *
            "received size $(size(output))"))
        size(output) == size(previous) || throw(DimensionMismatch(
            "Fourier convergence snapshot shape differs from its output"))
        for component in _fourier_component_indices(convergence, output)
            current_component = selectdim(output, 2, component)
            previous_component = selectdim(previous, 2, component)
            @inbounds for index in eachindex(current_component, previous_component)
                delta = current_component[index] - previous_component[index]
                abs(delta) <= tolerance * abs(current_component[index]) || return false
            end
        end
    end
    return true
end

"Update the pass counter and report whether the Fourier loop may stop."
function _fourier_convergence_step!(
    convergence::Union{IntensityConvergence,StokesConvergence},
    outputs::Tuple,
    snapshots::Tuple,
    consecutive_passes::Int,
    m::Int,
    m_max::Int)
    if m <= _fourier_guard_through_m(convergence, m_max)
        _update_fourier_snapshots!(snapshots, outputs)
        return false, 0
    end
    passed = _fourier_outputs_pass(convergence, outputs, snapshots)
    _update_fourier_snapshots!(snapshots, outputs)
    consecutive_passes = passed ? consecutive_passes + 1 : 0
    return consecutive_passes >= convergence.n_consecutive, consecutive_passes
end

"""
    _fourier_proxy_step!(convergence, accumulated, J₀⁻, ...)

Atmosphere-only convergence step used while constructing a reusable
`AtmosphereRTCache`. Postprocessing is intentionally skipped on that path, so
this helper extracts the current moment at the requested view nodes using the
same nearest-stream indices and azimuthal Stokes weights as production
postprocessing. The cache callback runs before this test, ensuring the retained
moment set and the accumulated proxy remain synchronized.

The proxy sees the atmosphere-only source. It is exact for Lambertian-family
surface replays (their `m>0` surface term is zero) and must be overridden to a
full series when the cache targets an azimuthally structured surface.
"""
function _fourier_proxy_step!(
    convergence::Union{IntensityConvergence,StokesConvergence},
    accumulated,
    J₀⁻,
    vza,
    vaz,
    qp_μ,
    pol_type,
    m::Int,
    weight,
    consecutive_passes::Int,
    m_max::Int)
    vza_info = _precompute_vza_weights(vza, vaz, qp_μ, pol_type, m, weight)
    source = Array(J₀⁻)
    passed = true
    @inbounds for i in eachindex(vza)
        istart, iend, w = vza_info[i]
        contribution = w * @view(source[istart:iend, 1, :])
        current = @view accumulated[i, :, :]
        current .+= contribution
        if m > _fourier_guard_through_m(convergence, m_max)
            for component in _fourier_component_indices(convergence, accumulated)
                for spectral_index in axes(current, 2)
                    delta = contribution[component, spectral_index]
                    new = current[component, spectral_index]
                    passed &= abs(delta) <= convergence.tolerance * abs(new)
                end
            end
        end
    end
    if m <= _fourier_guard_through_m(convergence, m_max)
        return false, 0
    end
    consecutive_passes = passed ? consecutive_passes + 1 : 0
    return consecutive_passes >= convergence.n_consecutive, consecutive_passes
end
