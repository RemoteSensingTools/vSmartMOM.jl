# Azimuthal Fourier-series convergence helpers.

"Highest Fourier order completed by the most recent forward/linearized solve."
const _LAST_FOURIER_M_USED = Ref{Int}(-1)

_fourier_convergence_active(::AllFourierMoments) = false
_fourier_convergence_active(::IntensityConvergence) = true
_fourier_convergence_active(::StokesConvergence) = true

"""
Highest order protected from convergence testing for this solve.

`convergence.min_m - 1` is the configured low-order safety guard (two for the
default `min_m=3`). Clamping it to the solver-resolved `m_max` means a series
with exact support below that guard is simply evaluated in full; convergence
never manufactures additional Fourier moments.
"""
@inline _fourier_guard_through_m(convergence, m_max::Int) =
    min(convergence.min_m - 1, m_max)

"Drop unavailable diagnostics while retaining a stable tuple of output arrays."
_fourier_outputs(outputs...) = Tuple(output for output in outputs if output !== nothing)

"Create the previous-partial-sum snapshots used to isolate the next m term."
_fourier_snapshots(outputs::Tuple) = map(copy, outputs)

"Copy current accumulated outputs into the reusable previous-moment snapshots."
function _update_fourier_snapshots!(snapshots::Tuple, outputs::Tuple)
    length(snapshots) == length(outputs) || throw(DimensionMismatch(
        "Fourier convergence snapshot/output counts differ"))
    foreach(copyto!, snapshots, outputs)
    return snapshots
end

"""
Return `true` when the latest contribution passes independently for every
Stokes component at every view and wavelength in every supplied output.
"""
function _fourier_stokes_passes(convergence::StokesConvergence,
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
        @inbounds for index in eachindex(output, previous)
            delta = output[index] - previous[index]
            abs(delta) <= tolerance * abs(output[index]) || return false
        end
    end
    return true
end

"""
Return `true` when the latest contribution passes the Stokes-I test at every
element of every supplied directional output. Directional radiance arrays use
the package layout `(view, stokes, spectral)`; the convergence test deliberately
does not inspect Q/U/V. The enclosing convergence strategy nevertheless
protects all available orders through `min(2, m_max)` before this test is
eligible under the default configuration.
"""
function _fourier_intensity_passes(convergence::IntensityConvergence,
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
        current_I = selectdim(output, 2, 1)
        previous_I = selectdim(previous, 2, 1)
        @inbounds for index in eachindex(current_I, previous_I)
            delta = current_I[index] - previous_I[index]
            abs(delta) <= tolerance * abs(current_I[index]) || return false
        end
    end
    return true
end

function _fourier_convergence_step!(convergence::StokesConvergence,
                                    outputs::Tuple,
                                    snapshots::Tuple,
                                    consecutive_passes::Int,
                                    m::Int,
                                    m_max::Int)
    if m <= _fourier_guard_through_m(convergence, m_max)
        _update_fourier_snapshots!(snapshots, outputs)
        return false, 0
    end
    passed = _fourier_stokes_passes(convergence, outputs, snapshots)
    _update_fourier_snapshots!(snapshots, outputs)
    consecutive_passes = passed ? consecutive_passes + 1 : 0
    return consecutive_passes >= convergence.n_consecutive, consecutive_passes
end

"Update the consecutive-pass counter and report whether the loop may stop."
function _fourier_convergence_step!(convergence::IntensityConvergence,
                                    outputs::Tuple,
                                    snapshots::Tuple,
                                    consecutive_passes::Int,
                                    m::Int,
                                    m_max::Int)
    if m <= _fourier_guard_through_m(convergence, m_max)
        _update_fourier_snapshots!(snapshots, outputs)
        return false, 0
    end
    passed = _fourier_intensity_passes(convergence, outputs, snapshots)
    _update_fourier_snapshots!(snapshots, outputs)
    consecutive_passes = passed ? consecutive_passes + 1 : 0
    return consecutive_passes >= convergence.n_consecutive, consecutive_passes
end
