# ==============================================================================
# Azimuthal Fourier-loop convergence — the VLIDORT accuracy test.
#
# See AbstractFourierConvergence / IntensityConvergence in types.jl for the
# strategy contract and LIDORT_CONVERGE (lidort_intensity.f90, LIDORT 3.7)
# for the reference semantics this mirrors:
#
#   moment m passes  ⇔  |ΔIₘ| ≤ tolerance·|I_accumulated|  at EVERY tested
#   output (zero contributions pass); n_consecutive passing moments end the
#   loop; a failing moment resets the counter (VLIDORT's TESTCONV).
# ==============================================================================

"""
Diagnostic: the highest Fourier moment actually computed by the most recent
`rt_run` on this task (VLIDORT's `FOURIER_SAVED`). Equals the stream-derived
`m_max` when the full series ran ([`AllFourierMoments`](@ref) or no early
convergence); smaller when [`IntensityConvergence`](@ref) terminated early.
"""
const _LAST_FOURIER_M_USED = Ref{Int}(-1)

# One tested output combination: pass ⇔ |Δ| ≤ tol·|new| (LIDORT's ACCUR test,
# with TAZM == 0 ⇒ pass because the inequality then reads 0 ≤ tol·|new|).
@inline _fc_pass(Δ, new, tol) = abs(Δ) ≤ tol * abs(new)

"""
    _fourier_full_passes(c, R_SFI, R_prev, T_SFI, T_prev) -> Bool

Monolithic-path test: the accumulators R_SFI/T_SFI (host arrays, shape
`(nVZA, pol_n, nSpec)`) already carry this moment's azimuth-weighted
contribution; the previous-moment snapshots give ΔIₘ by difference. Tests
Stokes-I (slot 1) at every view angle and spectral point, both directions —
LIDORT's "ALL directions AND ALL stream values AND ALL azimuths".
"""
function _fourier_full_passes(c::IntensityConvergence, R_SFI, R_prev, T_SFI, T_prev)
    tol = c.tolerance
    for (A, A0) in ((R_SFI, R_prev), (T_SFI, T_prev))
        @inbounds for s in axes(A, 3), v in axes(A, 1)
            _fc_pass(A[v, 1, s] - A0[v, 1, s], A[v, 1, s], tol) || return false
        end
    end
    return true
end

"""
    _fourier_proxy_passes!(c, I_acc, J₀⁻, vza, vaz, qp_μ, pol_type, m, weight) -> Bool

Atmosphere-only-path test (`stop_after_atmosphere = true`, i.e.
`rt_run_atmosphere`): postprocessing never runs there, so this moment's
TOA Stokes-I contribution at the user view angles is extracted directly from
the composite layer's upwelling source `J₀⁻` with the SAME nearest-stream
lookup and `weight·cos(m·Δφ)` azimuth factor postprocessing would apply
(`_precompute_vza_weights` IS the production helper). The contribution is
accumulated into the host buffer `I_acc (nVZA, nSpec)` and tested against it.
Exact for Lambertian-family surfaces (zero m > 0 surface term); see the
[`IntensityConvergence`](@ref) caveat for structured BRDFs.
"""
function _fourier_proxy_passes!(c::IntensityConvergence, I_acc, J₀⁻,
                                vza, vaz, qp_μ, pol_type, m, weight)
    tol = c.tolerance
    vza_info = _precompute_vza_weights(vza, vaz, qp_μ, pol_type, m, weight)
    pass = true
    @inbounds for i in eachindex(vza)
        istart, _, w = vza_info[i]
        # Stokes-I row of this view's stream block; w's I-diagonal entry is
        # weight·cos(m·Δφ) for both scalar and vector polarization types.
        wI = pol_type.n == 1 ? w : w.diag[1]
        Ji = Array(@view J₀⁻[istart, 1, :])          # (nSpec,) host copy
        for s in eachindex(Ji)
            Δ = wI * Ji[s]
            new = I_acc[i, s] + Δ
            I_acc[i, s] = new
            pass &= _fc_pass(Δ, new, tol)
        end
    end
    return pass
end
