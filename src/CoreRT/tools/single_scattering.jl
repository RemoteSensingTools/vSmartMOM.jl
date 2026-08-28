# =========================================================================
# Exact single-scattering engine — shared by the TMS correction and the
# standalone exact-SS reference (`rt_run_ss_exact`).
# =========================================================================
#
# SINGLE-DEFINITION POLICY (docs/dev_notes/forward_lin_parity_policy.md):
# every quantity here is computed in exactly one place and consumed by both
# the production TMS path and the exact-SS reference. Do NOT copy any of
# these formulas into solver code; call these functions.
#
# Design (Nakajima & Tanaka 1988, "Option B — never generate, then add"):
# the truncated direct-beam source is zeroed at the zero-weight view nodes
# in every Fourier mode (`_view_node_mask!`, applied right after the
# elemental step, before doubling — the beam term is the layer's initial
# j₀ at those rows, and zero-weight rows never feed other rows through the
# weighted operators, so the edit is strictly local to the view outputs).
# The exact single scattering is then added in real space after the
# Fourier sum (`tms_correction!`) — the exact kernel is never Fourier
# expanded (its azimuth series converges far too slowly for peaked phase
# functions; avoiding that is the entire point).
#
# N–T scaling without a materialized f: the memo's ω/(1−ωf) paired with
# δ-scaled Δτ′ is algebraically β_scat/Δτ′ (because Δτ′ = Δτ·(1−ωf)), so
# the accumulation consumes the solver's own scaled cumulative τ
# (`τ_sum_all`) plus UNSCALED scattering coefficients — no separate
# per-layer effective-f bookkeeping that could drift from the optics.

"""
    _ss_path_factor(τ′_above, Δτ′, μv, μ₀)

Path-length factor of the first-order source integrated through one layer,
observed at TOA: `exp(−τ′_above·(1/μv+1/μ₀)) · μ₀/(μv+μ₀) ·
(1 − exp(−Δτ′·(1/μv+1/μ₀)))`. Uses `expm1` so the Δτ′ → 0 limit is exact
to machine precision in Float32 and for ForwardDiff Duals (no Taylor
branch needed).
"""
@inline function _ss_path_factor(τ′_above, Δτ′, μv, μ₀)
    s = inv(μv) + inv(μ₀)
    return exp(-τ′_above * s) * (μ₀ / (μv + μ₀)) * (-expm1(-Δτ′ * s))
end

"""
    _exact_phase_columns(greek, phase_ν, phase_greek, μ₀, μv, Δϕ, ν_spec, n)
        -> Matrix (n × nSpec)

First column of the exact phase matrix `Z(μv, μ₀, Δϕ)` per spectral point,
for an unpolarized incident beam (the incident-side meridian rotation is
the identity on `[1,0,0,0]`, so the first column fully determines the SS
source). Evaluated once per stored spectral node via
`Scattering.phase_matrix_first_column` (which owns the scattering-angle
and rotation conventions — single definition) and interpolated linearly
in wavenumber, exactly mirroring how the truncated `phase_greek` nodes are
consumed. When `phase_ν === nothing` (single-node optics) the column is
constant over the band.
"""
function _exact_phase_columns(greek::Scattering.GreekCoefs,
                              phase_ν, phase_greek,
                              μ₀::FT, μv::FT, Δϕ::FT,
                              ν_spec::AbstractVector, n::Int) where {FT}
    cols = zeros(FT, n, length(ν_spec))
    if phase_ν === nothing || phase_greek === nothing || length(phase_ν) == 1
        c = Scattering.phase_matrix_first_column(greek, μ₀, μv, Δϕ, Val(n))
        for j in 1:n
            cols[j, :] .= FT(c[j])
        end
        return cols
    end
    node_cols = [Scattering.phase_matrix_first_column(g, μ₀, μv, Δϕ, Val(n))
                 for g in phase_greek]
    for (is, ν) in enumerate(ν_spec)
        # searchsortedfirst == findfirst(>=(ν)) on a sorted axis, without
        # allocating a closure per spectral point (this ran nSpec x nviews x
        # n_aer times).
        hi = min(searchsortedfirst(phase_ν, ν), length(phase_ν))
        hi = max(hi, 2)
        lo = hi - 1
        t = clamp((ν - phase_ν[lo]) / (phase_ν[hi] - phase_ν[lo]),
                  zero(FT), one(FT))
        for j in 1:n
            cols[j, is] = (1 - t) * FT(node_cols[lo][j]) + t * FT(node_cols[hi][j])
        end
    end
    return cols
end

"""
    _exact_ss_accumulate!(ΔI, model, iBand, τ_above, Δτ_layer; untruncated)

Core accumulation shared by [`tms_correction!`](@ref) and
[`rt_run_ss_exact`](@ref):

    ΔI[iv, :, s] += (F₀ / 4π) · Σ_layers path_factor(τ_above, Δτ, μv, μ₀) ·
                    Σ_scatterers β_scat(s, iz) · Zcol_scat(geom, s) / Δτ(s, iz)

- `τ_above`, `Δτ_layer` are `(nSpec, nZ)` cumulative-above and per-layer
  optical depths in whatever scaling the caller works in. TMS passes the
  solver's δ-SCALED depths (from `τ_sum_all`) — that IS the N–T scaling,
  since β_scat stays unscaled. The exact-SS reference passes UNSCALED
  depths (`f = 0` semantics).
- `untruncated = true` evaluates aerosol phase columns from the retained
  raw Greek sets (`phase_greek_raw`), falling back to the truncated set
  when no truncation occurred (exact == truncated by construction).
- Rayleigh is always exact (its expansion is complete at l = 2).
- Unpolarized solar source: callers must have verified `F₀[2:end] == 0`.
"""
function _exact_ss_accumulate!(ΔI::AbstractArray{FT,3}, model, iBand,
                               τ_above::AbstractMatrix, Δτ_layer::AbstractMatrix;
                               untruncated::Bool = true,
                               F₀_I::AbstractVector = ones(FT, size(τ_above, 1)),
                               view_selected::Union{Nothing,AbstractVector{Bool}} = nothing) where {FT}
    iB = iBand isa Integer ? Int(iBand) : Int(first(iBand))
    (; vza, vaz) = model.obs_geom
    # Exact μ₀ from the quadrature construction — the SAME value the beam
    # attenuation used (single-definition policy).
    μ₀ = FT(model.quad_points.μ₀)
    pol_n = polarization_type(model).n
    nSpec, nZ = size(Δτ_layer)
    ν_spec = model.atmosphere.spec_bands[iB]
    gr_src = model.optics.rayleigh.greek_rayleigh
    gr_ray = gr_src isa AbstractVector ? gr_src[iB] : gr_src

    τ_rayl = model.τ_rayl[iB]                 # (nSpec, nZ), scattering ≡ extinction
    τ_aer  = model.τ_aer[iB]                  # (nAer, nSpec, nZ), UNSCALED
    n_aer  = size(τ_aer, 1)
    aer_optics = model.optics.aerosols.aerosol_optics[iB]

    nS_rayl = size(τ_rayl, 1)
    nS_aer  = size(τ_aer, 2)

    # SPARSITY. In the layer-resolved aerosol injection every active LAYER is
    # registered as its own aerosol "mode", so `τ_aer` is diagonal in
    # (mode, layer): mode ia is non-zero in exactly one iz. Walking the full
    # (iz, s, ia) cross-product therefore spends ~98.6% of its iterations
    # multiplying a structural zero — at C90/OCO scale that is 319 M inner
    # iterations where 4.4 M suffice. The pattern depends only on τ_aer, not on
    # geometry, so it is built once here and reused for every view.
    aer_active = falses(n_aer, nZ)
    @inbounds for iz in 1:nZ, ia in 1:n_aer
        for s in 1:nS_aer
            if !iszero(τ_aer[ia, s, iz])
                aer_active[ia, iz] = true
                break
            end
        end
    end
    scale_s = Vector{FT}(undef, nSpec)

    for iv in eachindex(vza)
        view_selected !== nothing && !view_selected[iv] && continue
        μv = FT(cosd(vza[iv]))
        Δϕ = FT(deg2rad(vaz[iv]))

        # Exact phase columns at this geometry (n × nSpec), per scatterer.
        ray_col = _exact_phase_columns(gr_ray, nothing, nothing,
                                       μ₀, μv, Δϕ, ν_spec, pol_n)
        aer_cols = Vector{Matrix{FT}}(undef, n_aer)
        for ia in 1:n_aer
            ao = aer_optics[ia]
            g_nodes = untruncated && ao.phase_greek_raw !== nothing ?
                      ao.phase_greek_raw :
                      (ao.phase_greek !== nothing ? ao.phase_greek : nothing)
            g_single = untruncated && ao.phase_greek_raw !== nothing ?
                       ao.phase_greek_raw[1] : ao.greek_coefs
            aer_cols[ia] = _exact_phase_columns(g_single, ao.phase_ν, g_nodes,
                                                μ₀, μv, Δϕ, ν_spec, pol_n)
        end

        @inbounds for iz in 1:nZ
            # Per-(iz, s) geometry factor, computed ONCE and shared by the
            # Rayleigh and aerosol passes below (it used to be recomputed
            # inside a fused (iz, s) loop).
            for s in 1:nSpec
                Δτ = FT(Δτ_layer[s, iz])
                if Δτ <= zero(FT)
                    scale_s[s] = zero(FT)
                else
                    pf = _ss_path_factor(FT(τ_above[s, iz]), Δτ, μv, μ₀)
                    scale_s[s] = F₀_I[s] * pf / (Δτ * FT(4π))
                end
            end

            # β_rayleigh = τ_rayl (ω_ray ≡ 1).
            for s in 1:nSpec
                sc = scale_s[s]
                iszero(sc) && continue
                βr = FT(τ_rayl[min(s, nS_rayl), iz])
                for j in 1:pol_n
                    ΔI[iv, j, s] += sc * βr * ray_col[j, s]
                end
            end

            # β_aer = ω̃·τ_aer (unscaled), only for scatterers that actually
            # live in this layer — see `aer_active` above.
            for ia in 1:n_aer
                aer_active[ia, iz] || continue
                ω̃ = aer_optics[ia].ω̃
                ω̃_arr = ω̃ isa AbstractArray        # hoisted: was tested per (s, iz)
                ω̃_const = ω̃_arr ? zero(FT) : FT(ω̃)
                col = aer_cols[ia]
                for s in 1:nSpec
                    sc = scale_s[s]
                    iszero(sc) && continue
                    ω̃s = ω̃_arr ? FT(ω̃[s]) : ω̃_const
                    β = ω̃s * FT(τ_aer[ia, min(s, nS_aer), iz])
                    iszero(β) && continue
                    for j in 1:pol_n
                        ΔI[iv, j, s] += sc * β * col[j, s]
                    end
                end
            end
        end
    end
    return ΔI
end

"""
    _require_unpolarized_solar(F₀, what)

TMS/exact-SS currently assume an unpolarized incident beam (the exact
phase evaluation uses only the first Z column). Reject anything else
loudly rather than returning a silently incomplete polarized correction.
"""
function _require_unpolarized_solar(F₀::AbstractArray, what::String)
    size(F₀, 1) <= 1 && return nothing
    all(iszero, @view F₀[2:end, :]) || throw(ArgumentError(
        "$what requires an unpolarized solar source (F₀ Stokes components " *
        "Q/U/V must be zero); polarized-incident single scattering is not " *
        "implemented"))
    return nothing
end

"""
    _view_node_beam_mask(quad_points, arr_type, FT) -> device vector (NquadN)

Multiplicative mask that zeroes the direct-beam source at zero-weight
nodes: `0` where `wt_μN == 0`, `1` elsewhere. Zero-weight nodes are the
appended view directions (and, in the legacy embedded representation, the
solar node — whose beam source is only ever read back at that node's own
output, i.e. when it IS a view, so masking it is equally correct). True
quadrature nodes keep their beam term — they feed the diffuse field.
"""
function _view_node_beam_mask(quad_points, arr_type, ::Type{FT}) where {FT}
    m = FT.(.!iszero.(Array(quad_points.wt_μN)))
    return arr_type(m)
end

"""
    _mask_view_node_beam!(added_layer, mask)

Zero the upwelling beam source `j₀⁻` at masked (view) rows of a freshly
built elemental layer, BEFORE doubling: at that point `j₀⁻` at a view row
is exactly the truncated single-scattering source of the layer, and —
because zero-weight rows never feed other rows through the weighted
operators — removing it affects only the view outputs, leaving the
diffuse field bit-identical. The downwelling `j₀⁺` is left untouched: the
TMS correction currently replaces only the upwelling (TOA) single
scattering; BOA transmittance keeps the historical truncated SS.
"""
function _mask_view_node_beam!(added_layer, mask::AbstractVector)
    added_layer.j₀⁻ .*= reshape(mask, :, 1, 1)
    return nothing
end

"""
    _tms_correctable_views(model) -> Vector{Bool}

`true` for each view whose output node is a ZERO-WEIGHT quadrature node —
the nodes the beam mask acted on. A view whose zenith coincides with a
true (weighted) quadrature node cannot be TMS-corrected: its beam term
feeds the diffuse field and was deliberately left unmasked, so its
single scattering remains the truncated one (adding the exact SS there
would double count). Uses the SAME `nearest_point(qp_μ, cosd(vza))`
mapping as `postprocessing_vza!` (single-definition policy).
"""
function _tms_correctable_views(model)
    (; qp_μ, wt_μ) = model.quad_points
    qp_host = Array(qp_μ)
    wt_host = Array(wt_μ)
    sel = [iszero(wt_host[nearest_point(qp_host, cosd(v))])
           for v in model.obs_geom.vza]
    if !all(sel)
        @warn("TMSCorrection: $(count(!, sel)) of $(length(sel)) view zenith " *
              "angles coincide with weighted quadrature nodes; their beam " *
              "source feeds the diffuse field and cannot be masked, so their " *
              "single scattering remains truncated (no correction applied " *
              "there).", maxlog = 1)
    end
    return sel
end

"""
    tms_correction!(R_SFI, model, iBand, τ_sum_all, F₀)

Post-Fourier-sum exact single-scattering addition of the TMS correction
(the counterpart of the view-node beam mask; applied only at views whose
output node was actually masked — see [`_tms_correctable_views`](@ref)).
`τ_sum_all` is the solver's δ-SCALED cumulative optical depth,
`(nSpec, nZ+1)` — using it here is what implements the Nakajima–Tanaka
`ω/(1−ωf)` scaling (see the module header). Modifies
`R_SFI :: (nvza, nStokes, nSpec)` in place.
"""
function tms_correction!(R_SFI::AbstractArray{FT,3}, model, iBand,
                         τ_sum_all::AbstractMatrix, F₀::AbstractArray) where {FT}
    _require_unpolarized_solar(F₀, "TMSCorrection")
    nZ = size(τ_sum_all, 2) - 1
    τ_above = @view τ_sum_all[:, 1:nZ]
    Δτ = τ_sum_all[:, 2:end] .- τ_sum_all[:, 1:nZ]
    F₀_I = FT.(@view F₀[1, :])
    _exact_ss_accumulate!(R_SFI, model, iBand, τ_above, Δτ;
                          untruncated = true, F₀_I = F₀_I,
                          view_selected = _tms_correctable_views(model))
    return R_SFI
end

"""
    rt_run_ss_exact(model; i_band=1) -> Array (nvza × nStokes × nSpec)

Standalone exact (untruncated) single-scattering TOA radiance — the
analytic first-order reference. Uses the same engine as the TMS
correction with `f = 0` semantics: UNSCALED optical depths and the
untruncated phase functions, no Fourier loop at all
(O(nlayers × ngeom × nSpec)). This is the exact-SS reference the
Fourier-space `rt_run_ss` (which runs on truncated per-mode optics)
cannot be: the azimuth expansion of a peaked kernel converges too slowly.
"""
function rt_run_ss_exact(model; i_band::Integer = 1)
    FT = float_type(model)
    iB = Int(i_band)
    (; vza) = model.obs_geom
    pol_n = polarization_type(model).n
    nSpec = size(model.τ_abs[iB], 1)
    nZ    = size(model.τ_abs[iB], 2)

    # UNSCALED cumulative optical depth: gas + Rayleigh + raw aerosol.
    Δτ = Matrix{FT}(undef, nSpec, nZ)
    for iz in 1:nZ, s in 1:nSpec
        τa = zero(FT)
        for ia in 1:size(model.τ_aer[iB], 1)
            τa += FT(model.τ_aer[iB][ia, min(s, size(model.τ_aer[iB], 2)), iz])
        end
        Δτ[s, iz] = FT(model.τ_abs[iB][s, iz]) +
                    FT(model.τ_rayl[iB][min(s, size(model.τ_rayl[iB], 1)), iz]) + τa
    end
    τ_above = hcat(zeros(FT, nSpec), cumsum(Δτ, dims = 2)[:, 1:end-1])

    ΔI = zeros(FT, length(vza), pol_n, nSpec)
    _exact_ss_accumulate!(ΔI, model, iB, τ_above, Δτ; untruncated = true)
    return ΔI
end
