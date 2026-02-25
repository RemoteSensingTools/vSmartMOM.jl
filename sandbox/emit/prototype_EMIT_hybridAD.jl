# =================================================================
# EMIT Hybrid AD Prototype — ForwardDiff outer + Analytic RT inner
# =================================================================
#
# This script demonstrates the "hybrid AD" architecture:
#
#   1. Define an explicit state vector x
#   2. Wrap the RT pipeline in f(x) → spectrum
#   3. AD (ForwardDiff) handles the Mie computation — fast & exact
#   4. Analytic linearized RT handles the radiative transfer kernel
#   5. The chain rule connecting them is applied explicitly
#
# Key advantage: ForwardDiff AD Mie (~1.5x forward cost) replaces the
# expensive analytic Mie linearization (~30-70x forward cost), while
# keeping the efficient analytic RT Jacobian for the core solver.
# =================================================================

using Revise
using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.Scattering
using vSmartMOM.InelasticScattering
using Distributions
using Interpolations
using LinearAlgebra
using Printf

# ─────────────────────────────────────────────────────────────────
# Section 1: State vector definition
# ─────────────────────────────────────────────────────────────────

const IX_τ   = 1   # aerosol reference optical depth
const IX_nᵣ  = 2   # real part of refractive index
const IX_nᵢ  = 3   # imaginary part of refractive index
const IX_rₘ  = 4   # mean radius (μm)
const IX_σg  = 5   # geometric std dev of size distribution
const IX_p₀  = 6   # pressure peak of vertical profile (hPa)
const IX_σp  = 7   # pressure width of vertical profile (hPa)
const IX_alb = 8   # surface albedo (Lambertian, shared across bands)
const NX     = 8   # total state vector length

const PARAM_NAMES = ["τ_ref", "nᵣ", "nᵢ", "rₘ", "σ_g", "p₀", "σ_p", "albedo"]

# ─────────────────────────────────────────────────────────────────
# Section 2: AD-to-Analytic conversion helpers
# ─────────────────────────────────────────────────────────────────

"""
    ad_to_lin_aerosol(aero_ad, rₘ, σ_g) -> (AerosolOptics, linAerosolOptics)

Convert ForwardDiff AD Mie output into the analytic linearized format.

AD derivs columns: [rₘ, σ_g, nᵣ, nᵢ]  (from phase_function_autodiff.jl)
Analytic lin rows:  [nᵣ, nᵢ, rₘ, σ_g]  (convention in lin_model_from_parameters.jl)

Chain rule: AD computes d/drₘ and d/dσ_g, but the analytic code expects
d/d(log rₘ) and d/d(log σ_g). Multiply by rₘ and σ_g respectively.
"""
function ad_to_lin_aerosol(aero_ad::AerosolOptics, rₘ, σ_g)
    gc = aero_ad.greek_coefs
    L = length(gc.β)
    D = aero_ad.derivs  # [6L+2, 4], cols = [rₘ, σ_g, nᵣ, nᵢ]

    col_map = [3, 4, 1, 2]       # AD col → lin row reorder
    chain   = [1.0, 1.0, rₘ, σ_g] # chain rule for log-parameterization

    α̇ = zeros(4, L); β̇ = zeros(4, L); γ̇ = zeros(4, L)
    δ̇ = zeros(4, L); ϵ̇ = zeros(4, L); ζ̇ = zeros(4, L)
    ω̃̇ = zeros(4);    k̇ = zeros(4);    ḟᵗ = zeros(4)

    for lin_row in 1:4
        ad_col = col_map[lin_row]
        cf = chain[lin_row]
        α̇[lin_row, :] = D[0*L+1 : 1*L, ad_col] .* cf
        β̇[lin_row, :] = D[1*L+1 : 2*L, ad_col] .* cf
        γ̇[lin_row, :] = D[2*L+1 : 3*L, ad_col] .* cf
        δ̇[lin_row, :] = D[3*L+1 : 4*L, ad_col] .* cf
        ϵ̇[lin_row, :] = D[4*L+1 : 5*L, ad_col] .* cf
        ζ̇[lin_row, :] = D[5*L+1 : 6*L, ad_col] .* cf
        ω̃̇[lin_row]    = D[6*L+1, ad_col] * cf
        k̇[lin_row]     = D[6*L+2, ad_col] * cf
    end

    fwd = AerosolOptics(greek_coefs=gc, ω̃=aero_ad.ω̃, k=aero_ad.k, fᵗ=aero_ad.fᵗ)
    lin = linAerosolOptics(lin_greek_coefs=linGreekCoefs(α̇, β̇, γ̇, δ̇, ϵ̇, ζ̇),
                           ω̃̇=ω̃̇, k̇=k̇, ḟᵗ=ḟᵗ)
    return fwd, lin
end

"""
    interpolate_ad_aerosol(ad0, ad1, λ_grid, rₘ, σ_g)

Linearly interpolate two AD Mie results (at band endpoints) across the band.
Mirrors the interpolation in lin_model_from_parameters.jl.
"""
function interpolate_ad_aerosol(ad0::AerosolOptics, ad1::AerosolOptics,
                                λ_grid, rₘ, σ_g)
    fwd0, lin0 = ad_to_lin_aerosol(ad0, rₘ, σ_g)
    fwd1, lin1 = ad_to_lin_aerosol(ad1, rₘ, σ_g)

    Nl  = length(fwd1.greek_coefs.α)
    Nl_ = length(fwd0.greek_coefs.α)

    α = zeros(Nl); β = zeros(Nl); γ = zeros(Nl)
    δ = zeros(Nl); ϵ = zeros(Nl); ζ = zeros(Nl)
    α̇ = zeros(4, Nl); β̇ = zeros(4, Nl); γ̇ = zeros(4, Nl)
    δ̇ = zeros(4, Nl); ϵ̇ = zeros(4, Nl); ζ̇ = zeros(4, Nl)

    for (arr, f0f, f1f) in [
        (α, fwd0.greek_coefs.α, fwd1.greek_coefs.α),
        (β, fwd0.greek_coefs.β, fwd1.greek_coefs.β),
        (γ, fwd0.greek_coefs.γ, fwd1.greek_coefs.γ),
        (δ, fwd0.greek_coefs.δ, fwd1.greek_coefs.δ),
        (ϵ, fwd0.greek_coefs.ϵ, fwd1.greek_coefs.ϵ),
        (ζ, fwd0.greek_coefs.ζ, fwd1.greek_coefs.ζ)]
        arr[1:Nl_] .= 0.5 .* (f0f[1:Nl_] .+ f1f[1:Nl_])
        if Nl > Nl_; arr[Nl_+1:Nl] .= 0.5 .* f1f[Nl_+1:Nl]; end
    end
    for (darr, d0f, d1f) in [
        (α̇, lin0.lin_greek_coefs.α̇, lin1.lin_greek_coefs.α̇),
        (β̇, lin0.lin_greek_coefs.β̇, lin1.lin_greek_coefs.β̇),
        (γ̇, lin0.lin_greek_coefs.γ̇, lin1.lin_greek_coefs.γ̇),
        (δ̇, lin0.lin_greek_coefs.δ̇, lin1.lin_greek_coefs.δ̇),
        (ϵ̇, lin0.lin_greek_coefs.ϵ̇, lin1.lin_greek_coefs.ϵ̇),
        (ζ̇, lin0.lin_greek_coefs.ζ̇, lin1.lin_greek_coefs.ζ̇)]
        darr[:, 1:Nl_] .= 0.5 .* (d0f[:, 1:Nl_] .+ d1f[:, 1:Nl_])
        if Nl > Nl_; darr[:, Nl_+1:Nl] .= 0.5 .* d1f[:, Nl_+1:Nl]; end
    end

    ν_grid = [1e4 / λ_grid[1], 1e4 / λ_grid[end]]
    itp_kext = LinearInterpolation(ν_grid, [fwd0.k, fwd1.k])
    itp_ksca = LinearInterpolation(ν_grid, [fwd0.k * fwd0.ω̃, fwd1.k * fwd1.ω̃])

    nλ = length(λ_grid)
    k_vec = zeros(nλ); ω̃_vec = zeros(nλ); fᵗ_vec = zeros(nλ)
    for i in 1:nλ
        ν = 1e4 / λ_grid[i]
        k_vec[i] = itp_kext(ν)
        ω̃_vec[i] = itp_ksca(ν) / k_vec[i]
    end

    k̇_mat = zeros(4, nλ); ω̃̇_mat = zeros(4, nλ); ḟᵗ_mat = zeros(4, nλ)
    for ctr in 1:4
        dk_grid = [lin0.k̇[ctr], lin1.k̇[ctr]]
        dksca_grid = [lin0.k̇[ctr] * fwd0.ω̃ + fwd0.k * lin0.ω̃̇[ctr],
                      lin1.k̇[ctr] * fwd1.ω̃ + fwd1.k * lin1.ω̃̇[ctr]]
        itp_dk    = LinearInterpolation(ν_grid, dk_grid)
        itp_dksca = LinearInterpolation(ν_grid, dksca_grid)
        for i in 1:nλ
            ν = 1e4 / λ_grid[i]
            k̇_mat[ctr, i] = itp_dk(ν)
            ω̃̇_mat[ctr, i] = (itp_dksca(ν) - ω̃_vec[i] * k̇_mat[ctr, i]) / k_vec[i]
        end
    end

    fwd_out = AerosolOptics(greek_coefs=GreekCoefs(α, β, γ, δ, ϵ, ζ),
                            ω̃=ω̃_vec, k=k_vec, fᵗ=fᵗ_vec)
    lin_out = linAerosolOptics(lin_greek_coefs=linGreekCoefs(α̇, β̇, γ̇, δ̇, ϵ̇, ζ̇),
                               ω̃̇=ω̃̇_mat, k̇=k̇_mat, ḟᵗ=ḟᵗ_mat)
    return fwd_out, lin_out
end

# ─────────────────────────────────────────────────────────────────
# Section 3: Precomputed base — build once, reuse for all x
# ─────────────────────────────────────────────────────────────────

"""
    PrecomputedBase

Stores everything that does NOT depend on the aerosol/surface state vector:
gas absorption, Rayleigh, atmosphere profile, quadrature points, truncation.

The gas absorption computation is expensive (~30s per band); by precomputing
it once, each subsequent `update_model_from_x!` call takes only ~1s.
"""
struct PrecomputedBase{M, L}
    model::M
    lin_model::L
    params
    NAer::Int
    NGas::Int
    NSurf::Int
end

"""
    precompute_base(params) -> PrecomputedBase

Build the full analytic linearized model once. The expensive gas absorption
profiles are computed here and reused for all state vector evaluations.
"""
function precompute_base(params)
    p = deepcopy(params)
    model, lin_model = model_from_parameters(LinMode(), p)
    NAer  = length(p.scattering_params.rt_aerosols)
    NGas  = size(lin_model.τ̇_abs[1], 1)
    NSurf = length(p.spec_bands)
    return PrecomputedBase(model, lin_model, p, NAer, NGas, NSurf)
end

# ─────────────────────────────────────────────────────────────────
# Section 4: Update model in-place from state vector x
# ─────────────────────────────────────────────────────────────────

"""
    update_model_from_x!(base, x)

Update the precomputed model in-place with new aerosol/surface state
from x, using ForwardDiff AD for the Mie computation.

Only recomputes: Mie optics, τ_aer, τ̇_aer, surface albedo.
Does NOT recompute: gas absorption, Rayleigh, atmosphere profile.
"""
function update_model_from_x!(base::PrecomputedBase, x)
    FT = Float64
    τ_ref = x[IX_τ]; nᵣ = x[IX_nᵣ]; nᵢ = x[IX_nᵢ]
    rₘ = x[IX_rₘ]; σ_g = x[IX_σg]
    p₀ = x[IX_p₀]; σ_p = x[IX_σp]; alb = x[IX_alb]

    model = base.model
    lin_model = base.lin_model
    params = base.params

    # Update surface albedo
    for ib in 1:base.NSurf
        model.params.brdf[ib] = CoreRT.LambertianSurfaceScalar(alb)
    end

    n_bands = length(params.spec_bands)
    truncation_type = Scattering.δBGE{FT}(params.l_trunc, params.Δ_angle)

    sd = LogNormal(log(rₘ), log(σ_g))
    mie_aerosol = Aerosol(sd, nᵣ, nᵢ)

    # Reference extinction at λ_ref (analytic linearized, for k_ref chain rule)
    mie_model_ref = make_mie_model(params.scattering_params.decomp_type,
        mie_aerosol, params.scattering_params.λ_ref, params.polarization_type,
        truncation_type, params.scattering_params.r_max,
        params.scattering_params.nquad_radius)
    k_ref, k̇_ref = compute_ref_aerosol_extinction(LinMode(), mie_model_ref, FT)

    # Vertical profile + derivatives w.r.t. p₀ and σ_p
    τₚ, dτₚ_dp₀, dτₚ_dσp = CoreRT.getAerosolLayerOptProp(
        LinMode(), FT(1), p₀, σ_p, model.profile.p_half)

    for i_band in 1:n_bands
        curr_band_λ = FT(1e4) ./ params.spec_bands[i_band]

        # ── AD Mie computation ────────────────────────────────
        if length(curr_band_λ) == 1
            λ_mid = (maximum(curr_band_λ) + minimum(curr_band_λ)) / 2
            mie_m = make_mie_model(params.scattering_params.decomp_type,
                mie_aerosol, λ_mid, params.polarization_type, truncation_type,
                params.scattering_params.r_max, params.scattering_params.nquad_radius)
            aero_ad = compute_aerosol_optical_properties(mie_m; autodiff=true)
            aero_raw, lin_raw = ad_to_lin_aerosol(aero_ad, rₘ, σ_g)
        else
            mie_m0 = make_mie_model(params.scattering_params.decomp_type,
                mie_aerosol, curr_band_λ[1], params.polarization_type, truncation_type,
                params.scattering_params.r_max, params.scattering_params.nquad_radius)
            mie_m1 = make_mie_model(params.scattering_params.decomp_type,
                mie_aerosol, curr_band_λ[end], params.polarization_type, truncation_type,
                params.scattering_params.r_max, params.scattering_params.nquad_radius)
            aero_ad0 = compute_aerosol_optical_properties(mie_m0; autodiff=true)
            aero_ad1 = compute_aerosol_optical_properties(mie_m1; autodiff=true)
            aero_raw, lin_raw = interpolate_ad_aerosol(
                aero_ad0, aero_ad1, curr_band_λ, rₘ, σ_g)
        end

        # ── Truncation ────────────────────────────────────────
        if length(aero_raw.greek_coefs.β) > truncation_type.l_max
            model.aerosol_optics[i_band][1], lin_model.lin_aerosol_optics[i_band][1] =
                Scattering.truncate_phase(truncation_type, aero_raw, lin_raw;
                                          reportFit=false)
        else
            model.aerosol_optics[i_band][1] = aero_raw
            lin_model.lin_aerosol_optics[i_band][1] = lin_raw
        end

        # ── τ_aer and 7-column chain-rule derivatives ─────────
        k_band = model.aerosol_optics[i_band][1].k
        k̇_band = lin_model.lin_aerosol_optics[i_band][1].k̇

        # τ_aer = (τ_ref / k_ref) * k(λ) * τₚ
        model.τ_aer[i_band][1, :, :] = (τ_ref / k_ref) .* k_band .* τₚ'

        # d/d(τ_ref)
        lin_model.τ̇_aer[i_band][1, 1, :, :] .= (k_band ./ k_ref) .* τₚ'

        # d/d(nᵣ, nᵢ, rₘ, σ_g) via product rule on k/k_ref
        for ctr in 1:4
            lin_model.τ̇_aer[i_band][1, ctr+1, :, :] =
                ((τ_ref / k_ref) .* k̇_band[ctr, :] .-
                 (τ_ref / k_ref^2) .* k̇_ref[ctr] .* k_band) .* τₚ'
        end

        # d/d(p₀, σ_p) via vertical profile
        lin_model.τ̇_aer[i_band][1, 6, :, :] = (τ_ref / k_ref) .* k_band .* dτₚ_dp₀'
        lin_model.τ̇_aer[i_band][1, 7, :, :] = (τ_ref / k_ref) .* k_band .* dτₚ_dσp'
    end

    return nothing
end

# ─────────────────────────────────────────────────────────────────
# Section 5: f(x) wrapper — state vector → spectrum + Jacobian
# ─────────────────────────────────────────────────────────────────

"""
    hybrid_ad_run(base, x; i_band=1) -> (R, T, dR_dx, dT_dx)

Full forward model + Jacobian w.r.t. state vector x.

Returns:
  R, T      — reflectance/transmittance [nVZA, nStokes, nSpec]
  dR_dx     — Jacobian [NX, nVZA, nStokes, nSpec]
"""
function hybrid_ad_run(base::PrecomputedBase, x; i_band=1)
    update_model_from_x!(base, x)

    RS_type = InelasticScattering.noRS(
        fscattRayl  = [1.0],
        ϖ_Cabannes  = [1.0],
        bandSpecLim = UnitRange{Int64}[],
        iBand       = [1],
        F₀          = zeros(1, 1),
        SIF₀        = zeros(1, 1))

    R, T, Ṙ, Ṫ = CoreRT.rt_run_test(RS_type, base.model, base.lin_model,
                                      base.NAer, base.NGas, base.NSurf, i_band)

    # Remap: internal [aer×7 | gas×NGas | surf×NSurf] → state vector
    nParams, nVZA, nStokes, nSpec = size(Ṙ)
    dR_dx = zeros(NX, nVZA, nStokes, nSpec)
    dT_dx = zeros(NX, nVZA, nStokes, nSpec)

    dR_dx[1:7, :, :, :] = Ṙ[1:7, :, :, :]
    dT_dx[1:7, :, :, :] = Ṫ[1:7, :, :, :]

    # Albedo: sum across all per-band surface parameters
    surf_start = 7 * base.NAer + base.NGas + 1
    for ib in 1:base.NSurf
        dR_dx[IX_alb, :, :, :] .+= Ṙ[surf_start + ib - 1, :, :, :]
        dT_dx[IX_alb, :, :, :] .+= Ṫ[surf_start + ib - 1, :, :, :]
    end

    return R, T, dR_dx, dT_dx
end

"""
    forward_only(base, x; i_band=1) -> R

Forward radiance only. Used for finite-difference validation.
"""
function forward_only(base::PrecomputedBase, x; i_band=1)
    R, _, _, _ = hybrid_ad_run(base, x; i_band=i_band)
    return R
end

# ─────────────────────────────────────────────────────────────────
# Section 6: Finite-difference Jacobian utility
# ─────────────────────────────────────────────────────────────────

"""
    fd_jacobian(f, x0, ε_vec) -> Matrix [nSpec, NX]

Central-difference Jacobian with per-parameter step sizes.
"""
function fd_jacobian(f, x0, ε_vec)
    y0 = vec(f(x0))
    n = length(x0)
    m = length(y0)
    J = zeros(m, n)
    for i in 1:n
        xp = copy(x0); xp[i] += ε_vec[i]
        xm = copy(x0); xm[i] -= ε_vec[i]
        J[:, i] = (vec(f(xp)) .- vec(f(xm))) ./ (2ε_vec[i])
    end
    return J
end

# ═════════════════════════════════════════════════════════════════
# Section 7: Main execution
# ═════════════════════════════════════════════════════════════════

println("=" ^ 72)
println("  EMIT Hybrid AD Prototype")
println("=" ^ 72)

params = parameters_from_yaml(joinpath(@__DIR__, "..", "..",
    "test", "test_parameters", "ParamsEMIT_fast.yaml"))

c_aero = params.scattering_params.rt_aerosols[1]
x0 = Float64[
    c_aero.τ_ref,
    c_aero.aerosol.nᵣ,
    c_aero.aerosol.nᵢ,
    exp(c_aero.aerosol.size_distribution.μ),
    exp(c_aero.aerosol.size_distribution.σ),
    mean(c_aero.profile),
    std(c_aero.profile),
    params.brdf[1].albedo,
]

println("\nState vector x₀:")
for i in 1:NX
    @printf("  x[%d] %-6s = %g\n", i, PARAM_NAMES[i], x0[i])
end

# ═══════════════════════════════════════════════════════════════
# Precompute base (gas absorption, Rayleigh, atmosphere — once)
# ═══════════════════════════════════════════════════════════════
println("\n--- Precomputing base model (gas absorption, etc.) ---")
t_precompute = @elapsed begin
    base = precompute_base(params)
end
@printf("  Base model built in %.1f s\n", t_precompute)

# ═══════════════════════════════════════════════════════════════
# Run hybrid AD for band 1
# ═══════════════════════════════════════════════════════════════
println("\n--- Running hybrid AD for band 1 ---")
t_hybrid = @elapsed begin
    R, T, dR_dx, dT_dx = hybrid_ad_run(base, x0; i_band=1)
end
@printf("  Spectrum shape: (%s), Jacobian shape: (%s)\n",
    join(size(R), ", "), join(size(dR_dx), ", "))
@printf("  Hybrid AD time (Mie AD + RT): %.3f s\n", t_hybrid)

# ═══════════════════════════════════════════════════════════════
# Finite-difference validation (fast — reuses precomputed base)
# ═══════════════════════════════════════════════════════════════
println("\n--- Finite-difference validation (band 1, VZA 1, Stokes I) ---")
f_scalar = x -> forward_only(base, x; i_band=1)[1, 1, :]

ε_vec = [1e-5, 1e-5, 1e-7, 1e-6, 1e-5, 1e-1, 1e-1, 1e-5]
J_ad = dR_dx[1:NX, 1, 1, :]'   # [nSpec, NX]
J_fd = fd_jacobian(f_scalar, x0, ε_vec)

# The linearized RT stores derivs 4,5 w.r.t. log(rₘ) and log(σ_g),
# while FD perturbs rₘ and σ_g directly. Convert FD to log-space:
#   dR/d(log rₘ) = dR/drₘ * rₘ
J_fd[:, IX_rₘ] .*= x0[IX_rₘ]
J_fd[:, IX_σg] .*= x0[IX_σg]

println()
@printf("  %-6s  %12s  %12s  %10s  %s\n", "param", "max|AD|", "max|FD|", "rel_err", "note")
println("  " * "-" ^ 60)
for i in 1:NX
    ad_col = J_ad[:, i]
    fd_col = J_fd[:, i]
    max_ad = maximum(abs.(ad_col))
    max_fd = maximum(abs.(fd_col))
    denom  = max(max_ad, max_fd, 1e-30)
    rel_err = maximum(abs.(ad_col .- fd_col)) / denom
    note = (i == IX_rₘ || i == IX_σg) ? "  (log-param)" : ""
    @printf("  %-6s  %12.4e  %12.4e  %10.2e%s\n", PARAM_NAMES[i], max_ad, max_fd, rel_err, note)
end

# ═══════════════════════════════════════════════════════════════
# Compare with fully analytic linearized
# ═══════════════════════════════════════════════════════════════
println("\n--- Comparison with fully analytic linearized ---")

RS_an = InelasticScattering.noRS(
    fscattRayl=[1.0], ϖ_Cabannes=[1.0], bandSpecLim=UnitRange{Int64}[],
    iBand=[1], F₀=zeros(1,1), SIF₀=zeros(1,1))

# Use a fresh precomputed base (with analytic Mie from model_from_parameters)
base_an = precompute_base(params)
R_an, T_an, Ṙ_an, Ṫ_an = CoreRT.rt_run_test(
    RS_an, base_an.model, base_an.lin_model,
    base_an.NAer, base_an.NGas, base_an.NSurf, 1)

println()
@printf("  %-6s  %12s  %12s  %10s\n", "param", "max|hybrid|", "max|analytic|", "rel_err")
println("  " * "-" ^ 48)
for i in 1:7
    h_col = dR_dx[i, 1, 1, :]
    a_col = Ṙ_an[i, 1, 1, :]
    max_h = maximum(abs.(h_col))
    max_a = maximum(abs.(a_col))
    denom = max(max_h, max_a, 1e-30)
    rel_err = maximum(abs.(h_col .- a_col)) / denom
    @printf("  %-6s  %12.4e  %12.4e  %10.2e\n", PARAM_NAMES[i], max_h, max_a, rel_err)
end

# ═══════════════════════════════════════════════════════════════
# Timing comparison (compiled — Mie + RT only, no gas absorption)
# ═══════════════════════════════════════════════════════════════
println("\n--- Timing comparison (compiled, Mie + RT only) ---")
t2_hybrid = @elapsed hybrid_ad_run(base, x0; i_band=1)

# For analytic: rebuild Mie from scratch (model_from_parameters includes gas abs)
p_cpu = deepcopy(params); p_cpu.architecture = CPU()
t2_analytic_build = @elapsed model_from_parameters(LinMode(), p_cpu)
t2_analytic_rt = @elapsed CoreRT.rt_run_test(
    RS_an, base_an.model, base_an.lin_model,
    base_an.NAer, base_an.NGas, base_an.NSurf, 1)

@printf("  Hybrid AD (AD Mie + RT):       %.3f s\n", t2_hybrid)
@printf("  Analytic  (build + RT):        %.3f s + %.3f s = %.3f s\n",
    t2_analytic_build, t2_analytic_rt, t2_analytic_build + t2_analytic_rt)

if t2_hybrid > 0
    @printf("  Speedup (analytic/hybrid):     %.1fx\n",
        (t2_analytic_build + t2_analytic_rt) / t2_hybrid)
end

println("\n" * "=" ^ 72)
println("  Done.")
println("=" ^ 72)
